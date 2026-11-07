===============================================
         СОФТОВЫЕ УЛУЧШЕНИЯ
         Радиотелескоп "Водород-1420"
===============================================

===============================================
        1. АВТОСКАНИРОВАНИЕ НЕБА
===============================================

Скрипт auto_scanner.py выполняет полное сканирование неба:
- Проходит по сетке (az, alt)
- Для каждой точки делает FFT, ищет пик линии HI
- Сохраняет полные спектры для повторной обработки
- Вычисляет галактические координаты

ПАРАМЕТРЫ:
- line_half_hz: 200e3 (±200 кГц для поиска пика HI)
- integration_time: 300 сек (5 минут на точку)
- save_spectra: True (сохранение полных спектров)

===============================================
        2. СКРИПТ ЗАПИСИ СПЕКТРА
===============================================

Python скрипт для записи и сохранения спектра.

ТРЕБОВАНИЯ:
- Python 3.8+
- numpy, matplotlib
- pyrtlsdr (или SoapySDR)

УСТАНОВКА:
pip install numpy matplotlib pyrtlsdr

КОД:
-----------------------------------
import numpy as np
import matplotlib.pyplot as plt
from rtlsdr import RtlSdr

sdr = RtlSdr()

# Настройки (LO смещён на -300 кГц для избежания DC-пика)
sdr.sample_rate = 2.4e6  # 2.4 MS/s
sdr.center_freq = 1420105752  # 1420.105752 MHz (HI в +300 кГц от центра)
sdr.gain = 20  # LNA уже даёт ~32 дБ, начинать с 20

# Запись
samples = sdr.read_samples(256*1024)

# FFT
fft_data = np.fft.fft(samples)
fft_data = np.fft.fftshift(fft_data)
power = 10*np.log10(np.abs(fft_data)**2 + 1e-12)

# Частоты
freq = np.fft.fftshift(np.fft.fftfreq(len(power), 1/sdr.sample_rate)) / 1e6
freq += sdr.center_freq/1e6  # МГц

# Сохранение
np.savez('spectrum.npz', freq=freq, power=power)

# График
plt.figure(figsize=(10, 6))
plt.plot(freq, power)
plt.xlabel('Частота (МГц)')
plt.ylabel('Мощность (дБ)')
plt.grid(True)
plt.savefig('spectrum.png')
plt.show()
-----------------------------------

# RtlSdr закрывается после plt.show() (или используйте try/finally)

===============================================
        2. ПОСТРОЕНИЕ КАРТЫ НЕБА
===============================================

Скрипт для сканирования неба и построения тепловой карты.

ТРЕБУЕТСЯ:
- Arduino с сервоприводами
- RTL-SDR
- Python скрипт

АЛГОРИТМ:
1. Навести на точку (az, alt)
2. Измерить сигнал
3. Записать
4. Перейти к следующей точке
5. Построить карту

КОД (основная часть):
-----------------------------------
import os
import sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', 'software'))
import numpy as np
import matplotlib.pyplot as plt

# Импорт сканера с реальным измерением
from auto_scanner import RadioTelescope, SETTINGS, AZ_LIMITS, ALT_LIMITS

def _edges(values):
    """Границы ячеек pcolormesh по центрам."""
    v = np.asarray(values, dtype=float)
    if len(v) == 1:
        return np.array([v[0] - 0.5, v[0] + 0.5])
    mid = (v[:-1] + v[1:]) / 2
    return np.r_[v[0] - (v[1] - v[0]) / 2, mid, v[-1] + (v[-1] - v[-2]) / 2]

# Параметры сканирования
azimuths = list(range(0, 181, 15))  # серво: 0-180°
elevations = [10, 25, 40, 55, 70, 80]

data = np.zeros((len(elevations), len(azimuths)))

# Используем реальный сканер (или NO_SIGNAL для теста)
telescope = RadioTelescope(
    arduino_port=SETTINGS['serial_port'],
    use_sdr=True,
    use_arduino=True
)

try:
    for i, alt in enumerate(elevations):
        for j, az in enumerate(azimuths):
            # measure_point() возвращает (result_dB, ok, timestamp, spectrum_data)
            result, ok, timestamp, spectrum_data = telescope.measure_point(
                az, alt, integration_time=10
            )
            if ok:
                data[i, j] = result
                print(f"Az: {az}°, Alt: {alt}° -> {result:.3f} дБ")
            else:
                print(f"Az: {az}°, Alt: {alt}° -> ошибка наведения")

    # Построить карту с pcolormesh (правильные границы ячеек)
    plt.figure(figsize=(12, 8))
    plt.pcolormesh(_edges(azimuths), _edges(elevations), data,
                   cmap='hot', shading='flat')
    plt.colorbar(label='Пик / база (дБ)')
    plt.xlabel('Азимут серво (°)')
    plt.ylabel('Высота серво (°)')
    plt.title('Карта неба на 1420 МГц (az/alt; небо вращается во время скана)')
    plt.savefig('skymap.png', dpi=150)
    plt.show()
finally:
    telescope.close()
-----------------------------------

===============================================
        3. ИНТЕГРАЦИЯ СИГНАЛА
===============================================

Усреднение для улучшения отношения С/Ш.

МЕТОДЫ:

А) Простое усреднение
-----------------------------------
def average_signal(samples, n=100):
    result = np.zeros(len(samples[0]))
    for i in range(n):
        result += np.abs(np.fft.fft(samples[i]))**2  # мощность, не амплитуда
    return result / n
-----------------------------------

Б) экспоненциальное усреднение (онлайн)
-----------------------------------
alpha = 0.1  # коэффициент
smoothed = 0

for new_sample in stream:
    smoothed = alpha * new_sample + (1-alpha) * smoothed
-----------------------------------

В) Водопад (Waterfall)
- Накапливать FFT во времени
- Показывать как 2D график
- SDRAngel уже это умеет

===============================================
        4. ОНЛАЙН ТРАНСЛЯЦИЯ
===============================================

Отправка данных на сервер для отображения в браузере.

ВАРИАНТЫ:

А) Простой веб-сервер
- Python Flask
- JSON API с текущими данными
- Обновление через JS

Б) ThingSpeak / IoT
- Бесплатный сервис
- Отправка данных
- Графики

В) Discord webhook
- Отправка уведомлений
- Картинки

ПРИМЕР (ThingSpeak):
-----------------------------------
import requests

def send_to_thingspeak(signal):
    url = "https://api.thingspeak.com/update"
    data = {
        'api_key': 'ВАШ_API_KEY',
        'field1': signal
    }
    requests.post(url, data=data)
-----------------------------------

===============================================
        5. АВТООПОВЕЩЕНИЕ
===============================================

При обнаружении сильного сигнала — уведомление.

ТРИГГЕРЫ:
- Превышение порога
- Обнаружение пика
- Изменение фона

УВЕДОМЛЕНИЯ:
- Звук
- Email
- Telegram бот

ПРИМЕР (Telegram):
-----------------------------------
import requests

def notify_telegram(message):
    token = 'BOT_TOKEN'
    chat_id = 'CHAT_ID'
    url = f"https://api.telegram.org/bot{token}/sendMessage"
    requests.post(url, {'chat_id': chat_id, 'text': message})

# Использование
if peak_signal > threshold:
    notify_telegram(f"Обнаружен сигнал: {peak_signal} дБ")
-----------------------------------

===============================================
        6. АРХИВ И БАЗА ДАННЫХ
===============================================

Хранение наблюдений для анализа.

СТРУКТУРА:
observations/
  2024-01-15/
    spectrum_21-30.npz
    notes.txt
    weather.txt
  2024-01-16/
    ...

МЕТАДАННЫЕ:
- Дата и время
- Координаты наведения
- Погода
- Настройки приёмника
- Интеграция

===============================================
        7. КАЛЕНДАРЬ НАБЛЮДЕНИЙ
===============================================

Когда что наблюдать:

| Дата | Событие | Объект |
|------|---------|--------|
| Постоянно | Млечный путь | Водород |
| 2026-2027 | Юпитер | Радиоизлучение (в созвездии Льва) |
| 12-13 авг | Персеиды | Метеоры |
| 17-18 ноя | Леониды | Метеоры |
| 13-14 дек | Геминиды | Метеоры |

АВТОНАПОМИНАНИЕ:
- Скрипт проверяет дату
- Присылает уведомление
- Рассчитывает координаты

===============================================
        8. ИНТЕГРАЦИЯ С АСТРОНОМИЧЕСКИМИ ПРОГРАММАМИ
===============================================

ЭКСПОРТ В:
- KML для Google Earth
- FITS для FITSViewer
- CSV для Excel

ПРИМЕР (CSV — рекомендуется):
-----------------------------------
def save_csv(data, az_values, alt_values, time_tags, galactic_l, galactic_b, filename='observations.csv'):
    """Экспорт данных наблюдений в CSV.

    Галактические координаты нельзя отобразить в Google Earth напрямую —
    используйте CSV для анализа в Excel/Python.
    """
    import csv
    with open(filename, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['az', 'alt', 'signal_dB', 'time_utc', 'galactic_l', 'galactic_b'])
        for i, alt in enumerate(alt_values):
            for j, az in enumerate(az_values):
                l = galactic_l[i, j] if galactic_l is not None else ''
                b = galactic_b[i, j] if galactic_b is not None else ''
                l_str = f"{l:.1f}" if l else ''
                b_str = f"{b:.1f}" if b else ''
                writer.writerow([az, alt, data[i, j], time_tags[i, j], l_str, b_str])
    print(f"Сохранено: {filename}")

# Использование (после scan_sky() из auto_scanner.py):
# save_csv(data, az_values, alt_values, time_tags, galactic_l, galactic_b)
-----------------------------------

ПРИМЕР (KML — только для проверенных координат):
-----------------------------------
def save_kml(az, alt, signal):
    """KML имеет смысл только для реальных географических координат.

    Галактические координаты (l, b) нельзя положить в Google Earth.
    """
    # Пример: если вы преобразовали (l, b) в (lat, lon) через сторонний софт
    lat, lon = 55.7558, 37.6173  # пример: Москва
    kml = f'''<?xml version="1.0"?>
<kml>
  <Placemark>
    <name>Сигнал: {signal:.1f} дБ</name>
    <Point><coordinates>{lon},{lat},0</coordinates></Point>
  </Placemark>
</kml>'''
    with open('observation.kml', 'w') as f:
        f.write(kml)
-----------------------------------

===============================================
            ИТОГОВЫЙ СПИСОК
===============================================

| Улучшение | Сложность | Время |
|-----------|-----------|-------|
| Запись спектра | Простая | 1 час |
| Карта неба | Средняя | 4 часа |
| Интеграция | Простая | 1 час |
| Онлайн трансляция | Средняя | 2 часа |
| Автооповещение | Простая | 1 час |
| База данных | Средняя | 2 часа |
| Календарь | Простая | 1 час |
| Экспорт | Средняя | 2 часа |