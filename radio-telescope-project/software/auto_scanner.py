"""
АВТОСКАНИРОВАНИЕ НЕБА
Радиотелескоп "Водород-1420"

Требования: pip install numpy matplotlib pyrtlsdr pyserial astropy
Arduino: software/arduino_telescope.ino (протокол "SET <az> <alt>" -> "POS <az> <alt>")
Анализ спектра: software/hi_analysis.py

Алгоритм в каждой точке:
1. Интеграция: усреднение мощности FFT (окно Ханна) за integration_time секунд.
2. LO сдвинут на -300 кГц от линии HI, линия лежит на +100...+500 кГц от DC,
   DC-пик RTL-SDR в полосу поиска не попадает.
3. Линейная базовая линия по двум опорным полосам (-400...-200 и +800...+1000 кГц
   от LO), сглаживание, поиск пика, оценка значимости sigma.
4. Пик считается обнаруженным при sigma >= detection_sigma.

Ограничения:
- Угол серво НЕ равен географическому азимуту без калибровки по реальному
  источнику (Солнце, Луна). Галактические координаты считаются в предположении,
  что az/alt серво совпадают с настоящими.
- Скан по сетке az/alt занимает часы, небо за это время поворачивается.
  Карта az/alt это не снимок неба: для анализа используйте галактические
  координаты (l, b) и время каждой точки из файла .npz.
- Файл .npz содержит только обычные массивы (allow_pickle не нужен): карта data/sigma/
  galactic_l/galactic_b/time_tags и спектры spectra_az, spectra_alt, spectra_time,
  spectra_power (N x fft_size), spectra_freqs, spectra_sigma, spectra_detected,
  spectra_duration (накопленные отсчёты, с), spectra_wall_duration (реальное время, с).

Режимы: REAL (SDR + Arduino), TEST (test_mode: 6 точек по 3 с), NO_SIGNAL
(no_signal_test: искусственный шум, не требует SDR, Arduino и pyrtlsdr/pyserial).
"""

import os
import time
from datetime import datetime, timezone

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from hi_analysis import calculate_frequencies, band_slices, fit_baseline, find_line_peak

# ═══════════════════════════════════════════════════
#                    НАСТРОЙКИ
# ═══════════════════════════════════════════════════

SETTINGS = {
    # Arduino
    'serial_port': 'COM3',  # изменить!
    'serial_baud': 9600,

    # SDR
    'center_freq': 1420105752,  # LO на 300 кГц ниже линии HI (1420.405752 МГц)
    'sample_rate': 2.4e6,
    # Усиление тюнера. LNA уже даёт ~32 дБ: начните с 20 дБ и поднимайте,
    # пока шум не поднимется заметно над уровнем шума АЦП, но спектр не забит
    # паразитными линиями (признак перегрузки).
    'gain': 20,
    'fft_size': 8192,

    # Полосы (относительно центра SDR)
    'line_center_hz': 300000,  # линия HI на +300 кГц от LO
    'line_half_hz': 200e3,     # полоса поиска ±200 кГц вокруг линии (±42 км/с)
    'ref_offset_hz': 600e3,    # опорные полосы на ±600 кГц от центра линии
    'ref_half_hz': 100e3,      # половина ширины опорной полосы
    'skip_bins': 2,            # пропуск бинов вокруг DC

    # Статистика обнаружения
    'smooth_bins': 9,          # сглаживание (9 бинов ~ 2.6 кГц ~ 0.6 км/с)
    'detection_sigma': 6.0,    # порог обнаружения (на шуме sigma ~3, максимум за 300 прогонов 5.0)

    # Сканирование (серво: азимут 0-180°, высота 10-80°)
    'az_start': 0,
    'az_end': 181,
    'az_step': 15,
    'alt_values': [10, 25, 40, 55, 70, 80],

    'integration_time': 300,   # секунд на точку (5 минут)

    'output_dir': 'skymaps',

    # Тестовый режим (мало точек для быстрой проверки)
    'test_mode': False,
    # Проверка алгоритма на искусственном шуме (SDR не нужен, на выходе sigma ~ 3)
    'no_signal_test': False,

    # Координаты наблюдателя для галактических координат
    # ИСПРАВИТЬ: замените на фактическое место установки телескопа!
    'observer_lat': 55.7558,
    'observer_lon': 37.6173,
    'observer_height': 156,

    # Сохранение полных спектров для повторной обработки
    'save_spectra': True,
}

AZ_LIMITS = (0, 180)
ALT_LIMITS = (10, 80)
STEP_DELAY_S = 0.02  # delay(20) на градус в прошивке
SIDEREAL_DEG_PER_MIN = 0.25068  # вращение неба
TEST_INTEGRATION_S = 3          # test_mode: секунд на точку
BLOCK_FFTS = 64                 # FFT за одно чтение SDR (быстрее и без пропусков между буферами)


def move_timeout(old, new):
    steps = max(abs(new[0] - old[0]), abs(new[1] - old[1]))
    return steps * STEP_DELAY_S + 1.5


class RadioTelescope:
    def __init__(self, arduino_port=None, use_sdr=None, use_arduino=None):
        """Режимы: REAL (по умолчанию), NO_SIGNAL (SETTINGS['no_signal_test']).

        В NO_SIGNAL SDR и Arduino не подключаются, pyrtlsdr и pyserial не нужны.
        """
        no_signal = SETTINGS['no_signal_test']
        self.use_sdr = (not no_signal) if use_sdr is None else (use_sdr and not no_signal)
        use_arduino = (not no_signal) if use_arduino is None else (use_arduino and not no_signal)
        self.use_arduino = bool(use_arduino and arduino_port)

        self.window = np.hanning(SETTINGS['fft_size'])
        self.freqs = calculate_frequencies(SETTINGS['sample_rate'], SETTINGS['fft_size'])
        self.line_idx, self.ref_idx = band_slices(
            SETTINGS['sample_rate'], SETTINGS['fft_size'], SETTINGS)
        self.pos = (90, 45)  # стартовая позиция прошивки

        self.sdr = None
        if self.use_sdr:
            from rtlsdr import RtlSdr  # импорт здесь: в NO_SIGNAL библиотека не нужна
            print("Подключение к RTL-SDR...")
            self.sdr = RtlSdr()
            self.sdr.sample_rate = SETTINGS['sample_rate']
            self.sdr.center_freq = SETTINGS['center_freq']
            self.sdr.gain = SETTINGS['gain']
        else:
            print("RTL-SDR не используется (режим NO_SIGNAL)" if no_signal
                  else "RTL-SDR не используется")

        self.arduino = None
        if self.use_arduino:
            import serial
            print(f"Подключение к Arduino: {arduino_port}")
            self.arduino = serial.Serial(arduino_port, SETTINGS['serial_baud'], timeout=0.5)
            time.sleep(2)  # перезагрузка Arduino
        else:
            print("Arduino не подключена: наведение пропускается")

    def goto(self, az, alt):
        """Навести и дождаться POS. Возвращает True при подтверждении."""
        if not self.arduino:
            return True
        az = int(max(AZ_LIMITS[0], min(AZ_LIMITS[1], az)))
        alt = int(max(ALT_LIMITS[0], min(ALT_LIMITS[1], alt)))

        self.arduino.reset_input_buffer()
        self.arduino.write(f"SET {az} {alt}\n".encode())

        deadline = time.time() + move_timeout(self.pos, (az, alt))
        while time.time() < deadline:
            line = self.arduino.readline().decode(errors='ignore').strip()
            if line.startswith("POS"):
                self.pos = (az, alt)
                return True
        print("  -> Таймаут ожидания POS")
        return False

    def measure_point(self, az, alt, integration_time):
        """Наведение и интеграция.

        Возвращает (result, ok, timestamp, spectrum_data):
            result        отношение пик/база в дБ (nan при ошибке наведения)
            ok            наведение подтверждено
            timestamp     середина интервала измерения (UTC)
            spectrum_data словарь со спектром и статистикой (None при ошибке)
        """
        ok = self.goto(az, alt)
        if not ok:
            print("  -> Наведение не подтверждено, измерение пропущено")
            return np.nan, False, datetime.now(timezone.utc), None

        if self.sdr is None and not SETTINGS['no_signal_test']:
            raise RuntimeError("SDR не подключён: включите no_signal_test или подключите RTL-SDR")

        time.sleep(1)  # механическая стабилизация
        time_start = datetime.now(timezone.utc)
        t0 = time.monotonic()

        fft_size = SETTINGS['fft_size']
        num = max(1, int(integration_time * SETTINGS['sample_rate'] / fft_size))
        sample_time = num * fft_size / SETTINGS['sample_rate']  # длительность накопленных отсчётов
        acc = np.zeros(fft_size, dtype=np.float64)

        if SETTINGS['no_signal_test']:
            print("    [NO-SIGNAL TEST] искусственный шум")
        else:
            self.sdr.read_samples(fft_size)  # выбросить первый буфер
        done = 0
        next_report = 5000
        while done < num:
            m = min(BLOCK_FFTS, num - done)
            if SETTINGS['no_signal_test']:
                x = (np.random.randn(m, fft_size) + 1j * np.random.randn(m, fft_size)) / np.sqrt(2)
            else:
                x = np.asarray(self.sdr.read_samples(m * fft_size)).reshape(m, fft_size)
            spec = np.fft.fftshift(np.fft.fft(x * self.window, axis=1), axes=1)
            acc += np.sum(np.abs(spec) ** 2, axis=0)
            done += m
            if done >= next_report:
                print(f"    [{done}/{num}] FFT")
                next_report += 5000
        acc /= num

        wall_time = time.monotonic() - t0
        time_end = datetime.now(timezone.utc)
        time_avg = time_start + (time_end - time_start) / 2

        baseline = fit_baseline(self.freqs, acc, self.ref_idx)
        peak = find_line_peak(self.freqs, acc, self.line_idx, baseline, SETTINGS)
        result = 10 * np.log10(peak['peak_power'] / peak['baseline_power'])

        spectrum_data = {
            'freqs': self.freqs + SETTINGS['center_freq'],  # абсолютные частоты, Гц
            'power': acc,
            'baseline': baseline,
            'peak_abs_freq': SETTINGS['center_freq'] + peak['peak_freq'],
            'sample_integration_time': sample_time,  # накопленные отсчёты, с
            'wall_clock_time': wall_time,            # реальное время измерения, с
            **peak,
        }
        return result, ok, time_avg, spectrum_data

    def close(self):
        if self.sdr:
            self.sdr.close()
        if self.arduino:
            self.arduino.close()


def altaz_to_galactic(az, alt, time_obs):
    """Горизонтальные координаты -> галактические (l, b) в градусах.

    Возвращает (None, None), если astropy не установлен. Любая другая ошибка
    не прерывает скан: возвращается (None, None) с сообщением.
    """
    try:
        from astropy.coordinates import SkyCoord, AltAz, Galactic, EarthLocation
        import astropy.units as u
    except ImportError:
        return None, None
    try:
        location = EarthLocation(
            lat=SETTINGS['observer_lat'] * u.deg,
            lon=SETTINGS['observer_lon'] * u.deg,
            height=SETTINGS['observer_height'] * u.m)
        coord = SkyCoord(alt=alt * u.deg, az=az * u.deg,
                         frame=AltAz(obstime=time_obs, location=location))
        galactic = coord.transform_to(Galactic())
        return float(galactic.l.deg), float(galactic.b.deg)
    except Exception as e:  # нет данных IERS/эфемерид и т.п.
        print(f"  (галактические координаты недоступны: {e})")
        return None, None


def scan_grid():
    if SETTINGS['test_mode']:
        print("Тестовый режим: 3 азимута x 2 высоты = 6 точек")
        return [45, 90, 135], [30, 60]
    az_values = list(range(SETTINGS['az_start'], SETTINGS['az_end'], SETTINGS['az_step']))
    return az_values, list(SETTINGS['alt_values'])


def scan_sky(telescope):
    az_values, alt_values = scan_grid()
    integration = TEST_INTEGRATION_S if SETTINGS['test_mode'] else SETTINGS['integration_time']
    shape = (len(alt_values), len(az_values))
    data = np.full(shape, np.nan)
    sigma_map = np.full(shape, np.nan)
    time_tags = np.full(shape, '', dtype='U32')  # строки ISO, без pickle
    galactic_l = np.full(shape, np.nan)
    galactic_b = np.full(shape, np.nan)
    spectra_list = []

    total = len(az_values) * len(alt_values)
    est_min = total * (integration + 4) / 60
    print(f"Сканирование: {total} точек, ~{est_min:.0f} минут")
    if est_min > 60:
        print(f"ВНИМАНИЕ: за время скана небо повернётся на ~{est_min * SIDEREAL_DEG_PER_MIN:.0f}°. "
              "Карта az/alt не снимок неба; используйте (l, b) и время точек.")

    n = 0
    for i, alt in enumerate(alt_values):
        az_iter = list(enumerate(az_values))
        if i % 2:
            az_iter.reverse()
        for j, az in az_iter:
            n += 1
            print(f"[{n}/{total}] Az: {az}°, Alt: {alt}°")

            val, ok, time_avg, sd = telescope.measure_point(az, alt, integration)
            data[i, j] = val if ok else np.nan
            time_tags[i, j] = time_avg.isoformat()

            l, b = altaz_to_galactic(az, alt, time_avg)
            if l is not None:
                galactic_l[i, j], galactic_b[i, j] = l, b

            if sd is None:
                print("  → пропуск (позиция не подтверждена)")
                continue

            sigma_map[i, j] = sd['sigma']
            if SETTINGS['save_spectra']:
                spectra_list.append({'az': az, 'alt': alt,
                                     'time': time_avg.isoformat(), **sd})

            verdict = "ОБНАРУЖЕНА" if sd['detected'] else "не обнаружена"
            print(f"  → {val:.3f} дБ, σ={sd['sigma']:.1f}, линия {verdict}")
            print(f"     τ={sd['sample_integration_time']:.1f} с отсчётов, "
                  f"реально {sd['wall_clock_time']:.1f} с (план {integration} с)")
            if sd['detected']:
                print(f"     пик на {sd['peak_abs_freq'] / 1e6:.6f} МГц")
            if l is not None:
                print(f"     l={l:.1f}°, b={b:.1f}°  время: {time_avg.isoformat()}")
    return data, az_values, alt_values, time_tags, galactic_l, galactic_b, spectra_list, sigma_map


def _edges(values):
    """Границы ячеек pcolormesh по центрам (в т.ч. для неравномерной сетки)."""
    v = np.asarray(values, dtype=float)
    if len(v) == 1:
        return np.array([v[0] - 0.5, v[0] + 0.5])
    mid = (v[:-1] + v[1:]) / 2
    return np.r_[v[0] - (v[1] - v[0]) / 2, mid, v[-1] + (v[-1] - v[-2]) / 2]


def plot_skymap(data, filename, az_values, alt_values):
    """Карта отношения линия/база в координатах servo az/alt (не снимок неба)."""
    plt.figure(figsize=(14, 8))
    plt.pcolormesh(_edges(az_values), _edges(alt_values), data, cmap='hot', shading='flat')
    plt.colorbar(label='Пик / база (дБ)')
    plt.xlabel('Азимут серво (°)')
    plt.ylabel('Высота серво (°)')
    plt.title('Карта на 1420 МГц (относительная, az/alt; небо вращается во время скана)')
    plt.tight_layout()
    plt.savefig(filename, dpi=150)
    plt.close()
    print(f"Карта сохранена: {filename}")


def main():
    os.makedirs(SETTINGS['output_dir'], exist_ok=True)
    timestamp = time.strftime("%Y-%m-%d_%H-%M")
    telescope = RadioTelescope(arduino_port=SETTINGS['serial_port'])
    try:
        (data, az_vals, alt_vals, time_tags, gal_l, gal_b,
         spectra_list, sigma_map) = scan_sky(telescope)

        filename = f"{SETTINGS['output_dir']}/skymap_{timestamp}"
        save_dict = {
            'data': data, 'sigma': sigma_map,
            'az': np.array(az_vals), 'alt': np.array(alt_vals),
            'time_tags': time_tags, 'galactic_l': gal_l, 'galactic_b': gal_b,
        }
        if SETTINGS['save_spectra'] and spectra_list:
            save_dict.update({
                'spectra_az': np.array([r['az'] for r in spectra_list]),
                'spectra_alt': np.array([r['alt'] for r in spectra_list]),
                'spectra_time': np.array([r['time'] for r in spectra_list], dtype='U32'),
                'spectra_power': np.array([r['power'] for r in spectra_list]),
                'spectra_freqs': spectra_list[0]['freqs'],
                'spectra_sigma': np.array([r['sigma'] for r in spectra_list]),
                'spectra_detected': np.array([r['detected'] for r in spectra_list]),
                'spectra_duration': np.array([r['sample_integration_time'] for r in spectra_list]),
                'spectra_wall_duration': np.array([r['wall_clock_time'] for r in spectra_list]),
            })
            print(f"Сохранено {len(spectra_list)} спектров")

        np.savez(filename + '.npz', **save_dict)
        plot_skymap(data, filename + '.png', az_vals, alt_vals)
        print(f"\nГотово! Данные: {filename}.npz (np.load без allow_pickle)")

        valid = data[~np.isnan(data)]
        if len(valid) > 0:
            print(f"Отношение пик/база: среднее {np.mean(valid):.3f} дБ, "
                  f"σ {np.std(valid):.3f} дБ, мин {np.min(valid):.3f}, макс {np.max(valid):.3f}")
        n_det = int(np.sum(sigma_map >= SETTINGS['detection_sigma']))
        print(f"Точек с обнаруженной линией (σ >= {SETTINGS['detection_sigma']}): {n_det}")
    finally:
        telescope.close()


if __name__ == "__main__":
    main()
