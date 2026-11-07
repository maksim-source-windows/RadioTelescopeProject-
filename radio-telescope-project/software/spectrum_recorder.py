"""
ЗАПИСЬ СПЕКТРА ВОДОРОДА
Радиотелескоп "Водород-1420"

Требования: pip install numpy pyrtlsdr matplotlib
Анализ линии: software/hi_analysis.py (тот же, что в auto_scanner.py)
"""

import os
import time

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from rtlsdr import RtlSdr

from hi_analysis import calculate_frequencies, band_slices, fit_baseline, find_line_peak

HI_FREQ_HZ = 1420405752  # частота линии HI

SETTINGS = {
    'center_freq': 1420105752,  # LO на 300 кГц ниже линии (как в auto_scanner.py)
    'sample_rate': 2.4e6,
    'gain': 20,                 # LNA даёт ~32 дБ: начните с малого усиления тюнера
    'fft_size': 8192,
    # 100 усреднений = 0.34 с: только быстрый тест.
    # Для 5 минут нужно ~88000 (num_averages = 300 * sample_rate / fft_size).
    'num_averages': 100,
    'line_center_hz': 300000,
    'line_half_hz': 200e3,
    'ref_offset_hz': 600e3,
    'ref_half_hz': 100e3,
    'skip_bins': 2,
    'smooth_bins': 9,
    'detection_sigma': 6.0,
    'window': 'hann',
    'output_dir': 'recordings',
}


def get_window(size, name='hann'):
    """Окно для спектрального анализа"""
    if name == 'hann':
        return np.hanning(size)
    if name == 'hamming':
        return np.hamming(size)
    if name == 'blackman':
        return np.blackman(size)
    return np.ones(size)


def record_spectrum(sdr, num_averages=None):
    """Запись и усреднение спектра (мощность, не амплитуда)"""
    if num_averages is None:
        num_averages = SETTINGS['num_averages']

    fft_size = SETTINGS['fft_size']
    window = get_window(fft_size, SETTINGS['window'])
    result = np.zeros(fft_size, dtype=np.float64)

    sdr.read_samples(fft_size)  # выбросить первый буфер
    for i in range(num_averages):
        samples = sdr.read_samples(fft_size)
        result += np.abs(np.fft.fftshift(np.fft.fft(samples * window))) ** 2
        if i % 1000 == 0:
            print(f"  Усреднение {i + 1}/{num_averages}")
    return result / num_averages


def save_spectrum(freqs_abs, spectrum, timestamp, baseline=None):
    """Сохранение спектра (npz + png). freqs_abs: абсолютные частоты, Гц."""
    os.makedirs(SETTINGS['output_dir'], exist_ok=True)
    npz_file = f"{SETTINGS['output_dir']}/spectrum_{timestamp}.npz"
    np.savez(npz_file, freqs=freqs_abs, spectrum=spectrum,
             center_freq=SETTINGS['center_freq'],
             sample_rate=SETTINGS['sample_rate'])

    plt.figure(figsize=(12, 6))
    plt.plot(freqs_abs / 1e6, 10 * np.log10(spectrum + 1e-30), label='Спектр')
    if baseline is not None:
        plt.plot(freqs_abs / 1e6, 10 * np.log10(np.maximum(baseline, 1e-30)),
                 'r--', alpha=0.7, label='Базовая линия')
    plt.axvline(x=HI_FREQ_HZ / 1e6, color='green', linestyle='--', alpha=0.5,
                label='Линия HI 1420.4058 МГц')
    plt.xlabel('Частота (МГц)')
    plt.ylabel('Мощность (отн. ед., дБ)')
    plt.title(f'Спектр около 1420.4 МГц\n{timestamp}')
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    png_file = f"{SETTINGS['output_dir']}/spectrum_{timestamp}.png"
    plt.savefig(png_file, dpi=150)
    plt.close()
    print(f"Спектр сохранён: {npz_file} и {png_file}")
    return npz_file, png_file


def analyze(spectrum, freqs_rel):
    """Базовая линия + поиск пика. Возвращает (baseline, peak) или (None, None)."""
    fs, n = SETTINGS['sample_rate'], SETTINGS['fft_size']
    line_idx, ref_idx = band_slices(fs, n, SETTINGS)
    baseline = fit_baseline(freqs_rel, spectrum, ref_idx)
    return baseline, find_line_peak(freqs_rel, spectrum, line_idx, baseline, SETTINGS)


def main():
    print("=== ЗАПИСЬ СПЕКТРА ВОДОРОДА ===")
    print(f"Центральная частота: {SETTINGS['center_freq'] / 1e6:.6f} МГц")
    print(f"Полоса: {SETTINGS['sample_rate'] / 1e6:.1f} МГц, окно: {SETTINGS['window']}")
    print(f"Усиление: {SETTINGS['gain']} дБ")
    t_rec = SETTINGS['num_averages'] * SETTINGS['fft_size'] / SETTINGS['sample_rate']
    print(f"Время накопления: {t_rec:.1f} с")
    if t_rec < 60:
        print("ВНИМАНИЕ: накопление короче минуты, слабую линию это не покажет "
              "(для 5 минут num_averages ~ 88000).")

    sdr = RtlSdr()
    sdr.sample_rate = SETTINGS['sample_rate']
    sdr.center_freq = SETTINGS['center_freq']
    sdr.gain = SETTINGS['gain']

    freqs_rel = calculate_frequencies(SETTINGS['sample_rate'], SETTINGS['fft_size'])
    timestamp = time.strftime("%Y-%m-%d_%H-%M-%S")

    try:
        print(f"\nНачало записи в {timestamp}")
        spectrum = record_spectrum(sdr, SETTINGS['num_averages'])
        baseline, peak = analyze(spectrum, freqs_rel)
        save_spectrum(freqs_rel + SETTINGS['center_freq'], spectrum, timestamp, baseline)

        abs_f = SETTINGS['center_freq'] + peak['peak_freq']
        ratio_db = 10 * np.log10(peak['peak_power'] / peak['baseline_power'])
        print("\nПоиск пика в полосе HI:")
        print(f"  Частота пика: {abs_f / 1e6:.6f} МГц "
              f"(от центра SDR: {peak['peak_freq'] / 1e3:.1f} кГц)")
        print(f"  Скорость относительно покоя на частоте HI: "
              f"{(HI_FREQ_HZ - abs_f) / HI_FREQ_HZ * 299792.458:+.1f} км/с "
              "(без поправки на движение Земли/LSR)")
        print(f"  Пик/база: {ratio_db:.3f} дБ, значимость {peak['sigma']:.1f} σ "
              f"(порог {SETTINGS['detection_sigma']})")
        print("  ЛИНИЯ ОБНАРУЖЕНА" if peak['detected'] else "  Линия не обнаружена")
    except Exception as e:
        print(f"\nОшибка: {e}")
        import traceback
        traceback.print_exc()
    finally:
        sdr.close()
        print("SDR закрыт")


if __name__ == "__main__":
    main()
