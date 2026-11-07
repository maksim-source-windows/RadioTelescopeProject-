"""
АНАЛИЗ СПЕКТРА ЛИНИИ HI (общий для auto_scanner.py и spectrum_recorder.py)
Радиотелескоп "Водород-1420"

Схема: LO сдвинут от линии (линия на +line_center_hz от центра SDR), поэтому
DC-пик RTL-SDR не попадает в полосу поиска. Две опорные полосы по обе стороны
от линии дают линейную базовую линию (наклон полосы пропускания и дрейф усиления).

Статистика: спектр сглаживается скользящим средним (smooth_bins бинов), шум
оценивается по разностям отсчётов, разнесённых на smooth_bins (независимые окна;
не зависит от формы линии и наклона базовой линии). Пик считается обнаруженным
при sigma >= detection_sigma. На чистом шуме максимум по ~150 независимым
отсчётам даёт sigma ~ 3 (проверка на 300 синтетических спектрах: медиана 3.2,
максимум 5.0), поэтому порог 6 даёт ложные срабатывания очень редко.
"""

import numpy as np

REQUIRED_KEYS = ('line_center_hz', 'line_half_hz', 'ref_offset_hz',
                 'ref_half_hz', 'skip_bins', 'smooth_bins', 'detection_sigma')


def calculate_frequencies(fs, n):
    """Частоты FFT-бинов относительно центральной частоты SDR (после fftshift)."""
    return np.fft.fftshift(np.fft.fftfreq(n, 1 / fs))


def band_slices(fs, n, cfg):
    """Индексы бинов полосы линии и опорных полос (после fftshift).

    Возвращает (line_idx, ref_idx). Бросает ValueError, если опорная полоса
    выходит за пределы спектра: линейная базовая линия требует обеих сторон.
    """
    for key in REQUIRED_KEYS:
        if key not in cfg:
            raise KeyError(f"в настройках нет ключа '{key}'")

    freqs = calculate_frequencies(fs, n)
    bin_hz = fs / n
    c = n // 2
    lc = int(round(cfg['line_center_hz'] / bin_hz))
    lh = int(round(cfg['line_half_hz'] / bin_hz))
    line_idx = np.arange(c + lc - lh, c + lc + lh + 1)
    if line_idx[0] < 0 or line_idx[-1] >= n:
        raise ValueError("полоса поиска линии выходит за пределы спектра")

    # DC-пик: исключаем, если он вдруг попал в полосу поиска
    line_idx = line_idx[np.abs(freqs[line_idx]) > cfg['skip_bins'] * bin_hz]

    ro = int(round(cfg['ref_offset_hz'] / bin_hz))
    rh = int(round(cfg['ref_half_hz'] / bin_hz))
    ref_left = np.arange(c + lc - ro - rh, c + lc - ro + rh + 1)
    ref_right = np.arange(c + lc + ro - rh, c + lc + ro + rh + 1)
    if ref_left[0] < 0 or ref_right[-1] >= n:
        raise ValueError(
            "опорная полоса выходит за пределы спектра: уменьшите ref_offset_hz "
            "или сдвиньте line_center_hz")
    ref_idx = np.r_[ref_left, ref_right]
    return line_idx, ref_idx


def fit_baseline(freqs, power, ref_idx):
    """Линейная базовая линия по опорным полосам (МНК), массив длины power."""
    slope, intercept = np.polyfit(freqs[ref_idx], power[ref_idx], 1)
    return slope * freqs + intercept


def find_line_peak(freqs, power, line_idx, baseline, cfg):
    """Поиск пика линии со статистической оценкой.

    Возвращает словарь:
        peak_idx          индекс пика в полном спектре
        peak_freq         частота пика относительно центра SDR, Гц
        peak_power        сглаженная мощность в пике
        baseline_power    базовая линия в точке пика
        excess            сглаженное превышение над базовой линией в пике
        sigma             excess / шум (в единицах σ)
        integrated_power  сумма превышений в бинах, где сглаженный спектр > 2σ
        detected          sigma >= detection_sigma
    """
    k = max(1, int(cfg['smooth_bins']))
    excess = power[line_idx] - baseline[line_idx]
    if len(excess) < 3 * k:
        raise ValueError("полоса поиска слишком узкая для заданного smooth_bins")

    if k > 1:
        sm = np.convolve(excess, np.ones(k) / k, mode='valid')
        centers = line_idx[k // 2: k // 2 + len(sm)]
    else:
        sm = excess
        centers = line_idx

    # Шум по разностям независимых окон (устойчив к форме линии и наклону)
    d = sm[k:] - sm[:-k]
    noise = 1.4826 * np.median(np.abs(d - np.median(d))) / np.sqrt(2)
    noise = max(noise, 1e-30)

    p = int(np.argmax(sm))
    peak_idx = int(centers[p])
    base_at_peak = float(baseline[peak_idx])
    sigma = float(sm[p] / noise)
    significant = sm > 2 * noise
    return {
        'peak_idx': peak_idx,
        'peak_freq': float(freqs[peak_idx]),
        'peak_power': base_at_peak + float(sm[p]),
        'baseline_power': base_at_peak,
        'excess': float(sm[p]),
        'sigma': sigma,
        'integrated_power': float(np.sum(sm[significant])),
        'detected': bool(sigma >= cfg['detection_sigma']),
    }
