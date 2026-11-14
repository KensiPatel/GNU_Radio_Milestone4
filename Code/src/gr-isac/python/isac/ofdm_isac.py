#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright 2026 Group 05 ECE310.
#
# SPDX-License-Identifier: GPL-3.0-or-later
#
"""Signal processing of the OFDM ISAC chain (no GNU Radio code in this file).

ISAC = Integrated Sensing And Communication. One OFDM waveform does two jobs:
  * communication: the QPSK symbols on the subcarriers carry data bits,
  * sensing: the receiver knows what was sent, divides the echo by it, and
    reads distance (range) and speed (Doppler) of a target from the result.

Layout of one frame (all lengths in samples):
    n_symbols OFDM symbols, each = cp_len (cyclic prefix) + fft_size useful samples.
Active subcarriers: n_active of the fft_size, symmetric around DC, DC unused.
Pilots: every pilot_spacing-th active subcarrier carries a known symbol.

Sensing steps (see range_doppler_map):
    1. H[m, k] = Y[m, k] / X[m, k]     echo spectrum divided by sent spectrum
    2. IFFT over k (subcarriers)  -> delay, i.e. range
    3. FFT over m (symbols)       -> Doppler shift, i.e. speed
Theory (used by the tests and experiments):
    range resolution   dR = c / (2 B),            B = n_active * subcarrier spacing
    speed resolution   dv = lambda / (2 M T_sym)  (M symbols of length T_sym)
    max unambiguous range (cyclic prefix)  R_max = c * cp_len / (2 fs)
    max unambiguous speed                  v_max = lambda / (4 T_sym)
"""

import math
from dataclasses import dataclass

import numpy

C0 = 299792458.0  # speed of light, m/s


@dataclass(frozen=True)
class OfdmConfig:
    """All numbers that define the waveform and the processing."""

    fft_size: int = 256
    cp_len: int = 64
    n_symbols: int = 128
    n_active: int = 192
    pilot_spacing: int = 4
    sample_rate: float = 10e6
    carrier_freq: float = 5.8e9

    def __post_init__(self):
        if self.fft_size < 8 or self.fft_size & (self.fft_size - 1):
            raise ValueError("fft_size must be a power of two, at least 8")
        if not 0 < self.cp_len < self.fft_size:
            raise ValueError("cp_len must be between 1 and fft_size - 1")
        if self.n_symbols < 2:
            raise ValueError("n_symbols must be at least 2")
        if self.n_active < 2 or self.n_active % 2 or self.n_active >= self.fft_size:
            raise ValueError("n_active must be even and smaller than fft_size")
        if self.pilot_spacing < 0:
            raise ValueError("pilot_spacing must be 0 (no pilots) or positive")
        if self.sample_rate <= 0 or self.carrier_freq <= 0:
            raise ValueError("sample_rate and carrier_freq must be positive")

    # derived numbers ---------------------------------------------------
    @property
    def symbol_len(self):
        return self.fft_size + self.cp_len

    @property
    def frame_len(self):
        return self.n_symbols * self.symbol_len

    @property
    def symbol_time(self):
        return self.symbol_len / self.sample_rate

    @property
    def subcarrier_spacing(self):
        return self.sample_rate / self.fft_size

    @property
    def bandwidth(self):
        return self.n_active * self.subcarrier_spacing

    @property
    def wavelength(self):
        return C0 / self.carrier_freq

    @property
    def range_resolution(self):
        return C0 / (2.0 * self.bandwidth)

    @property
    def speed_resolution(self):
        return self.wavelength / (2.0 * self.n_symbols * self.symbol_time)

    @property
    def max_range(self):
        return C0 * self.cp_len / (2.0 * self.sample_rate)

    @property
    def max_speed(self):
        return self.wavelength / (4.0 * self.symbol_time)

    # subcarrier bookkeeping ------------------------------------------
    @property
    def active_bins(self):
        """FFT indices of the active subcarriers: +1..+n/2 and -n/2..-1 (DC unused)."""
        half = self.n_active // 2
        return numpy.concatenate((numpy.arange(1, half + 1), numpy.arange(self.fft_size - half, self.fft_size)))

    @property
    def pilot_bins(self):
        """FFT indices of the pilot subcarriers (a subset of the active ones)."""
        if self.pilot_spacing == 0:
            return numpy.zeros(0, dtype=int)
        return self.active_bins[:: self.pilot_spacing]

    @property
    def data_bins(self):
        return numpy.setdiff1d(self.active_bins, self.pilot_bins)

    @property
    def bits_per_frame(self):
        return 2 * self.n_symbols * self.data_bins.size


def _qpsk(bits):
    """Gray-mapped QPSK, unit power. bits has shape (..., 2)."""
    return ((1 - 2 * bits[..., 0]) + 1j * (1 - 2 * bits[..., 1])) / math.sqrt(2.0)


def frame_symbols(cfg, seed, frame_index):
    """Frequency-domain symbols X (n_symbols, fft_size) and the data bits of one frame.

    The random generator is seeded with (seed, frame_index), so a receiver that
    knows the seed can rebuild every frame (needed to compute the bit error rate).
    Returns (X, bits) with bits of shape (n_symbols, n_data_bins, 2).
    """
    rng = numpy.random.default_rng([int(seed), int(frame_index)])
    n_data = cfg.data_bins.size
    bits = rng.integers(0, 2, size=(cfg.n_symbols, n_data, 2))
    X = numpy.zeros((cfg.n_symbols, cfg.fft_size), dtype=numpy.complex128)
    X[:, cfg.data_bins] = _qpsk(bits)
    if cfg.pilot_bins.size:
        X[:, cfg.pilot_bins] = pilot_values(cfg)[None, :]
    return X, bits


def pilot_values(cfg):
    """Known pilot symbols (fixed pseudo-random QPSK, the same in every symbol)."""
    rng = numpy.random.default_rng(12345)
    return _qpsk(rng.integers(0, 2, size=(cfg.pilot_bins.size, 2)))


def modulate(X, cfg):
    """OFDM modulator: IFFT, cyclic prefix, scaling to unit average power.

    X: (n_symbols, fft_size). Returns complex64 array of length frame_len.
    """
    x = numpy.fft.ifft(X, axis=1) * math.sqrt(cfg.fft_size * cfg.fft_size / cfg.n_active)
    x = numpy.concatenate((x[:, -cfg.cp_len:], x), axis=1)
    return x.reshape(-1).astype(numpy.complex64)


def demodulate_symbols(rx, cfg, start=0):
    """Cut a frame into OFDM symbols and take the FFT of the useful part.

    rx: samples; the first symbol starts at index start. Returns (n_symbols, fft_size)
    with the same scaling as frame_symbols() (an undistorted frame gives X back).
    """
    rx = numpy.asarray(rx)
    need = start + cfg.frame_len
    if rx.size < need:
        raise ValueError("not enough samples: need %d, got %d" % (need, rx.size))
    sym = rx[start:need].reshape(cfg.n_symbols, cfg.symbol_len)[:, cfg.cp_len:]
    return numpy.fft.fft(sym, axis=1) * math.sqrt(cfg.n_active) / cfg.fft_size


# ---------------------------------------------------------------------------
# Sensing
# ---------------------------------------------------------------------------

def _window(n, use_window):
    return numpy.hanning(n + 1)[:-1] if use_window else numpy.ones(n)


def range_doppler_map(Y, X, cfg, use_window=True, range_pad=2, doppler_pad=2, latency_samples=0):
    """Complex range-Doppler map from received (Y) and sent (X) symbols.

    Y, X: (n_symbols, fft_size) from demodulate_symbols(). Only the active
    subcarriers are used. Returns array (n_doppler, n_range):
      axis 0: Doppler, fftshifted (index n_doppler // 2 = zero speed),
      axis 1: delay latency_samples .. cp_len samples (range_pad bins per sample), so
              column 0 is range 0. latency_samples is the fixed delay of the channel
              (the cyclic prefix only protects delays up to cp_len in total).
    """
    half = cfg.n_active // 2
    # Channel estimate on the contiguous band -half .. +half. The unused DC subcarrier is a
    # hole in the middle of the band that would spread energy over all ranges, so it is
    # filled from its two neighbours (the channel is smooth over one subcarrier).
    band = numpy.arange(-half, half + 1)
    H = numpy.zeros((cfg.n_symbols, band.size), dtype=numpy.complex128)
    nz = band != 0
    H[:, nz] = Y[:, band[nz] % cfg.fft_size] / X[:, band[nz] % cfg.fft_size]
    lo, hi = H[:, half - 1], H[:, half + 1]
    # amplitude: mean of the neighbours; phase: halfway between them (the delay makes the phase turn)
    H[:, half] = 0.5 * (numpy.abs(lo) + numpy.abs(hi)) * numpy.exp(
        1j * (numpy.angle(lo) + 0.5 * numpy.angle(hi * numpy.conj(lo))))
    H = H * _window(band.size, use_window)[None, :] * _window(cfg.n_symbols, use_window)[:, None]
    n_fast = cfg.fft_size * range_pad
    # IFFT over subcarriers gives the delay profile. Zero padding to n_fast points
    # gives range_pad profile points per sample of delay.
    padded = numpy.zeros((cfg.n_symbols, n_fast), dtype=numpy.complex128)
    padded[:, band % n_fast] = H
    profile = numpy.fft.ifft(padded, axis=1)[:, int(latency_samples) * range_pad: cfg.cp_len * range_pad]
    dop = numpy.fft.fft(profile, n=cfg.n_symbols * doppler_pad, axis=0)
    return numpy.fft.fftshift(dop, axes=0)


def axes(cfg, range_pad=2, doppler_pad=2, latency_samples=0):
    """Range axis (m) and speed axis (m/s, positive = target moving towards us)."""
    n_r = (cfg.cp_len - int(latency_samples)) * range_pad
    ranges = numpy.arange(n_r) * C0 / (2.0 * cfg.sample_rate * range_pad)
    n_d = cfg.n_symbols * doppler_pad
    doppler = (numpy.arange(n_d) - n_d // 2) / (n_d * cfg.symbol_time)  # Hz
    return ranges, doppler * cfg.wavelength / 2.0


def _box_sum(a, half_d, half_r):
    """Sum over a (2*half_d+1) x (2*half_r+1) box. Doppler axis wraps, range axis does not."""
    p = numpy.pad(a, ((half_d, half_d), (0, 0)), mode="wrap")
    p = numpy.pad(p, ((0, 0), (half_r, half_r)), mode="constant")
    c = numpy.zeros((p.shape[0] + 1, p.shape[1] + 1))
    c[1:, 1:] = p.cumsum(axis=0).cumsum(axis=1)
    h, w = 2 * half_d + 1, 2 * half_r + 1
    n_d, n_r = a.shape
    return c[h:h + n_d, w:w + n_r] - c[:n_d, w:w + n_r] - c[h:h + n_d, :n_r] + c[:n_d, :n_r]


def ca_cfar_2d(power, pfa, guard=(4, 4), train=(8, 8)):
    """Cell-averaging CFAR on a power map.

    For every cell the noise level is the mean of the training cells around it
    (a ring outside the guard cells). If the noise cells are independent and
    exponentially distributed (complex Gaussian noise), the threshold factor
    alpha = N * (pfa^(-1/N) - 1) gives the false-alarm probability pfa per cell.
    N is the number of training cells (fewer at the range edges).

    Args:
        power: (n_doppler, n_range) linear power.
        pfa: false-alarm probability per cell.
        guard, train: (Doppler, range) half-widths in cells.
    Returns (detections, noise_level): boolean array and the local noise estimate.
    """
    if not 0.0 < pfa < 1.0:
        raise ValueError("pfa must be between 0 and 1")
    gd, gr_ = guard
    td, tr = train
    hd, hr = gd + td, gr_ + tr
    if 2 * hd + 1 > power.shape[0]:
        raise ValueError("training window is larger than the Doppler axis")
    ones = numpy.ones_like(power)
    total = _box_sum(power, hd, hr) - _box_sum(power, gd, gr_)
    count = _box_sum(ones, hd, hr) - _box_sum(ones, gd, gr_)
    noise = total / numpy.maximum(count, 1.0)
    alpha = count * (pfa ** (-1.0 / numpy.maximum(count, 1.0)) - 1.0)
    return power > alpha * noise, noise


def _parabola(a, b, c):
    """Offset of the peak (in bins, -0.5..0.5) from three values around it."""
    den = a - 2.0 * b + c
    return 0.0 if den == 0.0 else 0.5 * (a - c) / den


def extract_targets(cmap, detections, noise, cfg, range_pad=2, doppler_pad=2):
    """Turn detected cells into a list of targets (range m, speed m/s, SNR dB).

    Neighbouring detected cells belong to one target: only local maxima (in a
    3 x 3 neighbourhood) are kept. The peak position is refined with a parabola
    through the log-magnitude of the peak and its two neighbours.
    """
    mag = numpy.abs(cmap)
    n_d, n_r = mag.shape
    ranges, speeds = axes(cfg, range_pad, doppler_pad, 0)
    bin_r = ranges[1] - ranges[0]
    bin_v = speeds[1] - speeds[0]
    logm = numpy.log(mag + 1e-30)
    out = []
    for i, j in zip(*numpy.nonzero(detections)):
        nb = mag[numpy.ix_([(i - 1) % n_d, i, (i + 1) % n_d], [max(j - 1, 0), j, min(j + 1, n_r - 1)])]
        if mag[i, j] < nb.max():
            continue
        dj = _parabola(logm[i, j - 1], logm[i, j], logm[i, j + 1]) if 0 < j < n_r - 1 else 0.0
        di = _parabola(logm[(i - 1) % n_d, j], logm[i, j], logm[(i + 1) % n_d, j])
        r = (j + dj) * bin_r
        v = speeds[i] + di * bin_v
        snr = 10.0 * math.log10(mag[i, j] ** 2 / max(noise[i, j], 1e-30))
        out.append((float(r), float(v), float(snr)))
    out.sort(key=lambda t: -t[2])
    return out


def sense_frame(rx, tx, cfg, pfa=1e-8, use_window=True, range_pad=2, doppler_pad=2,
                guard=(4, 4), train=(8, 8), latency_samples=0):
    """Complete sensing chain for one frame. Returns (targets, power_map_db).

    rx, tx: frame_len samples each (received echo and the known transmitted signal).
    """
    Y = demodulate_symbols(rx, cfg)
    X = demodulate_symbols(tx, cfg)
    cmap = range_doppler_map(Y, X, cfg, use_window, range_pad, doppler_pad, latency_samples)
    power = numpy.abs(cmap) ** 2
    det, noise = ca_cfar_2d(power, pfa, guard, train)
    targets = extract_targets(cmap, det, noise, cfg, range_pad, doppler_pad)
    return targets, 10.0 * numpy.log10(power + 1e-30)


# ---------------------------------------------------------------------------
# Communication side (the same waveform used as a data link)
# ---------------------------------------------------------------------------

def estimate_timing_cp(rx, cfg, max_offset=None):
    """Find where the first OFDM symbol starts, using the cyclic prefix.

    The cyclic prefix is a copy of the last cp_len samples of the symbol, so
    rx[n] and rx[n + fft_size] are equal inside the prefix. The offset with the
    largest summed correlation is the symbol start. Needs no knowledge of the data.
    """
    rx = numpy.asarray(rx)
    max_offset = cfg.cp_len if max_offset is None else int(max_offset)
    best, best_val = 0, -1.0
    for o in range(max_offset + 1):
        total = 0.0 + 0.0j
        for m in range(cfg.n_symbols - 1):
            a = rx[o + m * cfg.symbol_len: o + m * cfg.symbol_len + cfg.cp_len]
            b = rx[o + m * cfg.symbol_len + cfg.fft_size: o + m * cfg.symbol_len + cfg.fft_size + cfg.cp_len]
            total += numpy.sum(a * numpy.conj(b))
        if abs(total) > best_val:
            best, best_val = o, abs(total)
    return best


def equalize_and_decide(Y, cfg):
    """One-tap equaliser per subcarrier from the pilots, then hard QPSK decisions.

    For every OFDM symbol the channel at the pilots is Y / pilot; it is
    interpolated linearly (real and imaginary part) over all active subcarriers.
    Returns bits (n_symbols, n_data_bins, 2).
    """
    pil = cfg.pilot_bins
    if pil.size < 2:
        raise ValueError("pilots are needed for equalisation")

    def signed(bins):  # FFT index -> frequency index in -fft_size/2 .. fft_size/2
        return numpy.where(bins > cfg.fft_size // 2, bins - cfg.fft_size, bins)

    h_pil = Y[:, pil] / pilot_values(cfg)[None, :]
    order = numpy.argsort(signed(pil))
    f_pil = signed(pil)[order]
    f_act = signed(cfg.active_bins)
    h_full = numpy.zeros_like(Y)
    for m in range(Y.shape[0]):
        hp = h_pil[m, order]
        h_full[m, cfg.active_bins] = numpy.interp(f_act, f_pil, hp.real) + 1j * numpy.interp(f_act, f_pil, hp.imag)
    d = cfg.data_bins
    z = Y[:, d] / h_full[:, d]
    bits = numpy.empty(z.shape + (2,), dtype=int)
    bits[..., 0] = (z.real < 0).astype(int)
    bits[..., 1] = (z.imag < 0).astype(int)
    return bits


def bit_errors(rx, cfg, seed, frame_index, start):
    """Demodulate one received frame and count bit errors against the sent bits."""
    Y = demodulate_symbols(rx, cfg, start)
    bits = equalize_and_decide(Y, cfg)
    _, sent = frame_symbols(cfg, seed, frame_index)
    return int(numpy.sum(bits != sent)), int(sent.size)


# ---------------------------------------------------------------------------
# Channel model as a plain function (the block moving_target_channel_cc does the same, streaming)
# ---------------------------------------------------------------------------

def simulate_channel(x, targets, cfg, noise_power=0.0, rng=None, one_way=False, n_taps=33):
    """Apply delay, Doppler shift, gain and noise to the signal x (a numpy array).

    targets: list of (range_m, speed_mps, gain_db), see moving_target_channel_cc.
    The fixed latency (n_taps - 1) / 2 samples of the delay filter is included.
    """
    x = numpy.asarray(x, dtype=numpy.complex128)
    legs = 1.0 if one_way else 2.0
    idx = numpy.arange(x.size, dtype=numpy.float64)
    y = numpy.zeros(x.size, dtype=numpy.complex128)
    for r, v, g in targets:
        delay = legs * float(r) / C0 * cfg.sample_rate
        n_int = int(math.floor(delay))
        k = numpy.arange(n_taps) - (n_taps - 1) / 2.0 - (delay - n_int)
        taps = numpy.sinc(k) * numpy.hanning(n_taps + 2)[1:-1]
        delayed = numpy.zeros(x.size, dtype=numpy.complex128)
        full = numpy.convolve(x, taps)[: x.size]
        delayed[n_int:] = full[: x.size - n_int]
        f_dop = legs * float(v) * cfg.carrier_freq / C0 / cfg.sample_rate
        phase0 = -2.0 * math.pi * ((cfg.carrier_freq * legs * float(r) / C0) % 1.0)
        y += 10.0 ** (float(g) / 20.0) * delayed * numpy.exp(1j * (2.0 * math.pi * ((f_dop * idx) % 1.0) + phase0))
    if noise_power > 0.0:
        rng = numpy.random.default_rng() if rng is None else rng
        s = math.sqrt(noise_power / 2.0)
        y += s * (rng.standard_normal(x.size) + 1j * rng.standard_normal(x.size))
    return y
