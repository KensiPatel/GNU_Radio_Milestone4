#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright 2026 Group 05 ECE310.
#
# SPDX-License-Identifier: GPL-3.0-or-later
#
"""Energy detector for spectrum sensing.

The block turns a complex baseband stream into a power spectrum and decides,
for every frequency bin, whether the band is "occupied" or "free".

Processing per output item (fft_size * avg_len input samples):
  1. cut the input into avg_len frames of fft_size samples,
  2. multiply each frame by a Hann window and take its FFT,
  3. average |FFT|^2 over the avg_len frames (averaging lowers the noise
     variance, which is what makes weak signals detectable),
  4. normalise so that white noise of power sigma^2 gives a mean bin value of
     sigma^2 (the window power is divided out),
  5. compare every bin with a threshold.

New in M4 (all optional, defaults keep the M3 behaviour):
  * pfa: the threshold is calculated from theory for a chosen false-alarm
    probability per bin (CFAR), so no dB value has to be guessed,
  * decision_mode: the band decision can use the total power over all bins,
    which sees wideband signals that no single bin shows,
  * a message output "bands" that lists the occupied frequency bands.
"""

import math
from statistics import NormalDist

import numpy
import pmt
from gnuradio import gr

# Tiny value that keeps log10() finite when a bin is exactly zero.
_EPS = 1e-20


def energy_spectrum(samples, fft_size, avg_len):
    """Averaged, window-normalised power spectrum.

    Args:
        samples: complex array of length n_out * avg_len * fft_size.
        fft_size: FFT length (number of bins).
        avg_len: number of FFT frames averaged for one output.

    Returns:
        float64 array of shape (n_out, fft_size), linear power, fftshifted so
        that index fft_size // 2 is 0 Hz. White noise of power sigma^2 gives a
        mean of sigma^2 in every bin.
    """
    samples = numpy.asarray(samples)
    n_out = samples.size // (fft_size * avg_len)
    frames = samples[: n_out * avg_len * fft_size].reshape(n_out, avg_len, fft_size)
    window = numpy.hanning(fft_size + 1)[:-1]  # periodic Hann
    spectrum = numpy.fft.fft(frames * window, axis=2)
    power = numpy.abs(spectrum) ** 2 / numpy.sum(window ** 2)
    return numpy.fft.fftshift(power.mean(axis=1), axes=1)


def noise_floor_estimate(power, avg_len):
    """Robust noise power estimate from one averaged spectrum.

    Uses the median over bins, corrected for the fact that the median of an
    averaged exponential variable (a Gamma(avg_len) variable) is about
    (avg_len - 1/3) / avg_len times its mean. This works while the signal
    occupies less than half of the bins.
    """
    return numpy.median(power) / ((avg_len - 1.0 / 3.0) / avg_len)


def cfar_threshold_ratio(avg_len, pfa):
    """Threshold for one bin, as a multiple of the noise power.

    A noise-only bin of the averaged spectrum has the distribution of the mean
    of avg_len independent exponential variables, so
    P(bin > t * sigma^2) = Q(avg_len, avg_len * t), where Q is the regularised
    upper incomplete gamma function. This function returns the t for which
    this probability equals pfa (found by bisection).

    Example: avg_len = 1, pfa = 1e-3 gives -ln(1e-3) = 6.91.
    """
    if not 0.0 < pfa < 1.0:
        raise ValueError("pfa must be between 0 and 1")
    L = int(avg_len)

    def tail(x):  # Q(L, x) = exp(-x) * sum_{j<L} x^j / j!
        if x <= 0.0:
            return 1.0
        terms = [math.exp(j * math.log(x) - x - math.lgamma(j + 1.0)) for j in range(L)]
        return min(1.0, sum(terms))

    lo, hi = 0.0, 1.0
    while tail(L * hi) > pfa:
        hi *= 2.0
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if tail(L * mid) > pfa:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def total_power_threshold_ratio(fft_size, avg_len, pfa):
    """Threshold for the total power over all bins, as a multiple of the noise power.

    The mean over all bins of the averaged spectrum equals the Hann-weighted
    mean of |x|^2. For white noise its variance is
    sigma^4 * sum(w^4) / sum(w^2)^2 / avg_len, so it behaves like a Gamma
    variable with shape k = avg_len * sum(w^2)^2 / sum(w^4). The quantile of
    this Gamma variable is found with the Wilson-Hilferty approximation, which
    is accurate for large k (k is about 0.51 * fft_size * avg_len).
    """
    if not 0.0 < pfa < 1.0:
        raise ValueError("pfa must be between 0 and 1")
    w = numpy.hanning(int(fft_size) + 1)[:-1]
    k = int(avg_len) * numpy.sum(w ** 2) ** 2 / numpy.sum(w ** 4)
    z = NormalDist().inv_cdf(1.0 - pfa)
    return (1.0 - 1.0 / (9.0 * k) + z / (3.0 * math.sqrt(k))) ** 3


def occupied_bands(flags, fft_size):
    """Runs of occupied bins as (low, high) normalised frequencies.

    flags: 1-D array of fft_size booleans, fftshifted (index fft_size // 2 is
    0 Hz). Frequencies are in cycles per sample, from -0.5 to 0.5, and mark the
    centres of the first and last bin of each run.
    """
    flags = numpy.asarray(flags, dtype=bool)
    padded = numpy.concatenate(([False], flags, [False]))
    change = numpy.flatnonzero(padded[1:] != padded[:-1])
    starts, stops = change[0::2], change[1::2] - 1
    centre = fft_size // 2
    return [((a - centre) / fft_size, (b - centre) / fft_size) for a, b in zip(starts, stops)]


def detect(power, fft_size, avg_len, threshold_db=6.0, relative=True, pfa=0.0,
           noise_power=1.0, decision_mode=0, total_pfa=0.01):
    """Decision logic of the block, on already computed spectra.

    power: array (n_out, fft_size) from energy_spectrum().
    Returns (occupied, decision): boolean arrays of shape (n_out, fft_size) and (n_out,).
    decision_mode: 0 = any bin occupied, 1 = total power above its threshold,
    2 = either of the two.
    """
    power = numpy.atleast_2d(power)
    power_db = 10.0 * numpy.log10(power + _EPS)
    if pfa > 0.0:
        ratio = cfar_threshold_ratio(avg_len, pfa)
        if relative:
            base = numpy.array([noise_floor_estimate(row, avg_len) for row in power])[:, None]
        else:
            base = float(noise_power)
        occupied = power > base * ratio
    else:
        if relative:
            floor_db = numpy.array(
                [10.0 * numpy.log10(noise_floor_estimate(row, avg_len) + _EPS) for row in power]
            )
            threshold = floor_db[:, None] + threshold_db
        else:
            threshold = threshold_db
        occupied = power_db > threshold
    any_bin = occupied.any(axis=1)
    if decision_mode == 0:
        return occupied, any_bin
    total = power.mean(axis=1) > float(noise_power) * total_power_threshold_ratio(fft_size, avg_len, total_pfa)
    return occupied, (total if decision_mode == 1 else (total | any_bin))


class energy_detector_cf(gr.decim_block):
    """Spectrum-sensing energy detector.

    Input 0 : complex baseband samples.
    Output 0: power spectrum in dB, vector of fft_size floats, fftshifted.
    Output 1: occupancy flags, vector of fft_size bytes (1 = occupied).
    Output 2: band decision, one byte (see decision_mode).
    Message output "bands": for every output item, a list of (low, high)
    normalised frequencies of the occupied bands (cycles per sample).

    Args:
        fft_size: number of FFT bins (>= 2).
        avg_len: number of FFT frames averaged per output (>= 1).
        threshold_db: if relative is False, absolute threshold in dB (power
            1.0 = 0 dB). If relative is True, threshold above the estimated
            noise floor in dB.
        relative: True to measure the threshold from the estimated noise floor.
        pfa: if larger than 0, the per-bin threshold is calculated from theory
            for this false-alarm probability and threshold_db is ignored.
        noise_power: noise power (linear) that is assumed known. Used for the
            absolute CFAR threshold and for the total-power decision.
        decision_mode: 0 = any bin occupied, 1 = total power above threshold,
            2 = either.
        total_pfa: false-alarm probability of the total-power decision.
    """

    def __init__(self, fft_size=1024, avg_len=8, threshold_db=6.0, relative=True,
                 pfa=0.0, noise_power=1.0, decision_mode=0, total_pfa=0.01):
        fft_size = int(fft_size)
        avg_len = int(avg_len)
        if fft_size < 2:
            raise ValueError("fft_size must be at least 2")
        if avg_len < 1:
            raise ValueError("avg_len must be at least 1")
        if not 0.0 <= pfa < 1.0:
            raise ValueError("pfa must be 0 (off) or between 0 and 1")
        if noise_power <= 0.0:
            raise ValueError("noise_power must be positive")
        if int(decision_mode) not in (0, 1, 2):
            raise ValueError("decision_mode must be 0, 1 or 2")
        if not 0.0 < total_pfa < 1.0:
            raise ValueError("total_pfa must be between 0 and 1")
        gr.decim_block.__init__(
            self,
            name="energy_detector_cf",
            in_sig=[numpy.complex64],
            out_sig=[(numpy.float32, fft_size), (numpy.uint8, fft_size), numpy.uint8],
            decim=fft_size * avg_len,
        )
        self.fft_size = fft_size
        self.avg_len = avg_len
        self.threshold_db = float(threshold_db)
        self.relative = bool(relative)
        self.pfa = float(pfa)
        self.noise_power = float(noise_power)
        self.decision_mode = int(decision_mode)
        self.total_pfa = float(total_pfa)
        self.message_port_register_out(pmt.intern("bands"))

    def set_threshold_db(self, threshold_db):
        self.threshold_db = float(threshold_db)

    def set_relative(self, relative):
        self.relative = bool(relative)

    def set_pfa(self, pfa):
        self.pfa = float(pfa)

    def set_noise_power(self, noise_power):
        self.noise_power = float(noise_power)

    def set_decision_mode(self, decision_mode):
        self.decision_mode = int(decision_mode)

    def set_total_pfa(self, total_pfa):
        self.total_pfa = float(total_pfa)

    def work(self, input_items, output_items):
        n_out = len(output_items[0])
        power = energy_spectrum(input_items[0], self.fft_size, self.avg_len)
        occupied, decision = detect(
            power, self.fft_size, self.avg_len, self.threshold_db, self.relative,
            self.pfa, self.noise_power, self.decision_mode, self.total_pfa,
        )
        output_items[0][:] = 10.0 * numpy.log10(power + _EPS)
        output_items[1][:] = occupied
        output_items[2][:] = decision
        for row in occupied:
            bands = occupied_bands(row, self.fft_size)
            self.message_port_pub(
                pmt.intern("bands"),
                pmt.to_pmt([[float(a), float(b)] for a, b in bands]),
            )
        return n_out
