#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright 2026 Group 05 ECE310.
#
# SPDX-License-Identifier: GPL-3.0-or-later
#
"""Simulated channel with moving targets (for tests and demos without hardware)."""

import math

import numpy
from gnuradio import gr

from .ofdm_isac import C0


def fractional_delay_taps(delay, n_taps):
    """Windowed-sinc FIR that delays a signal by (n_taps - 1) / 2 + delay samples, 0 <= delay < 1."""
    k = numpy.arange(n_taps) - (n_taps - 1) / 2.0 - delay
    return numpy.sinc(k) * numpy.hanning(n_taps + 2)[1:-1]


class moving_target_channel_cc(gr.sync_block):
    """Echo channel: delayed, Doppler-shifted, attenuated copies of the input, plus noise.

    For every target i with range R_i (m), speed v_i (m/s, positive = moving towards
    the radar) and echo gain g_i (dB, relative to the input power):
        y[n] = sum_i 10^(g_i/20) * x[n - D_i] * exp(j 2 pi (f_i n / fs - fc * tau_i)) + noise
    with the echo delay tau_i = 2 R_i / c and the Doppler shift f_i = 2 v_i fc / c.
    With one_way=True the signal travels only once (a communication link):
    tau_i = R_i / c and f_i = v_i fc / c.

    The delay is done with a windowed-sinc FIR of n_taps taps, so it need not be a
    whole number of samples. This filter adds a fixed latency of (n_taps - 1) / 2
    samples to every path (property latency). A receiver has to subtract it.

    Args:
        targets: list of (range_m, speed_mps, gain_db) tuples.
        sample_rate: samples per second.
        carrier_freq: carrier frequency in Hz.
        noise_power: power of the complex Gaussian noise (linear, input power is 1.0).
        one_way: False for a radar echo, True for a one-way communication path.
        n_taps: odd number of FIR taps of the fractional delay (default 33).
        max_range: largest range (m) that set_targets() may use later.
        seed: seed of the noise generator.
    """

    def __init__(self, targets=((300.0, 20.0, -20.0),), sample_rate=10e6, carrier_freq=5.8e9,
                 noise_power=1e-3, one_way=False, n_taps=33, max_range=1000.0, seed=7):
        if sample_rate <= 0 or carrier_freq <= 0:
            raise ValueError("sample_rate and carrier_freq must be positive")
        if noise_power < 0:
            raise ValueError("noise_power must not be negative")
        if int(n_taps) < 3 or int(n_taps) % 2 == 0:
            raise ValueError("n_taps must be odd and at least 3")
        self.sample_rate = float(sample_rate)
        self.carrier_freq = float(carrier_freq)
        self.noise_power = float(noise_power)
        self.one_way = bool(one_way)
        self.n_taps = int(n_taps)
        legs = 1.0 if self.one_way else 2.0
        self._max_int_delay = int(math.ceil(legs * max(float(max_range), 0.0) / C0 * self.sample_rate)) + 1
        gr.sync_block.__init__(
            self, name="moving_target_channel_cc", in_sig=[numpy.complex64], out_sig=[numpy.complex64]
        )
        self.set_history(self._max_int_delay + self.n_taps)
        self._rng = numpy.random.default_rng(int(seed))
        self.set_targets(targets)

    @property
    def latency(self):
        """Fixed delay of the fractional-delay filter, in samples."""
        return (self.n_taps - 1) // 2

    def set_targets(self, targets):
        legs = 1.0 if self.one_way else 2.0
        paths = []
        for r, v, g in targets:
            if r < 0:
                raise ValueError("range must not be negative")
            delay = legs * float(r) / C0 * self.sample_rate  # samples
            n_int = int(math.floor(delay))
            if n_int > self._max_int_delay:
                raise ValueError("range %.0f m is larger than max_range" % r)
            taps = fractional_delay_taps(delay - n_int, self.n_taps)
            f_dop = legs * float(v) * self.carrier_freq / C0 / self.sample_rate  # cycles per sample
            phase0 = -2.0 * math.pi * ((self.carrier_freq * legs * float(r) / C0) % 1.0)
            paths.append((n_int, taps, f_dop, phase0, 10.0 ** (float(g) / 20.0)))
        self._paths = paths

    def set_noise_power(self, noise_power):
        if noise_power < 0:
            raise ValueError("noise_power must not be negative")
        self.noise_power = float(noise_power)

    def work(self, input_items, output_items):
        x = input_items[0]
        out = output_items[0]
        n = len(out)
        hist = self.history() - 1
        n0 = self.nitems_written(0)
        idx = numpy.arange(n0, n0 + n, dtype=numpy.float64)
        y = numpy.zeros(n, dtype=numpy.complex128)
        for n_int, taps, f_dop, phase0, amp in self._paths:
            # output[i] needs x[i - n_int - j] for j = 0 .. n_taps-1; in the buffer x index shifts by hist
            seg = x[hist - n_int - (self.n_taps - 1): hist - n_int + n].astype(numpy.complex128)
            delayed = numpy.convolve(seg, taps, mode="valid")[:n]
            rot = numpy.exp(1j * (2.0 * math.pi * ((f_dop * idx) % 1.0) + phase0))
            y += amp * delayed * rot
        if self.noise_power > 0.0:
            s = math.sqrt(self.noise_power / 2.0)
            y += s * (self._rng.standard_normal(n) + 1j * self._rng.standard_normal(n))
        out[:] = y
        return n
