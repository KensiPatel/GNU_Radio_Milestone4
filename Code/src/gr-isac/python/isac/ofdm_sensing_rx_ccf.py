#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright 2026 Group 05 ECE310.
#
# SPDX-License-Identifier: GPL-3.0-or-later
#
"""OFDM sensing receiver: range-Doppler map and target list."""

import numpy
import pmt
from gnuradio import gr

from .ofdm_isac import (OfdmConfig, axes, ca_cfar_2d, demodulate_symbols, extract_targets,
                        range_doppler_map)


class ofdm_sensing_rx_ccf(gr.decim_block):
    """Radar receiver of the ISAC chain (monostatic: it knows what was transmitted).

    Input 0: received signal (the echo), complex.
    Input 1: the transmitted signal, complex, time-aligned with the start of a frame.
    Output 0: range-Doppler map in dB, one vector of n_doppler * n_range floats per frame
        (row = Doppler, fftshifted, row n_doppler // 2 is zero speed; column = range).
    Message output "targets": for every frame a list of [range_m, speed_mps, snr_db]
        for every detected target (positive speed = moving towards the receiver).

    Args:
        fft_size, cp_len, n_symbols, n_active, pilot_spacing, sample_rate, carrier_freq:
            must be the same as in the transmitter, see OfdmConfig.
        pfa: false-alarm probability per map cell of the CA-CFAR detector.
        use_window: Hann window over subcarriers and symbols (lower side lobes).
        range_pad, doppler_pad: zero-padding factors of the two FFTs (finer map grid).
        guard_cells, train_cells: CFAR half-widths in cells, given as (Doppler, range).
        latency_samples: fixed delay of the channel in samples (the channel block
            reports it as its property "latency"). The map starts after it, so that
            column 0 is range 0. The largest range is then (cp_len - latency_samples)
            samples of delay, i.e. c * (cp_len - latency_samples) / (2 * sample_rate).
    """

    def __init__(self, fft_size=256, cp_len=64, n_symbols=128, n_active=192, pilot_spacing=4,
                 sample_rate=10e6, carrier_freq=5.8e9, pfa=1e-8, use_window=True,
                 range_pad=2, doppler_pad=2, guard_cells=(4, 4), train_cells=(8, 8), latency_samples=0):
        self.cfg = OfdmConfig(int(fft_size), int(cp_len), int(n_symbols), int(n_active), int(pilot_spacing),
                              float(sample_rate), float(carrier_freq))
        if not 0.0 < pfa < 1.0:
            raise ValueError("pfa must be between 0 and 1")
        if int(range_pad) < 1 or int(doppler_pad) < 1:
            raise ValueError("range_pad and doppler_pad must be at least 1")
        if int(latency_samples) < 0 or int(latency_samples) >= self.cfg.cp_len:
            raise ValueError("latency_samples must be between 0 and cp_len - 1")
        self.pfa = float(pfa)
        self.use_window = bool(use_window)
        self.range_pad = int(range_pad)
        self.doppler_pad = int(doppler_pad)
        self.guard_cells = tuple(int(g) for g in guard_cells)
        self.train_cells = tuple(int(t) for t in train_cells)
        if 2 * (self.guard_cells[0] + self.train_cells[0]) + 1 > self.cfg.n_symbols * self.doppler_pad:
            raise ValueError("training window is larger than the Doppler axis")
        self.latency_samples = int(latency_samples)
        self.n_doppler = self.cfg.n_symbols * self.doppler_pad
        self.n_range = (self.cfg.cp_len - self.latency_samples) * self.range_pad
        self.last_targets = []
        gr.decim_block.__init__(
            self,
            name="ofdm_sensing_rx_ccf",
            in_sig=[numpy.complex64, numpy.complex64],
            out_sig=[(numpy.float32, self.n_doppler * self.n_range)],
            decim=self.cfg.frame_len,
        )
        self.message_port_register_out(pmt.intern("targets"))

    def set_pfa(self, pfa):
        if not 0.0 < pfa < 1.0:
            raise ValueError("pfa must be between 0 and 1")
        self.pfa = float(pfa)

    def range_axis(self):
        return axes(self.cfg, self.range_pad, self.doppler_pad, self.latency_samples)[0]

    def speed_axis(self):
        return axes(self.cfg, self.range_pad, self.doppler_pad, self.latency_samples)[1]

    def process_frame(self, rx, tx):
        """Range-Doppler map (complex, shape n_doppler x n_range) and target list of one frame."""
        Y = demodulate_symbols(rx, self.cfg)
        X = demodulate_symbols(tx, self.cfg)
        cmap = range_doppler_map(Y, X, self.cfg, self.use_window, self.range_pad, self.doppler_pad, self.latency_samples)
        power = numpy.abs(cmap) ** 2
        det, noise = ca_cfar_2d(power, self.pfa, self.guard_cells, self.train_cells)
        targets = extract_targets(cmap, det, noise, self.cfg, self.range_pad, self.doppler_pad)
        return cmap, targets

    def work(self, input_items, output_items):
        n_frames = len(output_items[0])
        fl = self.cfg.frame_len
        for f in range(n_frames):
            rx = input_items[0][f * fl:(f + 1) * fl]
            tx = input_items[1][f * fl:(f + 1) * fl]
            cmap, targets = self.process_frame(rx, tx)
            output_items[0][f, :] = (10.0 * numpy.log10(numpy.abs(cmap) ** 2 + 1e-30)).reshape(-1)
            self.last_targets = targets
            self.message_port_pub(pmt.intern("targets"), pmt.to_pmt([[r, v, s] for r, v, s in targets]))
        return n_frames
