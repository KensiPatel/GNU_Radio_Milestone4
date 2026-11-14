#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright 2026 Group 05 ECE310.
#
# SPDX-License-Identifier: GPL-3.0-or-later
#
"""OFDM ISAC transmitter (source block, no input)."""

import numpy
from gnuradio import gr

from .ofdm_isac import OfdmConfig, frame_symbols, modulate


class ofdm_isac_tx_c(gr.sync_block):
    """Transmitter of the ISAC chain: a stream of OFDM frames carrying random QPSK data.

    Output 0: complex baseband samples, average power 1.0. Frame k starts at
    sample k * frame_len. Frame k holds the data bits of frame_symbols(cfg, seed, k),
    so any receiver that knows the seed can check the bits.

    Args:
        fft_size, cp_len, n_symbols, n_active, pilot_spacing: see OfdmConfig.
        sample_rate: samples per second (needed only for the physical units).
        carrier_freq: carrier frequency in Hz (needed only for the physical units).
        seed: seed of the data generator.
    """

    def __init__(self, fft_size=256, cp_len=64, n_symbols=128, n_active=192, pilot_spacing=4,
                 sample_rate=10e6, carrier_freq=5.8e9, seed=1):
        self.cfg = OfdmConfig(int(fft_size), int(cp_len), int(n_symbols), int(n_active), int(pilot_spacing),
                              float(sample_rate), float(carrier_freq))
        gr.sync_block.__init__(self, name="ofdm_isac_tx_c", in_sig=None, out_sig=[numpy.complex64])
        self.seed = int(seed)
        self._frame = 0
        self._buf = numpy.zeros(0, dtype=numpy.complex64)

    def frame(self, index):
        """Time-domain samples of frame number index (also usable without a flowgraph)."""
        X, _ = frame_symbols(self.cfg, self.seed, index)
        return modulate(X, self.cfg)

    def work(self, input_items, output_items):
        out = output_items[0]
        n = len(out)
        while self._buf.size < n:
            self._buf = numpy.concatenate((self._buf, self.frame(self._frame)))
            self._frame += 1
        out[:] = self._buf[:n]
        self._buf = self._buf[n:]
        return n
