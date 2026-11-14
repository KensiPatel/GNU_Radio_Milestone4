#!/usr/bin/env python3
# Copyright 2026 Group 05 ECE310.
# SPDX-License-Identifier: GPL-3.0-or-later
"""Terminal demo of the features added to the Energy Detector in M4.

Part 1: CFAR threshold. We give a false-alarm probability and the block
        calculates the threshold from theory.
Part 2: total-power decision. A signal that fills the whole band, only 15 %
        stronger than the noise, is invisible bin by bin but easy to see in the total.
Part 3: message output "bands". Two tones give two occupied bands.

Usage:  python3 demo_new_features.py
"""
import numpy as np
import pmt
from gnuradio import blocks, gr
from gnuradio.isac import energy_detector_cf
from gnuradio.isac.energy_detector_cf import cfar_threshold_ratio

FFT, AVG, OUTPUTS = 256, 8, 400
FS = 1_000_000                # pretend sample rate: 1 MS/s
rng = np.random.default_rng(3)


def noise(n_out, power=1.0):
    n = n_out * FFT * AVG
    s = np.sqrt(power / 2.0)
    return (s * (rng.standard_normal(n) + 1j * rng.standard_normal(n))).astype(np.complex64)


def tone(n, bin_index, amplitude):
    return (amplitude * np.exp(2j * np.pi * bin_index * np.arange(n) / FFT)).astype(np.complex64)


def run(x, with_messages=False, **kw):
    tb = gr.top_block()
    det = energy_detector_cf(FFT, AVG, **kw)
    dec = blocks.vector_sink_b(1)
    tb.connect(blocks.vector_source_c(x.tolist(), False), det)
    tb.connect((det, 0), blocks.null_sink(4 * FFT))
    tb.connect((det, 1), blocks.null_sink(FFT))
    tb.connect((det, 2), dec)
    dbg = blocks.message_debug()
    if with_messages:
        tb.msg_connect(det, "bands", dbg, "store")
    tb.run()
    return np.array(dec.data()), dbg


print("Part 1: CFAR threshold (noise power 1, false-alarm probability 1e-3 per bin)")
for L in (1, 8, 16):
    t = cfar_threshold_ratio(L, 1e-3)
    print(f"  {L:2d} averages: threshold = {t:5.2f} x noise power = {10 * np.log10(t):5.2f} dB")

print("\nPart 2: a wideband signal, 15 % more power than the noise (400 outputs each)")
common = dict(threshold_db=6.0, relative=False, pfa=1e-4, noise_power=1.0, total_pfa=0.01)
for label, x in (("noise only", noise(OUTPUTS)), ("noise + wideband signal", noise(OUTPUTS, 1.15))):
    d_bins, _ = run(x, decision_mode=0, **common)
    d_total, _ = run(x, decision_mode=1, **common)
    print(f"  {label:24s} says 'busy': per-bin test {100 * d_bins.mean():5.1f} %,"
          f" total power {100 * d_total.mean():5.1f} % of the outputs")

print("\nPart 3: message output 'bands' (two tones, at +40 and -60 bins)")
n_out = 3
x = noise(n_out) + tone(n_out * FFT * AVG, 40, 0.6) + tone(n_out * FFT * AVG, -60, 0.6)
_, dbg = run(x, with_messages=True, relative=True, pfa=1e-5)
msg = dbg.get_message(0)
for i in range(pmt.length(msg)):
    lo, hi = (pmt.to_double(pmt.vector_ref(pmt.vector_ref(msg, i), j)) for j in range(2))
    print(f"  occupied band {i + 1}: {lo * FS / 1e3:7.1f} kHz to {hi * FS / 1e3:7.1f} kHz")
