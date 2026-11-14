#!/usr/bin/env python3
# Copyright 2026 Group 05 ECE310.
# SPDX-License-Identifier: GPL-3.0-or-later
"""Terminal demo of the isac Energy Detector.

Makes two computer-generated test signals, runs each through the block inside a
real GNU Radio flowgraph, and prints what the block decided.

  Signal A: noise only            -> block should say "free"
  Signal B: noise + one tone      -> block should say "occupied" at the tone

Usage:  python3 demo_terminal.py
"""
import numpy as np
from gnuradio import blocks, gr
from gnuradio.isac import energy_detector_cf

FFT = 256            # number of frequency bins
AVG = 8              # FFT frames averaged per output
OUTPUTS = 5          # how many outputs we collect
FS = 1_000_000       # pretend sample rate: 1 MS/s, so bin width = FS/FFT
TONE_BIN = 32        # tone sits 32 bins above 0 Hz
THRESH_DB = 6.0      # threshold: 6 dB above the noise floor

rng = np.random.default_rng(1)
n = OUTPUTS * FFT * AVG


def noise():
    return ((rng.standard_normal(n) + 1j * rng.standard_normal(n)) / np.sqrt(2)).astype(np.complex64)


def tone(power):
    k = np.arange(n)
    return (np.sqrt(power) * np.exp(2j * np.pi * TONE_BIN * k / FFT)).astype(np.complex64)


def run(x):
    tb = gr.top_block()
    src = blocks.vector_source_c(x.tolist(), False)
    det = energy_detector_cf(FFT, AVG, THRESH_DB, True)
    s_psd, s_occ, s_dec = blocks.vector_sink_f(FFT), blocks.vector_sink_b(FFT), blocks.vector_sink_b(1)
    tb.connect(src, det)
    tb.connect((det, 0), s_psd)
    tb.connect((det, 1), s_occ)
    tb.connect((det, 2), s_dec)
    tb.run()
    m = len(s_dec.data())
    return (np.array(s_psd.data()).reshape(m, FFT),
            np.array(s_occ.data()).reshape(m, FFT),
            np.array(s_dec.data()))


def report(title, x, expect):
    psd, occ, dec = run(x)
    print("=" * 66)
    print(title)
    print("=" * 66)
    print(f"Block settings : FFT size {FFT}, {AVG} averages, threshold {THRESH_DB} dB above noise floor")
    print(f"Samples in     : {len(x)}  ->  outputs produced: {len(dec)}")
    print(f"Bin width      : {FS / FFT:.0f} Hz   (0 Hz is bin {FFT // 2})")
    print()
    print(" output | band decision | occupied bins (bin number = frequency)")
    for i in range(len(dec)):
        bins = np.nonzero(occ[i])[0]
        freqs = ", ".join(f"{(b - FFT // 2) * FS / FFT / 1e3:+.0f} kHz" for b in bins) or "none"
        verdict = "OCCUPIED" if dec[i] else "free"
        print(f"   {i + 1:4d} | {verdict:13s} | {list(bins) or '-'}  ({freqs})")
    peak = int(np.argmax(psd.mean(axis=0)))
    others = np.delete(psd.mean(axis=0), [peak - 1, peak, peak + 1])
    print()
    print(f"Strongest bin  : {peak}  ({(peak - FFT // 2) * FS / FFT / 1e3:+.0f} kHz), {psd.mean(axis=0)[peak]:.2f} dB")
    print(f"Average of other bins (noise level): {10 * np.log10((10 ** (others / 10)).mean()):.2f} dB   (true noise = 0 dB)")
    print(f"Expected       : {expect}")
    print()


report("TEST A: noise only", noise(), "all outputs 'free', no occupied bins")
report(f"TEST B: noise + tone at bin {FFT // 2 + TONE_BIN} (+{TONE_BIN * FS / FFT / 1e3:.0f} kHz)",
       noise() + tone(1.0),
       f"all outputs 'OCCUPIED', occupied bins around {FFT // 2 + TONE_BIN} "
       "(the window spreads one tone over up to 3 neighbouring bins)")
