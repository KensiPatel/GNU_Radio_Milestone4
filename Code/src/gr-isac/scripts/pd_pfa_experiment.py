#!/usr/bin/env python3
# Copyright 2026 Group 05 ECE310.
# SPDX-License-Identifier: GPL-3.0-or-later
"""Detection experiment for the isac Energy Detector.

Measures, by Monte Carlo simulation on synthetic data:
  * the false-alarm probability Pfa on noise-only bins, and
  * the detection probability Pd of one complex tone in white noise,
for several averaging lengths, and compares both with the closed-form theory.

Usage:  python3 pd_pfa_experiment.py OUTPUT_DIR
Needs GNU Radio (for the final block-in-flowgraph check), numpy, scipy, matplotlib.
"""
import csv
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import chi2, ncx2

from gnuradio import blocks, gr
from gnuradio.isac import energy_detector_cf
from gnuradio.isac.energy_detector_cf import energy_spectrum

N = 256              # FFT size
TONE_BIN = 32        # tone sits at +32 bins from DC
PFA_TARGET = 1e-3    # per-bin false-alarm probability used to set the threshold
TRIALS = 4000        # independent outputs per point
LS = [1, 4, 16]      # averaging lengths
SNR_DB = np.arange(-34, -9, 1.0)   # per-sample tone-to-noise power ratio
WIN = np.hanning(N + 1)[:-1]
GAIN = WIN.sum() ** 2 / (WIN ** 2).sum()   # processing gain of the windowed FFT (170.7 = 22.3 dB)


def cnoise(rng, shape):
    return (rng.standard_normal(shape) + 1j * rng.standard_normal(shape)) / np.sqrt(2)


def threshold(L):
    """Linear power threshold giving PFA_TARGET per bin when noise power is 1."""
    return chi2.isf(PFA_TARGET, 2 * L) / (2 * L)


def theory_pd(snr_lin, L):
    lam = 2 * L * snr_lin * GAIN
    return ncx2.sf(2 * L * threshold(L), 2 * L, lam)


def simulate(L, snr_db, rng):
    n = TRIALS * L * N
    k = np.arange(n)
    phase = rng.uniform(0, 2 * np.pi)
    x = cnoise(rng, n) + np.sqrt(10 ** (snr_db / 10)) * np.exp(1j * (2 * np.pi * TONE_BIN * k / N + phase))
    p = energy_spectrum(x, N, L)
    return p, x


def main(outdir):
    os.makedirs(outdir, exist_ok=True)
    rng = np.random.default_rng(2026)
    rows = []

    # 1) Pfa from noise-only data, all bins.
    pfa_rows = []
    for L in LS:
        n = TRIALS * L * N
        p = energy_spectrum(cnoise(rng, n), N, L)
        measured = float((p > threshold(L)).mean())
        pfa_rows.append((L, PFA_TARGET, measured))

    # 2) Pd versus SNR at the tone bin.
    pd_curves = {}
    for L in LS:
        meas, theo = [], []
        for snr in SNR_DB:
            p, _ = simulate(L, snr, rng)
            tone_idx = N // 2 + TONE_BIN
            meas.append(float((p[:, tone_idx] > threshold(L)).mean()))
            theo.append(float(theory_pd(10 ** (snr / 10), L)))
            rows.append((L, snr, meas[-1], theo[-1]))
        pd_curves[L] = (meas, theo)

    # 3) The block itself in a real flowgraph gives the same decisions as the core function.
    L, snr = 8, -24.0
    p_core, x = simulate(L, snr, rng)
    thr_db = 10 * np.log10(threshold(L))
    tb = gr.top_block()
    src = blocks.vector_source_c(x.astype(np.complex64).tolist(), False)
    det = energy_detector_cf(N, L, thr_db, False)
    sink = blocks.vector_sink_b(N)
    tb.connect(src, det)
    tb.connect((det, 1), sink)
    tb.connect((det, 0), blocks.null_sink(gr.sizeof_float * N))  # every output port must be connected
    tb.connect((det, 2), blocks.null_sink(gr.sizeof_char))
    tb.run()
    occ = np.array(sink.data()).reshape(-1, N)
    tone_idx = N // 2 + TONE_BIN
    block_pd = float(occ[:, tone_idx].mean())
    core_pd = float((p_core[:, tone_idx] > threshold(L)).mean())
    block_pfa = float(np.delete(occ, [tone_idx - 1, tone_idx, tone_idx + 1], axis=1).mean())

    with open(os.path.join(outdir, "pd_vs_snr.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["avg_len", "snr_db_per_sample", "pd_measured", "pd_theory"])
        w.writerows(rows)
    with open(os.path.join(outdir, "pfa_check.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["avg_len", "pfa_target", "pfa_measured"])
        w.writerows(pfa_rows)

    plt.figure(figsize=(6.4, 4.2))
    for L, c in zip(LS, ["tab:blue", "tab:orange", "tab:green"]):
        meas, theo = pd_curves[L]
        plt.plot(SNR_DB, theo, "-", color=c, label=f"theory, {L} avg")
        plt.plot(SNR_DB, meas, "o", color=c, ms=4, label=f"simulation, {L} avg")
    plt.axhline(0.9, color="gray", ls=":", lw=1)
    plt.xlabel("Tone-to-noise power per sample (dB)")
    plt.ylabel("Detection probability Pd")
    plt.title(f"Energy detector, FFT {N}, Pfa = {PFA_TARGET:g} per bin")
    plt.grid(alpha=0.3)
    plt.legend(fontsize=7, ncol=2)
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, "pd_vs_snr.png"), dpi=200)

    # Required SNR for Pd = 0.9, measured (interpolated) and theory
    summary = []
    for L in LS:
        meas, theo = pd_curves[L]

        def cross(v):
            v = np.array(v)
            i = np.argmax(v >= 0.9)
            if v[i] < 0.9:
                return float("nan")
            if i == 0:
                return float(SNR_DB[0])
            return float(np.interp(0.9, [v[i - 1], v[i]], [SNR_DB[i - 1], SNR_DB[i]]))
        summary.append((L, cross(meas), cross(theo)))

    lines = ["Energy detector experiment", f"FFT size {N}, trials per point {TRIALS}, seed 2026",
             f"Window gain (Hann) {10*np.log10(GAIN):.2f} dB", "", "False-alarm check (noise only, all bins):"]
    for L, t, m in pfa_rows:
        lines.append(f"  avg_len {L:2d}: target {t:.1e}, measured {m:.2e}")
    lines += ["", "SNR per sample needed for Pd = 0.9:"]
    for L, m, t in summary:
        lines.append(f"  avg_len {L:2d}: simulation {m:6.2f} dB, theory {t:6.2f} dB")
    lines += ["", "Block-in-flowgraph check (avg_len 8, SNR -24 dB, absolute threshold):",
              f"  Pd from block {block_pd:.4f}, Pd from core function {core_pd:.4f}, "
              f"Pfa on other bins from block {block_pfa:.2e}, theory Pd {float(theory_pd(10**(-2.4), 8)):.4f}"]
    text = "\n".join(lines)
    print(text)
    with open(os.path.join(outdir, "summary.txt"), "w") as f:
        f.write(text + "\n")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "results")
