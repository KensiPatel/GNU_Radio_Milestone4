#!/usr/bin/env python3
# Copyright 2026 Group 05 ECE310.
# SPDX-License-Identifier: GPL-3.0-or-later
"""M4 experiments for the isac Energy Detector.

Three experiments, all on synthetic data with fixed random seeds:

  bench   Throughput of the block inside a GNU Radio flowgraph (MS/s).
  noise   What happens to the false-alarm rate when the real noise power is not
          the power the threshold was designed for (noise uncertainty), for the
          absolute-threshold mode and the relative (median noise floor) mode.
  base    Comparison with two baselines on the same sample budget:
          (B) one long FFT without averaging, and
          (C) a plain time-domain energy detector, which is what the standard
          GNU Radio blocks (Complex to Mag^2, Moving Average, threshold) compute.

Usage:
    python3 m4_experiments.py OUTPUT_DIR [--only bench|noise|base|new]

Needs GNU Radio with the gr-isac module installed, numpy, scipy, matplotlib.
Running everything takes about 3 to 5 minutes on a laptop.
"""
import argparse
import csv
import os
import platform
import statistics
import time

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import chi2, ncx2  # noqa: F401  (ncx2 kept for readers who extend the script)

from gnuradio import blocks, gr
from gnuradio.isac import energy_detector_cf
from gnuradio.isac.energy_detector_cf import detect, energy_spectrum, noise_floor_estimate

# Colours: first three slots of a colour-blind-checked categorical palette.
BLUE, ORANGE, AQUA, INK, GREY = "#2a78d6", "#eb6834", "#1baf7a", "#0b0b0b", "#52514e"


def style_axes(ax):
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.grid(alpha=0.25, lw=0.6)
    ax.tick_params(labelsize=8)


def cnoise(rng, shape):
    """Complex white Gaussian noise of unit power."""
    return (rng.standard_normal(shape) + 1j * rng.standard_normal(shape)) / np.sqrt(2)


# --------------------------------------------------------------------------- bench
def run_flowgraph(n_samples, fft_size, avg_len, relative):
    tb = gr.top_block()
    src = blocks.null_source(gr.sizeof_gr_complex)
    head = blocks.head(gr.sizeof_gr_complex, n_samples)
    det = energy_detector_cf(fft_size, avg_len, 6.0, relative)
    tb.connect(src, head, det)
    tb.connect((det, 0), blocks.null_sink(gr.sizeof_float * fft_size))
    tb.connect((det, 1), blocks.null_sink(gr.sizeof_char * fft_size))
    tb.connect((det, 2), blocks.null_sink(gr.sizeof_char))
    t0 = time.perf_counter()
    tb.run()
    return time.perf_counter() - t0


def run_empty_flowgraph(n_samples):
    """Scheduler-only floor: the same source and head, no detector."""
    tb = gr.top_block()
    tb.connect(blocks.null_source(gr.sizeof_gr_complex),
               blocks.head(gr.sizeof_gr_complex, n_samples),
               blocks.null_sink(gr.sizeof_gr_complex))
    t0 = time.perf_counter()
    tb.run()
    return time.perf_counter() - t0


def exp_bench(outdir, repeats=3, total=2 ** 24):
    configs = [(256, 1, True), (1024, 8, False), (1024, 8, True), (4096, 8, True)]
    rows = []
    base_t = statistics.median(run_empty_flowgraph(total) for _ in range(repeats))
    rows.append(("no detector (source -> head -> sink)", "-", "-", "-", total / base_t / 1e6))
    for fft_size, avg_len, relative in configs:
        per_out = fft_size * avg_len
        n = (total // per_out) * per_out
        t = statistics.median(run_flowgraph(n, fft_size, avg_len, relative) for _ in range(repeats))
        rows.append(("energy_detector_cf", fft_size, avg_len,
                     "relative" if relative else "absolute", n / t / 1e6))
    with open(os.path.join(outdir, "bench_throughput.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["what", "fft_size", "avg_len", "threshold_mode", "throughput_Msps"])
        for r in rows:
            w.writerow(r[:4] + (f"{r[4]:.2f}",))
    lines = [f"Throughput benchmark: {total} samples per run, median of {repeats} runs.",
             "Input is complex64 from a Null Source (the speed of the FFT does not depend on the data).", ""]
    for r in rows:
        lines.append(f"  {r[0]:38s} fft {str(r[1]):>5s} avg {str(r[2]):>2s} {str(r[3]):9s} {r[4]:8.2f} MS/s")
    return lines, rows


# --------------------------------------------------------------------------- noise
def exp_noise(outdir, rng, n_fft=256, L=8, trials=4000, pfa_target=1e-3):
    lam = chi2.isf(pfa_target, 2 * L) / (2 * L)          # linear threshold for noise power 1
    rel_db = 10 * np.log10(lam)                           # the same number used as "dB above floor"
    deltas = [-2.0, -1.0, -0.5, 0.0, 0.5, 1.0, 2.0, 3.0]  # real noise power minus design noise power (dB)
    tone_snr_db = -18.0                                   # per-sample tone-to-nominal-noise ratio
    tone_bin = 32
    k = np.arange(trials * L * n_fft)
    tone = np.sqrt(10 ** (tone_snr_db / 10)) * np.exp(2j * np.pi * tone_bin * k / n_fft)
    base_noise = cnoise(rng, k.size)
    rows = []
    for d in deltas:
        sigma2 = 10 ** (d / 10)
        p_noise = energy_spectrum(np.sqrt(sigma2) * base_noise, n_fft, L)
        p_sig = energy_spectrum(np.sqrt(sigma2) * base_noise + tone, n_fft, L)
        floor_noise = np.array([noise_floor_estimate(r, L) for r in p_noise])
        floor_sig = np.array([noise_floor_estimate(r, L) for r in p_sig])
        thr_rel_noise = floor_noise[:, None] * lam
        thr_rel_sig = floor_sig[:, None] * lam
        ti = n_fft // 2 + tone_bin
        pfa_abs = float((p_noise > lam).mean())
        pfa_rel = float((p_noise > thr_rel_noise).mean())
        pd_abs = float((p_sig[:, ti] > lam).mean())
        pd_rel = float((p_sig[:, ti] > thr_rel_sig[:, 0]).mean())
        pfa_abs_theory = float(chi2.sf(2 * L * lam / sigma2, 2 * L))
        rows.append((d, pfa_abs, pfa_abs_theory, pfa_rel, pd_abs, pd_rel))

    # Check: the real block in a flowgraph flags exactly the same bins as the core functions.
    n_out = 500
    d = 1.0
    x = (np.sqrt(10 ** (d / 10)) * (cnoise(rng, n_out * L * n_fft))).astype(np.complex64)
    match = {}
    for relative, thr_db in ((False, rel_db), (True, rel_db)):
        tb = gr.top_block()
        det = energy_detector_cf(n_fft, L, thr_db, relative)
        s_occ = blocks.vector_sink_b(n_fft)
        tb.connect(blocks.vector_source_c(x.tolist(), False), det)
        tb.connect((det, 1), s_occ)
        tb.connect((det, 0), blocks.null_sink(gr.sizeof_float * n_fft))
        tb.connect((det, 2), blocks.null_sink(gr.sizeof_char))
        tb.run()
        occ = np.array(s_occ.data()).reshape(-1, n_fft).astype(bool)
        p = energy_spectrum(x, n_fft, L)
        if relative:
            fl = np.array([noise_floor_estimate(r, L) for r in p])
            ref = 10 * np.log10(p + 1e-20) > 10 * np.log10(fl + 1e-20)[:, None] + thr_db
        else:
            ref = 10 * np.log10(p + 1e-20) > thr_db
        match["relative" if relative else "absolute"] = (int((occ != ref).sum()), occ.size)

    with open(os.path.join(outdir, "noise_uncertainty.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["noise_power_error_db", "pfa_absolute_measured", "pfa_absolute_theory",
                    "pfa_relative_measured", "pd_absolute", "pd_relative"])
        for r in rows:
            w.writerow([f"{r[0]:+.1f}"] + [f"{v:.3e}" for v in r[1:]])

    fig, ax = plt.subplots(figsize=(6.4, 3.9))
    xs = [r[0] for r in rows]
    ax.semilogy(xs, [r[2] for r in rows], "--", color=GREY, lw=1.4, label="absolute: theory")
    ax.semilogy(xs, [max(r[1], 1e-5) for r in rows], "o", color=ORANGE, ms=6, label="absolute threshold: measured")
    ax.semilogy(xs, [max(r[3], 1e-5) for r in rows], "s-", color=BLUE, ms=6, lw=1.6, label="relative threshold: measured")
    ax.axhline(pfa_target, color=INK, lw=0.8, ls=":")
    ax.text(0.6, pfa_target * 0.55, "design target 1e-3", fontsize=8, color=GREY, va="top")
    ax.set_xlabel("Real noise power minus design noise power (dB)", fontsize=9)
    ax.set_ylabel("False-alarm probability per bin", fontsize=9)
    ax.set_title(f"Noise uncertainty, FFT {n_fft}, {L} averages", fontsize=10)
    ax.legend(fontsize=8, frameon=False, loc="upper left")
    style_axes(ax)
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "noise_uncertainty.png"), dpi=200)
    plt.close(fig)

    lines = [f"Noise-uncertainty experiment: FFT {n_fft}, {L} averages, {trials} outputs per point, "
             f"design Pfa {pfa_target:g}, tone at {tone_snr_db:.0f} dB per sample (designed noise power 1).",
             "  err dB | Pfa abs (meas) | Pfa abs (theory) | Pfa rel (meas) | Pd abs | Pd rel"]
    for r in rows:
        lines.append(f"  {r[0]:+5.1f}  |   {r[1]:.2e}     |   {r[2]:.2e}      |   {r[3]:.2e}     | {r[4]:.3f}  | {r[5]:.3f}")
    lines.append(f"Block in flowgraph vs core functions, differing flags / total flags: "
                 f"absolute {match['absolute'][0]} / {match['absolute'][1]}, "
                 f"relative {match['relative'][0]} / {match['relative'][1]}")
    return lines


# --------------------------------------------------------------------------- baselines
def required_snr(snr_db, pd, level=0.9):
    pd = np.asarray(pd)
    if pd.max() < level:
        return float("nan")
    i = int(np.argmax(pd >= level))
    if i == 0:
        return float(snr_db[0])
    return float(np.interp(level, [pd[i - 1], pd[i]], [snr_db[i - 1], snr_db[i]]))


def exp_base(outdir, rng, budget=2048, trials=4000, band_pfa=0.01):
    """Three detectors, same number of samples per decision (budget), same band-level Pfa.

    A: FFT 256 x 8 averages, per-bin threshold from the Sidak correction.
    B: FFT 2048 x 1 (no averaging), per-bin threshold from the Sidak correction.
    C: time-domain energy detector over all `budget` samples (one number for the whole band).
    The decision of A and B is "at least one bin above threshold" (the block's decision output).
    The noise power is assumed known (1.0). The tone sits half a bin off-centre, the worst case.
    """
    nA, LA = 256, budget // 256
    nB, LB = budget, 1

    def bin_pfa(nb):
        return 1 - (1 - band_pfa) ** (1.0 / nb)

    thrA = chi2.isf(bin_pfa(nA), 2 * LA) / (2 * LA)
    thrB = chi2.isf(bin_pfa(nB), 2 * LB) / (2 * LB)
    thrC = chi2.isf(band_pfa, 2 * budget) / (2 * budget)

    def decisions(x):
        pA = energy_spectrum(x, nA, LA)
        pB = energy_spectrum(x, nB, LB)
        a = (pA > thrA).any(axis=1)
        b = (pB > thrB).any(axis=1)
        c = (np.abs(x.reshape(-1, budget)) ** 2).mean(axis=1) > thrC
        return a.mean(), b.mean(), c.mean()

    # False-alarm calibration on noise only, 5 x trials outputs.
    pf = np.zeros(3)
    reps = 5
    for _ in range(reps):
        pf += np.array(decisions(cnoise(rng, trials * budget)))
    pf /= reps

    k = np.arange(trials * budget)
    # Tone offsets: half a bin for each FFT size (worst case for that detector).
    fA = 32.5 / nA
    fB = 260.5 / nB
    out = {"narrowband": [], "wideband": []}
    snr_nb = np.arange(-30, -5, 1.0)
    snr_wb = np.arange(-16, 5, 1.0)
    noise = cnoise(rng, trials * budget)
    for s_db in snr_nb:
        a_ = np.sqrt(10 ** (s_db / 10))
        ph = rng.uniform(0, 2 * np.pi)
        pA = energy_spectrum(noise + a_ * np.exp(1j * (2 * np.pi * fA * k + ph)), nA, LA)
        pB = energy_spectrum(noise + a_ * np.exp(1j * (2 * np.pi * fB * k + ph)), nB, LB)
        xC = noise + a_ * np.exp(1j * (2 * np.pi * fA * k + ph))
        c = (np.abs(xC.reshape(-1, budget)) ** 2).mean(axis=1) > thrC
        out["narrowband"].append((s_db, (pA > thrA).any(axis=1).mean(), (pB > thrB).any(axis=1).mean(), c.mean()))
    for s_db in snr_wb:
        x = noise * np.sqrt(1 + 10 ** (s_db / 10))  # a white, noise-like signal that fills the whole band
        out["wideband"].append((s_db,) + tuple(decisions(x)))

    names = ["A: FFT 256 x 8 avg", "B: FFT 2048, no avg", "C: time-domain energy"]
    colors = [BLUE, ORANGE, AQUA]
    markers = ["o", "s", "^"]
    with open(os.path.join(outdir, "baselines.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["signal", "snr_db_per_sample", "pd_A_fft256x8", "pd_B_fft2048", "pd_C_time_domain"])
        for sig, rws in out.items():
            for r in rws:
                w.writerow([sig, f"{r[0]:.0f}"] + [f"{v:.4f}" for v in r[1:]])

    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.8), sharey=True)
    for ax, (sig, ttl) in zip(axes, [("narrowband", "One narrowband tone"), ("wideband", "Wideband noise-like signal")]):
        arr = np.array(out[sig])
        for j in range(3):
            ax.plot(arr[:, 0], arr[:, 1 + j], marker=markers[j], color=colors[j], ms=4, lw=1.3, label=names[j])
        ax.axhline(0.9, color=GREY, lw=0.8, ls=":")
        ax.set_title(ttl, fontsize=10)
        ax.set_xlabel("Signal-to-noise power per sample (dB)", fontsize=9)
        style_axes(ax)
    axes[0].set_ylabel(f"Detection probability (band Pfa {band_pfa:g})", fontsize=9)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=3, fontsize=8, frameon=False)
    fig.tight_layout(rect=(0, 0.07, 1, 1))
    fig.savefig(os.path.join(outdir, "baselines.png"), dpi=200)
    plt.close(fig)

    lines = [f"Baseline comparison: {budget} samples per decision, {trials} decisions per point, "
             f"target band Pfa {band_pfa:g} (noise power assumed known).",
             f"  Measured band Pfa on noise only ({reps * trials} decisions): A {pf[0]:.4f}, B {pf[1]:.4f}, C {pf[2]:.4f}",
             "  SNR per sample needed for Pd = 0.9 (dB):"]
    for sig in ("narrowband", "wideband"):
        arr = np.array(out[sig])
        vals = [required_snr(arr[:, 0], arr[:, 1 + j]) for j in range(3)]
        lines.append(f"    {sig:10s}: A {vals[0]:6.2f}   B {vals[1]:6.2f}   C {vals[2]:6.2f}")
    return lines


def exp_new(outdir, rng, budget=2048, trials=4000, band_pfa=0.01):
    """Do the M4 features (CFAR threshold, total-power decision) fix the weaknesses found above?

    Same setting as exp_base: 2048 samples per decision, FFT 256 x 8 averages, band Pfa 0.01.
      M3 block : per-bin test only, Sidak threshold for the whole band.
      New block: decision_mode 2 (either). The band Pfa is shared: half for the per-bin test (CFAR,
                 Sidak) and half for the total-power test, so the sum stays about 0.01.
      Time dom.: one energy number over all samples (the best case for a wideband signal).
    The noise power is assumed known (1.0) in all three.
    Second part: the real noise power is wrong by -1 ... +1 dB (noise uncertainty).
    """
    n, L = 256, budget // 256
    half = band_pfa / 2

    def per_bin(pf):
        return 1 - (1 - pf) ** (1.0 / n)

    thrC = chi2.isf(band_pfa, 2 * budget) / (2 * budget)

    def three(x):
        p = energy_spectrum(x, n, L)
        _, old = detect(p, n, L, relative=False, pfa=per_bin(band_pfa), noise_power=1.0, decision_mode=0)
        _, new = detect(p, n, L, relative=False, pfa=per_bin(half), noise_power=1.0,
                        decision_mode=2, total_pfa=half)
        c = (np.abs(x.reshape(-1, budget)) ** 2).mean(axis=1) > thrC
        return old.mean(), new.mean(), c.mean()

    reps = 5
    pf = sum(np.array(three(cnoise(rng, trials * budget))) for _ in range(reps)) / reps

    k = np.arange(trials * budget)
    f0 = 32.5 / n
    noise = cnoise(rng, trials * budget)
    out = {"narrowband": [], "wideband": []}
    for s_db in np.arange(-30, -5, 1.0):
        a_ = np.sqrt(10 ** (s_db / 10))
        x = noise + a_ * np.exp(1j * (2 * np.pi * f0 * k + rng.uniform(0, 2 * np.pi)))
        out["narrowband"].append((s_db,) + three(x))
    for s_db in np.arange(-16, 5, 1.0):
        out["wideband"].append((s_db,) + three(noise * np.sqrt(1 + 10 ** (s_db / 10))))

    with open(os.path.join(outdir, "new_features.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["signal", "snr_db_per_sample", "pd_M3_block", "pd_new_block_either", "pd_time_domain"])
        for sig, rws in out.items():
            for r in rws:
                w.writerow([sig, f"{r[0]:.0f}"] + [f"{v:.4f}" for v in r[1:]])

    names = ["M3 block (per-bin test)", "New block (per-bin + total power)", "Time-domain energy"]
    colors, markers = [ORANGE, BLUE, AQUA], ["s", "o", "^"]
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.8), sharey=True)
    for ax, (sig, ttl) in zip(axes, [("narrowband", "One narrowband tone"), ("wideband", "Wideband noise-like signal")]):
        arr = np.array(out[sig])
        for j in range(3):
            ax.plot(arr[:, 0], arr[:, 1 + j], marker=markers[j], color=colors[j], ms=4, lw=1.3, label=names[j])
        ax.axhline(0.9, color=GREY, lw=0.8, ls=":")
        ax.set_title(ttl, fontsize=10)
        ax.set_xlabel("Signal-to-noise power per sample (dB)", fontsize=9)
        style_axes(ax)
    axes[0].set_ylabel(f"Detection probability (band Pfa {band_pfa:g})", fontsize=9)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=3, fontsize=7.5, frameon=False)
    fig.tight_layout(rect=(0, 0.08, 1, 1))
    fig.savefig(os.path.join(outdir, "new_features.png"), dpi=200)
    plt.close(fig)

    # Noise uncertainty on noise-only input: false-alarm rate when the real noise power is wrong.
    rows = []
    for e_db in (-1.0, -0.5, -0.25, 0.0, 0.25, 0.5, 1.0):
        x = cnoise(rng, trials * budget) * np.sqrt(10 ** (e_db / 10))
        p = energy_spectrum(x, n, L)
        _, bins = detect(p, n, L, relative=False, pfa=per_bin(band_pfa), noise_power=1.0, decision_mode=0)
        _, tot = detect(p, n, L, relative=False, pfa=per_bin(half), noise_power=1.0, decision_mode=1, total_pfa=half)
        _, rel = detect(p, n, L, relative=True, pfa=per_bin(band_pfa), decision_mode=0)
        rows.append((e_db, bins.mean(), rel.mean(), tot.mean()))
    with open(os.path.join(outdir, "new_noise_error.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["noise_error_db", "pfa_per_bin_test_known_noise", "pfa_per_bin_test_relative", "pfa_total_power"])
        for r in rows:
            w.writerow([f"{r[0]:+.2f}"] + [f"{v:.4f}" for v in r[1:]])

    # The block in a real flowgraph must give the same decisions as detect().
    x = cnoise(rng, 200 * budget) + 0.1 * np.exp(2j * np.pi * f0 * np.arange(200 * budget))
    tb = gr.top_block()
    det = energy_detector_cf(n, L, 6.0, False, per_bin(half), 1.0, 2, half)
    src = blocks.vector_source_c(x.astype(np.complex64).tolist(), False)
    s0, s1, s2 = blocks.null_sink(4 * n), blocks.null_sink(n), blocks.vector_sink_b(1)
    tb.connect(src, det)
    tb.connect((det, 0), s0)
    tb.connect((det, 1), s1)
    tb.connect((det, 2), s2)
    tb.run()
    _, ref = detect(energy_spectrum(x.astype(np.complex64), n, L), n, L, relative=False, pfa=per_bin(half),
                    noise_power=1.0, decision_mode=2, total_pfa=half)
    differ = int(np.sum(np.array(s2.data()).astype(bool) != ref))

    lines = [f"New features (M4): {budget} samples per decision, FFT {n} x {L}, band Pfa {band_pfa:g}, "
             f"{trials} decisions per point.",
             f"  Measured band Pfa on noise only ({reps * trials} decisions): M3 block {pf[0]:.4f}, "
             f"new block {pf[1]:.4f}, time-domain {pf[2]:.4f}",
             "  SNR per sample needed for Pd = 0.9 (dB):"]
    for sig in ("narrowband", "wideband"):
        arr = np.array(out[sig])
        v = [required_snr(arr[:, 0], arr[:, 1 + j]) for j in range(3)]
        lines.append(f"    {sig:10s}: M3 block {v[0]:6.2f}   new block {v[1]:6.2f}   time-domain {v[2]:6.2f}")
    lines += ["  False-alarm rate when the real noise power is wrong (noise only):",
              "    error dB | per-bin, known noise | per-bin, relative | total power"]
    for r in rows:
        lines.append(f"    {r[0]:+6.2f}   |       {r[1]:.4f}         |      {r[2]:.4f}       |   {r[3]:.4f}")
    lines.append(f"  Block in flowgraph against detect(): differing decisions {differ} / {len(ref)}")
    return lines


# --------------------------------------------------------------------------- main
def environment_lines():
    import scipy
    return [f"Machine: {platform.platform()}", f"Python {platform.python_version()}, numpy {np.__version__}, "
            f"scipy {scipy.__version__}, GNU Radio {gr.version()}"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("outdir")
    ap.add_argument("--only", choices=["bench", "noise", "base", "new"])
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    rng = np.random.default_rng(20261)
    lines = ["gr-isac M4 experiments", *environment_lines(), ""]
    if args.only in (None, "noise"):
        lines += exp_noise(args.outdir, rng) + [""]
    if args.only in (None, "base"):
        lines += exp_base(args.outdir, rng) + [""]
    if args.only in (None, "new"):
        lines += exp_new(args.outdir, rng) + [""]
    if args.only in (None, "bench"):
        b, _ = exp_bench(args.outdir)
        lines += b + [""]
    text = "\n".join(lines)
    print(text)
    name = "m4_summary.txt" if args.only is None else f"m4_summary_{args.only}.txt"
    with open(os.path.join(args.outdir, name), "w") as f:
        f.write(text + "\n")


if __name__ == "__main__":
    main()
