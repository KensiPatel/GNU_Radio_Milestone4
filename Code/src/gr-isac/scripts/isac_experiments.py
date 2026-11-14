#!/usr/bin/env python3
# Copyright 2026 Group 05 ECE310.
# SPDX-License-Identifier: GPL-3.0-or-later
"""Experiments for the OFDM ISAC chain (simulation only, fixed random seeds).

  map       Example range-Doppler map with two moving targets (figure).
  theory    Resolution and limits from the formulas, next to the measured ones:
            two-target resolution in range and in speed.
  accuracy  Detection probability, range and speed error against the echo SNR,
            for the OFDM method and for two correlation-radar baselines.
  pfa       Is the CFAR false-alarm rate really the requested one?
  ber       The same signal used as a data link: bit error rate against SNR, with and
            without an extra echo, against the QPSK theory.
  ici       Dynamic range of the map for fast targets (inter-carrier interference).
  speed     How fast the Python receiver is compared with real time.

Usage:  python3 isac_experiments.py OUTPUT_DIR [--only map|theory|accuracy|pfa|ber|ici|speed]
Needs GNU Radio with gr-isac installed, numpy, scipy, matplotlib. About 5 to 8 minutes.
"""
import argparse
import csv
import math
import os
import time

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import gamma

from gnuradio import blocks, gr
from gnuradio.isac import (moving_target_channel_cc, ofdm_isac, ofdm_isac_tx_c, ofdm_sensing_rx_ccf)

BLUE, ORANGE, AQUA, INK, GREY = "#2a78d6", "#eb6834", "#1baf7a", "#0b0b0b", "#52514e"
CFG = ofdm_isac.OfdmConfig()
LAT = 16  # latency of the 33-tap delay filter (samples)
NOISE = 1e-3  # noise power in the radar experiments; the transmitted power is 1.0
SEED_TX = 1
TX = ofdm_isac_tx_c(seed=SEED_TX)
TX0 = np.concatenate((TX.frame(0), TX.frame(1)[:200]))  # frame 0 plus a little extra for the baselines
C0 = ofdm_isac.C0


def style_axes(ax):
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.grid(alpha=0.25, lw=0.6)
    ax.tick_params(labelsize=8)


def write_csv(path, header, rows):
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(header)
        w.writerows(rows)


def echo(targets, rng, noise=NOISE, tx=TX0):
    return ofdm_isac.simulate_channel(tx, targets, CFG, noise, rng)


def sense(rx, pfa=1e-8, use_window=True):
    n = CFG.frame_len
    return ofdm_isac.sense_frame(rx[:n], TX0[:n], CFG, pfa=pfa, use_window=use_window, latency_samples=LAT)


# ------------------------------------------------------------------ baselines
def lag_correlation(rx, lags):
    """r[m, l] = sum over OFDM symbol m of rx[n + l] * conj(tx[n]), the matched filter of one symbol."""
    n, m, s = CFG.frame_len, CFG.n_symbols, CFG.symbol_len
    tx_sym = TX0[:n].reshape(m, s)
    out = np.empty((m, len(lags)), dtype=np.complex128)
    for k, lag in enumerate(lags):
        out[:, k] = np.sum(rx[lag:lag + n].reshape(m, s) * np.conj(tx_sym), axis=1)
    return out


LAGS = np.arange(LAT, CFG.cp_len)  # delays 0 .. cp_len - LAT samples of echo delay


def parabola_peak(p, i):
    if 0 < i < len(p) - 1:
        a, b, c = np.log(p[i - 1] + 1e-30), np.log(p[i] + 1e-30), np.log(p[i + 1] + 1e-30)
        den = a - 2 * b + c
        return i + (0.0 if den == 0 else 0.5 * (a - c) / den)
    return float(i)


def baseline_a(rx, pfa_cell):
    """Correlation radar, non-coherent: |matched filter of one symbol|^2 summed over all symbols."""
    r = lag_correlation(rx, LAGS)
    p = np.sum(np.abs(r) ** 2, axis=0)
    # noise known (best case for the baseline): each term is exponential with mean NOISE * symbol energy
    e_sym = np.sum(np.abs(TX0[:CFG.symbol_len]) ** 2)
    thr = NOISE * e_sym * gamma.isf(pfa_cell, a=CFG.n_symbols)
    i = int(np.argmax(p))
    return (p[i] > thr), parabola_peak(p, i) * C0 / (2 * CFG.sample_rate)


def baseline_b(rx, pfa_cell):
    """Correlation radar, coherent over the whole frame, no Doppler search."""
    r = np.sum(lag_correlation(rx, LAGS), axis=0)
    p = np.abs(r) ** 2
    e_frame = np.sum(np.abs(TX0[:CFG.frame_len]) ** 2)
    thr = NOISE * e_frame * (-math.log(pfa_cell))
    i = int(np.argmax(p))
    return (p[i] > thr), parabola_peak(p, i) * C0 / (2 * CFG.sample_rate)


# ------------------------------------------------------------------ map
def exp_map(out):
    rng = np.random.default_rng(11)
    targets = [(300.0, 20.0, -20.0), (450.0, -60.0, -24.0), (150.0, 0.0, -30.0)]
    rx = echo(targets, rng)
    found, _ = sense(rx)
    Y = ofdm_isac.demodulate_symbols(rx[:CFG.frame_len], CFG)
    X = ofdm_isac.demodulate_symbols(TX0[:CFG.frame_len], CFG)
    cm = ofdm_isac.range_doppler_map(Y, X, CFG, True, 2, 2, LAT)
    db = 10 * np.log10(np.abs(cm) ** 2 + 1e-30)
    ranges, speeds = ofdm_isac.axes(CFG, 2, 2, LAT)
    fig, ax = plt.subplots(figsize=(6.8, 4.0))
    im = ax.imshow(db, aspect="auto", origin="lower", cmap="viridis", vmin=db.max() - 70, vmax=db.max(),
                   extent=[ranges[0], ranges[-1], speeds[0], speeds[-1]])
    for (r, v, g), mark in zip(targets, ("o", "s", "^")):
        ax.plot(r, v, mark, mfc="none", mec="white", ms=13, mew=1.4)
    for r, v, s in found:
        ax.plot(r, v, "x", color=ORANGE, ms=7, mew=1.6)
    ax.set_xlabel("range (m)")
    ax.set_ylabel("speed (m/s), positive = towards the receiver")
    ax.set_title("Range-Doppler map, 3 targets (white = true, orange x = detected)", fontsize=9)
    cb = fig.colorbar(im, ax=ax)
    cb.set_label("dB", fontsize=8)
    cb.ax.tick_params(labelsize=7)
    ax.tick_params(labelsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(out, "range_doppler_map.png"), dpi=200)
    plt.close(fig)
    rows = [("true", t[0], t[1], t[2]) for t in targets] + [("detected", r, v, s) for r, v, s in found]
    write_csv(os.path.join(out, "map_example.csv"), ["kind", "range_m", "speed_mps", "gain_db_or_snr_db"], rows)
    print("map example: true", [(t[0], t[1]) for t in targets], "detected",
          [(round(r, 1), round(v, 1)) for r, v, _ in found])


# ------------------------------------------------------------------ theory / resolution
def exp_theory(out):
    c = CFG
    print("theory: B = %.2f MHz, range resolution c/2B = %.2f m, speed resolution = %.2f m/s" %
          (c.bandwidth / 1e6, c.range_resolution, c.speed_resolution))
    print("        max range (CP, minus filter latency) = %.0f m, max speed = %.0f m/s" %
          (C0 * (c.cp_len - LAT) / (2 * c.sample_rate), c.max_speed))
    # Two equally strong targets. The second one gets a random extra distance of 0 .. one half
    # wavelength (a few cm): this changes the carrier phase between the two echoes at random, as in
    # reality. Without it the phase difference would be the same in every trial.
    rng = np.random.default_rng(21)
    trials = 60
    fr = [0.5, 0.75, 1.0, 1.25, 1.5, 2.0, 2.5, 3.0]
    rows = []
    for win, wname in ((True, "Hann"), (False, "rectangular")):
        for axis in ("range", "speed"):
            for f in fr:
                ok = 0
                for _ in range(trials):
                    r0 = rng.uniform(100, 400)
                    v0 = rng.uniform(-100, 100)
                    jit = rng.uniform(0, c.wavelength / 2)
                    if axis == "range":
                        t = [(r0, v0, -20.0), (r0 + f * c.range_resolution + jit, v0, -20.0)]
                    else:
                        t = [(r0, v0, -20.0), (r0 + jit, v0 + f * c.speed_resolution, -20.0)]
                    found, _ = sense(echo(t, rng), use_window=win)
                    good = len(found) == 2
                    if good:
                        for tr, tv, _g in t:
                            if not any(abs(d[0] - tr) < 0.5 * c.range_resolution
                                       and abs(d[1] - tv) < 0.5 * c.speed_resolution for d in found):
                                good = False
                    ok += int(good)
                rows.append((wname, axis, f, ok / trials))
                print("  %-11s two targets, %s separation %.2f cells: both found in %3.0f %% of trials" %
                      (wname, axis, f, 100 * ok / trials))
    write_csv(os.path.join(out, "resolution.csv"),
              ["window", "axis", "separation_in_resolution_cells", "fraction_resolved"], rows)
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.4), sharey=True)
    for ax, axis in zip(axes, ("range", "speed")):
        for wname, col, mk in (("Hann", BLUE, "o"), ("rectangular", ORANGE, "s")):
            pts = [(r[2], r[3]) for r in rows if r[0] == wname and r[1] == axis]
            ax.plot([p[0] for p in pts], [p[1] for p in pts], marker=mk, color=col, lw=1.4, ms=5, label=wname + " window")
        ax.axvline(1.0, color=GREY, ls=":", lw=1)
        ax.text(1.04, 0.04, "formula: 1 cell", fontsize=7, color=GREY)
        ax.set_title("separation in %s" % axis, fontsize=9)
        ax.set_xlabel("separation (resolution cells from the formula)", fontsize=8)
        style_axes(ax)
    axes[0].set_ylabel("both targets found", fontsize=8)
    axes[0].legend(fontsize=7.5, frameon=False, loc="center right")
    fig.tight_layout()
    fig.savefig(os.path.join(out, "resolution.png"), dpi=200)
    plt.close(fig)


# ------------------------------------------------------------------ accuracy and Pd
def exp_accuracy(out):
    snrs = list(range(-38, -9, 2))  # echo SNR per sample before processing, dB (noise power 1e-3)
    trials = 60
    rng = np.random.default_rng(31)
    pfa_frame = CFG.n_symbols * 2 * (CFG.cp_len - LAT) * 2 * 1e-8
    pfa_cell_base = pfa_frame / len(LAGS)
    rows = []
    for speed_case in (0.0, 100.0):
        for snr in snrs:
            g = snr + 10 * math.log10(NOISE)
            hit_o = hit_o_rs = hit_a = hit_b = 0
            e_r, e_v, e_ra, e_rb = [], [], [], []
            for _ in range(trials):
                r0 = rng.uniform(60, 650)
                v0 = speed_case if speed_case == 0.0 else rng.uniform(80, 120) * rng.choice([-1, 1])
                rx = echo([(r0, v0, g)], rng)
                found, _ = sense(rx)
                near = [t for t in found if abs(t[0] - r0) < CFG.range_resolution]
                if near:
                    hit_o += 1
                    best = min(near, key=lambda t: abs(t[0] - r0))
                    e_r.append(best[0] - r0)
                    if abs(best[1] - v0) < CFG.speed_resolution:
                        hit_o_rs += 1
                        e_v.append(best[1] - v0)
                d, rb = baseline_a(rx, pfa_cell_base)
                if d and abs(rb - r0) < CFG.range_resolution:
                    hit_a += 1
                    e_ra.append(rb - r0)
                d, rb = baseline_b(rx, pfa_cell_base)
                if d and abs(rb - r0) < CFG.range_resolution:
                    hit_b += 1
                    e_rb.append(rb - r0)
            rms = lambda e: float(np.sqrt(np.mean(np.square(e)))) if e else float("nan")  # noqa: E731
            rows.append((speed_case, snr, hit_o / trials, hit_o_rs / trials, hit_a / trials, hit_b / trials,
                         rms(e_r), rms(e_v), rms(e_ra), rms(e_rb)))
            print("  speed %5.0f echo SNR %4d dB: Pd OFDM %.2f (range+speed %.2f)  corr.non-coh. %.2f  corr.coherent %.2f"
                  "  rms range %.1f m speed %.2f m/s" % (speed_case, snr, hit_o / trials, hit_o_rs / trials,
                                                         hit_a / trials, hit_b / trials, rms(e_r), rms(e_v)))
    write_csv(os.path.join(out, "accuracy_pd.csv"),
              ["target_speed_mps", "echo_snr_db", "pd_ofdm_range", "pd_ofdm_range_and_speed", "pd_corr_noncoherent",
               "pd_corr_coherent", "rms_range_ofdm_m", "rms_speed_ofdm_mps", "rms_range_corr_noncoh_m",
               "rms_range_corr_coh_m"], rows)
    arr = np.array(rows, dtype=float)
    fig, axes = plt.subplots(1, 2, figsize=(7.4, 3.6), sharey=True)
    for ax, sp, title in zip(axes, (0.0, 100.0), ("static target (0 m/s)", "moving target (about 100 m/s)")):
        d = arr[arr[:, 0] == sp]
        ax.plot(d[:, 1], d[:, 2], "o-", color=BLUE, ms=4, lw=1.4, label="OFDM range-Doppler (ours)")
        ax.plot(d[:, 1], d[:, 5], "^-", color=AQUA, ms=4, lw=1.4, label="correlation, coherent, no Doppler")
        ax.plot(d[:, 1], d[:, 4], "s-", color=ORANGE, ms=4, lw=1.4, label="correlation, non-coherent")
        ax.axhline(0.9, color=GREY, ls=":", lw=0.8)
        ax.set_title(title, fontsize=9)
        ax.set_xlabel("echo SNR per sample before processing (dB)", fontsize=8)
        style_axes(ax)
    axes[0].set_ylabel("probability of detection (range correct)", fontsize=8)
    axes[0].legend(fontsize=7, frameon=False, loc="lower right")
    fig.tight_layout()
    fig.savefig(os.path.join(out, "pd_vs_snr.png"), dpi=200)
    plt.close(fig)
    fig, axes = plt.subplots(1, 2, figsize=(7.4, 3.4))
    d = arr[arr[:, 0] == 100.0]
    axes[0].semilogy(d[:, 1], d[:, 6], "o-", color=BLUE, ms=4, lw=1.4)
    axes[0].axhline(CFG.range_resolution, color=GREY, ls=":", lw=1)
    axes[0].text(d[0, 1], CFG.range_resolution * 1.1, "one resolution cell", fontsize=7, color=GREY)
    axes[0].set_ylabel("rms range error (m)", fontsize=8)
    axes[1].semilogy(d[:, 1], d[:, 7], "o-", color=ORANGE, ms=4, lw=1.4)
    axes[1].axhline(CFG.speed_resolution, color=GREY, ls=":", lw=1)
    axes[1].text(d[0, 1], CFG.speed_resolution * 1.1, "one resolution cell", fontsize=7, color=GREY)
    axes[1].set_ylabel("rms speed error (m/s)", fontsize=8)
    for ax in axes:
        ax.set_xlabel("echo SNR per sample before processing (dB)", fontsize=8)
        style_axes(ax)
    fig.tight_layout()
    fig.savefig(os.path.join(out, "accuracy_vs_snr.png"), dpi=200)
    plt.close(fig)


# ------------------------------------------------------------------ pfa calibration
def exp_pfa(out):
    rng = np.random.default_rng(41)
    X, _ = ofdm_isac.frame_symbols(CFG, SEED_TX, 0)
    frames = 40
    rows = []
    variants = (("no window, no padding", False, 1), ("Hann, no padding", True, 1), ("Hann, padding 2 (default)", True, 2))
    for name, win, pad in variants:
        maps = []
        for _ in range(frames):
            Y = (rng.standard_normal(X.shape) + 1j * rng.standard_normal(X.shape)) * 0.1
            maps.append(np.abs(ofdm_isac.range_doppler_map(Y, X, CFG, win, pad, pad, LAT)) ** 2)
        for pfa in (1e-2, 1e-3, 1e-4):
            hits = cells = 0
            for p in maps:
                det, _ = ofdm_isac.ca_cfar_2d(p, pfa)
                hits += int(det.sum())
                cells += det.size
            rows.append((name, pfa, hits / cells, hits / cells / pfa))
            print("  %-26s requested %.0e measured %.2e  ratio %.2f" % (name, pfa, hits / cells, hits / cells / pfa))
    write_csv(os.path.join(out, "cfar_calibration.csv"), ["processing", "pfa_requested", "pfa_measured", "ratio"], rows)


# ------------------------------------------------------------------ ber
def run_ber(snr_db, extra_echo, n_frames, rng):
    """Returns (errors with pilot channel estimate, errors with perfect channel knowledge, bits)."""
    tx = np.concatenate([TX.frame(k) for k in range(n_frames + 1)])
    noise = 10 ** (-snr_db / 10.0)
    paths = [(180.0, 15.0, 0.0)] + ([(270.0, 20.0, -15.0)] if extra_echo else [])
    clean = ofdm_isac.simulate_channel(tx, paths, CFG, 0.0, rng, one_way=True)
    s = math.sqrt(noise / 2.0)
    y = clean + s * (rng.standard_normal(clean.size) + 1j * rng.standard_normal(clean.size))
    e_est = e_ideal = bits = 0
    d = CFG.data_bins
    for k in range(n_frames):
        sl = slice(k * CFG.frame_len, (k + 2) * CFG.frame_len)
        start = ofdm_isac.estimate_timing_cp(y[sl], CFG)
        e, b = ofdm_isac.bit_errors(y[sl], CFG, SEED_TX, k, start)
        e_est += e
        bits += b
        X, sent = ofdm_isac.frame_symbols(CFG, SEED_TX, k)
        Yc = ofdm_isac.demodulate_symbols(clean[sl], CFG, start)
        Yn = ofdm_isac.demodulate_symbols(y[sl], CFG, start)
        H = Yc[:, d] / X[:, d]
        z = Yn[:, d] / H
        dec = np.stack(((z.real < 0).astype(int), (z.imag < 0).astype(int)), axis=-1)
        e_ideal += int(np.sum(dec != sent))
    return e_est, e_ideal, bits


def q_func(x):
    return 0.5 * math.erfc(x / math.sqrt(2.0))


def exp_ber(out):
    rng = np.random.default_rng(51)
    snrs = [0, 2, 4, 6, 8, 10, 12]
    rows = []
    for s in snrs:
        r_los = run_ber(s, False, 10, rng)
        r_ech = run_ber(s, True, 10, rng)
        es_n0 = 10 ** ((s + 10 * math.log10(CFG.fft_size / CFG.n_active)) / 10.0)
        theory = q_func(math.sqrt(es_n0))
        rows.append((s, r_los[0] / r_los[2], r_ech[0] / r_ech[2], r_los[1] / r_los[2], theory))
        print("  SNR %2d dB: BER pilot-estimated %.2e | with echo -15 dB %.2e | perfect CSI %.2e | QPSK theory %.2e"
              % (s, rows[-1][1], rows[-1][2], rows[-1][3], theory))
    write_csv(os.path.join(out, "ber.csv"),
              ["snr_in_db", "ber_pilot_estimate", "ber_pilot_estimate_plus_echo", "ber_perfect_csi",
               "ber_qpsk_theory"], rows)
    arr = np.array(rows, dtype=float)
    floor = 3e-6
    fig, ax = plt.subplots(figsize=(5.8, 3.8))
    ax.semilogy(arr[:, 0], np.maximum(arr[:, 4], floor), "--", color=GREY, lw=1.4, label="QPSK theory")
    ax.semilogy(arr[:, 0], np.maximum(arr[:, 3], floor), "^-", color=AQUA, ms=4, lw=1.4,
                label="measured, perfect channel knowledge")
    ax.semilogy(arr[:, 0], np.maximum(arr[:, 1], floor), "o-", color=BLUE, ms=4, lw=1.4,
                label="measured, channel from pilots")
    ax.semilogy(arr[:, 0], np.maximum(arr[:, 2], floor), "s-", color=ORANGE, ms=4, lw=1.4,
                label="pilots + extra echo (-15 dB)")
    ax.set_xlabel("SNR per sample (dB)", fontsize=9)
    ax.set_ylabel("bit error rate", fontsize=9)
    ax.set_ylim(floor, 1)
    ax.legend(fontsize=7.5, frameon=False)
    style_axes(ax)
    fig.tight_layout()
    fig.savefig(os.path.join(out, "ber.png"), dpi=200)
    plt.close(fig)


# ------------------------------------------------------------------ ici
def exp_ici(out):
    rng = np.random.default_rng(61)
    rows = []
    for v in (0.0, 20.0, 60.0, 120.0, 250.0, 350.0):
        gaps = []
        for _ in range(6):
            rx = echo([(300.0, v, -10.0)], rng)
            Y = ofdm_isac.demodulate_symbols(rx[:CFG.frame_len], CFG)
            X = ofdm_isac.demodulate_symbols(TX0[:CFG.frame_len], CFG)
            db = 10 * np.log10(np.abs(ofdm_isac.range_doppler_map(Y, X, CFG, True, 2, 2, LAT)) ** 2 + 1e-30)
            gaps.append(db.max() - np.median(db))
        rows.append((v, float(np.mean(gaps))))
        print("  target speed %5.0f m/s: peak above map median = %.1f dB" % (v, np.mean(gaps)))
    write_csv(os.path.join(out, "dynamic_range_vs_speed.csv"), ["speed_mps", "peak_minus_median_db"], rows)


# ------------------------------------------------------------------ speed
def exp_speed(out):
    rng = np.random.default_rng(71)
    rx = echo([(300.0, 20.0, -20.0)], rng)
    t0 = time.perf_counter()
    reps = 20
    for _ in range(reps):
        sense(rx)
    per_frame = (time.perf_counter() - t0) / reps
    frame_time = CFG.frame_len / CFG.sample_rate
    # whole flowgraph
    n_frames = 20
    tb = gr.top_block()
    tx = ofdm_isac_tx_c(seed=SEED_TX)
    head = blocks.head(gr.sizeof_gr_complex, n_frames * CFG.frame_len)
    ch = moving_target_channel_cc(targets=[(300.0, 20.0, -20.0)], noise_power=NOISE)
    rxb = ofdm_sensing_rx_ccf(latency_samples=ch.latency)
    tb.connect(tx, head, ch, (rxb, 0))
    tb.connect(head, (rxb, 1))
    tb.connect(rxb, blocks.null_sink(gr.sizeof_float * rxb.n_doppler * rxb.n_range))
    t0 = time.perf_counter()
    tb.run()
    total = time.perf_counter() - t0
    rate = n_frames * CFG.frame_len / total
    print("  sensing core: %.1f ms per frame, the frame lasts %.2f ms at %.0f MS/s"
          " -> takes %.2f x the signal duration (1.0 = just real time)"
          % (per_frame * 1e3, frame_time * 1e3, CFG.sample_rate / 1e6, per_frame / frame_time))
    print("  whole flowgraph (tx, channel, sensing): %.2f MS/s -> takes %.1f x the signal duration at 10 MS/s"
          % (rate / 1e6, CFG.sample_rate / rate))
    write_csv(os.path.join(out, "speed.csv"), ["item", "value"],
              [("sensing_ms_per_frame", per_frame * 1e3), ("frame_duration_ms", frame_time * 1e3),
               ("flowgraph_msps", rate / 1e6)])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("outdir")
    ap.add_argument("--only", choices=["map", "theory", "accuracy", "pfa", "ber", "ici", "speed"])
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    steps = [("map", exp_map), ("theory", exp_theory), ("accuracy", exp_accuracy), ("pfa", exp_pfa),
             ("ber", exp_ber), ("ici", exp_ici), ("speed", exp_speed)]
    for name, fn in steps:
        if args.only in (None, name):
            print("== %s" % name)
            t0 = time.perf_counter()
            fn(args.outdir)
            print("   (%.0f s)" % (time.perf_counter() - t0))


if __name__ == "__main__":
    main()
