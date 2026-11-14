#!/usr/bin/env python3
# Copyright 2026 Group 05 ECE310.
# SPDX-License-Identifier: GPL-3.0-or-later
"""Terminal demo of the OFDM ISAC chain (no hardware, about 5 seconds).

Builds the flowgraph  ofdm_isac_tx_c -> moving_target_channel_cc -> ofdm_sensing_rx_ccf,
sends 3 frames through a channel with two moving targets, prints the detected targets next to
the true ones, then uses the same signal as a data link and prints the bit error rate.

Usage:  python3 demo_isac.py [OUTPUT_PNG]
"""
import sys

import numpy as np
import pmt
from gnuradio import blocks, gr
from gnuradio.isac import (moving_target_channel_cc, ofdm_isac, ofdm_isac_tx_c, ofdm_sensing_rx_ccf)

TARGETS = [(300.0, 20.0, -20.0), (500.0, -45.0, -26.0)]  # range m, speed m/s, echo gain dB


def main():
    png = sys.argv[1] if len(sys.argv) > 1 else None
    cfg = ofdm_isac.OfdmConfig()
    print("OFDM ISAC demo: %d subcarriers (%d used), %d symbols per frame, %.1f MHz bandwidth, %.1f GHz carrier"
          % (cfg.fft_size, cfg.n_active, cfg.n_symbols, cfg.bandwidth / 1e6, cfg.carrier_freq / 1e9))
    print("Theory: range resolution %.1f m, speed resolution %.2f m/s" % (cfg.range_resolution, cfg.speed_resolution))
    print("True targets (range m, speed m/s, echo gain dB):", TARGETS)

    n_frames = 3
    tb = gr.top_block()
    tx = ofdm_isac_tx_c(seed=1)
    head = blocks.head(gr.sizeof_gr_complex, n_frames * cfg.frame_len)
    ch = moving_target_channel_cc(targets=TARGETS, noise_power=1e-3)
    rx = ofdm_sensing_rx_ccf(latency_samples=ch.latency)
    sink = blocks.vector_sink_f(rx.n_doppler * rx.n_range)
    dbg = blocks.message_debug()
    tb.connect(tx, head, ch, (rx, 0))
    tb.connect(head, (rx, 1))
    tb.connect(rx, sink)
    tb.msg_connect(rx, "targets", dbg, "store")
    tb.run()

    for k in range(dbg.num_messages()):
        msg = pmt.to_python(dbg.get_message(k))
        print("frame %d: %d target(s) detected" % (k, len(msg)))
        for r, v, snr in msg:
            print("    range %7.1f m   speed %+7.2f m/s   SNR %5.1f dB" % (r, v, snr))

    # the same signal as a data link (one-way path, 180 m, 15 m/s, noise 10 dB below the signal)
    tx2 = np.concatenate([tx.frame(k) for k in range(3)])
    rng = np.random.default_rng(5)
    y = ofdm_isac.simulate_channel(tx2, [(180.0, 15.0, 0.0)], cfg, 0.1, rng, one_way=True)
    errors = bits = 0
    for k in range(2):
        seg = y[k * cfg.frame_len:(k + 2) * cfg.frame_len]
        start = ofdm_isac.estimate_timing_cp(seg, cfg)
        e, b = ofdm_isac.bit_errors(seg, cfg, 1, k, start)
        errors += e
        bits += b
    print("Communication at 10 dB SNR: %d bit errors in %d bits (BER %.1e)" % (errors, bits, errors / bits))

    if png:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        db = np.array(sink.data())[-rx.n_doppler * rx.n_range:].reshape(rx.n_doppler, rx.n_range)
        r, v = rx.range_axis(), rx.speed_axis()
        fig, ax = plt.subplots(figsize=(6.5, 3.8))
        ax.imshow(db, aspect="auto", origin="lower", cmap="viridis", vmin=db.max() - 70, vmax=db.max(),
                  extent=[r[0], r[-1], v[0], v[-1]])
        ax.set_xlabel("range (m)")
        ax.set_ylabel("speed (m/s)")
        fig.tight_layout()
        fig.savefig(png, dpi=150)
        print("map saved to", png)


if __name__ == "__main__":
    main()
