#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright 2026 Group 05 ECE310.
#
# SPDX-License-Identifier: GPL-3.0-or-later
#
"""QA tests of the OFDM ISAC chain: transmitter, moving-target channel, sensing receiver."""

import math

import numpy
import pmt
from gnuradio import blocks, gr, gr_unittest

try:
    from gnuradio.isac import (energy_detector_cf, moving_target_channel_cc, ofdm_isac,
                               ofdm_isac_tx_c, ofdm_sensing_rx_ccf)
except ImportError:
    import os
    import sys

    sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "build", "test_modules"))
    from gnuradio.isac import (energy_detector_cf, moving_target_channel_cc, ofdm_isac,
                               ofdm_isac_tx_c, ofdm_sensing_rx_ccf)

CFG = ofdm_isac.OfdmConfig()


def run_chain(targets, noise_power=1e-3, n_frames=2, seed=1, rx_kwargs=None, one_way=False):
    """tx -> channel -> sensing receiver; the receiver also gets a second copy of the tx stream."""
    n = CFG.frame_len * n_frames
    tb = gr.top_block()
    tx = ofdm_isac_tx_c(seed=seed)
    ref = ofdm_isac_tx_c(seed=seed)
    ch = moving_target_channel_cc(targets=targets, noise_power=noise_power, one_way=one_way)
    kw = dict(latency_samples=ch.latency)
    kw.update(rx_kwargs or {})
    rx = ofdm_sensing_rx_ccf(**kw)
    sink = blocks.vector_sink_f(rx.n_doppler * rx.n_range)
    dbg = blocks.message_debug()
    tb.connect(tx, blocks.head(gr.sizeof_gr_complex, n), ch, (rx, 0))
    tb.connect(ref, blocks.head(gr.sizeof_gr_complex, n), (rx, 1))
    tb.connect(rx, sink)
    tb.msg_connect(rx, "targets", dbg, "store")
    tb.run()
    rx.messages = dbg
    return rx, numpy.array(sink.data())


class qa_ofdm_isac(gr_unittest.TestCase):

    def test_001_theory_numbers(self):
        """Resolution and limits follow the textbook formulas."""
        c = CFG
        self.assertAlmostEqual(c.bandwidth, 192 * 10e6 / 256)
        self.assertAlmostEqual(c.range_resolution, ofdm_isac.C0 / (2 * c.bandwidth))
        self.assertAlmostEqual(c.speed_resolution, c.wavelength / (2 * 128 * c.symbol_time))
        self.assertAlmostEqual(c.max_range, ofdm_isac.C0 * 64 / (2 * 10e6))
        self.assertAlmostEqual(c.max_speed, c.wavelength / (4 * c.symbol_time))
        self.assertAlmostEqual(c.range_resolution, 19.99, places=1)
        self.assertAlmostEqual(c.speed_resolution, 6.31, places=1)

    def test_002_bad_arguments(self):
        """Wrong settings are refused."""
        with self.assertRaises(ValueError):
            ofdm_isac.OfdmConfig(fft_size=100)
        with self.assertRaises(ValueError):
            ofdm_isac.OfdmConfig(n_active=257)
        with self.assertRaises(ValueError):
            ofdm_isac.OfdmConfig(cp_len=0)
        with self.assertRaises(ValueError):
            ofdm_sensing_rx_ccf(pfa=0.0)
        with self.assertRaises(ValueError):
            ofdm_sensing_rx_ccf(range_pad=0)
        with self.assertRaises(ValueError):
            moving_target_channel_cc(n_taps=32)
        with self.assertRaises(ValueError):
            moving_target_channel_cc(noise_power=-1.0)
        with self.assertRaises(ValueError):
            moving_target_channel_cc(targets=[(5000.0, 0.0, -20.0)], max_range=1000.0)
        with self.assertRaises(ValueError):
            ofdm_isac.ca_cfar_2d(numpy.ones((300, 100)), 1.5)

    def test_003_frame_structure(self):
        """Unit power, cyclic prefix is a copy of the tail, demodulation gives the symbols back."""
        X, bits = ofdm_isac.frame_symbols(CFG, 3, 0)
        x = ofdm_isac.modulate(X, CFG)
        self.assertEqual(x.size, CFG.frame_len)
        self.assertAlmostEqual(float(numpy.mean(numpy.abs(x) ** 2)), 1.0, delta=0.02)
        sym = x.reshape(CFG.n_symbols, CFG.symbol_len)
        self.assertTrue(numpy.allclose(sym[:, :CFG.cp_len], sym[:, -CFG.cp_len:]))
        back = ofdm_isac.demodulate_symbols(x, CFG)
        self.assertLess(float(numpy.abs(back - X).max()), 1e-5)
        self.assertEqual(bits.shape, (CFG.n_symbols, CFG.data_bins.size, 2))
        self.assertEqual(float(numpy.abs(X[:, 0]).max()), 0.0)  # DC unused

    def test_004_tx_block(self):
        """The transmitter block streams frame 0, then frame 1, and the seed decides the data."""
        tb = gr.top_block()
        src = ofdm_isac_tx_c(seed=5)
        head = blocks.head(gr.sizeof_gr_complex, 2 * CFG.frame_len)
        sink = blocks.vector_sink_c()
        tb.connect(src, head, sink)
        tb.run()
        out = numpy.array(sink.data())
        self.assertEqual(out.size, 2 * CFG.frame_len)
        self.assertTrue(numpy.allclose(out[:CFG.frame_len], src.frame(0)))
        self.assertTrue(numpy.allclose(out[CFG.frame_len:], src.frame(1)))
        self.assertFalse(numpy.allclose(src.frame(0), src.frame(1)))
        self.assertFalse(numpy.allclose(src.frame(0), ofdm_isac_tx_c(seed=6).frame(0)))

    def test_005_channel_noise_power(self):
        """With no target the output is noise of the requested power."""
        tb = gr.top_block()
        src = blocks.null_source(gr.sizeof_gr_complex)
        ch = moving_target_channel_cc(targets=[], noise_power=0.25)
        sink = blocks.vector_sink_c()
        tb.connect(src, blocks.head(gr.sizeof_gr_complex, 200000), ch, sink)
        tb.run()
        out = numpy.array(sink.data())
        self.assertAlmostEqual(float(numpy.mean(numpy.abs(out) ** 2)), 0.25, delta=0.01)

    def test_006_channel_doppler_sign_and_gain(self):
        """A target moving towards us shifts a tone up by 2 v fc / c, and the echo gain is right."""
        v, g = 30.0, -10.0
        tb = gr.top_block()
        src = blocks.vector_source_c([1.0 + 0j] * 40000)
        ch = moving_target_channel_cc(targets=[(200.0, v, g)], noise_power=0.0)
        sink = blocks.vector_sink_c()
        tb.connect(src, ch, sink)
        tb.run()
        y = numpy.array(sink.data())[2000:34000]
        spec = numpy.abs(numpy.fft.fft(y * numpy.hanning(y.size), 1 << 20))
        f = numpy.fft.fftfreq(1 << 20, 1.0 / 10e6)[spec.argmax()]
        expected = 2 * v * CFG.carrier_freq / ofdm_isac.C0
        self.assertAlmostEqual(f, expected, delta=20.0)
        self.assertGreater(f, 0.0)
        self.assertAlmostEqual(10 * math.log10(float(numpy.mean(numpy.abs(y) ** 2))), g, delta=0.3)

    def test_007_channel_delay(self):
        """The echo of an impulse arrives at 2 R / c plus the filter latency, also between samples."""
        for rng_m, frac_ok in ((150.0, True), (157.5, True)):
            x = numpy.zeros(4000, dtype=numpy.complex64)
            x[1000] = 1.0
            tb = gr.top_block()
            ch = moving_target_channel_cc(targets=[(rng_m, 0.0, 0.0)], noise_power=0.0)
            sink = blocks.vector_sink_c()
            tb.connect(blocks.vector_source_c(x.tolist()), ch, sink)
            tb.run()
            y = numpy.abs(numpy.array(sink.data()))
            idx = numpy.arange(y.size)
            centre = float(numpy.sum(idx * y ** 2) / numpy.sum(y ** 2))
            expected = 1000 + 2 * rng_m / ofdm_isac.C0 * 10e6 + ch.latency
            self.assertAlmostEqual(centre, expected, delta=0.1)

    def test_008_single_target(self):
        """Range and speed of one target are found within a quarter of a resolution cell."""
        for r, v in ((300.0, 20.0), (517.3, -35.0), (120.0, 0.0)):
            rx, _ = run_chain([(r, v, -20.0)])
            self.assertEqual(len(rx.last_targets), 1, rx.last_targets)
            rr, vv, snr = rx.last_targets[0]
            self.assertAlmostEqual(rr, r, delta=CFG.range_resolution / 4)
            self.assertAlmostEqual(vv, v, delta=CFG.speed_resolution / 4)
            self.assertGreater(snr, 15.0)

    def test_009_two_targets_resolved(self):
        """Two targets 5 range cells apart are both reported; two inside one cell give one report."""
        rx, _ = run_chain([(300.0, 20.0, -20.0), (450.0, -60.0, -24.0)])
        found = sorted((round(t[0]), round(t[1])) for t in rx.last_targets)
        self.assertEqual(len(found), 2, rx.last_targets)
        self.assertAlmostEqual(found[0][0], 300, delta=5)
        self.assertAlmostEqual(found[1][0], 450, delta=5)
        rx, _ = run_chain([(300.0, 20.0, -20.0), (306.0, 20.0, -20.0)])
        self.assertEqual(len(rx.last_targets), 1, rx.last_targets)

    def test_010_noise_only_no_detection(self):
        """Noise alone gives no detection at pfa = 1e-7 (10 frames, about 330 000 cells)."""
        rx, _ = run_chain([], noise_power=1e-2, n_frames=10, rx_kwargs=dict(pfa=1e-7))
        # the receiver keeps only the last frame, so count over several runs of the core function
        total = 0
        rng = numpy.random.default_rng(11)
        for k in range(10):
            tx = ofdm_isac.modulate(ofdm_isac.frame_symbols(CFG, 1, k)[0], CFG)
            noise = (rng.standard_normal(tx.size) + 1j * rng.standard_normal(tx.size)) * 0.1
            targets, _ = ofdm_isac.sense_frame(noise, tx, CFG, pfa=1e-7)
            total += len(targets)
        self.assertEqual(len(rx.last_targets), 0)
        self.assertEqual(total, 0)

    def test_011_cfar_false_alarm_rate(self):
        """On independent exponential cells the CA-CFAR gives the requested pfa."""
        rng = numpy.random.default_rng(3)
        power = rng.exponential(1.0, size=(256, 128))
        for pfa in (1e-2, 1e-3):
            hits = 0
            cells = 0
            for _ in range(40):
                power = rng.exponential(1.0, size=(256, 128))
                det, _ = ofdm_isac.ca_cfar_2d(power, pfa)
                hits += int(det.sum())
                cells += det.size
            self.assertAlmostEqual(hits / cells / pfa, 1.0, delta=0.2)

    def test_012_map_output_and_message(self):
        """The block outputs one map per frame (zero speed in the middle row) and a targets message."""
        rx, data = run_chain([(240.0, 0.0, -15.0)], n_frames=2)
        self.assertEqual(data.size, 2 * rx.n_doppler * rx.n_range)
        cmap = data[-rx.n_doppler * rx.n_range:].reshape(rx.n_doppler, rx.n_range)
        i, j = numpy.unravel_index(cmap.argmax(), cmap.shape)
        self.assertEqual(i, rx.n_doppler // 2)
        self.assertAlmostEqual(rx.range_axis()[j], 240.0, delta=CFG.range_resolution / 2)
        self.assertAlmostEqual(rx.speed_axis()[i], 0.0)
        self.assertEqual(rx.messages.num_messages(), 2)  # one message per frame
        last = pmt.to_python(rx.messages.get_message(1))
        self.assertEqual(len(last), 1)
        self.assertAlmostEqual(last[0][0], 240.0, delta=CFG.range_resolution / 4)
        self.assertAlmostEqual(last[0][1], 0.0, delta=CFG.speed_resolution / 4)

    def test_013_timing_and_ber(self):
        """The same waveform works as a data link: CP timing is found and the bits are right."""
        n_frames = 3
        n = CFG.frame_len * (n_frames + 1)
        tb = gr.top_block()
        tx = ofdm_isac_tx_c(seed=9)
        ch = moving_target_channel_cc(targets=[(180.0, 15.0, -3.0)], noise_power=0.001, one_way=True)
        sink = blocks.vector_sink_c()
        tb.connect(tx, blocks.head(gr.sizeof_gr_complex, n), ch, sink)
        tb.run()
        y = numpy.array(sink.data())
        delay = 180.0 / ofdm_isac.C0 * 10e6 + ch.latency  # samples (6 + 16)
        start = ofdm_isac.estimate_timing_cp(y, CFG)
        self.assertAlmostEqual(start, delay, delta=1.0)
        errors = bits = 0
        for k in range(n_frames):
            e, b = ofdm_isac.bit_errors(y[k * CFG.frame_len:], CFG, 9, k, start)
            errors += e
            bits += b
        self.assertLess(errors / bits, 1e-3)

    def test_014_ber_noise_only_is_half(self):
        """Without a signal the bits are random: BER about 0.5 (the test cannot pass by accident)."""
        rng = numpy.random.default_rng(1)
        noise = rng.standard_normal(CFG.frame_len + 10) + 1j * rng.standard_normal(CFG.frame_len + 10)
        e, b = ofdm_isac.bit_errors(noise, CFG, 9, 0, 0)
        self.assertAlmostEqual(e / b, 0.5, delta=0.03)

    def test_015_energy_detector_sees_isac_signal(self):
        """The M4 total-power decision of the Energy Detector sees the wideband OFDM signal."""
        n = 1024 * 8 * 20
        decisions = {}
        for name, sig_gain in (("noise", -300.0), ("signal", -6.0)):
            tb = gr.top_block()
            tx = ofdm_isac_tx_c(seed=2)
            ch = moving_target_channel_cc(targets=[(0.0, 0.0, sig_gain)], noise_power=1.0, max_range=100.0)
            det = energy_detector_cf(1024, 8, 6.0, False, 1e-3, 1.0, 1, 0.01)
            d = blocks.vector_sink_b()
            tb.connect(tx, blocks.head(gr.sizeof_gr_complex, n + 200), ch, det)
            tb.connect((det, 0), blocks.null_sink(gr.sizeof_float * 1024))
            tb.connect((det, 1), blocks.null_sink(gr.sizeof_char * 1024))
            tb.connect((det, 2), d)
            tb.run()
            decisions[name] = numpy.mean(numpy.array(d.data()))
        self.assertLess(decisions["noise"], 0.2)
        self.assertGreater(decisions["signal"], 0.9)

    def test_016_channel_function_equals_block(self):
        """The plain-Python channel function (used by the experiments) equals the block."""
        tx = ofdm_isac_tx_c(seed=4)
        x = numpy.concatenate((tx.frame(0), tx.frame(1)))
        targets = [(157.5, -42.0, -8.0), (301.0, 15.0, -12.0)]
        tb = gr.top_block()
        ch = moving_target_channel_cc(targets=targets, noise_power=0.0)
        sink = blocks.vector_sink_c()
        tb.connect(blocks.vector_source_c(x.tolist()), ch, sink)
        tb.run()
        ref = ofdm_isac.simulate_channel(x, targets, CFG)
        self.assertLess(float(numpy.abs(numpy.array(sink.data()) - ref).max()), 1e-4)


if __name__ == "__main__":
    gr_unittest.run(qa_ofdm_isac)
