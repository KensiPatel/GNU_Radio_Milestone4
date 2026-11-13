#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright 2026 Group 05 ECE310.
#
# SPDX-License-Identifier: GPL-3.0-or-later
#

import numpy
from gnuradio import blocks, gr, gr_unittest

try:
    from gnuradio.isac import energy_detector_cf
    from gnuradio.isac.energy_detector_cf import (cfar_threshold_ratio, detect, energy_spectrum,
                                                  noise_floor_estimate)
except ImportError:  # running from the build tree
    import os
    import sys
    sys.path.append(os.path.join(os.path.dirname(__file__), "..", "..", "python"))
    from isac import energy_detector_cf
    from isac.energy_detector_cf import cfar_threshold_ratio, detect, energy_spectrum, noise_floor_estimate

FFT = 256
AVG = 8


def complex_noise(rng, n, power):
    """Complex white Gaussian noise with the given total power."""
    s = numpy.sqrt(power / 2.0)
    return (s * (rng.standard_normal(n) + 1j * rng.standard_normal(n))).astype(numpy.complex64)


def tone(n, bin_index, amplitude, fft_size=FFT):
    """Complex exponential centred exactly on an FFT bin (bin 0 = DC)."""
    k = numpy.arange(n)
    return (amplitude * numpy.exp(2j * numpy.pi * bin_index * k / fft_size)).astype(numpy.complex64)


def run_block(data, **kw):
    """Run data through the block in a real flowgraph; return the three outputs."""
    tb = gr.top_block()
    src = blocks.vector_source_c(data.tolist(), False)
    det = energy_detector_cf(**kw)
    fft_size = kw.get("fft_size", FFT)
    s_psd = blocks.vector_sink_f(fft_size)
    s_occ = blocks.vector_sink_b(fft_size)
    s_dec = blocks.vector_sink_b(1)
    tb.connect(src, det)
    tb.connect((det, 0), s_psd)
    tb.connect((det, 1), s_occ)
    tb.connect((det, 2), s_dec)
    tb.run()
    n = len(s_dec.data())
    psd = numpy.array(s_psd.data()).reshape(n, fft_size)
    occ = numpy.array(s_occ.data()).reshape(n, fft_size)
    dec = numpy.array(s_dec.data())
    return psd, occ, dec


class qa_energy_detector_cf(gr_unittest.TestCase):

    def test_001_instance(self):
        energy_detector_cf()

    def test_002_bad_arguments(self):
        with self.assertRaises(ValueError):
            energy_detector_cf(fft_size=1)
        with self.assertRaises(ValueError):
            energy_detector_cf(avg_len=0)

    def test_003_output_shapes(self):
        rng = numpy.random.default_rng(1)
        n_out = 5
        x = complex_noise(rng, n_out * FFT * AVG, 1.0)
        psd, occ, dec = run_block(x, fft_size=FFT, avg_len=AVG, threshold_db=10.0, relative=False)
        self.assertEqual(psd.shape, (n_out, FFT))
        self.assertEqual(occ.shape, (n_out, FFT))
        self.assertEqual(dec.shape, (n_out,))

    def test_004_noise_level_is_normalised(self):
        # White noise of power 2.0 must read about 3.01 dB in every bin.
        rng = numpy.random.default_rng(2)
        x = complex_noise(rng, 40 * FFT * AVG, 2.0)
        power = energy_spectrum(x, FFT, AVG)
        self.assertAlmostEqual(10 * numpy.log10(power.mean()), 3.01, delta=0.1)

    def test_005_tone_is_detected_in_the_right_bin(self):
        # Tone at +32 bins from DC, noise power 1.0, tone power 1.0 (0 dB SNR
        # per sample; the FFT gain lifts it far above the noise).
        rng = numpy.random.default_rng(3)
        n = 10 * FFT * AVG
        x = complex_noise(rng, n, 1.0) + tone(n, 32, 1.0)
        psd, occ, dec = run_block(x, fft_size=FFT, avg_len=AVG, threshold_db=6.0, relative=True)
        centre = FFT // 2  # fftshift: DC at index FFT/2
        self.assertTrue(numpy.all(occ[:, centre + 32] == 1))
        self.assertTrue(numpy.all(dec == 1))
        self.assertEqual(int(numpy.argmax(psd.mean(axis=0))), centre + 32)

    def test_006_noise_only_is_free_with_a_high_threshold(self):
        rng = numpy.random.default_rng(4)
        x = complex_noise(rng, 10 * FFT * AVG, 1.0)
        psd, occ, dec = run_block(x, fft_size=FFT, avg_len=AVG, threshold_db=6.0, relative=True)
        self.assertEqual(int(occ.sum()), 0)
        self.assertTrue(numpy.all(dec == 0))

    def test_007_absolute_threshold(self):
        # Noise power 1.0 (0 dB). A tone far above an absolute 10 dB threshold
        # is occupied, and the same noise alone stays free.
        rng = numpy.random.default_rng(5)
        n = 4 * FFT * AVG
        noise = complex_noise(rng, n, 1.0)
        _, _, dec_noise = run_block(noise, fft_size=FFT, avg_len=AVG, threshold_db=10.0, relative=False)
        _, _, dec_tone = run_block(noise + tone(n, -20, 1.0), fft_size=FFT, avg_len=AVG,
                                   threshold_db=10.0, relative=False)
        self.assertTrue(numpy.all(dec_noise == 0))
        self.assertTrue(numpy.all(dec_tone == 1))

    def test_008_noise_floor_estimate(self):
        rng = numpy.random.default_rng(6)
        x = complex_noise(rng, 30 * FFT * AVG, 4.0)
        power = energy_spectrum(x, FFT, AVG)
        estimates = [noise_floor_estimate(row, AVG) for row in power]
        self.assertAlmostEqual(10 * numpy.log10(numpy.mean(estimates)), 6.02, delta=0.2)

    def test_009_false_alarm_rate_matches_theory(self):
        # Threshold set for a per-bin false-alarm probability of 1e-2 on
        # noise power 1.0. The measured rate over many bins must be close.
        from scipy.stats import chi2
        pfa = 1e-2
        thr = chi2.isf(pfa, 2 * AVG) / (2 * AVG)
        rng = numpy.random.default_rng(7)
        x = complex_noise(rng, 100 * FFT * AVG, 1.0)
        psd, occ, dec = run_block(x, fft_size=FFT, avg_len=AVG,
                                  threshold_db=10 * numpy.log10(thr), relative=False)
        measured = occ.mean()
        self.assertAlmostEqual(measured, pfa, delta=0.4 * pfa)

    def test_010_matches_independent_reference(self):
        # A plain loop written with scipy (not the vectorised code used by the
        # block) must give the same spectrum in dB as the block output.
        from scipy.fft import fft, fftshift
        from scipy.signal import get_window
        rng = numpy.random.default_rng(8)
        n_out = 6
        x = complex_noise(rng, n_out * FFT * AVG, 1.0) + tone(n_out * FFT * AVG, 20, 0.3)
        psd, _, _ = run_block(x, fft_size=FFT, avg_len=AVG, threshold_db=6.0, relative=True)
        w = get_window("hann", FFT, fftbins=True)
        ref = numpy.zeros((n_out, FFT))
        for i in range(n_out):
            acc = numpy.zeros(FFT)
            for j in range(AVG):
                start = (i * AVG + j) * FFT
                frame = x[start:start + FFT].astype(numpy.complex128)
                acc += numpy.abs(fft(frame * w)) ** 2 / numpy.sum(w ** 2)
            ref[i] = 10 * numpy.log10(fftshift(acc / AVG) + 1e-20)
        self.assertLess(float(numpy.max(numpy.abs(psd - ref))), 1e-3)

    def test_011_relative_mode_ignores_input_gain(self):
        # The relative threshold follows the noise floor, so making the whole
        # input 30 dB louder must not change which bins are flagged.
        rng = numpy.random.default_rng(9)
        n = 8 * FFT * AVG
        x = complex_noise(rng, n, 1.0) + tone(n, -15, 0.5)
        _, occ_quiet, _ = run_block(x, fft_size=FFT, avg_len=AVG, threshold_db=8.0, relative=True)
        _, occ_loud, _ = run_block((x * 31.6228).astype(numpy.complex64), fft_size=FFT, avg_len=AVG,
                                   threshold_db=8.0, relative=True)
        self.assertTrue(numpy.array_equal(occ_quiet, occ_loud))
        self.assertTrue(numpy.all(occ_quiet[:, FFT // 2 - 15] == 1))

    def test_012_decision_is_any_occupied_bin(self):
        rng = numpy.random.default_rng(10)
        n = 12 * FFT * AVG
        x = complex_noise(rng, n, 1.0)
        x[: 6 * FFT * AVG] += tone(6 * FFT * AVG, 40, 0.5)   # tone in the first 6 outputs only
        _, occ, dec = run_block(x, fft_size=FFT, avg_len=AVG, threshold_db=8.0, relative=True)
        self.assertTrue(numpy.array_equal(dec, occ.any(axis=1).astype(dec.dtype)))
        self.assertTrue(numpy.all(dec[:6] == 1))
        self.assertTrue(numpy.all(dec[6:] == 0))

    def test_013_cfar_ratio_matches_theory(self):
        self.assertAlmostEqual(cfar_threshold_ratio(1, 1e-3), -numpy.log(1e-3), places=6)   # exponential tail
        self.assertAlmostEqual(10 * numpy.log10(cfar_threshold_ratio(8, 1e-3)), 3.897, places=2)
        self.assertAlmostEqual(10 * numpy.log10(cfar_threshold_ratio(16, 1e-3)), 2.906, places=2)
        self.assertLess(cfar_threshold_ratio(16, 1e-3), cfar_threshold_ratio(4, 1e-3))      # more averaging, lower threshold

    def test_014_cfar_false_alarm_rate_for_any_noise_level(self):
        # The same Pfa setting must give the target false-alarm rate whatever the noise power is.
        rng = numpy.random.default_rng(14)
        n_out = 1500
        for power in (0.2, 3.0):
            x = complex_noise(rng, n_out * FFT * AVG, power)
            _, occ, _ = run_block(x, fft_size=FFT, avg_len=AVG, relative=True, pfa=0.01)
            self.assertLess(abs(occ.mean() / 0.01 - 1.0), 0.15)
            _, occ, _ = run_block(x, fft_size=FFT, avg_len=AVG, relative=False, pfa=0.01, noise_power=power)
            self.assertLess(abs(occ.mean() / 0.01 - 1.0), 0.15)

    def test_015_cfar_still_finds_the_tone(self):
        rng = numpy.random.default_rng(15)
        n = 20 * FFT * AVG
        x = complex_noise(rng, n, 1.0) + tone(n, 30, 0.4)
        _, occ, dec = run_block(x, fft_size=FFT, avg_len=AVG, relative=True, pfa=1e-4)
        self.assertTrue(numpy.all(occ[:, FFT // 2 + 30] == 1))
        self.assertTrue(numpy.all(dec == 1))

    def test_016_total_power_sees_a_wideband_signal(self):
        # 15 % more noise power in every bin: no single bin shows it, the total power does.
        rng = numpy.random.default_rng(16)
        n_out = 400
        noise_only = complex_noise(rng, n_out * FFT * AVG, 1.0)
        louder = complex_noise(rng, n_out * FFT * AVG, 1.15)
        kw = dict(fft_size=FFT, avg_len=AVG, relative=False, pfa=1e-4, noise_power=1.0,
                  decision_mode=1, total_pfa=0.01)
        _, _, dec_noise = run_block(noise_only, **kw)
        _, _, dec_loud = run_block(louder, **kw)
        self.assertLess(dec_noise.mean(), 0.04)
        self.assertGreater(dec_loud.mean(), 0.95)
        kw["decision_mode"] = 0
        _, _, dec_bins = run_block(louder, **kw)
        self.assertLess(dec_bins.mean(), 0.2)

    def test_017_decision_mode_either(self):
        rng = numpy.random.default_rng(17)
        n = 50 * FFT * AVG
        x = complex_noise(rng, n, 1.0) + tone(n, 25, 0.5)
        p = energy_spectrum(x, FFT, AVG)
        occ, d_bins = detect(p, FFT, AVG, relative=False, pfa=1e-4, decision_mode=0)
        _, d_total = detect(p, FFT, AVG, relative=False, pfa=1e-4, decision_mode=1)
        _, d_either = detect(p, FFT, AVG, relative=False, pfa=1e-4, decision_mode=2)
        self.assertTrue(numpy.array_equal(d_either, d_bins | d_total))
        self.assertTrue(numpy.all(d_either))

    def test_018_message_output_lists_the_bands(self):
        rng = numpy.random.default_rng(18)
        n_out = 4
        n = n_out * FFT * AVG
        x = complex_noise(rng, n, 1.0) + tone(n, 40, 0.6) + tone(n, -60, 0.6)
        tb = gr.top_block()
        src = blocks.vector_source_c(x.tolist(), False)
        det = energy_detector_cf(fft_size=FFT, avg_len=AVG, relative=True, pfa=1e-5)
        sinks = [blocks.null_sink(4 * FFT), blocks.null_sink(FFT), blocks.null_sink(1)]
        dbg = blocks.message_debug()
        tb.connect(src, det)
        tb.connect((det, 0), sinks[0])
        tb.connect((det, 1), sinks[1])
        tb.connect((det, 2), sinks[2])
        tb.msg_connect(det, "bands", dbg, "store")
        tb.run()
        self.assertEqual(dbg.num_messages(), n_out)
        import pmt
        msg = dbg.get_message(0)
        bands = [[pmt.to_double(pmt.vector_ref(pmt.vector_ref(msg, i), j)) for j in range(2)]
                 for i in range(pmt.length(msg))]
        self.assertEqual(len(bands), 2)                       # two separate tones, two bands
        centres = sorted((lo + hi) / 2 for lo, hi in bands)
        self.assertAlmostEqual(centres[0], -60 / FFT, delta=2.0 / FFT)
        self.assertAlmostEqual(centres[1], 40 / FFT, delta=2.0 / FFT)

    def test_019_bad_new_arguments(self):
        for kw in (dict(pfa=1.5), dict(pfa=-0.1), dict(noise_power=0.0), dict(decision_mode=3), dict(total_pfa=0.0)):
            with self.assertRaises(ValueError):
                energy_detector_cf(fft_size=FFT, avg_len=AVG, **kw)


if __name__ == '__main__':
    gr_unittest.run(qa_energy_detector_cf)
