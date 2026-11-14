# Changelog

## M4, second part (14 November 2026): OFDM ISAC chain
New blocks (Python, simulation only):
- `ofdm_isac_tx_c`: OFDM transmitter with QPSK data and pilots (a frame is 128 symbols).
- `moving_target_channel_cc`: echo channel with moving targets (delay, Doppler, gain) and noise, also a one-way mode for communication.
- `ofdm_sensing_rx_ccf`: range-Doppler map, 2-D CA-CFAR detector, target list on a message port.
- `ofdm_isac.py`: library with the waveform, range-Doppler processing, CFAR, target extraction, CP timing, pilot equaliser, bit counting and the channel as a plain function.
- GRC descriptions for all three blocks, example `examples/isac_ofdm_demo.grc`.
Tests and tools:
- `qa_ofdm_isac.py`: 16 tests (35 with the Energy Detector tests).
- `scripts/isac_experiments.py`: example map, resolution, Pd and error against SNR with two correlation-radar baselines, CFAR calibration, bit error rate, dynamic range against speed, speed.
- `scripts/demo_isac.py`; `scripts/run_m4_all.sh` extended to 9 steps.
Found during the work (and fixed): the DC subcarrier hole produced ghost detections along the Doppler row of a target; it is now interpolated.
Known: CFAR false-alarm rate is 1.3 to 3 times the requested one with the default window; inter-carrier interference raises the map floor at high speed.

## M4 (13 November 2026)
Block changes (all optional; the defaults give exactly the M3 behaviour):
- CFAR threshold: new parameter `pfa` calculates the per-bin threshold from theory, in relative or absolute mode.
- Total-power decision: new parameters `decision_mode`, `noise_power`, `total_pfa`. Sees wideband signals (about 6 dB better than the M3 block for a band-filling signal).
- New optional message output `bands`: occupied bands as (low, high) frequencies.
- Reason: the M4 experiments showed that a hand-picked threshold is fragile and that the per-bin test is weak for wideband signals.

Tests and tools:
- QA tests: 19 in total (10 new): independent scipy reference, relative mode ignores input gain, decision equals "any bin occupied", CFAR values and false-alarm rate, total power, decision mode, `bands` message, bad arguments.
- Added `scripts/m4_experiments.py` (noise-level error, baselines, new features, speed) and `scripts/run_m4_all.sh`.
- Added `.flake8`; fixed two style warnings in `scripts/pd_pfa_experiment.py`.

## M3
- First version: `isac.energy_detector_cf` block, GRC description, example flowgraph, 9 QA tests, Pd/Pfa experiment.
