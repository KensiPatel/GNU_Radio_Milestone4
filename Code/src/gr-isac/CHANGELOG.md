# Changelog

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
