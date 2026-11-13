# gr-isac

A GNU Radio 3.10 out-of-tree (OOT) module by Group 05, ECE 310.
License: GPL-3.0-or-later. Written in Python, created with `gr_modtool`.

## What is in it

| Block | Status |
|-------|--------|
| `isac.energy_detector_cf` (GRC name: **Energy Detector**) | Implemented, tested |
| OFDM ISAC transmitter, sensing receiver, range-Doppler processor, target detector | Not started (see the GSoC 2027 proposal) |

### Energy Detector

Takes complex baseband samples and, for every output item, uses
`fft_size * avg_len` input samples:

1. split into `avg_len` frames of `fft_size` samples,
2. Hann window and FFT of each frame,
3. average `|FFT|^2` over the frames, normalised so white noise of power s² reads s² in every bin,
4. compare every bin with a threshold.

| Parameter | Meaning |
|-----------|---------|
| `fft_size` | number of bins, at least 2 |
| `avg_len` | frames averaged per output, at least 1. More averaging detects weaker signals but reacts more slowly |
| `threshold_db` | threshold in dB |
| `relative` | `True`: threshold above the estimated noise floor (median bin power; needs the signal to fill less than half the band). `False`: absolute level, where power 1.0 is 0 dB |
| `pfa` (new in M4) | `0` (default): use `threshold_db`. Larger than 0: the threshold is calculated from theory for this false-alarm probability per bin (CFAR), so no dB value has to be guessed |
| `noise_power` (new in M4) | Noise power (linear) that you know. Used by the absolute CFAR threshold and by the total-power decision. Default 1.0 |
| `decision_mode` (new in M4) | `0` (default): decision = any bin occupied. `1`: total power over all bins above its threshold (sees wideband signals). `2`: either |
| `total_pfa` (new in M4) | False-alarm probability of the total-power decision. Default 0.01 |

| Output | Type | Meaning |
|--------|------|---------|
| 0 `psd_db` | float vector, `fft_size` | power spectrum in dB, fft-shifted (0 Hz in the middle) |
| 1 `occupied` | byte vector, `fft_size` | 1 = bin above threshold |
| 2 `decision` | byte | the band decision (see `decision_mode`; default: at least one bin occupied) |
| `bands` (new in M4) | message, optional | for every output, a list of `(low, high)` frequencies in cycles per sample, one pair per run of occupied bins |

All three stream outputs must be connected (GNU Radio refuses to start otherwise); use a Null Sink for any you do not need.
A Hann-windowed tone can light up the bin next to it, so a single tone may flag three bins.

Which mode to use (measured in `../../../Results/m4_experiments/`): the per-bin test is best for narrow signals and, in
relative mode, does not care if the noise power changes. The total-power decision is about 6 dB more sensitive to a signal
that fills the whole band, but it only works if `noise_power` is correct to within a few tenths of a dB (0.25 dB too much
noise gave 25 % false alarms in our test). Use `decision_mode = 1 or 2` only when the noise power is known well.

## Build and install

```
sudo apt install gnuradio gnuradio-dev cmake build-essential pybind11-dev   # Ubuntu 24.04
cd gr-isac
mkdir build && cd build
cmake ..
make
make test          # runs python/isac/qa_energy_detector_cf.py
sudo make install
sudo ldconfig
```

Then open GNU Radio Companion, use **Reload Blocks**, and the block appears under **[isac]**.
`examples/energy_detector_demo.grc` is a ready flowgraph (a tone in noise, results written to files in /tmp).
It can run without a GUI: `grcc examples/energy_detector_demo.grc` then run the Python file it writes.

## Tests

`python/isac/qa_energy_detector_cf.py` has 19 tests (9 from M3, 10 added in M4). They run real flowgraphs and check output shapes,
the noise-level normalisation, tone detection in the right bin, no false alarms on noise with a high threshold,
absolute and relative thresholds, the noise-floor estimator, and a false-alarm rate against theory.
The M4 tests compare the block with an independent scipy implementation (agreement to 4e-7 dB),
check that the relative threshold ignores a change of input gain, and check that the decision output is 1 exactly when some bin is occupied.
The other new tests cover the CFAR threshold (values against theory, the false-alarm rate for any noise level, tone detection), the total-power decision on a wideband signal, decision mode "either", the `bands` message and bad arguments.

## See it run in the terminal

After installing the module, `python3 scripts/demo_terminal.py` generates two test signals (noise only, and noise plus one tone),
runs them through the block and prints the decisions. To see each test by name, run
`python3 python/isac/qa_energy_detector_cf.py -v`. Sample outputs are saved in `../../../Results/`.

## Experiment

`scripts/pd_pfa_experiment.py OUTPUT_DIR` simulates a tone in white noise and compares measured detection and
false-alarm probabilities with closed-form theory. Results used in the M3 report are in `M3/Results/` of the group repository; the M4 re-run is in `../../../Results/`.

## M4 experiments

`scripts/m4_experiments.py OUTPUT_DIR` runs four experiments (use `--only noise`, `--only base`, `--only new` or `--only bench` for one):
noise-level error (fixed against relative threshold), comparison with a time-domain energy detector and one long FFT,
the effect of the new CFAR and total-power features, and the speed of the block in a flowgraph. `bash scripts/run_m4_all.sh` runs everything (environment log, tests,
demo, flake8, M3 and M4 experiments) and takes 3 to 5 minutes. Set `PYTHON=/usr/bin/python3` if your GNU Radio uses the system Python.
Results and figures are in `../../../Results/`.

## Code style

`flake8` (configuration in `.flake8`) reports no problems.

## Layout

```
python/isac/energy_detector_cf.py      block code
python/isac/qa_energy_detector_cf.py   tests
grc/isac_energy_detector_cf.block.yml  GRC block description
examples/energy_detector_demo.grc      example flowgraph
scripts/pd_pfa_experiment.py           Monte Carlo experiment (M3)
scripts/m4_experiments.py              noise-level, baseline and speed experiments (M4)
scripts/run_m4_all.sh                  runs tests, demo, flake8 and all experiments
scripts/demo_terminal.py               prints the block's decisions in the terminal
```

## Limits

- Python block, simulated data only. Not tried with an SDR.
- The relative mode needs a noise-only part of the band; a signal over more than half the bins biases the floor.
- Only complex float input, one band, no GUI sink of our own (use the QT GUI sinks on `psd_db`).
- Not yet submitted to GNU Radio; the module lives in our group repository.
