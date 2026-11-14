# Results (M4)

All files were produced by running the code in `../Code/src/gr-isac` with `scripts/run_m4_all.sh`
(GNU Radio 3.10.9.2, Ubuntu 24.04, Python 3.12.3, 2-CPU cloud machine, 13 November 2026).
All data is simulated (complex tones and white Gaussian noise). Nothing comes from a real radio.

| File | What it is |
|------|------------|
| `environment.txt` | Date, operating system, versions |
| `qa_tests_log.txt` | 19 Energy Detector QA tests, all OK (9 from M3 and 10 added in M4) |
| `qa_isac_tests_log.txt` | 16 QA tests of the OFDM ISAC chain, all OK |
| `flake8_log.txt` | Code-style check: no problems |
| `demo_terminal_output.txt` | M3 demo: noise only reads "free", noise plus one tone reads "OCCUPIED" |
| `demo_new_features_output.txt` | M4 demo: CFAR thresholds, a wideband signal seen by the total-power decision but not bin by bin, and the `bands` message |
| `m3_reproduced/` | The M3 Pd/Pfa experiment run again: `summary.txt` is identical to the M3 file. Also the team member's own Ubuntu 22.04 run from M3 (`team_run_*`) |
| `m4_experiments/noise_uncertainty.csv/.png` | Pfa and Pd when the real noise power differs from the design value (-2 to +3 dB): fixed against relative threshold |
| `m4_experiments/baselines.csv/.png` | M3 method against a time-domain energy detector and one long FFT, same samples and same band false-alarm rate |
| `m4_experiments/new_features.csv/.png` | Effect of the M4 features (CFAR and total power): M3 block, new block and time-domain detector for a narrowband and a wideband signal |
| `m4_experiments/new_noise_error.csv` | False-alarm rate of the per-bin tests and of the total-power decision when the real noise power is wrong |
| `m4_experiments/bench_throughput.csv` | Speed in million samples per second |
| `m4_experiments/m4_summary.txt` | Text summary of all Energy Detector experiments |
| `demo_isac_output.txt`, `demo_isac_map.png` | Terminal demo of the ISAC chain: two targets detected at the right range and speed, and a bit error rate; the range-Doppler map |
| `isac_experiments_log.txt` | Console output of the ISAC experiments |
| `isac_experiments/range_doppler_map.png`, `map_example.csv` | Example range-Doppler map with three targets, true and detected values |
| `isac_experiments/resolution.csv/.png` | Two-target resolution in range and in speed, Hann and rectangular window |
| `isac_experiments/accuracy_pd.csv`, `pd_vs_snr.png`, `accuracy_vs_snr.png` | Detection probability and rms error against echo SNR, for the OFDM chain and two correlation-radar baselines, still and moving target |
| `isac_experiments/cfar_calibration.csv` | Requested against measured false-alarm rate of the CFAR detector |
| `isac_experiments/ber.csv/.png` | Bit error rate against SNR on the same signal: pilots, pilots plus echo, perfect channel knowledge, QPSK theory |
| `isac_experiments/dynamic_range_vs_speed.csv` | Map contrast against target speed (inter-carrier interference) |
| `isac_experiments/speed.csv` | Processing time per frame and flowgraph throughput (varies by 10 to 40 % between runs) |

Reproduce: install the module (see `../Code/src/gr-isac/README.md`), then `bash scripts/run_m4_all.sh OUTPUT_DIR` from the module folder.
The ISAC results (`isac_experiments/`, `qa_isac_tests_log.txt`, `demo_isac_*`) come from the cloud machine only; the member's laptop run in `nikkesh_ubuntu2204_run/` covers the Energy Detector part unless a newer log is added there.
The random generator is seeded, so the Pfa and Pd numbers repeat exactly; the speed numbers depend on the machine and vary by about 10 to 20 percent between runs.

`nikkesh_ubuntu2204_run/` holds the terminal output of the complete M4 run on a team member's laptop (Ubuntu 22.04, GNU Radio 3.10.1.1, Python 3.10.12), from a fresh clone of this repository: `terminal_output_build_and_run.txt` (install, build, ctest, run) and `terminal_output_results.txt` (environment, test log, flake8, demos, experiment summary). The member's own `~/m4_run` files can be added to this folder.
