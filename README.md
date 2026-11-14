# M4 — GNU Radio Spectrum Sensing and an OFDM ISAC Chain (Group 5)

**Course:** ECE 310 — Wireless Communications
**Milestone:** M4 — Final Consolidation, Documentation, Maintainer Feedback and Refined GSoC 2027 Proposal
**Category:** GSoC (Google Summer of Code) — AI Spectrum in GNU Radio
**Team:** Group 5 — Kensi Patel (AU2440102), Maitri Lad (AU2440146), Rishika Bothra (AU2440179), Nikkesh Parekh (AU2440214), Aryan Thakkar (AU2440248)

## Project overview

We build towards spectrum sensing and an integrated sensing and communication (ISAC) framework in GNU Radio.
The code is `gr-isac`, an out-of-tree module. M3 had one block, an Energy Detector. In M4 we tested it harder and improved it,
built a **simulated OFDM ISAC chain** (transmitter, moving-target channel, sensing receiver with range-Doppler map; one signal carries data and works as a radar),
measured where everything works and where it does not, worked on five upstream issues, listed all our pull requests with their real status,
and refined the GSoC 2027 proposal. Everything is simulation; nothing was tested on a radio.

## What is new since M3

**1. The block was improved** (all new options are optional; the defaults give exactly the M3 behaviour):

| Weakness found by our experiments | Change in the block | Measured result |
|---|---|---|
| A threshold in dB has to be guessed, and a fixed one fails when the noise level changes (1 dB more noise gives 13 times too many false alarms) | **CFAR threshold:** parameter `pfa` calculates the threshold from theory | Measured false-alarm rate follows the target for any noise power |
| Per-bin detection needs 8 dB more SNR than a time-domain detector for a signal that fills the band | **Total-power decision:** `decision_mode` 1 or 2 | 6.2 dB less SNR needed than the M3 block (-8.98 against -2.74 dB), 1.8 dB from the time-domain detector, 0.2 dB cost for a narrow tone. Limit: it needs the true noise power (0.25 dB too much noise gives 24 % false alarms) |
| No list of what is occupied | **Message output `bands`:** list of (low, high) frequencies | Tested: two tones give two bands at the right frequencies |

**2. A simulated OFDM ISAC chain (new in the second part of M4)** — three GNU Radio blocks (`ofdm_isac_tx_c`, `moving_target_channel_cc`, `ofdm_sensing_rx_ccf`) and a library, with GRC definitions, an example flowgraph (`examples/isac_ofdm_demo.grc`) and a terminal demo (`scripts/demo_isac.py`). Measured against theory and against two correlation-radar baselines (`Results/isac_experiments/`):

| Check | Result |
|---|---|
| Resolution from the formulas | 19.99 m in range, 6.31 m/s in speed (default settings); usable range 720 m |
| One target | within a quarter of a resolution cell (rms error 0.9 m and 0.27 m/s at -20 dB echo SNR) |
| Two targets | reliably separated only from about 2 cells apart with the Hann window (38 % at 1 cell) |
| Against correlation radars (Pd = 0.9) | still target: ours -26.4 dB, coherent correlation radar -32.3 dB (6 dB better than ours), non-coherent -26.5 dB. Target at 100 m/s: ours -26.5 dB, non-coherent -26.6 dB, coherent fails. Our chain adds the speed; it is not more sensitive for a still target |
| CFAR false-alarm rate | as requested without window and padding; 1.3 to 3 times the request with the default Hann window and padding |
| Data on the same signal | BER equals QPSK theory with perfect channel knowledge; about 1.5 dB worse with the channel estimated from pilots |
| Fast targets | map dynamic range falls from 63 dB (0 m/s) to 45 dB (350 m/s) |
| Speed | the Python flowgraph takes 2 to 3.5 times the duration of the signal at 10 MS/s |

**3. More testing:** 35 QA tests (Energy Detector 19, M3: 9; ISAC chain 16). One test compares the Energy Detector with an independent scipy implementation (agreement 4e-7 dB). `flake8` reports no problems.

**4. New Energy Detector experiments** (`Results/m4_experiments/`): noise-level error; comparison with a time-domain detector and one long FFT; effect of the new features; speed of the Python block (3.9 to 32.2 million samples/s on a 2-CPU machine). The M3 results were reproduced exactly.

**5. Upstream work:** five issues worked on in M4 (table below), eight pull requests in gnuradio/gnuradio over the project (M3: four), recorded with the real status and the maintainer feedback. **None is merged.**

**6. Documentation:** this README, `Code/README.md`, the module README and changelog, the report, and a refined GSoC 2027 proposal (v0.2, `Proposal/`).

## Issues worked on in this milestone

| Issue | What it is about (our summary) | Pull request | Outcome so far | Explains in the video |
|---|---|---|---|---|
| #8082 | filter: mmse interpolator/differentiator Python bindings do not accept arrays | #8238 | Closed 8 Nov (template sections missing), no code review | Maitri Lad |
| #8239 | ControlPort Thrift back-end does not build on main | #8244 | Open, labelled, no review yet | Kensi Patel |
| #6441 | blocktool fails confusingly when pygccxml/castxml is missing | #8267 | Changes requested 8 Nov; rewritten 9 Nov (about 470 to 79 lines); waiting for re-review | Nikkesh Parekh |
| #4554 | gr-wavelet blocks have no GRC definitions | #8268 | Closed 8 Nov (template sections missing); review also asked for fixes | Rishika Bothra |
| #8162 | GRC: Ctrl+A in variable/import text fields selects the whole flowgraph | none yet | In progress; not counted | — |

Earlier milestones: #7697 (PR #8224), #3110 (#8230), #8173 (#8232), #5672 (#8234); see the table below.

## Upstream pull requests (status as of 13 November 2026)

| PR | Issue | Area | Status |
|----|-------|------|--------|
| [#8234](https://github.com/gnuradio/gnuradio/pull/8234) | #5672 | GRC: canvas grid | Approved by a maintainer, not merged (status of 1 Nov) |
| [#8224](https://github.com/gnuradio/gnuradio/pull/8224) | #7697 | gr-network: QA tests for `tcp_sink` | Open, no review (1 Nov) |
| [#8230](https://github.com/gnuradio/gnuradio/pull/8230) | #3110 | gr-qtgui: remove unused Vector Sink markers | Open, no response; may overlap #8075 (1 Nov) |
| [#8232](https://github.com/gnuradio/gnuradio/pull/8232) | #8173 | GRC: guard against missing `bokeh_gui` option | Duplicate of #8202 (same fix, opened earlier by another contributor); no further work planned |
| [#8238](https://github.com/gnuradio/gnuradio/pull/8238) | #8082 | filter: mmse bindings accept arrays | Closed 8 Nov by a maintainer (template sections missing) |
| [#8244](https://github.com/gnuradio/gnuradio/pull/8244) | #8239 | cmake/runtime: Thrift build fixes | Open, labelled, no review yet |
| [#8267](https://github.com/gnuradio/gnuradio/pull/8267) | #6441 | blocktool: clear error when pygccxml/castxml is missing | Changes requested 8 Nov; rewrite posted 9 Nov; waiting for re-review |
| [#8268](https://github.com/gnuradio/gnuradio/pull/8268) | #4554 | wavelet: GRC blocks and example | Closed 8 Nov (template sections missing) |

In progress, not counted: issue #8162 (Ctrl+A in GRC text fields selects the whole flowgraph); the pull request is not opened yet.

Who explains which pull request in the video is listed in the report (section 10.1). That is a presentation assignment, not authorship; GitHub shows the real author of each pull request.

## Folder structure

```
M4_G5/
├── README.md                this file
├── Report/                  M4_G5_Report.pdf (and .docx)
├── Video/                   M4_G5_video_link.txt (link to the 20-minute video), video script (.docx, .pdf)
├── Proposal/                refined GSoC 2027 proposal v0.2 (.docx, .pdf)
├── Code/
│   ├── README.md            build and run steps
│   └── src/gr-isac/         the module: blocks, tests, GRC files, examples, scripts, README, CHANGELOG
├── Data/README.md           no external data; all data is simulated by the scripts
└── Results/                 logs, CSV files and figures, with README.md
```

## How to reproduce

See `Code/README.md`. In short: install GNU Radio 3.10, build and install `Code/src/gr-isac`, then run `bash scripts/run_m4_all.sh OUTPUT_DIR` in that folder (9 steps, about 5 to 8 minutes).

## Limitations

- The OFDM ISAC chain is a **simulation**: one antenna, point targets, no clutter, no self-interference (the receiver gets an exact copy of the transmitted signal), one-tap equaliser without coding, Python only. Usable range 720 m. Its false-alarm rate is 1.3 to 3 times the requested value, resolution in practice is about 2 cells, and fast targets raise the map floor. For a still target a coherent correlation radar is about 6 dB more sensitive than our chain.
- All results use simulated data. Nothing was tested on a real radio.
- The M4 experiments were run on one cloud machine (Ubuntu 24.04, GNU Radio 3.10.9.2); the complete M4 run was also repeated by a team member on Ubuntu 22.04 (GNU Radio 3.10.1.1, Python 3.10.12) from a fresh clone of this repository: all 19 Energy Detector tests passed and all Pfa and Pd numbers were identical (`Results/nikkesh_ubuntu2204_run/`). Only the speeds differ.
- The speed numbers depend on the machine.
- The total-power decision needs the noise power to be known to within a few tenths of a dB. A method to estimate it is not built yet.
- The ideas are not new: gr-inspector already contains an energy detector, and gr-ofdmradar and Braun's thesis cover OFDM radar. We did not run gr-inspector, gr-ofdmradar or gr-radar, so we make no claim to be better. See report section 2 for the literature comparison.
- No pull request is merged.

## Use of generative AI

Generative AI was used to help with code, experiments, and writing. The results come from experiments that were run, and every number in the report can be traced to a file in `Results/`. See the report (section 10.3).
