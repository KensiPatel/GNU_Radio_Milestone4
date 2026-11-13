# M4 — GNU Radio Spectrum Sensing and ISAC (Group 5)

**Course:** ECE 310 — Wireless Communications
**Milestone:** M4 — Final Consolidation, Documentation, Maintainer Feedback and Refined GSoC 2027 Proposal
**Category:** GSoC (Google Summer of Code) — AI Spectrum in GNU Radio
**Team:** Group 5 — Kensi Patel (AU2440102), Maitri Lad (AU2440146), Rishika Bothra (AU2440179), Nikkesh Parekh (AU2440214), Aryan Thakkar (AU2440248)

## Project overview

We build towards spectrum sensing and an integrated sensing and communication (ISAC) framework in GNU Radio.
The code is `gr-isac`, an out-of-tree module with an Energy Detector block (from M3). In M4 we tested it harder,
measured where it works and where it does not, listed all our upstream pull requests with their real status,
and refined the GSoC 2027 proposal.

## What is new since M3

**1. The block was improved** (all new options are optional; the defaults give exactly the M3 behaviour):

| Weakness found by our experiments | Change in the block | Measured result |
|---|---|---|
| A threshold in dB has to be guessed, and a fixed one fails when the noise level changes (1 dB more noise gives 13 times too many false alarms) | **CFAR threshold:** parameter `pfa` calculates the threshold from theory | Measured false-alarm rate follows the target for any noise power |
| Per-bin detection needs 8 dB more SNR than a time-domain detector for a signal that fills the band | **Total-power decision:** `decision_mode` 1 or 2 | 6.2 dB less SNR needed than the M3 block (-8.98 against -2.74 dB), 1.8 dB from the time-domain detector, 0.2 dB cost for a narrow tone. Limit: it needs the true noise power (0.25 dB too much noise gives 24 % false alarms) |
| No list of what is occupied | **Message output `bands`:** list of (low, high) frequencies | Tested: two tones give two bands at the right frequencies |

**2. More testing:** 19 QA tests (M3: 9). One test compares the block with an independent scipy implementation (agreement 4e-7 dB). `flake8` reports no problems.

**3. New experiments** (`Results/`): noise-level error; comparison with a time-domain detector and one long FFT; effect of the new features; speed of the Python block (3.9 to 32.2 million samples/s on a 2-CPU machine). The M3 results were reproduced exactly.

**4. Upstream work:** eight pull requests in gnuradio/gnuradio (M3: four), recorded below with the real status and the maintainer feedback. **None is merged.**

**5. Documentation:** this README, `Code/README.md`, the module README and changelog, the report, and a refined GSoC 2027 proposal (v0.2, `Proposal/`).

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
M4/
├── README.md                this file
├── Report/                  M4_G5_Report.pdf (and .docx)
├── Video/                   M4_G5_video_link.txt (link to the 20-minute video), script (PDF)
├── Proposal/                refined GSoC 2027 proposal v0.2 (.docx, .pdf)
├── Code/
│   ├── README.md            build and run steps
│   └── src/gr-isac/         the module: block, tests, GRC file, example, scripts, README, CHANGELOG
├── Data/README.md           no external data; all data is simulated by the scripts
└── Results/                 logs, CSV files and figures, with README.md
```

## How to reproduce

See `Code/README.md`. In short: install GNU Radio 3.10, build and install `Code/src/gr-isac`, then run `bash scripts/run_m4_all.sh` in that folder.

## Limitations

- Only the Energy Detector exists (now with the M4 features). The OFDM ISAC transmitter, sensing receiver and range-Doppler processing are planned in the proposal, not built.
- All results use simulated data. Nothing was tested on a real radio.
- The M4 experiments were run on one cloud machine (Ubuntu 24.04, GNU Radio 3.10.9.2); the M3 tests were also repeated on a team member's Ubuntu 22.04 computer (GNU Radio 3.10.1.1). The M4 scripts have not been run on that computer yet.
- The speed numbers depend on the machine.
- The total-power decision needs the noise power to be known to within a few tenths of a dB. A method to estimate it is not built yet.
- The idea is not new: gr-inspector already contains an energy detector. We do not claim novelty; the aim is a tested, documented block that fits a full ISAC chain.
- No pull request is merged.

## Use of generative AI

Generative AI was used to help with code, experiments, and writing. The results come from experiments that were run, and every number in the report can be traced to a file in `Results/`. See the report (section 10.3).
