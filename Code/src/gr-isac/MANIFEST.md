title: gr-isac — spectrum sensing and ISAC blocks for GNU Radio
brief: Energy Detector for spectrum sensing plus a simulated OFDM ISAC chain (transmitter, moving-target channel, range-Doppler sensing receiver); 35 QA tests; simulation only.
tags:
  - sdr
  - spectrum sensing
  - isac
  - ofdm
author:
  - Group 05, ECE 310 (Kensi Patel, Maitri Lad, Rishika Bothra, Nikkesh Parekh, Aryan Thakkar)
copyright_owner:
  - Group 05 ECE310
license: GPL-3.0-or-later
gr_supported_version: 3.10
repo: https://github.com/Aryanth29/g5_gsoc_spectrumsensing
---
gr-isac is a GNU Radio out-of-tree module. It holds an Energy Detector block (averaged FFT frames,
per-bin decision, optional CFAR threshold, total-power decision and a message output of occupied bands) and,
new in M4, a simulated OFDM ISAC chain: an OFDM transmitter that carries data, a channel with moving targets,
and a sensing receiver that makes a range-Doppler map and lists detected targets, using the same signal.
Everything is tested in simulation; nothing was tried on a real radio. See README.md for usage, tests and measured results.
