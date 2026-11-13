title: gr-isac — spectrum sensing and ISAC blocks for GNU Radio
brief: Energy Detector block for spectrum sensing with CFAR threshold, total-power decision and band messages (19 QA tests); first block of a planned ISAC (sensing + communication) module.
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
gr-isac is a GNU Radio out-of-tree module. It currently holds one block, the
Energy Detector, which averages FFT frames of a complex baseband stream and
decides for every frequency bin whether the band is occupied or free (optionally with a
theory-based CFAR threshold, a total-power decision and a message output listing the occupied bands). It is
the first step of a planned ISAC (integrated sensing and communication)
framework. See README.md for usage, tests and measured results.
