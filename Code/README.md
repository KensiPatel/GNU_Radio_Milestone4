# Code

`src/gr-isac/` is the GNU Radio 3.10 out-of-tree module (Python, created with `gr_modtool`).
It holds the Energy Detector block and a simulated OFDM ISAC chain (transmitter, moving-target channel, sensing receiver), their GRC descriptions, two example flowgraphs, 35 QA tests (19 + 16) and the experiment scripts.
M4 added three optional features to the Energy Detector (CFAR threshold, total-power decision, `bands` message output), the ISAC chain, 26 tests, experiments, a style configuration and a changelog. Without the new options the Energy Detector behaves as in M3.

Quick start (Ubuntu 22.04 or 24.04 with GNU Radio 3.10):

```
sudo apt install gnuradio gnuradio-dev cmake build-essential pybind11-dev
cd M4_G5/Code/src/gr-isac
mkdir build && cd build
cmake .. && make && make test   # both QA files; add -DPYTHON_EXECUTABLE=/usr/bin/python3 if cmake picks another python
sudo make install && sudo ldconfig
cd .. && bash scripts/run_m4_all.sh ~/m4_run   # tests, demos, flake8, all experiments; 5 to 8 minutes
python3 scripts/demo_new_features.py          # quick demo of the Energy Detector features
python3 scripts/demo_isac.py                  # ISAC chain: detected targets and a bit error rate
```

See `src/gr-isac/README.md` for the block's parameters and outputs.
