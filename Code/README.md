# Code

`src/gr-isac/` is the GNU Radio 3.10 out-of-tree module (Python, created with `gr_modtool`).
It holds the Energy Detector block, its GRC description, an example flowgraph, 19 QA tests and the experiment scripts.
M4 added three optional features to the block (CFAR threshold, total-power decision, `bands` message output), 10 tests (19 in total), experiments, a style configuration and a changelog. Without the new options the block behaves as in M3.

Quick start (Ubuntu 22.04 or 24.04 with GNU Radio 3.10):

```
sudo apt install gnuradio gnuradio-dev cmake build-essential pybind11-dev
cd Code/src/gr-isac
mkdir build && cd build
cmake .. && make && make test
sudo make install && sudo ldconfig
cd .. && bash scripts/run_m4_all.sh ~/m4_run   # tests, demos, flake8, all experiments; 4 to 6 minutes
python3 scripts/demo_new_features.py          # quick demo of the M4 features
```

See `src/gr-isac/README.md` for the block's parameters and outputs.
