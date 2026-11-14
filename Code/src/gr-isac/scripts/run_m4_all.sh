#!/usr/bin/env bash
# Runs every check of the gr-isac module and saves the logs in one folder.
# Usage (module must be installed, see README.md):  bash scripts/run_m4_all.sh OUTPUT_DIR
# Use PYTHON=/usr/bin/python3.12 bash scripts/run_m4_all.sh OUT   if python3 is not the GNU Radio python.
# The style check needs flake8 (pip install flake8); set FLAKE_PY if it lives in another python.
# Takes about 5 to 8 minutes.
set -u
PY=${PYTHON:-python3}
OUT=${1:-m4_run}
HERE=$(cd "$(dirname "$0")/.." && pwd)
mkdir -p "$OUT"

{
  echo "date (UTC): $(date -u)"
  uname -a
  $PY --version
  $PY -c "from gnuradio import gr; print('GNU Radio', gr.version())"
  $PY -c "import numpy, scipy; print('numpy', numpy.__version__, 'scipy', scipy.__version__)"
} > "$OUT/environment.txt" 2>&1

echo "1/9 QA tests, Energy Detector (19 tests)"; (cd "$HERE/python/isac" && $PY qa_energy_detector_cf.py -v) > "$OUT/qa_tests_log.txt" 2>&1
echo "2/9 QA tests, OFDM ISAC chain (16 tests)"; (cd "$HERE/python/isac" && $PY qa_ofdm_isac.py -v) > "$OUT/qa_isac_tests_log.txt" 2>&1
echo "3/9 terminal demo";       $PY "$HERE/scripts/demo_terminal.py" > "$OUT/demo_terminal_output.txt" 2>&1
echo "4/9 new-features demo";   $PY "$HERE/scripts/demo_new_features.py" > "$OUT/demo_new_features_output.txt" 2>&1
echo "5/9 ISAC demo";           $PY "$HERE/scripts/demo_isac.py" "$OUT/demo_isac_map.png" > "$OUT/demo_isac_output.txt" 2>&1
echo "6/9 style check";         (${FLAKE_PY:-$PY} -m flake8 --version && cd "$HERE" && ${FLAKE_PY:-$PY} -m flake8 python/isac/*.py scripts/*.py && echo "flake8: no problems") > "$OUT/flake8_log.txt" 2>&1
echo "7/9 M3 experiment";       $PY "$HERE/scripts/pd_pfa_experiment.py" "$OUT/m3_experiment" > "$OUT/m3_experiment_log.txt" 2>&1
echo "8/9 M4 experiments (Energy Detector)"; $PY "$HERE/scripts/m4_experiments.py" "$OUT/m4_experiments" > "$OUT/m4_experiments_log.txt" 2>&1
echo "9/9 M4 experiments (OFDM ISAC chain)"; $PY "$HERE/scripts/isac_experiments.py" "$OUT/isac_experiments" > "$OUT/isac_experiments_log.txt" 2>&1
echo "Done. Logs are in $OUT"
