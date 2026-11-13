#!/usr/bin/env bash
# Runs every check of the gr-isac module and saves the logs in one folder.
# Usage (module must be installed, see README.md):  bash scripts/run_m4_all.sh OUTPUT_DIR
# Use PYTHON=/usr/bin/python3.12 bash scripts/run_m4_all.sh OUT   if python3 is not the GNU Radio python.
# The style check needs flake8 (pip install flake8); set FLAKE_PY if it lives in another python.
# Takes about 4 to 6 minutes.
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

echo "1/6 QA tests";            (cd "$HERE/python/isac" && $PY qa_energy_detector_cf.py -v) > "$OUT/qa_tests_log.txt" 2>&1
echo "2/6 terminal demo";       $PY "$HERE/scripts/demo_terminal.py" > "$OUT/demo_terminal_output.txt" 2>&1
echo "3/6 new-features demo"; $PY "$HERE/scripts/demo_new_features.py" > "$OUT/demo_new_features_output.txt" 2>&1
echo "4/6 style check";         (${FLAKE_PY:-$PY} -m flake8 --version && cd "$HERE" && ${FLAKE_PY:-$PY} -m flake8 python/isac/*.py scripts/*.py && echo "flake8: no problems") > "$OUT/flake8_log.txt" 2>&1
echo "5/6 M3 experiment";       $PY "$HERE/scripts/pd_pfa_experiment.py" "$OUT/m3_experiment" > "$OUT/m3_experiment_log.txt" 2>&1
echo "6/6 M4 experiments";      $PY "$HERE/scripts/m4_experiments.py" "$OUT/m4_experiments" > "$OUT/m4_experiments_log.txt" 2>&1
echo "Done. Logs are in $OUT"
