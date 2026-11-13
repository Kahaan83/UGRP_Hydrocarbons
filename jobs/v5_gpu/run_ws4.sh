#!/usr/bin/env bash
# v5 WS4 GPU job pack (delta_ws4_mixture_ml).  Run from the repository root:   bash jobs/v5_gpu/run_ws4.sh [cuda|cpu]
# Writes results/v5/BASE/delta_ws4_mixture_ml_gpu/<cand>/.  NOT run in the authoring environment (no GPU there).
set -euo pipefail
DEV="${1:-cuda}"
PY="${PYTHON:-python}"
EPOCHS="${EPOCHS:-300}"          # declared GPU reading: 300 epochs; everything else as in the CPU reading
cd "$(dirname "$0")/../.."
export PYTHONPATH=src
$PY -c "import torch,sys; print('torch', torch.__version__, 'cuda', torch.cuda.is_available()); sys.exit(0 if ('$DEV'=='cpu' or torch.cuda.is_available()) else 1)"
for f in data/cache/pure_grid.npz data/unified/measurements_dedup.csv data/unified/compounds.csv data/real_fuels/v5/boehm2022_viscosity.csv \
         data/processed/v2/thermoml_clean_hc.csv results/v5/S1/viscosity_lab_residuals_v5.csv results/v4/L/density_lab_residuals.csv \
         results/v5/BASE/delta_ml_viscosity/rows_scored.csv results/v5/BASE/delta_ml_density/rows_scored.csv; do
  [ -f "$f" ] || { echo "missing input $f (copy the gitignored data folders too)"; exit 1; }
done
OUT=results/v5/BASE/delta_ws4_mixture_ml_gpu
for c in pds gnn st; do
  mkdir -p "$OUT/$c"
  # 1. dev out-of-fold (locked systems untouched)
  $PY scripts/v5/base/delta_ws4_mixture_ml.py --stage dev    --cands "$c" --device "$DEV" --epochs "$EPOCHS" --out "$OUT/$c" 2>&1 | tee "$OUT/$c/run_log_dev.txt"
  # 2. locked, once, final 5-seed ensemble on all dev training rows
  $PY scripts/v5/base/delta_ws4_mixture_ml.py --stage locked --cands "$c" --device "$DEV" --epochs "$EPOCHS" --out "$OUT/$c" 2>&1 | tee "$OUT/$c/run_log_locked.txt"
done
$PY -m pytest tests/test_v5_ws4_set_models.py -q
echo "done: compare $OUT/<cand>/readings.json with results/v5/BASE/delta_ws4_mixture_ml/<cand>/readings.json (CPU reading)"
