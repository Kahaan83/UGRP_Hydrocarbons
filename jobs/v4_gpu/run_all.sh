#!/usr/bin/env bash
# v4 Task 16 GPU job pack.  Run from the repository root:   bash jobs/v4_gpu/run_all.sh [cuda|cpu]
# Writes results/v4/E_gpu* and checkpoints runs/v4/gpu/ (gitignored).  Not run in the authoring environment (no GPU there).
set -euo pipefail
DEV="${1:-cuda}"
PY="${PYTHON:-python}"
cd "$(dirname "$0")/../.."
$PY -c "import torch,sys; print('torch', torch.__version__, 'cuda', torch.cuda.is_available()); sys.exit(0 if ('$DEV'=='cpu' or torch.cuda.is_available()) else 1)"
mkdir -p runs/v4/gpu results/v4/E_gpu_attention results/v4/E_gpu_sum
# 1. primary (attention pooling) incl. lab GroupKFold(system_id) evaluation; 5-seed deep ensemble
$PY scripts/v4/train_set_encoder.py --device "$DEV" --seeds 5 --epochs 300 --pools attention --lab-folds 5 \
    --out results/v4/E_gpu_attention --save-models runs/v4/gpu/attention
# 2. sum-pool ablation
$PY scripts/v4/train_set_encoder.py --device "$DEV" --seeds 5 --epochs 300 --pools sum --lab-folds 0 \
    --out results/v4/E_gpu_sum --save-models runs/v4/gpu/sum
# 3. tests
$PY -m pytest tests/test_v4_set_encoder.py -q
echo "done: compare results/v4/E_gpu_attention/promotion_decision.json with results/v4/E/promotion_decision.json (CPU run)"
