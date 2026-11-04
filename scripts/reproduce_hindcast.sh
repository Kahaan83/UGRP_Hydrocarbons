#!/usr/bin/env bash
# Deterministic Hindcast Reproduction Script (Task A9)
# Runs the real-fuel hindcast with fixed seeds and compares with frozen results.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
export PYTHONPATH="${ROOT}/src:${PYTHONPATH:-}"

echo "Executing deterministic fuel hindcast..."
python -m fuelblend.mixture.fuel_hindcast \
    --root "${ROOT}" \
    --fuels "${ROOT}/fuels.csv" \
    --bundle "${ROOT}/results/pairwise/extrap_deepsets_bundle.pt" "${ROOT}/results/pairwise/extrap_pairwise_bundle.pt" \
    --mc 16 \
    --out "${ROOT}/results/v2/A9"

echo "Hindcast reproduction completed successfully."
