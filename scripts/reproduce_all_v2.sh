#!/usr/bin/env bash
# reproduce_all_v2.sh - Shell script for one-command execution of UGRP Hydrocarbons v2 pipeline.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"

echo "========================================================================"
echo "UGRP Hydrocarbons v2 - Full Pipeline Reproduction"
echo "Root directory: ${ROOT_DIR}"
echo "========================================================================"

if [ -f "${ROOT_DIR}/.venv/Scripts/python.exe" ]; then
    PYTHON_BIN="${ROOT_DIR}/.venv/Scripts/python.exe"
elif [ -f "${ROOT_DIR}/.venv/bin/python" ]; then
    PYTHON_BIN="${ROOT_DIR}/.venv/bin/python"
else
    PYTHON_BIN="python"
fi

echo "Using Python: ${PYTHON_BIN}"
cd "${ROOT_DIR}"

"${PYTHON_BIN}" scripts/v2/reproduce_all.py "$@"
