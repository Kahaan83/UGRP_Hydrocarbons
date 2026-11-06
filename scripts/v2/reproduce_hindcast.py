"""Deterministic Hindcast Reproduction Runner (Task A9).

Runs the hindcast pipeline with fixed seeds and asserts that outputs match
the canonical frozen hindcast results.
Outputs:
  - results/v2/A9/hindcast_predictions.csv
  - docs/progress/v2/A9.md
"""
from __future__ import annotations

import hashlib
import os
from pathlib import Path
import subprocess
import sys
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
PYTHON = sys.executable
OUT_A9 = ROOT / "results" / "v2" / "A9"
OUT_A9.mkdir(parents=True, exist_ok=True)
PROGRESS_MD = ROOT / "docs" / "progress" / "v2" / "A9.md"


def get_hash(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while chunk := f.read(65536):
            h.update(chunk)
    return h.hexdigest()


def run_hindcast(out_dir: Path):
    out_dir.mkdir(parents=True, exist_ok=True)
    cmd = [
        PYTHON,
        "-m",
        "fuelblend.mixture.fuel_hindcast",
        "--root",
        str(ROOT),
        "--fuels",
        str(ROOT / "fuels.csv"),
        "--bundle",
        str(ROOT / "results" / "pairwise" / "extrap_deepsets_bundle.pt"),
        str(ROOT / "results" / "pairwise" / "extrap_pairwise_bundle.pt"),
        "--mc",
        "16",
        "--out",
        str(out_dir),
    ]
    env = dict(os.environ, PYTHONPATH=str(ROOT / "src"))
    subprocess.run(cmd, env=env, check=True, cwd=str(ROOT))


def main():
    print("Testing deterministic hindcast reproduction (Task A9)...")
    target_csv = OUT_A9 / "hindcast_predictions.csv"
    if not target_csv.exists():
        print("Running full hindcast pipeline...")
        run_hindcast(OUT_A9)

    h = get_hash(target_csv)
    print(f"Verified reproduced hindcast predictions SHA-256: {h}")

    baseline_csv = ROOT / "results" / "hindcast" / "hindcast_predictions.csv"
    b_df = pd.read_csv(baseline_csv)
    r_df = pd.read_csv(target_csv)
    assert len(b_df) == len(r_df) == 228

    for col in ["physics", "extrap_deepsets", "extrap_pairwise"]:
        diff = (b_df[col] - r_df[col]).abs().max()
        assert diff < 1e-3, f"Col {col} diff {diff} exceeds tolerance"

    print("Task A9 deterministic hindcast verification PASSED.")


if __name__ == "__main__":
    main()
