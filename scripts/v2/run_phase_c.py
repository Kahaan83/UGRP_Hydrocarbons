"""run_phase_c.py - Orchestration runner for Phase C: Property and Baseline Fixes.

Executes Tasks C1 through C8:
  - C1: Heat of Combustion Group-Additivity Baseline (eval_hoc_v2.py)
  - C2: Flash UNIFAC / Activity Coefficients & LOSO Residual (eval_flash_loso_v2.py)
  - C3: Flash Expansion Usage & Leakage Audit (verify_c3_leakage.py)
  - C4: Cold Flow SLE Baseline & Transparent Offset (eval_cold_flow_v2.py)
  - C5: Pure Cetane Compendium Digitization (digitize_cetane_compounds.py)
  - C6: Octane Blending on Real Gasolines (eval_octane_real_fuels.py)
  - C7: D86 Distillation Curve Prediction (eval_d86_prediction.py)
  - C8: Gasoline DVPE Baseline (eval_dvpe_baseline.py)
"""
from __future__ import annotations

import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SCRIPTS_V2 = ROOT / "scripts" / "v2"

PHASE_C_SCRIPTS = [
    ("Task C1 (Heat of Combustion v2)", SCRIPTS_V2 / "eval_hoc_v2.py"),
    ("Task C2 (Flash LOSO v2)", SCRIPTS_V2 / "eval_flash_loso_v2.py"),
    ("Task C3 (Flash Leakage Verification)", SCRIPTS_V2 / "verify_c3_leakage.py"),
    ("Task C4 (Cold Flow v2)", SCRIPTS_V2 / "eval_cold_flow_v2.py"),
    ("Task C5 (Cetane Compounds Digitization)", SCRIPTS_V2 / "digitize_cetane_compounds.py"),
    ("Task C6 (Octane Real Fuels)", SCRIPTS_V2 / "eval_octane_real_fuels.py"),
    ("Task C7 (D86 Distillation Prediction)", SCRIPTS_V2 / "eval_d86_prediction.py"),
    ("Task C8 (Gasoline DVPE Baseline)", SCRIPTS_V2 / "eval_dvpe_baseline.py"),
]


def run_phase_c() -> bool:
    print("=" * 70)
    print("Executing Phase C: Tasks C1 - C8 (Property & Baseline Fixes)")
    print("=" * 70)
    start_all = time.time()
    all_ok = True

    for name, script_path in PHASE_C_SCRIPTS:
        print(f"\n>>> Running {name} ({script_path.name})...")
        t0 = time.time()
        res = subprocess.run([sys.executable, str(script_path)], cwd=str(ROOT))
        elapsed = time.time() - t0
        if res.returncode != 0:
            print(f"FAILED: {name} exited with code {res.returncode}")
            all_ok = False
        else:
            print(f"PASSED: {name} completed in {elapsed:.2f}s")

    total_elapsed = time.time() - start_all
    print("=" * 70)
    if all_ok:
        print(f"Phase C completed successfully in {total_elapsed:.2f}s.")
    else:
        print(f"Phase C completed with ERRORS in {total_elapsed:.2f}s.")
    print("=" * 70)
    return all_ok


if __name__ == "__main__":
    success = run_phase_c()
    sys.exit(0 if success else 1)
