"""run_phase_d.py - Orchestration runner for Phase D: Evidence-Driven Model Strategy.

Executes Tasks D1 through D6:
  - D1: Real-Fuel Hindcast of Bilinear Pairwise Model (eval_bilinear_hindcast.py)
  - D2: Class-Level Interaction Learning (eval_class_interactions.py)
  - D3: Shrinkage Factor Sweep on Real Fuels (eval_shrinkage_sweep.py)
  - D4: Stacked Hierarchical Model (eval_stacked_model.py)
  - D5: Error vs Component Complexity N (eval_error_vs_n_components.py)
  - D6: Truncation & Lumping Sensitivity (eval_truncation_sensitivity.py)
"""
from __future__ import annotations

import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SCRIPTS_V2 = ROOT / "scripts" / "v2"

PHASE_D_SCRIPTS = [
    ("Task D1 (Bilinear Hindcast)", SCRIPTS_V2 / "eval_bilinear_hindcast.py"),
    ("Task D2 (Class Interactions)", SCRIPTS_V2 / "eval_class_interactions.py"),
    ("Task D3 (Shrinkage Sweep)", SCRIPTS_V2 / "eval_shrinkage_sweep.py"),
    ("Task D4 (Stacked Model)", SCRIPTS_V2 / "eval_stacked_model.py"),
    ("Task D5 (Error vs N Components)", SCRIPTS_V2 / "eval_error_vs_n_components.py"),
    ("Task D6 (Truncation Sensitivity)", SCRIPTS_V2 / "eval_truncation_sensitivity.py"),
]


def run_phase_d() -> bool:
    print("=" * 70)
    print("Executing Phase D: Tasks D1 - D6 (Evidence-Driven Model Strategy)")
    print("=" * 70)
    start_all = time.time()
    all_ok = True

    for name, script_path in PHASE_D_SCRIPTS:
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
        print(f"Phase D completed successfully in {total_elapsed:.2f}s.")
    else:
        print(f"Phase D completed with ERRORS in {total_elapsed:.2f}s.")
    print("=" * 70)
    return all_ok


if __name__ == "__main__":
    success = run_phase_d()
    sys.exit(0 if success else 1)
