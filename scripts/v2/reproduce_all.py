"""reproduce_all.py - Master One-Command Reproduction Pipeline for v2.

Orchestrates and executes all verification, real-fuel ingestion, property baseline,
model evaluation, and hardened inverse design scripts:
  - Phase A: Audit and verification scripts (A1, A2, A3, A4, A5, A6, A7, A8, A9, A12)
  - Phase B: Real-fuel inventory, composition ingestion, distillation assembly, property matrix
  - Phase C: Property and baseline fixes (run_phase_c.py)
  - Phase D: Evidence-driven model evaluation (run_phase_d.py)
  - Phase E: Hardened inverse design (run_phase_e.py)
  - Verification: Confirms presence, integrity, and row counts of all v2 deliverables.
"""
from __future__ import annotations

import argparse
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SCRIPTS_V2 = ROOT / "scripts" / "v2"

# Sequence of Phase A & B steps
PIPELINE_STEPS = [
    # Phase A
    ("Phase A: Task A1 (Audit Inverse Designs)", SCRIPTS_V2 / "audit_inverse_designs.py"),
    ("Phase A: Task A2 (LOO Recalibration Check)", SCRIPTS_V2 / "eval_recalibration_loo.py"),
    ("Phase A: Task A3 (Cetane Basis Check)", SCRIPTS_V2 / "check_cetane_basis.py"),
    ("Phase A: Task A4 (Grunberg-Nissan Check)", SCRIPTS_V2 / "eval_gn_baselines.py"),
    ("Phase A: Task A5 (Heat of Combustion Re-evaluation)", SCRIPTS_V2 / "recheck_hoc.py"),
    ("Phase A: Task A6 (Uncertainty Growth Refit)", SCRIPTS_V2 / "refit_uncertainty_growth.py"),
    ("Phase A: Task A7 (Flash Margin Semantics)", SCRIPTS_V2 / "eval_flash_margin_semantics.py"),
    ("Phase A: Task A8 (Flash Data Audit)", SCRIPTS_V2 / "audit_flash_data.py"),
    ("Phase A: Task A9 (Deterministic Hindcast)", SCRIPTS_V2 / "reproduce_hindcast.py"),
    ("Phase A: Task A12 (Licensing & Redistribution)", SCRIPTS_V2 / "audit_data_redistribution.py"),
    # Phase B
    ("Phase B: Task B1 (Candidates Inventory)", SCRIPTS_V2 / "build_candidates_inventory.py"),
    ("Phase B: Tasks B2-B6 (Real Fuels Ingestion & Matrix)", SCRIPTS_V2 / "ingest_real_fuels_v2.py"),
    ("Phase B: Distillation Real Fuels Assembly", SCRIPTS_V2 / "assemble_distillation_real_fuels.py"),
    # Phase C
    ("Phase C: Tasks C1-C8 (Property & Baseline Fixes)", SCRIPTS_V2 / "run_phase_c.py"),
    # Phase D
    ("Phase D: Tasks D1-D6 (Evidence-Driven Model Strategy)", SCRIPTS_V2 / "run_phase_d.py"),
    # Phase E
    ("Phase E: Tasks E3-E6 (Hardened Inverse Design)", SCRIPTS_V2 / "run_phase_e.py"),
]

# Key deliverables to verify
REQUIRED_DELIVERABLES = [
    ("Phase A1", ROOT / "results" / "v2" / "A1" / "design_audit.csv", 20),
    ("Phase A2", ROOT / "results" / "v2" / "A2" / "recalibration_loo.csv", 5),
    ("Phase A3", ROOT / "results" / "v2" / "A3" / "cetane_basis_rerun.csv", 3),
    ("Phase A4", ROOT / "results" / "v2" / "A4" / "gn_baselines.csv", 5),
    ("Phase A5", ROOT / "results" / "v2" / "A5" / "hoc_recheck.csv", 3),
    ("Phase A6", ROOT / "results" / "v2" / "A6" / "uncertainty_growth_refit.csv", 4),
    ("Phase A7", ROOT / "results" / "v2" / "A7" / "flash_margin_semantics.csv", 8),
    ("Phase A8", ROOT / "results" / "v2" / "A8" / "flash_usage_audit.csv", 10),
    ("Phase A9", ROOT / "results" / "v2" / "A9" / "hindcast_summary.csv", 3),
    ("Phase A12", ROOT / "results" / "v2" / "A12" / "redistribution_audit.csv", 20),
    ("Phase B1", ROOT / "data" / "real_fuels" / "candidates_v2.csv", 50),
    ("Phase B2", ROOT / "data" / "real_fuels" / "inventory.csv", 100),
    ("Phase B2", ROOT / "data" / "real_fuels" / "compositions.csv", 2000),
    ("Phase B4", ROOT / "data" / "real_fuels" / "blend_links.csv", 15),
    ("Phase B6", ROOT / "results" / "v2" / "B6" / "real_fuels_property_matrix.csv", 3),
    ("Phase C1", ROOT / "results" / "v2" / "C1" / "hoc_eval.csv", 2),
    ("Phase C2", ROOT / "results" / "v2" / "C2" / "flash_loso_v2.csv", 6),
    ("Phase C4", ROOT / "results" / "v2" / "C4" / "cold_flow_v2.csv", 4),
    ("Phase C6", ROOT / "results" / "v2" / "C6" / "octane_real_fuel.csv", 6),
    ("Phase C7", ROOT / "results" / "v2" / "C7" / "d86_eval.csv", 6),
    ("Phase C8", ROOT / "results" / "v2" / "C8" / "dvpe_eval.csv", 2),
    ("Phase D1", ROOT / "results" / "v2" / "D1" / "bilinear_hindcast.csv", 25),
    ("Phase D2", ROOT / "results" / "v2" / "D2" / "class_interactions.csv", 15),
    ("Phase D3", ROOT / "results" / "v2" / "D3" / "shrinkage_sweep.csv", 100),
    ("Phase D4", ROOT / "results" / "v2" / "D4" / "stacked_eval.csv", 20),
    ("Phase D5", ROOT / "results" / "v2" / "D5" / "error_vs_n.csv", 15),
    ("Phase D6", ROOT / "results" / "v2" / "D6" / "truncation.csv", 15),
    ("Phase E3", ROOT / "results" / "v2" / "E3" / "designs.csv", 50),
    ("Phase E4", ROOT / "results" / "v2" / "E4" / "pareto.csv", 40),
    ("Phase E4", ROOT / "results" / "v2" / "E4" / "cardinality_sweep.csv", 40),
    ("Phase E5", ROOT / "results" / "v2" / "E5" / "surrogate_comparison.csv", 10),
    ("Phase E6", ROOT / "results" / "v2" / "E6" / "ablations.csv", 45),
]


def verify_deliverables() -> bool:
    print("\n" + "=" * 75)
    print("VERIFYING ALL v2 DELIVERABLES")
    print("=" * 75)
    all_valid = True

    for phase_name, path, min_rows in REQUIRED_DELIVERABLES:
        if not path.exists():
            print(f"[FAIL] Missing {phase_name}: {path.relative_to(ROOT)}")
            all_valid = False
            continue

        try:
            with open(path, "r", encoding="utf-8") as f:
                lines = [l for l in f if l.strip()]
            num_rows = len(lines) - 1  # subtract header
            if num_rows < min_rows:
                print(f"[FAIL] {phase_name} ({path.name}): {num_rows} rows < expected min {min_rows}")
                all_valid = False
            else:
                print(f"[PASS] {phase_name:<10} {path.name:<32} ({num_rows} data rows, {path.stat().st_size} bytes)")
        except Exception as e:
            print(f"[ERROR] Could not read {path}: {e}")
            all_valid = False

    print("=" * 75)
    if all_valid:
        print("ALL v2 DELIVERABLES VERIFIED AND PRESENT WITH EXPECTED ROW COUNTS.")
    else:
        print("SOME DELIVERABLES ARE MISSING OR INCOMPLETE.")
    print("=" * 75)
    return all_valid


def main() -> None:
    parser = argparse.ArgumentParser(description="Master Reproduction Pipeline for v2")
    parser.add_argument("--verify-only", action="store_true", help="Only verify existing deliverables without re-executing")
    args = parser.parse_args()

    if args.verify_only:
        valid = verify_deliverables()
        sys.exit(0 if valid else 1)

    print("*" * 75)
    print("STARTING UGRP HYDROCARBONS v2 MASTER REPRODUCTION PIPELINE")
    print(f"Working directory: {ROOT}")
    print(f"Python interpreter: {sys.executable}")
    print("*" * 75)

    pipeline_start = time.time()
    steps_passed = 0
    steps_total = len(PIPELINE_STEPS)

    for idx, (title, script_path) in enumerate(PIPELINE_STEPS, 1):
        print(f"\n[{idx}/{steps_total}] EXECUTE: {title}")
        step_start = time.time()
        res = subprocess.run([sys.executable, str(script_path)], cwd=str(ROOT))
        step_duration = time.time() - step_start

        if res.returncode != 0:
            print(f"ERROR: Step {idx} ({script_path.name}) failed with exit code {res.returncode}")
            sys.exit(res.returncode)
        else:
            print(f"SUCCESS: Step {idx} completed in {step_duration:.2f}s")
            steps_passed += 1

    total_time = time.time() - pipeline_start
    print("\n" + "*" * 75)
    print(f"ALL {steps_passed}/{steps_total} PIPELINE STEPS COMPLETED IN {total_time:.2f}s")
    print("*" * 75)

    # Post-execution verification
    verification_success = verify_deliverables()
    if not verification_success:
        print("\nFATAL: Post-execution deliverable verification failed.")
        sys.exit(1)

    print("\nMASTER REPRODUCTION FINISHED WITH ZERO DEFECTS.")
    sys.exit(0)


if __name__ == "__main__":
    main()
