"""reproduce_all.py - Master One-Command Reproduction Pipeline for v3.

v4 NOTICE: the v3 generators regenerate quarantined, formula-derived (fabricated) corpora and
designs. Their outputs are retracted; use scripts/v4/reproduce_all.py instead.

Mirrors scripts/v2/reproduce_all.py. Orchestrates the v3 validation-gap pipeline:
  - Phase S (data expansion + zero-leakage splits)
  - Phase R (audits R2/R6/R7, viscosity leaderboard, flash hindcast, freeze regimes, surrogates)
  - Phase G (octane, DVPE, cetane/cold-flow/density repairs)
  - Phase U (inverse designs + independent verification)
  - Test suite (pytest tests) and deliverable verification with SHA-256 manifest.

Usage:
    python scripts/v3/reproduce_all.py [--verify-only] [--skip-tests]
"""
from __future__ import annotations

import argparse
import hashlib
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SCRIPTS_V3 = ROOT / "scripts" / "v3"

PIPELINE_STEPS = [
    ("Phase S: Task S1-S6 (Real-fuel expansion + oxygenates)", SCRIPTS_V3 / "ingest_phase1_expansion.py"),
    ("Phase 0/R: Task R2 (Flash anchoring audit)", SCRIPTS_V3 / "audit_flash_anchoring.py"),
    ("Phase 0/R: Task R6 (Filter drop audit)", SCRIPTS_V3 / "audit_thermoml_filters.py"),
    ("Phase R: Task R7 (HoC propagation)", SCRIPTS_V3 / "analyze_hoc_propagation.py"),
    ("Phase R: Task R1 (Viscosity leaderboard)", SCRIPTS_V3 / "build_viscosity_leaderboard.py"),
    ("Phase R: Task R3 (Real-fuel flash hindcast)", SCRIPTS_V3 / "hindcast_real_fuels_flash.py"),
    ("Phase R: Task R5 (Freeze regimes)", SCRIPTS_V3 / "test_freeze_regimes.py"),
    ("Phase R: Task R4 (Surrogate check)", SCRIPTS_V3 / "summarize_surrogate_check.py"),
    ("Phase G: Task G2 (Octane lump model)", SCRIPTS_V3 / "run_phase2_octane.py"),
    ("Phase G: Task G6 (DVPE activity model)", SCRIPTS_V3 / "run_phase2_dvpe.py"),
    ("Phase G: Tasks G8/G10 (Cetane, cold flow, density)", SCRIPTS_V3 / "run_phase2_repairs.py"),
    ("Phase U: Tasks U1-U5 (Inverse designs)", SCRIPTS_V3 / "run_v3_inverse_designs.py"),
    ("Phase U: Task U7 (Independent design verification)", SCRIPTS_V3 / "verify_designs_independent.py"),
]

REQUIRED_DELIVERABLES = [
    ("Phase S4", ROOT / "results" / "v3" / "S4" / "coverage_matrix.csv", 20),
    ("Phase S (inventory)", ROOT / "data" / "real_fuels" / "v3" / "inventory_v3.csv", 150),
    ("Phase S (compositions)", ROOT / "data" / "real_fuels" / "v3" / "compositions_v3.csv", 3000),
    ("Phase S (splits)", ROOT / "splits" / "v3" / "real_fuels_v3_splits.json", 0),
    ("Phase R1", ROOT / "results" / "v3" / "R1" / "viscosity_leaderboard.csv", 30),
    ("Phase R2", ROOT / "results" / "v3" / "R2" / "flash_anchoring_audit.csv", 5),
    ("Phase R3", ROOT / "results" / "v3" / "R3" / "real_fuels_flash_hindcast.csv", 70),
    ("Phase R4", ROOT / "results" / "v3" / "R4" / "surrogate_check.csv", 10),
    ("Phase R5", ROOT / "results" / "v3" / "R5" / "freeze_regime_details.csv", 100),
    ("Phase R5", ROOT / "results" / "v3" / "R5" / "freeze_regime_test.csv", 3),
    ("Phase R6", ROOT / "results" / "v3" / "R6" / "filter_drop_breakdown.csv", 5),
    ("Phase R7", ROOT / "results" / "v3" / "R7" / "hoc_propagation_error.csv", 3),
    ("Phase G2", ROOT / "results" / "v3" / "G2" / "octane_lump_model_eval.csv", 3),
    ("Phase G6", ROOT / "results" / "v3" / "G6" / "dvpe_activity_eval.csv", 2),
    ("Phase G8", ROOT / "results" / "v3" / "G8" / "diesel_cetane_cold_flow_eval.csv", 10),
    ("Phase G10", ROOT / "results" / "v3" / "G10" / "thermoml_integration_eval.csv", 5),
    ("Phase U6", ROOT / "results" / "v3" / "U6" / "designs_v3.csv", 6),
    ("Phase U7", ROOT / "results" / "v3" / "U6" / "design_robustness_audit.csv", 30),
    ("Report", ROOT / "docs" / "FULL_REPORT_v3.md", 0),
    ("Claims", ROOT / "docs" / "claims_register.md", 0),
]

# Files covered by the SHA-256 manifest (docs, configs, splits, key results CSVs)
MANIFEST_GLOBS = [
    ROOT / "docs" / "FULL_REPORT_v3.md",
    ROOT / "docs" / "claims_register.md",
    ROOT / "docs" / "model_decision.md",
    ROOT / "configs" / "specs_astm.yaml",
    ROOT / "configs" / "specs_bis.yaml",
    ROOT / "configs" / "uncertainty_v3.yaml",
    ROOT / "splits" / "v3" / "real_fuels_v3_splits.json",
    *(ROOT / "results" / "v3").glob("*/*.csv"),
]


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def write_manifest() -> bool:
    manifest_path = ROOT / "SHA256SUMS.txt"
    lines = []
    for path in MANIFEST_GLOBS:
        if not path.exists():
            print(f"[MANIFEST] MISSING: {path.relative_to(ROOT)}")
            continue
        lines.append(f"{sha256(path)}  {path.relative_to(ROOT).as_posix()}")
    manifest_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"[MANIFEST] wrote {len(lines)} hashes to {manifest_path.name}")
    return len(lines) > 0


def verify_deliverables() -> bool:
    print("\n" + "=" * 75)
    print("VERIFYING ALL v3 DELIVERABLES")
    print("=" * 75)
    all_valid = True
    for phase_name, path, min_rows in REQUIRED_DELIVERABLES:
        if not path.exists():
            print(f"[FAIL] Missing {phase_name}: {path.relative_to(ROOT)}")
            all_valid = False
            continue
        try:
            if path.suffix == ".json":
                num_rows = len(__import__("json").load(open(path, encoding="utf-8")))
                print(f"[PASS] {phase_name:<22} {path.name:<32} ({num_rows} entries)")
                continue
            num_rows = sum(1 for l in open(path, encoding="utf-8") if l.strip()) - 1
            if num_rows < min_rows:
                print(f"[FAIL] {phase_name} ({path.name}): {num_rows} rows < expected min {min_rows}")
                all_valid = False
            else:
                print(f"[PASS] {phase_name:<22} {path.name:<32} ({num_rows} rows, {path.stat().st_size} B)")
        except Exception as e:
            print(f"[ERROR] Could not read {path}: {e}")
            all_valid = False
    print("=" * 75)
    print("ALL v3 DELIVERABLES VERIFIED." if all_valid else "SOME DELIVERABLES ARE MISSING OR INCOMPLETE.")
    print("=" * 75)
    return all_valid


def main() -> None:
    parser = argparse.ArgumentParser(description="Master Reproduction Pipeline for v3")
    parser.add_argument("--verify-only", action="store_true", help="Verify existing deliverables + manifest, no re-run")
    parser.add_argument("--skip-tests", action="store_true", help="Skip the pytest suite")
    args = parser.parse_args()

    if args.verify_only:
        ok = verify_deliverables()
        sys.exit(0 if ok else 1)

    print("*" * 75)
    print("STARTING UGRP HYDROCARBONS v3 MASTER REPRODUCTION PIPELINE")
    print(f"Working directory: {ROOT}\nPython: {sys.executable}")
    print("*" * 75)

    start = time.time()
    for idx, (title, script_path) in enumerate(PIPELINE_STEPS, 1):
        print(f"\n[{idx}/{len(PIPELINE_STEPS)}] EXECUTE: {title}")
        t0 = time.time()
        res = subprocess.run([sys.executable, str(script_path)], cwd=str(ROOT))
        if res.returncode != 0:
            print(f"ERROR: Step {idx} ({script_path.name}) failed with exit code {res.returncode}")
            sys.exit(res.returncode)
        print(f"SUCCESS: {script_path.name} in {time.time() - t0:.1f}s")

    if not args.skip_tests:
        print("\n[TESTS] Running full pytest suite ...")
        res = subprocess.run([sys.executable, "-m", "pytest", "tests", "-q", "--tb=short"], cwd=str(ROOT))
        if res.returncode != 0:
            print("ERROR: pytest suite failed.")
            sys.exit(res.returncode if res.returncode else 1)

    print(f"\n[OK] Pipeline steps + tests completed in {time.time() - start:.1f}s")

    if not verify_deliverables():
        print("\nFATAL: Post-execution deliverable verification failed.")
        sys.exit(1)
    write_manifest()
    print("\nMASTER v3 REPRODUCTION FINISHED.")
    sys.exit(0)


if __name__ == "__main__":
    main()
