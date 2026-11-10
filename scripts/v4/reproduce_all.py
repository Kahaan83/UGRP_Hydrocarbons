"""reproduce_all.py - master one-command reproduction pipeline for v4 (mirrors scripts/v3/reproduce_all.py).

Runs every deterministic, cheap v4 build/eval script in dependency order, then pytest, then verifies
deliverable row counts and writes SHA256SUMS_v4.txt.

Steps that need unavailable or expensive resources are SKIP-GUARDED and reported as SKIPPED in the
summary (never silently dropped):
  * steps whose inputs live in the git-ignored nested repo ``extras/`` are skipped when it is absent;
  * ``train_set_encoder.py`` (CPU ~1600 s / GPU pack) only runs with ``--include-encoder``;
    its committed outputs in results/v4/E are verified and hashed but not regenerated;
  * ``e_posthoc_fixes.py`` is a one-off docs patch of the committed E outputs (not re-run);
  * ThermoML SLE re-extraction and any network download (fetch.py) are not part of v4 reproduction.

Usage:
    python scripts/v4/reproduce_all.py [--verify-only] [--skip-tests] [--include-encoder] [--manifest-out PATH]
"""
from __future__ import annotations

import argparse
import hashlib
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
S4 = ROOT / "scripts" / "v4"
EXTRAS = ROOT / "extras" / "fuelblend_data"

# (title, script, extra args, required paths [skip if any missing])
R2C = "data/real_fuels/v4/compositions_crc647_r2.csv"
STEPS = [
    ("P1  Task 2  identity check", "check_identities.py", [], []),
    ("P1  Task 1  provenance audit", "audit_provenance.py", [], []),
    ("P1  Task 3  quarantine + v4 tables", "quarantine.py", [], []),
    ("P1  Task 4  pure-constant audit", "audit_pure_constants.py", [], []),
    ("P2  Task 5  gasoline octane corpus", "build_gasoline_octane.py", [], [EXTRAS]),
    ("P2  Task 6  jet corpus", "build_jet_corpus.py", [], []),
    ("P2  Task 7  v4 splits", "build_splits_v4.py", [], [EXTRAS]),
    ("R2  review fix 1  realistic CRC 647 reconstruction", "rebuild_crc_reconstruction.py", [], []),
    ("F   Task 9  flash (ideal vs UNIFAC)", "eval_flash_v4.py", [], []),
    ("F   Task 19 flash input sensitivity (mesitylene / n-hexylbenzene)", "flash_input_sensitivity.py", [], []),
    ("B   Task 12 density/HoC/freezing", "eval_bulk_props_v4.py", [], [EXTRAS]),
    ("V   Task 11 viscosity nested CV", "build_viscosity_leaderboard_v4.py", [], []),
    ("R2  review fix 2  bulk props on r2 compositions", "eval_bulk_props_v4.py",
     ["--compositions", R2C, "--out", "results/v4/R2/B"], [EXTRAS]),
    ("R2  review fix 2  flash on r2 compositions", "eval_flash_v4.py",
     ["--compositions", R2C, "--out", "results/v4/R2/F"], []),
    ("R2  review fix 2  viscosity on r2 compositions", "build_viscosity_leaderboard_v4.py",
     ["--compositions", R2C, "--out", "results/v4/R2/V"], []),
    ("R2  review fix 2  saturates-split sensitivity", "recon_sensitivity.py", [], []),
    ("R2  review fix 2  v4-vs-r2 readings + README headline", "r2_readings.py", [], []),
    ("SLE Task 13 solid-liquid equilibrium", "eval_sle.py", [], []),
    ("SP  Task 14 sooting inputs", "build_sooting_inputs.py", [], [EXTRAS]),
    ("SP  Task 14 smoke point eval", "eval_smoke_point.py", [], []),
    ("O   Task 15 octane", "eval_octane_v4.py", [], []),
    ("O   Task 15 DVPE diagnostic", "eval_dvpe_v4.py", [], [EXTRAS]),
    ("L   review fix 3  known-composition lab residuals", "eval_lab_known_composition.py", [], []),
    ("U   Task 10 conformal margins", "build_uncertainty_v4.py", [], []),
    ("D2  review fix 4  designs with lab margins", "run_designs_v4.py", ["--out", "results/v4/D2"], []),
    ("D2  review fix 4  independent verification", "verify_designs_v4.py", ["--out", "results/v4/D2"], []),
]

NOT_RUN = [
    ("E   Task 16 set-encoder training (CPU ~1600 s / GPU)", "train_set_encoder.py",
     "expensive; committed outputs results/v4/E are verified + hashed. Use --include-encoder to retrain"),
    ("E   Task 16 post-hoc docs patch", "e_posthoc_fixes.py",
     "one-off patch of the committed E outputs; re-running without retraining adds nothing"),
    ("D   Task 17 designs + verification (results/v4/D)", "run_designs_v4.py / verify_designs_v4.py",
     "frozen: since review fix 4 the default engine uses the lab margins, so re-running into results/v4/D would overwrite the committed v4 "
     "evidence with D2-like content; results/v4/D is verified + hashed but not regenerated (D2 is regenerated instead)"),
    ("--  Antoine/Cox gate comparison", "compare_antoine_cox_gate.py",
     "separate uncommitted work, not a v4 plan deliverable"),
    ("--  ThermoML SLE re-extraction / downloads", "(extras thermoml_extract.py --sle, fetch.py)",
     "needs the git-ignored extras repo + 189 MB archive + network; committed thermoml_sle_hc.csv is used"),
]

# (label, path, min data rows)  -- min rows = observed count at v4 close
REQUIRED = [
    ("P1 ledger", "results/v4/P1/provenance_ledger.csv", 13546),
    ("P1 identity mismatches", "results/v4/P1/identity_mismatches.csv", 297),
    ("P1 pure-constant diffs", "results/v4/P1/pure_constant_diffs.csv", 483),
    ("P1 pure props", "data/processed/v4/pure_props_v4.csv", 411),
    ("P1 v4 inventory", "data/real_fuels/v4/inventory_v4.csv", 61),
    ("P1 v4 compositions", "data/real_fuels/v4/compositions_v4.csv", 2835),
    ("P1 quarantined compositions", "data/quarantine/compositions_quarantined.csv", 518),
    ("P2 gasoline octane", "data/real_fuels/v4/gasoline_octane_blends.csv", 210),
    ("P2 jet class compositions", "data/real_fuels/v4/jet_class_compositions.csv", 7),
    ("P2 jet properties", "data/real_fuels/v4/jet_properties_v4.csv", 2768),
    ("P2 coverage matrix", "results/v4/P2/coverage_matrix.csv", 1877),
    ("F flash lab", "results/v4/F/flash_lab_loso.csv", 126),
    ("F flash real", "results/v4/F/flash_real_fuels.csv", 58),
    ("F flash summary", "results/v4/F/flash_summary.csv", 15),
    ("F flash input sensitivity", "results/v4/F/flash_input_sensitivity.csv", 12),
    ("B density", "results/v4/B/density_v4.csv", 398),
    ("B hoc", "results/v4/B/hoc_v4.csv", 57),
    ("B freezing", "results/v4/B/freezing_v4.csv", 57),
    ("V viscosity leaderboard", "results/v4/V/viscosity_leaderboard_v4.csv", 63),
    ("V viscosity predictions", "results/v4/V/viscosity_predictions_v4.csv", 171),
    ("SLE eval", "results/v4/SLE/sle_eval.csv", 842),
    ("SP smoke point", "results/v4/SP/smoke_point_eval.csv", 17),
    ("O octane", "results/v4/O/octane_eval_v4.csv", 174),
    ("O dvpe", "results/v4/O/dvpe_eval_v4.csv", 14),
    ("U octane residuals", "results/v4/U/octane_oof_residuals.csv", 1048),
    ("U smoke residuals", "results/v4/U/smoke_point_oof_residuals.csv", 99),
    ("E encoder vs T11", "results/v4/E/viscosity_encoder_vs_task11.csv", 39),
    ("D designs", "results/v4/D/designs_v4.csv", 13),
    ("D constraints", "results/v4/D/constraints_v4.csv", 77),
    ("D robustness", "results/v4/D/robustness_v4.csv", 108),
    ("D verification", "results/v4/D/verification_v4.csv", 25),
    ("R2 reconstruction", "data/real_fuels/v4/compositions_crc647_r2.csv", 1),
    ("R2 bulk density", "results/v4/R2/B/density_v4.csv", 1),
    ("R2 flash real", "results/v4/R2/F/flash_real_fuels.csv", 1),
    ("R2 viscosity leaderboard", "results/v4/R2/V/viscosity_leaderboard_v4.csv", 1),
    ("R2 split sensitivity", "results/v4/R2/sat_split_sensitivity.csv", 1),
    ("R2 readings", "results/v4/R2/readings.csv", 1),
    ("L density lab residuals", "results/v4/L/density_lab_residuals.csv", 1),
    ("L viscosity lab residuals", "results/v4/L/viscosity_lab_residuals.csv", 1),
    ("D2 designs", "results/v4/D2/designs_v4.csv", 13),
    ("D2 constraints", "results/v4/D2/constraints_v4.csv", 1),
    ("D2 verification", "results/v4/D2/verification_v4.csv", 1),
    ("Lab splits", "splits/v4/lab_systems_splits.json", 0),
    ("Splits", "splits/v4/real_fuels_v4_splits.json", 0),
    ("Margins", "configs/uncertainty_v4.yaml", 0),
    ("Report", "docs/FULL_REPORT_v4.md", 0),
    ("Claims", "docs/claims_register.md", 0),
]

MANIFEST_FILES = [
    "docs/FULL_REPORT_v4.md", "docs/claims_register.md", "docs/model_decision.md",
    "configs/specs_astm.yaml", "configs/specs_bis.yaml", "configs/uncertainty_v4.yaml",
    "splits/v4/real_fuels_v4_splits.json", "splits/v4/lab_systems_splits.json",
]
MANIFEST_GLOBS = [
    ("results/v4", ["*/*.csv", "*/*.json", "*/*.md", "*/*/*.csv", "*/*/*.json"]),
    ("data/real_fuels/v4", ["*.csv"]),
    ("data/quarantine", ["*.csv"]),
    ("data/processed/v4", ["*.csv", "sooting/*.csv"]),
]
EXCLUDE_PARTS = {"antoine_cox_gate"}  # separate uncommitted user work


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def manifest_paths() -> list[Path]:
    out = {ROOT / p for p in MANIFEST_FILES}
    for base, pats in MANIFEST_GLOBS:
        for pat in pats:
            for p in (ROOT / base).glob(pat):
                if not (set(p.relative_to(ROOT).parts) & EXCLUDE_PARTS):
                    out.add(p)
    return sorted(out, key=lambda p: p.relative_to(ROOT).as_posix())


def write_manifest(out_path: Path) -> int:
    lines = []
    for p in manifest_paths():
        if not p.exists():
            print(f"[MANIFEST] MISSING: {p.relative_to(ROOT)}")
            continue
        lines.append(f"{sha256(p)}  {p.relative_to(ROOT).as_posix()}")
    out_path.write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")
    print(f"[MANIFEST] wrote {len(lines)} hashes to {out_path.name}")
    return len(lines)


def verify_deliverables() -> bool:
    print("\n" + "=" * 75 + "\nVERIFYING v4 DELIVERABLES\n" + "=" * 75)
    ok = True
    for label, rel, min_rows in REQUIRED:
        path = ROOT / rel
        if not path.exists():
            print(f"[FAIL] missing {label}: {rel}")
            ok = False
            continue
        if path.suffix == ".csv":
            n = sum(1 for line in open(path, encoding="utf-8") if line.strip()) - 1
            if n < min_rows:
                print(f"[FAIL] {label}: {n} rows < expected {min_rows} ({rel})")
                ok = False
            else:
                print(f"[PASS] {label:<28} {n:>6} rows  {rel}")
        else:
            print(f"[PASS] {label:<28} present        {rel}")
    print("ALL v4 DELIVERABLES VERIFIED." if ok else "SOME v4 DELIVERABLES ARE MISSING OR INCOMPLETE.")
    return ok


def main() -> None:
    ap = argparse.ArgumentParser(description="Master reproduction pipeline for v4")
    ap.add_argument("--verify-only", action="store_true", help="verify row counts + write manifest, no re-run")
    ap.add_argument("--skip-tests", action="store_true")
    ap.add_argument("--include-encoder", action="store_true", help="also retrain the Task 16 set encoder (slow)")
    ap.add_argument("--manifest-out", default=str(ROOT / "SHA256SUMS_v4.txt"))
    args = ap.parse_args()
    manifest_out = Path(args.manifest_out)

    if args.verify_only:
        ok = verify_deliverables()
        write_manifest(manifest_out)
        sys.exit(0 if ok else 1)

    print(f"v4 MASTER REPRODUCTION  cwd={ROOT}  python={sys.executable}")
    start = time.time()
    ran, skipped = [], []
    for idx, (title, script, extra, needs) in enumerate(STEPS, 1):
        missing = [str(n.relative_to(ROOT)) for n in needs if not n.exists()]
        if missing:
            print(f"\n[{idx}/{len(STEPS)}] SKIPPED {title}: missing {missing}")
            skipped.append((title, f"missing {missing}"))
            continue
        print(f"\n[{idx}/{len(STEPS)}] EXECUTE {title}")
        t0 = time.time()
        res = subprocess.run([sys.executable, str(S4 / script), *extra], cwd=str(ROOT))
        if res.returncode != 0:
            print(f"ERROR: {script} failed with exit code {res.returncode}")
            sys.exit(res.returncode)
        print(f"SUCCESS {script} in {time.time() - t0:.1f}s")
        ran.append((title, time.time() - t0))

    if args.include_encoder:
        print("\n[ENCODER] retraining set encoder (CPU)")
        res = subprocess.run([sys.executable, str(S4 / "train_set_encoder.py"), "--device", "cpu", "--seeds", "5",
                              "--epochs", "40", "--pools", "attention", "--lab-folds", "5"], cwd=str(ROOT))
        if res.returncode != 0:
            sys.exit(res.returncode)
    else:
        skipped.append((NOT_RUN[0][0], NOT_RUN[0][2]))
    for title, _script, why in NOT_RUN[1:]:
        skipped.append((title, why))

    if not args.skip_tests:
        print("\n[TESTS] pytest tests -q (a baseline of pre-existing failures exists: results/v4/baseline_pytest.txt)")
        res = subprocess.run([sys.executable, "-m", "pytest", "tests", "-q", "--tb=no", "-p", "no:cacheprovider"],
                             cwd=str(ROOT))
        print(f"[TESTS] pytest exit code {res.returncode} (nonzero is expected while the documented baseline failures remain)")

    print(f"\n[OK] {len(ran)} steps ran in {time.time() - start:.1f}s")
    print("\nSKIPPED / NOT RE-RUN (reported, not dropped):")
    for title, why in skipped:
        print(f"  - {title}: {why}")

    ok = verify_deliverables()
    write_manifest(manifest_out)
    if not ok:
        print("\nFATAL: deliverable verification failed.")
        sys.exit(1)
    print("\nMASTER v4 REPRODUCTION FINISHED.")


if __name__ == "__main__":
    main()
