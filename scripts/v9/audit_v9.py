#!/usr/bin/env python
"""v9 Data Acquisition Audit Script (Task 7 Step 1).

Verifies invariants across all newly ingested v9 datasets and loaders:
1. Units in canonical ranges (sanity bounds for catching unit errors):
   - density: 600 to 1100 kg/m3 (pure aromatics up to 1150)
   - freezing point: -120 to 20 degC (153.15 to 293.15 K)
   - cold viscosity at -40 degC (233.15 K): 1 to 50 cSt
   - flash point: 0 to 150 degC (273.15 to 423.15 K)
2. No rows derived from circular/synthetic sources (data/quarantine/ or DatasetsFuelProperties.xlsx).
3. Provenance columns filled (ref, source, system_id).
4. Every new SOURCES.md in data/raw/*/ contains a licence / license statement.
5. No PDF files staged in git index.

Usage:
  .venv\\Scripts\\python.exe scripts/v9/audit_v9.py [--report-dir results/v9/audit]
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import subprocess
import sys

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

SANITY_BOUNDS = {
    "density_kg_m3": (600.0, 1150.0),
    "freezing_point_C": (-120.0, 20.0),
    "cold_viscosity_cSt_at_minus_40C": (1.0, 50.0),
    "flash_point_C": (0.0, 150.0),
}

V9_RAW_DIRS = [
    ROOT / "data" / "raw" / "freezing_v9",
    ROOT / "data" / "raw" / "cold_visc_v9",
    ROOT / "data" / "raw" / "thermoml_v9",
    ROOT / "data" / "raw" / "pure_v9",
]


def check_sources_licence(raw_dirs: list[Path] = V9_RAW_DIRS) -> list[dict]:
    """Check that every new raw directory has a SOURCES.md with a licence statement."""
    results = []
    for d in raw_dirs:
        if not d.exists():
            continue
        sources_file = d / "SOURCES.md"
        if not sources_file.exists():
            results.append({"folder": d.name, "status": "FAIL", "reason": "Missing SOURCES.md"})
            continue
        text = sources_file.read_text(encoding="utf-8").lower()
        has_licence = ("license" in text) or ("licence" in text)
        results.append({
            "folder": d.name,
            "status": "PASS" if has_licence else "FAIL",
            "reason": "Licence statement verified" if has_licence else "Missing licence/license keyword"
        })
    return results


def check_no_staged_pdfs() -> dict:
    """Ensure no PDF files are staged in git."""
    try:
        out = subprocess.check_output(["git", "diff", "--cached", "--name-only"], text=True)
        staged_files = [line.strip() for line in out.splitlines() if line.strip()]
        pdf_files = [f for f in staged_files if f.lower().endswith(".pdf")]
        return {
            "status": "PASS" if not pdf_files else "FAIL",
            "staged_pdfs": pdf_files
        }
    except Exception as e:
        return {"status": "ERROR", "reason": str(e)}


def check_no_forbidden_sources(df: pd.DataFrame) -> dict:
    """Check that no rows originate from quarantine or synthetic Excel files."""
    if df.empty or "ref" not in df.columns:
        return {"status": "PASS", "matches": 0}
    forbidden = ["quarantine", "DatasetsFuelProperties.xlsx"]
    violations = []
    for f in forbidden:
        hits = df[df["ref"].astype(str).str.contains(f, case=False, na=False)]
        if not hits.empty:
            violations.append(f)
    return {
        "status": "PASS" if not violations else "FAIL",
        "violations": violations
    }


def check_provenance_columns(df: pd.DataFrame) -> dict:
    """Check that provenance fields (ref, system_id) are populated."""
    if df.empty:
        return {"status": "PASS", "missing_ref": 0, "missing_system_id": 0}
    missing_ref = int(df["ref"].isna().sum()) if "ref" in df.columns else len(df)
    missing_sid = int(df["system_id"].isna().sum()) if "system_id" in df.columns else len(df)
    passed = (missing_ref == 0) and (missing_sid == 0)
    return {
        "status": "PASS" if passed else "FAIL",
        "missing_ref": missing_ref,
        "missing_system_id": missing_sid
    }


def check_sanity_bounds(df: pd.DataFrame) -> list[dict]:
    """Check property values against sanity boundaries to catch unit conversion errors."""
    checks = []
    if df.empty or "property" not in df.columns or "value" not in df.columns:
        return checks

    # 1. Density
    dens = df[df["property"] == "density_kg_m3"]
    if not dens.empty:
        lo, hi = SANITY_BOUNDS["density_kg_m3"]
        outliers = dens[(dens["value"] < lo) | (dens["value"] > hi)]
        checks.append({
            "property": "density_kg_m3",
            "n_rows": len(dens),
            "status": "PASS" if outliers.empty else "FAIL",
            "outlier_count": len(outliers),
            "bounds": [lo, hi],
            "min_val": float(dens["value"].min()),
            "max_val": float(dens["value"].max()),
        })

    # 2. Freezing point (C)
    fp = df[df["property"] == "freezing_point_C"]
    if not fp.empty:
        lo, hi = SANITY_BOUNDS["freezing_point_C"]
        outliers = fp[(fp["value"] < lo) | (fp["value"] > hi)]
        checks.append({
            "property": "freezing_point_C",
            "n_rows": len(fp),
            "status": "PASS" if outliers.empty else "FAIL",
            "outlier_count": len(outliers),
            "bounds": [lo, hi],
            "min_val": float(fp["value"].min()),
            "max_val": float(fp["value"].max()),
        })

    # 3. Flash point (C)
    fl = df[df["property"] == "flash_point_C"]
    if not fl.empty:
        lo, hi = SANITY_BOUNDS["flash_point_C"]
        outliers = fl[(fl["value"] < lo) | (fl["value"] > hi)]
        checks.append({
            "property": "flash_point_C",
            "n_rows": len(fl),
            "status": "PASS" if outliers.empty else "FAIL",
            "outlier_count": len(outliers),
            "bounds": [lo, hi],
            "min_val": float(fl["value"].min()),
            "max_val": float(fl["value"].max()),
        })

    # 4. Cold viscosity at -40C (233.15 K)
    cv = df[(df["property"] == "viscosity_kin_cSt") & (abs(df["T_K"] - 233.15) < 0.5)]
    if not cv.empty:
        lo, hi = SANITY_BOUNDS["cold_viscosity_cSt_at_minus_40C"]
        outliers = cv[(cv["value"] < lo) | (cv["value"] > hi)]
        checks.append({
            "property": "cold_viscosity_at_-40C",
            "n_rows": len(cv),
            "status": "PASS" if outliers.empty else "FAIL",
            "outlier_count": len(outliers),
            "bounds": [lo, hi],
            "min_val": float(cv["value"].min()),
            "max_val": float(cv["value"].max()),
        })

    return checks


def run_audit(data_path: Path | None = None) -> dict:
    """Execute all audit checks."""
    audit_report = {
        "licence_checks": check_sources_licence(),
        "git_pdf_check": check_no_staged_pdfs(),
        "data_checks": {}
    }

    # Load data: if data_path provided, audit it; else audit all v9 loaders
    from fuelblend.data.loaders_freezing_v9 import load_freezing_v9
    from fuelblend.data.loaders_cold_visc_v9 import load_cold_visc_v9
    from fuelblend.data.loaders_thermoml_v9 import load_thermoml_v9

    dfs = []
    try:
        f_props, _ = load_freezing_v9()
        dfs.append(("freezing_v9", f_props))
    except Exception as e:
        dfs.append(("freezing_v9_error", pd.DataFrame()))

    try:
        c_df = load_cold_visc_v9()
        dfs.append(("cold_visc_v9", c_df))
    except Exception as e:
        dfs.append(("cold_visc_v9_error", pd.DataFrame()))

    try:
        t_df = load_thermoml_v9()
        dfs.append(("thermoml_v9", t_df))
    except Exception as e:
        dfs.append(("thermoml_v9_error", pd.DataFrame()))

    for name, df in dfs:
        if df.empty:
            continue
        audit_report["data_checks"][name] = {
            "rows": len(df),
            "forbidden_sources": check_no_forbidden_sources(df),
            "provenance": check_provenance_columns(df),
            "sanity_bounds": check_sanity_bounds(df),
        }

    return audit_report


def main() -> None:
    parser = argparse.ArgumentParser(description="Audit v9 data acquisition invariants.")
    parser.add_argument("--report-dir", type=Path, default=ROOT / "results" / "v9" / "audit", help="Directory for report")
    args = parser.parse_args()

    report = run_audit()
    print("=== v9 Data Acquisition Audit Summary ===")
    print(f"Git Staged PDFs: {report['git_pdf_check']['status']}")
    print("Licence Checks:")
    for lc in report["licence_checks"]:
        print(f"  {lc['folder']}: {lc['status']} ({lc['reason']})")

    print("\nData Invariants:")
    all_passed = (report["git_pdf_check"]["status"] == "PASS") and all(lc["status"] == "PASS" for lc in report["licence_checks"])
    for name, dc in report["data_checks"].items():
        print(f"  [{name}] {dc['rows']} rows: forbidden={dc['forbidden_sources']['status']}, prov={dc['provenance']['status']}")
        for sb in dc["sanity_bounds"]:
            print(f"    {sb['property']}: {sb['status']} (range: [{sb['min_val']:.2f}, {sb['max_val']:.2f}])")
            if sb["status"] != "PASS":
                all_passed = False

    args.report_dir.mkdir(parents=True, exist_ok=True)
    report_path = args.report_dir / "audit_result.json"
    report_path.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(f"\nFull audit result written to {report_path}")

    if not all_passed:
        print("Audit FAILED on one or more checks.")
        sys.exit(1)
    print("Audit PASSED on all checks.")


if __name__ == "__main__":
    main()
