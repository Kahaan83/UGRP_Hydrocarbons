"""test_v2_f_freeze.py - Verification suite for Phase F (Freeze v2 & Reproduction).

Validates:
  1. Existence and non-emptiness of reproduction scripts, Makefile, report, experiment plan, and manifest.
  2. Makefile contains 'reproduce' and 'test' targets.
  3. archive/v2_manifest.json is valid and checksums match on-disk files.
  4. FULL_REPORT_v2.md covers all mandatory sections and claims CLM-01 through CLM-11.
  5. experiment_plan.md defines the 3 blends, ASTM test standards, and 90% uncertainty margins.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def compute_sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while chunk := f.read(65536):
            h.update(chunk)
    return h.hexdigest()


def test_core_files_exist():
    required_files = [
        ROOT / "scripts" / "v2" / "reproduce_all.py",
        ROOT / "scripts" / "v2" / "run_phase_c.py",
        ROOT / "scripts" / "v2" / "run_phase_d.py",
        ROOT / "scripts" / "v2" / "run_phase_e.py",
        ROOT / "scripts" / "reproduce_all_v2.sh",
        ROOT / "Makefile",
        ROOT / "archive" / "v2_manifest.json",
        ROOT / "docs" / "FULL_REPORT_v2.md",
        ROOT / "docs" / "experiment_plan.md",
        ROOT / "docs" / "claims_register.md",
    ]
    for p in required_files:
        assert p.is_file(), f"Missing required file: {p.relative_to(ROOT)}"
        assert p.stat().st_size > 0, f"File is empty: {p.relative_to(ROOT)}"


def test_makefile_targets():
    makefile = ROOT / "Makefile"
    content = makefile.read_text(encoding="utf-8")
    assert "reproduce:" in content
    assert "test:" in content
    assert "reproduce_all.py" in content


def test_manifest_integrity():
    manifest_path = ROOT / "archive" / "v2_manifest.json"
    data = json.loads(manifest_path.read_text(encoding="utf-8"))

    assert data["tag"] == "v2-2026-11-05"
    assert "files" in data
    assert len(data["files"]) >= 50

    # Verify a random sample of critical files from the manifest
    critical_paths = [
        "results/v2/E3/designs.csv",
        "results/v2/A1/design_audit.csv",
        "results/v2/A2/recalibration_loo.csv",
        "results/v2/D1/bilinear_hindcast.csv",
        "configs/uncertainty_v2.yaml",
        "docs/FULL_REPORT_v2.md",
        "docs/experiment_plan.md",
    ]

    indexed_map = {entry["path"]: entry for entry in data["files"]}
    for cp in critical_paths:
        assert cp in indexed_map, f"Missing {cp} in manifest"
        entry = indexed_map[cp]
        disk_path = ROOT / cp
        assert disk_path.is_file()
        assert disk_path.stat().st_size == entry["bytes"]
        actual_sha = compute_sha256(disk_path)
        assert actual_sha == entry["sha256"], f"SHA mismatch for {cp}"


def test_full_report_structure_and_claims():
    report_path = ROOT / "docs" / "FULL_REPORT_v2.md"
    content = report_path.read_text(encoding="utf-8")

    # Check key sections
    assert "1. Executive Summary & The v1 Audit" in content
    assert "2. Physical & Thermodynamic Baselines" in content
    assert "3. Real-Fuel Validation Corpus" in content
    assert "4. Model Strategy Decided by Evidence" in content
    assert "5. Uncertainty Quantification & Conformal Safety Margins" in content
    assert "6. Hardened Inverse Design" in content
    assert "7. Master Claims Register Summary" in content
    assert "8. Reproduction, Codebase Architecture & Limitations" in content

    # Check all claims CLM-01 through CLM-11
    for i in range(1, 12):
        claim_tag = f"CLM-{i:02d}"
        assert claim_tag in content, f"Claim {claim_tag} missing from FULL_REPORT_v2.md"


def test_experiment_plan_structure():
    plan_path = ROOT / "docs" / "experiment_plan.md"
    content = plan_path.read_text(encoding="utf-8")

    # Check for the 3 blends
    assert "Blend 1: Jet A Minimal 3-Component Formulation" in content
    assert "Blend 2: Jet A Realistic 4-Component Formulation" in content
    assert "Blend 3: Diesel No. 2 Realistic 5-Component Formulation" in content

    # Check for required ASTM test methods
    assert "ASTM D445" in content
    assert "ASTM D4052" in content
    assert "ASTM D56" in content or "ASTM D93" in content
    assert "ASTM D5972" in content
    assert "ASTM D4809" in content
    assert "ASTM D86" in content
    assert "ASTM D6890" in content

    # Check for uncertainty and safety sections
    assert "90% Confidence Envelope" in content
    assert "Health, Safety, and Hazardous Waste Management" in content
    assert "500 mL" in content
