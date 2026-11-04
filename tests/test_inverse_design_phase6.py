"""Unit tests for Phase 6 Inverse Fuel Blend Design Deliverables (T6.1 - T6.6)."""
from pathlib import Path
import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parents[1]
CONFIGS_DIR = ROOT / "configs"
T6_DIR = ROOT / "results" / "t6"
DESIGNS_CSV = T6_DIR / "inverse_designs_three_families.csv"
MARGIN_CSV = T6_DIR / "flash_margin_sweep.csv"
ABLATION_CSV = T6_DIR / "inverse_design_ablations.csv"


def test_t61_yaml_configs():
    for fam in ["jet", "diesel", "gasoline"]:
        p = CONFIGS_DIR / f"specs_{fam}.yaml"
        assert p.exists(), f"Missing config: {p}"
        data = yaml.safe_load(p.read_text(encoding="utf-8"))
        assert data["fuel_family"] == fam
        assert "constraints" in data
        assert len(data["constraints"]) >= 4
        assert "uncertainty_inflation_factors" in data
        assert data["uncertainty_inflation_factors"]["viscosity_kin_cSt"] == 1.377


def test_t64_inverse_designs_feasibility():
    assert DESIGNS_CSV.exists(), f"Missing {DESIGNS_CSV}"
    df = pd.read_csv(DESIGNS_CSV)
    assert len(df) == 3
    assert set(df["fuel_family"]) == {"jet", "diesel", "gasoline"}
    # Acceptance criterion: every reported design is labeled FEASIBLE
    assert (df["feasibility_status"] == "FEASIBLE").all()
    assert (df["all_constraints_passed"] == True).all()


def test_t65_sparsity_cardinality():
    df = pd.read_csv(DESIGNS_CSV)
    # Cardinality must be sparse (<= 5 components) even though pool is 45 candidates
    for c in df["cardinality"]:
        assert 1 <= c <= 5


def test_t62_flash_margin_sweep():
    assert MARGIN_CSV.exists(), f"Missing {MARGIN_CSV}"
    df = pd.read_csv(MARGIN_CSV)
    assert len(df) >= 4
    assert 3.4 in df["flash_safety_margin_C"].values  # Principled 90% quantile


def test_t66_ablations():
    assert ABLATION_CSV.exists(), f"Missing {ABLATION_CSV}"
    df = pd.read_csv(ABLATION_CSV)
    assert len(df) == 4
    assert "Physics-Only Baseline Design" in df["ablation_experiment"].values
    assert "Learned Residual + Principled UQ (Full Engine)" in df["ablation_experiment"].values
