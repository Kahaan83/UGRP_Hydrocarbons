"""test_v2_e_inverse_design.py - Test Suite for Phase E (Tasks E1-E6).

Verifies:
  - Task E1: Hard constraint enforcement and deliberate infeasibility triggers.
  - Task E2: Calibrated uncertainty margin configuration schema and metadata.
  - Task E3: Physical feasibility of designed blends for Jet A, Diesel No. 2, and Gasoline.
  - Task E4: Pareto front mappings and candidate pool / cardinality sweeps.
  - Task E5: Surrogate hindcasting and side-by-side comparison deliverables.
  - Task E6: Model formulation and uncertainty margin ablation deliverables.
"""
from pathlib import Path
import pytest
import pandas as pd
import yaml

from fuelblend.design.inverse_design_v2 import (
    CANDIDATE_REGISTRY,
    compute_density_15C,
    compute_viscosity_kin_cSt,
    compute_flash_point_C,
    compute_freezing_point_wat_C,
    compute_heat_of_combustion_MJ_kg,
    compute_aromatic_content,
    compute_distillation_d86,
    compute_dvpe_kpa,
    compute_dcn,
    compute_octane,
    evaluate_blend,
)

ROOT = Path(__file__).resolve().parents[1]
CONFIGS_DIR = ROOT / "configs"
UNCERTAINTY_YAML = CONFIGS_DIR / "uncertainty_v2.yaml"
E3_CSV = ROOT / "results" / "v2" / "E3" / "designs.csv"
E4_PARETO_CSV = ROOT / "results" / "v2" / "E4" / "pareto.csv"
E4_CARD_CSV = ROOT / "results" / "v2" / "E4" / "cardinality_sweep.csv"
E5_CSV = ROOT / "results" / "v2" / "E5" / "surrogate_comparison.csv"
E6_CSV = ROOT / "results" / "v2" / "E6" / "ablations.csv"


# ===========================================================================
# TASK E1 TESTS: Constraint Enforcement & Infeasibility Triggers
# ===========================================================================
def test_e1_deliberate_infeasibility_jet_hexadecane_biphenyl():
    """Verify that Jet A with 74.5% n-hexadecane + 25.5% biphenyl fails freezing point and viscosity."""
    comp = [("n-hexadecane", 0.745), ("biphenyl", 0.255)]
    res = evaluate_blend(comp, "jet")
    assert not res["all_passed"], "Deliberately infeasible Jet A blend must fail!"

    failed_props = {c["property"] for c in res["constraints"] if not c["passed"]}
    assert "freezing_point_C" in failed_props, "Must fail freezing point (WAT far above -40 C)"
    assert res["properties"]["freezing_point_C"] > 0.0, "Hexadecane+biphenyl freezes above 0 C"


def test_e1_deliberate_infeasibility_jet_dodecane_tetralin():
    """Verify that Jet A with 74.5% n-dodecane + 25.5% tetralin fails freezing point and aromatics."""
    comp = [("n-dodecane", 0.745), ("tetralin", 0.255)]
    res = evaluate_blend(comp, "jet")
    assert not res["all_passed"], "Dodecane + tetralin must fail specification!"

    failed_props = {c["property"] for c in res["constraints"] if not c["passed"]}
    assert "freezing_point_C" in failed_props, "Dodecane+tetralin must fail freezing point <= -40 C"
    assert "aromatics_wt_pct" in failed_props, "25.5 wt% tetralin must fail aromatics <= 25.0 wt%"
    assert res["properties"]["freezing_point_C"] > -32.0, "Freezing point cannot reach -40 C"


def test_e1_deliberate_infeasibility_gasoline_polycyclics():
    """Verify that gasoline composed of heavy non-volatile polycyclics fails density and volatility."""
    comp = [("hexylbenzene", 0.25), ("biphenyl", 0.25), ("1-methylnaphthalene", 0.25), ("tetralin", 0.25)]
    res = evaluate_blend(comp, "gasoline")
    assert not res["all_passed"], "Heavy polycyclics must fail gasoline specification!"

    failed_props = {c["property"] for c in res["constraints"] if not c["passed"]}
    assert "density_kg_m3" in failed_props, "Polycyclics density must exceed 775 kg/m3"
    assert "vapor_pressure_kPa" in failed_props, "Heavy polycyclics must fail minimum vapor pressure"


def test_e1_all_eleven_properties_computed():
    """Verify that all 11 target physical properties are computed with valid finite numbers."""
    comp = [("isododecane", 0.58), ("n-dodecane", 0.12), ("decalin", 0.15), ("tetralin", 0.15)]
    res = evaluate_blend(comp, "jet")
    props = res["properties"]

    assert props["density_kg_m3"] > 0.0
    assert -100.0 < props["freezing_point_C"] < 50.0
    assert 0.0 < props["flash_point_C"] < 150.0
    assert 0.0 < props["viscosity_minus20C_recal_cSt"] < 50.0
    assert props["heat_of_combustion_net_MJ_kg"] > 40.0
    assert 0.0 <= props["aromatics_wt_pct"] <= 100.0
    assert 0.0 <= props["aromatics_vol_pct"] <= 100.0
    assert props["T10_D86_C"] > 0.0
    assert props["T90_D86_C"] > props["T10_D86_C"]
    assert props["dcn"] > 0.0
    assert props["ron"] > 0.0


# ===========================================================================
# TASK E2 TESTS: Uncertainty Safety Margins Schema & Values
# ===========================================================================
def test_e2_uncertainty_yaml_exists_and_valid():
    assert UNCERTAINTY_YAML.exists(), f"Missing {UNCERTAINTY_YAML}"
    with open(UNCERTAINTY_YAML, "r", encoding="utf-8") as f:
        data = yaml.safe_load(f)

    margins = data.get("margins", {})
    required_keys = [
        "flash_point_C", "viscosity_kin_cSt", "density_kg_m3", "freezing_point_C",
        "heat_of_combustion_net_MJ_kg", "dcn", "ron", "mon",
        "distillation_T10_C", "distillation_T90_C", "vapor_pressure_kPa"
    ]
    for k in required_keys:
        assert k in margins, f"Missing margin key: {k}"
        entry = margins[k]
        assert "value" in entry
        assert "target_coverage" in entry
        assert "source_task" in entry
        assert "source_file" in entry

    # Verify exact calibrated numbers
    assert margins["flash_point_C"]["value"] == 5.56
    assert margins["viscosity_kin_cSt"]["value"] == 1.94
    assert margins["density_kg_m3"]["value"] == 2.50
    assert margins["freezing_point_C"]["value"] == 6.50
    assert margins["heat_of_combustion_net_MJ_kg"]["value"] == 0.28
    assert margins["dcn"]["value"] == 5.8
    assert margins["ron"]["value"] == 8.5
    assert margins["mon"]["value"] == 7.8
    assert margins["distillation_T10_C"]["value"] == 6.8
    assert margins["distillation_T90_C"]["value"] == 7.2
    assert margins["vapor_pressure_kPa"]["value"] == 24.5


# ===========================================================================
# TASK E3 TESTS: Re-run Inverse Designs
# ===========================================================================
def test_e3_designs_csv_exists_and_valid():
    assert E3_CSV.exists(), f"Missing {E3_CSV}"
    df = pd.read_csv(E3_CSV)
    assert len(df) >= 50

    design_ids = df["design_id"].unique()
    assert "JET_V1_CLAIMED_FAIL" in design_ids
    assert "JET_V2_MINIMAL_3COMP" in design_ids
    assert "JET_V2_REALISTIC_4COMP" in design_ids
    assert "DIESEL_V1_CSV_FAIL" in design_ids
    assert "DIESEL_V2_REALISTIC_5COMP" in design_ids
    assert "GASOLINE_V1_FAIL" in design_ids
    assert "GASOLINE_V2_REALISTIC_5COMP" in design_ids


def test_e3_jet_v2_designs_pass_all_constraints():
    df = pd.read_csv(E3_CSV)

    # 3-component minimal Jet A must pass all constraints
    sub3 = df[df["design_id"] == "JET_V2_MINIMAL_3COMP"]
    assert sub3["overall_design_passed"].all(), "JET_V2_MINIMAL_3COMP must pass all constraints!"
    assert (sub3["effective_margin"] >= 0.0).all(), "All effective margins must be non-negative!"

    # 4-component realistic Jet A must pass all constraints
    sub4 = df[df["design_id"] == "JET_V2_REALISTIC_4COMP"]
    assert sub4["overall_design_passed"].all(), "JET_V2_REALISTIC_4COMP must pass all constraints!"
    assert (sub4["effective_margin"] >= 0.0).all(), "All effective margins must be non-negative!"


def test_e3_diesel_and_gasoline_v2_pass_all_constraints():
    df = pd.read_csv(E3_CSV)

    # Realistic 5-comp diesel must pass all constraints
    sub_d = df[df["design_id"] == "DIESEL_V2_REALISTIC_5COMP"]
    assert sub_d["overall_design_passed"].all(), "DIESEL_V2_REALISTIC_5COMP must pass all constraints!"
    assert (sub_d["effective_margin"] >= 0.0).all(), "All effective margins must be non-negative!"

    # Realistic 5-comp gasoline must pass all constraints
    sub_g = df[df["design_id"] == "GASOLINE_V2_REALISTIC_5COMP"]
    assert sub_g["overall_design_passed"].all(), "GASOLINE_V2_REALISTIC_5COMP must pass all constraints!"
    assert (sub_g["effective_margin"] >= 0.0).all(), "All effective margins must be non-negative!"


# ===========================================================================
# TASK E4 TESTS: Pareto Fronts & Complexity Sweeps
# ===========================================================================
def test_e4_pareto_csv_exists_and_valid():
    assert E4_PARETO_CSV.exists(), f"Missing {E4_PARETO_CSV}"
    df = pd.read_csv(E4_PARETO_CSV)
    assert len(df) >= 30

    tradeoffs = df["tradeoff"].unique()
    assert "Freezing Point vs Flash Point" in tradeoffs
    assert "DCN vs Cloud Point" in tradeoffs
    assert "RON vs Volatility (DVPE)" in tradeoffs


def test_e4_cardinality_sweep_csv_exists_and_valid():
    assert E4_CARD_CSV.exists(), f"Missing {E4_CARD_CSV}"
    df = pd.read_csv(E4_CARD_CSV)
    assert len(df) >= 30

    assert set(df["fuel_family"].unique()) == {"jet", "diesel", "gasoline"}
    assert set(df["cardinality_limit_M"].unique()) == {2, 3, 4, 6, 10}
    assert set(df["pool_name"].unique()) == {"K_15", "K_30", "K_60"}


# ===========================================================================
# TASK E5 TESTS: Published Surrogate Comparisons
# ===========================================================================
def test_e5_surrogate_comparison_csv_exists_and_valid():
    assert E5_CSV.exists(), f"Missing {E5_CSV}"
    df = pd.read_csv(E5_CSV)
    assert len(df) >= 12

    surr_ids = df["surrogate_id"].unique()
    # Check required surrogates from spec
    assert "SURR_Edwards_2comp" in surr_ids
    assert "SURR_POSF10325_3comp" in surr_ids
    assert "SURR_Dooley_2ndGen_4comp" in surr_ids
    assert "SURR_NIST_RP1_14comp" in surr_ids
    assert "SURR_IDEAL_4comp" in surr_ids
    assert "SURR_Mueller_4comp" in surr_ids
    assert "SURR_Peugeot_4comp" in surr_ids
    assert "SURR_TRF_70" in surr_ids
    assert "SURR_TRF_80" in surr_ids
    assert "SURR_TRF_90" in surr_ids
    assert "SURR_CoOptima_4comp" in surr_ids

    # Check side-by-side v2 designs
    assert "DESIGN_JetA_v2_4comp" in surr_ids
    assert "DESIGN_Diesel_v2_5comp" in surr_ids
    assert "DESIGN_Gasoline_v2_5comp" in surr_ids

    # Edwards 2-comp must fail freezing point and aromatics
    edw = df[df["surrogate_id"] == "SURR_Edwards_2comp"].iloc[0]
    assert not edw["all_constraints_passed"]
    assert "freezing_point_C" in edw["failed_constraints"]

    # Designed Jet A, Diesel, Gasoline must pass
    jet_d = df[df["surrogate_id"] == "DESIGN_JetA_v2_4comp"].iloc[0]
    assert jet_d["all_constraints_passed"]

    gas_d = df[df["surrogate_id"] == "DESIGN_Gasoline_v2_5comp"].iloc[0]
    assert gas_d["all_constraints_passed"]


# ===========================================================================
# TASK E6 TESTS: Ablation Study
# ===========================================================================
def test_e6_ablations_csv_exists_and_valid():
    assert E6_CSV.exists(), f"Missing {E6_CSV}"
    df = pd.read_csv(E6_CSV)
    assert len(df) >= 45

    models = df["model_type"].unique()
    assert "Raw_Physics" in models
    assert "Recalibrated_Physics" in models
    assert "Recalibrated_Plus_Shrunk_Residual" in models

    regimes = df["uncertainty_regime"].unique()
    assert "M0_Zero_Margin" in regimes
    assert "M1_v1_Heuristic_Margin" in regimes
    assert "M2_v2_Conformal_Margin" in regimes
