"""Unit tests for Phase 4 Many-Component Modeling Deliverables (T4.1 - T4.9)."""
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
T4_DIR = ROOT / "results" / "t4"


def test_t4_deliverables_exist():
    files = [
        "pair_coverage.csv",
        "truncation_study.csv",
        "reconstruction_sensitivity.csv",
        "error_vs_n_components.csv",
        "uncertainty_growth.csv",
        "physics_first_regularization.csv",
    ]
    for fn in files:
        p = T4_DIR / fn
        assert p.exists(), f"Missing deliverable: {p}"
        df = pd.read_csv(p)
        assert len(df) > 0, f"Empty deliverable: {p}"


def test_uncertainty_growth_factors():
    df = pd.read_csv(T4_DIR / "uncertainty_growth.csv")
    assert len(df) == 3
    for _, r in df.iterrows():
        assert r["fitted_growth_factor_g"] > 1.0  # Multicomponent error exceeds binary
        assert r["empirical_90pct_coverage_real_fuels"] > 0.90  # Valid coverage


def test_truncation_study_mass_monotonicity():
    df = pd.read_csv(T4_DIR / "truncation_study.csv")
    # Retained mass fraction must increase with number of retained components
    retained_mass = df["retained_mass_fraction"].tolist()
    assert retained_mass == sorted(retained_mass)
    assert retained_mass[-1] > 0.95  # >95% mass retained at top 40 components
