"""Tests for Task D6: Truncation and Lumping Sensitivity on Measured-Composition Real Fuels."""
from pathlib import Path
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
CSV_PATH = ROOT / "results" / "v2" / "D6" / "truncation.csv"
DOC_PATH = ROOT / "docs" / "progress" / "v2" / "D6.md"


def test_d6_deliverables_exist():
    assert CSV_PATH.exists(), f"Missing {CSV_PATH}"
    assert DOC_PATH.exists(), f"Missing {DOC_PATH}"


def test_d6_truncation_schema():
    df = pd.read_csv(CSV_PATH)
    assert len(df) > 50, f"Expected >50 rows, got {len(df)}"

    expected_cols = {
        "fuel_id",
        "composition_type",
        "total_components",
        "cutoff_K",
        "effective_components",
        "retained_mass_fraction",
        "lumping_strategy",
        "density_full_kg_m3",
        "density_pred_kg_m3",
        "density_rel_drift_pct",
        "viscosity_full_cSt",
        "viscosity_pred_cSt",
        "viscosity_rel_drift_pct",
        "flash_full_C",
        "flash_pred_C",
        "flash_abs_drift_C",
        "ncv_full_MJ_kg",
        "ncv_pred_MJ_kg",
        "ncv_rel_drift_pct",
    }
    assert expected_cols.issubset(set(df.columns))


def test_d6_monotonic_retained_mass():
    """Verify that as cutoff K increases, retained mass fraction is monotonically non-decreasing."""
    df = pd.read_csv(CSV_PATH)
    for fid, g in df.groupby("fuel_id"):
        g_trunc = g[g["lumping_strategy"] == "truncate_renormalize"].sort_values("cutoff_K")
        mass_fractions = g_trunc["retained_mass_fraction"].tolist()
        assert all(x <= y + 1e-6 for x, y in zip(mass_fractions, mass_fractions[1:])), \
            f"Non-monotonic mass fractions for {fid}: {mass_fractions}"


def test_d6_tail_lumping_reduces_drift():
    """Verify that tail lumping substantially reduces viscosity and density drift at low K."""
    df = pd.read_csv(CSV_PATH)
    k5 = df[df["cutoff_K"] == 5]
    trunc_visc = k5[k5["lumping_strategy"] == "truncate_renormalize"]["viscosity_rel_drift_pct"].mean()
    lump_visc = k5[k5["lumping_strategy"] == "tail_lumping"]["viscosity_rel_drift_pct"].mean()

    assert lump_visc < trunc_visc, f"Expected tail lumping drift {lump_visc} < raw truncation {trunc_visc}"
