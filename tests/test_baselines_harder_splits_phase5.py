"""Unit tests for Phase 5 Literature Baselines, Harder Splits, and Multi-Seed Runs (T5.3 - T5.6)."""
import json
from pathlib import Path
import pandas as pd
import numpy as np

from fuelblend.baselines.literature_baselines import (
    grunberg_nissan_viscosity,
    kendall_monroe_viscosity,
    redlich_kister_density,
)

ROOT = Path(__file__).resolve().parents[1]
T5_DIR = ROOT / "results" / "t5"
SPLITS_DIR = ROOT / "splits"


def test_literature_mixing_rules():
    # Equimolar binary of octane (eta=0.54 mPa*s, rho=703 kg/m3, MW=114.2)
    # and toluene (eta=0.59 mPa*s, rho=867 kg/m3, MW=92.14)
    x = [0.5, 0.5]
    eta = [0.54, 0.59]
    rho = [703.0, 867.0]
    mw = [114.2, 92.14]

    eta_gn = grunberg_nissan_viscosity(x, eta, G_interaction=0.0)
    assert 0.54 < eta_gn < 0.59

    eta_km = kendall_monroe_viscosity(x, eta)
    assert 0.54 < eta_km < 0.59

    rho_rk = redlich_kister_density(x, mw, rho, A0=0.5, A1=0.0)
    assert 703.0 < rho_rk < 867.0


def test_harder_splits_json():
    p = SPLITS_DIR / "harder_splits.json"
    assert p.exists()
    data = json.loads(p.read_text(encoding="utf-8"))
    assert "unseen_n_components" in data
    assert "unseen_temperature_ranges" in data
    assert "unseen_component_pairs" in data

    # Verify split separation
    train_n2 = set(data["unseen_n_components"]["train_system_ids"])
    test_n3 = set(data["unseen_n_components"]["test_system_ids"])
    assert len(train_n2.intersection(test_n3)) == 0


def test_multiseed_ensemble_csv():
    p = T5_DIR / "multiseed_ensemble_eval.csv"
    assert p.exists()
    df = pd.read_csv(p)
    assert len(df) == 7  # 5 seeds + mean + std
    assert "MEAN (5 seeds)" in df["seed"].values
    assert "SPREAD (std dev)" in df["seed"].values
