import numpy as np

from fuelblend.v8.ml_constants import ml_constant, oof_constants


def test_ml_layer_returns_ml_status():
    v, st = ml_constant("CCCCCCCCCCCC", "Tm")      # dodecane
    assert st.startswith("ml_") and np.isfinite(v)


def test_out_of_domain_oxygenate():
    assert ml_constant("CCO", "rho15")[1] == "ml_out_of_domain"


def test_invalid_smiles():
    assert ml_constant("not_a_smiles", "Tm")[1] == "invalid_smiles"


def test_oof_has_no_in_sample_rows():
    df = oof_constants("Tm")
    assert df["fold"].notna().all() and df["inchikey"].is_unique
