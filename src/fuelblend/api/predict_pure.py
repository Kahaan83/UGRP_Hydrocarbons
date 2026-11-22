"""Pure-compound property prediction with the best method per property (see `best_methods`).

Fitted models live in a gitignored cache (`data/cache/pure_models/`). Random forests fit on first use (seconds, fixed
seed). The D-MPNN ensemble is fitted by `scripts/v7/final/fit_pure_models.py` (GPU recommended). Training rows are the
Track B dev pool with the locked ysi_holdout labels masked. Errors are scaffold-CV MAE on dev rows (optimistic); there is
no conformal margin for pure compounds.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from fuelblend.api import best_methods as BM

ROOT = Path(__file__).resolve().parents[3]
CACHE = ROOT / "data/cache/pure_models"
RF_TREES, RF_SEED = 500, 0
_MODELS: dict = {}          # in-process cache of loaded models, keyed by file path


def _rf_path(prop: str) -> Path:
    return CACHE / f"rf_{prop}.joblib"


def fit_rf(prop: str, comp=None):
    """Fit and cache the RF for one property on all dev rows."""
    import joblib
    from sklearn.ensemble import RandomForestRegressor
    from fuelblend.features.rdkit_desc import rdkit_2d
    from fuelblend.v7 import runner
    if comp is None:
        comp, _, _ = runner.load_track_b_pool()
    y = comp[prop].values
    has = ~np.isnan(y)
    X = rdkit_2d(comp.smiles[has].tolist())
    cols = X.columns[(X.nunique() > 1).values]
    X = X[cols]
    med = X.median()
    rf = RandomForestRegressor(n_estimators=RF_TREES, n_jobs=4, random_state=RF_SEED).fit(X.fillna(med), y[has])
    CACHE.mkdir(parents=True, exist_ok=True)
    joblib.dump(dict(rf=rf, cols=list(cols), med=med.to_dict(), n_train=int(has.sum())), _rf_path(prop), compress=3)
    return _rf_path(prop)


def _predict_rf(prop: str, smiles: list) -> np.ndarray:
    import joblib
    import pandas as pd
    from fuelblend.features.rdkit_desc import rdkit_2d
    if not _rf_path(prop).exists():
        fit_rf(prop)
    key = str(_rf_path(prop))
    if key not in _MODELS:
        _MODELS[key] = joblib.load(_rf_path(prop))
    b = _MODELS[key]
    X = rdkit_2d(smiles)[b["cols"]]
    return b["rf"].predict(X.fillna(pd.Series(b["med"])))


def _predict_dmpnn(prop: str, smiles: list) -> np.ndarray:
    import torch
    from fuelblend.ml import dmpnn as D
    p, meta = CACHE / "dmpnn_long.pt", CACHE / "dmpnn_long.json"
    if not p.exists() or not meta.exists():
        raise FileNotFoundError("D-MPNN weights missing; run scripts/v7/final/fit_pure_models.py "
                                "(GPU venv, FUELBLEND_DEVICE=cuda)")
    m = json.loads(meta.read_text())
    if str(p) not in _MODELS:
        nets = []
        for sd in torch.load(p, map_location=D.DEVICE):
            net = D.DMPNN(len(m["tasks"]), **m["arch"]).to(D.DEVICE)
            net.load_state_dict(sd)
            net.eval()
            nets.append(net)
        _MODELS[str(p)] = nets
    nets = _MODELS[str(p)]
    j = m["tasks"].index(prop)
    pred, _ = D.predict(nets, [D.featurize(s) for s in smiles])
    return pred[:, j] * m["sd"][j] + m["mu"][j]


def predict_pure(smiles, props=None) -> dict:
    """{smiles: {prop: {value, unit, method, evidence_mae, significance, status}}} for the best method of each property.

    props defaults to all ten. A property whose model is unavailable gets status 'model_missing' and no value.
    Unparseable SMILES give status 'invalid_smiles'."""
    from rdkit import Chem
    smi = [smiles] if isinstance(smiles, str) else list(smiles)
    props = list(props) if props else list(BM.PURE)
    ok = [Chem.MolFromSmiles(s) is not None for s in smi]
    good = [s for s, o in zip(smi, ok) if o]
    out = {s: {} for s in smi}
    for prop in props:
        e = BM.best_method(prop)
        base = dict(unit=e["unit"], method=e["method"], evidence_mae=BM.evidence_mae(prop),
                    evidence="scaffold 5-fold CV, dev rows (optimistic)", significance=e["significance"])
        vals, err = None, None
        if good:
            try:
                vals = (_predict_rf if e["method"] == "rf_rdkit2d" else _predict_dmpnn)(prop, good)
            except FileNotFoundError as ex:
                err = str(ex)
        it = iter(vals) if vals is not None else None
        for s, o in zip(smi, ok):
            if not o:
                out[s][prop] = dict(base, value=float("nan"), status="invalid_smiles")
            elif it is None:
                out[s][prop] = dict(base, value=float("nan"), status="model_missing", reason=err)
            else:
                out[s][prop] = dict(base, value=float(next(it)), status="ok")
    return out
