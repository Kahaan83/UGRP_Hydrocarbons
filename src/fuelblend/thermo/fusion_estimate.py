"""Tm / dSfus estimator for hydrocarbons (v5 delta_ws3_pure_gaps Part B).  Estimates only; never replaces a measured value.

``family_of`` is the grouped-CV key declared in results/v5/BASE/delta_ws3_pure_gaps/README.md (fixed before any fit).
"""
from __future__ import annotations

import warnings

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors, GraphDescriptors, rdMolDescriptors

RDLogger.DisableLog("rdApp.*")
SOURCE = "ml_estimated"
FAMILIES = ("n_alkane", "iso_alkane", "acyclic_unsaturated", "cycloalkane_mono", "cycloalkane_poly", "cycloalkene",
            "benzene_mono", "aromatic_naphthenic", "aromatic_poly")
LGB_PARAMS = dict(n_estimators=400, learning_rate=0.03, num_leaves=15, min_child_samples=5, random_state=0, verbose=-1)


def _mol(smiles):
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        raise ValueError(f"bad SMILES {smiles!r}")
    return m


def family_of(smiles: str) -> str:
    m = _mol(smiles)
    ri = m.GetRingInfo()
    arom_rings = [r for r in ri.AtomRings() if all(m.GetAtomWithIdx(i).GetIsAromatic() for i in r)]
    sat_rings = [r for r in ri.AtomRings() if any(not m.GetAtomWithIdx(i).GetIsAromatic() for i in r)]   # non-aromatic ring (may share atoms)
    unsat_nonarom = any(b.GetBondType() != Chem.BondType.SINGLE and not b.GetIsAromatic() for b in m.GetBonds())
    if arom_rings:
        if len(arom_rings) >= 2:
            return "aromatic_poly"
        return "aromatic_naphthenic" if sat_rings else "benzene_mono"
    if ri.NumRings() == 0:
        if unsat_nonarom:
            return "acyclic_unsaturated"
        branched = any(a.GetSymbol() == "C" and a.GetDegree() > 2 for a in m.GetAtoms())
        return "iso_alkane" if branched else "n_alkane"
    if unsat_nonarom:
        return "cycloalkene"
    return "cycloalkane_mono" if ri.NumRings() == 1 else "cycloalkane_poly"


def n_carbon(smiles: str) -> int:
    return sum(a.GetSymbol() == "C" for a in _mol(smiles).GetAtoms())


def _joback(m):
    try:
        from thermo.group_contribution.joback import Joback
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            J = Joback(m)
            tm, hf = J.Tm(J.counts), J.Hfus(J.counts)
        return (float(tm) if tm is not None else np.nan, float(hf) / 1000.0 if hf is not None else np.nan)
    except Exception:  # noqa: BLE001 - Joback cannot fragment some structures
        return np.nan, np.nan


def _wiener(m):
    d = Chem.GetDistanceMatrix(m)
    return float(d.sum() / 2.0)


def features(smiles_list) -> pd.DataFrame:
    rows = []
    for s in smiles_list:
        m = _mol(s)
        ranks = list(Chem.CanonicalRankAtoms(m, breakTies=False))
        nh = m.GetNumHeavyAtoms()
        jt, jh = _joback(m)
        rows.append(dict(
            n_carbon=n_carbon(s), mw=Descriptors.MolWt(m), n_rings=rdMolDescriptors.CalcNumRings(m),
            n_arom_rings=rdMolDescriptors.CalcNumAromaticRings(m), frac_sp3=rdMolDescriptors.CalcFractionCSP3(m),
            n_rot=rdMolDescriptors.CalcNumRotatableBonds(m),
            n_ch3=sum(a.GetSymbol() == "C" and a.GetTotalNumHs() == 3 for a in m.GetAtoms()),
            n_quat=sum(a.GetSymbol() == "C" and a.GetDegree() == 4 for a in m.GetAtoms()),
            wiener=_wiener(m), kappa1=GraphDescriptors.Kappa1(m), kappa2=GraphDescriptors.Kappa2(m), kappa3=GraphDescriptors.Kappa3(m),
            balaban=GraphDescriptors.BalabanJ(m), symmetry=len(set(ranks)) / nh if nh else np.nan,
            joback_tm=jt, joback_dhfus=jh))
    return pd.DataFrame(rows)


def _model():
    import lightgbm as lgb
    return lgb.LGBMRegressor(**LGB_PARAMS)


def oof_predict(X: pd.DataFrame, y: np.ndarray, groups) -> np.ndarray:
    """Leave-one-group-out out-of-fold predictions (every group is one fold)."""
    groups = np.asarray(groups)
    pred = np.full(len(y), np.nan)
    for g in np.unique(groups):
        te = groups == g
        if te.all():
            continue
        mdl = _model().fit(X[~te], y[~te])
        pred[te] = mdl.predict(X[te])
    return pred


def fit_predict(X: pd.DataFrame, y: np.ndarray, Xnew: pd.DataFrame) -> np.ndarray:
    return _model().fit(X, y).predict(Xnew)


def fill(measured: dict, estimates: dict) -> dict:
    """Merge: a measured entry always wins; estimates only fill keys that have no measured entry."""
    out = {k: dict(v) for k, v in estimates.items() if k not in measured}
    out.update({k: dict(v) for k, v in measured.items()})
    return out
