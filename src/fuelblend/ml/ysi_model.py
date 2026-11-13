"""Learned group-contribution model of pure-compound yield sooting index (YSI), v5 delta_ml_ysi.

Ridge regression on 16 RDKit group counts (StandardScaler inside, so the model stays linear in the counts and
extrapolates linearly in carbon number), raw-YSI target, sample weight 1/(YSI+50)^2.  Declared in
results/v5/BASE/delta_ml_ysi/README.md.  ``predict_ysi`` NEVER overwrites a measured value: compounds in the 224-compound
table are returned as ``measured``; the rest are ``ml_predicted`` (inside the strict tier) or ``ml_predicted_extrapolated``.
"""
from __future__ import annotations

import json
from functools import lru_cache
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors
from rdkit.Chem.Scaffolds import MurckoScaffold

ROOT = Path(__file__).resolve().parents[3]
MODEL_DIR = ROOT / "models" / "v5" / "ml_ysi"
YSI_CSV = ROOT / "data" / "processed" / "v4" / "sooting" / "ysi_hc.csv"

FEATURES = ["n_C", "n_H", "H_per_C", "MW", "arom_C", "arom_rings", "aliph_rings", "fused_ring_atoms", "sp2_nonarom_C", "sp_C",
            "CH3", "tertiary_CH", "quaternary_C", "benzylic_C", "rot_bonds", "DBE"]
ALPHAS = (0.1, 1.0, 10.0)
CLASSES = ["paraffin", "iso", "other", "cyclo", "mono_arom", "poly_arom"]
STRICT_CARBON_MARGIN = 2


def canon(smiles: str) -> str:
    """Canonical SMILES without stereo; raises ValueError if unparseable."""
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        raise ValueError(f"unparseable SMILES {smiles!r}")
    return Chem.MolToSmiles(m, isomericSmiles=False)


def _is_hc(m) -> bool:
    return all(a.GetSymbol() in ("C", "H") for a in m.GetAtoms())


def n_carbon(m) -> int:
    return sum(a.GetSymbol() == "C" for a in m.GetAtoms())


def features(smiles: str) -> np.ndarray:
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        raise ValueError(f"unparseable SMILES {smiles!r}")
    mh = Chem.AddHs(m)
    nC = n_carbon(m)
    nH = sum(a.GetSymbol() == "H" for a in mh.GetAtoms())
    ri = m.GetRingInfo()
    n_arom = rdMolDescriptors.CalcNumAromaticRings(m)
    n_aliph = rdMolDescriptors.CalcNumRings(m) - n_arom
    ring_count = np.array([ri.NumAtomRings(a.GetIdx()) for a in m.GetAtoms()])
    arom_C = sp2 = sp = ch3 = tert = quat = benz = 0
    for a in m.GetAtoms():
        if a.GetSymbol() != "C":
            continue
        hyb = a.GetHybridization()
        nb = a.GetDegree()
        if a.GetIsAromatic():
            arom_C += 1
        elif hyb == Chem.HybridizationType.SP2:
            sp2 += 1
        elif hyb == Chem.HybridizationType.SP:
            sp += 1
        if hyb == Chem.HybridizationType.SP3 and not a.GetIsAromatic():
            if nb == 1:
                ch3 += 1
            if nb == 3:
                tert += 1
            if nb == 4:
                quat += 1
            if any(n.GetIsAromatic() for n in a.GetNeighbors()):
                benz += 1
    dbe = (2 * nC + 2 - nH) / 2.0
    return np.array([nC, nH, nH / max(nC, 1), Descriptors.MolWt(m), arom_C, n_arom, n_aliph, int((ring_count >= 2).sum()), sp2, sp,
                     ch3, tert, quat, benz, rdMolDescriptors.CalcNumRotatableBonds(m), dbe], float)


def mol_class(smiles: str) -> str:
    m = Chem.MolFromSmiles(smiles)
    n_arom = rdMolDescriptors.CalcNumAromaticRings(m)
    if n_arom >= 2:
        return "poly_arom"
    if n_arom == 1:
        return "mono_arom"
    if rdMolDescriptors.CalcNumRings(m) > 0:
        return "cyclo"
    unsat = any(b.GetBondType() != Chem.BondType.SINGLE for b in m.GetBonds())
    if unsat:
        return "other"
    branched = any(a.GetDegree() >= 3 for a in m.GetAtoms())
    return "iso" if branched else "paraffin"


def is_saturated_acyclic_or_cyclic(smiles: str) -> bool:
    """No double/triple bond, no aromatic ring."""
    m = Chem.MolFromSmiles(smiles)
    return all(b.GetBondType() == Chem.BondType.SINGLE for b in m.GetBonds())


def is_nonaromatic_ring(smiles: str) -> bool:
    m = Chem.MolFromSmiles(smiles)
    return rdMolDescriptors.CalcNumRings(m) > 0 and rdMolDescriptors.CalcNumAromaticRings(m) == 0


def group_key(smiles: str) -> str:
    """Generic Murcko scaffold for ring compounds, ``acyclic_C{n}`` otherwise."""
    m = Chem.MolFromSmiles(smiles)
    if rdMolDescriptors.CalcNumRings(m) == 0:
        return f"acyclic_C{n_carbon(m)}"
    sc = MurckoScaffold.GetScaffoldForMol(m)
    return Chem.MolToSmiles(MurckoScaffold.MakeScaffoldGeneric(sc))


def scaffold_generic(smiles: str) -> str | None:
    m = Chem.MolFromSmiles(smiles)
    if rdMolDescriptors.CalcNumRings(m) == 0:
        return None
    return Chem.MolToSmiles(MurckoScaffold.MakeScaffoldGeneric(MurckoScaffold.GetScaffoldForMol(m)))


def load_training(path: Path = YSI_CSV) -> pd.DataFrame:
    """Canonical (stereo-free) SMILES, duplicates collapsed to the median YSI."""
    d = pd.read_csv(path)
    d["smiles_c"] = [canon(s) for s in d.smiles]
    g = d.groupby("smiles_c", as_index=False).agg(ysi=("ysi", "median"), name=("name", "first"))
    return g.rename(columns={"smiles_c": "smiles"})


def weights(y) -> np.ndarray:
    return 1.0 / (np.asarray(y, float) + 50.0) ** 2


# ---------------------------------------------------------------------------------------------- fitting
class RidgeGC:
    """Standardised weighted ridge; stores plain arrays so it can be saved as JSON."""

    def __init__(self, alpha: float = 1.0):
        self.alpha = alpha

    def fit(self, X, y, w):
        X, y, w = np.asarray(X, float), np.asarray(y, float), np.asarray(w, float)
        self.mean_ = X.mean(0)
        sd = X.std(0)
        self.scale_ = np.where(sd > 0, sd, 1.0)
        Z = (X - self.mean_) / self.scale_
        sw = w / w.sum()
        ym = float(np.dot(sw, y))
        zm = sw @ Z
        Zc, yc = Z - zm, y - ym
        A = Zc.T @ (Zc * sw[:, None]) + self.alpha * np.eye(Z.shape[1]) / len(y)
        self.coef_ = np.linalg.solve(A, Zc.T @ (yc * sw))
        self.intercept_ = ym - float(zm @ self.coef_)
        return self

    def predict(self, X):
        Z = (np.asarray(X, float) - self.mean_) / self.scale_
        return np.clip(Z @ self.coef_ + self.intercept_, 0.0, None)

    def to_dict(self):
        return dict(alpha=self.alpha, mean=self.mean_.tolist(), scale=self.scale_.tolist(), coef=self.coef_.tolist(), intercept=self.intercept_)

    @classmethod
    def from_dict(cls, d):
        o = cls(d["alpha"])
        o.mean_, o.scale_, o.coef_, o.intercept_ = np.array(d["mean"]), np.array(d["scale"]), np.array(d["coef"]), float(d["intercept"])
        return o


def select_alpha(X, y, groups, n_inner: int = 4) -> float:
    """Inner GroupKFold(4), criterion = weighted MSE with the fit weights."""
    from sklearn.model_selection import GroupKFold
    w = weights(y)
    groups = np.asarray(groups)
    k = min(n_inner, len(np.unique(groups)))
    best, best_s = ALPHAS[0], np.inf
    for a in ALPHAS:
        s = 0.0
        for tr, te in GroupKFold(k).split(X, y, groups):
            m = RidgeGC(a).fit(X[tr], y[tr], w[tr])
            s += float(np.sum(w[te] * (m.predict(X[te]) - y[te]) ** 2))
        if s < best_s - 1e-15:
            best, best_s = a, s
    return best


def fit_ridge(X, y, groups) -> RidgeGC:
    return RidgeGC(select_alpha(X, y, groups)).fit(X, y, weights(y))


# ---------------------------------------------------------------------------------------------- inference
@lru_cache(maxsize=1)
def _load():
    d = json.loads((MODEL_DIR / "model.json").read_text(encoding="utf-8"))
    return RidgeGC.from_dict(d["ridge"]), d


def tier(smiles: str, meta: dict) -> str:
    """strict | extrapolated for a compound without measured YSI."""
    cls = mol_class(smiles)
    nC = n_carbon(Chem.MolFromSmiles(smiles))
    ok = nC <= meta["class_max_nC"].get(cls, -99) + STRICT_CARBON_MARGIN
    sc = scaffold_generic(smiles)
    if sc is not None and sc not in set(meta["train_scaffolds"]):
        ok = False
    return "strict" if ok else "extrapolated"


def predict_ysi(smiles_list) -> pd.DataFrame:
    """DataFrame[smiles, ysi, lo, hi, source].  Hydrocarbons only; a non-hydrocarbon gets NaN and source ``out_of_domain``.

    source: ``measured`` (value from the 224-compound table, never modified) | ``ml_predicted`` | ``ml_predicted_extrapolated``.
    lo/hi: split-conformal 90 % interval (README, R-MLY-4); for extrapolated fills the larger of the OOF and extrapolation-test quantiles
    (descriptive, NOT calibrated for carbon-number extrapolation).  Unparseable SMILES raise ValueError."""
    model, meta = _load()
    tab = meta["measured"]
    rows = []
    for s in smiles_list:
        c = canon(s)
        m = Chem.MolFromSmiles(c)
        if not _is_hc(m):
            rows.append(dict(smiles=s, ysi=np.nan, lo=np.nan, hi=np.nan, source="out_of_domain"))
            continue
        if c in tab:
            v = float(tab[c])
            rows.append(dict(smiles=s, ysi=v, lo=v, hi=v, source="measured"))
            continue
        v = float(model.predict(features(c)[None, :])[0])
        ext = tier(c, meta) == "extrapolated"
        q = max(meta["q_oof"], meta["q_extrap"]) if ext else meta["q_oof"]
        rows.append(dict(smiles=s, ysi=v, lo=max(0.0, v - q * (v + 50.0)), hi=v + q * (v + 50.0),
                         source="ml_predicted_extrapolated" if ext else "ml_predicted"))
    return pd.DataFrame(rows, columns=["smiles", "ysi", "lo", "hi", "source"])
