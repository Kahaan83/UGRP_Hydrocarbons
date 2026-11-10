"""v4 Task 15: octane models for Co-Optima bRON blends (base-fuel lump + one measured blendstock).

Every model maps one blend row to a predicted RON (or MON):
  * ``f``       blendstock VOLUME fraction (measured, exact),
  * ``base_on`` the base gasoline's measured RON (MON) -- a model INPUT (ruling R17),
  * blendstock descriptors (lump class, structure) -- NOT the blendstock's own bRON, which is the
    quantity being predicted and unknown for a held-out blendstock (ruling R19).

All hyper-parameter selection happens inside ``fit(train)`` using only ``train`` (ruling R6).
"""
from __future__ import annotations

import json
import re
from functools import lru_cache

import numpy as np
import pandas as pd
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import ConstantKernel, RBF, WhiteKernel

from fuelblend.models import octane_ghosh_lumps as G

# Assumed bulk properties of the (compositionless) base gasoline, used ONLY by linear-by-mole.
# Typed typical values, not measured here; disclosed in the README.
BASE_MW_ASSUMED = 100.0
BASE_RHO_ASSUMED = 745.0     # kg/m3 at 15 C

FEATURES = ["f", "base_ron", "base_mon", "base_sens", "mw", "rho", "n_c", "o_mass_frac",
            "unsat", "ring", "arom"]


@lru_cache(maxsize=None)
def blendstock_descriptors(cas: str | None, name: str = "") -> dict:
    """Structure descriptors of a blendstock from its CAS (chemicals/RDKit); NaN for named mixtures."""
    out = dict(smiles=None, mw=np.nan, rho=np.nan, n_c=np.nan, o_mass_frac=np.nan,
               unsat=np.nan, ring=np.nan, arom=np.nan)
    if cas:
        from chemicals.identifiers import search_chemical
        from rdkit import Chem
        from rdkit.Chem import Descriptors

        smi = search_chemical(cas).smiles
        m = Chem.MolFromSmiles(smi)
        mw = Descriptors.MolWt(m)
        n_o = sum(a.GetSymbol() == "O" for a in m.GetAtoms())
        out.update(smiles=smi, mw=mw, n_c=sum(a.GetSymbol() == "C" for a in m.GetAtoms()),
                   o_mass_frac=16.0 * n_o / mw,
                   unsat=float(m.HasSubstructMatch(Chem.MolFromSmarts("C=C"))),
                   ring=float(m.GetRingInfo().NumRings() > 0),
                   arom=float(any(a.GetIsAromatic() for a in m.GetAtoms())))
        try:
            from thermo import Chemical
            out["rho"] = float(Chemical(cas, T=288.15).rho)
        except Exception:                                  # density optional (linear-by-mole only)
            pass
    return out


def prepare_blends(df: pd.DataFrame, prop: str) -> tuple[pd.DataFrame, dict]:
    """Rows usable as ``prop`` targets, as a flat model table, plus an exclusion report.

    Uses the per-property flag ``target_ok_<prop>`` (never the OR column ``target_ok``); a row
    must also carry the measured value and the matching base value.  Rows with an unresolved
    (``frac=None``) component are dropped and counted, never imputed.
    """
    assert prop in ("ron", "mon")
    flag = f"target_ok_{prop}"
    rep = dict(n_rows=len(df), n_flag_false=int((~df[flag].astype(bool)).sum()))
    d = df[df[flag].astype(bool)]
    rep["n_flag_true"] = len(d)
    rows, n_none, n_shape, n_nan = [], 0, 0, 0
    for r in d.to_dict("records"):
        comps = json.loads(r["components"])
        if any(c["frac"] is None for c in comps):
            n_none += 1
            continue
        lump = [c for c in comps if str(c["inchikey"]).startswith("LUMP:BRON_BASE")]
        stock = [c for c in comps if c not in lump]      # named mixtures (ketone/furan) stay blendstocks
        base_on = r[f"base_{prop}"]
        if lump:                                           # lump id carries BOTH base values: BRON_BASE_RON90.3_MON84.7
            mm = re.search(r"RON([\d.]+)_MON([\d.]+)", lump[0]["inchikey"])
            if mm:
                r = dict(r, base_ron=float(mm.group(1)), base_mon=float(mm.group(2)))
                base_on = r[f"base_{prop}"]
        if len(stock) != 1 or len(lump) != 1:
            n_shape += 1
            continue
        if pd.isna(r[prop]) or pd.isna(base_on) or pd.isna(r["base_ron"]) or pd.isna(r["base_mon"]):
            n_nan += 1
            continue
        group = r["blend_id"].rsplit(":", 1)[0]            # BRON:<CAS> or BRON:LUMP:<NAME>
        cas = None if group.startswith("BRON:LUMP:") else group.split(":", 1)[1]
        ds = blendstock_descriptors(cas, stock[0]["name"])
        rows.append(dict(
            blend_id=r["blend_id"], group=group, cas=cas, name=stock[0]["name"], f=float(stock[0]["frac"]),
            base_on=float(base_on), base_ron=float(r["base_ron"]), base_mon=float(r["base_mon"]),
            base_sens=float(r["base_ron"] - r["base_mon"]), y=float(r[prop]),
            lump=G.blendstock_lump(ds["smiles"], stock[0]["name"]), **ds))
    rep.update(n_dropped_frac_none=n_none, n_dropped_shape=n_shape, n_dropped_nan=n_nan, n_used=len(rows))
    return pd.DataFrame(rows), rep


def linear_by_volume_pure(components: list[tuple[str, float]], prop: str = "ron") -> float:
    """Linear-by-volume octane of a blend of components with literature lump priors.

    ``components`` = [(smiles, volume_fraction), ...]; for an n-heptane/iso-octane PRF this is
    exactly the vol % iso-octane (the definition of the PRF scale).
    """
    tot = sum(v for _, v in components)
    return float(sum(v / tot * G.lump_prior(G.blendstock_lump(s), prop) for s, v in components))


# ------------------------------------------------------------------------------- models
class Model:
    name = "?"

    def fit(self, train):
        self.selection_groups_ = set(train["group"])
        return self


class ConstantMean(Model):
    name = "constant_mean"

    def fit(self, train):
        super().fit(train)
        self.m = float(train["y"].mean())
        return self

    def predict(self, d):
        return np.full(len(d), self.m)


class BaseUnchanged(Model):
    name = "base_unchanged"

    def predict(self, d):
        return d["base_on"].to_numpy(float)


class LinearVolumePrior(Model):
    name = "linear_volume_prior"

    def __init__(self, prop):
        self.prop = prop

    def predict(self, d):
        f = d["f"].to_numpy(float)
        pri = np.array([G.lump_prior(l, self.prop) for l in d["lump"]])
        return (1 - f) * d["base_on"].to_numpy(float) + f * pri


class LinearMolePrior(Model):
    """Mole-fraction linear blending with lump priors; needs blendstock MW and density (NaN otherwise)."""
    name = "linear_mole_prior"

    def __init__(self, prop):
        self.prop = prop

    def predict(self, d):
        f = d["f"].to_numpy(float)
        nb = f * d["rho"].to_numpy(float) / d["mw"].to_numpy(float)
        na = (1 - f) * BASE_RHO_ASSUMED / BASE_MW_ASSUMED
        x = nb / (nb + na)
        pri = np.array([G.lump_prior(l, self.prop) for l in d["lump"]])
        return (1 - x) * d["base_on"].to_numpy(float) + x * pri


def _bov(d):
    return (d["y"].to_numpy(float) - (1 - d["f"].to_numpy(float)) * d["base_on"].to_numpy(float)) / d["f"].to_numpy(float)


class LinearVolumePooledBOV(Model):
    """Linear-by-volume with ONE blending octane value (mean of training BOVs) for every blendstock."""
    name = "linear_volume_pooled_bov"

    def fit(self, train):
        super().fit(train)
        self.bov = float(np.mean(_bov(train)))
        return self

    def predict(self, d):
        f = d["f"].to_numpy(float)
        return (1 - f) * d["base_on"].to_numpy(float) + f * self.bov


class GhoshLumps(Model):
    name = "ghosh_lumps"

    def __init__(self, prop):
        self.prop = prop

    def fit(self, train):
        self.m = G.fit(train, self.prop)
        self.selection_groups_ = self.m.selection_groups_
        return self

    def predict(self, d):
        return self.m.predict(d)


def _xy(d):
    return d[FEATURES].to_numpy(float)


class _Residual(Model):
    """Learner on the residual of the pooled-BOV linear baseline."""

    def _resid(self, d, bov):
        f = d["f"].to_numpy(float)
        return d["y"].to_numpy(float) - ((1 - f) * d["base_on"].to_numpy(float) + f * bov)

    def _base(self, d):
        f = d["f"].to_numpy(float)
        return (1 - f) * d["base_on"].to_numpy(float) + f * self.bov


class GBMResidual(_Residual):
    name = "gbm_resid"
    GRID = [dict(max_iter=n, max_depth=dp, learning_rate=0.1) for n in (20, 60) for dp in (1, 2)]

    def __init__(self, params=None):
        self.params = params

    def _new(self, p):
        return HistGradientBoostingRegressor(min_samples_leaf=5, random_state=0, **p)

    def fit(self, train):
        super().fit(train)
        train = train.reset_index(drop=True)
        self.bov = float(np.mean(_bov(train)))
        # re-selected on every fit() so a reused instance can never carry a previous fold's choice
        self.params_ = self.params if self.params is not None else             self.GRID[int(np.argmin(self._select(train)))]
        self.m = self._new(self.params_).fit(_xy(train), self._resid(train, self.bov))
        return self

    def _select(self, train):
        # inner grouped CV on `train` only; each candidate is a fixed-parameter GBMResidual
        return G.grouped_inner_scores(lambda p: GBMResidual(p), self.GRID, train)

    def predict(self, d):
        return self._base(d) + self.m.predict(_xy(d))


class GPResidual(_Residual):
    name = "gp_resid"

    def fit(self, train):
        super().fit(train)
        train = train.reset_index(drop=True)
        self.bov = float(np.mean(_bov(train)))
        X = _xy(train)
        self.mu = np.nanmean(X, axis=0)
        self.mu = np.where(np.isnan(self.mu), 0.0, self.mu)
        Xi = np.where(np.isnan(X), self.mu, X)
        self.sd = Xi.std(axis=0)
        self.sd[self.sd == 0] = 1.0
        Z = (Xi - self.mu) / self.sd
        k = ConstantKernel(4.0, (1e-2, 1e3)) * RBF(np.ones(Z.shape[1]) * 3.0, (0.5, 1e2)) + \
            WhiteKernel(1.0, (1e-3, 1e2))
        # kernel hyper-parameters come from marginal likelihood on the TRAINING rows only
        self.m = GaussianProcessRegressor(k, normalize_y=False, n_restarts_optimizer=0, random_state=0)
        import warnings
        from sklearn.exceptions import ConvergenceWarning
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', ConvergenceWarning)   # length-scale at bound = feature ignored
            self.m.fit(Z, self._resid(train, self.bov))
        return self

    def predict(self, d):
        X = _xy(d)
        Z = (np.where(np.isnan(X), self.mu, X) - self.mu) / self.sd
        return self._base(d) + self.m.predict(Z)


def make_models(prop: str) -> list[Model]:
    return [ConstantMean(), BaseUnchanged(), LinearVolumePrior(prop), LinearMolePrior(prop),
            LinearVolumePooledBOV(), GhoshLumps(prop), GBMResidual(), GPResidual()]


def fit_predict(model: Model, train: pd.DataFrame, test: pd.DataFrame) -> np.ndarray:
    """Fit on ``train`` only and predict ``test`` (test labels are never read)."""
    t = test.drop(columns=["y"]) if "y" in test.columns else test
    return np.asarray(model.fit(train).predict(t), float)


# ------------------------------------------------------------------------------- metrics
def metrics(y, p) -> dict:
    y, p = np.asarray(y, float), np.asarray(p, float)
    ok = np.isfinite(p)
    e = p[ok] - y[ok]
    if not len(e):
        return dict(n=0, mae=np.nan, rmse=np.nan, bias=np.nan, median_abs=np.nan)
    return dict(n=int(ok.sum()), mae=float(np.abs(e).mean()), rmse=float(np.sqrt((e ** 2).mean())),
                bias=float(e.mean()), median_abs=float(np.median(np.abs(e))))


def cluster_boot(groups, a, b=None, n_boot=4000, seed=0):
    """Cluster bootstrap over groups of mean |e_a| (or mean |e_a| - mean |e_b|, paired). Rows with a
    non-finite error in either series are dropped first.  Returns (point, lo, hi, n_groups)."""
    g = np.asarray(groups)
    a = np.abs(np.asarray(a, float))
    ok = np.isfinite(a)
    if b is not None:
        b = np.abs(np.asarray(b, float))
        ok &= np.isfinite(b)
    g, a = g[ok], a[ok]
    v = a - (b[ok] if b is not None else 0.0)
    ug = np.unique(g)
    S = np.array([v[g == u].sum() for u in ug])
    N = np.array([(g == u).sum() for u in ug], float)
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, len(ug), (n_boot, len(ug)))
    est = S[idx].sum(1) / N[idx].sum(1)
    return float(S.sum() / N.sum()), float(np.percentile(est, 2.5)), float(np.percentile(est, 97.5)), len(ug)
