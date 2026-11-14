"""v5 WS7 (results/v5/BASE/delta_ws7_design_realism): real blendstocks, D86 constraints and a blendstock optimiser.

* A ``Blendstock`` is a fixed recipe of registry compounds: a real fuel's GCxGC class x carbon-number bins mapped with the T2 route
  (``composition_route.kind_bins`` -> ``map_bins`` with ``reconstruction_r3.SubReps``), or a single pure compound (verified table name).
* ``BlendPool`` expands top-level fractions ``u`` to compound mass fractions ``w = A @ u``.  Every engine call (property models, recipe margins,
  graded penalties, extrapolation, class mix, by-construction aromatics / naphthalenes, blocked check) runs on the compound pool, so nothing in
  ``engine_v4`` / ``library_route`` needs to know about blendstocks.  A pool of one blendstock at weight 1 is the rep-level recipe exactly.
* D86 models (declared in the README): Tb from the default flash-registry Antoine at 101.325 kPa, ``baselines.distillation`` cuts on volume
  fractions (pure-grid rho at 288.15 K); FBP = TBP cut at 99 vol % (no offset); slope quantities as differences.  Distillation domain: monotone D86
  curve and TBP T90-T10 inside the margin-fit fuels' range; outside, the recipe margin is None (``outside_distillation_domain``).
* ``optimise``: multi-start SLSQP in ``u`` with the graded engine objective, an optional uncertainty penalty (alternate-model spread) and hard
  constraints (extrapolation, synthetic cap, distillation domain).
"""
from __future__ import annotations

import dataclasses
import functools
import json
import math
import warnings
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Optional

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import inchi as rdinchi
from scipy.optimize import minimize

from fuelblend.baselines.distillation import predict_d86_from_composition, predict_tbp_cuts
from fuelblend.baselines.flash_unifac import _lookup
from fuelblend.design import engine_v4 as E
from fuelblend.design import library_route as LR

ROOT = Path(__file__).resolve().parents[3]
V5 = ROOT / "data" / "real_fuels" / "v5"
DELTA = ROOT / "results" / "v5" / "BASE" / "delta_ws7_design_realism"
D86_MARGINS = DELTA / "d86_margins.json"
P_ATM_KPA = 101.325
T15_K = 288.15
FBP_CUT_PCT = 99.0
LAMBDA_UNC = 0.25
SYNTHETIC_CAP = 0.50
REPRESENTED_REASON = ("represented_blendstock_margin_not_validated: recipe contains a GCxGC-represented blendstock; margins were validated on "
                      "known compositions only (real-fuel margins not adopted, WS6)")
DOMAIN_REASON = "outside_distillation_domain"

# D86 quantity -> (constraint canonical name, limit side: 'hi' = maximum limit, 'lo' = minimum limit)
D86_PROPS = {"T10": ("distillation_T10", "hi"), "FBP": ("distillation_FBP", "hi"),
             "T50_T10": ("distillation_T50_T10", "lo"), "T90_T10": ("distillation_T90_T10", "lo")}

# declared blendstocks (README R-W7-1): name -> (source, fuel_id, pathway label)
BLENDSTOCKS = {
    "jetA_posf10325": ("ajftd", 2, "Jet A base (AJFTD GCxGC fuel 2, label 10325, POSF-10325)"),
    "hefa_came": ("fuellib", "FUELLIB_hefa-came", "HEFA-SPK (FuelLib HEFA-CAME)"),
    "atj_posf11498": ("fuellib", "FUELLIB_posf11498", "ATJ-SPK (FuelLib C-1 POSF11498 ATJ)"),
}
BLOCKED_BLENDSTOCKS = {"SAK": "no named SAK composition in data/real_fuels or data/raw (AJFTD fuels 47/139 are ~97 % aromatic but unlabelled)"}


# ------------------------------------------------------------------------------------------------ registry / reps
@functools.lru_cache(maxsize=1)
def _reg_prov_reps():
    from fuelblend.data.reconstruction_r3 import SubReps
    from fuelblend.mixture.common import CompoundRegistry, PurePropertyProvider
    reg = CompoundRegistry(ROOT / "data" / "unified" / "compounds.csv")
    prov = PurePropertyProvider(reg, ROOT / "data" / "cache" / "pure_grid.npz", verbose=False)
    reps = SubReps(reg, prov, vp_lookup=lambda s: _lookup(s, None))
    return reg, prov, reps


def _ik(smiles: str) -> str:
    return rdinchi.MolToInchiKey(Chem.MolFromSmiles(smiles))


@dataclass
class Blendstock:
    name: str
    label: str
    source: str
    fuel_id: object
    smiles: list                       # compound SMILES (or table names for pure compounds)
    w: np.ndarray                      # mass fractions over ``smiles`` (sum 1)
    unmapped_fraction: float = 0.0
    represented: bool = True           # True = GCxGC bins mapped to representatives (not the fuel's own molecules)
    synthetic: bool = False            # counts toward the D7566 synthetic-blendstock cap

    def class_fractions(self) -> dict:
        from fuelblend.eval.composition_route import hc_class
        out: dict = {}
        for s, wi in zip(self.smiles, self.w):
            m = Chem.MolFromSmiles(_table_smiles().get(s, s))
            k = hc_class(m)
            out[k] = out.get(k, 0.0) + float(wi)
        return out


@functools.lru_cache(maxsize=1)
def _table_smiles() -> dict:
    d = pd.read_csv(LR.PURE_PROPS_V5).drop_duplicates("compound")
    return dict(zip(d.compound, d.smiles))


def _table_ik() -> dict:
    """InChIKey -> table name; when several table names share an InChIKey the one with the most ``verified`` rows wins (ties: name order)."""
    pp = pd.read_csv(LR.PURE_PROPS_V5)
    nver = pp[pp.status == "verified"].groupby("compound").size()
    best: dict = {}
    for n, s in sorted(_table_smiles().items()):
        k = _ik(s)
        if k not in best or nver.get(n, 0) > nver.get(best[k], 0):
            best[k] = n
    return best


def load_blendstock(name: str) -> Blendstock:
    from fuelblend.eval.composition_route import kind_bins, map_bins
    if name not in BLENDSTOCKS:
        raise KeyError(f"{name}: not a declared blendstock ({sorted(BLENDSTOCKS)}); blocked: {BLOCKED_BLENDSTOCKS}")
    src, fid, label = BLENDSTOCKS[name]
    reg, _prov, reps = _reg_prov_reps()
    if src == "ajftd":
        df = pd.read_csv(V5 / "ajftd_gcxgc_composition.csv")
        sub = df[df.fuel_id == fid]
    else:
        df = pd.read_csv(V5 / "fuellib_gcxgc_compositions.csv")
        sub = df[df.fuel_id == fid]
    if sub.empty:
        raise ValueError(f"{name}: no composition rows")
    m = map_bins(kind_bins(src, sub), reps)
    return Blendstock(name=name, label=label, source=src, fuel_id=fid, smiles=[reg.smiles[i] for i in m.idx], w=np.asarray(m.w, float),
                      unmapped_fraction=float(m.unmapped_fraction), represented=True, synthetic=name != "jetA_posf10325")


def pure_component(table_name: str) -> Blendstock:
    if table_name not in _table_smiles():
        raise KeyError(f"{table_name}: not in pure_props_v5 (pure components must be verified-table names)")
    return Blendstock(name=table_name, label=f"pure {table_name}", source="pure_props_v5", fuel_id=None, smiles=[table_name],
                      w=np.array([1.0]), represented=False)


class BlendPool:
    """Top-level members -> compound pool ids (table names win over SMILES for the same InChIKey) and matrix A (n_comp x n_top)."""

    def __init__(self, members: list):
        self.members = list(members)
        if len({m.name for m in self.members}) != len(self.members):
            raise ValueError("duplicate blendstock names")
        tab = _table_smiles()
        tab_ik = _table_ik()
        for m in self.members:                       # a pure member's own table name wins for its InChIKey
            if not m.represented:
                tab_ik[_ik(tab[m.smiles[0]])] = m.smiles[0]
        ids, pos = [], {}
        cols = []
        for m in self.members:
            col = {}
            for s, wi in zip(m.smiles, m.w):
                ik = _ik(tab.get(s, s))
                if ik not in pos:
                    pos[ik] = len(ids)
                    ids.append(tab_ik.get(ik, s))
                col[pos[ik]] = col.get(pos[ik], 0.0) + float(wi)
            cols.append(col)
        self.comp_ids = ids
        self.A = np.zeros((len(ids), len(self.members)))
        for j, col in enumerate(cols):
            for i, v in col.items():
                self.A[i, j] = v
        self.names = [m.name for m in self.members]

    def expand(self, u) -> np.ndarray:
        u = np.clip(np.asarray(u, float), 0.0, None)
        u = u / u.sum()
        return self.A @ u

    def represented_mass(self, u) -> float:
        u = np.clip(np.asarray(u, float), 0.0, None)
        u = u / u.sum()
        return float(sum(ui for ui, m in zip(u, self.members) if m.represented))

    def synthetic_mass(self, u) -> float:
        u = np.clip(np.asarray(u, float), 0.0, None)
        u = u / u.sum()
        return float(sum(ui for ui, m in zip(u, self.members) if m.synthetic))


# ------------------------------------------------------------------------------------------------ D86
def tb_C(smiles: str) -> float:
    """Normal boiling point (C) from the default flash-registry Antoine (kPa, K form); NaN when absent."""
    r = _lookup(smiles, None)
    if r is None:
        return float("nan")
    A, B, C = r[0]
    t = B / (A - math.log10(P_ATM_KPA)) - C
    return float(t - 273.15) if 250.0 <= t <= 800.0 else float("nan")


def d86_quantities(tbs, vol) -> dict:
    """Declared D86 model -> dict T10, T50, T90, FBP, T50_T10, T90_T10, tbp_spread (TBP T90 - T10), monotone."""
    tbs, vol = np.asarray(tbs, float), np.asarray(vol, float)
    keep = vol > 0
    d = predict_d86_from_composition(tbs[keep], list(vol[keep]))
    fbp = predict_tbp_cuts(tbs[keep], vol[keep], cut_percents=(FBP_CUT_PCT,))[f"T{int(FBP_CUT_PCT)}_TBP_C"]
    q = dict(T10=d["T10_D86_C"], T50=d["T50_D86_C"], T90=d["T90_D86_C"], FBP=float(fbp))
    q["T50_T10"] = q["T50"] - q["T10"]
    q["T90_T10"] = q["T90"] - q["T10"]
    q["tbp_spread"] = d["T90_TBP_C"] - d["T10_TBP_C"]
    q["monotone"] = bool(q["T10"] <= q["T50"] <= q["T90"] <= q["FBP"])
    return q


def load_d86_margins(path: Path = D86_MARGINS) -> dict:
    """{'margins': {q: margin or 'inf'}, 'domain': {'tbp_spread_min', 'tbp_spread_max'}}; absent file -> every margin inf, empty domain."""
    if not Path(path).exists():
        return dict(margins={q: "inf" for q in D86_PROPS}, domain=dict(tbp_spread_min=float("nan"), tbp_spread_max=float("nan")))
    return json.load(open(path, encoding="utf-8"))


def domain_distance(q: dict, dom: dict) -> float:
    """>= 0; 0 inside the declared distillation domain.  Units: C / 100 (optimiser penalty; the judge only tests > 0)."""
    lo, hi = float(dom.get("tbp_spread_min", float("nan"))), float(dom.get("tbp_spread_max", float("nan")))
    if not (math.isfinite(lo) and math.isfinite(hi)):
        return 1.0                                                  # no fitted domain: never inside
    s = q["tbp_spread"]
    d = max(lo - s, 0.0) + max(s - hi, 0.0)
    d += max(q["T10"] - q["T50"], 0.0) + max(q["T50"] - q["T90"], 0.0) + max(q["T90"] - q["FBP"], 0.0)
    return float(d) / 100.0


def add_d86_models(models: dict, d86: Optional[dict] = None) -> dict:
    """Copy of route ``models`` (compound pool) with the four D86 PropertyModels; margins as recipe ``margin_fn`` (None outside the domain
    or when the R-W7-2 margin is inf); optimiser-only graded stand-in = the margin (0 when inf) plus the domain distance penalty."""
    d86 = load_d86_margins() if d86 is None else d86
    L = models["library"]
    cache: dict = {}

    def quantities(pool, w):
        key = (tuple(pool), np.asarray(w, float).tobytes())
        q = cache.get(key)
        if q is None:
            if len(cache) > 4096:
                cache.clear()
            w = np.asarray(w, float)
            smi = [L.smiles[n] for n in pool]
            tbs = np.array([tb_C(s) for s in smi])
            if not np.all(np.isfinite(tbs[w > 0])):
                raise E.ModelUnavailable("no Antoine Tb for a recipe member")
            rho = np.array([_rho15(L, n) for n in pool])
            if not np.all(np.isfinite(rho[w > 0])):
                raise E.ModelUnavailable("no grid rho at 288.15 K for a recipe member")
            vol = np.where(w > 0, w / np.where(np.isfinite(rho), rho, 1.0), 0.0)
            q = cache[key] = d86_quantities(tbs, vol / vol.sum())
        return q

    props = dict(models["props"])
    for qn, (canon, side) in D86_PROPS.items():
        raw = d86["margins"].get(qn, "inf")
        m_fit = None if (isinstance(raw, str) or not math.isfinite(float(raw))) else float(raw)

        def fn(pool, w, T_C, _q=qn):
            return float(quantities(pool, w)[_q])

        def compute(pool, w, T_C, _m=m_fit):
            q = quantities(pool, w)
            if domain_distance(q, d86["domain"]) > 0:
                return None, (f"{DOMAIN_REASON} (TBP T90-T10 {q['tbp_spread']:.1f} C, domain {d86['domain'].get('tbp_spread_min')}-"
                              f"{d86['domain'].get('tbp_spread_max')}; D86 monotone {q['monotone']})")
            if _m is None:
                return None, "no valid D86 margin (R-W7-2 gate failed: margin inf)"
            return _m, None

        def pen(pool, w, T_C):
            return domain_distance(quantities(pool, w), d86["domain"])

        props[canon] = E.PropertyModel(
            name=f"antoine_tb_d86:{qn}", fn=fn, margin_fn=LR._RecipeMargin(compute), margin_ln=False,
            margin_direction="under" if side == "hi" else "over", margin_source=f"delta_ws7_design_realism/d86_margins.json:{qn}",
            margin_stand_in=float(m_fit) if m_fit is not None else 0.0, margin_penalty_fn=pen,
            notes="Tb from the flash-registry Antoine (extrapolated to Tb); Engler offsets of baselines/distillation.py")
    out = dict(models)
    out["props"] = props
    return out


def _rho15(L, name) -> float:
    j = L.reg.ik2i.get(L.ik(name))
    if j is None or not L.provider.valid[j]:
        return float("nan")
    _, rho = L.provider.lookup(np.array([[j]]), T15_K)
    return float(rho[0][0])


# ------------------------------------------------------------------------------------------------ margin fit (R-W7-2)
def conformal_margin(r, alpha: float = 0.10) -> float:
    """ceil((n+1)(1-alpha))-th smallest residual (unsafe-side residuals), floored at 0; inf when n is too small for that rank."""
    r = np.sort(np.asarray([x for x in r if np.isfinite(x)], float))
    n = len(r)
    k = int(math.ceil((n + 1) * (1 - alpha)))
    if n == 0 or k > n:
        return float("inf")
    return float(max(r[k - 1], 0.0))


def fit_d86_margins(df: pd.DataFrame, alpha: float = 0.10, min_total: int = 30, min_per_source: int = 10) -> dict:
    """df: one row per fuel with columns source, pred_<q>, meas_<q> for q in T10, T50, T90, FBP, T50_T10, T90_T10.

    Unsafe-side residual: meas - pred for maximum limits (T10, FBP) and pred - meas for minimum limits (T50_T10, T90_T10); T50/T90 both sides.
    Certifying margin = max(pooled, per-source margins of sources with >= ``min_per_source`` fuels).  Gate: >= ``min_total`` fuels and every
    one of the sources AV-33-22 ('av3322') and AJFTD ('ajftd') with >= ``min_per_source`` fuels, else inf."""
    out = {}
    sides = {"T10": "under", "FBP": "under", "T50_T10": "over", "T90_T10": "over", "T50": "under", "T90": "under"}
    for q, dirn in sides.items():
        d = df[np.isfinite(df[f"pred_{q}"]) & np.isfinite(df[f"meas_{q}"])]
        r = (d[f"meas_{q}"] - d[f"pred_{q}"]) if dirn == "under" else (d[f"pred_{q}"] - d[f"meas_{q}"])
        n_src = d.source.value_counts().to_dict()
        pooled = conformal_margin(r, alpha)
        per = {s: conformal_margin(r[d.source == s], alpha) for s, k in n_src.items() if k >= min_per_source}
        gate = bool(len(d) >= min_total and all(n_src.get(s, 0) >= min_per_source for s in ("av3322", "ajftd")))
        cert = max([pooled] + list(per.values())) if gate else float("inf")
        err = d[f"pred_{q}"] - d[f"meas_{q}"]
        const = d[f"meas_{q}"] - d[f"meas_{q}"].mean()
        out[q] = dict(direction=dirn, n=int(len(d)), n_by_source={k: int(v) for k, v in n_src.items()}, gate_passed=gate,
                      pooled=pooled, per_source=per, certifying_margin=cert, bias=float(err.mean()) if len(d) else float("nan"),
                      mae=float(err.abs().mean()) if len(d) else float("nan"),
                      mae_constant=float(const.abs().mean()) if len(d) else float("nan"))
    return out


# ------------------------------------------------------------------------------------------------ uncertainty penalty
def alternate_spread_fn(models: dict, spec: dict) -> Callable:
    """U(pool, w) = max over spec constraints whose property has alternates of |alternate - primary| / spec scale (declared, optimiser only)."""
    items = []
    for key, c in spec["constraints"].items():
        prop, T_C = E.canonical_property(key)
        alts = models.get("alternates", {}).get(prop, [])
        pm = models["props"].get(prop)
        if pm is None or not alts:
            continue
        items.append((pm, alts, T_C, E._scale(c.get("min"), c.get("max"))))
    memo: dict = {}

    def U(pool, w):
        k = np.asarray(w, float).tobytes()
        if k in memo:
            return memo[k]
        if len(memo) > 4096:
            memo.clear()
        best = 0.0
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            for pm, alts, T_C, sc in items:
                try:
                    p = float(pm.fn(pool, w, T_C))
                except Exception:
                    continue
                for a in alts:
                    try:
                        v = float(a.fn(pool, w, T_C))
                    except Exception:
                        continue
                    if math.isfinite(p) and math.isfinite(v):
                        best = max(best, abs(v - p) / sc)
        memo[k] = best
        return best
    return U


# ------------------------------------------------------------------------------------------------ optimiser
def optimise(vector: Callable, n: int, n_starts: int = 10, seed: int = 42, hard: Optional[Callable] = None, n_hard_tail: int = 0,
             penalty: Optional[Callable] = None, lam: float = 0.0, max_iter: int = 120) -> dict:
    """Maximise min(vector(u)) - lam * penalty(u) over the simplex in ``u``.

    ``vector(u)`` -> engine slacks (its last ``n_hard_tail`` entries, the extrapolation sides, are ALSO hard constraints >= 0);
    ``hard(u)`` -> array that must be >= 0 (synthetic cap, distillation domain).  Returns dict(u, score, feasible, n_converged, n_starts)."""
    rng = np.random.default_rng(seed)
    starts = [np.full(n, 1.0 / n)] + [rng.dirichlet(np.ones(n)) for _ in range(max(n_starts - 1, 0))]

    def norm(z):
        u = np.clip(z, 0.0, None)
        s = u.sum()
        return u / s if s > 0 else np.full(n, 1.0 / n)

    def obj(u):
        v = vector(u)
        return float(np.min(v)) - (lam * float(penalty(u)) if (penalty is not None and lam) else 0.0)

    def hard_vals(u):
        parts = []
        if n_hard_tail:
            parts.append(np.asarray(vector(u), float)[-n_hard_tail:])
        if hard is not None:
            parts.append(np.asarray(hard(u), float))
        return np.concatenate(parts) if parts else np.zeros(0)

    def feasible(u):
        h = hard_vals(u)
        return bool(h.size == 0 or np.all(h >= -1e-9))

    best = dict(u=starts[0], score=-np.inf, feasible=False)
    n_ok = 0
    for u0 in starts:
        z0 = np.append(u0, obj(u0))
        cons = [{"type": "eq", "fun": lambda z: np.sum(z[:-1]) - 1.0},
                {"type": "ineq", "fun": lambda z: np.asarray(vector(norm(z[:-1])), float)
                 - (lam * float(penalty(norm(z[:-1]))) if (penalty is not None and lam) else 0.0) - z[-1]}]
        if n_hard_tail or hard is not None:
            cons.append({"type": "ineq", "fun": lambda z: hard_vals(norm(z[:-1]))})
        try:
            res = minimize(lambda z: -z[-1], z0, method="SLSQP", bounds=[(0, 1)] * n + [(None, None)], constraints=cons,
                           options={"maxiter": max_iter, "ftol": 1e-9})
            u = norm(res.x[:-1])
            n_ok += int(res.success)
        except Exception:
            u = u0
        for cand in (u, u0):
            f, s = feasible(cand), obj(cand)
            if (f and not best["feasible"]) or (f == best["feasible"] and s > best["score"]):
                best = dict(u=cand, score=s, feasible=f)
    best.update(n_converged=n_ok, n_starts=len(starts))
    return best


# ------------------------------------------------------------------------------------------------ judge
def judge(spec_name: str, models: dict, w_comp: np.ndarray, margins_path, synthetic_blend: bool, alt_limits: bool, recipe_margins: bool,
          d86: bool, represented_mass: float = 0.0) -> E.DesignResult:
    """``library_route.evaluate_recipe`` with the d86 flag (same body), then the represented-blendstock rule."""
    pool = list(models["pool"])
    blocked = E.blocked_result(spec_name, pool, models, synthetic_blend=synthetic_blend)
    if blocked is not None:
        return blocked
    w = np.asarray(w_comp, float)
    w = w / w.sum()
    spec = E.load_spec(spec_name, synthetic_blend=synthetic_blend, alt_limits=alt_limits, d86=d86)
    margins = E.load_margins(margins_path)
    use = LR.with_recipe_margins(models, margins) if recipe_margins else models
    cons = E.evaluate_blend(spec, pool, w, use, margins)
    res = LR.annotate(E.DesignResult(weights=dict(zip(pool, map(float, w))), predictions={c.key: c.prediction for c in cons},
                                     margins={c.key: c.applied_margin for c in cons}, passed=False, extrapolating=False,
                                     binding_constraints=[], constraints=cons, spec_name=spec_name, pool=pool), margins_path,
                      models=models, alt_limits=alt_limits, recipe_margins=recipe_margins)
    apply_represented_rule(res, represented_mass, models)
    return res


def apply_represented_rule(res: E.DesignResult, represented_mass: float, models: Optional[dict] = None) -> None:
    """Declared rule: with represented blendstock mass > 0 every model-predicted (not by-construction) constraint is at most ``uncertified``."""
    if represented_mass <= 0:
        return
    props = (models or {}).get("props", {})
    for c in res.constraints:
        pm = props.get(c.prop)
        if c.status in (E.STATUS_CONSTRUCTION, E.STATUS_FAIL) or (pm is not None and pm.by_construction):
            continue
        if c.status != E.STATUS_UNCERTIFIED:
            c.status = E.STATUS_UNCERTIFIED
        if REPRESENTED_REASON not in c.reason:
            c.reason = (c.reason + "; " if c.reason else "") + REPRESENTED_REASON
    res.passed = bool(all(c.status in E.OK_STATUSES for c in res.constraints) and not res.extrapolating)


def design_pool(spec_name: str, bp: BlendPool, models: dict, margins_path, n_starts: int = 10, seed: int = 42, synthetic_blend: bool = True,
                alt_limits: bool = True, recipe_margins: bool = True, d86: bool = True, lam: float = LAMBDA_UNC,
                synthetic_cap: Optional[float] = None, max_iter: int = 120) -> tuple:
    """-> (DesignResult on the compound pool, top-level weights dict, optimiser info)."""
    pool = list(models["pool"])
    blocked = E.blocked_result(spec_name, pool, models, seed, synthetic_blend)
    if blocked is not None:
        return blocked, {n: float("nan") for n in bp.names}, dict(blocked=True)
    spec = E.load_spec(spec_name, synthetic_blend=synthetic_blend, alt_limits=alt_limits, d86=d86)
    margins = E.load_margins(margins_path)
    use = LR.with_recipe_margins(models, margins) if recipe_margins else models
    sides = E._Sides(spec, pool, use, margins, graded=bool(recipe_margins))
    U = alternate_spread_fn(models, spec) if lam else None
    dom_pms = [use["props"][c] for c, _ in D86_PROPS.values() if d86 and c in use["props"]]

    def vector(u):
        return sides.vector(bp.expand(u))

    def hard(u):
        h = []
        if synthetic_cap is not None:
            h.append(synthetic_cap - bp.synthetic_mass(u))
        if dom_pms:
            try:
                h.append(-float(dom_pms[0].margin_penalty_fn(pool, bp.expand(u), None)))
            except Exception:
                h.append(-1.0)
        return np.array(h, float)

    has_hard = synthetic_cap is not None or bool(dom_pms)
    info = optimise(vector, len(bp.names), n_starts=n_starts, seed=seed, hard=hard if has_hard else None,
                    n_hard_tail=len(sides.extrap_sets), penalty=(lambda u: U(pool, bp.expand(u))) if U else None, lam=lam, max_iter=max_iter)
    u = np.where(info["u"] < 1e-6, 0.0, info["u"])
    u = u / u.sum()
    res = judge(spec_name, models, bp.expand(u), margins_path, synthetic_blend, alt_limits, recipe_margins, d86, bp.represented_mass(u))
    res.optimiser = dict(res.optimiser, best_objective=float(info["score"]), hard_feasible=bool(info["feasible"]),
                         n_converged=int(info["n_converged"]), n_starts=int(info["n_starts"]), seed=int(seed), lam=float(lam),
                         min_slack=float(np.nanmin([c.slack for c in res.constraints])) if res.constraints else float("nan"),
                         uncertainty_U=float(U(pool, bp.expand(u))) if U else float("nan"))
    return res, dict(zip(bp.names, map(float, u))), info


def build_models(bp: BlendPool) -> dict:
    """Route models on the compound pool, in the BlendPool's compound order (``models['pool']`` is the resolved name list).

    Own grid and flash-registry data are required; Tm/dHfus are NOT (many GCxGC representatives, e.g. iso-paraffins, have none): a member
    without them makes freezing ``model unavailable`` (uncertified) instead of refusing the pool.  ``models['sle_gaps']`` lists them."""
    models = LR.build_route_models(bp.comp_ids, require=("grid", "flash"))
    models["sle_gaps"] = sorted(n for n, g in models["library"].own_data_gaps(("sle",)).items())
    if len(models["pool"]) != len(bp.comp_ids):
        raise ValueError("route pool size differs from the blend pool compound list")
    return models
