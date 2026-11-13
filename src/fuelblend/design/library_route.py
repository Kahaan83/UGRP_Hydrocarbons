"""v5 Task 6: inverse design on a SMILES / InChIKey pool, with an honest status for every constant.

``RouteLibrary`` is a ``models_v4.Library`` whose pool may be given as table names, SMILES or InChIKeys.  A SMILES or
InChIKey that matches a compound of the verified constants table uses the table name and the verified constants (identical
design to the name-based pool).  A compound that is not in the table gets route constants, each with a status:

    MW                  RDKit                                 verified
    rho15               pure-property grid at 288.15 K        route_grid
    Tm, dHfus           ``build_tmap`` (Coutinho table)       coutinho_table
    flash               anchored Antoine / LFL registry       route_anchored
    freezing            the API function: ideal SvL in mole fraction on the Coutinho tmap for every member
                        (so every member's Tm/dHfus is recorded as coutinho_table and freezing is uncertified)
    RON, MON, DCN, Psat37_8   none                            absent  (ModelUnavailable)

``v()`` returns route values too and records every non-verified one in ``nonverified_used``.  ``build_route_models`` turns
that record into ``margin_missing_reason`` so the unchanged engine reports the constraint as ``uncertified``.
Flash is computed on the LFL overlay registry (U1), the same registry the API's designated flash uses.  Disclosure: kinematic
viscosity of a non-table member comes from the pure-property grid, which is fitted and unpinned; the same grid rho is flagged
``route_grid`` in density but is used silently inside viscosity (no verified viscosity constant exists, so it stays certifiable).
No representatives are used inside a pool: a member without its own data is refused by name.  ``ensure()`` is never called.
"""
from __future__ import annotations

import dataclasses
import math
import re
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, inchi as rdinchi

from fuelblend.baselines import flash_unifac as fu
from fuelblend.design import engine_v4 as E
from fuelblend.design import models_v4 as M
from fuelblend.design.engine_v4 import ModelUnavailable

ROOT = Path(__file__).resolve().parents[3]
PURE_PROPS_V5 = ROOT / "data" / "processed" / "v5" / "pure_props_v5.csv"
_IK_RE = re.compile(r"^[A-Z]{14}-[A-Z]{10}-[A-Z]$")
T15_K = 288.15
ROUTE_PROPS = {"rho15": "route_grid", "Tm": "coutinho_table", "dHfus": "coutinho_table", "flash": "route_anchored"}


class RouteLibrary(M.Library):
    def __init__(self, pool_ids, pure_props=PURE_PROPS_V5):
        super().__init__(pure_props)
        from fuelblend.api import predict_route as pr
        self._pr = pr
        self._ctx = pr.get_ctx()
        self.provider = pr.ReadOnlyProvider(self.reg, pr.GRID_NPZ)     # frozen grid: never computes, never writes
        self._table_names = set(self.smiles)
        self._table_ik = {self.ik(n): n for n in self._table_names}
        self._route: dict = {}                    # (name, prop) -> (value, status) for compounds outside the table
        self.nonverified_used: set = set()
        self.pool = [self._resolve(p) for p in pool_ids]
        if len(set(self.pool)) != len(self.pool):
            raise ValueError("duplicate compounds in the pool")

    # ------------------------------------------------------------------ pool id -> name
    def _resolve(self, pid: str) -> str:
        pid = str(pid).strip()
        if pid in self._table_names:
            return pid
        reg = self._ctx.reg
        if _IK_RE.match(pid):
            ik = pid
            if ik in self._table_ik:
                return self._table_ik[ik]
            i = reg.ik2i.get(ik)
            if i is None:
                raise ModelUnavailable(f"{pid}: InChIKey not in the compound registry, so no own data")
            smi = reg.smiles[i]
        else:
            mol = Chem.MolFromSmiles(pid)
            if mol is None:
                raise ValueError(f"unparseable pool id: {pid!r}")
            ik = rdinchi.MolToInchiKey(mol)
            if ik in self._table_ik:
                return self._table_ik[ik]
            smi = Chem.MolToSmiles(mol)
        name = Chem.MolToSmiles(Chem.MolFromSmiles(smi))
        if name not in self.smiles:
            self.smiles[name] = smi
            self._ik[name] = ik
            self._build_route_constants(name, ik, smi)
        return name

    def _build_route_constants(self, name, ik, smi):
        mol = Chem.MolFromSmiles(smi)
        self._route[(name, "MW")] = (float(Descriptors.MolWt(mol)), "verified")
        j = self.reg.ik2i.get(ik)
        if j is not None and self.provider.valid[j]:
            _, rho = self.provider.lookup(np.array([[j]]), T15_K)
            if np.isfinite(rho[0][0]):
                self._route[(name, "rho15")] = (float(rho[0][0]), "route_grid")
        t = self._pr._tmap().get(ik) or {}
        tm, dh = t.get("tm_K", float("nan")), t.get("dh_fus_J_mol", float("nan"))
        if np.isfinite(tm) and np.isfinite(dh):
            self._route[(name, "Tm")] = (float(tm) - 273.15, "coutinho_table")
            self._route[(name, "dHfus")] = (float(dh) / 1000.0, "coutinho_table")
        if fu._lookup(smi, None) is not None:
            self._route[(name, "flash")] = (float("nan"), "route_anchored")

    # ------------------------------------------------------------------ constants
    def is_table(self, name) -> bool:
        return name in self._table_names

    def status(self, name, prop) -> str:
        if self.is_table(name):
            return self._tab.get((name, prop), (float("nan"), "absent"))[1]
        return self._route.get((name, prop), (float("nan"), "absent"))[1]

    def touch(self, name, prop):
        """Record the use of a route constant that the model reads without calling ``v`` (the flash registry)."""
        if not self.is_table(name):
            st = self.status(name, prop)
            if st == "absent":
                raise ModelUnavailable(f"{name}.{prop} status=absent")
            if st != "verified":
                self.nonverified_used.add(f"{name}.{prop}={st}")

    def v(self, name, prop):
        if self.is_table(name):
            return super().v(name, prop)
        val, st = self._route.get((name, prop), (float("nan"), "absent"))
        if st == "absent" or not math.isfinite(float(val)):
            raise ModelUnavailable(f"{name}.{prop} status=absent")
        if st != "verified":
            self.nonverified_used.add(f"{name}.{prop}={st}")
        return float(val)

    def own_data_gaps(self, require=("grid", "flash", "sle")) -> dict:
        """{member: [missing data kinds]} for members that lack own data (no representatives are used in a pool)."""
        gaps = {}
        for n in self.pool:
            j = self.idx(n)
            miss = []
            if "grid" in require and not (j is not None and self.provider.valid[j]):
                miss.append("grid (viscosity/density)")
            if "flash" in require and fu._lookup(self.smiles[n], None) is None:
                miss.append("flash registry")
            if "sle" in require:
                ok = all(np.isfinite(float(self._tab.get((n, p), (float("nan"),))[0])) for p in ("Tm", "dHfus")) \
                    if self.is_table(n) else ((n, "Tm") in self._route)
                if not ok:
                    miss.append("Tm/dHfus")
            if miss:
                gaps[n] = miss
        return gaps


# ------------------------------------------------------------------------------------------------ models
def _probe(L: RouteLibrary, pm: E.PropertyModel) -> set:
    L.nonverified_used.clear()
    n = len(L.pool)
    try:
        pm.fn(list(L.pool), np.full(n, 1.0 / n), -20.0)
    except ModelUnavailable:
        pass
    except Exception as e:                      # never leave a model certifiable when its probe could not run
        used = set(L.nonverified_used)
        L.nonverified_used.clear()
        return used, f"probe failed: {type(e).__name__}: {e}"
    used = set(L.nonverified_used)
    L.nonverified_used.clear()
    return used, ""


def _flag(L: RouteLibrary, pm: E.PropertyModel) -> E.PropertyModel:
    if pm.name.startswith("le_chatelier"):                         # flash on the LFL overlay registry (U1), like the API
        unifac = "unifac" in pm.name

        def fn(pool, w, T_C, _u=unifac):
            for n in pool:
                L.touch(n, "flash")
            reg = L._pr._overlay_ctx().flash_reg
            smi = [L.smiles[n] for n in pool]
            if any(fu._lookup(x, reg) is None for x in smi):
                raise ModelUnavailable("no Antoine/LFL entry in the flash registry")
            x = M._mole_frac(L, pool, w)
            with np.errstate(all="ignore"), warnings.catch_warnings():
                warnings.simplefilter("ignore")
                d = fu.flash_point_unifac_detail(smi, x, registry=reg, use_unifac=_u)
            if d["status"] != "ok" or not np.isfinite(d["flash_C"]):
                raise ModelUnavailable(f"flash solver status {d['status']}")
            return float(d["flash_C"])
        pm = dataclasses.replace(pm, fn=fn)
    used, err = _probe(L, pm)
    if not used and not err:
        return pm
    why = err if err else "non-verified constant(s): " + ", ".join(sorted(used))
    if err:                                      # the failing model must report unavailable, not raise
        inner = pm.fn

        def fn(pool, w, T_C, _inner=inner, _why=err):
            raise ModelUnavailable(_why)
        pm = dataclasses.replace(pm, fn=fn)
    return dataclasses.replace(pm, margin_keys=(), margin_keys_by_side=None, margin_keys_fn=None,
                               margin_missing_reason=why, by_construction=False, T_valid_C=None)


def _freezing_api_fn(L: RouteLibrary):
    """Design freezing = the function ``predict_fuel`` uses: ideal Schroeder-van Laar in mole fraction (RDKit MW) on the
    Coutinho ``build_tmap`` constants, for every pool member (D7, review C1).  The constants are route constants, so each
    member is recorded as ``coutinho_table`` and the model is uncertified (freezing has no valid margin anyway)."""
    from fuelblend.eval.composition_route import ideal_sle_C

    def fn(pool, w, T_C):
        pool = list(pool)
        tmap = L._pr._tmap()
        iks = [L.ik(n) for n in pool]
        for n, ik in zip(pool, iks):
            if not tmap.get(ik):
                raise ModelUnavailable(f"{n}.Tm status=absent (no Coutinho table entry)")
            L.nonverified_used.add(f"{n}.Tm={'coutinho_table'}")
            L.nonverified_used.add(f"{n}.dHfus={'coutinho_table'}")
        w = np.asarray(w, float)
        mw = np.array([L._ctx.mw(L.smiles[n]) for n in pool])
        x = (w / w.sum()) / mw
        x = x / x.sum()
        v = ideal_sle_C(x, iks, tmap)
        if not np.isfinite(v):
            raise ModelUnavailable("no crystallising component with Tm and dHfus")
        return float(v)
    return fn


def build_route_models(pool_ids, pure_props=PURE_PROPS_V5, require=("grid", "flash", "sle")) -> dict:
    """Models for ``engine_v4.design`` on a SMILES / InChIKey / name pool.  The resolved pool is ``models['pool']``."""
    L = RouteLibrary(pool_ids, pure_props)
    gaps = L.own_data_gaps(require)
    if gaps:
        n, m = next(iter(gaps.items()))
        raise ModelUnavailable(f"{n}: no own data for {', '.join(m)}; no representatives are used inside a pool "
                               f"(all gaps: {gaps})")
    models = M.build_models(L)
    props, alts = models["props"], models["alternates"]
    mole = next(p for p in alts["freezing"] if p.name == "ideal_sle_mole_fraction")
    alts["freezing"] = [p for p in alts["freezing"] if p.name != "ideal_sle_mole_fraction"] + [props["freezing"]]
    props["freezing"] = dataclasses.replace(mole, fn=_freezing_api_fn(L), certified_kinds={"jet"})
    models["props"] = {k: _flag(L, p) for k, p in props.items()}
    models["alternates"] = {k: [_flag(L, p) for p in v] for k, v in alts.items()}
    models["pool"] = list(L.pool)
    models["library"] = L
    return models


LAB_MARGIN_NOTE = "lab_margin_set_by_one_system"
SMOKE_MARGIN_KEY = "smoke_point_mm.lab_mixtures.bulk_regression"


def default_n_starts(pool) -> int:
    return 8 if len(pool) > 20 else 30


SMOKE_MARGIN_REASON = "smoke margin 7.8 mm from binaries and surrogates (user decision), applied inside the regression domain only"


def _rejudge_smoke(result: E.DesignResult, margins: dict) -> None:
    """Smoke margin, design side, DOMAIN-GATED (user decision 2026-11-13 plus controller ruling).

    The engine (unchanged) judges smoke without a margin and tags a recipe outside the regression domain (aromatics 8-23 vol %, result
    19-25 mm) with the OutOfDomainWarning reason, which keeps it ``uncertified`` with margin None.  The optimiser therefore also sees no
    smoke margin (a margin is recipe-dependent through the domain, the engine's margin is not).  Here, after the engine, a smoke constraint
    that is IN domain is re-judged with the v5 margin: slack = (pred - margin - limit) / scale; ``fail`` if negative.  Out of domain it is left
    exactly as the engine produced it."""
    e = margins.get(SMOKE_MARGIN_KEY)
    m, unit, direction, _ = E.read_applied_margin(e)
    if m is None:
        return
    for c in result.constraints:
        if c.prop != "smoke_point" or not math.isfinite(c.prediction) or "OutOfDomainWarning" in c.reason:
            continue
        if c.lo is None or c.hi is not None:                       # smoke is a minimum limit: over-prediction is the unsafe side
            continue
        sc = E._scale(c.lo, c.hi)
        c.applied_margin, c.margin_keys = m, (SMOKE_MARGIN_KEY,)
        c.slack = (c.prediction - m - c.lo) / sc
        if c.slack < 0 or c.prediction < c.lo:
            c.status = E.STATUS_FAIL
        c.reason = (c.reason + "; " if c.reason else "") + SMOKE_MARGIN_REASON + ("; limit violated after applying the margin" if c.status == E.STATUS_FAIL else "")
        result.margins[c.key] = m
    result.passed = result.passed and all(c.status in E.OK_STATUSES for c in result.constraints)


def annotate(result: E.DesignResult, margins_path=None) -> E.DesignResult:
    """Every flash constraint carrying a finite margin gets ``lab_margin_set_by_one_system`` in its reason (BASE Task 10:
    the lab margin is the group score of the single mesitylene system).  In-domain smoke constraints get the v5 smoke margin."""
    for c in result.constraints:
        has_margin = c.applied_margin is not None and math.isfinite(c.applied_margin)
        if c.prop == "flash" and has_margin and LAB_MARGIN_NOTE not in c.reason:
            c.reason = (c.reason + "; " if c.reason else "") + LAB_MARGIN_NOTE
    _rejudge_smoke(result, E.load_margins(margins_path))
    return result


def design_routed(spec_name, models, n_starts=None, seed=42, margins_path=None) -> E.DesignResult:
    pool = models["pool"]
    return annotate(E.design(spec_name, pool, models, n_starts=n_starts or default_n_starts(pool), seed=seed,
                             margins_path=margins_path), margins_path)


# ------------------------------------------------------------------------------------------------ forward check
@dataclasses.dataclass
class _Row:
    constraint: str
    design_pred: float
    forward_pred: float
    tol: float
    cause: str = ""
    informational: bool = False
    limit_key: str = ""                   # constraint key whose spec limit this row is judged against ("" = no limit judgement)


CHECK_COLUMNS = ["constraint", "design_pred", "forward_pred", "abs_diff", "tol", "pass", "cause",
                 "limit", "lo", "hi", "design_meets", "api_meets"]


def _name_smiles() -> dict:
    d = pd.read_csv(PURE_PROPS_V5).drop_duplicates("compound")
    return dict(zip(d.compound, d.smiles))


def evaluate_recipe(spec_name: str, models: dict, weights: dict, margins_path=None) -> E.DesignResult:
    """A ``DesignResult`` for a given recipe (no optimisation), scored by the engine with ``models``."""
    pool = list(models["pool"]) if "pool" in models else list(weights)
    w = np.array([weights.get(n, 0.0) for n in pool], float)
    w = w / w.sum()
    spec = E.load_spec(spec_name)
    cons = E.evaluate_blend(spec, pool, w, models, E.load_margins(margins_path))
    return annotate(E.DesignResult(weights=dict(zip(pool, map(float, w))), predictions={c.key: c.prediction for c in cons},
                          margins={c.key: c.applied_margin for c in cons}, passed=False, extrapolating=False,
                          binding_constraints=[], constraints=cons, spec_name=spec_name, pool=pool), margins_path)


def _fnum(x) -> float:
    try:
        return float("nan") if x is None else float(x)
    except (TypeError, ValueError):
        return float("nan")


def _meets(v, lo, hi) -> bool:
    """Point-estimate verdict against the spec limit(s); a missing or non-finite value never meets."""
    v = _fnum(v)
    if not math.isfinite(v):
        return False
    return not ((lo is not None and v < lo) or (hi is not None and v > hi))


def _limit_text(lo, hi) -> str:
    if lo is not None and hi is not None:
        return f"{lo:g}..{hi:g}"
    return f">= {lo:g}" if lo is not None else (f"<= {hi:g}" if hi is not None else "")


def _api_value(f: dict, key: str) -> tuple:
    """(value or nan, cause suffix) of an API property dict; a None value gives the status/reason as the cause."""
    v = f.get(key)
    if v is None:
        return float("nan"), f"API returned no property {key}"
    val = v.get("value")
    if val is None or not math.isfinite(_fnum(val)):
        why = "; ".join(v.get("domain_reasons") or []) or str(v.get("status", ""))
        return float("nan"), f"API value is None (status {v.get('status')}): {why}"
    return float(val), ""


def design_forward_check(result: E.DesignResult, T_C_visc: float | None = None) -> pd.DataFrame:
    """Re-predict ``result.weights`` through ``predict_fuel`` and compare with the design-side predictions.

    Two kinds of rows.  (1) For EVERY constrained property: ``design_pred`` (the engine's point estimate), ``forward_pred`` (the
    API value), the spec ``limit``, ``design_meets`` and ``api_meets`` (each point estimate against the limit, margins
    ignored).  A design meets its limits only if ``design_meets`` and ``api_meets`` are both true for every such row
    (``limits_met_by_both``); a forward value of None gives ``api_meets`` False with the cause.  (2) Like-for-like rows with a
    tolerance (R-API-4): flash 0.01 C against the API designated (overlay) flash, viscosity 1e-6 relative, net heat exact
    (1e-9) for the composition-H model.  ``pass`` is judged on those tolerance rows; the other rows are reported."""
    from fuelblend.api import predict_route as pr
    names = _name_smiles()
    pool = list(result.pool) or list(result.weights)
    smi = {n: names.get(n, n) for n in pool}
    comp = {smi[n]: float(result.weights[n]) for n in pool if result.weights[n] > 0}
    cons = {c.key: c for c in result.constraints}
    L = RouteLibrary(pool)
    kv = next((k for k in cons if k.startswith("kin_visc_")), None)
    if kv is not None:
        mo = re.match(r"^kin_visc_(m?)(\d+)C_cSt$", kv)
        T_spec = (-1.0 if mo.group(1) else 1.0) * float(mo.group(2))
        if T_C_visc is not None and abs(T_C_visc - T_spec) > 1e-9:
            raise ValueError(f"T_C_visc={T_C_visc} does not match the spec constraint {kv} ({T_spec} C)")
        T_C_visc = T_spec
    elif T_C_visc is None:
        T_C_visc = -20.0
    f_vis = pr.predict_fuel(comp, T_K=T_C_visc + 273.15, fuel_type="jet")
    f_15 = pr.predict_fuel(comp, T_K=T15_K, fuel_type="jet")
    rows: list[_Row] = []

    def dp(key):
        c = cons.get(key)
        return _fnum(c.prediction) if c is not None else float("nan")

    # (1) every constrained property against the API value and the spec limit
    api_key = {"density_15C_kg_m3": ("density_kg_m3", f_15), "flash_point_C": ("flash_point_C", f_15),
               "freezing_point_C": ("freezing_point_C", f_15), "nhoc_MJ_kg": ("heat_net_MJ_kg", f_15),
               "smoke_point_mm": ("smoke_point_mm", f_15), "aromatics_wt_pct": ("aromatics_wt_pct", f_15),
               "aromatics_vol_pct": ("aromatics_vol_pct", f_15)}
    if kv:
        api_key[kv] = ("viscosity_kin_cSt", f_vis)
    notes = {
        "density_15C_kg_m3": "design: ideal volume mixing of the table rho15; API: pure-property grid rho at 288.15 K",
        "flash_point_C": "design and API both use the LFL-overlay registry (U1)",
        "freezing_point_C": "design and API both use ideal Schroeder-van Laar in mole fraction on the Coutinho tmap constants",
        "nhoc_MJ_kg": "design primary model uses constant H = 13.8 wt%; the API uses composition H",
        "smoke_point_mm": "design and API both use the AJFTD bulk regression on aromatics vol % (API aromatics from the grid rho at 288.15 K)",
        "aromatics_wt_pct": "class mass fraction, by construction",
        "aromatics_vol_pct": "design: table rho15 volume fractions; API: grid rho at 288.15 K",
    }
    if kv:
        notes[kv] = "kinematic viscosity at the spec temperature, Refutas on the pure-property grid"
    for key in cons:
        if key not in api_key:
            continue
        ak, f = api_key[key]
        fv, why = _api_value(f, ak)
        rows.append(_Row(key, dp(key), fv, float("nan"), notes.get(key, "") + ("; " + why if why else ""),
                         informational=True, limit_key=key))

    # (2) like-for-like tolerance rows
    if "flash_point_C" in cons:
        fl = f_15["flash_point_C"]
        rows.append(_Row("flash_point_C[tol]", dp("flash_point_C"), _fnum(fl.get("value")), 0.01,
                         "design and API both use the LFL-overlay registry (U1)"))
        v4d = (fl.get("variants") or {}).get("v4_designated")
        rows.append(_Row("flash_point_C[v4_designated]", dp("flash_point_C"), _fnum(v4d), float("nan"),
                         "informational: plain registry (no LFL overlay); the diff is the overlay effect", True))
    if kv:
        rows.append(_Row(kv + "[tol]", dp(kv), _api_value(f_vis, "viscosity_kin_cSt")[0], 1e-6))
    if "nhoc_MJ_kg" in cons:
        # The API uses the composition-H group-additivity model; the design primary uses constant H.  Two like-for-like rows
        # use the design's composition-H model: as designed (table MW, 2 decimals) and with RDKit MW, the API's MW.
        w = np.array([result.weights[n] for n in pool], float)
        fv = _api_value(f_15, "heat_net_MJ_kg")[0]
        fn = M.make_fns(L)["nhoc_comp_h"]
        rows.append(_Row("nhoc_MJ_kg[compH_table_MW]", fn(pool, w, None), fv, float("nan"),
                         "same composition-H model; the design reads listed MW (2 decimals) from pure_props_v5, the API reads RDKit MW", True))
        Lr = RouteLibrary(pool)
        for n in pool:
            Lr._tab[(n, "MW")] = (float(Descriptors.MolWt(Chem.MolFromSmiles(smi[n]))), "verified")
        rows.append(_Row("nhoc_MJ_kg[compH_rdkit_MW]", M.make_fns(Lr)["nhoc_comp_h"](pool, w, None), fv, 1e-9,
                         "design-side composition-H group-additivity model with RDKit MW (alternate model; primary is constant H)"))
    out = []
    for r in rows:
        fin = math.isfinite(_fnum(r.forward_pred)) and math.isfinite(_fnum(r.design_pred))
        d = abs(r.design_pred - r.forward_pred) if fin else float("nan")
        if r.informational or not math.isfinite(r.tol):
            ok = "reported"
        elif r.constraint.startswith("kin_visc_"):
            ok = bool(fin and d <= r.tol * abs(r.forward_pred))
        else:
            ok = bool(fin and d <= r.tol)
        row = dict(constraint=r.constraint, design_pred=_fnum(r.design_pred), forward_pred=_fnum(r.forward_pred), abs_diff=d,
                   tol=r.tol, cause=r.cause, limit="", lo=np.nan, hi=np.nan, design_meets=None, api_meets=None)
        row["pass"] = ok
        if r.limit_key:
            c = cons[r.limit_key]
            row.update(limit=_limit_text(c.lo, c.hi), lo=_fnum(c.lo), hi=_fnum(c.hi),
                       design_meets=_meets(r.design_pred, c.lo, c.hi), api_meets=_meets(r.forward_pred, c.lo, c.hi))
        out.append(row)
    return pd.DataFrame(out, columns=CHECK_COLUMNS)


def limits_met_by_both(fc: pd.DataFrame) -> bool:
    """True iff every constrained property meets its limit at the design point estimate AND at the API value."""
    lim = fc[fc["limit"] != ""]
    return bool(len(lim) > 0 and lim["design_meets"].astype(bool).all() and lim["api_meets"].astype(bool).all())
