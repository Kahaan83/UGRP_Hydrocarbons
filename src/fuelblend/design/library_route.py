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

import copy
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
from fuelblend.design import smoke_margin_v5 as SM
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


FREEZING_MARGIN_KEY = "freezing_point_C.lab_liquidus.svl_verified"


def _freezing_api_fn(L: RouteLibrary):
    """Design freezing = the function ``predict_fuel`` uses: ideal Schroeder-van Laar in mole fraction (RDKit MW) on the designated freezing
    constants (``predict_route._tmap_freezing``: tm_dhfus_verified.csv where present, else the Coutinho entry).  A member whose constants come from the
    verified csv is not recorded as non-verified; a member on the Coutinho entry is recorded ``coutinho_table`` (freezing then stays uncertified)."""
    from fuelblend.eval.composition_route import VERIFIED_SOURCE, ideal_sle_C

    def fn(pool, w, T_C):
        pool = list(pool)
        tmap = L._pr._tmap_freezing()
        iks = [L.ik(n) for n in pool]
        for n, ik in zip(pool, iks):
            t = tmap.get(ik)
            if not t:
                raise ModelUnavailable(f"{n}.Tm status=absent (no verified or Coutinho table entry)")
            if t.get("source") != VERIFIED_SOURCE:
                L.nonverified_used.add(f"{n}.Tm={t.get('source', 'coutinho_table')}")
                L.nonverified_used.add(f"{n}.dHfus={t.get('source', 'coutinho_table')}")
        w = np.asarray(w, float)
        mw = np.array([L._ctx.mw(L.smiles[n]) for n in pool])
        x = (w / w.sum()) / mw
        x = x / x.sum()
        v = ideal_sle_C(x, iks, tmap)
        if not np.isfinite(v):
            raise ModelUnavailable("no crystallising component with Tm and dHfus")
        return float(v)
    return fn


VISC_ML_NAME = "refutas+ml_pair_gn"
VISC_ML_REASON = ("designated ML pair-residual viscosity (as predict_fuel) with its R-MLV-6 ln margin "
                  "(results/v5/BASE/delta_ml_viscosity/readings.json); judged after the optimiser")


def visc_ml_status(smiles, w, T_C: float) -> dict:
    """The API's ML viscosity decision for a known recipe: ``{'in_domain', 'g', 'margin', 'reason'}``.

    Same rule as ``predict_route._add_viscosity_ml`` / ``_promote_ml_viscosity`` for a known composition: members with weight > 0,
    mole fractions from RDKit MW (``ctx.mw``), ``viscosity_pair.domain_status`` plus the 233-393 K training range; ``g`` is the ln
    correction and ``margin`` the R-MLV-6 margin of the temperature bin.  Out of domain: ``g`` 0, ``margin`` None."""
    from fuelblend.api import predict_route as pr
    T_K = float(T_C) + 273.15
    w = np.asarray(w, float)
    keep = w > 0
    smi = [s for s, k in zip(smiles, keep) if k]
    out = dict(in_domain=False, g=0.0, margin=None, reason="")
    if not smi:
        out["reason"] = "empty recipe"
        return out
    try:
        from fuelblend.ml import viscosity_pair as vp
        ctx = pr.get_ctx()
        x = w[keep] / np.array([ctx.mw(s) for s in smi])
        x = x / x.sum()
        dom = vp.domain_status(smi, T_K)
        if dom["in_domain"] and not (pr.ML_VISC_T_RANGE_K[0] <= T_K <= pr.ML_VISC_T_RANGE_K[1]):
            dom = dict(in_domain=False, reason=f"T outside the training range {pr.ML_VISC_T_RANGE_K[0]}-{pr.ML_VISC_T_RANGE_K[1]} K")
        if not dom["in_domain"]:
            out["reason"] = dom["reason"]
            return out
        out.update(in_domain=True, g=float(vp.ln_correction(smi, x, T_K)), margin=float(pr._ml_margin_for(T_K)))
    except Exception as e:                                      # model files missing or feature failure: Refutas stays designated
        out["reason"] = f"ml model unavailable ({type(e).__name__})"
    return out


_VSTAT_CAP = 4096


def _visc_status_cached(L: RouteLibrary, pool, w, T_C: float) -> dict:
    """``visc_ml_status`` for a pool recipe, memoised per (members, weights, T) on the library (same values; the value function
    and the recipe margin function of a design both ask for it at every optimiser evaluation)."""
    cache = L.__dict__.setdefault("_vstat", {})
    key = (tuple(pool), np.asarray(w, float).tobytes(), float(T_C))
    st = cache.get(key)
    if st is None:
        if len(cache) >= _VSTAT_CAP:
            cache.clear()
        st = cache[key] = visc_ml_status([L.smiles[n] for n in pool], w, T_C)
    return st


def _visc_designated_fn(L: RouteLibrary, refutas_fn):
    """Design viscosity = the API designated value: Refutas x exp(ML ln correction) where the ML model is in domain, else Refutas."""
    def fn(pool, w, T_C):
        nu = refutas_fn(pool, w, T_C)
        st = _visc_status_cached(L, pool, w, T_C)
        return float(nu * math.exp(st["g"])) if st["in_domain"] else nu
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
    # the yaml key is attached unconditionally: a margins file without it gives the engine's "no such entry" reason (uncertified)
    props["freezing"] = dataclasses.replace(mole, fn=_freezing_api_fn(L), certified_kinds={"jet"},
                                            margin_keys=(FREEZING_MARGIN_KEY,), margin_missing_reason="")
    # viscosity: same value as predict_fuel (ML pair residual designated in its domain, cc5f5ef); the optimiser runs with the
    # Refutas yaml margin of the T bin, the final recipe is re-judged with the ML margin in ``annotate`` (``_rejudge_viscosity``)
    props["viscosity"] = dataclasses.replace(props["viscosity"], fn=_visc_designated_fn(L, props["viscosity"].fn))
    models["props"] = {k: _flag(L, p) for k, p in props.items()}
    models["alternates"] = {k: [_flag(L, p) for p in v] for k, v in alts.items()}
    models["pool"] = list(L.pool)
    models["library"] = L
    return models


LAB_MARGIN_NOTE = "lab_margin_set_by_one_system"
SMOKE_MARGIN_KEY = "smoke_point_mm.lab_mixtures.bulk_regression"
SMOKE_V5_KEY = "smoke_point_mm.ajftd_v5_in_domain.bulk_regression"        # real-fuel margin (Task 1.4); absent until the user decides
SMOKE_OOD_REASON = "outside the smoke regression domain (aromatics 8-23 vol %, result 19-25 mm)"


def smoke_recipe_margin(pool, w, classes: dict, margins: dict, in_domain: bool) -> tuple:
    """Smoke margin of one recipe -> (margin | None, yaml key | None, reason | None).

    Outside the regression domain: None.  If ``SMOKE_V5_KEY`` is PRESENT in ``margins`` it REPLACES the lab key (not max over both)
    and applies only inside the real-fuel class-mix range (``smoke_margin_v5`` gate; outside: None, reason
    ``class_mix_outside_real_fuel_range``); a present but unusable entry gives None, never a fallback to the lab key.  Key absent:
    the lab-mixture margin (7.8 mm) exactly as before, with no class-mix gate."""
    if not in_domain:
        return None, None, SMOKE_OOD_REASON
    if SMOKE_V5_KEY in margins:
        m = E.read_applied_margin(margins[SMOKE_V5_KEY])[0]
        if m is None:
            return None, SMOKE_V5_KEY, "real-fuel smoke margin entry has no valid applied_margin"
        ok, why = SM.class_mix_ok_for_recipe(pool, w, classes)
        return (m, SMOKE_V5_KEY, None) if ok else (None, SMOKE_V5_KEY, why)
    m = E.read_applied_margin(margins.get(SMOKE_MARGIN_KEY))[0]
    return (m, SMOKE_MARGIN_KEY, None) if m is not None else (None, SMOKE_MARGIN_KEY, "smoke point has no valid yaml margin")


class _RecipeMargin:
    """Callable ``margin_fn(pool, w, T_C) -> margin | None`` with a memo (one entry per weight vector, FIFO-bounded) and a
    ``why`` giving the reason of a None.  ``compute(pool, w, T_C) -> (margin | None, reason | None)``."""

    def __init__(self, compute, cap: int = 4096):
        self._compute, self._cap, self._memo = compute, cap, {}

    def _get(self, pool, w, T_C):
        key = (tuple(pool), np.asarray(w, float).tobytes(), None if T_C is None else float(T_C))
        hit = self._memo.get(key)
        if hit is None:
            if len(self._memo) >= self._cap:
                self._memo.clear()
            hit = self._memo[key] = self._compute(list(pool), np.asarray(w, float), T_C)
        return hit

    def __call__(self, pool, w, T_C):
        return self._get(pool, w, T_C)[0]

    def why(self, pool, w, T_C):
        return self._get(pool, w, T_C)[1] or ""


def _smoke_margin_fn(smoke_fn, classes: dict, margins: dict) -> _RecipeMargin:
    def compute(pool, w, T_C):
        with warnings.catch_warnings(record=True) as rec:
            warnings.simplefilter("always")
            smoke_fn(pool, w, T_C)
        in_dom = not any(r.category.__name__.endswith("OutOfDomainWarning") for r in rec)
        m, _key, why = smoke_recipe_margin(pool, w, classes, margins, in_dom)
        return m, why
    return _RecipeMargin(compute)


def _visc_margin_fn(L: RouteLibrary, margins: dict) -> _RecipeMargin:
    def compute(pool, w, T_C):
        if T_C is None:
            return None, "viscosity margin needs a temperature"
        st = _visc_status_cached(L, pool, w, T_C)
        if st["in_domain"]:
            m = st["margin"]
            return (float(m), None) if m is not None and math.isfinite(float(m)) else (None, "ML viscosity margin unavailable")
        key = f"viscosity_kin_cSt.lab_mixtures.raw_refutas.{E.viscosity_t_bin(T_C)}"
        m, _u, _d, why = E.read_applied_margin(margins.get(key))
        return m, (f"{key}: {why}" if m is None else None)
    return _RecipeMargin(compute)


def with_recipe_margins(models: dict, margins: dict) -> dict:
    """A copy of ``models`` whose viscosity and smoke models carry a recipe-dependent ``margin_fn`` (Task 1.5):
    viscosity = ML R-MLV-6 ln margin where the ML model is in domain, else the Refutas lab margin of the T bin; smoke = the
    in-domain smoke margin (``smoke_recipe_margin``).  A model that the route flagged (non-verified constants, ``margin_missing_reason``)
    keeps its uncertified status.  ``models`` itself is not modified."""
    props = dict(models["props"])
    L, classes = models["library"], models.get("classes", {})
    pv = props.get("viscosity")
    if pv is not None and not pv.margin_missing_reason:
        props["viscosity"] = dataclasses.replace(
            pv, margin_fn=_visc_margin_fn(L, margins), margin_ln=True, margin_direction="two_sided",
            margin_source="recipe_margin_fn:viscosity (ML R-MLV-6 ln margin in domain, else Refutas lab margin of the T bin)")
    ps = props.get("smoke_point")
    if ps is not None and not ps.margin_missing_reason:
        key = SMOKE_V5_KEY if SMOKE_V5_KEY in margins else SMOKE_MARGIN_KEY
        props["smoke_point"] = dataclasses.replace(
            ps, margin_fn=_smoke_margin_fn(ps.fn, classes, margins), margin_ln=False, margin_direction="over",
            margin_source=f"recipe_margin_fn:smoke_point ({key}, in domain only)")
    out = dict(models)
    out["props"] = props
    return out


def default_n_starts(pool) -> int:
    return 8 if len(pool) > 20 else 30


SMOKE_MARGIN_REASON = "smoke margin 7.8 mm from binaries and surrogates (user decision), applied inside the regression domain only"


def _rejudge_smoke(result: E.DesignResult, margins: dict, models: dict | None = None, alt_limits: bool = False) -> None:
    """Smoke margin, design side, DOMAIN-GATED (user decision 2026-11-13 plus controller ruling).

    The engine without ``margin_fn`` judges smoke without a margin and tags a recipe outside the regression domain (aromatics 8-23
    vol %, result 19-25 mm) with the OutOfDomainWarning reason, which keeps it ``uncertified`` with margin None.  Here, after the engine,
    a smoke constraint that is IN domain is re-judged with ``smoke_recipe_margin`` (lab key, or the real-fuel key when present, with its
    class-mix gate; the gate needs ``models``): slack = (pred - margin - limit) / scale; ``fail`` if negative.  With ``alt_limits`` the
    slack is max(primary, min(alternative, naphthalenes)) like the optimiser; the naphthalene value comes from ``models`` (unknown ->
    primary only).  Out of domain it is left exactly as the engine produced it.  With ``margin_fn`` models this function is only the
    consistency check (``annotate`` runs it on a copy)."""
    spec_c = {}
    if alt_limits:
        try:
            spec_c = E.load_spec(result.spec_name, alt_limits=True)["constraints"]
        except Exception:
            spec_c = {}
    pool = list(result.pool) or list(result.weights)
    w = np.array([result.weights.get(n, 0.0) for n in pool], float)
    classes = (models or {}).get("classes", {})
    for c in result.constraints:
        if c.prop != "smoke_point" or not math.isfinite(c.prediction) or "OutOfDomainWarning" in c.reason:
            continue
        if c.lo is None or c.hi is not None:                       # smoke is a minimum limit: over-prediction is the unsafe side
            continue
        if SMOKE_V5_KEY in margins and models is None:             # the class-mix gate needs the recipe classes: leave the engine's judgement
            continue
        m, key, _why = smoke_recipe_margin(pool, w, classes, margins, True)
        if m is None:
            continue
        sc = E._scale(c.lo, c.hi)
        c.applied_margin, c.margin_keys = m, (key,)
        slack = (c.prediction - m - c.lo) / sc
        alt = spec_c.get(c.key) or {}
        primary_fail = slack < 0 or c.prediction < c.lo
        branch = None                                              # None = no alternative branch; "fail" / "unknown" / "met"
        if alt_limits and alt.get("alt_min") is not None and alt.get("alt_requires") and models is not None:
            a_slack = (c.prediction - m - float(alt["alt_min"])) / E._scale(alt["alt_min"], None)
            vals, unknown = [a_slack], False
            for rk, rv in alt["alt_requires"].items():
                mo = re.match(r"^(.*)_(max|min)$", rk)
                pm2 = models["props"].get(E.canonical_property(mo.group(1))[0]) if mo else None
                # a requirement is judged without margin only when its model is by_construction (margin 0 legitimately, as in the
                # engine); any other requirement has no margin here, so it is "not established" (never silently margin 0)
                if pm2 is None or not pm2.by_construction:
                    unknown = True
                    continue
                try:
                    q = float(pm2.fn(pool, w, None))
                except Exception:
                    q = float("nan")
                if not math.isfinite(q):
                    unknown = True
                    continue
                vals.append((float(rv) - q) / E._scale(None, rv) if mo.group(2) == "max" else (q - float(rv)) / E._scale(rv, None))
            branch = "fail" if min(vals) < 0 else "unknown" if unknown else "met"
            if not unknown:                                        # same number the optimiser / engine report: max(primary, alt combo)
                slack = max(slack, min(vals))
        c.slack = slack
        if primary_fail:
            if branch in (None, "fail"):
                c.status = E.STATUS_FAIL
            else:                                                  # the alternative branch is met (or not established): at most uncertified
                verified = str(alt.get("verification", "")).strip().lower() == "verified"
                c.status = E.STATUS_PASS if (branch == "met" and verified) else E.STATUS_UNCERTIFIED
        elif slack < 0:
            c.status = E.STATUS_FAIL
        c.reason = (c.reason + "; " if c.reason else "") + SMOKE_MARGIN_REASON + ("; limit violated after applying the margin" if c.status == E.STATUS_FAIL else "")
        result.margins[c.key] = m
    result.passed = result.passed and all(c.status in E.OK_STATUSES for c in result.constraints)


def _rejudge_viscosity(result: E.DesignResult) -> None:
    """Viscosity certification on the API designated value (user decision "Make it main", cc5f5ef).

    In the ML domain the prediction already is the ML value (``_visc_designated_fn``); here its margin becomes the R-MLV-6 ln margin of
    the T bin (two-sided ln-ratio): slack recomputed, ``fail`` if the limit is violated after the margin, else ``pass`` unless an engine
    reason other than the Refutas margin blocks it.  Out of the ML domain the engine's Refutas judgement is left as it is (a note says why)."""
    names = _name_smiles()
    pool = list(result.pool) or list(result.weights)
    smi = [names.get(n, n) for n in pool]
    w = np.array([result.weights.get(n, 0.0) for n in pool], float)
    for c in result.constraints:
        if c.prop != "viscosity" or not math.isfinite(c.prediction):
            continue
        _, T_C = E.canonical_property(c.key)
        st = visc_ml_status(smi, w, T_C)
        if not st["in_domain"]:
            note = f"ml_not_designated: Refutas kept ({st['reason']})"
            c.reason = (c.reason + "; " if c.reason else "") + note
            continue
        m = st["margin"]
        sc = E._scale(c.lo, c.hi)
        slacks = []
        if c.lo is not None:
            slacks.append((c.prediction * math.exp(-m) - c.lo) / sc)
        if c.hi is not None:
            slacks.append((c.hi - c.prediction * math.exp(m)) / sc)
        nominal = not ((c.lo is not None and c.prediction < c.lo) or (c.hi is not None and c.prediction > c.hi))
        keep = [r for r in (c.reason.split("; ") if c.reason else [])
                if r and "raw_refutas" not in r and not r.startswith("limit violated")]
        c.model, c.applied_margin, c.margin_keys = VISC_ML_NAME, m, ("results/v5/BASE/delta_ml_viscosity/readings.json:R-MLV-6",)
        c.slack = float(min(slacks)) if slacks else float("nan")
        if not nominal or c.slack < 0:
            c.status = E.STATUS_FAIL
            keep.append("limit violated after applying the ML margin")
        elif any(not r.startswith("margin validated on") for r in keep):
            c.status = E.STATUS_UNCERTIFIED
        elif c.status not in (E.STATUS_SCOPE,):
            c.status = E.STATUS_PASS
        c.reason = "; ".join([VISC_ML_REASON] + keep)
        result.margins[c.key] = m
        result.predictions[c.key] = c.prediction
    result.passed = bool(all(c.status in E.OK_STATUSES for c in result.constraints) and not result.extrapolating)


def _consistency(result: E.DesignResult, chk: E.DesignResult) -> list:
    """Differences between the engine's judgement of viscosity / smoke (``margin_fn`` path) and the independent post-hoc rejudge ``chk``:
    slack, margin and status (both ways) for viscosity and smoke.  NB ``chk`` starts as a copy of the engine result, so this is only a
    partial check; the independent check is the tests' comparison of the alt-on status against the flags-off legacy path."""
    def same(a, b):
        a = float("nan") if a is None else float(a)
        b = float("nan") if b is None else float(b)
        return (math.isnan(a) and math.isnan(b)) or abs(a - b) <= 1e-9
    out = []
    mine = {c.key: c for c in result.constraints}
    for c2 in chk.constraints:
        c1 = mine[c2.key]
        if c1.prop not in ("viscosity", "smoke_point"):
            continue
        if not same(c1.slack, c2.slack):
            out.append(f"{c1.key}: slack {c1.slack} vs rejudge {c2.slack}")
        if not same(c1.applied_margin, c2.applied_margin):
            out.append(f"{c1.key}: margin {c1.applied_margin} vs rejudge {c2.applied_margin}")
        if c1.status != c2.status:
            out.append(f"{c1.key}: status {c1.status} vs rejudge {c2.status}")
    return out


def annotate(result: E.DesignResult, margins_path=None, models: dict | None = None, alt_limits: bool = False,
             recipe_margins: bool = False) -> E.DesignResult:
    """Every flash constraint carrying a finite margin gets ``lab_margin_set_by_one_system`` in its reason (BASE Task 10:
    the lab margin is the group score of the single mesitylene system).  In-domain smoke constraints get the v5 smoke margin.

    ``recipe_margins`` (the engine already judged viscosity and smoke with ``margin_fn``, see ``with_recipe_margins``): the rejudges run
    on a copy as a consistency check and any difference is listed in ``result.optimiser['margin_fn_consistency']`` (empty = consistent);
    the engine's judgement is kept.  ``models`` / ``alt_limits`` make the smoke rejudge class-mix- and alternative-limit-aware."""
    for c in result.constraints:
        has_margin = c.applied_margin is not None and math.isfinite(c.applied_margin)
        if c.prop == "flash" and has_margin and LAB_MARGIN_NOTE not in c.reason:
            c.reason = (c.reason + "; " if c.reason else "") + LAB_MARGIN_NOTE
    margins = E.load_margins(margins_path)
    if recipe_margins:
        chk = copy.deepcopy(result)
        _rejudge_viscosity(chk)
        _rejudge_smoke(chk, margins, models, alt_limits)
        result.optimiser = dict(result.optimiser, margin_fn_consistency=_consistency(result, chk))
        result.passed = bool(all(c.status in E.OK_STATUSES for c in result.constraints) and not result.extrapolating)
        return result
    _rejudge_viscosity(result)
    _rejudge_smoke(result, margins, models, alt_limits)
    return result


BLOCKED_AROMATICS, BLOCKED_AROMATICS_REASON, blocked_result = E.BLOCKED_AROMATICS, E.BLOCKED_AROMATICS_REASON, E.blocked_result


def design_routed(spec_name, models, n_starts=None, seed=42, margins_path=None, synthetic_blend: bool = False,
                  alt_limits: bool = False, recipe_margins: bool = False) -> E.DesignResult:
    """``synthetic_blend``: add the D7566 aromatics minimum (8 vol %); a pool with no aromatic member then cannot meet it
    and the result is category ``blocked_aromatics_min_infeasible`` (no optimisation run).  ``alt_limits``: honour the
    D1655 alternative smoke limit (spec flag ``use_alt_limits``; the caller's spec is not mutated, it is re-read).
    ``recipe_margins`` (default off = results identical to before): viscosity and smoke carry recipe-dependent ``margin_fn`` margins
    (``with_recipe_margins``), so the optimiser sees them; ``annotate`` then only checks consistency."""
    pool = models["pool"]
    blocked = E.blocked_result(spec_name, pool, models, seed, synthetic_blend)
    if blocked is not None:
        return blocked
    use = with_recipe_margins(models, E.load_margins(margins_path)) if recipe_margins else models
    return annotate(E.design(spec_name, pool, use, n_starts=n_starts or default_n_starts(pool), seed=seed,
                             margins_path=margins_path, synthetic_blend=synthetic_blend, alt_limits=alt_limits),
                    margins_path, models=models, alt_limits=alt_limits, recipe_margins=recipe_margins)


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


def evaluate_recipe(spec_name: str, models: dict, weights: dict, margins_path=None, synthetic_blend: bool = False,
                    alt_limits: bool = False, recipe_margins: bool = False) -> E.DesignResult:
    """A ``DesignResult`` for a given recipe (no optimisation), scored by the engine with ``models``."""
    pool = list(models["pool"]) if "pool" in models else list(weights)
    blocked = E.blocked_result(spec_name, pool, models, synthetic_blend=synthetic_blend)
    if blocked is not None:
        return blocked
    w = np.array([weights.get(n, 0.0) for n in pool], float)
    w = w / w.sum()
    spec = E.load_spec(spec_name, synthetic_blend=synthetic_blend, alt_limits=alt_limits)
    margins = E.load_margins(margins_path)
    use = with_recipe_margins(models, margins) if recipe_margins else models
    cons = E.evaluate_blend(spec, pool, w, use, margins)
    return annotate(E.DesignResult(weights=dict(zip(pool, map(float, w))), predictions={c.key: c.prediction for c in cons},
                          margins={c.key: c.applied_margin for c in cons}, passed=False, extrapolating=False,
                          binding_constraints=[], constraints=cons, spec_name=spec_name, pool=pool), margins_path,
                    models=models, alt_limits=alt_limits, recipe_margins=recipe_margins)


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
        "freezing_point_C": "design and API both use ideal Schroeder-van Laar in mole fraction on the designated freezing constants (verified csv, else Coutinho)",
        "nhoc_MJ_kg": "design primary model uses constant H = 13.8 wt%; the API uses composition H",
        "smoke_point_mm": "design and API both use the AJFTD bulk regression on aromatics vol % (API aromatics from the grid rho at 288.15 K)",
        "aromatics_wt_pct": "class mass fraction, by construction",
        "aromatics_vol_pct": "design: table rho15 volume fractions; API: grid rho at 288.15 K",
    }
    if kv:
        notes[kv] = ("kinematic viscosity at the spec temperature, API designated value on both sides: Refutas x ML pair residual in the "
                     "ML domain, else Refutas on the pure-property grid")
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
        var = f_vis["viscosity_kin_cSt"].get("variants") or {}
        rr = var.get("refutas_raw")
        if rr is not None:                                  # ML designated: show the Refutas value it replaced
            rows.append(_Row(kv + "[refutas_raw]", dp(kv), _fnum(rr.get("value")), float("nan"),
                             "informational: API Refutas variant (replaced by the designated ML value in its domain)", True))
        else:
            mlv = var.get("ml_pair_residual") or {}
            rows.append(_Row(kv + "[ml_pair_residual]", dp(kv), _fnum(mlv.get("value")), float("nan"),
                             "informational: ML not designated here; "
                             + ("; ".join(mlv.get("warnings") or []) or str(mlv.get("reason") or "no ML value")), True))
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
