"""v5 Task 6: inverse design on a SMILES / InChIKey pool, with an honest status for every constant.

``RouteLibrary`` is a ``models_v4.Library`` whose pool may be given as table names, SMILES or InChIKeys.  A SMILES or
InChIKey that matches a compound of the verified constants table uses the table name and the verified constants (identical
design to the name-based pool).  A compound that is not in the table gets route constants, each with a status:

    MW                  RDKit                                 verified
    rho15               pure-property grid at 288.15 K        route_grid
    Tm, dHfus           ``build_tmap`` (Coutinho table)       coutinho_table
    flash               anchored Antoine / LFL registry       route_anchored
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
    props["freezing"] = dataclasses.replace(mole, certified_kinds={"jet"})
    models["props"] = {k: _flag(L, p) for k, p in props.items()}
    models["alternates"] = {k: [_flag(L, p) for p in v] for k, v in alts.items()}
    models["pool"] = list(L.pool)
    models["library"] = L
    return models


LAB_MARGIN_NOTE = "lab_margin_set_by_one_system"


def default_n_starts(pool) -> int:
    return 8 if len(pool) > 20 else 30


def annotate(result: E.DesignResult) -> E.DesignResult:
    """Every flash constraint carrying a finite margin gets ``lab_margin_set_by_one_system`` in its reason (BASE Task 10:
    the lab margin is the group score of the single mesitylene system)."""
    for c in result.constraints:
        has_margin = c.applied_margin is not None and math.isfinite(c.applied_margin)
        if c.prop == "flash" and has_margin and LAB_MARGIN_NOTE not in c.reason:
            c.reason = (c.reason + "; " if c.reason else "") + LAB_MARGIN_NOTE
    return result


def design_routed(spec_name, models, n_starts=None, seed=42, margins_path=None) -> E.DesignResult:
    pool = models["pool"]
    return annotate(E.design(spec_name, pool, models, n_starts=n_starts or default_n_starts(pool), seed=seed,
                             margins_path=margins_path))


# ------------------------------------------------------------------------------------------------ forward check
@dataclasses.dataclass
class _Row:
    constraint: str
    design_pred: float
    forward_pred: float
    tol: float
    cause: str = ""
    informational: bool = False


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
                          binding_constraints=[], constraints=cons, spec_name=spec_name, pool=pool))


def design_forward_check(result: E.DesignResult, T_C_visc: float | None = None) -> pd.DataFrame:
    """Re-predict ``result.weights`` through ``predict_fuel`` and compare with the design-side predictions.

    Tolerances (R-API-4): flash 0.01 C against the API designated (overlay) flash, viscosity 1e-6
    relative, net heat exact (1e-9, like-for-like composition-H model).  Density and freezing differences are reported
    with the cause, not judged."""
    from fuelblend.api import predict_route as pr
    from fuelblend.baselines.heat_of_combustion_v2 import mixture_heat_of_combustion_v2
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
        return float(c.prediction) if c is not None else float("nan")

    if "density_15C_kg_m3" in cons:
        rows.append(_Row("density_15C_kg_m3", dp("density_15C_kg_m3"), f_15["density_kg_m3"]["value"], float("nan"),
                         "design: ideal volume mixing of the table rho15; forward: pure-property grid rho at 288.15 K "
                         "(table rho15 is a listed value, the grid is a fitted curve)", True))
    if "flash_point_C" in cons:
        fl = f_15["flash_point_C"]
        rows.append(_Row("flash_point_C", dp("flash_point_C"), fl["value"], 0.01,
                         "design and API both use the LFL-overlay registry (U1)"))
        rows.append(_Row("flash_point_C[v4_designated]", dp("flash_point_C"), fl["variants"]["v4_designated"], float("nan"),
                         "informational: plain registry (no LFL overlay); the diff is the overlay effect", True))
    if "freezing_point_C" in cons:
        rows.append(_Row("freezing_point_C", dp("freezing_point_C"), f_15["freezing_point_C"]["value"], float("nan"),
                         "design: Tm/dHfus from pure_props_v5 (listed values) in mole-fraction ideal SvL; forward: Tm/dHfus "
                         "from build_tmap (Coutinho table). Same equation, different constants", True))
    if kv:
        rows.append(_Row(kv, dp(kv), f_vis["viscosity_kin_cSt"]["value"], 1e-6))
    if "nhoc_MJ_kg" in cons:
        # The API uses the composition-H group-additivity model; the design primary uses constant H.  Two like-for-like rows
        # use the design's composition-H model: as designed (table MW, 2 decimals) and with RDKit MW, the API's MW.
        w = np.array([result.weights[n] for n in pool], float)
        fn = M.make_fns(L)["nhoc_comp_h"]
        rows.append(_Row("nhoc_MJ_kg[compH_table_MW]", fn(pool, w, None), f_15["heat_net_MJ_kg"]["value"], float("nan"),
                         "same composition-H model; the design reads listed MW (2 decimals) from pure_props_v5, the API reads RDKit MW", True))
        Lr = RouteLibrary(pool)
        for n in pool:
            Lr._tab[(n, "MW")] = (float(Descriptors.MolWt(Chem.MolFromSmiles(smi[n]))), "verified")
        rows.append(_Row("nhoc_MJ_kg", M.make_fns(Lr)["nhoc_comp_h"](pool, w, None), f_15["heat_net_MJ_kg"]["value"], 1e-9,
                         "design-side composition-H group-additivity model with RDKit MW (alternate model; primary is constant H)"))
        rows.append(_Row("nhoc_MJ_kg[primary_constH]", dp("nhoc_MJ_kg"), f_15["heat_net_MJ_kg"]["value"], float("nan"),
                         "design primary model uses constant H = 13.8 wt%; the API uses composition H", True))
    for key, fk in (("aromatics_wt_pct", "aromatics_wt_pct"), ("aromatics_vol_pct", "aromatics_vol_pct")):
        if key in cons:
            rows.append(_Row(key, dp(key), f_15[fk]["value"], float("nan"),
                             "design: table rho15 volume fractions; forward: grid rho at 288.15 K" if "vol" in key else
                             "class mass fraction, by construction", True))
    out = []
    for r in rows:
        d = abs(r.design_pred - r.forward_pred) if (r.forward_pred is not None and np.isfinite(r.design_pred)) else float("nan")
        if r.informational:
            ok = "reported"
        elif r.constraint.startswith("kin_visc_"):
            ok = bool(d <= r.tol * abs(r.forward_pred))
        else:
            ok = bool(d <= r.tol)
        out.append(dict(constraint=r.constraint, design_pred=r.design_pred, forward_pred=r.forward_pred, abs_diff=d,
                        tol=r.tol, **{"pass": ok}, cause=r.cause))
    return pd.DataFrame(out, columns=["constraint", "design_pred", "forward_pred", "abs_diff", "tol", "pass", "cause"])
