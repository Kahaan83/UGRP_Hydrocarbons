"""v5 API route: route context (registry + frozen provider grid + SubReps) and per-property component mapping.

Mirrors the BASE scoring route (``fuelblend.eval.composition_route``): a component missing own data for a property is
replaced, for that property only, by ``SubReps.get``; non-hydrocarbons and components with no representative are unmapped.
The provider is read-only: ``ensure()`` never runs and the npz is never written.
"""
from __future__ import annotations

import importlib.util
import math
import re
import warnings
from dataclasses import asdict, dataclass
from functools import lru_cache
from pathlib import Path
from typing import Callable, Optional

import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors, inchi

from fuelblend.baselines.flash_unifac import _lookup, flash_point_unifac_detail
from fuelblend.data.reconstruction_r3 import SubReps, subkind
from fuelblend.eval.composition_route import NEAREST_ANY, RouteCtx, build_tmap, refutas_nu_safe
from fuelblend.mixture.common import CompoundRegistry, PurePropertyProvider

RDLogger.DisableLog("rdApp.*")

ROOT = Path(__file__).resolve().parents[3]
GRID_NPZ = ROOT / "data" / "cache" / "pure_grid.npz"
COMPOUNDS_CSV = ROOT / "data" / "unified" / "compounds.csv"
_IK_RE = re.compile(r"^[A-Z]{14}-[A-Z]{10}-[A-Z]$")


class ReadOnlyProvider(PurePropertyProvider):
    """Frozen pure-property grid: loads the npz, never computes or writes."""

    def __init__(self, reg, cache_path=GRID_NPZ):
        if not Path(cache_path).exists():
            raise FileNotFoundError(f"pure-property grid not found: {cache_path}")
        super().__init__(reg, cache_path, verbose=False)

    def ensure(self, comp_idx):                       # noqa: D401 - deliberate no-op
        return None


@lru_cache(maxsize=1)
def get_ctx() -> RouteCtx:
    reg = CompoundRegistry(COMPOUNDS_CSV)
    prov = ReadOnlyProvider(reg, GRID_NPZ)
    reps = SubReps(reg, prov, vp_lookup=lambda s: _lookup(s, None))
    ctx = RouteCtx(reg, prov)
    ctx.reps = reps
    return ctx


@dataclass
class Component:
    ik: str
    smiles: str
    w: float
    reg_idx: Optional[int]
    subkind: Optional[str]
    n_carbon: int
    is_hc: bool


def _mol_info(smi: str):
    m = Chem.MolFromSmiles(smi)
    can = Chem.MolToSmiles(m)
    nc = sum(a.GetSymbol() == "C" for a in m.GetAtoms())
    hc = nc > 0 and all(a.GetSymbol() in ("C", "H") for a in m.GetAtoms())
    return m, can, nc, hc


def resolve(composition, mole_fractions: Optional[bool] = None) -> list[Component]:
    """composition: {id: weight} or [(id, weight), ...]; id is an InChIKey or SMILES.  Weights are mass fractions unless
    ``mole_fractions`` is truthy.  Duplicates are merged by InChIKey, zero weights dropped, weights normalised."""
    ctx = get_ctx()
    reg = ctx.reg
    items = list(composition.items()) if hasattr(composition, "items") else list(composition)
    acc: dict[str, list] = {}
    for key, wt in items:
        wt = float(wt)
        if not np.isfinite(wt) or wt < 0:
            raise ValueError(f"invalid weight {wt!r} for {key!r}")
        key = str(key).strip()
        if _IK_RE.match(key):
            i = reg.ik2i.get(key)
            if i is None:
                raise ValueError(f"InChIKey not in registry: {key}")
            smi = reg.smiles[i]
        else:
            m = Chem.MolFromSmiles(key)
            if m is None:
                raise ValueError(f"unparseable SMILES: {key!r}")
            smi = Chem.MolToSmiles(m)
            key = inchi.MolToInchiKey(m)
            i = reg.ik2i.get(key)
            if i is not None:
                smi = reg.smiles[i]
        if key in acc:
            acc[key][1] += wt
        else:
            acc[key] = [smi, wt, i]
    out = []
    for ik, (smi, wt, i) in acc.items():
        if mole_fractions:
            wt = wt * ctx.mw(smi)
        if wt > 0:
            out.append([ik, smi, wt, i])
    tot = sum(o[2] for o in out)
    if not out or tot <= 0:
        raise ValueError("composition has no positive weight")
    comps = []
    for ik, smi, wt, i in out:
        _, can, nc, hc = _mol_info(smi)
        comps.append(Component(ik=ik, smiles=smi, w=wt / tot, reg_idx=i, subkind=subkind(smi) if hc else None,
                               n_carbon=nc, is_hc=hc))
    return comps


@dataclass
class PropMap:
    idx: np.ndarray
    w: np.ndarray
    own_mass: float
    substituted_mass: float
    unmapped_mass: float
    substitutions: list


def map_for(comps: list, has_own: Callable[[int], bool], reps) -> PropMap:
    """Per-property mapping.  Own data wins; else a SubReps representative (gap <= 3, 99 for NEAREST_ANY kinds); else unmapped.
    Weights are renormalised over the mapped mass; representatives are aggregated as in ``map_bins``."""
    agg: dict[int, float] = {}
    own = sub = lost = 0.0
    subs = []
    for c in comps:
        if not c.is_hc:
            lost += c.w
            continue
        if c.reg_idx is not None and has_own(c.reg_idx):
            agg[c.reg_idx] = agg.get(c.reg_idx, 0.0) + c.w
            own += c.w
            continue
        rep = None
        if c.subkind is not None:
            rep = reps.get(c.subkind, c.n_carbon, max_gap=99 if c.subkind in NEAREST_ANY else 3)
            if rep is not None and not has_own(rep):      # a representative must itself carry own data for this property
                rep = None
        if rep is None:
            lost += c.w
            continue
        agg[rep] = agg.get(rep, 0.0) + c.w
        sub += c.w
        rep_nc = _mol_info(reps.smiles[rep])[2]
        subs.append((c.ik, reps.ik[rep], abs(rep_nc - c.n_carbon)))
    idx = np.array(sorted(agg), dtype=np.int64)
    w = np.array([agg[i] for i in idx], float)
    w = w / w.sum() if w.size and w.sum() > 0 else w
    return PropMap(idx=idx, w=w, own_mass=own, substituted_mass=sub, unmapped_mass=lost, substitutions=subs)


@lru_cache(maxsize=1)
def _tmap() -> dict:
    ctx = get_ctx()
    ok = [i for i in range(len(ctx.reg.ik)) if ctx.reg.usable[i]]
    return build_tmap(ctx.reg, ok)


@lru_cache(maxsize=1)
def _ysi_smiles() -> frozenset:
    from fuelblend.baselines.sooting import ysi_table      # data/processed/v4/sooting/ysi_hc.csv
    return frozenset(ysi_table())


@lru_cache(maxsize=1)
def _overlay_ctx() -> RouteCtx:
    """Designated-flash context (BASE Task 10, U1): tabulated LFL values replace still-estimated ones.  Same call that produced
    ``pred_lfl_overlay``; the script is loaded by path and not edited."""
    base = get_ctx()
    spec = importlib.util.spec_from_file_location("cells_t2_flash_api", ROOT / "scripts" / "v5" / "base" / "cells_t2_flash.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    ctx = RouteCtx(base.reg, base.prov, lfl_override=mod.lfl_overlay_dict())
    ctx.reps = base.reps
    return ctx


def predicates(ctx: Optional[RouteCtx] = None) -> dict:
    """Own-data predicates over registry indices: HAS_GRID, HAS_FLASH, HAS_SLE, HAS_YSI."""
    ctx = ctx or get_ctx()
    reg, prov = ctx.reg, ctx.prov

    def has_sle(i):
        t = _tmap().get(reg.ik[i])
        return bool(t) and bool(np.isfinite(t.get("tm_K", np.nan)) and np.isfinite(t.get("dh_fus_J_mol", np.nan)))

    def has_ysi(i):
        m = Chem.MolFromSmiles(reg.smiles[i])
        return m is not None and Chem.MolToSmiles(m) in _ysi_smiles()

    return {"HAS_GRID": lambda i: bool(prov.valid[i]),
            "HAS_FLASH": lambda i: _lookup(reg.smiles[i], None) is not None,
            "HAS_SLE": has_sle, "HAS_YSI": has_ysi}


# ====================================================================================================================
# predict_fuel: viscosity, density and flash on the BASE scoring route
# ====================================================================================================================
LOW_COVERAGE_MASS = 0.10                 # unmapped mass above this gives value None / status low_coverage
MESITYLENE_IK = inchi.MolToInchiKey(Chem.MolFromSmiles("Cc1cc(C)cc(C)c1"))
FLASH_METHOD = "ideal_le_chatelier_anchored"
_UNITS = {"viscosity_kin_cSt": "cSt", "viscosity_dyn_mPas": "mPa.s", "density_kg_m3": "kg/m3", "flash_point_C": "degC"}


def _coverage(pm: PropMap) -> dict:
    return {"own_mass": float(pm.own_mass), "substituted_mass": float(pm.substituted_mass),
            "unmapped_mass": float(pm.unmapped_mass),
            "substitutions": [{"input": a, "representative": b, "carbon_gap": int(g)} for a, b, g in pm.substitutions]}


def _new_prop(prop: str, method: str, pm: PropMap) -> dict:
    res = {"value": None, "std": None, "unit": _UNITS[prop], "method": method, "status": "no_evaluation",
           "coverage": _coverage(pm), "lo": math.nan, "hi": math.nan, "margin": None, "in_domain": False,
           "indicative_interval": None, "variants": {}, "warnings": []}
    cov = res["coverage"]
    if cov["substituted_mass"] > 0:
        res["warnings"].append(f"substituted_components: {len(pm.substitutions)} component(s), mass {cov['substituted_mass']:.3f}, "
                               "represented by a nearest same-kind compound for this property")
    if 0 < cov["unmapped_mass"] <= LOW_COVERAGE_MASS:
        res["warnings"].append(f"unmapped_mass_dropped: mass {cov['unmapped_mass']:.3f} has no data for this property; renormalised")
    elif cov["unmapped_mass"] > LOW_COVERAGE_MASS:
        res["warnings"].append(f"low_coverage: unmapped mass {cov['unmapped_mass']:.3f} > {LOW_COVERAGE_MASS}")
    return res


def _no_value(res: dict, status: str, reason: str) -> dict:
    res.update(value=None, method="none", status=status, model_name="none", interval=[math.nan, math.nan],
               domain_reasons=[reason], tuple=(None, math.nan, math.nan, "none", False))
    res["warnings"].append(f"{status}: {reason}")
    return res


def _finish(res: dict, yaml_prop: str, iks: list, T_K: float, fuel_type, kind: str) -> dict:
    """Attach lo/hi/in_domain.  The lab/yaml margins apply only when nothing was substituted or unmapped."""
    from fuelblend.api import predict_blend as pb
    if res["value"] is None:
        return res
    cov = res["coverage"]
    if cov["substituted_mass"] > 1e-12 or cov["unmapped_mass"] > 1e-12:
        res.update(model_name=res["method"], interval=[math.nan, math.nan], domain_reasons=["substituted_or_unmapped_mass"],
                   tuple=(res["value"], math.nan, math.nan, res["method"], False))
        return res
    return pb._attach_uncertainty(yaml_prop, res, iks, T_K, fuel_type, kind)


def predict_fuel(composition, T_K: float = 298.15, fuel_type: Optional[str] = None, mole_fractions=None,
                 composition_kind: str = "known") -> dict:
    """Viscosity (kinematic and dynamic), density and flash point of a hydrocarbon blend on the BASE scoring route.

    composition: {id: weight} or [(id, weight), ...]; id = InChIKey or SMILES; weights are mass fractions (``mole_fractions``,
    a list aligned with the composition order, replaces them with mole fractions).  Returns
    ``{T_K, warnings, components, <property>: {value, unit, method, status, coverage, lo, hi, margin, in_domain,
    indicative_interval, variants, warnings, ...}}``.  A property with no yaml margin, or with substituted/unmapped mass,
    has NaN lo/hi and in_domain False.  Unparseable SMILES raise ValueError."""
    from fuelblend.api import predict_blend as pb
    if composition_kind not in pb.COMPOSITION_KINDS:
        raise ValueError(f"composition_kind must be one of {pb.COMPOSITION_KINDS}, got {composition_kind!r}")
    T_K = float(T_K)
    items = list(composition.items()) if hasattr(composition, "items") else list(composition)
    if mole_fractions is not None:
        if len(mole_fractions) != len(items):
            raise ValueError("mole_fractions must be aligned with the composition")
        comps = resolve([(k, float(m)) for (k, _), m in zip(items, mole_fractions)], mole_fractions=True)
    else:
        comps = resolve(items)
    ctx = get_ctx()
    P = predicates(ctx)
    iks = [c.ik for c in comps]
    top_warn = [f"non_hc_input: {c.ik} is not a hydrocarbon; unmapped for physics, margins off" for c in comps if not c.is_hc]
    out = {"T_K": T_K, "warnings": top_warn, "components": [asdict(c) for c in comps]}

    # ---- viscosity and density: frozen provider grid, Refutas / ideal volume
    pm = map_for(comps, P["HAS_GRID"], ctx.reps)
    recal = (composition_kind == "reconstructed" and pb.VISC_RECAL_T_C[0] + 273.15 <= T_K <= pb.VISC_RECAL_T_C[1] + 273.15)
    vmethod = "refutas_raw" if composition_kind == "known" else ("refutas_recalibrated" if recal else "refutas")
    props = {"viscosity_kin_cSt": _new_prop("viscosity_kin_cSt", vmethod, pm),
             "viscosity_dyn_mPas": _new_prop("viscosity_dyn_mPas", vmethod, pm),
             "density_kg_m3": _new_prop("density_kg_m3", "dens_idvol", pm)}
    if not (ctx.prov.T[0] <= T_K <= ctx.prov.T[-1]):
        for r in props.values():
            _no_value(r, "T_out_of_grid", f"T={T_K:.2f} K outside the pure-property grid [{ctx.prov.T[0]:.0f}, {ctx.prov.T[-1]:.0f}] K")
    elif pm.unmapped_mass > LOW_COVERAGE_MASS or pm.idx.size == 0:
        for r in props.values():
            _no_value(r, "low_coverage", f"unmapped mass {pm.unmapped_mass:.3f} > {LOW_COVERAGE_MASS}")
    else:
        nu, rho = refutas_nu_safe(ctx.prov, pm.idx, pm.w, T_K)
        if recal:
            b0, b1 = pb._recal_beta()
            nu = nu * math.exp(b0 + b1 * 1000.0 / T_K)
        elif composition_kind == "reconstructed":
            props["viscosity_kin_cSt"]["warnings"].append(
                f"recalibration_not_applied: T={T_K - 273.15:.1f} C outside {pb.VISC_RECAL_T_C} C; raw Refutas, no certified margin")
        for name, v in (("viscosity_kin_cSt", nu), ("viscosity_dyn_mPas", nu * rho / 1000.0), ("density_kg_m3", rho)):
            if np.isfinite(v):
                props[name].update(value=float(v), status="ok")
            else:
                _no_value(props[name], "non_finite", "Refutas / ideal-volume result is not finite")
    for name, r in props.items():
        yp = "viscosity_dyn_mPas" if name.startswith("viscosity") else name
        out[name] = _finish(r, yp, iks, T_K, fuel_type, composition_kind)

    # ---- flash point: ideal Le Chatelier on the anchored Antoine / LFL registry (does not depend on T)
    pf = map_for(comps, P["HAS_FLASH"], ctx.reps)
    fl = _new_prop("flash_point_C", FLASH_METHOD, pf)
    if pf.unmapped_mass > LOW_COVERAGE_MASS or pf.idx.size == 0:
        _no_value(fl, "low_coverage", f"unmapped mass {pf.unmapped_mass:.3f} > {LOW_COVERAGE_MASS}")
    else:
        smi = [ctx.reg.smiles[i] for i in pf.idx]
        x = pf.w / np.array([ctx.mw(s) for s in smi])
        x = x / x.sum()

        def _flash(c):
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                d = flash_point_unifac_detail(smi, x, registry=c.flash_reg, use_unifac=False)
            return (d["flash_C"] if d["status"] == "ok" and np.isfinite(d["flash_C"]) else None), d["status"]

        val, st = _flash(_overlay_ctx())
        plain, _ = _flash(ctx)
        fl["variants"]["v4_designated"] = plain
        if val is None:
            _no_value(fl, "non_finite", f"flash solver status {st}")
        else:
            fl.update(value=float(val), status="ok")
            if not (0.0 <= val <= 100.0):
                fl["warnings"].append("antoine_nominal_window_exceeded: flash point outside 0-100 C; the anchored Antoine "
                                      "window (273.15-373.15 K) is nominal")
    if MESITYLENE_IK in iks:
        fl["warnings"].append("antoine_suspect:mesitylene")
    out["flash_point_C"] = _finish(fl, "flash_point_C", iks, T_K, fuel_type, composition_kind)
    if np.isfinite(fl["lo"]) and np.isfinite(fl["hi"]):
        fl["warnings"].append("lab_margin_set_by_one_system")
    return out
