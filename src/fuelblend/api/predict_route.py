"""v5 API route: route context (registry + frozen provider grid + SubReps) and per-property component mapping.

Mirrors the BASE scoring route (``fuelblend.eval.composition_route``): a component missing own data for a property is
replaced, for that property only, by ``SubReps.get``; non-hydrocarbons and components with no representative are unmapped.
The provider is read-only: ``ensure()`` never runs and the npz is never written.
"""
from __future__ import annotations

import hashlib
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
from fuelblend.eval.composition_route import (NEAREST_ANY, VERIFIED_SOURCE, RouteCtx, build_tmap, build_tmap_designated,
                                              build_tmap_verified, refutas_nu_safe)
from fuelblend.mixture.common import CompoundRegistry, PurePropertyProvider

RDLogger.DisableLog("rdApp.*")

ROOT = Path(__file__).resolve().parents[3]
GRID_NPZ = ROOT / "data" / "cache" / "pure_grid.npz"
COMPOUNDS_CSV = ROOT / "data" / "unified" / "compounds.csv"
_IK_RE = re.compile(r"^[A-Z]{14}-[A-Z]{10}-[A-Z]$")
GRID_SHA256 = "b931e969926698dff52d69ff51dd99a990c24ff12eca545aa98e970afd3e7daf"      # pinned frozen pure-property grid


@lru_cache(maxsize=1)
def grid_sha_status() -> tuple:
    """(matches the pinned GRID_SHA256, actual sha256 hex).  Hashed once per process; never raises on a mismatch."""
    h = hashlib.sha256()
    with open(GRID_NPZ, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest() == GRID_SHA256, h.hexdigest()


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
def _tmap_freezing() -> dict:
    """Freezing constants: tm_dhfus_verified.csv where present, else the Coutinho entry (tagged ``source``)."""
    return build_tmap_designated(_tmap(), build_tmap_verified())


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
_UNITS = {"viscosity_kin_cSt": "cSt", "viscosity_dyn_mPas": "mPa.s", "density_kg_m3": "kg/m3", "flash_point_C": "degC",
          "freezing_point_C": "degC", "heat_net_MJ_kg": "MJ/kg", "smoke_point_mm": "mm", "aromatics_wt_pct": "wt%",
          "aromatics_vol_pct": "vol%"}


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


FLASH_MARGINS_YAML = ROOT / "configs" / "uncertainty_v5.yaml"     # generated by scripts/v5/base/build_uncertainty_v5_q2.py (user decision Q2)
_flash_margins_cache: Optional[dict] = None


def _v5_margins() -> dict:
    global _flash_margins_cache
    if _flash_margins_cache is None:
        import yaml
        _flash_margins_cache = yaml.safe_load(FLASH_MARGINS_YAML.read_text(encoding="utf-8"))["margins"]
    return _flash_margins_cache


SMOKE_MARGIN_KEY = "smoke_point_mm.lab_mixtures.bulk_regression"
SMOKE_MARGIN_NOTE = "smoke_margin_from_binaries_and_surrogates (user decision)"          # lab key 7.8 mm (used only when the v5 key is absent)
SMOKE_MARGIN_NOTE_V5 = "smoke_margin_from_real_fuel_ajftd_v5_in_domain (user decision R-SC-2b)"
SMOKE_V5_KEY = "smoke_point_mm.ajftd_v5_in_domain.bulk_regression"        # real-fuel margin; replaces the lab key when present in the yaml


@lru_cache(maxsize=1)
def _registry_classes() -> list:
    """Library class (n_paraffin / iso_paraffin / cycloalkane / aromatic / ...) per registry index (same file and order as the design side)."""
    from fuelblend.design import models_v4 as M
    from fuelblend.mixture.cold_flow import compound_meta
    return compound_meta(M.COMPOUNDS)[1]


def _smoke_class_mix_ok(comps) -> tuple:
    """Class-mix gate of the real-fuel smoke margin on an API composition (same basis as the design side)."""
    from fuelblend.design import smoke_margin_v5 as SM
    cls = _registry_classes()
    pool = [c.ik for c in comps]
    classes = {c.ik: (cls[c.reg_idx] if c.reg_idx is not None else "unknown") for c in comps}
    return SM.class_mix_ok_for_recipe(pool, [c.w for c in comps], classes)


def _flash_margin_v5(res: dict) -> None:
    """Flash only: replace the v4 yaml margin by the one in configs/uncertainty_v5.yaml (same keys, max of applied_margin).
    Other properties and the legacy predict_blend paths keep uncertainty_v4.yaml.  Domain reasons from _attach_uncertainty stay."""
    from fuelblend.api import predict_blend as pb
    keys = pb.MARGIN_KEYS.get(("flash_point_C", res["method"]))
    if keys is None or res.get("value") is None:
        return
    margin, unit, direction, reason = pb._applied_margin(keys, _v5_margins())
    if margin is None:
        return
    lo, hi = pb.interval_from_margin(res["value"], margin, unit)
    res.update(lo=lo, hi=hi, interval=[lo, hi], margin=margin, margin_unit=unit, margin_direction=direction)
    if direction == "over":
        res["safe_bound"] = lo
    res["tuple"] = (res["value"], lo, hi, res["method"], res["in_domain"])


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
    res = pb._attach_uncertainty(yaml_prop, res, iks, T_K, fuel_type, kind)
    if yaml_prop == "flash_point_C":
        _flash_margin_v5(res)
    return res


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
    top_warn = []
    sha_ok, sha = grid_sha_status()
    if not sha_ok:
        top_warn.append(f"grid_sha_mismatch: data/cache/pure_grid.npz sha256 {sha[:12]}... differs from the pinned "
                        f"{GRID_SHA256[:12]}...; viscosity, density and design numbers may differ from the BASE evidence")
    top_warn += [f"non_hc_input: {c.ik} is not a hydrocarbon; unmapped for physics, margins off" for c in comps if not c.is_hc]
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
    exceeded = False
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
                exceeded = True
                fl["warnings"].append("antoine_nominal_window_exceeded: flash point outside 0-100 C; the anchored Antoine "
                                      "window (273.15-373.15 K) is nominal")
    if MESITYLENE_IK in iks or any(r == MESITYLENE_IK for _, r, _ in pf.substitutions):   # input, or entering as a representative
        fl["warnings"].append("antoine_suspect:mesitylene")
    out["flash_point_C"] = _finish(fl, "flash_point_C", iks, T_K, fuel_type, composition_kind)
    if exceeded:                                       # outside the model's nominal Antoine window: no interval, not in domain
        nan = math.nan
        fl.update(lo=nan, hi=nan, interval=[nan, nan], margin=None, in_domain=False,
                  domain_reasons=list(fl.get("domain_reasons") or []) + ["antoine_nominal_window_exceeded"],
                  tuple=(fl["value"], nan, nan, fl["method"], False))
    if np.isfinite(fl["lo"]) and np.isfinite(fl["hi"]):
        fl["warnings"].append("lab_margin_set_by_one_system")
        fl["warnings"].append("flash_margin_excludes_mesitylene (user decision Q2)")
    _add_viscosity_k1(out["viscosity_kin_cSt"], T_K)
    _add_viscosity_ml(out, ctx, pm, T_K, composition_kind)
    _promote_ml_viscosity(out)
    _add_density_ml(out, ctx, pm, T_K, composition_kind)
    _add_task3_props(out, comps, ctx, fuel_type, T_K, composition_kind)
    _attach_indicative(out, fuel_type, T_K)
    return out


# ====================================================================================================================
# freezing, net heat, smoke point, aromatics, and the K1 / F1x labelled variants
# ====================================================================================================================
SMOKE_YSI_MIN = 0.90
HEAT_YAML_PROP = "heat_net_MJ_kg"                      # predict_blend.MARGIN_KEYS maps (this, "ga_comp_H") to the crc647 yaml key
K1_TOL_K = 0.01
HAS_ANY_REGISTRY = lambda i: True                       # noqa: E731 - the scoring route's map_bins does not ask for property data


@lru_cache(maxsize=1)
def _pure_data():
    """Group-additivity GCV model and table (copy of the ``PureData`` set-up in scripts/v5/base/cells_t2_bulk.py)."""
    import pandas as pd
    from fuelblend.baselines.heat_of_combustion_v2 import GroupAdditivityGCVModel
    g = pd.read_csv(ROOT / "data" / "processed" / "v2" / "gcv_v2_MJ_kg.csv").set_index("inchikey")
    tier_a = g[g.tier == "A"].reset_index()
    mw = {k: Descriptors.MolWt(Chem.MolFromSmiles(s)) for k, s in zip(tier_a.inchikey, tier_a.smiles)}
    model = GroupAdditivityGCVModel()
    model.fit_from_tier_a(tier_a.assign(mw=tier_a.inchikey.map(mw)))
    return model, g


_COMP_CACHE: dict = {}


def _comp_info(ctx, i: int) -> dict:
    """Per registry compound: smiles, mw, GCV (MJ/kg), carbon count, YSI availability.  Same logic as ``PureData.comp_info``."""
    if i not in _COMP_CACHE:
        model, g = _pure_data()
        ik, smi = ctx.reg.ik[i], ctx.reg.smiles[i]
        m = Chem.MolFromSmiles(smi)
        mw = Descriptors.MolWt(m)
        gcv = float(g.loc[ik, "value"]) if ik in g.index else model.predict_gcv(smi, mw)
        nc = float(sum(a.GetSymbol() == "C" for a in m.GetAtoms()))
        _COMP_CACHE[i] = dict(smiles=smi, mw=mw, gcv=gcv, nc=nc, ysi_ok=Chem.MolToSmiles(m) in _ysi_smiles())
    return _COMP_CACHE[i]


@lru_cache(maxsize=1)
def _k1() -> dict:
    """K1 per-temperature ln offsets fitted once on all designated viscosity rows of BASE rows_t2_bulk.csv (in-sample on T2)."""
    import pandas as pd
    from fuelblend.eval.cold_offset import fit_k1
    d = pd.read_csv(ROOT / "results" / "v5" / "BASE" / "rows_t2_bulk.csv")
    v = d[(d.property == "viscosity_kin_cSt") & (d.variant == "designated") & (d.status == "ok") & (~d.low_coverage.astype(bool))]
    return fit_k1(v.T_K, np.log(v.measured) - np.log(v.pred))


def _add_viscosity_k1(res: dict, T_K: float) -> None:
    """Variant only: no margin, only within K1_TOL_K of a fitted temperature, labelled in_sample_T2 (D3b)."""
    if res.get("value") is None or res.get("method") != "refutas_raw":      # K1 was fitted on raw Refutas only (no double correction)
        return
    for tf, off in _k1().items():
        if abs(T_K - tf) <= K1_TOL_K:
            res["variants"]["k1_cold_offset_in_sample"] = {
                "value": float(res["value"] * math.exp(off)), "offset_ln": float(off), "T_fit_K": float(tf), "label": "in_sample_T2",
                "note": "per-temperature ln offset fitted on designated raw-Refutas T2 viscosity rows (status ok, not low-coverage); in-sample, no margin"}
            return


ML_VISC_T_RANGE_K = (233.15, 393.15)        # delta_ml_viscosity README: correction defined over the 233-393 K training range


@lru_cache(maxsize=1)
def _ml_visc_margins() -> dict:
    """R-MLV-6 ML margins (ln, alpha 0.10) read from the committed delta_ml_viscosity readings."""
    import json
    d = json.loads((ROOT / "results" / "v5" / "BASE" / "delta_ml_viscosity" / "readings.json").read_text())["R-MLV-6"]["ml"]
    return {k: float(v["margin"]) for k, v in d.items()}


def _ml_margin_for(T_K: float) -> float:
    m = _ml_visc_margins()
    # same bins as the R-MLV-6 scoring (eval_lab_known_composition_v5.t_bin, engine_v4.viscosity_t_bin): 263.15 and 283.15 K go to the warmer bin
    return m["T_lt_263K"] if T_K < 263.15 else (m["T_263_283K"] if T_K < 283.15 else m["T_ge_283K"])


def _add_viscosity_ml(out: dict, ctx, pm: PropMap, T_K: float, composition_kind: str) -> None:
    """Variant ml_pair_residual (ln eta = ln eta_Refutas + sum x_i x_j g_ij, mole fractions).  Variant only: designated stays Refutas.
    R-MLV-1..4 passed, so a margin is attached in domain; reconstructed real fuels get the value with in_domain False, margin None."""
    res = out["viscosity_kin_cSt"]
    if res.get("value") is None:
        return
    dyn = out["viscosity_dyn_mPas"]
    base = {"method": "refutas+ml_pair_gn", "margin_unit": "ln", "training_scope": "lab binaries, 264 dev systems",
            "source": "results/v5/BASE/delta_ml_viscosity"}

    def _none(reason):
        v = dict(base, value=None, ln_correction=None, status="out_of_range", in_domain=False, margin=None, reason=reason,
                 warnings=[f"ml_pair_residual_out_of_range: {reason}"])
        res["variants"]["ml_pair_residual"] = v
        dyn["variants"]["ml_pair_residual"] = dict(v)

    if pm.substituted_mass > 0 or pm.unmapped_mass > 0 or pm.idx.size == 0:
        return _none("component without its own grid data (substituted or unmapped mass)")
    try:
        from fuelblend.ml import viscosity_pair as vp
        smi = [ctx.reg.smiles[i] for i in pm.idx]
        x = pm.w / np.array([ctx.mw(s) for s in smi])
        x = x / x.sum()
        dom = vp.domain_status(smi, T_K)
        if dom["in_domain"] and not (ML_VISC_T_RANGE_K[0] <= T_K <= ML_VISC_T_RANGE_K[1]):
            dom = dict(in_domain=False, reason=f"T outside the training range {ML_VISC_T_RANGE_K[0]}-{ML_VISC_T_RANGE_K[1]} K")
        if not dom["in_domain"]:
            return _none(dom["reason"])
        g = vp.ln_correction(smi, x, T_K)
    except Exception as e:                                      # model files missing or feature failure: variant unavailable
        return _none(f"ml model unavailable ({type(e).__name__})")
    known = composition_kind == "known"
    m = _ml_margin_for(T_K) if known else None
    w = [] if known else ["ml_pair_residual_reconstructed_fuel: R-MLV-5 transfer was mixed; in_domain False, no margin"]
    for r in (res, dyn):
        r["variants"]["ml_pair_residual"] = dict(base, value=float(r["value"] * math.exp(g)), ln_correction=float(g), status="ok",
                                                 in_domain=known, margin=m, warnings=list(w))


ML_DENSITY_READINGS = ROOT / "results" / "v5" / "BASE" / "delta_ml_density" / "readings.json"


def _add_density_ml(out: dict, ctx, pm: PropMap, T_K: float, composition_kind: str) -> None:
    """Variant ml_pair_excess_volume (V = V_ideal + sum x_i x_j g_ij).  Variant only: designated stays dens_idvol, because
    R-MLD-3 (ternaries) failed.  Margin (R-MLD-4, locked, kg/m3) only for known binaries in domain; 3+ components get the
    value with in_domain False and no margin."""
    res = out["density_kg_m3"]
    if res.get("value") is None:
        return
    base = {"method": "ideal_volume+ml_pair_excess_volume", "margin_unit": "kg/m3", "source": "results/v5/BASE/delta_ml_density"}

    def _none(reason):
        res["variants"]["ml_pair_excess_volume"] = dict(base, value=None, status="out_of_range", in_domain=False, margin=None,
                                                        warnings=[f"ml_pair_excess_volume_out_of_range: {reason}"])

    if pm.substituted_mass > 0 or pm.unmapped_mass > 0 or pm.idx.size == 0:
        return _none("component without its own grid data (substituted or unmapped mass)")
    try:
        import json as _json
        from fuelblend.ml import density_pair as dp
        smi = [ctx.reg.smiles[i] for i in pm.idx]
        dom = dp.domain_status(smi)
        if not dom["in_domain"]:
            return _none(dom["reason"])
        x = pm.w / np.array([ctx.mw(s) for s in smi])
        r = dp.density_ml(smi, x / x.sum(), T_K)
        m = float(_json.loads(ML_DENSITY_READINGS.read_text())["R-MLD-4"]["ml"]["locked_margin_kg_m3"])
    except Exception as e:                                      # model files missing or feature failure
        return _none(f"ml model unavailable ({type(e).__name__})")
    ok = composition_kind == "known" and len(smi) == 2
    w = [] if ok else ["ml_pair_excess_volume_no_margin: margin only for known binaries (R-MLD-3 ternaries failed; real fuels not tested)"]
    v = float(res["value"] * r["rho_ml"] / r["rho_ideal"])     # same correction ratio applied on the API's ideal value
    lo, hi = (v - m, v + m) if ok else (math.nan, math.nan)
    res["variants"]["ml_pair_excess_volume"] = dict(base, value=v, excess_volume_cm3_mol=r["excess_volume_cm3_mol"], status="ok",
                                                    in_domain=ok, margin=m if ok else None, lo=lo, hi=hi, warnings=w)


ML_DESIGNATED_NOTE ="designated_ml_pair_residual (user decision 2026-11-13, after R-MLV-1..4 passed)"


def _promote_ml_viscosity(out: dict) -> None:
    """User decision 2026-11-13 ("Make it main"): where the ml_pair_residual variant is ok and in domain (known composition,
    own grid data, 233-393 K, class pairs seen in training) it becomes the designated viscosity with its own conformal ln
    margin (R-MLV-6).  The Refutas result is kept as variants['refutas_raw'].  Elsewhere designated stays Refutas, with a
    warning naming why the ML value was not used."""
    from fuelblend.api import predict_blend as pb
    for name in ("viscosity_kin_cSt", "viscosity_dyn_mPas"):
        r = out[name]
        ml = r["variants"].get("ml_pair_residual")
        if r.get("value") is None or ml is None:
            continue
        if not (ml.get("status") == "ok" and ml.get("in_domain") and ml.get("margin") is not None):
            why = ml.get("reason") or "; ".join(ml.get("warnings") or []) or "not in the ML domain"
            r["warnings"].append(f"ml_not_designated: Refutas kept ({why})")
            continue
        r["variants"]["refutas_raw"] = {k: r.get(k) for k in ("value", "method", "lo", "hi", "margin", "margin_unit",
                                                              "in_domain", "interval")}
        v, m = ml["value"], ml["margin"]
        lo, hi = pb.interval_from_margin(v, m, "ln-ratio")
        r.update(value=v, method="refutas+ml_pair_gn", model_name="refutas+ml_pair_gn", margin=m, margin_unit="ln-ratio",
                 margin_source="results/v5/BASE/delta_ml_viscosity/readings.json (R-MLV-6)", lo=lo, hi=hi,
                 interval=[lo, hi], in_domain=True, uncertified=False,
                 tuple=(v, lo, hi, "refutas+ml_pair_gn", True))
        r["warnings"].append(ML_DESIGNATED_NOTE)


def _uncertified(res: dict, reason: str) -> dict:
    """No yaml margin: NaN lo/hi, in_domain False, and an explicit uncertified flag."""
    v = res["value"]
    res.update(model_name=res["method"], interval=[math.nan, math.nan], domain_reasons=[reason], uncertified=True,
               tuple=(v, math.nan, math.nan, res["method"], False))
    return res


def _mapped_input_weights(comps: list, pm: PropMap) -> dict:
    """input inchikey -> (mass fraction over the mapped mass, rep inchikey or None for own); unmapped inputs are absent."""
    sub = {a: b for a, b, _ in pm.substitutions}
    tot = 1.0 - pm.unmapped_mass
    out = {}
    for c in comps:
        if not c.is_hc:
            continue
        if c.reg_idx is not None:
            out[c.ik] = (c.w / tot, None)
        elif c.ik in sub:
            out[c.ik] = (c.w / tot, sub[c.ik])
    return out


FREEZING_MARGIN_KEY = "freezing_point_C.lab_liquidus.svl_verified"


def _freezing_margin_v5(prop: dict, iks: list, tmap: dict, composition_kind: str) -> bool:
    """v5 delta_freeze_margin: the liquidus margin (direction under, added to the prediction: safe_bound = hi) applies only to a known
    composition whose every mapped component has a verified Tm / dHfus and nothing was substituted or unmapped.  Returns True if applied."""
    e = _v5_margins().get(FREEZING_MARGIN_KEY)
    cov = prop["coverage"]
    if (e is None or prop.get("value") is None or composition_kind != "known" or cov["substituted_mass"] > 1e-12 or cov["unmapped_mass"] > 1e-12
            or not all((tmap.get(ik) or {}).get("source") == VERIFIED_SOURCE for ik in iks)):
        return False
    from fuelblend.api import predict_blend as pb
    margin, unit, direction, _ = pb._applied_margin([FREEZING_MARGIN_KEY], _v5_margins())
    if margin is None or direction != "under":
        return False
    lo, hi = pb.interval_from_margin(prop["value"], margin, unit)
    prop.update(model_name=prop["method"], interval=[lo, hi], lo=lo, hi=hi, margin=margin, margin_unit=unit, margin_direction=direction,
                safe_bound=hi, in_domain=True, domain_reasons=[], uncertified=False,
                tuple=(prop["value"], lo, hi, prop["method"], True))
    return True


def _freezing(prop: dict, comps, ctx, pm: PropMap, fuel_type, composition_kind: str = "known") -> None:
    from fuelblend.eval import freezing_f1 as ff
    from fuelblend.eval.composition_route import ideal_sle_C
    reg = ctx.reg
    iks = [reg.ik[int(i)] for i in pm.idx]
    x = pm.w / np.array([ctx.mw(reg.smiles[int(i)]) for i in pm.idx])
    x = x / x.sum()
    tmap = _tmap_freezing()                              # verified constants where present, else the Coutinho entry (tagged)

    def has_sle(ik):
        t = tmap.get(ik)
        return bool(t) and bool(np.isfinite(t.get("tm_K", np.nan)) and np.isfinite(t.get("dh_fus_J_mol", np.nan)))

    no_sle = float(sum(w for w, ik in zip(pm.w, iks) if not has_sle(ik)))
    prop["coverage"]["sle_missing_mapped_mass"] = no_sle   # mapped mass (own or representative) with no Tm / dHfus; counts in x
    v = ideal_sle_C(x, iks, tmap)
    v_old = ideal_sle_C(x, iks, _tmap())                  # all-Coutinho constants, the designated value before delta_freeze_margin
    prop["variants"]["coutinho_constants"] = {
        "value": float(v_old) if np.isfinite(v_old) else None,
        "label": "ideal SvL in mole fraction on the Coutinho Tm / dHfus map only (value before delta_freeze_margin); uncertified"}
    if not np.isfinite(v):
        _no_value(prop, "non_finite", "no component with Tm / dHfus constants above x = 1e-4")
        prop["uncertified"] = True
    else:
        prop.update(value=float(v), status="ok")
        if _freezing_margin_v5(prop, iks, tmap, composition_kind):
            prop["warnings"].append("freezing_margin_from_lab_liquidus_binaries (delta_freeze_margin): verified constants for every component")
        else:
            _uncertified(prop, "uncertified: no freezing margin applies (needs a known composition, verified Tm / dHfus for every component, "
                               "nothing substituted or unmapped)")
    prop["label"] = ("ideal SvL in mole fraction on verified Tm / dHfus where available (else the Coutinho entry); "
                     + ("lab-liquidus margin applied" if prop.get("uncertified") is False else "uncertified"))
    if prop.get("uncertified") is not False:
        prop["warnings"].append("freezing_uncertified")
    if no_sle > 0:
        prop["warnings"].append(f"components_without_sle_constants: mass {no_sle:.3f} is liquid solvent only")
    if str(fuel_type).lower() != "jet":                    # F1x variant: jet fuels with n-alkanes only
        return
    mw_map = _mapped_input_weights(comps, pm)
    mol_na, other = {}, 0.0
    for c in comps:
        if c.ik not in mw_map:
            continue
        w, rep = mw_map[c.ik]
        if c.subkind == "n_paraffin":
            mol_na[c.n_carbon] = mol_na.get(c.n_carbon, 0.0) + w / ff.nalkane_mw(c.n_carbon)
        else:
            smi = c.smiles if rep is None else reg.smiles[reg.ik2i[rep]]
            other += w / ctx.mw(smi)
    if not mol_na:
        return
    tot = sum(mol_na.values()) + other
    T, arg = ff.t_first(mol_na, tot, ff.coutinho_consts())
    if not np.isfinite(T):
        return
    a, b = ff.fit_f1x()
    tc = T - 273.15
    prop["variants"]["f1x_recalibrated"] = {
        "value": float(a + b * tc), "T_first_C": float(tc), "first_species_nC": arg, "a": a, "b": b,
        "logo_mae_C": ff.F1X_LOGO_MAE,
        "label": "a + b*T_first, complete Coutinho n-alkane table, a and b fitted on 49 S5 headline groups",
        "score_note": "only the out-of-sample pooled LOGO MAE (logo_mae_C, n=49 groups) may be quoted; uncertified"}


def _heat(prop: dict, ctx, pm: PropMap, fuel_type) -> None:
    from fuelblend.baselines.heat_of_combustion_v2 import mixture_heat_of_combustion_v2
    if str(fuel_type).lower() != "jet":
        prop["warnings"].append("heat_net_jet_only: pass fuel_type='jet'")
        _no_value(prop, "no_evaluation", "heat of combustion is evaluated for jet fuels only")
        return
    ci = [_comp_info(ctx, int(i)) for i in pm.idx]
    r = mixture_heat_of_combustion_v2(w_mass=pm.w, comp_gcv_MJ_kg=[c["gcv"] for c in ci], comp_mw_g_mol=[c["mw"] for c in ci],
                                      comp_n_carbon=[c["nc"] for c in ci], use_standard_jet_h=False)["ncv_MJ_kg"]
    if np.isfinite(r):
        prop.update(value=float(r), status="ok")
    else:
        _no_value(prop, "non_finite", "group-additivity heat is not finite")


def _smoke_ysi_variant(prop: dict, comps, ctx, pm: PropMap) -> None:
    """YSI route as a labelled variant only (no evaluable fuel in BASE: YSI coverage < 0.90 on all 39 T2 fuels); gated by own-coverage."""
    from fuelblend.baselines.sooting import mixture_mw, smoke_point_from_ysi, ysi_mix
    ci = [_comp_info(ctx, int(i)) for i in pm.idx]
    ok = np.array([c["ysi_ok"] for c in ci])
    cov_mapped = float(pm.w[ok].sum()) if ok.size else 0.0         # mapped mass (own or representative) with YSI
    prop["coverage"]["ysi_mapped_mass"] = cov_mapped
    tot = 1.0 - pm.unmapped_mass
    cov = float(sum(c.w for c in comps if c.is_hc and c.reg_idx is not None and _comp_info(ctx, c.reg_idx)["ysi_ok"]) / tot) if tot > 0 else 0.0
    prop["coverage"]["ysi_own_mass"] = cov                          # input compounds that themselves carry YSI; this gates
    if cov < SMOKE_YSI_MIN:
        prop["variants"]["ysi_route"] = {"value": None, "status": "no_evaluation", "ysi_own_mass": cov,
                                         "note": f"YSI own-coverage {cov:.3f} < {SMOKE_YSI_MIN}"}
        return
    w = pm.w[ok] / pm.w[ok].sum()
    smi = [c["smiles"] for c, o in zip(ci, ok) if o]
    mw = np.array([c["mw"] for c, o in zip(ci, ok) if o])
    x = (w / mw) / (w / mw).sum()
    prop["variants"]["ysi_route"] = {"value": float(smoke_point_from_ysi(ysi_mix(smi, x), mixture_mw(smi, x))), "status": "ok",
                                     "ysi_own_mass": cov, "note": "no measured skill; variant only"}


def _smoke(prop: dict, comps, ctx, pm: PropMap, arom_vol: dict) -> None:
    """Designated smoke point = the v4 AJFTD bulk regression on aromatics vol % (the model the design side uses; user decision 2026-11-13).
    Margin 7.8 mm from configs/uncertainty_v5.yaml only inside the domain (user decision 2026-11-13); uncertified and warns outside the regression domain (aromatics 8-23 vol %, result 19-25 mm)."""
    import warnings as _w
    from fuelblend.baselines.sooting import OutOfDomainWarning, smoke_point_from_bulk
    _smoke_ysi_variant(prop, comps, ctx, pm)
    a = arom_vol.get("value")
    if a is None:
        _no_value(prop, "no_evaluation", "aromatics vol % unavailable (low coverage on the density grid)")
        prop["method"] = "ajftd_bulk_regression"
        _uncertified(prop, "aromatics vol % unavailable")
        return
    with _w.catch_warnings(record=True) as rec:
        _w.simplefilter("always")
        v = float(smoke_point_from_bulk(float(a)))
    for r in rec:
        if issubclass(r.category, OutOfDomainWarning):
            prop["warnings"].append(f"out_of_domain: {r.message}")
    prop.update(value=v, status="ok", method="ajftd_bulk_regression")
    prop["inputs"] = {"aromatics_vol_pct": float(a)}
    cov = prop["coverage"]
    in_dom = not any(w.startswith("out_of_domain") for w in prop["warnings"])
    margins = _v5_margins()
    if SMOKE_V5_KEY in margins:                         # real-fuel margin present: it REPLACES the lab key, inside the class-mix range only
        e = margins[SMOKE_V5_KEY]
        gate_ok, gate_why = _smoke_class_mix_ok(comps)
    else:
        e, gate_ok, gate_why = margins.get(SMOKE_MARGIN_KEY), True, None
    m = None if e is None else e.get("applied_margin")
    if (not in_dom or cov["substituted_mass"] > 1e-12 or cov["unmapped_mass"] > 1e-12 or m is None
            or isinstance(m, str) or not math.isfinite(float(m)) or not gate_ok):
        _uncertified(prop, (gate_why or "smoke point has no valid yaml margin") if in_dom else "outside the smoke regression domain (aromatics 8-23 vol %, result 19-25 mm)")
        return
    m = float(m)
    lo, hi = v - m, v + m                   # over-prediction is unsafe for a minimum smoke point: lo is the safe bound
    prop.update(model_name=prop["method"], lo=lo, hi=hi, interval=[lo, hi], margin=m, margin_unit=e["unit"],
                margin_direction=e["direction"], safe_bound=lo, in_domain=True, domain_reasons=[],
                tuple=(v, lo, hi, prop["method"], True))
    prop["warnings"].append(SMOKE_MARGIN_NOTE_V5 if SMOKE_V5_KEY in margins else SMOKE_MARGIN_NOTE)


def _is_aromatic(smi: str) -> bool:
    return any(a.GetIsAromatic() for a in Chem.MolFromSmiles(smi).GetAtoms())


def _aromatics(out: dict, comps, ctx, P) -> None:
    hc = sum(c.w for c in comps if c.is_hc)                         # same HC-only basis as the volume share
    wt = 100.0 * sum(c.w for c in comps if c.is_hc and _is_aromatic(c.smiles)) / hc if hc > 0 else 0.0
    r = _new_prop("aromatics_wt_pct", "by_construction", map_for(comps, HAS_ANY_REGISTRY, ctx.reps))
    r.update(value=float(wt), status="ok")
    _uncertified(r, "aromatics content is by construction (no model, no margin)")
    out["aromatics_wt_pct"] = r
    pg = map_for(comps, P["HAS_GRID"], ctx.reps)
    v = _new_prop("aromatics_vol_pct", "by_construction", pg)
    if pg.unmapped_mass > LOW_COVERAGE_MASS or pg.idx.size == 0:
        _no_value(v, "low_coverage", f"unmapped mass {pg.unmapped_mass:.3f} > {LOW_COVERAGE_MASS}")
    else:
        _, rho = ctx.prov.lookup(pg.idx[None, :], 288.15)
        vol = pg.w / rho[0]
        flag = np.array([_is_aromatic(ctx.reg.smiles[int(i)]) for i in pg.idx])
        v.update(value=float(100.0 * vol[flag].sum() / vol.sum()), status="ok")
        v["warnings"].append("volume basis: pure-component grid density at 288.15 K")
        _uncertified(v, "aromatics content is by construction (no model, no margin)")
    out["aromatics_vol_pct"] = v


def _add_task3_props(out: dict, comps, ctx, fuel_type, T_K: float, composition_kind: str) -> None:
    P = predicates(ctx)
    iks = [c.ik for c in comps]
    pm = map_for(comps, HAS_ANY_REGISTRY, ctx.reps)             # the scoring route's property-independent mapping
    props = {k: _new_prop(k, m, pm) for k, m in (("freezing_point_C", "ideal_svl_mole_fraction"), ("heat_net_MJ_kg", "ga_comp_H"),
                                                  ("smoke_point_mm", "ajftd_bulk_regression"))}
    if pm.unmapped_mass > LOW_COVERAGE_MASS or pm.idx.size == 0:
        for r in props.values():
            _no_value(r, "low_coverage", f"unmapped mass {pm.unmapped_mass:.3f} > {LOW_COVERAGE_MASS}")
            r["uncertified"] = True
    else:
        _freezing(props["freezing_point_C"], comps, ctx, pm, fuel_type, composition_kind)
        _heat(props["heat_net_MJ_kg"], ctx, pm, fuel_type)
    out["freezing_point_C"] = props["freezing_point_C"]
    out["heat_net_MJ_kg"] = _finish(props["heat_net_MJ_kg"], HEAT_YAML_PROP, iks, T_K, fuel_type, composition_kind)
    out["smoke_point_mm"] = props["smoke_point_mm"]
    _aromatics(out, comps, ctx, P)
    if not (pm.unmapped_mass > LOW_COVERAGE_MASS or pm.idx.size == 0):
        _smoke(props["smoke_point_mm"], comps, ctx, pm, out["aromatics_vol_pct"])


# ====================================================================================================================
# indicative T2 interval (D6): descriptive, separate from lab margins, never read by design
# ====================================================================================================================
INDICATIVE_KIND = "empirical_T2_descriptive_not_held_out"
INDICATIVE_MAX_UNMAPPED = 0.10
_BULK_ROWS_NAME = {"heat_net_MJ_kg": "heat_of_combustion_net_MJ_kg"}      # API key -> BASE rows_t2_bulk.csv property
_LN_PROPS = ("viscosity_kin_cSt",)                                         # residuals in ln space; interval multiplicative


@lru_cache(maxsize=1)
def _indicative_tables() -> dict:
    """Frozen BASE T2 residual tables, built once.  {prop: [(T_K or nan, q90_abs, n_groups), ...]}.

    Flash: |pred_designated - measured| over the chosen row of each group in rows_t2_flash.csv.  Note: the API's designated flash
    is the LFL overlay (U1), so this interval comes from the v4-designated residuals; the label stays
    empirical_T2_descriptive_not_held_out.  Bulk: designated variant, status ok, not low_coverage, one row per group and BASE T bin;
    viscosity error is ln(pred/measured), the others pred - measured."""
    import pandas as pd
    base = ROOT / "results" / "v5" / "BASE"
    tabs: dict = {}
    f = pd.read_csv(base / "rows_t2_flash.csv")
    f = f[(f.status == "ok") & (~f.low_coverage.astype(bool)) & f.chosen.astype(bool)].drop_duplicates("group")
    e = np.abs(f.pred_designated.to_numpy(float) - f.measured.to_numpy(float))
    if len(f):
        tabs["flash_point_C"] = [(math.nan, float(np.quantile(e, 0.9)), int(len(f)))]
    b = pd.read_csv(base / "rows_t2_bulk.csv")
    b = b[(b.variant == "designated") & (b.status == "ok") & (~b.low_coverage.astype(bool))]
    for api_key in ("viscosity_kin_cSt", "density_kg_m3", "heat_net_MJ_kg", "smoke_point_mm"):
        d = b[b.property == _BULK_ROWS_NAME.get(api_key, api_key)]
        d = d.sort_values("chosen", ascending=False, kind="stable").drop_duplicates(["T_K", "group"])
        if api_key in _LN_PROPS:
            d = d[(d.measured > 0) & (d.pred > 0)]
            err = np.log(d.pred.to_numpy(float)) - np.log(d.measured.to_numpy(float))
        else:
            err = d.pred.to_numpy(float) - d.measured.to_numpy(float)
        if len(d) == 0:
            continue
        rows = []
        for tk, idx in d.reset_index(drop=True).groupby(d.T_K.fillna(-1.0).reset_index(drop=True)).groups.items():
            ee = np.abs(err[list(idx)])
            if len(ee) == 0:
                continue
            assert len(ee) == d.reset_index(drop=True).loc[list(idx)].group.nunique()     # one row per group in each bin
            rows.append((math.nan if tk < 0 else float(tk), float(np.quantile(ee, 0.9)), int(len(ee))))
        if rows:
            tabs[api_key] = rows
    return tabs


INDICATIVE_MIN_GROUPS = 5
INDICATIVE_MAX_BIN_GAP_K = 10.0


def _indicative_interval(prop: str, T_K: float, value: float) -> Optional[dict]:
    """90 % abs-error quantile of the frozen T2 residuals, nearest BASE T bin.  None when no table exists for ``prop``.
    Gating on fuel type and unmapped mass is done by the caller.  Never feeds lo/hi."""
    rows = _indicative_tables().get(prop)
    if not rows or value is None or not np.isfinite(value):
        return None
    withT = [r for r in rows if np.isfinite(r[0])]
    tb = None
    if withT:
        tb, q, n = min(withT, key=lambda r: abs(r[0] - T_K))
        if abs(tb - T_K) > INDICATIVE_MAX_BIN_GAP_K:      # the nearest evidence bin is too far in temperature
            return None
    else:
        q, n = rows[0][1:]
    if n < INDICATIVE_MIN_GROUPS:                          # a quantile of fewer than 5 groups is not an interval
        return None
    lo, hi = (value * math.exp(-q), value * math.exp(q)) if prop in _LN_PROPS else (value - q, value + q)
    res = {"lo": float(lo), "hi": float(hi), "q90_abs": float(q), "n_groups": int(n), "kind": INDICATIVE_KIND}
    if tb is not None:
        res["T_bin_K"] = float(tb)
    if prop == "flash_point_C":
        res["residual_basis"] = "v4_designated_flash (pred_designated); the API's designated flash is the LFL overlay"
    return res


def _attach_indicative(out: dict, fuel_type, T_K: float) -> None:
    if str(fuel_type).lower() != "jet":
        return
    for k, r in out.items():
        if not isinstance(r, dict) or "indicative_interval" not in r or r.get("value") is None:
            continue
        if r["coverage"]["unmapped_mass"] > INDICATIVE_MAX_UNMAPPED:
            continue
        if r.get("method") == "refutas_recalibrated":      # the table holds raw-Refutas residuals; do not attach it to a recalibrated value
            continue
        r["indicative_interval"] = _indicative_interval(k, T_K, r["value"])
