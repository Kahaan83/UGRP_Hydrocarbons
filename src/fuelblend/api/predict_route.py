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
    _add_viscosity_k1(out["viscosity_kin_cSt"], T_K)
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
            res["variants"]["k1_cold_offset_loso"] = {
                "value": float(res["value"] * math.exp(off)), "offset_ln": float(off), "T_fit_K": float(tf), "label": "in_sample_T2",
                "note": "per-temperature ln offset fitted on designated raw-Refutas T2 viscosity rows (status ok, not low-coverage); in-sample, no margin"}
            return


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


def _freezing(prop: dict, comps, ctx, pm: PropMap, fuel_type) -> None:
    from fuelblend.eval import freezing_f1 as ff
    from fuelblend.eval.composition_route import ideal_sle_C
    reg = ctx.reg
    iks = [reg.ik[int(i)] for i in pm.idx]
    x = pm.w / np.array([ctx.mw(reg.smiles[int(i)]) for i in pm.idx])
    x = x / x.sum()
    tmap = _tmap()

    def has_sle(ik):
        t = tmap.get(ik)
        return bool(t) and bool(np.isfinite(t.get("tm_K", np.nan)) and np.isfinite(t.get("dh_fus_J_mol", np.nan)))

    no_sle = float(sum(w for w, ik in zip(pm.w, iks) if not has_sle(ik)))
    prop["coverage"]["sle_missing_mapped_mass"] = no_sle   # mapped mass (own or representative) with no Tm / dHfus; counts in x
    v = ideal_sle_C(x, iks, tmap)
    if not np.isfinite(v):
        _no_value(prop, "non_finite", "no component with Tm / dHfus constants above x = 1e-4")
        prop["uncertified"] = True
    else:
        prop.update(value=float(v), status="ok")
        _uncertified(prop, "uncertified: no skill shown over a constant; ideal SvL bias about +8 C; no yaml margin")
    prop["label"] = "ideal SvL in mole fraction; no skill shown, bias about +8 C; uncertified"
    prop["warnings"].append("freezing_uncertified_no_skill_shown")
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


def _smoke(prop: dict, comps, ctx, pm: PropMap) -> None:
    from fuelblend.baselines.sooting import mixture_mw, smoke_point_from_ysi, ysi_mix
    ci = [_comp_info(ctx, int(i)) for i in pm.idx]
    ok = np.array([c["ysi_ok"] for c in ci])
    cov_mapped = float(pm.w[ok].sum()) if ok.size else 0.0         # mapped mass (own or representative) with YSI
    prop["coverage"]["ysi_mapped_mass"] = cov_mapped
    tot = 1.0 - pm.unmapped_mass
    cov = float(sum(c.w for c in comps if c.is_hc and c.reg_idx is not None and _comp_info(ctx, c.reg_idx)["ysi_ok"]) / tot) if tot > 0 else 0.0
    prop["coverage"]["ysi_own_mass"] = cov                          # input compounds that themselves carry YSI; this gates
    if cov < SMOKE_YSI_MIN:
        prop.update(value=None, status="no_evaluation", method="ysi_route")
        prop["warnings"].append(f"no_evaluation: YSI own-coverage {cov:.3f} < {SMOKE_YSI_MIN}")
        _uncertified(prop, f"YSI coverage {cov:.3f} < {SMOKE_YSI_MIN}")
        return
    w = pm.w[ok] / pm.w[ok].sum()
    smi = [c["smiles"] for c, o in zip(ci, ok) if o]
    mw = np.array([c["mw"] for c, o in zip(ci, ok) if o])
    x = (w / mw) / (w / mw).sum()
    prop.update(value=float(smoke_point_from_ysi(ysi_mix(smi, x), mixture_mw(smi, x))), status="ok")
    _uncertified(prop, "smoke point has no valid yaml margin")


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
                                                  ("smoke_point_mm", "ysi_route"))}
    if pm.unmapped_mass > LOW_COVERAGE_MASS or pm.idx.size == 0:
        for r in props.values():
            _no_value(r, "low_coverage", f"unmapped mass {pm.unmapped_mass:.3f} > {LOW_COVERAGE_MASS}")
            r["uncertified"] = True
    else:
        _freezing(props["freezing_point_C"], comps, ctx, pm, fuel_type)
        _heat(props["heat_net_MJ_kg"], ctx, pm, fuel_type)
        _smoke(props["smoke_point_mm"], comps, ctx, pm)
    out["freezing_point_C"] = props["freezing_point_C"]
    out["heat_net_MJ_kg"] = _finish(props["heat_net_MJ_kg"], HEAT_YAML_PROP, iks, T_K, fuel_type, composition_kind)
    out["smoke_point_mm"] = props["smoke_point_mm"]
    _aromatics(out, comps, ctx, P)


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
    tabs["flash_point_C"] = [(math.nan, float(np.quantile(e, 0.9)), int(len(f)))]
    b = pd.read_csv(base / "rows_t2_bulk.csv")
    b = b[(b.variant == "designated") & (b.status == "ok") & (~b.low_coverage.astype(bool))]
    for api_key in ("viscosity_kin_cSt", "density_kg_m3", "heat_net_MJ_kg", "smoke_point_mm"):
        d = b[b.property == _BULK_ROWS_NAME.get(api_key, api_key)].drop_duplicates(["T_K", "group"])
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
            rows.append((math.nan if tk < 0 else float(tk), float(np.quantile(ee, 0.9)), int(len(ee))))
        tabs[api_key] = rows
    return tabs


def _indicative_interval(prop: str, T_K: float, value: float) -> Optional[dict]:
    """90 % abs-error quantile of the frozen T2 residuals, nearest BASE T bin.  None when no table exists for ``prop``.
    Gating on fuel type and unmapped mass is done by the caller.  Never feeds lo/hi."""
    rows = _indicative_tables().get(prop)
    if not rows or value is None or not np.isfinite(value):
        return None
    withT = [r for r in rows if np.isfinite(r[0])]
    q, n = (min(withT, key=lambda r: abs(r[0] - T_K))[1:] if withT else rows[0][1:])
    lo, hi = (value * math.exp(-q), value * math.exp(q)) if prop in _LN_PROPS else (value - q, value + q)
    return {"lo": float(lo), "hi": float(hi), "q90_abs": float(q), "n_groups": int(n), "kind": INDICATIVE_KIND}


def _attach_indicative(out: dict, fuel_type, T_K: float) -> None:
    if str(fuel_type).lower() != "jet":
        return
    for k, r in out.items():
        if not isinstance(r, dict) or "indicative_interval" not in r or r.get("value") is None:
            continue
        if r["coverage"]["unmapped_mass"] > INDICATIVE_MAX_UNMAPPED:
            continue
        r["indicative_interval"] = _indicative_interval(k, T_K, r["value"])
