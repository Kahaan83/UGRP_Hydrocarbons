"""v5 BASE: one route from a class x carbon-number composition to predicted properties.

Replaces the two parallel mapping paths (G1 ``map_fuel`` and S2 ``rep_for``) with one property-independent mapping, and wraps the
property calls S3 used (``predict_comp``) behind a context object that can carry in-memory input overrides (LFL, Antoine) so that
variants never touch a file on disk.  Function bodies are moved from scripts/v5/{eval_ajftd_gcxgc_v5,s3_composition_vs_physics}.py
unchanged except where noted (the Refutas guard).
"""
from __future__ import annotations

import re
import warnings
from dataclasses import dataclass, field
from typing import Optional

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors

from fuelblend.baselines.flash_unifac import UnifacFallbackWarning, _default_registry, flash_point_unifac_detail
from fuelblend.thermo.coutinho_wax import R_GAS, get_pure_thermal_properties

RDLogger.DisableLog("rdApp.*")

COVERAGE_MIN = 0.90
NEAREST_ANY = ("cyclo_di", "indan_tetralin", "aromatic_di")

# source class label -> SubReps kind (None = no realistic representative exists)
_AJFTD_LIKE = {"n_paraffin": "n_paraffin", "iso_paraffin": "iso_paraffin", "monocycloparaffin": "cyclo_mono",
               "alkylbenzene": "alkylbenzene", "alkene": "alkene_alkyne", "dicycloparaffin": "cyclo_di",
               "cycloaromatic": "indan_tetralin", "diaromatic": "aromatic_di", "tricycloparaffin": None, "diene": None}
CLASS_MAPS: dict[str, dict[str, Optional[str]]] = {
    "ajftd": dict(_AJFTD_LIKE),
    "av3322": {"n-alkane": "n_paraffin", "iso-alkane": "iso_paraffin", "monocycloalkane": "cyclo_mono", "dicycloalkane": "cyclo_di",
               "alkylbenzene": "alkylbenzene", "cycloaromatic": "indan_tetralin", "diaromatic": "aromatic_di",
               "alkene": "alkene_alkyne", "tricycloalkane": None},
    # FuelLib classes are read from the GC bin label (see fuellib_class); these are the label stems
    "fuellib": {"Isoparaffin": "iso_paraffin", "Monocycloparaffin": "cyclo_mono", "Dicycloparaffin": "cyclo_di",
                "Tricycloparaffin": None, "Diaromatic": "aromatic_di", "Cycloaromatic": "indan_tetralin", "Alkene": "alkene_alkyne"},
    # FACE diesel tables carry the same hc_class labels as the AJFTD GCxGC composition (what G1 map_fuel reads)
    "face_diesel": dict(_AJFTD_LIKE),
}


def hc_class(m) -> str:
    if any(a.GetIsAromatic() for a in m.GetAtoms()):
        return "aromatic"
    if any(b.GetBondType() != Chem.BondType.SINGLE for b in m.GetBonds()):
        return "alkene_alkyne"
    if m.GetRingInfo().NumRings() > 0:
        return "cycloalkane"
    return "iso_paraffin" if any(a.GetSymbol() == "C" and a.GetDegree() > 2 for a in m.GetAtoms()) else "n_paraffin"


def fuellib_class(gc_bin: str) -> Optional[str]:
    """FuelLib GC bin label -> SubReps kind (None = no representative).  Same rules as scripts/v5/s2_t2_flash.fuellib_kind."""
    b = gc_bin
    if re.match(r"n-C\d+$", b):
        return "n_paraffin"
    for key, kind in CLASS_MAPS["fuellib"].items():
        if key in b:
            return kind
    if b == "Toluene" or re.match(r"C\d+-Benzene$", b):
        return "alkylbenzene"
    raise ValueError(f"unrecognised FuelLib gc_bin {b!r}")


def kind_bins(source: str, df: pd.DataFrame) -> pd.DataFrame:
    """One fuel's rows in a source's own layout -> bins (hc_class = SubReps kind or 'no_rep:<label>', n_carbon, mass_pct)."""
    if source == "av3322":
        lab, cn, m = df["group"], df["carbon_number"], df["mass_pct"]
        kinds = [CLASS_MAPS["av3322"][g] for g in lab]
    elif source == "fuellib":
        lab, cn, m = df["gc_bin"], df["n_carbon"], df["fraction"]
        kinds = [fuellib_class(g) for g in lab]
    elif source in ("ajftd", "face_diesel"):
        lab, cn, m = df["hc_class"], df["n_carbon"], df["mass_pct"]
        kinds = [CLASS_MAPS[source][g] for g in lab]
    else:
        raise ValueError(source)
    hc = [k if k is not None else f"no_rep:{l}" for k, l in zip(kinds, lab)]
    return pd.DataFrame({"hc_class": hc, "n_carbon": np.asarray(cn, int), "mass_pct": np.asarray(m, float)})


@dataclass
class MappedFuel:
    idx: np.ndarray
    w: np.ndarray
    unmapped_fraction: float
    unmapped_by_class: dict

    @property
    def low_coverage(self) -> bool:
        return self.unmapped_fraction > 1.0 - COVERAGE_MIN


def map_bins(bins: pd.DataFrame, reps) -> MappedFuel:
    """bins: (hc_class = SubReps kind, n_carbon, mass_pct) of one fuel.  A bin is unmapped only when ``reps.get`` returns None
    (property-independent: nothing about Antoine, LFL or the property being scored enters here)."""
    agg: dict[int, float] = {}
    lost: dict[str, float] = {}
    tot = float(bins.mass_pct.sum())
    for r in bins.itertuples():
        kind = str(r.hc_class)
        i = None
        if not kind.startswith("no_rep"):
            i = reps.get(kind, int(r.n_carbon), max_gap=99 if kind in NEAREST_ANY else 3)
        if i is None:
            lost[kind] = lost.get(kind, 0.0) + r.mass_pct / tot
        else:
            agg[i] = agg.get(i, 0.0) + r.mass_pct
    mapped = sum(agg.values())
    idx = np.array(sorted(agg), dtype=np.int64)
    w = np.array([agg[i] for i in idx], float)
    w = w / w.sum() if w.size and w.sum() > 0 else w
    return MappedFuel(idx=idx, w=w, unmapped_fraction=1.0 - mapped / tot, unmapped_by_class=lost)


def refutas_nu_safe(provider, idx, w, T_K) -> tuple[float, float]:
    """Kinematic Refutas viscosity (cSt) and ideal-volume density (kg/m3).  Components with non-positive mass are dropped first
    (a zero-mass propane representative has nu + 0.8 < 1 and its NaN would poison 0 * NaN in the VBN sum); if a remaining component
    still has nu_i + 0.8 <= 1 the viscosity is NaN, never a silently wrong number."""
    idx, w = np.asarray(idx, np.int64), np.asarray(w, float)
    keep = w > 0
    idx, w = idx[keep], w[keep] / w[keep].sum()
    ln_eta, rho = provider.lookup(idx[None, :], T_K)
    ln_eta, rho = ln_eta[0], rho[0]
    rho_mix = float(1.0 / (w / rho).sum())
    nu_i = 1000.0 * np.exp(ln_eta) / rho
    if not np.all(np.isfinite(nu_i)) or np.any(nu_i + 0.8 <= 1.0):
        return float("nan"), rho_mix
    vbn_i = 14.534 * np.log(np.log(nu_i + 0.8)) + 10.975
    vbn = float((w * vbn_i).sum())
    return float(np.exp(np.exp((vbn - 10.975) / 14.534)) - 0.8), rho_mix


def ideal_sle_C(x, iks, tmap) -> float:
    best = np.nan
    for xi, ik in zip(x, iks):
        t = tmap.get(ik)
        if t and np.isfinite(t["tm_K"]) and np.isfinite(t["dh_fus_J_mol"]) and xi > 1e-4:
            inv = 1.0 / t["tm_K"] - R_GAS / t["dh_fus_J_mol"] * np.log(xi)
            if inv > 0:
                T = 1.0 / inv - 273.15
                best = T if not np.isfinite(best) else max(best, T)
    return best


def build_tmap(reg, idx) -> dict:
    """{inchikey: {tm_K, dh_fus_J_mol, ...}} for the registry compounds at ``idx`` (same call S3 made)."""
    idx = sorted(set(int(i) for i in idx))
    info = pd.DataFrame(index=[reg.ik[i] for i in idx])
    info["smiles"] = [reg.smiles[i] for i in idx]
    mols = [Chem.MolFromSmiles(s) for s in info.smiles]
    info["cls"] = [hc_class(m) for m in mols]
    info["n_carbon"] = [float(sum(a.GetSymbol() == "C" for a in m.GetAtoms())) for m in mols]
    return get_pure_thermal_properties(info)


@dataclass
class RouteCtx:
    """reg: CompoundRegistry; prov: PurePropertyProvider; tmap: {inchikey: {tm_K, dh_fus_J_mol}}.
    Overrides are keyed by InChIKey and applied to an in-memory copy of the flash registry only."""
    reg: object
    prov: object
    tmap: dict = field(default_factory=dict)
    lfl_override: Optional[dict] = None
    antoine_override: Optional[dict] = None
    flash_reg: pd.DataFrame = field(init=False, repr=False, default=None)

    def __post_init__(self):
        fr = _default_registry().copy()
        for ik, v in (self.lfl_override or {}).items():
            fr.loc[fr.inchikey == ik, "lfl_vol_pct"] = float(v)
        for ik, abc in (self.antoine_override or {}).items():
            for col, val in zip(("antoine_A", "antoine_B", "antoine_C"), abc):
                fr.loc[fr.inchikey == ik, col] = float(val)
        self.flash_reg = fr
        self._mw: dict[str, float] = {}

    def mw(self, smi: str) -> float:
        if smi not in self._mw:
            self._mw[smi] = Descriptors.MolWt(Chem.MolFromSmiles(smi))
        return self._mw[smi]


def predict_properties(ctx: RouteCtx, fuel: MappedFuel, targets: pd.DataFrame) -> pd.DataFrame:
    """targets has columns property, T_K.  Adds pred and status; status 'ok' implies a finite prediction.
    flash: ideal Le Chatelier; viscosity: Refutas; density: ideal volume; freezing: ideal SvL."""
    keep = fuel.w > 0
    idx, w = fuel.idx[keep], fuel.w[keep]
    w = w / w.sum()
    iks = [ctx.reg.ik[i] for i in idx]
    smi = [ctx.reg.smiles[i] for i in idx]
    x = w / np.array([ctx.mw(s) for s in smi])
    x = x / x.sum()
    preds, stats = [], []
    for r in targets.itertuples():
        pred, status = np.nan, "ok"
        if r.property == "flash_point_C":
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", UnifacFallbackWarning)
                d = flash_point_unifac_detail(smi, x, registry=ctx.flash_reg, use_unifac=False)
            pred, status = d["flash_C"], d["status"]
        elif r.property in ("viscosity_kin_cSt", "density_kg_m3"):
            nu, rho = refutas_nu_safe(ctx.prov, idx, w, float(r.T_K))
            pred = nu if r.property == "viscosity_kin_cSt" else rho
        elif r.property == "freezing_point_C":
            pred = ideal_sle_C(x, iks, ctx.tmap)
        else:
            status = "unsupported_property"
        pred = float(pred) if pred is not None else np.nan
        if not np.isfinite(pred) and status == "ok":
            status = "non_finite"
        preds.append(pred)
        stats.append(status)
    out = targets.reset_index(drop=True).copy()
    out["pred"], out["status"] = preds, stats
    return out
