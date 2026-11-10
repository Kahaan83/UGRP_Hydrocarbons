"""v4 review fix 3: residuals of pure-component physics on KNOWN-composition lab mixtures.

Why: the CRC647 margins (results/v4/U) mix reconstruction error with physics error.  Lab mixtures have measured
compositions, so their residuals isolate the physics.  Two datasets:
  * density   : ThermoML liquid mass density of hydrocarbon mixtures (n_comp >= 2), ideal volume mixing
                rho = 1 / sum(w_i / rho_i(T)), pure rho_i(T) from eval_bulk_props_v4.PureDensity TIER 1 ONLY (measured fits).
  * viscosity : build_viscosity_leaderboard_v4.build_lab selection (hydrocarbon binary/multi dynamic viscosity,
                non-conflicting, measurements_dedup.csv) WITHOUT its LAB_T_RANGE filter; raw Refutas + ideal-volume density.
Residual convention: predicted - measured (density kg/m3; viscosity ln(pred/meas), = minus build_lab's ln(meas/pred)).

Lab split: ALL structurally-eligible system ids of both datasets (sorted, before any pure-property validity filter and
before any residual exists) -> 20 % locked with np.random.default_rng(20261111) -> splits/v4/lab_systems_splits.json.
Margins use dev systems only (scripts/v4/build_uncertainty_v4.py::lab_known).

Run:  .venv/Scripts/python.exe scripts/v4/eval_lab_known_composition.py
"""
from __future__ import annotations

import functools
import hashlib
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

THERMOML = ROOT / "data" / "processed" / "v2" / "thermoml_clean_hc.csv"
MEAS = ROOT / "data" / "unified" / "measurements_dedup.csv"
SPLIT_PATH = ROOT / "splits" / "v4" / "lab_systems_splits.json"
OUT = ROOT / "results" / "v4" / "L"
SEED = 20261111
LOCK_FRAC = 0.20
T_REF_K = 288.71            # same reference as eval_bulk_props_v4.PureDensity
BINS = ("T_lt_263K", "T_263_283K", "T_ge_283K")


def t_bin(T_K: float) -> str:
    """263.15 and 283.15 K belong to the WARMER bin."""
    if T_K < 263.15:
        return "T_lt_263K"
    if T_K < 283.15:
        return "T_263_283K"
    return "T_ge_283K"


def _load(name: str, rel: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / rel)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


def system_id(iks) -> str:
    return hashlib.md5("|".join(sorted(iks)).encode()).hexdigest()[:12]


# ------------------------------------------------------------------------------------------------ density
@functools.lru_cache(maxsize=None)
def _inchi_info(inchi: str):
    """-> (inchikey, molar mass, is_hydrocarbon).  Hydrocarbon: only C/H atoms, neutral, >= 1 C."""
    m = Chem.MolFromInchi(inchi)
    if m is None:
        return None
    syms = {a.GetSymbol() for a in m.GetAtoms()}
    hc = syms <= {"C", "H"} and "C" in syms and all(a.GetFormalCharge() == 0 for a in m.GetAtoms())
    return Chem.inchi.InchiToInchiKey(inchi), float(Descriptors.MolWt(m)), bool(hc)


def prepare_density_records(df: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    """Structural selection only (no pure-density check): -> one de-duplicated row per (system, T, composition)."""
    t = df[df.prop.str.contains("Mass density") & (df.prop_phase == "Liquid") & (df.n_comp >= 2)
           & (df.P_kPa.isna() | (df.P_kPa <= 110))]
    st = dict(raw_rows=int(len(t)), skipped_non_hc=0, skipped_comp_type=0, skipped_lt2_components=0)
    recs = []
    for r in t.itertuples():
        if r.comp_type not in ("mass fraction", "mole fraction"):
            st["skipped_comp_type"] += 1
            continue
        inchis, vals = json.loads(r.inchis), np.array(json.loads(r.comp_values), float)
        info = [_inchi_info(i) for i in inchis]
        if any(i is None or not i[2] for i in info):
            st["skipped_non_hc"] += 1
            continue
        keep = vals > 0                                           # a 0.0 endpoint is the pure liquid, not a mixture
        if keep.sum() < 2 or not np.all(np.isfinite(vals)) or not np.isfinite(r.T_K):
            st["skipped_lt2_components"] += 1
            continue
        iks = [info[i][0] for i in range(len(info)) if keep[i]]
        mw = np.array([info[i][1] for i in range(len(info)) if keep[i]])
        v = vals[keep]
        w = v if r.comp_type == "mass fraction" else v * mw / np.sum(v * mw)
        w = w / w.sum()
        order = np.argsort(iks)
        iks, w = [iks[i] for i in order], w[order]
        recs.append(dict(system_id=system_id(iks), iks=tuple(iks), w=tuple(np.round(w, 3)), T_K=float(r.T_K),
                         T_key=round(float(r.T_K) * 2) / 2, value=float(r.value), w_full=tuple(w)))
    d = pd.DataFrame(recs)
    if d.empty:
        return d, st
    g = d.groupby(["system_id", "iks", "w", "T_key"], sort=True)
    out = g.agg(T_K=("T_K", "first"), measured=("value", "median"), n_obs=("value", "size"), w_full=("w_full", "first")).reset_index()
    out["n_comp"] = out["iks"].map(len)
    st["after_dedup_records"] = int(len(out))
    return out, st


def build_density_rows(df: pd.DataFrame, pdm, split_role: dict | None = None) -> tuple[pd.DataFrame, dict]:
    """Ideal-volume density residuals.  A system is EXCLUDED (and counted) when any of its components has no tier-1
    (measured) pure-density fit at all; a record is dropped (and counted) when a component's tier-1 fit does not cover
    its T (+-30 K extrapolation, PureDensity's own rule).  Nothing is filled from a class regression."""
    d, st = prepare_density_records(df)
    st = dict(st)
    st.update(n_excluded_no_pure_density=0, n_records_excluded_no_pure_density=0, n_records_dropped_T_outside_tier1=0,
              n_systems_structural=int(d.system_id.nunique()) if len(d) else 0)
    if d.empty:
        return pd.DataFrame(columns=["system_id", "split_role", "T_K", "t_bin", "n_comp", "predicted", "measured", "resid"]), st
    bad = {s for s, iks in zip(d.system_id, d.iks) if any(k not in pdm.t1 for k in iks)}
    st["n_excluded_no_pure_density"] = len(bad)
    st["n_records_excluded_no_pure_density"] = int(d.system_id.isin(bad).sum())
    rows = []
    for r in d[~d.system_id.isin(bad)].itertuples():
        rho = []
        for k in r.iks:
            val, tier = pdm.rho(k, "unused", 10, r.T_K)
            rho.append(val if tier == 1 else np.nan)
        if not np.all(np.isfinite(rho)):
            st["n_records_dropped_T_outside_tier1"] += 1
            continue
        pred = 1.0 / float(np.sum(np.array(r.w_full) / np.array(rho)))
        rows.append(dict(system_id=r.system_id, split_role=(split_role or {}).get(r.system_id, "unassigned"), T_K=r.T_K,
                         t_bin=t_bin(r.T_K), n_comp=r.n_comp, predicted=pred, measured=r.measured, resid=pred - r.measured))
    out = pd.DataFrame(rows, columns=["system_id", "split_role", "T_K", "t_bin", "n_comp", "predicted", "measured", "resid"])
    st["n_systems_scored"] = int(out.system_id.nunique())
    st["n_systems_no_valid_T"] = len(set(d.system_id) - bad - set(out.system_id))
    return out, st


# ---------------------------------------------------------------------------------------------- viscosity
def viscosity_structural_ids(meas: pd.DataFrame) -> list[str]:
    w = meas[meas.kind.isin(["binary", "multi"]) & meas.is_hc.astype(bool) & (meas.property == "viscosity_dyn_mPas")
             & ~meas.conflict.astype(bool)]
    return sorted(set(w.system_id))


# ------------------------------------------------------------------------------------------------- split
def make_split(ids) -> dict:
    ids = sorted(set(ids))
    n_lock = int(round(LOCK_FRAC * len(ids)))
    rng = np.random.default_rng(SEED)
    locked = sorted(rng.choice(np.array(ids), size=n_lock, replace=False).tolist())
    ls = set(locked)
    return dict(seed=SEED, locked_fraction=LOCK_FRAC, n_systems=len(ids), n_locked=n_lock,
                definition=("union of structurally-eligible system ids of the ThermoML density mixtures and the build_lab "
                            "viscosity mixtures (before any pure-property validity filter and before any residual); sorted; "
                            "locked = rng.choice(ids, round(0.2 n)) with np.random.default_rng(20261111)"),
                dev=[i for i in ids if i not in ls], locked=locked)


def write_split(ids) -> dict:
    sp = make_split(ids)
    if SPLIT_PATH.exists():
        old = json.loads(SPLIT_PATH.read_text())
        if old != sp:
            raise SystemExit(f"{SPLIT_PATH} exists and differs from the recomputed split; refusing to overwrite (inputs changed?)")
        return old
    SPLIT_PATH.parent.mkdir(parents=True, exist_ok=True)
    SPLIT_PATH.write_text(json.dumps(sp, indent=1) + "\n")
    return sp


# --------------------------------------------------------------------------------------------------- main
def summarize(df: pd.DataFrame, unit: str) -> list[str]:
    L = ["| slice | split | records | systems | MAE | bias (pred - meas) |", "|---|---|---|---|---|---|"]
    for name, sel in [("all T", np.ones(len(df), bool))] + [(b, (df.t_bin == b).to_numpy()) for b in BINS]:
        for role in ("dev", "locked", "all"):
            s = df[sel & ((df.split_role == role).to_numpy() if role != "all" else True)]
            if len(s):
                L.append(f"| {name} | {role} | {len(s)} | {s.system_id.nunique()} | {s.resid.abs().mean():.4g} {unit} | {s.resid.mean():+.4g} {unit} |")
            else:
                L.append(f"| {name} | {role} | 0 | 0 | n/a | n/a |")
    return L


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    thermo = pd.read_csv(THERMOML, low_memory=False)
    meas = pd.read_csv(MEAS, low_memory=False)

    # 1. split FIRST (ids only; no residual exists yet)
    dens_struct, dens_st0 = prepare_density_records(thermo)
    visc_ids = viscosity_structural_ids(meas)
    sp = write_split(set(dens_struct.system_id) | set(visc_ids))
    role = {**{s: "dev" for s in sp["dev"]}, **{s: "locked" for s in sp["locked"]}}
    print(f"split: {sp['n_systems']} systems, {sp['n_locked']} locked -> {SPLIT_PATH.relative_to(ROOT)}", flush=True)

    # 2. density
    ev = _load("eval_bulk_props_v4", "scripts/v4/eval_bulk_props_v4.py")
    pdm = ev.PureDensity(pd.DataFrame(columns=["inchikey_norm", "cls", "n_carbon"]))
    dens, dst = build_density_rows(thermo, pdm, role)
    dens.to_csv(OUT / "density_lab_residuals.csv", index=False, float_format="%.8g")

    # 3. viscosity (build_lab, LAB_T_RANGE filter removed)
    lb = _load("build_viscosity_leaderboard_v4", "scripts/v4/build_viscosity_leaderboard_v4.py")
    from fuelblend.mixture.cold_flow import compound_meta
    from fuelblend.mixture.common import CompoundRegistry, PurePropertyProvider
    reg = CompoundRegistry(ROOT / "data" / "unified" / "compounds.csv")
    _, cls, nc = compound_meta(ROOT / "data" / "unified" / "compounds.csv")
    provider = PurePropertyProvider(reg, ROOT / "data" / "cache" / "pure_grid.npz", verbose=False)
    lb.LAB_T_RANGE = (0.0, 1.0e9)
    lab = lb.build_lab(reg, provider, cls, nc)
    visc = pd.DataFrame(dict(system_id=lab.system_id, split_role=lab.system_id.map(role).fillna("unassigned"), T_K=lab.T_K,
                             t_bin=[t_bin(x) for x in lab.T_K], n_comp=lab.n_comp, predicted=lab.raw_refutas_mPas,
                             measured=lab.measured_mPas, resid=-lab.resid))      # build_lab resid = ln meas - ln pred
    assert np.allclose(visc.resid, np.log(visc.predicted) - np.log(visc.measured))
    visc.to_csv(OUT / "viscosity_lab_residuals.csv", index=False, float_format="%.8g")
    vst = dict(n_systems_structural=len(visc_ids), n_systems_scored=int(visc.system_id.nunique()),
               n_systems_excluded_no_valid_pure_property=len(visc_ids) - int(visc.system_id.nunique()))

    # 4. README
    assert not (set(dens.system_id) | set(visc.system_id)) - set(role), "residual system missing from split"
    L = ["# v4 review fix 3: known-composition lab residuals (density and raw Refutas viscosity)", "",
         "Generated by `scripts/v4/eval_lab_known_composition.py`. Residual = predicted - measured. Lab mixtures have measured "
         "compositions, so these residuals isolate the physics from the CRC647 reconstruction error.", "",
         "## Lab split", f"- {sp['n_systems']} structurally-eligible systems (union of density and viscosity lab systems), "
         f"{sp['n_locked']} locked (20 %, `np.random.default_rng({SEED})`), written to `splits/v4/lab_systems_splits.json` BEFORE any residual was computed; "
         "the script refuses to overwrite a differing split.",
         "- Margins (`configs/uncertainty_v4.yaml`) use dev systems only; locked systems feed a counts-only coverage check.",
         "- The independent unit for margins is the SYSTEM (distinct InChIKey set), as for the existing lab flash entries.", "",
         "## Density: ideal volume mixing", "",
         f"- Source `data/processed/v2/thermoml_clean_hc.csv`: liquid `Mass density`, n_comp >= 2, P missing or <= 110 kPa; {dst['raw_rows']} raw rows; "
         f"skipped: {dst['skipped_non_hc']} non-hydrocarbon, {dst['skipped_comp_type']} with a composition basis other than mass/mole fraction, "
         f"{dst['skipped_lt2_components']} with < 2 components of non-zero fraction; {dst['after_dedup_records']} records after de-duplication "
         "(median over system, T to 0.5 K, composition to 3 decimals).",
         f"- Pure densities: `eval_bulk_props_v4.PureDensity` TIER 1 only (measured ThermoML linear fits). **Excluded for missing pure density: "
         f"{dst['n_excluded_no_pure_density']} systems ({dst['n_records_excluded_no_pure_density']} records)** (a component has no tier-1 fit; "
         f"nothing is filled from a class regression). A further {dst['n_records_dropped_T_outside_tier1']} records were dropped because a "
         f"tier-1 fit does not cover the record's T (+-30 K), leaving {dst['n_systems_scored']} scored systems "
         f"({dst['n_systems_no_valid_T']} further systems lost entirely to that T rule).",
         f"- Scored: {len(dens)} records, {dens.system_id.nunique()} systems.", ""]
    L += summarize(dens, "kg/m3")
    L += ["", "## Viscosity: raw Refutas (ideal-volume density), ln-ratio", "",
          "- Selection = `build_viscosity_leaderboard_v4.build_lab` (hydrocarbon binary/multi dynamic viscosity from `measurements_dedup.csv`, "
          "no conflict flag, every component with valid pure eta/rho) with its `LAB_T_RANGE` (253-340 K) filter removed so that the cold bin is populated if data exist.",
          f"- {vst['n_systems_structural']} structurally eligible systems, {vst['n_systems_scored']} scored; "
          f"{vst['n_systems_excluded_no_valid_pure_property']} excluded because a component lacks valid pure viscosity/density "
          f"(for viscosity this is the registry/provider validity of `build_lab`, not the ThermoML tier-1 density fit).",
          f"- Scored: {len(visc)} records. `predicted`/`measured` are dynamic viscosity in mPa s; the ln-ratio is identical for kinematic viscosity at the same density.", ""]
    L += summarize(visc, "ln")
    L += ["", "### Temperature bins", "- `t_bin`: `T_lt_263K` (T < 263.15 K), `T_263_283K` (263.15 <= T < 283.15), `T_ge_283K`; the edges belong to the warmer bin.",
          "- The cold bin is thin (see the systems column above): when fewer than 9 independent dev systems exist, the conformal margin is `inf` and "
          "no substitute is given.", ""]
    L += ["## Caveats",
          "- Lab data are mostly binaries of n-alkanes/aromatics at 293-323 K; transfer to 100-component jet fuels is not demonstrated.",
          "- Density bias comes from ideal volume mixing (excess volumes neglected); the over/under margins are one-sided on that bias.",
          "- Raw Refutas is validated here on lab binaries, whereas Task 11 found it collapses on CRC647 fuels at -40 C; the lab margin cannot be read as a fuel-level guarantee."]
    (OUT / "README.md").write_text("\n".join(L) + "\n", encoding="utf-8")
    print(f"density {len(dens)} rows / {dens.system_id.nunique()} systems; viscosity {len(visc)} rows / {visc.system_id.nunique()} systems", flush=True)
    print(json.dumps(dst), json.dumps(vst), flush=True)


if __name__ == "__main__":
    main()
