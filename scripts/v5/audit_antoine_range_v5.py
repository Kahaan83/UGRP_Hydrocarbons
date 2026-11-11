"""v5 S2 Question B (R-B1): which pure compounds are used at a predicted flash temperature outside the validity range of the
stored Antoine set?

Records audited (predicted flash T taken from existing v4 outputs, nothing re-fitted):
  * 126 lab mixtures      results/v4/F/flash_lab_loso.csv           (components: comp_smiles)
  * 57 CRC 647 fuels + POSF-10325   results/v4/R3/F/flash_real_fuels.csv (r3 compositions: data/real_fuels/v4/compositions_crc647_r3.csv;
                                    POSF-10325 components from eval_flash_v4.posf_items())
Predicted flash T = ideal_C (the UNIFAC T is also reported; flag is on the ideal prediction, the primary model).

Stored validity range per compound (first rule that applies):
  1. the compound is in data/antoine_sources.csv with a range -> range_T_lo_K .. range_T_hi_K          (basis "antoine_sources")
  2. registry antoine_source == refit_measured -> data range T_min..T_max of data/processed/antoine_refit.csv (basis "refit_data_range")
  3. registry antoine_source == anchored_thermo_ambient -> 273.15..373.15 K, the anchoring window      (basis "anchored_window")
  4. otherwise no range is stored: reported as range_unknown, NOT counted as flagged                     (basis "none")
Reuses scripts/v4/audit_pure_constants.py (canon, _resolve) by import and flash.check_antoine_validity; v4 files are not edited.

Outputs (results/v5/S2/): antoine_range_flags.csv (one row per compound), antoine_range_flags_records.csv (one row per record x compound).
"""
from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from fuelblend.baselines.flash import check_antoine_validity, psat_antoine  # noqa: E402
from fuelblend.baselines.flash_unifac import _default_registry  # noqa: E402

RDLogger.DisableLog("rdApp.*")
OUT = ROOT / "results" / "v5" / "S2"


def _load(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(spec)
    sys.modules[name] = m
    spec.loader.exec_module(m)
    return m


apc = _load(ROOT / "scripts/v4/audit_pure_constants.py", "audit_pure_constants_v4")
ev4 = _load(ROOT / "scripts/v4/eval_flash_v4.py", "eval_flash_v4_mod")


def stored_ranges() -> dict:
    reg = _default_registry()
    src = pd.read_csv(ROOT / "data/antoine_sources.csv")
    src = src.assign(can=src.smiles.map(apc.canon)).set_index("can")
    refit = pd.read_csv(ROOT / "data/processed/antoine_refit.csv").set_index("inchikey")
    out = {}
    for can, r in reg.iterrows():
        lo = hi = np.nan
        basis = "none"
        if can in src.index and pd.notna(src.loc[can].range_T_lo_K):
            s = src.loc[can]
            lo, hi, basis = float(s.range_T_lo_K), float(s.range_T_hi_K), "antoine_sources"
        elif r.antoine_source == "refit_measured" and r.inchikey in refit.index and pd.notna(refit.loc[r.inchikey].T_min):
            lo, hi, basis = float(refit.loc[r.inchikey].T_min), float(refit.loc[r.inchikey].T_max), "refit_data_range"
        elif r.antoine_source == "anchored_thermo_ambient":
            lo, hi, basis = float(r.antoine_Tlo), float(r.antoine_Thi), "anchored_window"
        out[can] = dict(lo=lo, hi=hi, basis=basis, source=r.antoine_source, inchikey=r.inchikey,
                        abc=(r.antoine_A, r.antoine_B, r.antoine_C))
    return out


def collect_records() -> pd.DataFrame:
    rows = []
    lab = pd.read_csv(ROOT / "results/v4/F/flash_lab_loso.csv")
    for i, r in enumerate(lab.itertuples()):
        for s in json.loads(r.comp_smiles):
            rows.append(dict(record=f"lab:{i}", system=r.system_id, set="lab_mixtures", T_ideal_C=r.ideal_C, T_unifac_C=r.unifac_C,
                             smiles=apc.canon(s)))
    rf = pd.read_csv(ROOT / "results/v4/R3/F/flash_real_fuels.csv")
    comp = pd.read_csv(ROOT / "data/real_fuels/v4/compositions_crc647_r3.csv")
    comp = comp[~comp.is_lumped.astype(bool) & (comp.inchikey != "LUMP_NO_REALISTIC_REP")]
    posf = ev4.posf_items()[0]
    for r in rf.itertuples():
        if r.fuel_id == "POSF-10325":
            smis = posf
        else:
            smis = list(comp[comp.system_id == f"fuel:{r.fuel_id}"].smiles)
        for s in sorted(set(apc.canon(x) for x in smis)):
            rows.append(dict(record=f"fuel:{r.fuel_id}", system=r.fuel_id, set="real_fuels_r3", T_ideal_C=r.ideal_C,
                             T_unifac_C=r.unifac_C, smiles=s))
    return pd.DataFrame(rows)


def flag_outside(T_flash_C: float, lo_K: float, hi_K: float):
    """(outside_stored_range, extrapolation_K).  An unknown range (NaN) is never flagged."""
    if np.isnan(lo_K) or np.isnan(hi_K):
        return False, np.nan
    ok, d = check_antoine_validity(T_flash_C, [lo_K], [hi_K])
    return (not ok), d


def thermo_ratio(smiles: str, abc, T_K: float):
    """stored-Antoine Psat / thermo-library Psat at T_K (an independent reference correlation, not a measurement)."""
    try:
        import thermo
        chem = apc._resolve(smiles)
        if chem is None:
            return np.nan
        p = thermo.Chemical(chem.CASs, T=float(T_K)).Psat
        if not p or not np.isfinite(p) or p <= 0:
            return np.nan
        return psat_antoine(T_K, *abc) / (p / 1000.0)
    except Exception:
        return np.nan


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    rng = stored_ranges()
    rec = collect_records()
    info = rec.smiles.map(lambda s: rng.get(s))
    rec["range_lo_K"] = [i["lo"] if i else np.nan for i in info]
    rec["range_hi_K"] = [i["hi"] if i else np.nan for i in info]
    rec["range_basis"] = [i["basis"] if i else "not_in_registry" for i in info]
    flags, dist = [], []
    for r in rec.itertuples():
        f, d = flag_outside(r.T_ideal_C, r.range_lo_K, r.range_hi_K)
        flags.append(f)
        dist.append(d)
    rec["outside_range"] = flags
    rec["extrapolation_K"] = dist
    rec["range_unknown"] = rec.range_lo_K.isna()
    rec.to_csv(OUT / "antoine_range_flags_records.csv", index=False, float_format="%.6g")

    pv = pd.read_csv(ROOT / "data/processed/v4/pure_props_v4.csv")
    pv["can"] = pv.smiles.map(apc.canon)
    pv = pv[pv.prop == "Psat37_8"].set_index("can")
    tm = pd.read_csv(ROOT / "results/v5/thermoml_reextract/thermoml_hc_v5.csv", low_memory=False)
    tm = tm[(tm.n_comp == 1) & tm.prop.str.contains("Vapor", na=False) & (tm.prop_phase == "Liquid")]

    out = []
    for smi, g in rec.groupby("smiles"):
        ri = rng.get(smi)
        flagged = g[g.outside_range]
        Tk = g.T_ideal_C + 273.15
        abc = ri["abc"] if ri else (np.nan,) * 3
        rep_T = float(np.median(flagged.T_ideal_C + 273.15)) if len(flagged) else float(np.median(Tk))
        mol = Chem.MolFromSmiles(smi)
        inchi = Chem.MolToInchi(mol)
        meas = tm[tm.inchis.str.contains(inchi, regex=False, na=False)]
        tq = (flagged.T_ideal_C + 273.15) if len(flagged) else Tk
        near = meas[(meas.T_K >= tq.min() - 5) & (meas.T_K <= tq.max() + 5)] if len(meas) else meas
        v4row = pv.loc[smi] if smi in pv.index else None
        out.append(dict(
            smiles=smi, inchikey=ri["inchikey"] if ri else "", antoine_source=ri["source"] if ri else "not_in_registry",
            range_basis=ri["basis"] if ri else "none", range_lo_K=ri["lo"] if ri else np.nan, range_hi_K=ri["hi"] if ri else np.nan,
            range_unknown=bool(ri is None or np.isnan(ri["lo"])), n_records_using=g.record.nunique(),
            n_records_outside=flagged.record.nunique(), n_systems_outside=flagged.system.nunique(),
            n_lab_records_outside=flagged[flagged.set == "lab_mixtures"].record.nunique(),
            n_fuel_records_outside=flagged[flagged.set == "real_fuels_r3"].record.nunique(),
            flagged=bool(len(flagged) > 0),
            T_flash_min_K=float((tq).min()), T_flash_max_K=float((tq).max()),
            max_extrapolation_K=float(flagged.extrapolation_K.max()) if len(flagged) else 0.0,
            v4_psat37_8_status=v4row.status if v4row is not None else "no_v4_row",
            stored_over_thermo_ref_at_median_flash_T=thermo_ratio(smi, abc, rep_T) if ri else np.nan,
            median_flash_T_K=rep_T,
            n_thermoml_liquid_psat_points=len(meas),
            thermoml_T_range_K=(f"{meas.T_K.min():.1f}-{meas.T_K.max():.1f}" if len(meas) else ""),
            n_thermoml_points_within_5K_of_flagged_T=len(near),
            independent_psat_check_at_flash_T=bool(len(near) > 0),
            thermo_library_reference_ratio_available=bool(np.isfinite(thermo_ratio(smi, abc, rep_T))) if ri else False))
    res = pd.DataFrame(out).sort_values(["flagged", "n_records_outside"], ascending=False)
    res.to_csv(OUT / "antoine_range_flags.csv", index=False, float_format="%.6g")

    fl = res[res.flagged]
    print(f"records audited: {rec.record.nunique()} ({rec[rec.set=='lab_mixtures'].record.nunique()} lab + "
          f"{rec[rec.set=='real_fuels_r3'].record.nunique()} fuels); compounds: {len(res)}")
    print(res.groupby("range_basis").agg(n=("smiles", "size"), flagged=("flagged", "sum")))
    print(f"FLAGGED compounds: {len(fl)}; range_unknown compounds: {int(res.range_unknown.sum())}")
    print("records with >=1 flagged compound:", rec[rec.outside_range].record.nunique(),
          "systems:", rec[rec.outside_range].system.nunique())
    pd.set_option("display.width", 250)
    pd.set_option("display.max_columns", 40)
    print(fl[["smiles", "antoine_source", "range_basis", "range_lo_K", "range_hi_K", "n_records_outside", "n_systems_outside",
              "max_extrapolation_K", "stored_over_thermo_ref_at_median_flash_T", "independent_psat_check_at_flash_T"]].to_string())


if __name__ == "__main__":
    main()
