"""NREL FuelLib loader (v5): GCxGC bin compositions and digitised Edwards 2017/2020 properties.

Reads the untrusted-data copy in ``data/raw/fuellib/`` (see its PROVENANCE.md; BSD-3) and writes

* ``data/real_fuels/v5/fuellib_gcxgc_compositions.csv``  one row per (fuel, GC bin) with mass/mole fraction,
  SMILES/InChIKey of the bin's *reference compound*, class and carbon number.  The first nine columns match
  ``data/real_fuels/v4/compositions_*.csv``; ``is_lumped`` is True for bins (a bin stands for all isomers of
  its class x carbon number, represented by ONE reference molecule) and False for the single-compound fuels.
* ``data/real_fuels/v5/fuellib_properties.csv``  tier-B (secondary, digitised/transcribed by the FuelLib authors)
  properties, long format, columns of ``data/real_fuels/v4/jet_properties_v4.csv`` plus extras.
* ``data/real_fuels/v5/fuellib_build_log.json``

Double counting: the POSF/NJFCP property values are the same measurements as Edwards 2020 (AFRL-TR-2020-0017,
AD1093317) and Edwards 2017 (AIAA 2017-0146).  Every such row carries ``double_count_group =
'AFRL-TR-2020-0017'``; a consumer that also loads AFRL-TR-2020-0017 / CRC AV-27-18 / v4 POSF-10325 rows must
drop one copy (prefer the primary source).  Pure-compound property files (NIST WebBook copies) are NOT loaded.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import pandas as pd
import yaml
from rdkit import Chem
from rdkit.Chem import Descriptors

REPO = Path(__file__).resolve().parents[3]
RAW = REPO / "data" / "raw" / "fuellib"
OUT = REPO / "data" / "real_fuels" / "v5"
PROVENANCE = "fuellib_gcxgc_tierB"
PSI_TO_KPA = 6.894757293168

COMP_COLS = ["system_id", "inchikey", "smiles", "fraction", "basis", "is_lumped", "unexplained_fraction",
             "provenance", "verified", "fuel_id", "fuel_name", "category", "gc_bin", "reference_compound",
             "hc_class", "n_carbon", "mw", "mass_fraction_norm", "mole_fraction", "weight_pct_raw",
             "sum_weight_pct", "source_ref", "v4_overlap"]
PROP_COLS = ["fuel_id", "property", "value", "unit", "T_C", "method", "role", "source_file", "source_page",
             "provenance", "raw_value", "raw_unit", "fuel_name", "data_tier", "double_count_group",
             "digitisation", "T_note", "hefa_wt_pct", "v4_overlap"]
TARGETS = {"density_kg_m3", "viscosity_kin_cSt", "flash_point_C", "freezing_point_C",
           "net_heat_of_combustion_MJ_kg"}

# fuels whose property files are Edwards 2017/2020 (or NJFCP SwRI) data; others are NIST WebBook copies
PROPERTY_FUELS = ["posf10264", "posf10289", "posf10325", "posf11498"]
V4_OVERLAP = {"posf10325": "POSF-10325"}   # v4 loaders_jet: GC fuel POSF-10325 (= NJFCP A-2 = CRC AV-27-18 code 10325)
NJFCP_GROUP = {"posf10264", "posf10289", "posf10325", "posf11498"}


def _hc_class(gc_bin: str) -> str:
    b = gc_bin.lower()
    if "isoparaffin" in b:
        return "iso_paraffin"
    if b.startswith("n-c"):
        return "n_paraffin"
    if "cycloparaffin" in b:
        return "cycloalkane"
    if "alkene" in b:
        return "alkene_alkyne"
    if "benzene" in b or "aromatic" in b or b == "toluene":
        return "aromatic"
    return "other"


def _read_gc(path: Path) -> pd.DataFrame:
    d = pd.read_csv(path)
    d = d.rename(columns={"Weight %": "wt_pct", "GC-Bin": "gc_bin", "SMILES": "smiles",
                          "Reference Compound": "reference_compound"})
    if "reference_compound" not in d:
        d["reference_compound"] = ""
    d["wt_pct"] = pd.to_numeric(d["wt_pct"], errors="raise")
    return d[["gc_bin", "reference_compound", "smiles", "wt_pct"]]


def load_metadata(raw: Path = RAW) -> dict:
    return yaml.safe_load((raw / "fuel_metadata.yaml").read_text(encoding="utf-8"))["fuels"]


def build_compositions(raw: Path = RAW, log: dict | None = None) -> pd.DataFrame:
    log = {} if log is None else log
    meta = load_metadata(raw)
    rows = []
    log.setdefault("gc_fuels", {})
    for f in sorted((raw / "gcData").glob("*.csv")):
        slug = f.stem
        d = _read_gc(f)
        if d["gc_bin"].duplicated().any():
            raise ValueError(f"{f.name}: duplicated GC bins")
        m = meta.get(slug, {})
        total = float(d["wt_pct"].sum())
        nz = d[d["wt_pct"] > 0].copy()
        mols = [Chem.MolFromSmiles(s) for s in nz["smiles"]]
        if any(x is None for x in mols):
            raise ValueError(f"{f.name}: unparsable SMILES")
        nz["mw"] = [Descriptors.MolWt(x) for x in mols]
        nz["inchikey"] = [Chem.MolToInchiKey(x) for x in mols]
        nz["n_carbon"] = [sum(a.GetSymbol() == "C" for a in x.GetAtoms()) for x in mols]
        nz["hydrocarbon"] = [all(a.GetSymbol() in ("C", "H") for a in x.GetAtoms()) for x in mols]
        if not nz["hydrocarbon"].all():
            raise ValueError(f"{f.name}: non-hydrocarbon reference compound")
        wn = nz["wt_pct"] / nz["wt_pct"].sum()
        mol_ = (wn / nz["mw"]) / (wn / nz["mw"]).sum()
        simple = len(d) == 1 or m.get("category") == "Simple"
        for (_, r), w, x in zip(nz.iterrows(), wn, mol_):
            rows.append({
                "system_id": f"fuel:FUELLIB_{slug}", "inchikey": r["inchikey"], "smiles": r["smiles"],
                "fraction": r["wt_pct"] / 100.0, "basis": "mass", "is_lumped": not simple,
                "unexplained_fraction": max(0.0, 1.0 - total / 100.0), "provenance": PROVENANCE, "verified": False,
                "fuel_id": f"FUELLIB_{slug}", "fuel_name": m.get("name", slug), "category": m.get("category", ""),
                "gc_bin": r["gc_bin"], "reference_compound": r["reference_compound"],
                "hc_class": _hc_class(r["gc_bin"]), "n_carbon": int(r["n_carbon"]), "mw": float(r["mw"]),
                "mass_fraction_norm": float(w), "mole_fraction": float(x), "weight_pct_raw": float(r["wt_pct"]),
                "sum_weight_pct": total, "source_ref": f"{m.get('source', '')}; {m.get('reference', '')}",
                "v4_overlap": V4_OVERLAP.get(slug, ""),
            })
        log["gc_fuels"][slug] = {"n_bins_file": int(len(d)), "n_bins_nonzero": int(len(nz)),
                                 "sum_weight_pct": round(total, 3)}
    return pd.DataFrame(rows, columns=COMP_COLS)


# ---------------------------------------------------------------- properties
def _std(prop: str, val: float, unit: str):
    """FuelLib property row -> (canonical name, value, unit) or None."""
    if prop == "Density":
        return "density_kg_m3", val * 1000.0 if unit == "g/cm^3" else math.nan, "kg/m3"
    if prop == "Viscosity":
        return "viscosity_kin_cSt", val if unit == "mm^2/s" else math.nan, "cSt"
    if prop == "FlashPoint":
        return "flash_point_C", val, "C"
    if prop == "FreezePoint":
        return "freezing_point_C", val, "C"
    if prop == "DerivedCetaneNumber":
        return "dcn", val, "-"
    if prop == "HeatOfCombustion":
        return "net_heat_of_combustion_MJ_kg", val if unit == "MJ/kg" else math.nan, "MJ/kg"
    if prop == "SmokePoint":
        return "smoke_point_mm", val, "mm"
    if prop == "HeatCapacity":
        return "heat_capacity_kJ_kg_K", val if unit == "kJ/kg/K" else math.nan, "kJ/kg/K"
    if prop == "ThermalConductivity":
        return "thermal_conductivity_W_m_K", val if unit == "W/m/K" else math.nan, "W/m/K"
    if prop == "SurfaceTension":
        return "surface_tension_N_m", (val if unit == "N/m" else val * 1e-3 if unit == "dynes/cm" else math.nan), "N/m"
    if prop == "VaporPressure":
        return "vapor_pressure_kPa", (val if unit == "kPa" else val * PSI_TO_KPA if unit == "psi" else math.nan), "kPa"
    return None


def _digitisation(source: str) -> str:
    s = source or ""
    return "figure_digitised" if "Fig" in s else "table_transcribed"


def build_properties(raw: Path = RAW, log: dict | None = None, fuels=PROPERTY_FUELS) -> pd.DataFrame:
    log = {} if log is None else log
    meta = load_metadata(raw)
    out, dropped = [], 0
    for slug in fuels:
        p = raw / "propertiesData" / f"{slug}.csv"
        d = pd.read_csv(p, header=0, names=["Temp", "Temp_Units", "Property", "Property_Units",
                                            "Property_Value", "Source"])
        for _, r in d.iterrows():
            if str(r["Temp_Units"]).strip().lower() != "celsius":
                raise ValueError(f"{p.name}: unexpected temperature unit {r['Temp_Units']}")
            prop, unit = r["Property"], r["Property_Units"]
            pv = float(r["Property_Value"])
            T = r["Temp"]
            if prop == "Distillation":
                # Temp = boiling temperature (C), Property_Value = volume % recovered
                name, val, u, TC = f"distillation_T{int(round(pv))}_C", float(T), "C", math.nan
                raw_val, raw_unit = float(T), "celsius"
            else:
                std = _std(prop, pv, unit)
                if std is None or std[1] is None or (isinstance(std[1], float) and math.isnan(std[1])):
                    dropped += 1
                    continue
                name, val, u = std
                TC = float(T) if pd.notna(T) else math.nan
                raw_val, raw_unit = pv, unit
            src = str(r["Source"]) if pd.notna(r["Source"]) else "NJFCP (SwRI)"
            out.append({
                "fuel_id": f"FUELLIB_{slug}", "property": name, "value": val, "unit": u, "T_C": TC,
                "method": src, "role": "target" if name in TARGETS else "descriptor",
                "source_file": f"data/raw/fuellib/propertiesData/{slug}.csv", "source_page": src,
                "provenance": PROVENANCE, "raw_value": raw_val, "raw_unit": raw_unit,
                "fuel_name": meta.get(slug, {}).get("name", slug), "data_tier": "B",
                "double_count_group": "AFRL-TR-2020-0017" if slug in NJFCP_GROUP else "",
                "digitisation": _digitisation(src), "T_note": "", "hefa_wt_pct": math.nan,
                "v4_overlap": V4_OVERLAP.get(slug, ""),
            })
    # HEFA / Jet-A blends (mass %), Vozka et al. 2018 via FuelLib; temperature is NOT stated in the file
    bp = raw / "propertiesData" / "hefa-jet-a-blends.csv"
    b = pd.read_csv(bp, header=0, skiprows=[1]) if bp.exists() else pd.DataFrame()
    for _, r in b.iterrows():
        w = float(r["HEFA_concentration"])
        for feed in ("TALL", "CAME", "MFAT"):
            for pcol, name, mult, u, ru in (("Viscosity", "viscosity_kin_cSt", 1.0, "cSt", "mm^2/s"),
                                            ("Density", "density_kg_m3", 1000.0, "kg/m3", "g/cm^3")):
                v = float(r[f"{pcol}_{feed}"])
                out.append({
                    "fuel_id": f"FUELLIB_jet-a__hefa-{feed.lower()}__w{int(w):03d}", "property": name,
                    "value": v * mult, "unit": u, "T_C": math.nan,
                    "method": "FuelLib blend table (Vozka 2018 attribution)", "role": "target",
                    "source_file": "data/raw/fuellib/propertiesData/hefa-jet-a-blends.csv", "source_page": "",
                    "provenance": PROVENANCE, "raw_value": v, "raw_unit": ru,
                    "fuel_name": f"Jet-A + {w:g} wt% HEFA-{feed}", "data_tier": "B", "double_count_group": "",
                    "digitisation": "table_transcribed", "T_note": "temperature not stated in FuelLib file",
                    "hefa_wt_pct": w, "v4_overlap": "",
                })
    log["properties_dropped_unconvertible"] = dropped
    return pd.DataFrame(out, columns=PROP_COLS)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    log: dict = {}
    comp = build_compositions(RAW, log)
    prop = build_properties(RAW, log)
    comp.to_csv(OUT / "fuellib_gcxgc_compositions.csv", index=False)
    prop.to_csv(OUT / "fuellib_properties.csv", index=False)
    log["n_composition_rows"] = int(len(comp))
    log["n_property_rows"] = int(len(prop))
    (OUT / "fuellib_build_log.json").write_text(json.dumps(log, indent=1), encoding="utf-8")
    print(f"compositions: {len(comp)} rows, {comp['fuel_id'].nunique()} fuels; properties: {len(prop)} rows")


if __name__ == "__main__":
    main()
