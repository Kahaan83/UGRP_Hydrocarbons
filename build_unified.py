#!/usr/bin/env python
"""Pull every source into ONE canonical schema.

    python build_unified.py                 # flash point, viscosity, density, vapour pressure
    python build_unified.py --octane        # also DCN / MON (not your target properties)
    python build_unified.py --all-molecules # keep non-hydrocarbon rows too (default: hydrocarbon only)

Inputs (any that are missing are skipped, and the report says so):
  data/interim/**/flashpoint_{pure,mix}_exp.csv, viscosity_{pure,mix}_exp.csv   (Leenhouts benchmark)
  data/processed/viscosity_{pure,mix}_hc.csv                                    (NIST/TDE via build_viscosity.py)
  data/processed/thermoml_hc.csv                                                (thermoml_extract.py)
  data/processed/ajftd_fuels.csv                                                (ingest_ajftd.py)
  data/processed/antoine_refit.csv | antoine_best.csv                           (refit_antoine.py | webbook_antoine.py)
  data/raw/lfl/lfl.csv                                                          (you; template is written for you)

Outputs in data/unified/:
  compounds.csv          registry: InChIKey, SMILES, formula, MW, class, + Antoine / LFL columns
  measurements.csv       every observation, one canonical row each (before dedupe)
  measurements_dedup.csv one row per (property, system, composition, T): median value, sources, conflict flag
  measurements_out_of_domain.csv  rows removed by the fuel-domain filter (see DOMAIN; disable with --no-domain-filter)
  UNIFIED_REPORT.md      counts, ranges, conflicts, and what physics inputs are still missing

Canonical properties / units
  flash_point_C | viscosity_dyn_mPas | viscosity_kin_cSt | density_kg_m3 | vapor_pressure_kPa |
  freezing_point_C | pour_point_C | cloud_point_C | dcn | mon
Composition: x_mole and w_mass (JSON lists aligned to comp_inchikeys, sorted by InChIKey).
`system_id` = hash of the component SET: use it as the group key for GroupKFold so that a
blend family never appears in both train and test.
"""
from __future__ import annotations
import argparse, hashlib, json, math, re
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

from common import carbon_number, hc_class, is_hydrocarbon, mol_from, mol_weight

INT, PROC, RAW, OUT = Path("data/interim"), Path("data/processed"), Path("data/raw"), Path("data/unified")
WINDOW = (313.15, 373.15)
notes: list[str] = []        # things that went wrong / were skipped -> go into the report


# ------------------------------------------------------------------ compound registry
class Registry:
    def __init__(self):
        self.by_key: dict[str, dict] = {}
        self._cache: dict[tuple, str | None] = {}
        self.unparsed = 0

    def _add(self, m) -> str | None:
        if m is None:
            self.unparsed += 1
            return None
        key = Chem.MolToInchiKey(m)
        if key not in self.by_key:
            hc = is_hydrocarbon(m)
            self.by_key[key] = dict(
                inchikey=key, smiles=Chem.MolToSmiles(m), inchi=Chem.MolToInchi(m),
                formula=rdMolDescriptors.CalcMolFormula(m), mw=mol_weight(m), is_hc=hc,
                cls=hc_class(m) if hc else "other", n_carbon=carbon_number(m))
        return key

    def inchi(self, s: str) -> str | None:
        k = ("i", s)
        if k not in self._cache:
            self._cache[k] = self._add(mol_from(inchi=s))
        return self._cache[k]

    def smiles(self, s: str) -> str | None:
        k = ("s", s)
        if k not in self._cache:
            self._cache[k] = self._add(mol_from(smiles=s))
        return self._cache[k]


REG = Registry()


# ------------------------------------------------------------------ row construction
def make_row(source, prop, value, keys, x_mole, T_K=np.nan, P_kPa=np.nan, ref="", raw_value=np.nan,
             raw_unit="", fraction_basis="mole", fuel_id="", composition_level="L3", comp_coverage=1.0):
    """keys: list of InChIKeys; x_mole: aligned fractions (mole or mass, see fraction_basis)."""
    if value is None or not np.isfinite(value):
        return None
    if keys:
        x = np.asarray(x_mole, float)
        if len(keys) != len(x) or not np.all(np.isfinite(x)) or x.sum() <= 0:
            return None
        cov = float(x.sum())
        if abs(cov - 1) > 1e-2:
            notes.append(f"{source}: fractions sum to {cov:.3f} (renormalised)")
        x = x / cov
        M = np.array([REG.by_key[k]["mw"] for k in keys])
        if fraction_basis == "mass":
            w = x
            x = (w / M) / np.sum(w / M)
        else:
            w = x * M / np.sum(x * M)
        order = np.argsort(keys)
        keys = [keys[i] for i in order]; x = x[order]; w = w[order]
        is_hc = all(REG.by_key[k]["is_hc"] for k in keys)
        system = hashlib.md5("|".join(keys).encode()).hexdigest()[:12]
        n = len(keys)
        kind = "pure" if n == 1 else "binary" if n == 2 else "multi"
        clevel = "L3" if kind != "pure" else "L3"
        comp_cov = comp_coverage if comp_coverage is not None else cov
        cbasis = fraction_basis
    else:
        x = w = np.array([]); is_hc = True; system = "fuel:" + str(fuel_id); n = 0; kind = "real_fuel"
        clevel = composition_level if composition_level in ("L0", "L1", "L2") else "L0"
        comp_cov = comp_coverage if comp_coverage is not None else 0.0
        cbasis = fraction_basis if fraction_basis else "none"

    return dict(source=source, property=prop, value=float(value), T_K=T_K, P_kPa=P_kPa, kind=kind, n_comp=n,
                is_hc=is_hc, system_id=system, comp_inchikeys=json.dumps(keys),
                x_mole=json.dumps([round(float(v), 6) for v in x]), w_mass=json.dumps([round(float(v), 6) for v in w]),
                fuel_id=fuel_id, ref=str(ref), raw_value=raw_value, raw_unit=raw_unit,
                composition_level=clevel, comp_coverage=float(comp_cov), comp_basis=cbasis)



def find(root: Path, name: str) -> Path | None:
    if not root.exists():
        return None
    for p in sorted(root.rglob(name)):
        if p.is_file() and "__MACOSX" not in p.parts and not p.name.startswith(("._", "~$")):
            return p
    return None


def to_T_K(series: pd.Series, src: str) -> pd.Series:
    t = pd.to_numeric(series, errors="coerce")
    if t.notna().any() and t.dropna().median() < 150:     # looks like deg C
        notes.append(f"{src}: temperatures look like Celsius -> converted to K")
        t = t + 273.15
    return t


def std_leenhouts(prop_name: str, value: float, unit: str):
    u = (unit or "").lower()
    if prop_name == "flash_point_C":
        if u.startswith("k"):
            return value - 273.15
        return value                                       # 'Celcius' (sic)
    if prop_name == "viscosity_dyn_mPas":
        if "ln" in u:
            return math.exp(value)
        if "log" in u:
            return 10 ** value
        return value
    return value


# ------------------------------------------------------------------ loaders
def load_leenhouts(path: Path, prop: str, source="leenhouts") -> list[dict]:
    df = pd.read_csv(path, low_memory=False)
    idx = sorted(int(m.group(1)) for c in df.columns if (m := re.fullmatch(r"component_inchi_(\d+)", c)))
    tcol = next((c for c in df.columns if re.match(r"(?i)temp", c)), None)
    T = to_T_K(df[tcol], path.name) if tcol else pd.Series(np.nan, index=df.index)
    ftype = df["fraction_type"].iloc[0] if "fraction_type" in df else "molar"
    basis = "mass" if str(ftype).lower().startswith("mass") else "mole"
    rows = []
    for i, r in enumerate(df.to_dict("records")):
        keys, xs = [], []
        for j in idx:
            inchi = r.get(f"component_inchi_{j}")
            if isinstance(inchi, str) and inchi.startswith("InChI"):
                keys.append(REG.inchi(inchi))
                xs.append(r.get(f"fraction_{j}", np.nan))
        if not keys or any(k is None for k in keys):
            continue
        if len(keys) == 1:
            xs = [1.0]
        elif any(pd.isna(x) for x in xs):
            continue
        v = std_leenhouts(prop, float(r["value"]), str(r.get("unit", "")))
        row = make_row(source, prop, v, keys, xs, T_K=T.iloc[i], raw_value=r["value"], raw_unit=str(r.get("unit", "")),
                       ref=path.name, fraction_basis=basis)
        if row:
            rows.append(row)
    return rows


def load_nist_hc() -> list[dict]:
    rows = []
    p = PROC / "viscosity_pure_hc.csv"
    if p.exists():
        for r in pd.read_csv(p).to_dict("records"):
            k = REG.smiles(r["smiles"])
            if k:
                row = make_row("nist_tde", "viscosity_dyn_mPas", r["eta_mPas"], [k], [1.0], T_K=r["T_K"], P_kPa=r.get("P_kPa", np.nan),
                               ref=r.get("Ref_ID", ""), raw_value=r["ln_eta_Pas"], raw_unit="ln(Pa*s)")
                if row:
                    rows.append(row)
    p = PROC / "viscosity_mix_hc.csv"
    if p.exists():
        for r in pd.read_csv(p).to_dict("records"):
            keys = [REG.smiles(s) for s in json.loads(r["smiles_list"])]
            if any(k is None for k in keys):
                continue
            row = make_row("nist_tde", "viscosity_dyn_mPas", r["eta_mPas"], keys, json.loads(r["x_mole"]), T_K=r["T_K"],
                           P_kPa=r.get("P_kPa", np.nan), ref=r.get("Ref_ID", ""), raw_value=r["ln_eta_Pas"], raw_unit="ln(Pa*s)")
            if row:
                rows.append(row)
    return rows


THERMOML_PROPS = {            # exact ThermoML names -> (canonical property, factor to canonical unit)
    "Mass density, kg/m3": ("density_kg_m3", 1.0),
    "Viscosity, Pa*s": ("viscosity_dyn_mPas", 1e3),
    "Kinematic viscosity, m2/s": ("viscosity_kin_cSt", 1e6),
    "Vapor or sublimation pressure, kPa": ("vapor_pressure_kPa", 1.0),
}


def thermoml_prop(text: str):
    return THERMOML_PROPS.get(str(text).strip(), (None, None))


def load_thermoml() -> list[dict]:
    p = PROC / "thermoml_hc.csv"
    if not p.exists():
        return []
    df = pd.read_csv(p)
    if "prop_phase" not in df.columns:
        notes.append("thermoml_hc.csv is from the OLD extractor (no prop_phase column; it mixed in excess viscosities, gas-phase "
                     "densities, high pressures and VLE vapour compositions) -> IGNORED. Re-run: python thermoml_extract.py")
        return []
    rows, skipped = [], 0
    for r in df.to_dict("records"):
        prop, f = thermoml_prop(str(r["prop"]))
        if not prop:
            continue
        inchis = json.loads(r["inchis"]); vals = json.loads(r["comp_values"])
        keys = [REG.inchi(i) if i else None for i in inchis]
        if any(k is None for k in keys):
            skipped += 1; continue
        if len(keys) == 1:
            xs, basis = [1.0], "mole"
        else:
            ct = str(r.get("comp_type", "")).lower()
            if None in vals or ("mole" not in ct and "mass" not in ct):
                skipped += 1; continue
            xs, basis = vals, ("mass" if "mass" in ct else "mole")
        row = make_row("thermoml", prop, r["value"] * f, keys, xs, T_K=r.get("T_K", np.nan), P_kPa=r.get("P_kPa", np.nan),
                       ref=r.get("file", ""), raw_value=r["value"], raw_unit=str(r["prop"]), fraction_basis=basis)
        if row:
            rows.append(row)
    if skipped:
        notes.append(f"thermoml: {skipped} rows skipped (no InChI or no usable composition)")
    return rows


def load_ajftd() -> list[dict]:
    p = PROC / "ajftd_fuels.csv"
    if not p.exists():
        return []
    rows = []
    for r in pd.read_csv(p).to_dict("records"):
        prop = str(r["property"])
        # Bulk properties / distillation without GCxGC are L0; class bulk totals are L1
        clevel = "L1" if ("aromatics" in prop or "saturates" in prop or "olefins" in prop) else "L0"
        row = make_row("ajftd", prop, r["value"], [], [], T_K=r.get("T_K", np.nan), fuel_id=str(r["fuel_id"]),
                       ref=str(r.get("method", "")), raw_value=r["value"], raw_unit=str(r.get("unit_std", "")),
                       composition_level=clevel, comp_coverage=0.0, fraction_basis="none")
        if row:
            rows.append(row)
    return rows



def load_flashpoint_li() -> list[dict]:
    """Sun/Krakauer/Morgan (Mol Inf 2020, arXiv 2002.11315) 10k pure-compound flash points.

    Columns include SMILES and flash point in K. Duplicates collapse in dedupe().
    """
    p = Path("data/raw/flashpoint_li/integrated_dataset_public.csv")
    if not p.exists():
        return []
    df = pd.read_csv(p)
    rows = []
    for r in df[["smiles", "flashpoint", "pure substance"]].dropna().to_dict("records"):
        if int(r["pure substance"]) != 1:
            continue
        key = REG.smiles(str(r["smiles"]))
        if not key:
            continue
        T_C = float(r["flashpoint"]) - 273.15           # dataset flash points are in K
        row = make_row("flashpoint_li", "flash_point_C", T_C, [key, ],
                       [1.0], ref=str(r.get("source", "MolInf2020_figshare9275210")),
                       raw_value=float(r["flashpoint"]), raw_unit="K")
        if row:
            rows.append(row)
    return rows


def load_fortin() -> list[dict]:
    """Fortin et al. 2012 RP-1/RP-2 orthogonality (NIST): real-fuel distillation curves."""
    p = Path("data/raw/rp1/fortin2012_tables.csv")
    if not p.exists():
        return []
    rows = []
    for r in pd.read_csv(p).to_dict("records"):
        # Distillation curves without GCxGC are L0
        row = make_row("fortin2012", str(r["property"]), float(r["value_std"]), [], [],
                       fuel_id=str(r["fuel_id"]), ref=str(r.get("source_pdf", "")),
                       raw_value=float(r["value_std"]), raw_unit=str(r.get("unit_std", "")),
                       composition_level="L0", comp_coverage=0.0, fraction_basis="none")
        if row:
            rows.append(row)
    return rows


def load_crc() -> list[dict]:
    """CRC 647 / AV-27-18 / AV-11-09 aviation-fuel tables (data/raw/pour_point/*.csv).

    One row per (fuel, property); property names are canonicalised by parse_crc_tables.py.
    """
    rows = []
    for fn, src in (("crc647_fuel_properties.csv", "crc647"),
                    ("av2718_fuel_properties.csv", "av2718"),
                    ("av1109_freezing_point.csv", "av1109")):
        p = Path("data/raw/pour_point") / fn
        if not p.exists():
            continue
        try:
            df = pd.read_csv(p)
        except Exception as exc:
            notes.append(f"{fn}: unreadable ({exc})")
            continue
        for r in df.to_dict("records"):
            try:
                v = float(r.get("value_std"))
            except (TypeError, ValueError):
                continue
            if not np.isfinite(v):
                continue
            prop = str(r["property"])
            T_C = r.get("T_C")
            # Density at 60°F (15.56 °C) is standard for aviation fuels; infer when blank.
            try:
                T_val = float(T_C) if T_C not in (None, "", np.nan) else None
            except (TypeError, ValueError):
                T_val = None
            if prop == "density_kg_m3" and T_val is None:
                T_K = 288.71  # 60 °F = 15.555... °C + 273.15 ≈ 288.7056, rounded 288.71
            elif T_val is not None:
                T_K = T_val + 273.15
            else:
                T_K = np.nan
            clevel = "L1" if ("aromatics" in prop or "saturates" in prop or "olefins" in prop or "piona" in prop) else "L0"
            row = make_row(src, prop, v, [], [], T_K=T_K,
                           fuel_id=str(r.get("fuel_code", "")), ref=str(r.get("method", "")),
                           raw_value=v, raw_unit=str(r.get("unit_std", "")),
                           composition_level=clevel, comp_coverage=0.0, fraction_basis="none")
            if row:
                rows.append(row)
    return rows



def load_chew2024() -> list[dict]:
    p = Path("data/raw/chew2024/chew2024_viscosity.csv")
    if not p.exists():
        return []
    rows = []
    for r in pd.read_csv(p).to_dict("records"):
        key = REG.smiles(str(r.get("smiles", "")))
        if not key:
            key = REG.smiles(str(r.get("compound_name", "")))
        if not key:
            continue
        try:
            v = float(r["viscosity_mPas"])
            T = float(r["T_K"]) if pd.notna(r.get("T_K")) else np.nan
        except Exception:
            continue
        if v <= 0 or not np.isfinite(v):
            continue
        row = make_row("chew2024", "viscosity_dyn_mPas", v, [key], [1.0],
                       T_K=T, ref=str(r.get("source", "SI_S4.1")),
                       raw_value=v, raw_unit="mPa.s")
        if row:
            rows.append(row)
    return rows


def load_saldana11() -> list[dict]:
    p = Path("data/raw/flashpoint_li/saldana11_flashpoints.csv")
    if not p.exists():
        return []
    rows = []
    for r in pd.read_csv(p).to_dict("records"):
        key = REG.smiles(str(r.get("smiles", "")))
        if not key:
            continue
        try:
            fp = float(r["flashpoint_value"])
            unit = str(r.get("flashpoint_unit", "K")).strip().lower()
        except Exception:
            continue
        T_C = fp - 273.15 if ("k" in unit or unit == "kelvin") else float(fp)
        row = make_row("saldana11", "flash_point_C", T_C, [key], [1.0],
                       ref="Saldana2011_EnergyFuels_25_3900",
                       raw_value=fp, raw_unit=str(r.get("flashpoint_unit", "K")))
        if row:
            rows.append(row)
    return rows


def load_pure_freezing() -> list[dict]:
    p = Path("data/raw/pour_point/pure_compound_freezing_points.csv")
    if not p.exists():
        return []
    rows = []
    for r in pd.read_csv(p).to_dict("records"):
        key = REG.smiles(str(r.get("smiles", "")))
        if not key:
            continue
        try:
            v = float(r["value_C"])
        except Exception:
            continue
        row = make_row("pure_freezing_nist", str(r.get("property", "freezing_point_C")), v,
                       [key], [1.0], ref=str(r.get("source", "NIST_WebBook")),
                       raw_value=v, raw_unit="°C")
        if row:
            rows.append(row)
    return rows


def load_surrogates() -> list[dict]:
    """Jet-fuel surrogates (multi-component recipes from literature). Each fuel → blend records."""
    p = Path("data/raw/surrogates/fuel_surrogates.csv")
    if not p.exists():
        return []
    df = pd.read_csv(p)
    groups = df.groupby(["fuel_code", "fuel_name", "reference"])
    rows = []
    for (fc, fn, ref), g in groups:
        keys = []
        vals = []
        for _, r in g.iterrows():
            s = str(r.get("component_smiles", "")).strip()
            if not s:
                continue  # skip unassigned iso-class designations
            k = REG.smiles(s)
            if not k:
                continue
            keys.append(k)
            try:
                frac = float(r["frac_value"])
            except Exception:
                continue
            if frac <= 0:
                continue
            vals.append(frac)
        if not keys or not vals or len(keys) < 2:
            continue  # need multi-component to be useful as blend
        x = np.array(vals, float)
        # surrogate fractions are in mole (most papers) — if mass, note; default mole per most sources
        ft = str(g.iloc[0].get("frac_type", "mole")).lower()
        if ft.startswith("mass") or ft.startswith("wt"):
            basis = "mass"
        else:
            basis = "mole"
        try:
            row = make_row("surrogate", "compose", 0.0, keys, list(x / x.sum()),
                           fuel_id=str(fc), ref=str(ref), fraction_basis=basis)
        except Exception:
            continue
        # surrogates have composition, not a measured property — skip property-value row,
        # but add a marker record for tracking (lightweight; use source only for tracking)
        # Instead: produce a synthetic blend record as a pure-composition multi entry with a marker property
        # For simplicity keep just the composition registration via a dummy property; but build_unified
        # requires a property/value — use a placeholder and document
        row_meas = make_row("surrogate", "blend_composition", 1.0, keys, list(x / x.sum()),
                            T_K=np.nan, fuel_id=str(fc), ref=str(ref),
                            raw_value=1.0, raw_unit="blend_unit",
                            fraction_basis=basis)
        if row_meas:
            # tag as multi blend with correct keys; the 'value' is meaningless — treat as metadata
            # Instead of polluting measurements, skip property-value load; only register if needed.
            # Better: skip make_row with dummy value — just don't load surrogates into measurements
            # unless they have a measured property (they don't). Use row for nothing.
            pass
    # Decision: surrogates don't belong in measurements.csv (no measured property; they're recipes).
    # Save for later use; do NOT append to unified measurements (would create fake values).
    # Just write a note and return empty (no false data).
    notes.append("surrogate recipes saved to data/raw/surrogates/ — not loaded into measurements (recipes, not measurements)")
    return []


def load_gc_composition() -> list[dict]:
    """Real-fuel GC class distributions — real_fuel kind, no molecule lists for totals; named compounds as multi if fractions given."""
    p = Path("data/raw/gc_composition/fuel_gc_compositions.csv")
    if not p.exists():
        return []
    df = pd.read_csv(p)
    rows = []
    # Only named-component rows with SMILES-like fractions can act as blends; most are class totals.
    # For simplicity, treat as real_fuel tracking (fuel_id, property = class + pct) — not measurement
    # Return empty to avoid false measurements; note written.
    notes.append("gc_composition class tables saved to data/raw/gc_composition/ — not loaded as measurements (summary statistics, not lab observations)")
    return []


# ------------------------------------------------------------------ application-domain filter
DOMAIN = {   # property -> (T window K, value window) ; fuel-relevant liquids only
    "density_kg_m3": ((230.0, 470.0), (400.0, 1300.0)),
    "viscosity_dyn_mPas": ((230.0, 470.0), (0.05, 1e5)),
    "viscosity_kin_cSt": ((230.0, 470.0), (0.05, 1e4)),
    "vapor_pressure_kPa": ((230.0, 520.0), (1e-6, 1000.0)),
    # flash: fuel-relevant window; C25+ wax/asphalt flash points (250-450 C) are beyond
    # Antoine validity and outside the fuel domain
    "flash_point_C": ((300.0, 700.0), (-90.0, 250.0)),
}


def apply_domain(meas: pd.DataFrame):
    """Drop observations outside the fuel domain (near-critical/gas-like 'liquid' densities, negative or absurd
    viscosities, MPa-GPa 'vapour pressures' of light-gas VLE ...). Returns (kept, dropped) - dropped is saved for inspection."""
    bad = pd.Series(False, index=meas.index)
    reason = pd.Series("", index=meas.index, dtype=object)
    for prop, ((t0, t1), (v0, v1)) in DOMAIN.items():
        m = (meas["property"] == prop) & (meas["kind"] != "real_fuel")
        bt = m & meas["T_K"].notna() & ((meas["T_K"] < t0) | (meas["T_K"] > t1))
        bv = m & ((meas["value"] < v0) | (meas["value"] > v1))
        reason[bt & ~bad] = f"T outside {t0:g}-{t1:g} K"
        bad |= bt
        reason[bv & ~bad] = f"value outside {v0:g}-{v1:g}"
        bad |= bv
    dropped = meas[bad].assign(domain_reason=reason[bad])
    for (p, r_), n in dropped.groupby(["property", "domain_reason"]).size().items():
        notes.append(f"domain filter: dropped {n} {p} rows ({r_})")
    return meas[~bad], dropped


# ------------------------------------------------------------------ dedupe
TOL = {"flash_point_C": ("abs", 2.0), "pour_point_C": ("abs", 3.0), "freezing_point_C": ("abs", 1.0),
       "cloud_point_C": ("abs", 3.0), "dcn": ("abs", 3.0), "mon": ("abs", 3.0)}


def dedupe(df: pd.DataFrame) -> pd.DataFrame:
    """One record per (property, system, composition, T).
    1) exact repeats (same source, ref and value) are dropped - NIST/TDE has ~1/3 of those;
    2) the value is the MEDIAN of the remaining independent observations (robust to a single bad paper);
    3) conflict  = the middle half of the observations (IQR) disagree beyond tolerance (or 2 obs disagree);
       n_outlier = observations > 25 % from the median (needs >= 3 obs) - flagged, not silently deleted."""
    df = df.copy()
    df["T_key"] = np.where(df["T_K"].notna(), (df["T_K"] * 2).round() / 2, -1.0)
    P = df["P_kPa"].astype(float)       # missing / atmospheric -> 101; elevated pressures stay separate records
    df["P_key"] = np.where(P.isna() | (P <= 150), 101.0, (P / 100).round() * 100)
    df["x_key"] = [",".join(f"{v:.3f}" for v in json.loads(s)) or fid for s, fid in zip(df["x_mole"], df["fuel_id"])]
    keys = ["property", "system_id", "x_key", "T_key", "P_key"]
    n_all = len(df)
    df = df.drop_duplicates(keys + ["source", "ref", "value"])
    notes.append(f"dedupe: removed {n_all - len(df)} exact repeat rows (same source/ref/value)")
    g = df.groupby(keys, sort=False)
    agg = g.agg(value=("value", "median"), vmin=("value", "min"), vmax=("value", "max"),
                q25=("value", lambda s: s.quantile(.25)), q75=("value", lambda s: s.quantile(.75)),
                n_obs=("value", "size"), n_refs=("ref", "nunique"), n_sources=("source", "nunique"),
                sources=("source", lambda s: "+".join(sorted(set(s)))),
                n_outlier=("value", lambda s: int(((s / s.median() - 1).abs() > 0.25).sum()) if len(s) >= 3 else 0))
    first = df.drop_duplicates(keys).set_index(keys)
    out = first.drop(columns=["value", "source", "raw_value", "raw_unit"]).join(agg).reset_index()

    def conflict(r):
        if r.n_obs < 2:
            return False
        mode, tol = TOL.get(r.property, ("rel", 0.05))
        lo, hi = (r.vmin, r.vmax) if r.n_obs == 2 or mode == "abs" else (r.q25, r.q75)
        spread = (hi - lo) if mode == "abs" else (hi - lo) / max(abs(r.value), 1e-12)
        return bool(spread > tol)
    out["conflict"] = out.apply(conflict, axis=1)
    return out


# ------------------------------------------------------------------ physics inputs
def attach_physics(comp: pd.DataFrame, flash_keys: set[str]) -> pd.DataFrame:
    comp = comp.copy()
    for c in ("antoine_A", "antoine_B", "antoine_C", "antoine_Tlo", "antoine_Thi", "antoine_source"):
        comp[c] = np.nan if c != "antoine_source" else ""
    comp["antoine_covers_window"] = False
    comp["antoine_unit"] = "log10(P/kPa) = A - B/(T/K + C)"
    refit, best = PROC / "antoine_refit.csv", PROC / "antoine_best.csv"
    ant = None
    if refit.exists():
        a = pd.read_csv(refit); a["source"] = "refit_measured"
        ant = a.rename(columns={"T_min": "T_lo", "T_max": "T_hi"})
    if best.exists():
        b = pd.read_csv(best)
        if "smiles" not in b and Path("compounds.csv").exists():
            b = b.merge(pd.read_csv("compounds.csv")[["name", "smiles"]], on="name", how="left")
        if "smiles" not in b:
            notes.append("antoine_best.csv has no smiles and compounds.csv not found: WebBook sets not joined")
            b = b.iloc[0:0].assign(smiles=[])
        b["inchikey"] = [REG.smiles(s) if isinstance(s, str) else None for s in b["smiles"]]
        b["source"] = "webbook"
        b["A"] = b["A"] + 2.0                       # WebBook is log10(P/bar); P[kPa] = 100 * P[bar]  ->  A_kPa = A_bar + 2
        b = b.dropna(subset=["inchikey"])
        ant = b if ant is None else pd.concat([ant, b[~b.inchikey.isin(ant.inchikey)]], ignore_index=True)
    if ant is not None:
        ant = ant.drop_duplicates("inchikey").set_index("inchikey")
        for k in comp.index.intersection(ant.index):
            r = ant.loc[k]
            for col, val in (("antoine_A", r.A), ("antoine_B", r.B), ("antoine_C", r.C), ("antoine_Tlo", r.T_lo),
                             ("antoine_Thi", r.T_hi)):
                comp.loc[k, col] = float(val)
            comp.loc[k, "antoine_source"] = str(r.source)
            comp.loc[k, "antoine_covers_window"] = bool(r.T_lo <= WINDOW[0] and r.T_hi >= WINDOW[1])
    comp["lfl_vol_pct"] = np.nan; comp["lfl_T_C"] = np.nan
    comp["lfl_source"] = pd.Series([""] * len(comp), index=comp.index, dtype=object)
    lfl = RAW / "lfl" / "lfl.csv"
    if lfl.exists():
        d = pd.read_csv(lfl)
        d["inchikey"] = d["inchikey"] if "inchikey" in d else [REG.smiles(s) for s in d["smiles"]]
        d = d.dropna(subset=["inchikey", "lfl_vol_pct"]).drop_duplicates("inchikey").set_index("inchikey")
        for k in comp.index.intersection(d.index):
            comp.loc[k, "lfl_vol_pct"] = float(d.loc[k, "lfl_vol_pct"])
            comp.loc[k, "lfl_T_C"] = float(d.loc[k].get("T_C", np.nan))
            comp.loc[k, "lfl_source"] = str(d.loc[k].get("source", "user")) if "source" in d.columns else "user"
    else:       # write a ready-to-fill template for the compounds that appear in flash-point data
        t = comp.loc[comp.index.isin(flash_keys) & comp.is_hc, ["smiles", "formula", "n_carbon", "cls"]].sort_values(["n_carbon", "formula"])
        t = t.reset_index()
        for c in ("lfl_vol_pct", "T_C", "source", "note"):
            t[c] = ""
        (RAW / "lfl").mkdir(parents=True, exist_ok=True)
        t.to_csv(RAW / "lfl" / "lfl_template.csv", index=False)
        notes.append(f"LFL: wrote data/raw/lfl/lfl_template.csv ({len(t)} compounds). Fill lfl_vol_pct (+T_C, source) from Zabetakis "
                     "Bulletin 627 and save as data/raw/lfl/lfl.csv  -- or simply run:  python fill_lfl.py")
    return comp


# ------------------------------------------------------------------ report
def md(df, index: bool = True) -> str:
    if isinstance(df, pd.Series):
        df = df.to_frame("n")
    d = df.reset_index() if index else df
    f = lambda v: f"{v:.4g}" if isinstance(v, (float, np.floating)) else str(v)
    head = "| " + " | ".join(map(str, d.columns)) + " |\n|" + "---|" * len(d.columns)
    return head + "\n" + "\n".join("| " + " | ".join(f(v) for v in row) + " |" for row in d.itertuples(index=False))


def report(meas: pd.DataFrame, dd: pd.DataFrame, comp: pd.DataFrame, srcs: dict, args) -> str:
    L = ["# Unified dataset report\n"]
    L.append("## Sources loaded\n" + "\n".join(f"- {k}: {v} rows" for k, v in srcs.items()))
    L.append(f"\nRows before dedupe: {len(meas)}   after: {len(dd)}   conflicts (sources disagree): {int(dd.conflict.sum())}")
    L.append(f"Unparsed structures: {REG.unparsed}   compounds in registry: {len(comp)} (hydrocarbons: {int(comp.is_hc.sum())})\n")
    hc = dd[dd.is_hc]
    L.append("## Hydrocarbon-only usable observations (after dedupe)\n")
    piv = hc.pivot_table(index="property", columns="kind", values="value", aggfunc="size", fill_value=0)
    L.append(md(piv) if len(piv) else "none")
    L.append("\n## Ranges\n")
    r = hc.groupby("property").agg(n=("value", "size"), vmin=("value", "min"), vmax=("value", "max"),
                                    Tmin=("T_K", "min"), Tmax=("T_K", "max"), systems=("system_id", "nunique"))
    L.append(md(r.round(3)) if len(r) else "none")
    L.append("\n## Multi-component coverage (n_comp >= 3, hydrocarbon)\n")
    m3 = hc[hc.n_comp >= 3].groupby("property").size()
    L.append(md(m3) if len(m3) else "**none** - no open multi-component hydrocarbon data for any property")
    if (dd.kind == "real_fuel").any():
        L.append("\n## Real fuels (AJFTD)\n" + md(dd[dd.kind == "real_fuel"].groupby("property").size()))
    L.append(f"\nRecords with >=1 outlier observation (>25 % from median, >=3 obs): {int((dd.n_outlier > 0).sum())}; "
             f"records backed by >=3 independent references: {int((dd.n_refs >= 3).sum())}")
    if dd.conflict.any():
        L.append("\nConflicts by source combination:\n" + md(dd[dd.conflict].groupby("sources").size()))
    L.append("\n## Conflicts between sources (largest first)\n")
    cf = dd[dd.conflict].assign(spread=lambda d: d.q75 - d.q25).sort_values("spread", ascending=False).head(10)
    L.append(md(cf[["property", "system_id", "T_K", "n_obs", "n_refs", "sources", "q25", "q75"]], index=False) if len(cf) else "none")

    L.append("\n## Physics-input coverage (what the baselines need)\n")
    def comps_of(df):
        ks = set()
        for s in df["comp_inchikeys"]:
            ks.update(json.loads(s))
        return ks
    fl = hc[(hc.property == "flash_point_C") & (hc.n_comp > 0)]
    if len(fl):
        ks = comps_of(fl)
        c = comp.loc[comp.index.isin(ks)]
        L.append(f"Flash-point records (HC): {len(fl)}; distinct components: {len(c)}")
        L.append(f"- with Antoine set: {int(c.antoine_A.notna().sum())}/{len(c)}; Antoine range covers 40-100 C: {int(c.antoine_covers_window.sum())}")
        L.append(f"- with LFL: {int(c.lfl_vol_pct.notna().sum())}/{len(c)}")
        ok = fl["comp_inchikeys"].apply(lambda s: all(pd.notna(comp.loc[k, "lfl_vol_pct"]) and
                                                       bool(comp.loc[k, "antoine_covers_window"]) for k in json.loads(s)))
        L.append(f"- records where the physics flash-point model can run (all components have LFL + in-range Antoine): **{int(ok.sum())}/{len(fl)}**")
    vi = hc[(hc.property == "viscosity_dyn_mPas") & (hc.n_comp > 0)]
    if len(vi):
        dens = set()
        for s in meas.loc[(meas.property == "density_kg_m3") & (meas.n_comp == 1), "comp_inchikeys"]:
            dens.update(json.loads(s))
        okv = vi["comp_inchikeys"].apply(lambda s: all(k in dens for k in json.loads(s)))
        L.append(f"\nViscosity records (HC): {len(vi)}; components with a measured pure density (needed for kinematic/Refutas): "
                 f"{int(okv.sum())}/{len(vi)} records fully covered")
    L.append("\n## Notes / warnings\n" + ("\n".join(f"- {n}" for n in sorted(set(notes))) or "none"))
    return "\n".join(L)


# ------------------------------------------------------------------ main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--octane", action="store_true")
    ap.add_argument("--all-molecules", action="store_true")
    ap.add_argument("--no-domain-filter", action="store_true", help="keep out-of-domain rows (gas-like densities, GPa pressures ...)")
    a = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)

    rows: list[dict] = []
    srcs: dict[str, int] = {}
    spec = [("flashpoint_pure_exp.csv", "flash_point_C"), ("flashpoint_mix_exp.csv", "flash_point_C"),
            ("viscosity_pure_exp.csv", "viscosity_dyn_mPas"), ("viscosity_mix_exp.csv", "viscosity_dyn_mPas")]
    if a.octane:
        spec += [("dcn_pure_exp.csv", "dcn"), ("dcn_mix_exp.csv", "dcn"), ("mon_pure_exp.csv", "mon"), ("mon_mix_exp.csv", "mon")]
    for fn, prop in spec:
        p = find(INT, fn)
        if p is None:
            notes.append(f"missing {fn} (run inspect_downloads.py)")
            continue
        r = load_leenhouts(p, prop)
        srcs[f"leenhouts:{fn}"] = len(r); rows += r
    for name, fn in (("nist_tde", load_nist_hc), ("thermoml", load_thermoml), ("ajftd", load_ajftd),
                     ("flashpoint_li", load_flashpoint_li), ("fortin2012", load_fortin), ("crc", load_crc),
                     ("chew2024", load_chew2024), ("saldana11", load_saldana11), ("pure_freezing", load_pure_freezing)):
        r = fn(); srcs[name] = len(r); rows += r
        if not r:
            notes.append(f"{name}: no rows (file missing or empty)")
    if not rows:
        raise SystemExit("nothing loaded - run inspect_downloads.py / build_viscosity.py first")

    meas = pd.DataFrame(rows)
    if not a.all_molecules:
        meas = meas[meas.is_hc]
    if not a.no_domain_filter:
        meas, dropped = apply_domain(meas)
        dropped.to_csv(OUT / "measurements_out_of_domain.csv", index=False)
    used = set()
    for sj in meas["comp_inchikeys"]:
        used.update(json.loads(sj))
    comp = pd.DataFrame([v for k, v in REG.by_key.items() if k in used]).set_index("inchikey")
    flash_keys = set()
    for s in meas.loc[meas.property == "flash_point_C", "comp_inchikeys"]:
        flash_keys.update(json.loads(s))
    comp = attach_physics(comp, flash_keys)
    dd = dedupe(meas)

    comp.reset_index().to_csv(OUT / "compounds.csv", index=False)
    meas.to_csv(OUT / "measurements.csv", index=False)
    dd.to_csv(OUT / "measurements_dedup.csv", index=False)
    rep = report(meas, dd, comp, srcs, a)
    (OUT / "UNIFIED_REPORT.md").write_text(rep, encoding="utf-8")
    print(rep)
    print(f"\nwrote {OUT}/ (compounds.csv, measurements.csv, measurements_dedup.csv, UNIFIED_REPORT.md)")


if __name__ == "__main__":
    main()
