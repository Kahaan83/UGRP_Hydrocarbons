#!/usr/bin/env python
"""Ingest the AJFTD export(s) saved in data/raw/ajftd/  (.csv / .xlsx / .xls).

    python ingest_ajftd.py --describe   # writes data/interim/ajftd_description.txt  <- upload THIS to me first
    python ingest_ajftd.py              # best-effort mapping to the canonical schema

Outputs (data/processed/):
  ajftd_fuels.csv         long: fuel_id, category, property, value, unit_std, T_K, method, src_file, src_col
  ajftd_composition.csv   wide: fuel_id + GC / GCxGC class columns (n-/iso-/cyclo-paraffin, aromatics, olefins ...)
  data/interim/ajftd_unmapped.txt   columns / property names that were NOT mapped (so nothing disappears silently)

FORMATS: (1) the output of YOUR ajftd_scraper.py (columns fuel_id, metric, Value, Unit, temperature_value,
volume_evaporated_value, Test Method ...) is recognised automatically and parsed by `ingest_scraper`;
(2) generic WIDE / LONG exports are handled by best-effort regex rules.
"No Test" / "No search" / "Not Found" / blank values are dropped (and counted in the printout).

STATUS (generic layouts): the real export format was not available when this was written. Two layouts are handled
(WIDE: one row per fuel, one column per property;  LONG: one row per measurement) and the matching
rules are plain regular expressions in PROPS / COMP below - when --describe shows your actual
headers, edit those lists (or send me the description and I will).
"""
from __future__ import annotations
import re, sys
from pathlib import Path

import numpy as np
import pandas as pd

SRC, PROC, INT = Path("data/raw/ajftd"), Path("data/processed"), Path("data/interim")

# canonical property, header regex, default unit, kind ('temp' -> values are temperatures)
PROPS = [
    ("flash_point_C", r"flash", "C", "temp"),
    ("freezing_point_C", r"freez", "C", "temp"),
    ("pour_point_C", r"pour", "C", "temp"),
    ("cloud_point_C", r"cloud", "C", "temp"),
    ("viscosity_kin_cSt", r"kinematic|viscos(?!.*dynamic)", "cSt", "visc"),
    ("density_kg_m3", r"densit", "kg/m3", "dens"),
    ("net_heat_of_combustion_MJ_kg", r"net.*(heat|energy)|heat of comb", "MJ/kg", "num"),
    ("hydrogen_mass_pct", r"hydrogen", "%", "num"),
    ("smoke_point_mm", r"smoke", "mm", "num"),
    ("distillation_ibp_C", r"\bibp\b|initial boiling", "C", "temp"),
    ("distillation_t10_C", r"\bt ?10\b|10 ?% ?(rec|dist|evap)", "C", "temp"),
    ("distillation_t50_C", r"\bt ?50\b|50 ?% ?(rec|dist|evap)", "C", "temp"),
    ("distillation_t90_C", r"\bt ?90\b|90 ?% ?(rec|dist|evap)", "C", "temp"),
    ("distillation_fbp_C", r"\bfbp\b|final boiling|end ?point", "C", "temp"),
]
COMP = r"n[- ]?paraffin|normal paraffin|iso[- ]?paraffin|cyclo[- ]?paraffin|naphthen|cycloalkane|aromatic|olefin|alkene|alkane|paraffin|gc ?x ?gc|carbon number"
ID_PAT = r"(fuel.*(id|name|code)|sample.*(id|name)|posf|^id$|^name$|^fuel$)"
CAT_PAT = r"(categor|fuel ?type|class|feedstock|pathway|description|origin|blend)"
LONG_PROP = r"(property|test|analysis|characteristic|parameter|determination)"
LONG_VAL = r"(^value$|result|measurement)"


def read_tables() -> dict[str, pd.DataFrame]:
    out = {}
    for p in sorted(SRC.rglob("*")):
        if not p.is_file() or p.name.startswith(("._", "~$")) or "__MACOSX" in p.parts:
            continue
        try:
            if p.suffix.lower() in (".csv", ".tsv", ".txt"):
                sep = "\t" if p.suffix.lower() == ".tsv" else None
                try:
                    out[p.name] = pd.read_csv(p, sep=sep, engine="python", encoding="utf-8-sig")
                except UnicodeDecodeError:
                    out[p.name] = pd.read_csv(p, sep=sep, engine="python", encoding="latin-1")
            elif p.suffix.lower() in (".xlsx", ".xls"):
                for sh, df in pd.read_excel(p, sheet_name=None).items():
                    out[f"{p.name}::{sh}"] = df
        except Exception as e:  # noqa: BLE001
            print(f"could not read {p.name}: {type(e).__name__}: {str(e)[:100]}")
    return out


# ------------------------------------------------------------------ describe
def describe(tables):
    L = ["AJFTD export description (structure + a few sample values; delete anything you must not share)\n"]
    for name, df in tables.items():
        L.append(f"\n=== {name}: {df.shape[0]} rows x {df.shape[1]} cols")
        for c in df.columns:
            s = df[c].dropna()
            uniq = s.astype(str).unique()
            sample = list(uniq[:25]) if (len(uniq) <= 25 and s.dtype == object) else list(s.astype(str).head(3))
            L.append(f"- {str(c)!r:60} dtype={str(df[c].dtype):8} non-null={len(s):6} unique={len(uniq):6} sample={sample}")
    INT.mkdir(parents=True, exist_ok=True)
    (INT / "ajftd_description.txt").write_text("\n".join(L), encoding="utf-8")
    print("\n".join(L[:60]))
    print(f"\nwrote {INT/'ajftd_description.txt'}  - upload this file")


# ------------------------------------------------------------------ helpers
def find_col(cols, pat):
    return next((c for c in cols if re.search(pat, str(c), re.I)), None)


def num(v):
    if isinstance(v, str):
        v = re.sub(r"[<>≤≥~]", "", v).replace(",", "").strip()
    return pd.to_numeric(v, errors="coerce")


def temp_from(text: str):
    """'Viscosity @ -20 C', 'density at 15°C', 'at 313.15 K' -> Kelvin."""
    m = re.search(r"(-?\d+(?:\.\d+)?)\s*(?:°|º|deg(?:rees)?\.?\s*)?\s*([CFK])\b", str(text))
    if not m:
        return np.nan
    v, u = float(m.group(1)), m.group(2).upper()
    return v + 273.15 if u == "C" else (v - 32) * 5 / 9 + 273.15 if u == "F" else v


def unit_from(text: str, default: str) -> str:
    m = re.findall(r"\(([^()]*)\)", str(text))
    return m[-1].strip() if m else default


def convert(prop: str, kind: str, v: float, unit: str):
    u = unit.lower().replace(" ", "").replace("²", "2").replace("³", "3")
    if kind == "temp":
        if u in ("f", "°f", "degf"):
            return (v - 32) * 5 / 9
        if u in ("k", "kelvin"):
            return v - 273.15
        return v
    if kind == "visc":
        if "m2/s" in u and "mm2" not in u:
            return v * 1e6
        return v                                            # cSt == mm2/s
    if kind == "dens":
        if u in ("g/ml", "g/cm3", "g/cc", "gm/ml", "kg/l"):
            return v * 1000
        return v
    return v


def map_prop(text: str):
    t = str(text)
    if re.search(COMP, t, re.I):
        return None
    for prop, pat, unit, kind in PROPS:
        if re.search(pat, t, re.I):
            return prop, unit, kind
    return None


# ------------------------------------------------------------------ layouts
def ingest_wide(name, df, rows, comp_rows, unmapped):
    idc = find_col(df.columns, ID_PAT) or df.columns[0]
    catc = find_col([c for c in df.columns if c != idc], CAT_PAT)
    for c in df.columns:
        if c in (idc, catc):
            continue
        hit = map_prop(c)
        if hit is None:
            if re.search(COMP, str(c), re.I):
                comp_rows.setdefault(name, df[[idc]].copy().rename(columns={idc: "fuel_id"}))[f"comp::{c}"] = df[c].map(num)
            else:
                nn = int(df[c].notna().sum())
                if nn:
                    unmapped.append(f"{name}: column {c!r} ({nn} non-null)")
            continue
        prop, dunit, kind = hit
        unit = unit_from(c, dunit)
        T = temp_from(c) if kind in ("visc", "dens") else np.nan
        for fid, cat, v in zip(df[idc], df[catc] if catc else [""] * len(df), df[c].map(num)):
            if pd.notna(v):
                rows.append(dict(fuel_id=fid, category=cat, property=prop, value=convert(prop, kind, float(v), unit),
                                 unit_std=dunit, T_K=T, method="", src_file=name, src_col=str(c)))


def ingest_long(name, df, rows, unmapped):
    idc = find_col(df.columns, ID_PAT) or df.columns[0]
    pc, vc = find_col(df.columns, LONG_PROP), find_col(df.columns, LONG_VAL)
    catc = find_col([c for c in df.columns if c not in (idc, pc, vc)], CAT_PAT)
    uc, mc = find_col(df.columns, r"^unit"), find_col(df.columns, r"method|astm|standard")
    tc = find_col(df.columns, r"temp|condition")
    miss: dict[str, int] = {}
    for r in df.to_dict("records"):
        hit = map_prop(r[pc])
        v = num(r[vc])
        if hit is None or pd.isna(v):
            if hit is None and pd.notna(v):
                miss[str(r[pc])] = miss.get(str(r[pc]), 0) + 1
            continue
        prop, dunit, kind = hit
        unit = str(r[uc]) if uc and pd.notna(r[uc]) else unit_from(r[pc], dunit)
        T = temp_from(f"{r[pc]} {r[tc]}") if kind in ("visc", "dens") else np.nan
        if kind in ("visc", "dens") and np.isnan(T) and tc and pd.notna(r[tc]):
            t = num(r[tc]); T = t + 273.15 if pd.notna(t) else np.nan      # bare number in a temperature column: deg C
        rows.append(dict(fuel_id=r[idc], category=r[catc] if catc else "", property=prop, value=convert(prop, kind, float(v), unit),
                         unit_std=dunit, T_K=T, method=r[mc] if mc else "", src_file=name, src_col=str(r[pc])))
    unmapped += [f"{name}: property {k!r} ({n} values)" for k, n in sorted(miss.items(), key=lambda kv: -kv[1])]


# ------------------------------------------------------------------ scraper format (ajftd_scraper.py output)
NOT_A_VALUE = {"", "no test", "no search", "not found", "nan", "none", "n/a", "na", "-"}
SCRAPER_METRICS = {   # normalised metric label -> (canonical property, default unit, kind)
    "density": ("density_kg_m3", "kg/m3", "dens"),
    "viscosity_kinematic": ("viscosity_kin_cSt", "cSt", "visc"),
    "kinematic_viscosity": ("viscosity_kin_cSt", "cSt", "visc"),
    "viscosity": ("viscosity_dyn_mPas", "mPa*s", "num"),
    "flash_point": ("flash_point_C", "C", "temp"),
    "pour_point": ("pour_point_C", "C", "temp"),
    "freezing_point": ("freezing_point_C", "C", "temp"),
    "cloud_point": ("cloud_point_C", "C", "temp"),
    "smoke_point": ("smoke_point_mm", "mm", "num"),
    "distillation_loss": ("distillation_loss_pct", "%", "num"),
    "distillation_recovery": ("distillation_recovery_pct", "%", "num"),
    "distillation_residue": ("distillation_residue_pct", "%", "num"),
    "thermal_stability_spent_fuel_volume": ("thermal_stability_spent_fuel_mL", "mL", "num"),
}
GENERIC_NAMES = {"alternative jet fuels", "", "nan", "none"}


def is_scraper_format(df) -> bool:
    cols = {str(c).strip().lower() for c in df.columns}
    return {"fuel_id", "metric", "value", "unit"} <= cols


def api_to_density_kg_m3(api: float) -> float:
    """API gravity (60 F) -> density at 15.56 C: rho = 141.5/(API+131.5) * 999.016 kg/m3 (D1298 definition)."""
    return 141.5 / (api + 131.5) * 999.016


def ingest_scraper(name, df, rows, comp_rows, unmapped, report):
    df = df.copy()
    df.columns = [str(c).strip() for c in df.columns]
    lc = {c.lower(): c for c in df.columns}
    g = lambda r, k: r.get(lc.get(k.lower(), "\0"))
    comp: dict[str, dict] = {}
    for r in df.to_dict("records"):
        fid = str(g(r, "fuel_id"))
        label = g(r, "fuel_label") or g(r, "fuel_name")
        label = fid if str(label).strip().lower() in GENERIC_NAMES else str(label)
        raw = str(g(r, "Value")).strip()
        metric = re.sub(r"[^a-z0-9]+", "_", re.sub(r"\([^)]*\)", " ", str(g(r, "metric")).lower())).strip("_")
        report[metric]["rows"] += 1
        if raw.lower() in NOT_A_VALUE:
            report[metric]["no_value"] += 1
            continue
        qual = raw[0] if raw[:1] in "<>" else ""
        v = num(raw)
        if pd.isna(v):
            unmapped.append(f"{name}: non-numeric Value {raw!r} for metric {metric!r} (fuel {fid})")
            continue
        unit = str(g(r, "Unit") or "").strip()
        meth, date = str(g(r, "Test Method") or ""), str(g(r, "Test Date") or "")
        t_raw = num(g(r, "temperature_value")) if g(r, "temperature_value") is not None else np.nan
        t_unit = str(g(r, "temperature_unit") or "C").strip().upper()[:1] or "C"
        T = (t_raw + 273.15 if t_unit == "C" else (t_raw - 32) * 5 / 9 + 273.15 if t_unit == "F" else t_raw) if pd.notna(t_raw) else np.nan
        base = dict(fuel_id=fid, fuel_label=label, category=g(r, "category"), method=meth, test_date=date, qualifier=qual,
                    info=g(r, "info"), src_file=name, src_col=metric)
        out = None
        if metric == "density_api_by_meter":                          # API gravity -> density at 15.56 C (derived)
            out = dict(property="density_kg_m3", value=api_to_density_kg_m3(v), unit_std="kg/m3", T_K=288.706, derived="from API gravity")
        elif metric == "distillation":
            vol = num(g(r, "volume_evaporated_value"))
            if pd.notna(vol):
                out = dict(property=f"distillation_T{vol:g}_C", value=convert("", "temp", v, unit or "C"), unit_std="C", T_K=np.nan)
        elif metric in SCRAPER_METRICS:
            prop, dunit, kind = SCRAPER_METRICS[metric]
            out = dict(property=prop, value=convert(prop, kind, v, unit or dunit), unit_std=dunit,
                       T_K=T if kind in ("visc", "dens") else np.nan)
        else:                                                           # composition / contaminants / acidity / sulfur ...
            col = f"comp::{metric}[{unit}]" if unit else f"comp::{metric}"
            comp.setdefault(fid, {"fuel_id": fid, "fuel_label": label}).setdefault(col, []).append(v)
            report[metric]["numeric"] += 1
            continue
        if out is None:
            unmapped.append(f"{name}: distillation row without volume % (fuel {fid})")
            continue
        rows.append({**base, **out})
        report[metric]["numeric"] += 1
    if comp:
        wide = []
        for fid, d in comp.items():
            row = {k: (float(np.median(v)) if isinstance(v, list) else v) for k, v in d.items()}
            wide.append(row)
        comp_rows[name] = pd.DataFrame(wide)


def main():
    tables = read_tables()
    if not tables:
        sys.exit(f"No files found in {SRC}/ - save the AJFTD export(s) (csv/xlsx) there first.")
    if "--describe" in sys.argv:
        return describe(tables)
    rows, comp_rows, unmapped = [], {}, []
    import collections
    report: dict = collections.defaultdict(lambda: collections.Counter())
    for name, df in tables.items():
        df = df.dropna(how="all").dropna(axis=1, how="all")
        if df.empty:
            continue
        if is_scraper_format(df):
            print(f"{name}: SCRAPER layout (ajftd_scraper.py output)")
            ingest_scraper(name, df, rows, comp_rows, unmapped, report)
        elif find_col(df.columns, LONG_PROP) and find_col(df.columns, LONG_VAL):
            print(f"{name}: LONG layout"); ingest_long(name, df, rows, unmapped)
        else:
            print(f"{name}: WIDE layout"); ingest_wide(name, df, rows, comp_rows, unmapped)
    PROC.mkdir(parents=True, exist_ok=True); INT.mkdir(parents=True, exist_ok=True)
    out = pd.DataFrame(rows)
    out.to_csv(PROC / "ajftd_fuels.csv", index=False)
    if comp_rows:
        pd.concat(comp_rows.values(), ignore_index=True).to_csv(PROC / "ajftd_composition.csv", index=False)
    (INT / "ajftd_unmapped.txt").write_text("\n".join(unmapped) or "nothing unmapped", encoding="utf-8")
    print(f"\nfuels: {out.fuel_id.nunique() if len(out) else 0}   measurements: {len(out)}")
    if len(out):
        print(out.groupby("property").agg(n=("value", "size"), fuels=("fuel_id", "nunique"), vmin=("value", "min"), vmax=("value", "max")).round(3).to_string())
    if report:
        print("\ncoverage per metric (rows on the site vs rows with a real number):")
        for m_, c in sorted(report.items(), key=lambda kv: -kv[1]["rows"]):
            print(f"  {m_:48} rows={c['rows']:5}  numeric={c['numeric']:5}  no value ('No Test' etc.)={c['no_value']:5}")
    print(f"unmapped items: {len(unmapped)} (see data/interim/ajftd_unmapped.txt)")


if __name__ == "__main__":
    main()
