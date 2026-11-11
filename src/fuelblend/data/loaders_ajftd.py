"""AJFTD (Alternative Jet Fuel Test Database) bulk-descriptor loader (v5, task D1).

Reads the on-disk scrape ``extras/fuelblend_data/data/raw/ajftd/ajftd_all_metrics.csv`` (one row per
site table cell) and writes

* ``data/real_fuels/v5/ajftd_fuels_long.csv``  - numeric rows only, one per measurement
* ``data/real_fuels/v5/ajftd_fuels_wide.csv``  - one row per fuel (median per property)
* ``data/real_fuels/v5/ajftd_build_log.json``  - drop counts and coverage

Provenance class is always ``real_fuel_bulk_descriptors``: these are batch-level property tests of real
fuels with NO molecule list; hydrocarbon-type composition exists only for a handful of fuels (see the
D1 note in ``docs/ajftd_v5_note.md``).

Encodings found in the raw file (verified):
* kinematic viscosity temperature is in ``temperature_value`` (deg C, text; ``temperature_unit`` = C).
  It is ``No Test``/blank for ~590 numeric rows; those are the spec-table entry of the test with method
  ``310`` (DEF STAN 91-091 table number) or ``D445`` whose limit is ``8 Max`` (= at -20 C, the Jet A / Jet A-1
  specification condition).  They are loaded at T = 253.15 K with ``T_source = 'assumed_spec_-20C'``.
* method codes: ``D445``/``D 445`` (ASTM), ``310``/``310IPA`` (spec-table id, not an ASTM number),
  ``D2532`` (-40 C), ``D7042``.  ``measurement_index`` repeats the same value for replicates/spec rows;
  exact duplicates are dropped.
* density ``Unit`` kg/m3 with temperature 15 C; some values are g/mL (<2) -> x1000.
* ``density_api_by_meter`` rows with Unit ``API`` are API gravity; rows with Unit ``kg/m3`` are in fact
  net heat of combustion in kJ/kg (values ~44 000) - re-routed, flagged ``mislabelled_metric``.
"""
from __future__ import annotations

import json
import re
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[3]
RAW = REPO / "extras" / "fuelblend_data" / "data" / "raw" / "ajftd" / "ajftd_all_metrics.csv"
OUT = REPO / "data" / "real_fuels" / "v5"
PROVENANCE = "real_fuel_bulk_descriptors"

NOT_A_VALUE = {"", "no test", "no search", "not found", "nan", "none", "n/a", "na", "-", "see below"}
# plausibility windows (value dropped, counted in the log, if outside)
PLAUSIBLE = {
    "flash_point_C": (-10.0, 200.0),
    "freezing_point_C": (-110.0, 0.0),
    "viscosity_kin_cSt": (0.2, 40.0),
    "aromatics_vol_pct": (0.0, 100.0),
    "density_kg_m3": (600.0, 1000.0),
    "heat_of_combustion_net_MJ_kg": (40.0, 46.0),
}
LONG_COLS = ["fuel_id", "fuel_label", "source_group", "test_year", "test_date", "category", "property", "value",
             "unit", "T_K", "T_source", "method", "qualifier", "info", "derived", "provenance_class"]


def api_to_density_kg_m3(api: float) -> float:
    """API gravity (60 F) -> density at 15.56 C, kg/m3 (D1298 definition)."""
    return 141.5 / (api + 131.5) * 999.016


def parse_value(raw) -> tuple[float, str]:
    """-> (number or NaN, qualifier '<' / '>' / '').  Non-values ('No Test' ...) give NaN."""
    s = str(raw).strip()
    if s.lower() in NOT_A_VALUE:
        return np.nan, ""
    q = s[0] if s[:1] in "<>" else ""
    try:
        return float(re.sub(r"[<>≤≥~,\s]", "", s)), q
    except ValueError:
        return np.nan, q


def parse_temp_c(raw) -> float:
    """Temperature text -> deg C.  '-20' -> -20 ; '-20°C (-4°F),,,-40' -> -20 (first number) ; 'No Test' -> NaN."""
    m = re.match(r"\s*(-?\d+(?:\.\d+)?)", str(raw))
    return float(m.group(1)) if m else np.nan


def source_group(label: str) -> str:
    """Normalised prefix of the fuel label (the site's naming convention encodes the supplier/programme)."""
    s = str(label).strip()
    u = s.upper()
    if u.startswith("JF KENNEDY"):
        return "JFK"
    if u.startswith("JETSCREEN"):
        return "JETSCREEN"
    if u.startswith("CRC"):
        return "CRC"
    if re.match(r"HR\d", u):
        return "HR"
    if u.startswith("JAA"):
        return "JAA"
    if re.match(r"JA\d", u) or u.startswith("JA1"):
        return "JA"
    if s.startswith("("):
        return "PAREN_BATCH"
    if re.fullmatch(r"\d+", s):
        return "NUMERIC_ID"
    if re.match(r"\d{7}-\d", s) or re.match(r"\d{3}-\d{7}-\d", s):
        return "LAB_REPORT"
    return "OTHER"


def _test_year(date) -> float:
    m = re.search(r"/(\d{4})\b", str(date))
    return float(m.group(1)) if m else np.nan


def _unit_ok(unit: str, allowed: tuple[str, ...]) -> bool:
    return str(unit).strip().lower().replace(" ", "") in allowed or str(unit).strip() in ("", "nan")


def build_long(raw: pd.DataFrame, log: dict | None = None) -> pd.DataFrame:
    """Raw scrape table -> canonical long table (numeric rows only, exact duplicates dropped)."""
    log = {} if log is None else log
    dropped: dict[str, int] = log.setdefault("dropped", {})
    rows: list[dict] = []

    def drop(reason):
        dropped[reason] = dropped.get(reason, 0) + 1

    for r in raw.to_dict("records"):
        metric = str(r["metric"])
        wanted = metric in ("flash_point", "freezing_point", "viscosity_kinematic", "aromatics_total",
                            "distillation", "density", "density_api_by_meter", "net_heat_of_combustion")
        if not wanted:
            continue
        v, q = parse_value(r["Value"])
        if np.isnan(v):
            continue
        unit = str(r.get("Unit", "")).strip()
        meth = "" if pd.isna(r.get("Test Method")) else str(r.get("Test Method")).strip()
        base = dict(fuel_id=r["fuel_id"], fuel_label=r["fuel_label"], source_group=source_group(r["fuel_label"]),
                    test_year=_test_year(r.get("Test Date")), test_date=str(r.get("Test Date", "")),
                    category=r["category"], method=meth, qualifier=q,
                    info="" if pd.isna(r.get("info")) else str(r.get("info")),
                    derived="", T_K=np.nan, T_source="", provenance_class=PROVENANCE)
        tc = parse_temp_c(r.get("temperature_value"))
        out = None
        if metric == "flash_point":
            if unit.upper() not in ("C", "", "NAN"):
                drop("flash_unit_not_C"); continue
            out = dict(property="flash_point_C", value=v, unit="C")
        elif metric == "freezing_point":
            out = dict(property="freezing_point_C", value=v, unit="C")
        elif metric == "viscosity_kinematic":
            if np.isnan(tc):
                T_source, tc = "assumed_spec_-20C", -20.0
            elif str(r.get("temperature_value")).strip().lstrip("-").replace(".", "", 1).isdigit():
                T_source = "reported"
            else:
                T_source = "reported_parsed_text"
            base["T_K"] = tc + 273.15
            base["T_source"] = T_source
            out = dict(property="viscosity_kin_cSt", value=v, unit="cSt")      # mm2/s == cSt
        elif metric == "aromatics_total":
            if "mass" in unit.lower():
                out = dict(property="aromatics_mass_pct", value=v, unit="mass%")
            else:
                out = dict(property="aromatics_vol_pct", value=v, unit="vol%")
        elif metric == "distillation":
            vol = pd.to_numeric(r.get("volume_evaporated_value"), errors="coerce")
            if pd.isna(vol):
                drop("distillation_no_volume"); continue
            if meth and not re.search(r"(?:^|\D)(86|200[AB]?)(?:\D|$)|200IP", meth.replace(" ", "")):
                drop("distillation_non_D86_method"); continue
            out = dict(property=f"distillation_T{vol:g}_C", value=v, unit="C")
        elif metric == "density":
            if v < 2.0:
                v *= 1000.0
                base["info"] = (base["info"] + " g/mL->kg/m3").strip()
            if np.isnan(tc):
                base["T_K"], base["T_source"] = 288.15, "assumed_15C"
            else:
                base["T_K"], base["T_source"] = tc + 273.15, "reported"
            out = dict(property="density_kg_m3", value=v, unit="kg/m3")
        elif metric == "density_api_by_meter":
            if unit.upper() == "API":
                out = dict(property="density_kg_m3", value=api_to_density_kg_m3(v), unit="kg/m3")
                base.update(T_K=288.706, T_source="API_def_60F", derived="from API gravity")
            else:                                  # mislabelled heat of combustion in kJ/kg
                out = dict(property="heat_of_combustion_net_MJ_kg", value=v / 1000.0, unit="MJ/kg")
                base["derived"] = "mislabelled_metric(kJ/kg in density_api_by_meter)"
        elif metric == "net_heat_of_combustion":
            if v > 1000:
                v /= 1000.0
                base["derived"] = "kJ/kg->MJ/kg"
            out = dict(property="heat_of_combustion_net_MJ_kg", value=v, unit="MJ/kg")
        if out is None:
            continue
        lo_hi = PLAUSIBLE.get(out["property"])
        if lo_hi and not (lo_hi[0] <= out["value"] <= lo_hi[1]):
            drop(f"implausible_{out['property']}"); continue
        rows.append({**base, **out})

    df = pd.DataFrame(rows, columns=LONG_COLS)
    n0 = len(df)
    df["value"] = df["value"].round(6)
    df = df.drop_duplicates(subset=[c for c in LONG_COLS if c not in ()]).reset_index(drop=True)
    dropped["exact_duplicate_rows"] = n0 - len(df)
    return df[LONG_COLS]


def build_wide(long: pd.DataFrame) -> pd.DataFrame:
    """One row per fuel: median per property (censored '<'/'>' values excluded).  Viscosity columns are
    per temperature (-20 C, 40 C); density is at 15 C only (direct preferred over API-derived)."""
    L = long[long["qualifier"] == ""].copy()
    meta = long.drop_duplicates("fuel_id").set_index("fuel_id")[["fuel_label", "source_group", "test_year"]]
    cols = {}
    simple = L[~L["property"].isin(["viscosity_kin_cSt", "density_kg_m3"])]
    for p, g in simple.groupby("property"):
        cols[p] = g.groupby("fuel_id")["value"].median()
    v = L[L["property"] == "viscosity_kin_cSt"].copy()
    v["Tc"] = (v["T_K"] - 273.15).round(0).astype(int)
    for tc, name in [(-20, "viscosity_kin_cSt_m20C"), (-40, "viscosity_kin_cSt_m40C"), (40, "viscosity_kin_cSt_40C")]:
        g = v[v["Tc"] == tc]
        if len(g):
            cols[name] = g.groupby("fuel_id")["value"].median()
    d = L[(L["property"] == "density_kg_m3") & (L["T_K"].sub(288.15).abs() < 1.0)]
    direct = d[d["derived"] == ""].groupby("fuel_id")["value"].median()
    derived = d[d["derived"] != ""].groupby("fuel_id")["value"].median()
    dens = direct.combine_first(derived)
    cols["density_kg_m3_15C"] = dens
    cols["density_is_derived"] = (~dens.index.isin(direct.index)).astype(bool)
    cols["density_is_derived"] = pd.Series(cols["density_is_derived"], index=dens.index)
    wide = meta.join(pd.DataFrame(cols), how="left")
    wide["provenance_class"] = PROVENANCE
    return wide.reset_index().sort_values("fuel_id").reset_index(drop=True)


CORE5 = ["flash_point_C", "freezing_point_C", "viscosity_kin_cSt", "aromatics_vol_pct", "distillation"]


def n_core5(long: pd.DataFrame) -> int:
    """Fuels having flash + freezing + viscosity + aromatics (vol%) + at least one D86 point."""
    ok = None
    for p in CORE5:
        m = long["property"].str.startswith("distillation_T") if p == "distillation" else long["property"] == p
        s = set(long.loc[m, "fuel_id"])
        ok = s if ok is None else ok & s
    return len(ok)


def composition_summary(raw: pd.DataFrame) -> dict:
    """Which hydrocarbon-type composition rows carry numbers?  (Composition Metrics, numeric only.)"""
    c = raw[raw["category"] == "Composition Metrics"].copy()
    c["v"] = pd.to_numeric(c["Value"], errors="coerce")
    out = {"composition_rows_total": int(len(c)), "composition_rows_numeric": int(c["v"].notna().sum())}
    for m in ["aromatics_total", "saturates_total", "monoaromatics", "diaromatics", "olefins_total",
              "paraffins_total", "cycloalkanes_total", "cycloparaffins", "naphthalenes", "h_c_ratio",
              "hydrogen_content", "average_molecular_weight"]:
        s = c[(c["metric"] == m) & c["v"].notna()]
        out[m] = {"rows": int(len(s)), "fuels": int(s["fuel_id"].nunique())}
    sets = {m: set(c[(c["metric"] == m) & c["v"].notna()]["fuel_id"])
            for m in ["saturates_total", "aromatics_total", "paraffins_total", "cycloalkanes_total", "olefins_total"]}
    out["fuels_sat+arom"] = len(sets["saturates_total"] & sets["aromatics_total"])
    out["fuels_paraffin+cycloalkane+arom"] = len(sets["paraffins_total"] & sets["cycloalkanes_total"] & sets["aromatics_total"])
    return out


def build_ajftd_v5(raw_path: Path = RAW, out_dir: Path = OUT) -> dict:
    raw = pd.read_csv(raw_path, encoding="utf-8-sig", low_memory=False)
    log: dict = {"raw_rows": int(len(raw)), "raw_fuels": int(raw["fuel_id"].nunique())}
    long = build_long(raw, log)
    wide = build_wide(long)
    out_dir.mkdir(parents=True, exist_ok=True)
    long.to_csv(out_dir / "ajftd_fuels_long.csv", index=False)
    wide.to_csv(out_dir / "ajftd_fuels_wide.csv", index=False)
    log["long_rows"] = int(len(long))
    log["rows_per_property"] = {k: int(v) for k, v in long.groupby("property").size().items()}
    log["fuels_per_property"] = {k: int(v) for k, v in long.groupby("property")["fuel_id"].nunique().items()}
    log["fuels_with_numeric"] = int(long["fuel_id"].nunique())
    log["fuels_core5"] = n_core5(long)
    log["source_group_fuels"] = {k: int(v) for k, v in wide["source_group"].value_counts().items()}
    log["test_year_range"] = [float(long["test_year"].min()), float(long["test_year"].max())]
    log["composition"] = composition_summary(raw)
    (out_dir / "ajftd_build_log.json").write_text(json.dumps(log, indent=2), encoding="utf-8")
    return log


def main() -> None:
    print(json.dumps(build_ajftd_v5(), indent=2))


if __name__ == "__main__":
    main()
