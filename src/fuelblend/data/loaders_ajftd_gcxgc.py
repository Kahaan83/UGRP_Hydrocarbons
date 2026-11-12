"""AJFTD GCxGC hydrocarbon-type compositions (class x carbon number, weight %), v5 data audit.

Source: extras/fuelblend_data/data/raw/ajftd/ajftd_fuel_details.csv, the per-fuel detail tables saved by the scraper
(long format: fuel_id, fuel_label, table_index, row, column, value).  Table 0 of a fuel page is the composition tree; a
leaf row has a `Metric` such as `n-C10`, `C11-isoparaffins`, `C3-benzene (C09)` and a numeric `Weight %`.
docs/ajftd_v5_note.md looked only at ajftd_all_metrics.csv and therefore reported "no hydrocarbon-type composition".

Outputs (data/real_fuels/v5/):
  ajftd_gcxgc_composition.csv : one row per fuel x class x carbon number (weight %, as printed; not renormalised)
  ajftd_gcxgc_fuels.csv       : one row per fuel: leaf sum, class totals, number of bins, which bulk properties exist
  ajftd_gcxgc_build_log.json
Provenance class: real_fuel_gcxgc_bins (measured class x carbon-number bins; no molecule list, isomers unresolved).
AJFTD is login-gated: do not publish the raw rows.
"""
from __future__ import annotations

import json
import re
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
RAW = ROOT / "extras" / "fuelblend_data" / "data" / "raw" / "ajftd" / "ajftd_fuel_details.csv"
OUT = ROOT / "data" / "real_fuels" / "v5"
WIDE = OUT / "ajftd_fuels_wide.csv"
CLOSURE = (98.0, 102.0)                      # leaf weight % sum accepted as a complete composition
TARGETS = ("flash_point_C", "freezing_point_C", "viscosity_kin_cSt", "density_kg_m3", "heat_of_combustion_net_MJ_kg")

# order matters: the first keyword found decides the class
_CLASSES = (("tricyclo", "tricycloparaffin"), ("dicyclo", "dicycloparaffin"), ("monocyclo", "monocycloparaffin"),
            ("cycloaromatic", "cycloaromatic"), ("indene", "cycloaromatic"), ("diaromatic", "diaromatic"),
            ("naphthalene", "diaromatic"), ("decalin", "dicycloparaffin"), ("benzene", "alkylbenzene"),
            ("toluene", "alkylbenzene"), ("diene", "diene"), ("alkene", "alkene"), ("iso", "iso_paraffin"))
_CN = re.compile(r"C(\d{2})(\+)?")


def parse_metric(metric: str) -> tuple[str, int, bool] | None:
    """`C3-benzene (C09)` -> ('alkylbenzene', 9, False); `n-C10` -> ('n_paraffin', 10, False); a class total -> None."""
    m = str(metric).strip()
    low = m.lower()
    paren = re.search(r"\(C(\d{2})(\+)?\)", m)             # alkylbenzenes carry the true carbon number in brackets
    hit = paren or _CN.search(m)
    if not hit:
        return None
    cn, open_ended = int(hit.group(1)), bool(hit.group(2)) or "lower" in low or "+" in m
    if re.fullmatch(r"n-C\d{2}(\+| (and|&) lower)?", m):
        return "n_paraffin", cn, open_ended
    for key, cls in _CLASSES:
        if key in low:
            return cls, cn, open_ended
    return None


def load(raw: Path = RAW) -> tuple[pd.DataFrame, list[str]]:
    d = pd.read_csv(raw, low_memory=False)
    t = d[(d.table_index == 0) & d["column"].isin(["Metric", "Weight %"])]
    w = t.pivot_table(index=["fuel_id", "fuel_label", "row"], columns="column", values="value", aggfunc="first").reset_index()
    if "Weight %" not in w or "Metric" not in w:
        return pd.DataFrame(), []
    w["mass_pct"] = pd.to_numeric(w["Weight %"], errors="coerce")
    w = w[w.mass_pct.notna() & w.Metric.notna() & (w.Metric.astype(str).str.len() < 60)]
    rows, unparsed = [], set()
    for r in w.itertuples():
        p = parse_metric(r.Metric)
        if p is None:
            unparsed.add(str(r.Metric))
            continue
        rows.append(dict(fuel_id=r.fuel_id, fuel_label=r.fuel_label, hc_class=p[0], n_carbon=p[1], cn_open_ended=p[2],
                         mass_pct=float(r.mass_pct), metric_raw=r.Metric, provenance="real_fuel_gcxgc_bins"))
    comp = pd.DataFrame(rows)
    if len(comp):                                          # a page can repeat a bin: keep one value, they must agree
        comp = comp.groupby(["fuel_id", "fuel_label", "hc_class", "n_carbon", "cn_open_ended", "metric_raw", "provenance"],
                            as_index=False).mass_pct.first()
    return comp, sorted(unparsed)


def fuel_table(comp: pd.DataFrame, wide: Path = WIDE) -> pd.DataFrame:
    cls = comp.pivot_table(index=["fuel_id", "fuel_label"], columns="hc_class", values="mass_pct", aggfunc="sum", fill_value=0.0)
    f = cls.add_prefix("pct_").reset_index()
    f["leaf_sum_pct"] = cls.sum(axis=1).to_numpy()
    f["n_bins"] = comp.groupby(["fuel_id", "fuel_label"]).size().to_numpy()
    f["composition_complete"] = f.leaf_sum_pct.between(*CLOSURE)
    if wide.exists():
        p = pd.read_csv(wide, low_memory=False)
        key = "fuel_id" if "fuel_id" in p else p.columns[0]
        for tgt in TARGETS:
            cols = [c for c in p.columns if c == tgt or c.startswith(tgt)]
            has = p[cols].notna().any(axis=1) if cols else pd.Series(False, index=p.index)
            f[f"has_{tgt}"] = f.fuel_id.astype(str).map(dict(zip(p[key].astype(str), has))).fillna(False).astype(bool)
    return f


def main() -> None:
    comp, unparsed = load()
    f = fuel_table(comp)
    OUT.mkdir(parents=True, exist_ok=True)
    comp.to_csv(OUT / "ajftd_gcxgc_composition.csv", index=False)
    f.to_csv(OUT / "ajftd_gcxgc_fuels.csv", index=False)
    ok = f[f.composition_complete]
    log = dict(source=str(RAW.relative_to(ROOT)).replace("\\", "/"), fuels_with_bins=int(len(f)), fuels_complete=int(len(ok)),
               closure_window_pct=list(CLOSURE), bins=int(len(comp)), unparsed_metric_labels=unparsed,
               complete_with_property={t: int(ok[f"has_{t}"].sum()) for t in TARGETS if f"has_{t}" in ok},
               note="weight % as printed by AJFTD (AFRL GCxGC); bins are class x carbon number, isomers unresolved")
    (OUT / "ajftd_gcxgc_build_log.json").write_text(json.dumps(log, indent=1), encoding="utf-8")
    print(json.dumps(log, indent=1))


if __name__ == "__main__":
    main()
