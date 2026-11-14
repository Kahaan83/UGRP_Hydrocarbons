"""AJFTD smoke point + naphthalenes per fuel (v5 delta_smoke_cert 1.1).

Reads the raw scrape CSV (outside the repo, never committed); rules are declared in
results/v5/BASE/delta_smoke_cert/README.md.  One row per fuel_id; smoke = median of uncensored rows,
censored-only fuels keep NaN + qualifier; naphthalenes only from vol% D1840 rows (class total, not the
single compound 'naphthalene', no mass% conversion).
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from fuelblend.data.loaders_ajftd import parse_value, source_group

REPO = Path(__file__).resolve().parents[3]
RAW_DEFAULT = REPO.parent / "extras" / "fuelblend_data" / "data" / "raw" / "ajftd" / "ajftd_all_metrics.csv"
AROM_DEFAULT = REPO / "data" / "real_fuels" / "v5" / "ajftd_fuels_long.csv"
OUT_COLS = ["fuel_id", "source_group", "smoke_point_mm", "smoke_qualifier", "naphthalenes_vol_pct",
            "aromatics_vol_pct", "method_smoke"]


def _norm_method(m) -> str:
    return str(m).replace(" ", "") if pd.notna(m) else ""


def load_ajftd_smoke_naph(raw_path=RAW_DEFAULT, aromatics_path=AROM_DEFAULT) -> pd.DataFrame:
    raw = pd.read_csv(raw_path, encoding="utf-8-sig", dtype=str, keep_default_na=False)
    raw["fuel_id"] = raw["fuel_id"].astype(int)
    labels = raw.drop_duplicates("fuel_id").set_index("fuel_id")["fuel_label"]

    s = raw[raw["metric"] == "smoke_point"].copy()
    parsed = s["Value"].map(parse_value)
    s["num"] = [p[0] for p in parsed]
    s["q"] = [p[1] for p in parsed]
    s["method"] = s["Test Method"].map(_norm_method)
    keep = s[s["num"].notna() | (s["q"] != "")]  # drops No Test / Not Found / blank
    rows = {}
    for fid, g in keep.groupby("fuel_id"):
        unc = g[(g["q"] == "") & g["num"].notna()]
        if len(unc):
            val, q = float(unc["num"].median()), ""
        else:
            val, q = np.nan, g["q"].iloc[0]
        meth = g["method"].replace("", np.nan).dropna()
        rows[fid] = {"smoke_point_mm": val, "smoke_qualifier": q,
                     "method_smoke": meth.mode().sort_values().iloc[0] if len(meth) else ""}
    sm = pd.DataFrame.from_dict(rows, orient="index")

    n = raw[(raw["metric"] == "naphthalenes") & (raw["Unit"].str.strip() == "vol%")
            & (raw["Test Method"].map(_norm_method) == "D1840")].copy()
    n["num"] = [parse_value(v)[0] for v in n["Value"]]
    nap = n.dropna(subset=["num"]).groupby("fuel_id")["num"].median().rename("naphthalenes_vol_pct")

    ar = pd.read_csv(aromatics_path, usecols=["fuel_id", "property", "value"])
    ar = ar[ar["property"] == "aromatics_vol_pct"].groupby("fuel_id")["value"].median().rename("aromatics_vol_pct")

    out = sm.join(nap, how="outer").join(ar, how="outer")
    out.index.name = "fuel_id"
    out = out.reset_index()
    out["source_group"] = out["fuel_id"].map(labels).map(source_group)
    out["smoke_qualifier"] = out["smoke_qualifier"].fillna("")
    out["method_smoke"] = out["method_smoke"].fillna("")
    return out[OUT_COLS].sort_values("fuel_id").reset_index(drop=True)
