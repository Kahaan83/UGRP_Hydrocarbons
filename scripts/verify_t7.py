#!/usr/bin/env python
"""Check T7 file passes unified schema."""
import pandas as pd

df = pd.read_csv("data/processed/t7_leenhouts_flash_mix.csv")
print("T7 rows:", len(df))
print("columns:", list(df.columns))
for c in ["source", "property", "value", "T_K", "P_kPa", "kind", "n_comp",
          "is_hc", "system_id", "comp_inchikeys", "x_mole", "w_mass",
          "fuel_id", "ref", "raw_value", "raw_unit"]:
    print(f"  {c}: {'OK' if c in df.columns else 'MISSING'}")

if "property" in df.columns:
    print("property values:", df["property"].unique()[:10].tolist())
if "kind" in df.columns:
    print("kind values:", df["kind"].unique()[:10].tolist())
if "system_id" in df.columns:
    print("system_id count:", df["system_id"].nunique())
    print("sample system_ids:", df["system_id"].head(3).tolist())
if "T_K" in df.columns:
    print("T_K stats: min=", df["T_K"].min(), "max=", df["T_K"].max(), "nans=", df["T_K"].isna().sum())
if "value" in df.columns:
    bad = df["value"].isna().sum() + (~pd.isfinite(df["value"])).sum()
    print("value bad:", bad)
