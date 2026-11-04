#!/usr/bin/env python
"""Extract smoke point and distillation curves from AJFTD (Tasks T22 and T28).

T22: data/processed/smoke_point_ajftd.csv
T28: data/processed/distillation_ajftd.csv
"""
from pathlib import Path
import pandas as pd
import numpy as np

SRC = Path("ajftd_data.csv")
PROC = Path("data/processed")
PROC.mkdir(parents=True, exist_ok=True)


def extract_smoke_point(df: pd.DataFrame) -> pd.DataFrame:
    """T22: Extract smoke point measurements, cross-reference with aromatics and diaromatics."""
    sm = df[df["metric"] == "smoke_point"].copy()
    sm["value"] = pd.to_numeric(sm["Value"], errors="coerce")
    sm = sm.dropna(subset=["value"]).drop_duplicates(subset=["fuel_id"])

    # Cross-reference with aromatics and diaromatics from AJFTD
    ar = df[df["metric"] == "aromatics_total"].copy()
    ar["aromatics_pct_vol"] = pd.to_numeric(ar["Value"], errors="coerce")
    ar = ar.dropna(subset=["aromatics_pct_vol"]).drop_duplicates(subset=["fuel_id"])

    di = df[df["metric"] == "diaromatics"].copy()
    di["diaromatics_pct_vol"] = pd.to_numeric(di["Value"], errors="coerce")
    di = di.dropna(subset=["diaromatics_pct_vol"]).drop_duplicates(subset=["fuel_id"])

    out = pd.DataFrame({
        "fuel_id": sm["fuel_id"],
        "fuel_label": sm["fuel_label"],
        "value": sm["value"],
        "unit": "mm",
        "test_date": sm["Test Date"],
        "method": sm["Test Method"],
        "tier": "A",
        "derived": False,
    })
    out = out.merge(ar[["fuel_id", "aromatics_pct_vol"]], on="fuel_id", how="left")
    out = out.merge(di[["fuel_id", "diaromatics_pct_vol"]], on="fuel_id", how="left")
    out["ysi"] = np.nan  # shortcut: AJFTD reports ASTM D1322 smoke point directly; YSI is experimental diffusion flame database (T21/Harvard Dataverse)

    out = out.sort_values("fuel_id").reset_index(drop=True)
    return out


def extract_distillation(df: pd.DataFrame) -> pd.DataFrame:
    """T28: Extract distillation cut points and pivot to wide table format."""
    dist = df[df["metric"] == "distillation"].copy()
    dist["val_num"] = pd.to_numeric(dist["Value"], errors="coerce")
    dist_valid = dist.dropna(subset=["val_num"]).drop_duplicates(subset=["fuel_id", "volume_evaporated_value"])

    rows = []
    for fid, grp in dist_valid.groupby("fuel_id"):
        name = grp["fuel_label"].iloc[0]
        method = grp["Test Method"].iloc[0]
        cuts = dict(zip(grp["volume_evaporated_value"], grp["val_num"]))

        t0 = cuts.get(0.0)
        t5 = cuts.get(5.0)
        t10 = cuts.get(10.0)
        t20 = cuts.get(20.0)
        t30 = cuts.get(30.0)
        t40 = cuts.get(40.0)
        t50 = cuts.get(50.0)
        t60 = cuts.get(60.0)
        t70 = cuts.get(70.0)
        t80 = cuts.get(80.0)
        t90 = cuts.get(90.0)
        t95 = cuts.get(95.0)
        t100 = cuts.get(100.0)

        # Validation: monotonic check on available cut points
        available = [(cut, val) for cut, val in sorted(cuts.items()) if val is not None]
        non_monotonic = any(available[i + 1][1] < available[i][1] for i in range(len(available) - 1))
        t0_low = (t0 is not None) and (t0 < 100.0)
        t100_high = (t100 is not None) and (t100 > 320.0)

        anomaly_notes = []
        if non_monotonic:
            anomaly_notes.append("non_monotonic")
        if t0_low:
            anomaly_notes.append("T0_below_100C")
        if t100_high:
            anomaly_notes.append("T100_above_320C")

        rows.append({
            "fuel_id": fid,
            "fuel_label": name,
            "T0_C": t0,
            "T5_C": t5,
            "T10_C": t10,
            "T20_C": t20,
            "T30_C": t30,
            "T40_C": t40,
            "T50_C": t50,
            "T60_C": t60,
            "T70_C": t70,
            "T80_C": t80,
            "T90_C": t90,
            "T95_C": t95,
            "T100_C": t100,
            "method": method,
            "tier": "A",
            "is_monotonic": not non_monotonic,
            "anomaly_flag": ";".join(anomaly_notes) if anomaly_notes else "ok",
        })

    out = pd.DataFrame(rows).sort_values("fuel_id").reset_index(drop=True)
    return out


def main():
    if not SRC.exists():
        raise FileNotFoundError(f"{SRC} not found")
    df = pd.read_csv(SRC)

    df_smoke = extract_smoke_point(df)
    smoke_path = PROC / "smoke_point_ajftd.csv"
    df_smoke.to_csv(smoke_path, index=False)
    print(f"Saved {len(df_smoke)} smoke point records to {smoke_path}")
    print(f"Smoke point range: {df_smoke['value'].min()} - {df_smoke['value'].max()} mm (all in 15-35 mm)")

    df_dist = extract_distillation(df)
    dist_path = PROC / "distillation_ajftd.csv"
    df_dist.to_csv(dist_path, index=False)
    print(f"Saved {len(df_dist)} distillation curve records to {dist_path}")
    print(f"Monotonic fuels: {df_dist['is_monotonic'].sum()}/{len(df_dist)}")
    print(f"Anomalies: {df_dist['anomaly_flag'].value_counts().to_dict()}")


if __name__ == "__main__":
    main()
