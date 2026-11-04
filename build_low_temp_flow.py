"""Extract and consolidate real-fuel and blend low-temperature flow properties.

Combines:
- CRC aviation fuels (CRC Report 647, AV-11-09)
- AJFTD alternative & conventional jet fuels
- NASA NTRS 19860018789 (Partially Frozen Hydrocarbon Fuels, 8 fuels)
- Kirgina et al. 2021/2019 (Straight-run diesel fuels, pure and cold-flow improved)

Outputs:
  data/processed/low_temp_flow_properties_fuels.csv
"""

import os
import pandas as pd
import numpy as np

def build_low_temp_flow_properties():
    records = []

    # 1. CRC fuels (crc_fuels.csv)
    crc_path = os.path.join("data", "processed", "crc_fuels.csv")
    if os.path.exists(crc_path):
        crc = pd.read_csv(crc_path, dtype={"fuel_id": str})
        crc_lt = crc[crc["property"].isin(["pour_point_C", "freezing_point_C", "cloud_point_C"])].copy()
        for _, row in crc_lt.iterrows():
            raw_val = str(row["value_raw"]).strip()
            qualifier = "<" if "<" in raw_val else (">" if ">" in raw_val else np.nan)
            records.append({
                "fuel_id": str(row["fuel_id"]),
                "fuel_name": str(row["fuel_name"]),
                "source": str(row["source"]),
                "property": str(row["property"]),
                "value": float(row["value"]),
                "unit": "C",
                "qualifier": qualifier,
                "method": str(row["method"]) if pd.notna(row["method"]) else np.nan,
                "wax_content_pct": np.nan,
                "n_paraffin_content_pct": np.nan,
                "tier": str(row["tier"]) if pd.notna(row["tier"]) else "real_fuel",
                "derived": "False",
                "note": f"CRC aviation fuel; raw: {raw_val}"
            })

    # 2. AJFTD fuels (ajftd_fuels.csv)
    ajftd_path = os.path.join("data", "processed", "ajftd_fuels.csv")
    if os.path.exists(ajftd_path):
        ajftd = pd.read_csv(ajftd_path)
        ajftd_lt = ajftd[ajftd["property"].isin(["pour_point_C", "freezing_point_C", "cloud_point_C"])].copy()
        for _, row in ajftd_lt.iterrows():
            qualifier = row["qualifier"] if pd.notna(row["qualifier"]) else np.nan
            records.append({
                "fuel_id": f"AJFTD_{row['fuel_id']}",
                "fuel_name": str(row["fuel_label"]),
                "source": "ajftd",
                "property": str(row["property"]),
                "value": float(row["value"]),
                "unit": "C",
                "qualifier": qualifier,
                "method": str(row["method"]) if pd.notna(row["method"]) else np.nan,
                "wax_content_pct": np.nan,
                "n_paraffin_content_pct": np.nan,
                "tier": "real_fuel",
                "derived": str(row["derived"]) if pd.notna(row["derived"]) else "False",
                "note": "AJFTD database entry"
            })

    # 3. NASA NTRS 19860018789 (Cotterell et al. 1986)
    # 8 fuels with n-alkane wt%, ASTM mp, ASTM fp, Equil mp, DSC mp, Equil pour pt + wt% solid at pour pt
    nasa_fuels = [
        {
            "id": "NASA_LFP-3", "name": "LFP-3 (paraffinic diesel)", "n_paraffin": 18.4,
            "props": [
                ("freezing_point_C", -19.3, "ASTM D2386", 4.7, "ASTM freezing point"),
                ("freezing_point_C", -16.5, "ASTM D2386 (melting)", 4.7, "ASTM melting point"),
                ("freezing_point_C", -17.6, "Static equilibrium melting", 4.7, "Equilibrium melting point"),
                ("freezing_point_C", -16.8, "DSC (5 C/min)", 4.7, "DSC melting point"),
                ("pour_point_C", -27.0, "Static equilibrium pour (1 C/min)", 4.7, "Equilibrium pour point (wt% solid=4.7)")
            ]
        },
        {
            "id": "NASA_LFP-5", "name": "LFP-5 (paraffinic distillate)", "n_paraffin": 13.9,
            "props": [
                ("freezing_point_C", -30.5, "ASTM D2386", 3.4, "ASTM freezing point"),
                ("freezing_point_C", -28.2, "ASTM D2386 (melting)", 3.4, "ASTM melting point"),
                ("freezing_point_C", -28.7, "Static equilibrium melting", 3.4, "Equilibrium melting point"),
                ("freezing_point_C", -28.9, "DSC (5 C/min)", 3.4, "DSC melting point"),
                ("pour_point_C", -36.0, "Static equilibrium pour (1 C/min)", 3.4, "Equilibrium pour point (wt% solid=3.4)")
            ]
        },
        {
            "id": "NASA_Fuel_7", "name": "Fuel #7 (LFP-5 + 0.1% pour depressant)", "n_paraffin": 13.9,
            "props": [
                ("freezing_point_C", -27.4, "Static equilibrium melting", 7.7, "Equilibrium melting point"),
                ("freezing_point_C", -28.7, "DSC (5 C/min)", 7.7, "DSC melting point"),
                ("pour_point_C", -52.0, "Static equilibrium pour (1 C/min)", 7.7, "Equilibrium pour point (wt% solid=7.7, 0.1% additive)")
            ]
        },
        {
            "id": "NASA_LFPA-3", "name": "LFPA-3 (aviation turbine blend)", "n_paraffin": 22.2,
            "props": [
                ("freezing_point_C", -49.7, "ASTM D2386", 3.5, "ASTM freezing point"),
                ("freezing_point_C", -48.6, "ASTM D2386 (melting)", 3.5, "ASTM melting point"),
                ("freezing_point_C", -44.4, "DSC (5 C/min)", 3.5, "DSC melting point"),
                ("pour_point_C", -47.0, "Static equilibrium pour (1 C/min)", 3.5, "Equilibrium pour point (wt% solid=3.5)")
            ]
        },
        {
            "id": "NASA_Shale_II_JP-5", "name": "Shale II JP-5 (Paraho shale kerosene)", "n_paraffin": 36.0,
            "props": [
                ("freezing_point_C", -32.2, "ASTM D2386", 3.0, "ASTM freezing point"),
                ("freezing_point_C", -27.6, "ASTM D2386 (melting)", 3.0, "ASTM melting point"),
                ("freezing_point_C", -49.8, "Static equilibrium melting", 3.0, "Equilibrium melting point"),
                ("freezing_point_C", -50.0, "DSC (5 C/min)", 3.0, "DSC melting point"),
                ("pour_point_C", -53.0, "Static equilibrium pour (1 C/min)", 3.0, "Equilibrium pour point (wt% solid=3.0)")
            ]
        },
        {
            "id": "NASA_NAPC-2", "name": "NAPC-2 (modified JP-5)", "n_paraffin": 17.3,
            "props": [
                ("freezing_point_C", -24.0, "ASTM D2386", 3.2, "ASTM freezing point"),
                ("freezing_point_C", -26.1, "Static equilibrium melting", 3.2, "Equilibrium melting point"),
                ("freezing_point_C", -25.9, "DSC (5 C/min)", 3.2, "DSC melting point"),
                ("pour_point_C", -30.0, "Static equilibrium pour (1 C/min)", 3.2, "Equilibrium pour point (wt% solid=3.2)")
            ]
        },
        {
            "id": "NASA_NAPC-3", "name": "NAPC-3 (modified JP-5)", "n_paraffin": 19.2,
            "props": [
                ("freezing_point_C", -34.0, "ASTM D2386", 2.3, "ASTM freezing point"),
                ("freezing_point_C", -35.7, "Static equilibrium melting", 2.3, "Equilibrium melting point"),
                ("freezing_point_C", -36.3, "DSC (5 C/min)", 2.3, "DSC melting point"),
                ("pour_point_C", -41.0, "Static equilibrium pour (1 C/min)", 2.3, "Equilibrium pour point (wt% solid=2.3)")
            ]
        },
        {
            "id": "NASA_NAPC-5", "name": "NAPC-5 (low aromatic JP-5)", "n_paraffin": 16.0,
            "props": [
                ("freezing_point_C", -50.0, "ASTM D2386", 2.8, "ASTM freezing point"),
                ("freezing_point_C", -50.0, "Static equilibrium melting", 2.8, "Equilibrium melting point"),
                ("freezing_point_C", -50.0, "DSC (5 C/min)", 2.8, "DSC melting point"),
                ("pour_point_C", -57.0, "Static equilibrium pour (1 C/min)", 2.8, "Equilibrium pour point (wt% solid=2.8)")
            ]
        }
    ]

    for f in nasa_fuels:
        for prop, val, method, wax, note in f["props"]:
            records.append({
                "fuel_id": f["id"],
                "fuel_name": f["name"],
                "source": "NASA NTRS 19860018789",
                "property": prop,
                "value": float(val),
                "unit": "C",
                "qualifier": np.nan,
                "method": method,
                "wax_content_pct": wax,
                "n_paraffin_content_pct": f["n_paraffin"],
                "tier": "real_fuel",
                "derived": "False",
                "note": note
            })

    # 4. Kirgina et al. 2021 / 2019 (Straight-run diesel fuels)
    # From Kirgina et al., Oil & Gas Sci. Technol. 76 (2021) 40 / Pet. Coal 61(1) 2019
    # 5 straight-run diesel fuel samples with n-paraffin content, Cloud Point (Tcp), CFPP, Pour Point (Tpp)
    # pure and with cold flow improver Additive C
    kirgina_samples = [
        {"id": "Kirgina_DF1", "name": "Straight-run diesel sample 1", "n_paraffin": 20.77, "cp": -15.0, "cfpp": -21.0, "pp": -35.0},
        {"id": "Kirgina_DF2", "name": "Straight-run diesel sample 2", "n_paraffin": 19.08, "cp": -12.0, "cfpp": -22.0, "pp": -32.0},
        {"id": "Kirgina_DF3", "name": "Straight-run diesel sample 3", "n_paraffin": 19.61, "cp": -3.0, "cfpp": -4.0, "pp": -18.0},
        {"id": "Kirgina_DF4", "name": "Straight-run diesel sample 4", "n_paraffin": 22.05, "cp": -12.0, "cfpp": -24.0, "pp": -45.0},
        {"id": "Kirgina_DF5", "name": "Straight-run diesel sample 5", "n_paraffin": 19.88, "cp": -13.0, "cfpp": -17.0, "pp": -29.0},
        # Additive C blends (1C1 - 5C1)
        {"id": "Kirgina_DF1_AddC", "name": "Straight-run diesel sample 1 + Additive C", "n_paraffin": 20.77, "cp": -16.0, "cfpp": -32.0, "pp": -35.0},
        {"id": "Kirgina_DF2_AddC", "name": "Straight-run diesel sample 2 + Additive C", "n_paraffin": 19.08, "cp": -17.0, "cfpp": -30.0, "pp": -39.0},
        {"id": "Kirgina_DF3_AddC", "name": "Straight-run diesel sample 3 + Additive C", "n_paraffin": 19.61, "cp": -4.0, "cfpp": -17.0, "pp": -28.0},
        {"id": "Kirgina_DF4_AddC", "name": "Straight-run diesel sample 4 + Additive C", "n_paraffin": 22.05, "cp": -12.0, "cfpp": -26.0, "pp": -50.0},
        {"id": "Kirgina_DF5_AddC", "name": "Straight-run diesel sample 5 + Additive C", "n_paraffin": 19.88, "cp": -15.0, "cfpp": -28.0, "pp": -39.0},
    ]

    for s in kirgina_samples:
        is_add = "AddC" in s["id"]
        note_suffix = " (additized with cold-flow improver C)" if is_add else " (pure straight-run diesel)"
        # Cloud point
        records.append({
            "fuel_id": s["id"],
            "fuel_name": s["name"],
            "source": "Kirgina et al. 2019/2021",
            "property": "cloud_point_C",
            "value": s["cp"],
            "unit": "C",
            "qualifier": np.nan,
            "method": "GOST 5066 / ASTM D2500",
            "wax_content_pct": np.nan,
            "n_paraffin_content_pct": s["n_paraffin"],
            "tier": "real_fuel",
            "derived": "False",
            "note": f"Cloud point (Tcp){note_suffix}"
        })
        # CFPP
        records.append({
            "fuel_id": s["id"],
            "fuel_name": s["name"],
            "source": "Kirgina et al. 2019/2021",
            "property": "cfpp_C",
            "value": s["cfpp"],
            "unit": "C",
            "qualifier": np.nan,
            "method": "EN 116 / ASTM D6371",
            "wax_content_pct": np.nan,
            "n_paraffin_content_pct": s["n_paraffin"],
            "tier": "real_fuel",
            "derived": "False",
            "note": f"Cold filter plugging point (CFPP){note_suffix}"
        })
        # Pour point
        records.append({
            "fuel_id": s["id"],
            "fuel_name": s["name"],
            "source": "Kirgina et al. 2019/2021",
            "property": "pour_point_C",
            "value": s["pp"],
            "unit": "C",
            "qualifier": np.nan,
            "method": "GOST 20287 / ASTM D97",
            "wax_content_pct": np.nan,
            "n_paraffin_content_pct": s["n_paraffin"],
            "tier": "real_fuel",
            "derived": "False",
            "note": f"Pour point (Tpp){note_suffix}"
        })

    # Assemble DataFrame
    cols = [
        "fuel_id", "fuel_name", "source", "property", "value",
        "unit", "qualifier", "method", "wax_content_pct",
        "n_paraffin_content_pct", "tier", "derived", "note"
    ]
    df = pd.DataFrame(records, columns=cols)

    out_path = os.path.join("data", "processed", "low_temp_flow_properties_fuels.csv")
    df.to_csv(out_path, index=False)
    print(f"Saved {len(df)} records to {out_path}")
    print("Breakdown by property:")
    print(df["property"].value_counts())
    print("\nBreakdown by source:")
    print(df["source"].value_counts())
    return df

if __name__ == "__main__":
    build_low_temp_flow_properties()
