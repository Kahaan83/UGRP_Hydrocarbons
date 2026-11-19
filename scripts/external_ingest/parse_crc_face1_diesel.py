#!/usr/bin/env python3
"""parse_crc_face1_diesel.py - Parse properties and 2-D GC composition for 9 FACE diesel fuels.

Source: CRC Report No. FACE-1 "Chemical and Physical Properties of the Fuels for Advanced
Combustion Engines (FACE) Research Diesel Fuels" (Coordinating Research Council, 2014).
File: data/raw/diesel/CRC_FACE-1.pdf

Outputs:
  - data/interim/external_staging/crc_face1_diesel_properties.csv
  - data/interim/external_staging/crc_face1_diesel_composition_table52.csv

NOTE: Staged output ONLY. DOES NOT alter active canonical datasets or production pipeline.
"""

from pathlib import Path
import re
import pandas as pd
import pdfplumber

RAW_PDF = Path("data/raw/diesel/CRC_FACE-1.pdf")
STAGING_DIR = Path("data/interim/external_staging")

FUELS = ["FD1A", "FD2A", "FD3A", "FD4A", "FD5A", "FD6A", "FD7A", "FD8A", "FD9A"]


def parse_table52_composition(pdf: pdfplumber.PDF) -> pd.DataFrame:
    """Parse Table 5.2: Quantitative Group Type 2-D GC-FID Separation (Page 80, index 79)."""
    p80 = pdf.pages[79]
    txt = p80.extract_text() or ""
    rows = []
    # Lines like: FD1A 6.78 59.06 52.28 12.17 8.08 19.97 28.05 99.28
    # Notice FD5A is labelled 'Face5' in Table 5.2
    for line in txt.splitlines():
        line = line.strip()
        m = re.match(r"^(FD[1-9]A|Face5)\s+([\d\.]+)\s+([\d\.]+)\s+([\d\.]+)\s+([\d\.]+)\s+([\d\.]+)\s+([\d\.]+)\s+([\d\.]+)\s+([\d\.]+)$", line)
        if m:
            fuel_label, n_paraffins, paraffins_tot, isoparaffins, cycloalkanes, known_aro, unk_aro, tot_aro, tot = m.groups()
            norm_fuel = "FD5A" if fuel_label == "Face5" else fuel_label
            rows.append({
                "fuel_id": norm_fuel,
                "n_paraffins_pct": float(n_paraffins),
                "isoparaffins_pct": float(isoparaffins),
                "cycloalkanes_pct": float(cycloalkanes),
                "total_aromatics_pct": float(tot_aro),
                "total_recovered_pct": float(tot),
                "method": "2-D GC-FID group type separation",
                "source_table": "Table 5.2",
                "provenance": "CRC FACE-1 Report (2014)",
            })

    df = pd.DataFrame(rows).sort_values("fuel_id").reset_index(drop=True)
    return df


def parse_table31_dcn(pdf: pdfplumber.PDF) -> dict[str, float]:
    """Parse Table 3.1: ASTM D6890 Derived Cetane Number (Page 48, index 47)."""
    p48 = pdf.pages[47]
    txt = p48.extract_text() or ""
    dcn_map = {}
    for line in txt.splitlines():
        m = re.match(r"^(FD[1-9]A)\s+[\d\.]+\s+[\d\.]+\s+Eq\.\s+\([34]\.[12]\)\s+([\d\.]+)", line.strip())
        if m:
            fuel_id, dcn_str = m.groups()
            dcn_map[fuel_id] = float(dcn_str)
    return dcn_map


def parse_table22_properties(pdf: pdfplumber.PDF) -> pd.DataFrame:
    """Parse Table 2.2: Summary Data for ASTM Standard Tests (Pages 37-43)."""
    # Specific values compiled directly from Table 2.2 consensus/benchmark lab results:
    # 1. Cetane Number ASTM D613 (Engine Method, ORNL/SwRI):
    cn_d613 = {"FD1A": 30.7, "FD2A": 28.7, "FD3A": 30.7, "FD4A": 28.5, "FD5A": 55.0, "FD6A": 54.1, "FD7A": 45.9, "FD8A": 49.0, "FD9A": 43.5}
    # 2. Kinematic Viscosity at 40 °C in cSt (ASTM D445, ORNL/SwRI):
    visc_40 = {"FD1A": 1.741, "FD2A": 1.960, "FD3A": 1.320, "FD4A": 1.785, "FD5A": 1.795, "FD6A": 3.097, "FD7A": 1.759, "FD8A": 3.278, "FD9A": 2.107}
    # 3. Cloud Point in °C (ASTM D2500, ORNL/SwRI):
    cloud = {"FD1A": -27.0, "FD2A": -51.0, "FD3A": -35.0, "FD4A": -50.0, "FD5A": -29.0, "FD6A": -25.0, "FD7A": -21.0, "FD8A": -20.0, "FD9A": -30.0}
    # 4. Flash Point closed cup in °F (ASTM D93, ORNL/SwRI):
    flash_f = {"FD1A": 135.0, "FD2A": 128.0, "FD3A": 141.0, "FD4A": 139.0, "FD5A": 131.0, "FD6A": 166.0, "FD7A": 149.0, "FD8A": 158.0, "FD9A": 140.0}
    # 5. Specific Gravity at 15.56 °C (ASTM D4052, ORNL/SwRI):
    sg = {"FD1A": 0.8084, "FD2A": 0.8037, "FD3A": 0.8400, "FD4A": 0.8359, "FD5A": 0.8086, "FD6A": 0.8412, "FD7A": 0.8375, "FD8A": 0.8694, "FD9A": 0.8464}
    # 6. Net Heat of Combustion in BTU/lb (ASTM D240):
    hoc_btu = {"FD1A": 18402.0, "FD2A": 18553.1, "FD3A": 18120.0, "FD4A": 18269.4, "FD5A": 18442.6, "FD6A": 18399.4, "FD7A": 18211.3, "FD8A": 18140.9, "FD9A": 18256.5}

    dcn_map = parse_table31_dcn(pdf)

    rows = []
    for f in FUELS:
        # Cetane D613
        rows.append({
            "fuel_id": f,
            "property": "cetane_number_D613",
            "value": cn_d613[f],
            "unit": "-",
            "T_K": float("nan"),
            "method": "ASTM D613",
            "lab": "ORNL/SwRI",
            "source_table": "Table 2.2",
        })
        # DCN D6890
        if f in dcn_map:
            rows.append({
                "fuel_id": f,
                "property": "dcn",
                "value": dcn_map[f],
                "unit": "-",
                "T_K": float("nan"),
                "method": "ASTM D6890 (IQT)",
                "lab": "NREL/CRC",
                "source_table": "Table 3.1",
            })
        # Viscosity
        rows.append({
            "fuel_id": f,
            "property": "viscosity_kin_cSt",
            "value": visc_40[f],
            "unit": "cSt",
            "T_K": 313.15,
            "method": "ASTM D445",
            "lab": "ORNL/SwRI",
            "source_table": "Table 2.2",
        })
        # Cloud point
        rows.append({
            "fuel_id": f,
            "property": "cloud_point_C",
            "value": cloud[f],
            "unit": "C",
            "T_K": float("nan"),
            "method": "ASTM D2500",
            "lab": "ORNL/SwRI",
            "source_table": "Table 2.2",
        })
        # Flash point: convert F to C
        fp_c = round((flash_f[f] - 32.0) * 5.0 / 9.0, 2)
        rows.append({
            "fuel_id": f,
            "property": "flash_point_C",
            "value": fp_c,
            "unit": "C",
            "T_K": float("nan"),
            "method": "ASTM D93 (closed cup)",
            "lab": "ORNL/SwRI",
            "source_table": "Table 2.2",
        })
        # Density: SG * 999.016 kg/m3
        dens = round(sg[f] * 999.016, 2)
        rows.append({
            "fuel_id": f,
            "property": "density_kg_m3",
            "value": dens,
            "unit": "kg/m3",
            "T_K": 288.71,
            "method": "ASTM D4052",
            "lab": "ORNL/SwRI",
            "source_table": "Table 2.2",
        })
        # Heat of combustion: BTU/lb to MJ/kg (1 BTU/lb = 0.002326 MJ/kg)
        hoc_mj = round(hoc_btu[f] * 0.002326, 3)
        rows.append({
            "fuel_id": f,
            "property": "heat_of_combustion_net_MJ_kg",
            "value": hoc_mj,
            "unit": "MJ/kg",
            "T_K": 298.15,
            "method": "ASTM D240",
            "lab": "ORNL/SwRI",
            "source_table": "Table 2.2",
        })

    df = pd.DataFrame(rows)
    return df


def main():
    if not RAW_PDF.exists():
        print(f"ERROR: Raw PDF not found at {RAW_PDF}")
        return

    STAGING_DIR.mkdir(parents=True, exist_ok=True)
    with pdfplumber.open(RAW_PDF) as pdf:
        df_comp = parse_table52_composition(pdf)
        out_comp = STAGING_DIR / "crc_face1_diesel_composition_table52.csv"
        df_comp.to_csv(out_comp, index=False)
        print(f"Saved Table 5.2 composition: {len(df_comp)} fuels -> {out_comp}")

        df_props = parse_table22_properties(pdf)
        out_props = STAGING_DIR / "crc_face1_diesel_properties.csv"
        df_props.to_csv(out_props, index=False)
        print(f"Saved properties: {len(df_props)} rows -> {out_props}")

    print("CRC FACE-1 diesel data parsed and staged successfully!")


if __name__ == "__main__":
    main()
