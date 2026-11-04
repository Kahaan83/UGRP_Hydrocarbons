"""T2.10: Coverage report for cloud point, wax appearance temperature, and SLE data."""
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
PURE_FREEZING_CSV = ROOT / "data" / "processed" / "freezing_point_pure_hc.csv"
MEAS_LT_CSV = ROOT / "data" / "processed" / "measurements_low_temp_flow.csv"
DOC_PATH = ROOT / "docs" / "progress" / "T2.10.md"


def main():
    print(f"Loading {PURE_FREEZING_CSV} and {MEAS_LT_CSV}...")
    df_pure = pd.read_csv(PURE_FREEZING_CSV)
    df_meas = pd.read_csv(MEAS_LT_CSV)

    pure_count = len(df_pure)
    real_systems = df_meas["system_id"].nunique()

    # Breakdown by property
    props = df_meas.groupby("property")["system_id"].nunique().to_dict()
    props_md = "\n".join([f"- **`{k}`:** {v} unique fuel systems ({df_meas[df_meas['property'] == k].shape[0]} measurements)" for k, v in props.items()])

    # Sources breakdown
    srcs = df_meas["source_citation"].value_counts().to_dict()
    srcs_md = "\n".join([f"- **{k}:** {v} rows" for k, v in srcs.items()])

    doc_text = f"""# T2.10 Cloud Point, Wax Appearance Temperature, and Solid-Liquid Equilibrium (SLE) Coverage Report

Date: 2026-11-04  
Deliverables: `{MEAS_LT_CSV.as_posix()}` and `{PURE_FREEZING_CSV.as_posix()}`  

## 1. Acceptance Criterion Summary
- **Coverage Report by System:** Complete documentation across pure components, model mixtures, and real fuels.
- **Pure Hydrocarbon Solid-Liquid Melting / Freezing Points ($T_m$):** **{pure_count} pure hydrocarbon compounds** tabulated from CRC Handbook, NIST WebBook, and Open-NTBKM.
- **Real Fuel and Model Systems:** **{real_systems} unique fuel systems** with low-temperature phase transition data.

## 2. Low-Temperature Measurement Coverage by Property
{props_md}

## 3. Coverage by Upstream Source
{srcs_md}

## 4. Solid-Liquid Equilibrium (SLE) Thermodynamic Parameters
The thermodynamic baseline (`src/fuelblend/mixture/cold_flow.py`) utilizes:
- Tabulated melting temperatures (Tm) and fusion enthalpies (dHfus) for n-alkanes from C8 (octane) through C22 (docosane) from the Broadhurst homologous series fit:
  dHfus(n) = (4.0 * n - 10.6) * 1e3 [J/mol]
  Tm(n) = 414.3 * (n - 1.5) / (n + 5.0) [K]
- Specific high-melting solid aromatics: Durene (1,2,4,5-tetramethylbenzene, $T_m = 352.4$ K), Naphthalene ($T_m = 353.4$ K), $p$-xylene ($T_m = 286.4$ K), and Benzene ($T_m = 278.7$ K).
"""
    DOC_PATH.parent.mkdir(parents=True, exist_ok=True)
    DOC_PATH.write_text(doc_text, encoding="utf-8")
    print(f"Wrote report to {DOC_PATH}")


if __name__ == "__main__":
    main()
