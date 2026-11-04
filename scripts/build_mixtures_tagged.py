"""Tag all 510 existing mixture systems with fuel_family and real_fuel (T1.4)."""
import json
from pathlib import Path
import pandas as pd

from fuelblend.features.fuel_family import classify_mixture

ROOT = Path(__file__).resolve().parents[1]
MODEL_TABLE_PATH = ROOT / "data" / "processed" / "model_table.csv"
COMPOUNDS_PATH = ROOT / "data" / "unified" / "compounds.csv"
OUT_PATH = ROOT / "data" / "processed" / "mixtures_tagged.csv"


def main():
    print(f"Loading {MODEL_TABLE_PATH}...")
    df = pd.read_csv(MODEL_TABLE_PATH, low_memory=False)
    comp_df = pd.read_csv(COMPOUNDS_PATH).set_index("inchikey")

    # Select mixture systems (binary and multi)
    mix_df = df[df["kind"].isin(["binary", "multi"])].copy()
    unique_systems = mix_df["system_id"].unique()
    print(f"Found {len(unique_systems)} unique mixture systems in model_table.")
    assert len(unique_systems) == 510, f"Expected exactly 510 mixture systems, found {len(unique_systems)}"

    tagged_rows = []
    for sys_id in unique_systems:
        sys_rows = mix_df[mix_df["system_id"] == sys_id]
        first = sys_rows.iloc[0]
        n_comp = int(first["n_comp"])
        try:
            keys = json.loads(first["comp_inchikeys"])
        except Exception:
            keys = []

        c_nums = [comp_df.loc[k, "n_carbon"] for k in keys if k in comp_df.index and pd.notna(comp_df.loc[k, "n_carbon"])]
        c_nums = [int(c) for c in c_nums]

        props = sorted(list(sys_rows["property"].unique()))
        sources = "+".join(sorted(list(set("+".join(sys_rows["sources"].dropna().astype(str)).split("+")))))

        res = classify_mixture(
            system_id=sys_id,
            n_comp=n_comp,
            comp_inchikeys=keys,
            properties=props,
            comp_c_nums=c_nums,
            source=sources,
        )

        min_c = min(c_nums) if c_nums else None
        max_c = max(c_nums) if c_nums else None
        avg_c = round(sum(c_nums) / len(c_nums), 2) if c_nums else None

        tagged_rows.append({
            "system_id": sys_id,
            "n_comp": n_comp,
            "fuel_family": res["fuel_family"],
            "real_fuel": res["real_fuel"],
            "classification_reason": res["reason"],
            "properties": ";".join(props),
            "sources": sources,
            "min_carbon": min_c,
            "max_carbon": max_c,
            "avg_carbon": avg_c,
            "comp_inchikeys": first["comp_inchikeys"],
        })

    out_df = pd.DataFrame(tagged_rows)
    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    out_df.to_csv(OUT_PATH, index=False)
    print(f"Wrote {len(out_df)} tagged mixtures to {OUT_PATH}")
    print("\nBreakdown by fuel_family:")
    print(out_df["fuel_family"].value_counts().to_string())
    print("\nBreakdown by real_fuel:")
    print(out_df["real_fuel"].value_counts().to_string())


if __name__ == "__main__":
    main()
