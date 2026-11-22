"""Data inventory for v9 acquisition (counts per target from committed files).

Targets:
1. Real-fuel freezing: fuels with measured freezing AND class x carbon-number composition.
2. Mixture systems: systems per property and n_comp in measurements_dedup.csv, split by T < 263 K.
3. Pure constants: has_antoine/has_lfl/has_viscosity/has_density for the 7 candidates in candidate_gaps.csv.
4. Cold viscosity: real fuels with viscosity at 253.15 K and 233.15 K.
"""
from __future__ import annotations

import argparse
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]


def count_real_fuel_freezing(real_fuels_dir: Path, freezing_lit_dir: Path) -> dict[str, int]:
    # 1. AJFTD GCxGC
    gc_fuels_file = real_fuels_dir / "ajftd_gcxgc_fuels.csv"
    if gc_fuels_file.exists():
        gc_fuels = pd.read_csv(gc_fuels_file)
        ajftd_count = int(gc_fuels["has_freezing_point_C"].sum())
    else:
        ajftd_count = 0

    # 2. CRC AV-33-22
    crc_props_file = real_fuels_dir / "crc_av3322_properties.csv"
    crc_comp_file = real_fuels_dir / "crc_av3322_composition.csv"
    if crc_props_file.exists() and crc_comp_file.exists():
        crc_props = pd.read_csv(crc_props_file)
        crc_comp = pd.read_csv(crc_comp_file)
        fids_with_comp = set(crc_comp["fuel_id"].unique())
        crc_fz = crc_props[
            (crc_props["property"] == "freezing_point_C")
            & (crc_props["fuel_id"].isin(fids_with_comp))
        ]
        crc_count = int(crc_fz["fuel_id"].nunique())
    else:
        crc_count = 0

    # 3. freezing_lit (Vozka 2018)
    vozka_props_file = freezing_lit_dir / "vozka2018_properties.csv"
    if vozka_props_file.exists():
        vozka_props = pd.read_csv(vozka_props_file)
        vozka_fz = vozka_props[vozka_props["property"] == "freezing_point_C"]
        vozka_count = int(vozka_fz["sample_id"].nunique())
    else:
        vozka_count = 0

    return {
        "ajftd_gcxgc": ajftd_count,
        "crc_av3322": crc_count,
        "vozka2018": vozka_count,
        "total_freezing_fuels": ajftd_count + crc_count + vozka_count,
    }


def count_mixture_systems(dedup_path: Path) -> list[dict]:
    df = pd.read_csv(dedup_path, low_memory=False)
    mix = df[df["n_comp"] > 1]

    records = []
    # Total mixture systems across all properties
    total_mix_systems = int(mix["system_id"].nunique())
    records.append({
        "target": "target_2_mixtures",
        "metric": "all_mixture_systems_total",
        "count": total_mix_systems,
    })

    # Per property breakdown
    for prop in sorted(mix["property"].unique()):
        sub = mix[mix["property"] == prop]
        bin_sys = int(sub[sub["n_comp"] == 2]["system_id"].nunique())
        mul_sys = int(sub[sub["n_comp"] > 2]["system_id"].nunique())
        cold_sys = int(sub[sub["T_K"] < 263.15]["system_id"].nunique())

        records.append({
            "target": "target_2_mixtures",
            "metric": f"{prop}_binary_systems",
            "count": bin_sys,
        })
        records.append({
            "target": "target_2_mixtures",
            "metric": f"{prop}_multi_systems",
            "count": mul_sys,
        })
        records.append({
            "target": "target_2_mixtures",
            "metric": f"{prop}_cold_systems_lt263K",
            "count": cold_sys,
        })
    return records


def count_pure_candidate_constants(gaps_path: Path) -> list[dict]:
    df = pd.read_csv(gaps_path)
    records = []
    for const in ["has_antoine", "has_lfl", "has_viscosity", "has_density"]:
        records.append({
            "target": "target_3_pure_candidates",
            "metric": f"candidates_{const}",
            "count": int(df[const].sum()),
        })
    all_four = (
        df["has_antoine"]
        & df["has_lfl"]
        & df["has_viscosity"]
        & df["has_density"]
    ).sum()
    records.append({
        "target": "target_3_pure_candidates",
        "metric": "candidates_all_four_constants",
        "count": int(all_four),
    })
    records.append({
        "target": "target_3_pure_candidates",
        "metric": "candidates_total_evaluated",
        "count": len(df),
    })
    return records


def count_cold_viscosity_fuels(dedup_path: Path) -> list[dict]:
    df = pd.read_csv(dedup_path, low_memory=False)
    rf = df[df["kind"] == "real_fuel"]
    v = rf[rf["property"].str.startswith("viscosity")]

    v_253 = v[(v["T_K"] >= 253.0) & (v["T_K"] <= 253.3)]
    v_233 = v[(v["T_K"] >= 233.0) & (v["T_K"] <= 233.3)]

    n_253 = int(v_253["system_id"].nunique())
    n_233 = int(v_233["system_id"].nunique())
    both = set(v_253["system_id"]).intersection(set(v_233["system_id"]))

    return [
        {
            "target": "target_4_cold_viscosity",
            "metric": "fuels_viscosity_253K_minus20C",
            "count": n_253,
        },
        {
            "target": "target_4_cold_viscosity",
            "metric": "fuels_viscosity_233K_minus40C",
            "count": n_233,
        },
        {
            "target": "target_4_cold_viscosity",
            "metric": "fuels_viscosity_both_temps",
            "count": len(both),
        },
    ]


def run_inventory(dedup_path: Path, out_csv: Path) -> pd.DataFrame:
    records = []

    # Target 1: Freezing
    real_fuels_dir = ROOT / "data" / "real_fuels" / "v5"
    freezing_lit_dir = ROOT / "data" / "raw" / "freezing_lit"
    fz_counts = count_real_fuel_freezing(real_fuels_dir, freezing_lit_dir)
    for k, v in fz_counts.items():
        records.append({
            "target": "target_1_freezing",
            "metric": k,
            "count": v,
        })

    # Target 2: Mixtures
    records.extend(count_mixture_systems(dedup_path))

    # Target 3: Pure candidates
    gaps_path = ROOT / "results" / "v5" / "BASE" / "delta_representatives" / "candidate_gaps.csv"
    if gaps_path.exists():
        records.extend(count_pure_candidate_constants(gaps_path))

    # Target 4: Cold viscosity
    records.extend(count_cold_viscosity_fuels(dedup_path))

    inv_df = pd.DataFrame(records)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    inv_df.to_csv(out_csv, index=False)
    return inv_df


def write_readme(inv_df: pd.DataFrame, out_readme: Path):
    lines = [
        "# Data Inventory Before v9 Acquisition",
        "",
        "Counts per target from committed files at the start of v9 (baseline state).",
        "",
        "| Target | Metric | Count |",
        "|---|---|---|",
    ]
    for _, row in inv_df.iterrows():
        lines.append(f"| {row['target']} | {row['metric']} | {row['count']} |")
    lines.append("")
    out_readme.write_text("\n".join(lines), encoding="utf-8")


def main():
    parser = argparse.ArgumentParser(description="Inventory data on disk")
    parser.add_argument(
        "--dedup-path",
        type=Path,
        default=ROOT / "data" / "unified" / "measurements_dedup.csv",
        help="Path to measurements_dedup.csv",
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=ROOT / "results" / "v9" / "inventory" / "before.csv",
        help="Path to output CSV",
    )
    parser.add_argument(
        "--readme",
        type=Path,
        default=ROOT / "results" / "v9" / "inventory" / "README.md",
        help="Path to output README.md",
    )
    args = parser.parse_args()

    inv_df = run_inventory(args.dedup_path, args.out)
    if args.readme:
        write_readme(inv_df, args.readme)
    print(f"Inventory written to {args.out} ({len(inv_df)} rows)")


if __name__ == "__main__":
    main()
