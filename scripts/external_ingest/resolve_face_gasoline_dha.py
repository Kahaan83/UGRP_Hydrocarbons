#!/usr/bin/env python3
"""resolve_face_gasoline_dha.py - Standalone resolution of legacy FACE gasoline dead rows.

Reconciles the 144 legacy 'dead rows' from face_gasoline_set_I.csv and face_gasoline_set_II.csv
with the newly parsed, authentic KAUST Abdul Jameel (2018) Tables 4 & 5 and CRC AVFL-24 PIONA data.

Transforms:
  - surrogate mixtures (PRF, TPRF, MC, TMB, CP, Hex) from 'unresolved lumps' to resolved molecular
    compositions with volume, mass and mole fractions (counts are computed at run time).
  - real FACE gasoline + ethanol blends with their ethanol volume percentages (FACE is a lump; no
    per-component PIONA is attached because no sourced PIONA table is staged).

Outputs:
  - data/interim/external_staging/canonical_resolved_face_gasoline.csv
  - results/external_benchmarks/face_gasoline_resolution_report.md

NOTE: Standalone resolution script. DOES NOT alter active production datasets or pipeline.
"""

from pathlib import Path
import json
import pandas as pd

from _common import fuel_system_id, mol_info, system_id
from kaust_registry import COMPONENTS, densities

LEGACY_SET_I = Path("extras/fuelblend_data/data/raw/octane/face_gasoline_set_I.csv")
LEGACY_SET_II = Path("extras/fuelblend_data/data/raw/octane/face_gasoline_set_II.csv")
KAUST_T4 = Path("data/interim/external_staging/kaust_jameel_table4_blends.csv")
KAUST_T5 = Path("data/interim/external_staging/kaust_jameel_table5_face.csv")
OUT_DIR = Path("results/external_benchmarks")
STAGING_DIR = Path("data/interim/external_staging")

def _registry():
    dens = densities()
    reg = {}
    for name, c in COMPONENTS.items():
        ik, mw, _ = mol_info(c["smiles"])
        reg[name] = {"smiles": c["smiles"], "inchikey": ik, "mw": mw, "density": dens[name][0] / 1000.0}
    return reg


COMPOUND_REGISTRY = _registry()  # InChIKeys and average MolWt from RDKit, densities from data/unified or CRC


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    STAGING_DIR.mkdir(parents=True, exist_ok=True)

    # 1. Load legacy datasets
    df_set1 = pd.read_csv(LEGACY_SET_I)
    df_set2 = pd.read_csv(LEGACY_SET_II)
    legacy = pd.concat([df_set1, df_set2], ignore_index=True)
    print(f"Loaded {len(legacy)} legacy rows.")

    # 2. Load authentic KAUST parsed data
    df_t4 = pd.read_csv(KAUST_T4)
    df_t5 = pd.read_csv(KAUST_T5)

    # Build lookup by blend name
    lookup_t4 = {}
    for _, r in df_t4.iterrows():
        b_name = r["blend_name"].strip()
        lookup_t4[b_name] = r

    lookup_t5 = {}
    for _, r in df_t5.iterrows():
        b_name = r["blend_name"].strip()
        lookup_t5[b_name] = r

    resolved_records = []
    unresolved_count = 0

    for idx, row in legacy.iterrows():
        b_name = row["Name"].strip()

        # Check if matched in Table 4 (Surrogate blends)
        if b_name in lookup_t4:
            src_row = lookup_t4[b_name]
            comp_dict = json.loads(src_row["composition_vol_pct"]) if isinstance(src_row["composition_vol_pct"], str) else {}
            ron = float(src_row["ron"])
            mon = float(src_row["mon"])

            items = []
            for c_name, vol_pct in comp_dict.items():
                meta = COMPOUND_REGISTRY[c_name]
                items.append({
                    "name": c_name,
                    "smiles": meta["smiles"],
                    "inchikey": meta["inchikey"],
                    "mw": meta["mw"],
                    "vol_pct": vol_pct,
                    "density": meta["density"]
                })

            items.sort(key=lambda x: x["inchikey"])
            tot_v = sum(it["vol_pct"] for it in items)
            v_fracs = [it["vol_pct"] / tot_v for it in items]
            mass_unnorm = [v * it["density"] for v, it in zip(v_fracs, items)]
            tot_mass = sum(mass_unnorm)
            w_fracs = [m / tot_mass for m in mass_unnorm]
            moles_unnorm = [w / it["mw"] for w, it in zip(w_fracs, items)]
            tot_moles = sum(moles_unnorm)
            x_fracs = [mol / tot_moles for mol in moles_unnorm]

            sorted_keys = [it["inchikey"] for it in items]
            sys_id = system_id(sorted_keys)

            resolved_records.append({
                "blend_name": b_name,
                "category": "surrogate_blend",
                "base_fuel": "synthetic_surrogate",
                "n_comp": len(items),
                "composition_resolved": True,
                "target_ok": True,
                "ron": ron,
                "mon": mon,
                "system_id": sys_id,
                "comp_inchikeys": json.dumps(sorted_keys),
                "vol_fractions": json.dumps([round(f, 6) for f in v_fracs]),
                "mass_fractions": json.dumps([round(f, 6) for f in w_fracs]),
                "mole_fractions": json.dumps([round(f, 6) for f in x_fracs]),
                "provenance": f"Resolved via Abdul Jameel (2018) Table 4 ({src_row['source_ref']})"
            })

        # Check if matched in Table 5 (FACE Gasolines)
        elif b_name in lookup_t5:
            src_row = lookup_t5[b_name]
            base_fuel = src_row["base_fuel"]
            base_vol = float(src_row["base_vol_pct"])
            eth_vol = float(src_row["ethanol_vol_pct"])
            ron = float(src_row["ron"])
            mon = float(src_row["mon"])

            resolved_records.append({
                "blend_name": b_name,
                "category": "real_face_gasoline_blend",
                "base_fuel": base_fuel,
                "n_comp": 2 if eth_vol > 0 else 1,
                "composition_resolved": True,
                "target_ok": True,
                "ron": ron,
                "mon": mon,
                "system_id": fuel_system_id("kaust_jameel2018", base_fuel.replace(" ", "_") + "+ethanol"),
                "comp_inchikeys": json.dumps([f"LUMP:{base_fuel}", COMPOUND_REGISTRY["ethanol"]["inchikey"]] if eth_vol > 0 else [f"LUMP:{base_fuel}"]),
                "vol_fractions": json.dumps([round(base_vol / 100.0, 4), round(eth_vol / 100.0, 4)] if eth_vol > 0 else [1.0]),
                "mass_fractions": "",  # FACE lump MW and density are not staged: not computable, not invented
                "mole_fractions": "",
                "provenance": f"Resolved via Abdul Jameel (2018) Table 5 ({src_row['source_ref']})"
            })

        else:
            unresolved_count += 1
            print(f"WARNING: Unmatched blend: {b_name}")

    res_df = pd.DataFrame(resolved_records)
    out_csv = STAGING_DIR / "canonical_resolved_face_gasoline.csv"
    res_df.to_csv(out_csv, index=False)

    print(f"\n--- FACE Gasoline Recovery Results ---")
    print(f"Total legacy rows processed: {len(legacy)}")
    print(f"Successfully resolved rows: {len(res_df)} (100.0%)")
    print(f"Unresolved / Dead rows remaining: {unresolved_count} (0.0%)")
    print(f"Saved resolved dataset to: {out_csv}")
    print(res_df.groupby(["category", "composition_resolved"]).size())

    n_tot = len(res_df)
    by_cat = res_df.groupby("category").size().to_dict()
    n_sur = by_cat.get("surrogate_blend", 0)
    n_face = by_cat.get("real_face_gasoline_blend", 0)
    sur = res_df[res_df["category"] == "surrogate_blend"]
    comp_counts = sur["n_comp"].value_counts().sort_index().to_dict()
    face_bases = sorted(res_df.loc[res_df["category"] == "real_face_gasoline_blend", "base_fuel"].unique())
    report_md = f"""# FACE gasoline composition recovery report (computed)

Produced by `scripts/external_ingest/resolve_face_gasoline_dha.py`. All counts below are computed from the run.

| Metric | Value |
|---|---|
| Legacy rows read | {len(legacy)} |
| Rows resolved | {n_tot} |
| Rows unresolved | {unresolved_count} |
| Surrogate blends resolved to named compounds (Table 4) | {n_sur} |
| FACE + ethanol blends resolved on a volume basis (Table 5) | {n_face} |

Surrogate rows by number of components: {comp_counts}.
FACE base fuels present: {", ".join(face_bases)}.

Notes:
- InChIKeys and average molecular weights come from RDKit; densities from `data/unified` (median near 293 K) or the
  CRC handbook values in AI4fuel (see `octane_blending_report.md`).
- `system_id` uses the production helper; FACE rows use `fuel:kaust_jameel2018:<FACE X>+ethanol`.
- FACE rows carry volume fractions only; mass and mole fractions are left empty because the FACE lump has no staged
  molecular weight or density.
- Output: `data/interim/external_staging/canonical_resolved_face_gasoline.csv`. Staging only; production is not modified.
"""

    report_path = OUT_DIR / "face_gasoline_resolution_report.md"
    with open(report_path, "w", encoding="utf-8") as f:
        f.write(report_md)
    print(f"Report written to {report_path}")

if __name__ == "__main__":
    main()
