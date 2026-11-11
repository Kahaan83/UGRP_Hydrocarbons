"""Build a model table: one row per (system_id, T, P, property) from measurements_dedup.csv."""
import json
import pandas as pd
import numpy as np

DEDUP = "data/unified/measurements_dedup.csv"

def _json_loads(s):
    if isinstance(s, (list, dict)):
        return s
    if pd.isna(s):
        return []
    try:
        return json.loads(s)
    except Exception:
        return []

def build_model_table(path=DEDUP):
    df = pd.read_csv(path, low_memory=False)
    if "is_hc" in df.columns:
        df = df[df["is_hc"] == True].copy()
    
    # Rule R1: pure compounds are always eligible; blend rows are eligible only if composition is available (L3 or L2).
    # L1 (class totals) and L0 (no composition) are excluded from the model table and saved to excluded_by_R1.csv.
    clevel = df["composition_level"] if "composition_level" in df.columns else pd.Series("L3", index=df.index)
    kind = df["kind"]
    
    # Eligibility mask:
    # pure -> always True
    # blend (binary, multi, real_fuel) -> True iff composition_level in ('L2', 'L3')
    is_pure = (kind == "pure")
    has_comp = clevel.isin(["L2", "L3"])
    eligible_mask = is_pure | has_comp
    
    excluded_df = df[~eligible_mask].copy()
    excluded_df.to_csv("data/processed/excluded_by_R1.csv", index=False)
    print(f"Rule R1: routed {len(excluded_df)} L1/L0 blend rows to data/processed/excluded_by_R1.csv")
    
    df = df[eligible_mask].copy()
    clevel = clevel[eligible_mask]

    mt = pd.DataFrame({
        "property": df["property"],
        "system_id": df["system_id"],
        "T_K": df["T_K"],
        "P_kPa": df["P_kPa"],
        "kind": df["kind"],
        "n_comp": df["n_comp"],
        "composition_level": clevel,
        "comp_coverage": df.get("comp_coverage", 1.0),
        "comp_basis": df.get("comp_basis", "mole"),
        "value": df["value"],
        "vmin": df["vmin"],
        "vmax": df["vmax"],
        "n_obs": df["n_obs"],
        "n_refs": df["n_refs"],
        "n_sources": df["n_sources"],
        "sources": df["sources"],
        "n_outlier": df["n_outlier"],
        "conflict": df["conflict"],
        "comp_inchikeys": df["comp_inchikeys"].apply(_json_loads).apply(json.dumps),
        "x_mole": df["x_mole"].apply(_json_loads).apply(json.dumps),
        "w_mass": df["w_mass"].apply(_json_loads).apply(json.dumps),
    })
    mt["quality_weight"] = 1.0 / np.sqrt(df["n_obs"].clip(lower=1).astype(float))
    mt["tier"] = df["sources"].apply(_tier_of_sources)
    mt["method"] = "unknown"  # filled per-source by loaders later
    mt["domain"] = mt["property"].map(_domain)
    return mt


def _domain(prop):
    p = str(prop)
    if "viscosity" in p:
        return "viscosity"
    if "density" in p:
        return "density"
    if "flash" in p:
        return "flash"
    if "pour" in p or "freez" in p or "cloud" in p:
        return "cold_flow"
    if p in ("dcn", "mon"):
        return "combustion"
    if "gcv" in p or "heat" in p:
        return "gcv"
    if "smoke" in p:
        return "smoke"
    if "vapor" in p:
        return "vapor_pressure"
    return "other"

# Tier rules (ANALYSIS_PLAN §4): A = measured, method+conditions stated, peer-reviewed or gov;
# B = handbook/compilation or derived-with-stated-rule; C = estimates/ML-derived (never train).
TIER_OF_SOURCE = {
    "leenhouts": "A", "nist_tde": "A", "thermoml": "A", "ajftd": "A",
    "fortin2012": "A", "pure_freezing_nist": "A", "saldana11": "A",
    "crc647": "B", "av2718": "B", "av1109": "B", "chew2024": "B",
    "flashpoint_li": "A", "pubchem": "C", "lfl_estimated": "C",
    # added 2026-11-09; B because of known lab issues (Ha & Lee pure dodecane reads ~9 C high; Pires WSU flash and density are biased)
    "ha_lee_2016": "B", "pires2018": "B", "vozka2019": "B",
}

def _tier_of_sources(sources):
    tiers = [TIER_OF_SOURCE.get(s, "unknown") for s in str(sources).split("+")]
    order = {"C": 0, "unknown": 1, "B": 2, "A": 3}
    return min(tiers, key=lambda t: order.get(t, 1))


if __name__ == "__main__":
    mt = build_model_table()
    mt.to_csv("data/processed/model_table.csv", index=False)
    print("rows:", len(mt))
    print(mt["domain"].value_counts().to_string())
