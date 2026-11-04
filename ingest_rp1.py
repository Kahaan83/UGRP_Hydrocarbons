# -*- coding: utf-8 -*-
"""Ingest NIST RP-1 (NISTIR 6646) characterization, surrogates, and properties.

Task T11:
- data/processed/rp1_characterization.csv: 33 components from chemical characterization
  (NISTIR 6646 Tables 3, 4, 5) with mass %, normalized mass fractions, mole fractions,
  chemical formula, MW, SMILES, InChIKey, and hydrocarbon classification. Flag tier: A.
- data/processed/rp1_surrogates.csv: 14-component, 3-component, and 2-component surrogate
  blend specifications with exact fractions, formulas, SMILES, InChIKeys. Flag tier: A.
- data/processed/rp1_properties.csv: 1,480 property measurement rows covering density,
  kinematic viscosity, viscosity*density, and thermal conductivity with canonical units.
"""

from pathlib import Path
import pandas as pd
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

from common import mol_from, hc_class, is_hydrocarbon, mol_weight

ROOT = Path(__file__).resolve().parent
PROCESSED = ROOT / "data" / "processed"
PROCESSED.mkdir(parents=True, exist_ok=True)

# -----------------------------------------------------------------------------
# 1. Chemical characterization (33 components from Tables 3, 4, 5 of NISTIR 6646)
# -----------------------------------------------------------------------------
TIER1_COMPS = [
    ("1", "2,2-dimethylbutane", "000075-83-2", "CCC(C)(C)C", 2.375, 36, "Tier 1 (>=2 mass%)"),
    ("2", "3-methyldecane", "013151-34-3", "CCCCCCCC(C)CC", 3.985, 36, "Tier 1 (>=2 mass%)"),
    ("3", "3-ethyl-4,4-dimethyl-2-pentene", "053907-59-8", "CC=C(CC)C(C)(C)C", 2.726, 36, "Tier 1 (>=2 mass%)"),
    ("4", "2,9-dimethyldecane", "001002-17-1", "CCC(C)CCCCCC(C)C", 6.280, 36, "Tier 1 (>=2 mass%)"),
    ("5", "2-methyl-cis-decalin", "1000152-47-3", "CC1CCC2CCCCC2C1", 3.970, 36, "Tier 1 (>=2 mass%)"),
    ("6", "decahydro-2-methylnaphthalene", "002958-76-1", "CC1CCC2CCCCC2C1", 2.574, 36, "Tier 1 (>=2 mass%)"),
    ("7", "cis-syn-1-methyl-decalin", "1000158-89-1", "CC1CCC2CCCCC2C1", 4.652, 36, "Tier 1 (>=2 mass%)"),
    ("8a", "1-hexyl-3-methylcyclopentane or cyclododecane (coeluting peak 8a)", "061142-68-5 / 000294-62-2", "CCCCCCC1CCCC(C)1", 5.099, 36, "Tier 1 (>=2 mass%)"),
    ("9", "1-dodecene", "000112-41-4", "C=CCCCCCCCCCC", 5.995, 36, "Tier 1 (>=2 mass%)"),
    ("10", "2-methylundecane", "007045-71-8", "CCCCCCCCCC(C)C", 3.124, 36, "Tier 1 (>=2 mass%)"),
    ("11", "3-methylundecane", "001002-43-3", "CCCCCCCCC(C)CC", 2.839, 36, "Tier 1 (>=2 mass%)"),
    ("12", "2,2-dimethyl-deca-3,5-diene", "055638-50-1", "CCCCC=CC=CC(C)(C)C", 2.735, 36, "Tier 1 (>=2 mass%)"),
    ("13", "methylcyclododecane", "NA", "CC1CCCCCCCCCCC1", 3.580, 36, "Tier 1 (>=2 mass%)"),
    ("14", "dodecane", "000112-40-3", "CCCCCCCCCCCC", 5.327, 36, "Tier 1 (>=2 mass%)"),
    ("15", "2,7,10-trimethyldodecane", "074645-98-0", "CCC(C)CCC(C)CCCCC(C)C", 3.765, 36, "Tier 1 (>=2 mass%)"),
]

TIER2_COMPS = [
    ("a", "2,7-dimethyloctane or 2-methylnonane (coeluting/isomer-ambiguous)", "001072-16-8 / 000871-83-0", "CC(C)CCCCC(C)C", 1.329, 37, "Tier 2 (>=1 mass%)"),
    ("b", "cyclodecene", "003717-12-0", "C1=CCCCCCCCC1", 1.610, 37, "Tier 2 (>=1 mass%)"),
    ("c", "cis-decahydronaphthalene (cis-decalin)", "108746-01-6", "C1CCC2CCCCC2C1", 1.174, 37, "Tier 2 (>=1 mass%)"),
    ("d", "?-1,9-dodecadiene", "1000245-71-0", "C=CCCCCCCCC=CC", 1.754, 37, "Tier 2 (>=1 mass%)"),
    ("e", "4-methyl-4-undecene", "061142-40-3", "CCCCCCC=C(C)CCC", 1.663, 37, "Tier 2 (>=1 mass%)"),
    ("f", "x-tridecene (peak f)", "NA", "C=CCCCCCCCCCCC", 1.115, 37, "Tier 2 (>=1 mass%)"),
    ("g", "1-tridecene", "111270-56-1", "C=CCCCCCCCCCCC", 1.241, 37, "Tier 2 (>=1 mass%)"),
    ("h", "heptylcyclohexane", "005617-41-4", "CCCCCCCC1CCCCC1", 1.429, 37, "Tier 2 (>=1 mass%)"),
    ("i", "x-tridecene (peak i)", "NA", "C=CCCCCCCCCCCC", 1.305, 37, "Tier 2 (>=1 mass%)"),
    ("j", "x-tridecene (peak j)", "NA", "C=CCCCCCCCCCCC", 1.977, 37, "Tier 2 (>=1 mass%)"),
    ("k", "2,10-dimethylundecane", "017301-27-8", "CCC(C)CCCCCCC(C)C", 1.507, 37, "Tier 2 (>=1 mass%)"),
    ("l", "x-methyltridecane", "NA", "CCCCCCCCCCCC(C)C", 1.494, 37, "Tier 2 (>=1 mass%)"),
    ("m", "2-methylnaphthalene", "000091-57-6", "Cc1ccc2ccccc2c1", 1.249, 37, "Tier 2 (>=1 mass%)"),
    ("aa", "tridecane", "000629-50-5", "CCCCCCCCCCCCC", 1.080, 37, "Tier 2 (>=1 mass%)"),
]

TIER_LIGHT_COMPS = [
    ("1aa", "methane", "107902-82-8", "C", 0.0, 38, "Tier Light (trace)"),
    ("1a", "nonane", "000111-84-2", "CCCCCCCCC", 0.179, 38, "Tier Light (<1 mass%)"),
    ("1b", "1,3,5-trimethylcyclohexane", "001795-26-2", "CC1CC(C)CC(C)C1", 0.654, 38, "Tier Light (<1 mass%)"),
    ("1d", "2-methyldecane", "006975-98-0", "CCCCCCCCC(C)C", 0.817, 38, "Tier Light (<1 mass%)"),
]

def build_characterization() -> pd.DataFrame:
    all_comps = TIER1_COMPS + TIER2_COMPS + TIER_LIGHT_COMPS
    records = []
    for peak, name, cas, smi, m_pct, pg, tier in all_comps:
        m = mol_from(smi)
        assert m is not None, f"Failed to parse SMILES {smi} for {name}"
        assert is_hydrocarbon(m), f"Non-hydrocarbon {name}"
        records.append({
            "fuel": "RP-1",
            "tier_flag": "A",
            "peak_id": peak,
            "component_name": name,
            "cas": cas,
            "formula": rdMolDescriptors.CalcMolFormula(m),
            "smiles": smi,
            "inchikey": Chem.MolToInchiKey(m),
            "mw": round(mol_weight(m), 3),
            "hc_class": hc_class(m),
            "tier_description": tier,
            "mass_percent_raw": m_pct,
            "mass_fraction_raw": m_pct / 100.0,
            "source": "NISTIR 6646 (2007)",
            "page": pg,
        })
    df = pd.DataFrame(records)
    total_raw_mass = df["mass_fraction_raw"].sum()
    df["mass_fraction_norm"] = df["mass_fraction_raw"] / total_raw_mass
    moles_raw = df["mass_fraction_raw"] / df["mw"]
    df["mole_fraction"] = moles_raw / moles_raw.sum()
    return df

# -----------------------------------------------------------------------------
# 2. Surrogate definitions (14-comp, 3-comp, 2-comp)
# -----------------------------------------------------------------------------
def build_surrogates() -> pd.DataFrame:
    surrogates = []

    # 14-component surrogate from Table 2, p34
    t2_data = [
        ("3-ethyl-4,4-dimethyl-2-pentene", "53907-59-8", "CC=C(CC)C(C)(C)C", 9.98),
        ("cyclodecene", "3618-12-0", "C1=CCCCCCCCC1", 2.11),
        ("2-methylnonane", "871-83-0", "CCCCCCCC(C)C", 2.32),
        ("2-methylnaphthalene", "91-57-6", "Cc1ccc2ccccc2c1", 5.10),
        ("2-methyldecalin", "2958-76-1", "CC1CCC2CCCCC2C1", 22.35),
        ("3-methyldecane", "13151-34-3", "CCCCCCCC(C)CC", 10.84),
        ("1-dodecene", "112-41-4", "C=CCCCCCCCCCC", 2.64),
        ("cyclododecane", "294-62-2", "C1CCCCCCCCCCC1", 4.27),
        ("4-methyl-4-undecene", "61142-40-3", "CCCCCCC=C(C)CCC", 10.45),
        ("n-dodecane", "112-40-3", "CCCCCCCCCCCC", 1.93),
        ("heptylcyclohexane", "5617-41-4", "CCCCCCCC1CCCCC1", 14.22),
        ("1-tridecene", "2437-56-1", "C=CCCCCCCCCCCC", 1.45),
        ("2,7,10-trimethyldodecane", "74645-98-0", "CCC(C)CCC(C)CCCCC(C)C", 10.38),
        ("n-hexadecane", "544-76-3", "CCCCCCCCCCCCCCCC", 1.95),
    ]
    for name, cas, smi, mol_pct in t2_data:
        surrogates.append({
            "surrogate_id": "RP1_NIST_14comp",
            "surrogate_name": "NISTIR 6646 14-Component Surrogate (Table 2)",
            "target_fuel": "RP-1",
            "tier_flag": "A",
            "n_components": 14,
            "component_name": name,
            "cas": cas,
            "smiles": smi,
            "frac_type": "mole_percent",
            "frac_value": mol_pct,
            "mole_fraction": mol_pct / 99.99,
            "reference": "Magee et al. NISTIR 6646 (2007), Table 2, p34",
        })

    # 3-component surrogate (MacDonald / Edwards & Maurice / LLNL)
    s3_data = [
        ("methylcyclohexane", "108-87-2", "CC1CCCCC1", 59.0),
        ("n-dodecane", "112-40-3", "CCCCCCCCCCCC", 32.0),
        ("iso-cetane (2,2,4,4,6,8,8-heptamethylnonane)", "4390-04-9", "CC(C)(C)CC(C)(C)CC(C)CC(C)(C)C", 9.0),
    ]
    for name, cas, smi, mol_pct in s3_data:
        surrogates.append({
            "surrogate_id": "RP1_Edwards_3comp",
            "surrogate_name": "Edwards/Stanford 3-Component RP-1 Pyrolysis Surrogate",
            "target_fuel": "RP-1",
            "tier_flag": "A",
            "n_components": 3,
            "component_name": name,
            "cas": cas,
            "smiles": smi,
            "frac_type": "mole_percent",
            "frac_value": mol_pct,
            "mole_fraction": mol_pct / 100.0,
            "reference": "MacDonald et al. (Stanford/LLNL) / Edwards & Maurice (2001)",
        })

    # 2-component surrogate (binary baseline)
    s2_data = [
        ("methylcyclohexane", "108-87-2", "CC1CCCCC1", 50.0),
        ("n-dodecane", "112-40-3", "CCCCCCCCCCCC", 50.0),
    ]
    for name, cas, smi, mol_pct in s2_data:
        surrogates.append({
            "surrogate_id": "RP1_Binary_2comp",
            "surrogate_name": "RP-1 2-Component Binary Baseline (MCH/n-dodecane 50/50)",
            "target_fuel": "RP-1",
            "tier_flag": "A",
            "n_components": 2,
            "component_name": name,
            "cas": cas,
            "smiles": smi,
            "frac_type": "mole_percent",
            "frac_value": mol_pct,
            "mole_fraction": mol_pct / 100.0,
            "reference": "Edwards & Maurice (2001) / Combustion Literature Binary Baseline",
        })

    df = pd.DataFrame(surrogates)
    fml_list, mw_list, ik_list, cls_list = [], [], [], []
    for _, r in df.iterrows():
        m = mol_from(r["smiles"])
        assert m is not None, f"Failed to parse SMILES {r['smiles']}"
        assert is_hydrocarbon(m), f"Non-hydrocarbon {r['component_name']}"
        fml_list.append(rdMolDescriptors.CalcMolFormula(m))
        mw_list.append(round(mol_weight(m), 3))
        ik_list.append(Chem.MolToInchiKey(m))
        cls_list.append(hc_class(m))

    df["formula"] = fml_list
    df["mw"] = mw_list
    df["inchikey"] = ik_list
    df["hc_class"] = cls_list

    # Compute mass fraction within each surrogate
    mass_fracs = []
    for sid, group in df.groupby("surrogate_id"):
        moles = group["mole_fraction"]
        masses = moles * group["mw"]
        mass_frac = masses / masses.sum()
        for idx in group.index:
            df.loc[idx, "mass_fraction"] = mass_frac.loc[idx]

    return df

# -----------------------------------------------------------------------------
# 3. Property data (already extracted in data/raw/rp1/rp1_blends_properties.csv)
# -----------------------------------------------------------------------------
def build_properties() -> pd.DataFrame:
    raw_csv = ROOT / "data" / "raw" / "rp1" / "rp1_blends_properties.csv"
    assert raw_csv.exists(), f"Raw properties file missing: {raw_csv}"
    df = pd.read_csv(raw_csv)
    df["tier_flag"] = "A"
    # Ensure T_K is present alongside T_C
    df["T_K"] = (df["T_C"] + 273.15).round(3)
    cols = [
        "blend_id", "blend_name", "tier_flag", "property", "value_std",
        "unit_std", "T_C", "T_K", "method", "source_pdf", "page"
    ]
    return df[cols]

def main():
    df_char = build_characterization()
    out_char = PROCESSED / "rp1_characterization.csv"
    df_char.to_csv(out_char, index=False)
    print(f"Wrote {len(df_char)} rows to {out_char}")

    df_surr = build_surrogates()
    out_surr = PROCESSED / "rp1_surrogates.csv"
    df_surr.to_csv(out_surr, index=False)
    print(f"Wrote {len(df_surr)} rows to {out_surr}")

    df_props = build_properties()
    out_props = PROCESSED / "rp1_properties.csv"
    df_props.to_csv(out_props, index=False)
    print(f"Wrote {len(df_props)} rows to {out_props}")

    # Integrity verification
    assert len(df_char) == 33, f"Expected 33 components, got {len(df_char)}"
    assert set(df_surr["surrogate_id"].unique()) == {"RP1_NIST_14comp", "RP1_Edwards_3comp", "RP1_Binary_2comp"}
    assert len(df_props) == 1480, f"Expected 1480 property rows, got {len(df_props)}"
    print("All assertions and self-checks passed successfully.")

if __name__ == "__main__":
    main()
