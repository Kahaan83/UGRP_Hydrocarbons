"""v5 delta_ws3_pure_gaps, post-result annotation (no tier, gate or value changes).

The declared tier rule (README) used tm_dhfus_verified.csv and ThermoML only. data/processed/v5/pure_props_v5.csv also holds Tm / dHfus rows
with status 'verified' (inverse_design_v2.CANDIDATE_REGISTRY, purchasable_compounds.csv, WebBook). This script adds, without changing any tier,
columns tm_pure_props_v5 / dh_pure_props_v5 ("<status>: <value>") to gap_table.csv and fusion_ml_estimated.csv, so a reader can see where an
ML estimate sits next to a value listed elsewhere in the repo.
Run after delta_ws3_fusion.py: PYTHONPATH=src .venv/Scripts/python.exe scripts/v5/base/delta_ws3_annotate_pure_props.py
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd
from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")
ROOT = Path(__file__).resolve().parents[3]
PP = ROOT / "data/processed/v5/pure_props_v5.csv"
FILES = [ROOT / "results/v5/BASE/delta_ws3_pure_gaps/gap_table.csv", ROOT / "data/processed/v5/fusion_ml_estimated.csv"]


def ik(s):
    m = Chem.MolFromSmiles(s) if isinstance(s, str) else None
    return Chem.MolToInchiKey(m) if m is not None else None


def lookup() -> dict:
    p = pd.read_csv(PP)
    p = p[p.prop.isin(["Tm", "dHfus"])].copy()
    p["inchikey"] = p.smiles.map(ik)
    out = {}
    for r in p.itertuples():
        unit = "C" if r.prop == "Tm" else "kJ/mol"
        out.setdefault((r.inchikey, r.prop), []).append(f"{r.status}: {r.value:g} {unit}")
    return {k: "; ".join(v) for k, v in out.items()}


def main():
    lk = lookup()
    for f in FILES:
        d = pd.read_csv(f)
        d["tm_pure_props_v5"] = d.inchikey.map(lambda k: lk.get((k, "Tm"), ""))
        d["dh_pure_props_v5"] = d.inchikey.map(lambda k: lk.get((k, "dHfus"), ""))
        d.to_csv(f, index=False)
        cols = ["label", "tm_pure_props_v5", "dh_pure_props_v5"]
        if "source" in d:
            hit = d[(d.tm_pure_props_v5.str.contains("verified")) | (d.dh_pure_props_v5.str.contains("verified"))]
            print(f.name, "ml rows with a 'verified' pure_props_v5 value:", len(hit))
            print(hit[cols + ["tm_K", "dhfus_kJ_mol"]].to_string())
        else:
            print(f.name, "rows with any pure_props_v5 value:", int(((d.tm_pure_props_v5 != "") | (d.dh_pure_props_v5 != "")).sum()), "of", len(d))


if __name__ == "__main__":
    main()
