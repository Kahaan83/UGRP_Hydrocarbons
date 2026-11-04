#!/usr/bin/env python
"""Hydrocarbon-only viscosity table (pure + mixtures) from the NIST/TDE data
shipped in CheMixHub (datasets/nist-logV).

UNIT NOTE (verified): the 'value' column is labelled 'Log viscosity / cP' but is
actually ln(eta / Pa.s). Check: hexane 298.15 K -> exp(-8.1258)=0.296 mPa.s,
water 0.891, toluene 0.553, hexadecane 3.03 mPa.s - all match literature.

Output (data/processed/):
  viscosity_pure_hc.csv    one row per (compound, T, P)
  viscosity_mix_hc.csv     one row per mixture state; compositions as JSON lists
  compounds_hc.csv         every hydrocarbon found: smiles, class, C number, MW

Dynamic viscosity only. Refutas needs KINEMATIC viscosity -> divide by density
(rho, kg/m3) at the same T: nu[cSt] = eta[mPa.s] / rho[g/cm3].
"""
from __future__ import annotations
import ast, json, sys
from pathlib import Path

import numpy as np
import pandas as pd

from common import hc_class, carbon_number, is_hydrocarbon, mol_from, mol_weight, mole_to_mass

SRC = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(
    "data/raw/git/chemixhub/datasets/nist-logV/processed_data")
OUT = Path("data/processed")
OUT.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(SRC / "processed_NISTlogV.csv")
cmp = pd.read_csv(SRC / "compounds.csv")

# ---- hydrocarbon compound table
rows = {}
for cid, smi in zip(cmp.compound_id, cmp.smiles):
    m = mol_from(smiles=smi)
    if is_hydrocarbon(m):
        rows[int(cid)] = dict(compound_id=int(cid), smiles=smi, cls=hc_class(m),
                              n_carbon=carbon_number(m), mw=mol_weight(m))
comp = pd.DataFrame(rows.values())
comp.to_csv(OUT / "compounds_hc.csv", index=False)
hc_ids = set(rows)
print(f"hydrocarbons: {len(comp)} of {len(cmp)} compounds")

# ---- keep rows whose every component is a hydrocarbon
df["ids"] = df.cmp_ids.apply(ast.literal_eval)
df["x"] = df.cmp_mole_fractions.apply(ast.literal_eval)
df = df[df.ids.apply(lambda l: all(i in hc_ids for i in l))].copy()

df["ln_eta_Pas"] = df["value"]                      # see UNIT NOTE
df["eta_mPas"] = np.exp(df["ln_eta_Pas"]) * 1e3
df["n_comp"] = df.ids.apply(len)
df = df[(df.eta_mPas > 0) & np.isfinite(df.eta_mPas)]


def to_mass(r):
    M = [rows[i]["mw"] for i in r.ids]
    return json.dumps([round(float(w), 6) for w in mole_to_mass(r.x, M)])


pure = df[df.n_comp == 1].copy()
pure["compound_id"] = pure.ids.apply(lambda l: l[0])
pure = pure.merge(comp, on="compound_id")
pure[["compound_id", "smiles", "cls", "n_carbon", "mw", "T", "P", "eta_mPas",
      "ln_eta_Pas", "Visc_Unc", "Ref_ID"]].rename(columns={"T": "T_K", "P": "P_kPa"}) \
    .to_csv(OUT / "viscosity_pure_hc.csv", index=False)

mix = df[df.n_comp > 1].copy()
mix["smiles_list"] = mix.ids.apply(lambda l: json.dumps([rows[i]["smiles"] for i in l]))
mix["x_mole"] = mix.x.apply(json.dumps)
mix["w_mass"] = mix.apply(to_mass, axis=1)
mix["classes"] = mix.ids.apply(lambda l: json.dumps([rows[i]["cls"] for i in l]))
mix[["smiles_list", "x_mole", "w_mass", "classes", "n_comp", "T", "P", "eta_mPas",
     "ln_eta_Pas", "Visc_Unc", "Ref_ID"]].rename(columns={"T": "T_K", "P": "P_kPa"}) \
    .to_csv(OUT / "viscosity_mix_hc.csv", index=False)

print(f"pure rows: {len(pure)}  (compounds: {pure.compound_id.nunique()})")
print(f"mixture rows: {len(mix)}  by n_comp: {mix.n_comp.value_counts().sort_index().to_dict()}")
print(f"T range (K): {df['T'].min():.1f} - {df['T'].max():.1f}")
