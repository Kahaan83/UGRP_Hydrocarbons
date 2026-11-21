"""Track B step 1: assemble the pure-compound multi-task pool (one row per compound x target, median of sources).

Sources: data/unified pure rows (HC only), AI4fuel staging (pure, HC and oxygenates), v5 ThermoML fusion table, YSI table.
Rules: oxygenates (is_hc False) keep only the octane heads (ron, mon). Tb and LHV are NOT trained: no Tb source in the pool,
and AI4fuel heat_of_combustion values span -22100..-199 with unverified units (MJ/kg label), so LHV is held back.
"""
import json
from pathlib import Path
import numpy as np, pandas as pd
from rdkit import Chem, RDLogger
RDLogger.DisableLog("rdApp.*")
R = Path(__file__).resolve().parents[3]
OUT = R / "results/v6/gnn"
OCT = {"ron", "mon"}

def canon(s):
    m = Chem.MolFromSmiles(s) if isinstance(s, str) else None
    return Chem.MolToSmiles(m) if m else None

def iks(s):
    try: return json.loads(s)
    except Exception: return []

rows = []
# ---- AI4fuel
a = pd.read_csv(R / "data/interim/external_staging/canonical_ai4fuel_measurements.csv", low_memory=False)
ac = pd.read_csv(R / "data/interim/external_staging/canonical_ai4fuel_compounds.csv")
smi = dict(zip(ac.inchikey, ac.smiles))
amap = {"melting_point_C": ("tm_K", lambda v: v + 273.15, (-200, 400)), "flash_point_C": ("flash_C", lambda v: v, (-120, 300)),
        "density_kg_m3": ("density", lambda v: v, (500, 1500)), "dcn": ("dcn", lambda v: v, (0, 120)),
        "cetane_number": ("cetane", lambda v: v, (-20, 120)), "ron": ("ron", lambda v: v, (-20, 130)),
        "mon": ("mon", lambda v: v, (-20, 130))}
for _, r in a.iterrows():
    p = r["property"]
    k = iks(r["comp_inchikeys"])
    if len(k) != 1: continue
    if p == "viscosity_dyn_mPas":
        if abs(r["T_K"] - 293.15) > 1 or not r["value"] > 0: continue
        t, v = "logvisc_mPas_293", np.log10(r["value"])
    elif p in amap:
        t, f, (lo, hi) = amap[p]
        if not lo <= r["value"] <= hi: continue
        v = f(r["value"])
    else: continue
    rows.append((k[0], smi.get(k[0]), t, v, bool(r["is_hc"]), "ai4fuel"))
# ---- unified pure
u = pd.read_csv(R / "data/unified/measurements_dedup.csv", low_memory=False)
uc = pd.read_csv(R / "data/unified/compounds.csv", usecols=["inchikey", "smiles", "is_hc"])
usmi = dict(zip(uc.inchikey, uc.smiles))
u = u[(u.kind == "pure") & u.is_hc]
for _, r in u.iterrows():
    p = r["property"]; k = iks(r["comp_inchikeys"])
    if len(k) != 1: continue
    if p == "density_kg_m3" and abs(r["T_K"] - 293.15) <= 5: t, v = "density", r["value"]
    elif p == "flash_point_C": t, v = "flash_C", r["value"]
    elif p == "dcn": t, v = "dcn", r["value"]
    elif p == "mon": t, v = "mon", r["value"]
    elif p == "viscosity_dyn_mPas" and abs(r["T_K"] - 293.15) <= 1 and r["value"] > 0: t, v = "logvisc_mPas_293", np.log10(r["value"])
    else: continue
    rows.append((k[0], usmi.get(k[0]), t, v, True, "unified"))
# ---- ThermoML fusion (Tm, dHfus)
f = pd.read_csv(R / "data/processed/v5/thermoml_fusion_pure_hc_by_compound.csv")
for _, r in f.iterrows():
    if r["tm_K"] == r["tm_K"]: rows.append((r["inchikey"], r["smiles"], "tm_K", r["tm_K"], True, "thermoml"))
    if r["dhfus_kJ_mol"] == r["dhfus_kJ_mol"]: rows.append((r["inchikey"], r["smiles"], "dhfus_kJ_mol", r["dhfus_kJ_mol"], True, "thermoml"))
# ---- YSI
y = pd.read_csv(R / "data/processed/v4/sooting/ysi_hc.csv")
for _, r in y.iterrows(): rows.append((r["inchikey"], r["smiles"], "ysi", r["ysi"], True, "ysi_hc"))

d = pd.DataFrame(rows, columns=["inchikey", "smiles", "target", "value", "is_hc", "source"])
d["smiles"] = d.smiles.map(canon)
d = d.dropna(subset=["smiles"])
d = d[d.is_hc | d.target.isin(OCT)]
g = d.groupby(["inchikey", "target"]).agg(smiles=("smiles", "first"), value=("value", "median"), n=("value", "size"),
                                         spread=("value", lambda x: x.max() - x.min()), is_hc=("is_hc", "all"),
                                         sources=("source", lambda x: "+".join(sorted(set(x))))).reset_index()
g.to_csv(OUT / "pool.csv", index=False)
print(g.groupby("target").agg(n=("value", "size"), n_oxy=("is_hc", lambda x: int((~x).sum())), med_spread=("spread", "median")).to_string())
print("compounds", g.inchikey.nunique())
