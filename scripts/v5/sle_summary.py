"""Counts for the ThermoML SLE extraction (results/v5/BASE/data_thermoml_sle/README.md).  No scoring."""
import json
import sys

import pandas as pd
from rdkit import Chem

SRC = sys.argv[1] if len(sys.argv) > 1 else "data/processed/v5/thermoml_sle_hc.csv"
d = pd.read_csv(SRC)
cmp_ = pd.read_csv("data/unified/compounds.csv", usecols=["inchikey"])
known = set(cmp_.inchikey)


def feats(smi):
    m = Chem.MolFromSmiles(smi)
    nc = sum(a.GetSymbol() == "C" for a in m.GetAtoms())
    arom = any(a.GetIsAromatic() for a in m.GetAtoms())
    unsat = any(b.GetBondTypeAsDouble() > 1 and not b.GetIsAromatic() for b in m.GetBonds())
    ring = m.GetRingInfo().NumRings() > 0
    branched = any(a.GetDegree() > 2 for a in m.GetAtoms())
    jet = 8 <= nc <= 16 and (arom or not unsat)
    nalk = nc >= 10 and not ring and not unsat and not branched
    return jet, nalk


rows = []
for sid, g in d.groupby("system_id"):
    smi = json.loads(g.smiles.iloc[0])
    iks = json.loads(g.inchikeys.iloc[0])
    f = [feats(s) for s in smi]
    rows.append(dict(system_id=sid, n_comp=len(smi), jet=any(x[0] for x in f), nalk10=any(x[1] for x in f),
                     all_known=all(k in known for k in iks), has_liquidus=(g.point_kind == "liquidus").any(),
                     n_rows=len(g), names=" + ".join(json.loads(g.names.iloc[0])),
                     nalk10_crystal=False))
s = pd.DataFrame(rows)
print("rows", len(d), "| by kind", d.point_kind.value_counts().to_dict(), "| files", d.file.nunique(), "| dois", d.doi.nunique())
print("rows by n_comp", d.n_comp.value_counts().sort_index().to_dict())
print("systems", len(s), "| by n_comp", s.n_comp.value_counts().sort_index().to_dict(), "| >=3 comp", int((s.n_comp >= 3).sum()))
print("systems with liquidus rows", int(s.has_liquidus.sum()))
print("systems with a jet-range compound (C8-C16 sat./cyclo/aromatic)", int(s.jet.sum()), "| rows", int(d.system_id.isin(s[s.jet].system_id).sum()))
print("systems with an n-alkane C10+", int(s.nalk10.sum()), "| rows", int(d.system_id.isin(s[s.nalk10].system_id).sum()))
print("systems with all components in unified/compounds.csv", int(s.all_known.sum()), "| rows", int(d.system_id.isin(s[s.all_known].system_id).sum()))
print("  ... and with a jet-range compound", int((s.all_known & s.jet).sum()), "| and n-alkane C10+", int((s.all_known & s.nalk10).sum()))
print("T_K range (liquidus)", d[d.point_kind == "liquidus"].T_K.min(), d[d.point_kind == "liquidus"].T_K.max())
print("systems with >=3 components:", s[s.n_comp >= 3].names.tolist())
print("n-alkane C10+ systems:", s[s.nalk10].names.tolist())
print("jet-range, all-known systems:", s[s.jet & s.all_known].names.tolist())
