"""Assemble every open multicomponent (n_comp >= 3) liquidus / wax-disappearance row into the canonical measurement schema and count
scorable systems per class mix (verified vs ThermoML-measured Tm/dHfus). Sourcing only: no gate re-run, no margin, no model scoring.
Run: PYTHONPATH=src .venv/Scripts/python.exe scripts/v5/base/ws2_multicomp_data_ingest.py
"""
from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "src"))
spec = importlib.util.spec_from_file_location("ws2", ROOT / "scripts/v5/base/delta_ws2_freeze_multicomp.py")
W = importlib.util.module_from_spec(spec)
spec.loader.exec_module(W)
OUT = ROOT / "results/v5/BASE/delta_ws2_multicomp_data"
GATE = 8


def main():
    d, dropped = W.dedup(W.load_rows())
    tier, ver = W.constants()
    mws = {}
    recs = []
    for r in d.itertuples():
        m = [Descriptors.MolWt(Chem.MolFromSmiles(s)) for s in r.smiles]
        x = np.array(r.x)
        w = x * np.array(m)
        w = w / w.sum()
        order = np.argsort(r.inchikeys)
        T_C = r.T_K - 273.15
        recs.append(dict(
            source={"thermoml": "thermoml", "oakley2017": "literature", "ascent065A": "literature", "shariatrad2022": "literature",
                    "dauphin_softx2018": "literature"}[r.source],
            property="freezing_point_C", value=round(T_C, 3), T_K=round(r.T_K, 3), P_kPa=101.3, kind="multi", n_comp=r.n_comp,
            is_hc=True, system_id=r.system_id, comp_inchikeys=json.dumps([r.inchikeys[i] for i in order]),
            x_mole=json.dumps([round(float(x[i]), 6) for i in order]), w_mass=json.dumps([round(float(w[i]), 6) for i in order]),
            fuel_id="", ref=f"{r.source}:{r.doi}:{r.sample}", raw_value=round(T_C, 3), raw_unit="C",
            method=r.method, class_mix=W.FM.class_mix(r.smiles), smiles=json.dumps([r.smiles[i] for i in order])))
    m = pd.DataFrame(recs)
    # exact duplicate guard (system, T 0.5 K bin, x 3 dp)
    m["_k"] = m.system_id + "|" + (m.T_K / 0.5).round().astype(str) + "|" + m.x_mole.map(lambda s: str([round(v, 3) for v in json.loads(s)]))
    n0 = len(m)
    m = m.drop_duplicates("_k").drop(columns="_k").reset_index(drop=True)
    m.to_csv(OUT / "measurements_multicomp_liquidus.csv", index=False)

    rows = []
    for sid, g in m.groupby("system_id"):
        iks = json.loads(g.comp_inchikeys.iloc[0])
        smi = json.loads(g.smiles.iloc[0])
        miss_any = [s for s, k in zip(smi, iks) if k not in tier]
        miss_ver = [s for s, k in zip(smi, iks) if k not in ver]
        nc = [sum(a.GetSymbol() == "C" for a in Chem.MolFromSmiles(s).GetAtoms()) for s in smi]
        rows.append(dict(system_id=sid, class_mix=g.class_mix.iloc[0], method=g.method.iloc[0], source=g.ref.iloc[0].split(":")[0],
                         n_comp=len(iks), n_rows=len(g), max_carbon=max(nc), jet_range_le_C20=max(nc) <= 20,
                         scorable_verified_or_thermoml=not miss_any, scorable_verified_only=not miss_ver,
                         missing_any=json.dumps(miss_any), missing_verified=json.dumps(miss_ver)))
    s = pd.DataFrame(rows).sort_values(["class_mix", "source"])
    s.to_csv(OUT / "systems.csv", index=False)

    out = []
    for cm in ["alkane_only", "alkane_naphthene", "with_aromatic"]:
        t = s[s.class_mix == cm]
        nw = t[t.method == "wax_appearance"]
        out.append(dict(class_mix=cm, systems_in_set=len(t), of_which_wax_appearance_C18_C36=len(nw),
                        scorable_verified_or_thermoml=int(t.scorable_verified_or_thermoml.sum()),
                        scorable_verified_only=int(t.scorable_verified_only.sum()),
                        scorable_excl_wax_appearance=int(t[t.method != "wax_appearance"].scorable_verified_or_thermoml.sum()),
                        ceiling_if_all_constants_curated=len(t), ceiling_excl_wax_appearance=len(t) - len(nw),
                        gate=GATE, gap_to_gate_now=GATE - int(t.scorable_verified_or_thermoml.sum()),
                        gate_reachable_from_open_data=len(t) >= GATE and len(t) - len(nw) >= GATE))
    c = pd.DataFrame(out)
    c.to_csv(OUT / "counts_by_class_mix.csv", index=False)
    allmiss = sorted({x for v in s.missing_any for x in json.loads(v)})
    print(f"rows {n0}->{len(m)}; dedup-vs-thermoml dropped {len(dropped)}")
    print(c.to_string(index=False))
    print("constants missing (any tier):", allmiss)


if __name__ == "__main__":
    main()
