#!/usr/bin/env python3
"""ingest_kaust_jameel.py - Stage KAUST octane data (Abdul Jameel et al. 2018) in canonical schema.

Reads the parsed Tables 3, 4, 5 from data/interim/external_staging/. Writes (staging only):
  - kaust_jameel_canonical_measurements.csv
  - kaust_jameel_compounds.csv

Rules:
  - InChIKeys come from RDKit via the curated SMILES table in kaust_registry.py.
  - system_id: production helper (fuelblend.data.thermoml_sle.system_id) for named-compound systems;
    fuel:<source>:<sample> for FACE gasoline lumps.
  - Oxygenates are allowed for octane only; is_hc tags them.
  - Blend compositions are NOT renormalised silently: the volume-sum is recorded and checked.
Does not touch data/unified/.
"""
from __future__ import annotations

import json

import pandas as pd

from _common import STAGING_DIR, UNIFIED_COMP, fuel_system_id, mol_info, system_id
from kaust_registry import COMPONENTS, TABLE3_SMILES, densities, load_tables, pure_props

SOURCE = "kaust_jameel2018"
TIER = "literature_compiled"
P_KPA = 101.325


def _kind(n):
    return "pure" if n == 1 else ("binary" if n == 2 else "multi")


def _row(prop, val, kind, n_comp, is_hc, sid, iks, x, w, fuel_id, ref, family, extra=None):
    r = {"source": SOURCE, "property": prop, "value": float(val), "T_K": float("nan"), "P_kPa": P_KPA,
         "kind": kind, "n_comp": n_comp, "is_hc": bool(is_hc), "system_id": sid,
         "comp_inchikeys": json.dumps(iks), "x_mole": json.dumps(x), "w_mass": json.dumps(w),
         "fuel_id": fuel_id, "ref": ref, "raw_value": float(val), "raw_unit": "ON",
         "provenance_tier": TIER, "family": family}
    if extra:
        r.update(extra)
    return r


def main():
    t3, t4, t5 = load_tables()
    dens = densities()
    keys, mw, hc = {}, {}, {}
    for name, c in COMPONENTS.items():
        ik, m, h = mol_info(c["smiles"])
        keys[name], mw[name], hc[name] = ik, m, h

    rows, comp_rows = [], {}
    prod_iks = set(pd.read_csv(UNIFIED_COMP, usecols=["inchikey"])["inchikey"])

    # ---- Table 3: pure compounds (all 128 mapped via RDKit) ----
    for r in t3.itertuples():
        smi = TABLE3_SMILES[int(r.s_no)]
        ik, m, h = mol_info(smi)
        sid = system_id([ik])
        comp_rows[ik] = {"inchikey": ik, "smiles": smi, "name": r.name, "mw_avg": m, "is_hc": h,
                         "table3_s_no": int(r.s_no), "in_production_compounds": ik in prod_iks}
        for prop, val in (("ron", r.ron), ("mon", r.mon)):
            rows.append(_row(prop, val, "pure", 1, h, sid, [ik], [1.0], [1.0], "",
                             f"Jameel2018:Table3:{r.source_ref}", ik))

    # ---- Table 4: blends; single-component (+0) rows are pure ----
    comp_by_name = {k: v for k, v in keys.items()}
    sum_log = []
    for r in t4.itertuples():
        comp = json.loads(r.composition_vol_pct)
        unknown = [c for c in comp if c not in COMPONENTS]
        if unknown:
            raise ValueError(f"{r.blend_name}: unmapped components {unknown}")
        tot = sum(comp.values())
        sum_log.append((r.blend_name, tot))
        if abs(tot - 100.0) > 1.5:
            raise ValueError(f"{r.blend_name}: volume percentages sum to {tot}, refusing to renormalise")
        items = sorted(((comp_by_name[c], c, v / tot) for c, v in comp.items()))
        iks = [i[0] for i in items]
        wr = [i[2] * dens[i[1]][0] for i in items]
        w = [x / sum(wr) for x in wr]
        nr = [wi / mw[i[1]] for wi, i in zip(w, items)]
        x = [n / sum(nr) for n in nr]
        non_eth = sorted(c for c in comp if c != "ethanol")
        family = "+".join(non_eth)
        is_hc = all(hc[c] for c in comp)
        n = len(items)
        extra = {"vol_sum_pct": round(tot, 4), "vol_frac": json.dumps([round(i[2], 6) for i in items])}
        for prop, val in (("ron", r.ron), ("mon", r.mon)):
            rows.append(_row(prop, val, _kind(n), n, is_hc, system_id(iks), iks, [round(v, 6) for v in x],
                             [round(v, 6) for v in w], str(r.blend_name),
                             f"Jameel2018:Table4:{r.source_ref}", family, extra))
        for c in comp:
            ik = comp_by_name[c]
            comp_rows.setdefault(ik, {"inchikey": ik, "smiles": COMPONENTS[c]["smiles"], "name": c,
                                      "mw_avg": mw[c], "is_hc": hc[c], "table3_s_no": COMPONENTS[c]["t3"],
                                      "in_production_compounds": ik in prod_iks})

    # ---- Table 5: FACE gasoline + ethanol (volume basis; FACE is a lump, no per-component mass/mole) ----
    eth = keys["ethanol"]
    for r in t5.itertuples():
        base = r.base_fuel.replace(" ", "_")
        sid = fuel_system_id(SOURCE, f"{base}+ethanol")
        has_eth = r.ethanol_vol_pct > 0
        for prop, val in (("ron", r.ron), ("mon", r.mon)):
            rows.append(_row(prop, val, "binary" if has_eth else "pure", 2 if has_eth else 1, not has_eth,
                             sid, [eth] if has_eth else [], [], [], r.blend_name,
                             f"Jameel2018:Table5:{r.source_ref}", f"{base}",
                             {"vol_sum_pct": r.base_vol_pct + r.ethanol_vol_pct,
                              "vol_frac": json.dumps([round(r.ethanol_vol_pct / 100.0, 4)]),
                              "ethanol_vol_pct": r.ethanol_vol_pct}))

    df = pd.DataFrame(rows)
    df.to_csv(STAGING_DIR / "kaust_jameel_canonical_measurements.csv", index=False)
    pd.DataFrame(list(comp_rows.values())).to_csv(STAGING_DIR / "kaust_jameel_compounds.csv", index=False)

    print(f"Staged {len(df)} KAUST octane rows; {len(comp_rows)} compounds")
    print(df.groupby(["kind"]).size())
    print("Table 4 volume-sum range:", min(s for _, s in sum_log), max(s for _, s in sum_log))
    print("Densities:", {k: (round(v[0], 1), v[1]) for k, v in dens.items()})
    print("Pure octane props:", pure_props(t3, t4))


if __name__ == "__main__":
    main()
