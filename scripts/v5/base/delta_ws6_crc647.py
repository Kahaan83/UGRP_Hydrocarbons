"""WS6 application tier: models trained on the 145 measured-GCxGC fuels (split-G training set, all fuels) applied to the
CRC 647 jet fuels whose compositions are reconstructed (r3, measured hydrocarbon types + reconstructed carbon split; NOT GCxGC).
Label: fold_b / dev / exploratory.  CRC 647 never enters training.  Declared in README (application tier).
Cells: flash (v4 F ideal_C), ln viscosity at -20 and -40 C (raw Refutas from v4 R3 V).  Writes crc647_application.csv / .json.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem

sys.path.insert(0, str(Path(__file__).resolve().parent))
import delta_ws6_gcxgc as W  # noqa: E402

ROOT, OUT = W.ROOT, W.OUT


def classify(smi: str) -> tuple[str, int]:
    m = Chem.MolFromSmiles(smi)
    nC = sum(a.GetSymbol() == "C" for a in m.GetAtoms())
    ri = m.GetRingInfo()
    arom_rings = sum(all(m.GetAtomWithIdx(i).GetIsAromatic() for i in r) for r in ri.AtomRings())
    rings = ri.NumRings()
    if arom_rings >= 2:
        return "diaromatic", nC
    if arom_rings == 1:
        return ("alkylbenzene" if rings == 1 else "cycloaromatic"), nC
    if rings >= 3:
        return "tricycloparaffin", nC
    if rings == 2:
        return "dicycloparaffin", nC
    if rings == 1:
        return "monocycloparaffin", nC
    if any(b.GetBondType() != Chem.BondType.SINGLE for b in m.GetBonds()):
        return "alkene", nC
    branched = any(a.GetSymbol() == "C" and a.GetDegree() > 2 for a in m.GetAtoms())
    return ("iso_paraffin" if branched else "n_paraffin"), nC


def crc_features() -> pd.DataFrame:
    c = pd.read_csv(ROOT / "data/real_fuels/v4/compositions_crc647_r3.csv")
    c = c[c.smiles.notna()].copy()
    c["fuel_id"] = c.system_id.str.replace("fuel:", "", regex=False)
    rows = {}
    for fid, g in c.groupby("fuel_id"):
        v = np.zeros(len(W.FEATS))
        for r in g.itertuples():
            k, n = classify(r.smiles)
            v[W.CLASSES.index(k) * len(W.CN_BINS) + W.CN_BINS.index(W.cn_bin(n))] += r.fraction
        rows[fid] = v / v.sum()
    return pd.DataFrame.from_dict(rows, orient="index", columns=W.FEATS)


def main() -> None:
    t = pd.read_csv(OUT / "fuel_table.csv")
    F = pd.read_csv(OUT / "features.csv", index_col="entity")
    Fc = crc_features()
    fl = pd.read_csv(ROOT / "results/v4/R3/F/flash_real_fuels.csv")
    fl = fl[fl.family == "crc647"].set_index("fuel_id")
    vs = pd.read_csv(ROOT / "results/v4/R3/V/viscosity_predictions_v4.csv")
    rows = []
    cells = {
        "flash": [(f, float(fl.loc[f].measured_C), float(fl.loc[f].ideal_C), bool(fl.loc[f].get("statistics_allowed", True))) for f in fl.index if f in Fc.index and np.isfinite(fl.loc[f].ideal_C)],
    }
    for cell, T in (("visc_253", 253.15), ("visc_233", 233.15)):
        s = vs[(vs.T_K.round(2) == T)]
        cells[cell] = [(r.fuel_id, float(np.log(r.y)), float(np.log(r["pred_Raw Refutas"])), bool(~r.locked)) for _, r in s.iterrows() if r.fuel_id in Fc.index]
    for cell, items in cells.items():
        d = t[t.cell == cell]
        X, y, p = F.loc[d.entity].to_numpy(float), d.y.to_numpy(float), d.p.to_numpy(float)
        Xt = Fc.loc[[i[0] for i in items]].to_numpy(float)
        yt = np.array([i[1] for i in items]); pt = np.array([i[2] for i in items])
        pr = W.fit_predict(X, y, p, Xt, pt, True)
        for m, v in pr.items():
            e = v - yt
            rows.append(dict(cell=cell, model=m, n=len(items), mae=float(np.abs(e).mean()), bias=float(e.mean()),
                             over_gt5C=float((e > 5).mean()) if cell == "flash" else np.nan,
                             under_gt_0p4154=float((-e > 0.4154).mean()) if cell.startswith("visc") else np.nan))
    r = pd.DataFrame(rows)
    r.to_csv(OUT / "crc647_application.csv", index=False)
    print(r.round(3).to_string())


if __name__ == "__main__":
    main()
