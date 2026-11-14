"""WS6 post-hoc control (NOT declared in the README; labelled post hoc): C2s = physics + per-source offset fitted in-fold,
compared with the declared M1 on split G, same folds/repeats as delta_ws6_gcxgc.run_G (same seeds).  Tests whether M1's split-G gain is
just source identity (AV-33-22 and AJFTD sit in different composition regions and carry different measurement offsets).
Writes posthoc_source_offset.csv.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import delta_ws6_gcxgc as W  # noqa: E402


def main() -> None:
    t = pd.read_csv(W.OUT / "fuel_table.csv")
    oof = pd.read_csv(W.OUT / "oof_G.csv")
    m1 = oof[oof.model == "M1"].assign(ae=lambda d: d.err.abs()).groupby(["cell", "entity"]).ae.mean()
    rows = []
    for cell in W.CELLS:
        d = t[t.cell == cell].reset_index(drop=True)
        y, p, grp, fam = d.y.to_numpy(float), d.p.to_numpy(float), d.group.to_numpy(), d.family.to_numpy()
        ug = np.unique(grp)
        ae = np.zeros((W.N_REP, len(d)))
        for rep in range(W.N_REP):
            rng = np.random.default_rng(1000 + rep)
            fold_of = dict(zip(rng.permutation(ug), np.arange(len(ug)) % W.K_FOLDS))
            fo = np.array([fold_of[g] for g in grp])
            for f in range(W.K_FOLDS):
                te = fo == f
                if te.sum() == 0 or (~te).sum() < 8:
                    continue
                r = y - p
                for s in np.unique(fam[te]):
                    sel = te & (fam == s)
                    tr = (~te) & (fam == s)
                    off = r[tr].mean() if tr.sum() >= 3 else r[~te].mean()
                    ae[rep, sel] = np.abs(p[sel] + off - y[sel])
        c2s = ae.mean(0)
        mm = m1.loc[cell].reindex(d.entity).to_numpy()
        diff, lo, hi, p1, _ = W.boot(c2s - mm)
        rows.append(dict(cell=cell, n=len(d), mae_C2s=c2s.mean(), mae_M1=mm.mean(), d_C2s_minus_M1=diff, lo=lo, hi=hi, p1=p1))
    r = pd.DataFrame(rows)
    r.to_csv(W.OUT / "posthoc_source_offset.csv", index=False)
    print(r.round(3).to_string())


if __name__ == "__main__":
    main()
