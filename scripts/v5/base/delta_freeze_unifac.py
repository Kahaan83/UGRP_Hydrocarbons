"""delta_freeze_unifac (results/v5/BASE/delta_freeze_unifac/README.md): UNIFAC-SvL on scope S_F (R-FU-1) and all verified ThermoML (R-FU-2).

    python scripts/v5/base/delta_freeze_unifac.py
"""
from __future__ import annotations

import importlib.util
import json
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "src"))
warnings.simplefilter("ignore")
OUT = ROOT / "results" / "v5" / "BASE" / "delta_freeze_unifac"

from fuelblend.baselines.flash_unifac import make_unifac_gamma_fn  # noqa: E402
from fuelblend.eval import freeze_multicomp as FM  # noqa: E402
from fuelblend.eval.composition_route import build_tmap_verified  # noqa: E402


def _load(name, rel):
    s = importlib.util.spec_from_file_location(name, ROOT / rel)
    m = importlib.util.module_from_spec(s)
    s.loader.exec_module(m)
    return m


FMG = _load("fm_margin", "scripts/v5/base/delta_freeze_margin.py")
OLDM = FMG.OLDM


def main():
    tmap = build_tmap_verified()
    sle = pd.read_csv(FMG.SLE_CSV)
    sle = sle[sle.point_kind == "liquidus"].reset_index(drop=True)
    rows = FMG.scope_rows(tmap)                     # same row order as sle (one row per liquidus point)
    assert len(rows) == len(sle)
    pu, fbs = [], []
    for r, s in zip(rows.itertuples(), sle.itertuples()):
        if not r.all_verified:
            pu.append(np.nan)
            fbs.append(False)
            continue
        smi, iks = json.loads(s.smiles), json.loads(s.inchikeys)
        x = [float(v) for v in json.loads(s.x_liquid)]
        x = list(np.asarray(x) / np.sum(x))
        cc = [tmap[k] for k in iks]
        cc = [(float(c["tm_K"]), float(c["dh_fus_J_mol"])) for c in cc]
        gfn, fb = make_unifac_gamma_fn(smi)
        pu.append(FM.svl_activity_C(x, cc, gfn))
        fbs.append(bool(fb))
    rows["pred_unifac"] = pu
    rows["unifac_fallback"] = fbs
    rows.to_csv(OUT / "rows.csv", index=False, float_format="%.5f")
    S = []
    sf = rows[rows.in_SF].copy()
    allv = rows[rows.all_verified & rows.all_hc].copy()
    for scope, d in (("S_F", sf), ("thermoml_all_verified_any_T_any_nC", allv)):
        for model, col in (("svl_verified", "pred_verified"), ("svl_unifac", "pred_unifac")):
            S.append(OLDM.summarize(d, scope, model, col, "under"))
    summ = pd.DataFrame(S)
    summ.to_csv(OUT / "summary.csv", index=False, float_format="%.4f")
    pd.set_option("display.width", 250)
    print(summ[["scope", "model", "n_scored", "n_systems", "mae", "bias", "margin_alpha10", "frac_unsafe", "max_group_unsafe"]].to_string())
    sf["r_v"] = sf.pred_verified - sf.measured
    sf["r_u"] = sf.pred_unifac - sf.measured
    g = sf.groupby(["group", "sample"]).agg(n=("r_v", "size"), worst_ideal=("r_v", "min"), worst_unifac=("r_u", "min"),
                                             mae_ideal=("r_v", lambda v: v.abs().mean()), mae_unifac=("r_u", lambda v: v.abs().mean()),
                                             fallback=("unifac_fallback", "any")).sort_values("worst_ideal")
    g.to_csv(OUT / "per_system_S_F.csv", float_format="%.4f")
    print(g.to_string())
    a = summ[(summ.scope == "S_F") & (summ.model == "svl_verified")].iloc[0]
    b = summ[(summ.scope == "S_F") & (summ.model == "svl_unifac")].iloc[0]
    pref = bool(b.margin_alpha10 < a.margin_alpha10 and b.mae <= a.mae and not bool(sf.unifac_fallback.any()))
    rd = dict(R_FU_1=dict(ideal_margin=float(a.margin_alpha10), unifac_margin=float(b.margin_alpha10), ideal_mae=float(a.mae),
                          unifac_mae=float(b.mae), any_fallback=bool(sf.unifac_fallback.any()), unifac_preferred=pref))
    json.dump(rd, open(OUT / "readings.json", "w"), indent=1)
    print(json.dumps(rd, indent=1))


if __name__ == "__main__":
    main()
