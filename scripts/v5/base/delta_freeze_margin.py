"""v5 delta_freeze_margin: freezing (ideal SvL) scored with verified Tm / dHfus constants; margin feasibility.

Declared readings: results/v5/BASE/delta_freeze_margin/README.md (committed before this script existed).
Run: PYTHONPATH=src .venv/Scripts/python.exe scripts/v5/base/delta_freeze_margin.py
Writes rows.csv, summary.csv, readings.json.  Writes NOTHING to configs/ (the yaml is written by build_uncertainty_v5_q2.py only if the rule passes).
"""
from __future__ import annotations

import importlib.util
import json
import sys
import warnings
from pathlib import Path
from unittest import mock

import numpy as np
import pandas as pd
from rdkit import Chem

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "src"))
OUT = ROOT / "results" / "v5" / "BASE" / "delta_freeze_margin"
OLD = ROOT / "results" / "v5" / "BASE" / "delta_newdata_freeze_smoke"

from fuelblend.api import predict_route as pr                                     # noqa: E402
from fuelblend.eval.composition_route import build_tmap_verified, ideal_sle_C      # noqa: E402

_s = importlib.util.spec_from_file_location("old_fs", ROOT / "scripts/v5/base/delta_newdata_freeze_smoke.py")
OLDM = importlib.util.module_from_spec(_s)
_s.loader.exec_module(OLDM)

SLE_CSV = ROOT / "data" / "processed" / "v5" / "thermoml_sle_hc.csv"


def _nc(smi):
    return sum(a.GetSymbol() == "C" for a in Chem.MolFromSmiles(smi).GetAtoms())


def scope_rows(tmap):
    d = pd.read_csv(SLE_CSV)
    d = d[d.point_kind == "liquidus"].reset_index(drop=True)
    out = []
    for r in d.itertuples():
        smi, iks = json.loads(r.smiles), json.loads(r.inchikeys)
        x = [float(v) for v in json.loads(r.x_liquid)]
        ncs = [_nc(s) for s in smi]
        all_hc = bool(r.is_hc)
        all_c8_16 = all(8 <= n <= 16 for n in ncs)
        all_verified = all(k in tmap for k in iks)
        cold = float(r.T_K) <= 273.15
        in_sf = bool(all_hc and all_c8_16 and all_verified and cold)
        pred = float(ideal_sle_C(x, iks, tmap)) if all_verified else np.nan
        out.append(dict(system_id=r.system_id, group=r.system_id, doi=r.doi, sample=" + ".join(json.loads(r.names)), T_K=float(r.T_K),
                        measured=float(r.T_K) - 273.15, n_comp=int(r.n_comp), all_hc=all_hc, all_c8_c16=all_c8_16, all_verified=all_verified,
                        cold=cold, in_SF=in_sf, pred_verified=pred, censored=False,
                        status=("ok" if np.isfinite(pred) else ("not_all_verified" if not all_verified else "non_finite"))))
    return pd.DataFrame(out)


def add_old(rows):
    old = pd.read_csv(OLD / "rows_freezing.csv")
    old = old[old.dataset == "thermoml"][["system_id", "T_K", "sample", "pred_svl"]].rename(columns={"pred_svl": "pred_old"})
    return rows.merge(old, on=["system_id", "T_K", "sample"], how="left")


def summ(d, scope, model, pred_col):
    return OLDM.summarize(d, scope, model, pred_col, "under")


def descriptive(tmap):
    """Vozka and ASCENT with tmap_verified (verified-only) through predict_fuel."""
    with mock.patch.object(pr, "_tmap", lambda: tmap):
        asc = OLDM.ascent_rows()
        voz = OLDM.vozka_rows()
    return asc, voz


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    tmap = build_tmap_verified()
    rows = add_old(scope_rows(tmap))
    sf = rows[rows.in_SF].copy()
    S = []
    S.append(summ(sf, "S_F", "svl_verified", "pred_verified"))
    sf_old = sf[np.isfinite(sf.pred_old.astype(float))]
    S.append(summ(sf_old, "S_F_rows_old_finite", "svl_old_coutinho", "pred_old"))
    S.append(summ(sf_old, "S_F_rows_old_finite", "svl_verified", "pred_verified"))
    allv = rows[rows.all_verified & rows.all_hc]
    S.append(summ(allv, "thermoml_all_verified_any_T_any_nC", "svl_verified", "pred_verified"))   # descriptive
    ocold = rows[rows.all_hc & rows.all_c8_c16 & rows.cold & np.isfinite(rows.pred_old.astype(float))]
    S.append(summ(ocold, "old_jet_cold_scorable_rows", "svl_old_coutinho", "pred_old"))
    asc, voz = descriptive(tmap)
    asc = asc[asc.scope_kind == "ascent_mixtures"]
    S.append(summ(asc, "ascent_mixtures", "svl_verified_only", "pred_svl"))
    S.append(summ(voz, "vozka", "svl_verified_only", "pred_svl"))
    S.append(summ(voz[np.isfinite(voz.pred_f1x.astype(float))], "vozka", "f1x_variant_unchanged", "pred_f1x"))
    summary = pd.DataFrame(S)

    r = summary[(summary.scope == "S_F") & (summary.model == "svl_verified")].iloc[0].to_dict()
    boot = OLDM.bootstrap_diff(sf[np.isfinite(sf.pred_verified)], "pred_verified") if r["n_scored"] else None
    fm1 = bool(r["n_systems"] >= 9)
    fm2 = bool(r["n_scored"] > 0 and r["mae"] < r["mae_const"])
    fm3 = bool(np.isfinite(r["margin_alpha10"]))
    readings = {
        "S_F": dict(n_rows_scope=int(r["n_rows_scope"]), n_scored=int(r["n_scored"]), n_unscored=int(r["n_unscored"]), n_systems=int(r["n_systems"]),
                    systems=sorted(sf[np.isfinite(sf.pred_verified)].group.unique().tolist()),
                    bias_C=float(r["bias"]), max_group_under_prediction_C=float(r["max_group_unsafe"]),
                    measured_range_C=[float(r["measured_min"]), float(r["measured_max"])]),
        "R-FM-1": dict(n_systems=int(r["n_systems"]), threshold=9, passed=fm1),
        "R-FM-2": dict(mae_C=float(r["mae"]), const_mae_C=float(r["mae_const"]), macro_mae_C=float(r["mae_macro"]), macro_const_C=float(r["mae_const_macro"]),
                       bootstrap95_const_minus_svl=boot, interval_excludes_zero=bool(boot and boot[0] > 0), passed=fm2),
        "R-FM-3": dict(margin_C=float(r["margin_alpha10"]), finite=fm3, passed=fm3, direction="under", alpha=0.10),
        "rule_write_margin": bool(fm1 and fm2 and fm3),
    }
    rows.drop(columns=["censored"]).to_csv(OUT / "rows.csv", index=False)
    summary.to_csv(OUT / "summary.csv", index=False)
    (OUT / "readings.json").write_text(json.dumps(readings, indent=2, default=float), encoding="utf-8")
    pd.set_option("display.width", 250, "display.max_columns", 40)
    print(summary[["scope", "model", "n_rows_scope", "n_scored", "n_unscored", "n_systems", "mae", "bias", "mae_const", "margin_alpha10", "max_group_unsafe"]].to_string())
    print(json.dumps(readings, indent=2, default=float))


if __name__ == "__main__":
    warnings.simplefilter("ignore")
    main()
