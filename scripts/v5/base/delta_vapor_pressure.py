"""v5 BASE Task 16a: AV-33-22 vapour pressure against Raoult + Antoine on the mapped measured composition.

Outputs (results/v5/BASE/delta_vapor_pressure/): rows_vp.csv, delta.csv, flash_vs_vp.csv, readings.json
"""
from __future__ import annotations

import json
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "src"))
from fuelblend.eval.cellkit import frame, make_row  # noqa: E402

BASE = ROOT / "results" / "v5" / "BASE"
OUT = BASE / "delta_vapor_pressure"
V5 = ROOT / "data" / "real_fuels" / "v5"
COV_MIN = 0.90


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    warnings.simplefilter("ignore")
    from fuelblend.baselines.flash_unifac import _lookup
    from fuelblend.data.reconstruction_r3 import SubReps
    from fuelblend.eval.composition_route import RouteCtx, kind_bins, map_bins
    from fuelblend.mixture.common import CompoundRegistry, PurePropertyProvider
    reg = CompoundRegistry(ROOT / "data/unified/compounds.csv")
    prov = PurePropertyProvider(reg, ROOT / "data/cache/pure_grid.npz", verbose=False)
    reps = SubReps(reg, prov, vp_lookup=lambda s: _lookup(s, None))
    ctx = RouteCtx(reg, prov)
    ant = ctx.flash_reg.drop_duplicates("inchikey").set_index("inchikey")[["antoine_A", "antoine_B", "antoine_C"]]
    comp = pd.read_csv(V5 / "crc_av3322_composition.csv")
    props = pd.read_csv(V5 / "crc_av3322_properties.csv")
    vp = props[(props.property == "vapor_pressure_kPa") & props.qualifier.isna()]
    fu = pd.read_csv(BASE / "t2_fuels.csv")
    fu = fu[fu.source == "av3322"].set_index("fuel_id")
    rows = []
    for fid, g in vp.groupby("fuel_id"):
        mf = map_bins(kind_bins("av3322", comp[comp.fuel_id == fid]), reps)
        iks = [reg.ik[i] for i in mf.idx]
        x = mf.w / np.array([ctx.mw(reg.smiles[i]) for i in mf.idx])
        x = x / x.sum()
        a = ant.reindex(iks).to_numpy(float)
        have = np.isfinite(a).all(axis=1)
        cov = float(x[have].sum())
        for r in g.itertuples():
            pred = float(np.sum(x[have] * 10 ** (a[have, 0] - a[have, 1] / (r.T_K + a[have, 2])))) if cov >= COV_MIN else np.nan
            status = "ok" if cov >= COV_MIN and r.value > 0 else ("antoine_coverage_lt_0.90" if cov < COV_MIN else "measured_le_0")
            rows.append(dict(fuel_id=fid, entity=fu.loc[fid, "entity"], group=fu.loc[fid, "group"], chosen=bool(fu.loc[fid, "chosen_for_group"]),
                             T_K=round(float(r.T_K), 2), measured=r.value, pred=pred, antoine_mole_coverage=cov, status=status))
    d = pd.DataFrame(rows)
    d["err_ln"] = np.where(d.status == "ok", np.log(d.pred / d.measured), np.nan)
    d.to_csv(OUT / "rows_vp.csv", index=False, float_format="%.5g")
    ok = d[(d.status == "ok") & d.chosen]
    cells = [make_row("vapor_pressure_kPa", "T2_measured", "av3322", "raoult_antoine", "designated", "recomputed_v5", "descriptive",
                      "data/real_fuels/v5/crc_av3322_properties.csv (D6378)", err=g.err_ln, groups=g.group,
                      extra=dict(T_K=T, baseline_mae=np.nan, baseline_bias=np.nan, baseline_unsafe_rate=np.nan,
                                 n_excluded_status=int(((d.T_K == T) & (d.status != "ok")).sum())))
             for T, g in ok.groupby("T_K")]
    frame(cells).to_csv(OUT / "delta.csv", index=False)
    fl = pd.read_csv(BASE / "rows_t2_flash.csv")
    fl = fl[(fl.source == "av3322") & fl.chosen.astype(bool) & (fl.status == "ok")][["entity", "measured", "pred_designated"]]
    fl["flash_err"] = fl.pred_designated - fl.measured
    m = fl.merge(ok[ok.T_K == 313.15][["entity", "err_ln"]], on="entity")
    m.to_csv(OUT / "flash_vs_vp.csv", index=False, float_format="%.5g")
    rho = float(spearmanr(m.flash_err, m.err_ln).correlation) if len(m) > 3 else np.nan
    res = dict(n_rows_ok=int((d.status == "ok").sum()), n_rows_excluded=int((d.status != "ok").sum()), n_pairs=len(m), spearman_flash_vs_lnvp313=rho,
               reading_pass=bool(rho > 0.5), bias_ln_by_T={f"{c['T_K']:.2f}": c["bias"] for c in cells}, flash_bias_av3322=float(fl.flash_err.mean()))
    (OUT / "readings.json").write_text(json.dumps(res, indent=2))
    print(json.dumps(res, indent=2))


if __name__ == "__main__":
    main()
