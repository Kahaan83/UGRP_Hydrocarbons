"""v5 BASE Task 13: pooled T2 residual model (physics + ridge on composition features), leave-source-out.

Outputs (results/v5/BASE/delta_t2_residual/): delta.csv, readings.json
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
from fuelblend.eval.cellkit import frame, make_row  # noqa: E402
from fuelblend.eval.composition_route import kind_bins  # noqa: E402
from fuelblend.mixture.t2_residual import bin_features, fit_predict_residual  # noqa: E402

BASE = ROOT / "results" / "v5" / "BASE"
OUT = BASE / "delta_t2_residual"
V5 = ROOT / "data" / "real_fuels" / "v5"
S1 = {"T_lt_263K": 0.4153576080066075, "T_263_283K": 0.6084381907110359}
DIRECTIONS = [("av3322", ["ajftd_gcxgc"]), ("ajftd_gcxgc", ["av3322"]), ("fuellib", ["av3322", "ajftd_gcxgc"])]   # (test, train sources)

_s = importlib.util.spec_from_file_location("cells_t2_bulk_res", ROOT / "scripts/v5/base/cells_t2_bulk.py")
bulk = importlib.util.module_from_spec(_s)
_s.loader.exec_module(bulk)


def features() -> pd.DataFrame:
    fu = pd.read_csv(BASE / "t2_fuels.csv")
    comps = {"av3322": pd.read_csv(V5 / "crc_av3322_composition.csv"), "fuellib": pd.read_csv(V5 / "fuellib_gcxgc_compositions.csv"),
             "ajftd": pd.read_csv(V5 / "ajftd_gcxgc_composition.csv")}
    X = {}
    for r in fu.itertuples():
        c = comps[r.source]
        fid = int(r.fuel_id) if r.source == "ajftd" else r.fuel_id
        X[r.entity] = bin_features(kind_bins(r.source, c[c.fuel_id == fid]))
    return pd.DataFrame.from_dict(X, orient="index")


def targets() -> dict[str, pd.DataFrame]:
    """name -> rows with entity, source, group, chosen, measured, pred, log (bool), T_K."""
    out = {}
    fl = pd.read_csv(BASE / "rows_t2_flash.csv")
    fl = fl[(fl.status == "ok") & ~fl.low_coverage.astype(bool)].rename(columns={"pred_designated": "pred"})
    out["flash_point_C"] = fl.assign(T_K=np.nan, log=False)
    bk = pd.read_csv(BASE / "rows_t2_bulk.csv")
    bk = bk[(bk.variant == "designated") & (bk.status == "ok") & ~bk.low_coverage.astype(bool)]
    for name, prop, T, lg in (("density_kg_m3@288.15", "density_kg_m3", 288.15, False), ("viscosity_kin_cSt@313.15", "viscosity_kin_cSt", 313.15, True),
                              ("viscosity_kin_cSt@253.15", "viscosity_kin_cSt", 253.15, True), ("viscosity_kin_cSt@233.15", "viscosity_kin_cSt", 233.15, True)):
        out[name] = bk[(bk.property == prop) & (np.round(bk.T_K, 2) == T)].assign(log=lg)
    out["heat_of_combustion_net_MJ_kg"] = bk[bk.property == "heat_of_combustion_net_MJ_kg"].assign(log=False)
    s5 = pd.read_csv(ROOT / "results/v5/S5/s5_predictions.csv")
    fu = pd.read_csv(BASE / "t2_fuels.csv").set_index("entity")
    s5 = s5[(s5.role == "headline") & s5.entity.isin(fu.index) & s5.F1.notna() & s5.freezing_C.notna()].copy()
    src = {"av3322": "av3322", "fuellib": "fuellib", "ajftd": "ajftd_gcxgc"}
    s5["source"] = [src[fu.loc[e, "source"]] for e in s5.entity]
    s5["chosen"] = [bool(fu.loc[e, "chosen_for_group"]) for e in s5.entity]
    out["freezing_point_C"] = s5.rename(columns={"freezing_C": "measured", "F1": "pred"}).assign(T_K=np.nan, log=False)
    return out


def err_of(d, pred):
    return np.log(pred) - np.log(d.measured) if d.log.iloc[0] else pred - d.measured


def unsafe_for(name, T):
    if name.startswith("flash"):
        return lambda e: e > 5.0
    if name.startswith("visc") and T < 283.15:
        m = S1["T_lt_263K" if T < 263.15 else "T_263_283K"]
        return lambda e, m=m: e < -m
    return None


def run_direction(d: pd.DataFrame, X: pd.DataFrame, test_src: str, train_srcs: list[str]) -> pd.DataFrame | None:
    te = bulk.one_per_group(d[d.source == test_src])
    tr = d[d.source.isin(train_srcs) & ~d.group.isin(set(te.group)) & d.entity.isin(X.index)]
    tr = pd.concat([bulk.one_per_group(g) for _, g in tr.groupby("source")]) if len(tr) else tr
    te = te[te.entity.isin(X.index)]
    if len(tr) < 5 or te.empty:
        return None
    assert not set(tr.group) & set(te.group)
    r_tr = -err_of(tr, tr.pred)                                                    # residual = measured - physics (ln for viscosity)
    rhat = fit_predict_residual(X.loc[tr.entity].to_numpy(), r_tr.to_numpy(), tr.group.to_numpy(), X.loc[te.entity].to_numpy())
    te = te.copy()
    te["err_raw"] = err_of(te, te.pred)
    te["err_model"] = te.err_raw + rhat
    te["n_train_groups"] = tr.group.nunique()
    return te


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    warnings.simplefilter("ignore")
    X = features()
    rows, res = [], {}
    for name, d in targets().items():
        T = float(d.T_K.iloc[0]) if len(d) and np.isfinite(d.T_K.iloc[0]) else np.nan
        prop = name.split("@")[0]
        ok = []
        for test_src, train_srcs in DIRECTIONS:
            te = run_direction(d, X, test_src, train_srcs)
            if te is None:
                rows.append(dict(make_row(prop, "T2_measured", test_src, "physics+ridge_residual", "t2_residual", "no_evaluation", "exploratory", ""),
                                 reason="fewer than 5 training fuels or no test fuels", T_K=T, train_sources="+".join(train_srcs)))
                if test_src != "fuellib":
                    ok.append(False)
                continue
            un = unsafe_for(name, T)
            raw = make_row(prop, "T2_measured", test_src, "raw", "", "recomputed_v5", "exploratory", "", err=te.err_raw, groups=te.group, unsafe=un)
            r = make_row(prop, "T2_measured", test_src, "physics+ridge_residual", "t2_residual", "recomputed_v5", "exploratory",
                         "rows_t2_flash.csv / rows_t2_bulk.csv / S5 s5_predictions.csv", err=te.err_model, groups=te.group, unsafe=un,
                         err_trivial=te.err_raw, trivial_model="raw_physics",
                         extra=dict(T_K=T, train_sources="+".join(train_srcs), n_train_groups=int(te.n_train_groups.iloc[0]),
                                    baseline_mae=raw["mae"], baseline_bias=raw["bias"], baseline_unsafe_rate=raw["unsafe_rate"]))
            rows.append(r)
            if test_src != "fuellib":
                no_worse = not (np.isfinite(r["unsafe_rate"]) and np.isfinite(raw["unsafe_rate"]) and r["unsafe_rate"] > raw["unsafe_rate"] + 1e-12)
                ok.append(r["verdict"] == "beats_trivial" and no_worse)
        res[name] = dict(reading_pass=bool(ok) and all(ok))
    d = frame(rows)
    d.to_csv(OUT / "delta.csv", index=False)
    (OUT / "readings.json").write_text(json.dumps(res, indent=2))
    print(d[["property", "T_K", "set", "status", "n_groups", "n_train_groups", "baseline_mae", "mae", "baseline_bias", "bias",
             "baseline_unsafe_rate", "unsafe_rate", "ci_lo", "ci_hi", "verdict"]].to_string())
    print(json.dumps(res, indent=2))


if __name__ == "__main__":
    main()
