"""Task 14: limiting-species freezing baseline. Predicted freezing point = max over species of the ideal-solubility temperature T_i with
ln x_i = -(dHfus_i / R)(1/T - 1/Tm_i), i.e. 1/T_i = 1/Tm_i - R ln x_i / dHfus_i. Constants are the designated ones used by baseline (b) of v6.
Dev rows only (the 52 dev entities of results/v6/freezing/features_preds.csv); no fitting, so LOGO and in-sample coincide.
  python scripts/v7/freezing/baseline_d.py"""
import sys, json
from pathlib import Path
import numpy as np, pandas as pd

ROOT = Path(__file__).resolve().parents[3]
R = 8.314462618


def limiting_species_C(x, consts):
    """x: mole fractions; consts: per species (tm_K, dh_fus_J_mol) or None (does not crystallise). Returns degC, NaN if no species crystallises."""
    best = -np.inf
    for xi, c in zip(np.asarray(x, float), consts):
        if c is None or xi <= 0:
            continue
        tm, dh = c
        best = max(best, 1.0 / (1.0 / tm - R * np.log(xi) / dh))
    return float(best - 273.15) if np.isfinite(best) else float("nan")


def main():
    sys.path.insert(0, str(ROOT / "scripts/v6/freezing")); sys.path.insert(0, str(ROOT / "src"))
    import common as C
    from fuelblend.api.predict_route import _tmap_freezing
    ds = C.dataset()
    feat = pd.read_csv(ROOT / "results/v6/freezing/features_preds.csv", dtype={"fuel_id": str})
    dev = feat[feat.role == "dev"].reset_index(drop=True)
    ents = C.entities(ds)
    g1 = C.load("g1_v7bd", "scripts/v5/eval_ajftd_gcxgc_v5.py")
    d5 = C.load("d5_v7bd", "scripts/v5/base/delta_freeze_t2_unifac.py")
    tdes = _tmap_freezing()
    pred = []
    for e in dev.entity:
        x, iks, smi, _ = d5.s5_inputs(ents[e], g1)
        consts = []
        for ik in iks:
            t = tdes.get(ik)
            ok = bool(t) and np.isfinite(t.get("tm_K", np.nan)) and np.isfinite(t.get("dh_fus_J_mol", np.nan))
            consts.append((float(t["tm_K"]), float(t["dh_fus_J_mol"])) if ok else None)
        pred.append(limiting_species_C(x, consts))
    dev["pred_limiting"] = pred
    ref = pd.read_csv(ROOT / "results/v6/freezing/dev_oof_predictions.csv")
    assert (ref.entity.to_numpy() == dev.entity.to_numpy()).all() and np.allclose(ref.y, dev.freezing_C)
    y = dev.freezing_C.to_numpy()
    cols = {"a_constant": ref.a_constant, "b_unifac_svl": ref.b_unifac_svl, "c1_weighted_tm": ref.c1_weighted_tm, "c2_ridge": ref.c2_ridge,
            "limiting_species": pd.Series(pred)}
    ok = np.isfinite(pred)
    rows = []
    rng = np.random.default_rng(0)
    idx = rng.integers(0, len(y), (4000, len(y)))
    e_const = np.abs(cols["a_constant"].to_numpy() - y)
    for m, p in cols.items():
        p = p.to_numpy(float)
        e = np.abs(p - y)
        d = e_const - e
        est = d[idx].mean(1)
        rows.append(dict(model=m, n=len(y), n_nonfinite=int((~np.isfinite(p)).sum()), mae=float(np.nanmean(e)), bias=float(np.nanmean(p - y)),
                         mean_over_pred=float(np.nanmean(np.clip(p - y, 0, None))), frac_over=float(np.nanmean(p > y)),
                         p90_signed=float(np.nanpercentile(p - y, 90)), gain_vs_constant=float(np.nanmean(d)),
                         ci_lo=float(np.nanpercentile(est, 5)), ci_hi=float(np.nanpercentile(est, 95))))
    out = ROOT / "results/v7/freezing/baseline_d"
    out.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(out / "summary.csv", index=False, float_format="%.5g")
    dev[["entity", "freezing_C", "pred_limiting"]].to_csv(out / "predictions.csv", index=False, float_format="%.5g")
    print(pd.DataFrame(rows).round(3).to_string())


if __name__ == "__main__":
    main()
