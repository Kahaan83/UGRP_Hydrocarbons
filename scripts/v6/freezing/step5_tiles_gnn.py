"""Track A post-hoc step 5 (exploratory, DEV ONLY): fill the missing isomer constants of the registry with Track B predictions and
rerun baseline (e) isomer-averaged tiles. The Track A holdout was scored once and is NOT used here.

Fill order for each registry compound lacking Tm and/or dHfus in the designated tmap: (1) measured value in the Track B pool,
(2) estimate. Variants: e2 = Track B D-MPNN for Tm and dHfus; e3 = v5 fusion_estimate (fitted on the pool Tm) for Tm, D-MPNN for dHfus.
Estimates never replace designated values. No fuel label enters the pool (pure compounds only), so there is no target leakage."""
import sys, json
from pathlib import Path
import numpy as np, pandas as pd, torch
sys.path.insert(0, str(Path(__file__).resolve().parent))
import common as C                                        # noqa: E402
import step3_models as S3                                 # noqa: E402
from fuelblend.eval import composition_route as CR        # noqa: E402
from fuelblend.ml import dmpnn as D, molpool as MP        # noqa: E402
from fuelblend.thermo import fusion_estimate as F         # noqa: E402

R = C.ROOT
OUTD = C.OUT / "tiles_gnn"
OUTD.mkdir(exist_ok=True)


def main():
    from fuelblend.api.predict_route import _tmap_freezing
    from fuelblend.baselines import flash_unifac as fu
    from fuelblend.data.tile_isomers import get_registry_and_reps, tile_mean_tmap
    tdes = _tmap_freezing()
    reg, reps = get_registry_and_reps()
    pool = pd.read_csv(R / "results/v6/gnn/pool.csv")
    meas = {t: pool[pool.target == t].set_index("inchikey").value.to_dict() for t in ("tm_K", "dhfus_kJ_mol")}
    ck = torch.load(R / "results/v6/mixture/encoder.pt", map_location="cpu", weights_only=False)
    enc = MP.load_encoder(R / "results/v6/mixture/encoder.pt")
    jt, jd = ck["tasks"].index("tm_K"), ck["tasks"].index("dhfus_kJ_mol")
    smi = dict(zip(reg.ik, reg.smiles))

    def has(t, ik, key):
        return ik in t and np.isfinite(t[ik].get(key, np.nan))

    need = [ik for ik in reg.ik if not (has(tdes, ik, "tm_K") and has(tdes, ik, "dh_fus_J_mol"))]
    b = D.collate([D.featurize(smi[k]) for k in need])
    with torch.no_grad():
        z = torch.stack([m(b) for m in enc]).mean(0).cpu().numpy()
    g_tm, g_dh = dict(zip(need, z[:, jt] * ck["sd"][jt] + ck["mu"][jt])), dict(zip(need, z[:, jd] * ck["sd"][jd] + ck["mu"][jd]))
    tr_ik = list(meas["tm_K"])
    tr_smi = pd.read_csv(R / "results/v6/gnn/pool.csv").drop_duplicates("inchikey").set_index("inchikey").smiles
    f_tm = dict(zip(need, np.asarray(F.fit_predict(F.features([tr_smi[k] for k in tr_ik]), np.array([meas["tm_K"][k] for k in tr_ik]),
                                                   F.features([smi[k] for k in need])))))
    counts = {}

    def build(tm_est):
        t = dict(tdes)
        n = dict(designated=0, pool_measured=0, estimated=0)
        for ik in reg.ik:
            e = dict(t.get(ik) or {})
            for key, mk, est, scale in (("tm_K", "tm_K", tm_est, 1.0), ("dh_fus_J_mol", "dhfus_kJ_mol", g_dh, 1000.0)):
                if has(tdes, ik, key):
                    n["designated"] += 1
                elif ik in meas[mk]:
                    e[key] = float(meas[mk][ik]) * scale; n["pool_measured"] += 1
                else:
                    e[key] = float(est[ik]) * scale; n["estimated"] += 1
            if e.get("dh_fus_J_mol", 1) <= 0:
                e["dh_fus_J_mol"] = np.nan
            t[ik] = e
        counts.update({f"{id(tm_est)}": n})
        return t, n

    ds = pd.read_csv(C.OUT / "features_preds.csv", dtype={"fuel_id": str})
    ents = C.entities(ds)
    s5 = C.load("s5_v6c", "scripts/v5/s5_freezing_nalkane.py")
    g1 = C.load("g1_v6c", "scripts/v5/eval_ajftd_gcxgc_v5.py")
    dev = ds[ds.role == "dev"].reset_index(drop=True)
    y = dev.freezing_C.to_numpy()
    oof = pd.read_csv(C.OUT / "dev_oof_predictions.csv")
    const = oof["a_constant"].to_numpy()
    rows, info = [], {}
    for name, est in (("e2_tiles_gnn", g_tm), ("e3_tiles_fusion_tm", f_tm)):
        t, n = build(est)
        tm_t = tile_mean_tmap(reg, reps, t)
        changed = sum(1 for k in tm_t if tm_t[k] is not t.get(k) and tm_t[k].get("n_isomers"))
        p = np.array([C.unifac_designated(ents[e], s5, g1, CR, fu, tm_t) for e in dev.entity])
        info[name] = dict(fills=n, tiles_changed=changed)
        dev[f"pred_{name}"] = p
        for tag, mask in (("all", np.ones(len(y), bool)), ("no_fuel64", (dev.fuel_id != "64").to_numpy())):
            ok = mask & np.isfinite(p)
            rows.append(dict(model=name, subset=tag, n_nonfinite=int((mask & ~np.isfinite(p)).sum()), **S3.stats(y[ok], p[ok], const[ok])))
    for name, col in (("b_unifac_svl", "pred_b"), ("e_isomer_tiles_v6", "pred_e")):
        for tag, mask in (("all", np.ones(len(y), bool)), ("no_fuel64", (dev.fuel_id != "64").to_numpy())):
            rows.append(dict(model=name, subset=tag, n_nonfinite=0, **S3.stats(y[mask], dev[col].to_numpy(float)[mask], const[mask])))
    r = pd.DataFrame(rows)
    r.to_csv(OUTD / "dev_summary.csv", index=False, float_format="%.4f")
    dev[["entity", "fuel_id", "freezing_C", "pred_b", "pred_e", "pred_e2_tiles_gnn", "pred_e3_tiles_fusion_tm"]].to_csv(OUTD / "dev_preds.csv", index=False, float_format="%.4f")
    json.dump(info, open(OUTD / "fill_info.json", "w"), indent=1)
    pd.set_option("display.width", 200)
    print(json.dumps(info)); print(r[["model", "subset", "n", "n_nonfinite", "mae", "bias", "frac_under", "ci_lo", "ci_hi"]].round(2).to_string())


if __name__ == "__main__":
    main()
