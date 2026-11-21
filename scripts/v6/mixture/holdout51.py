"""Track C replication test on the 51 Leenhouts rows dropped for low pure-MON coverage (declared in results/v6/mixture/README.md)."""
import sys, json, importlib.util
from pathlib import Path
import numpy as np, pandas as pd, torch
R = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(R / "src")); sys.path.insert(0, str(R / "scripts/external_ingest")); sys.path.insert(0, str(R / "scripts/v6/mixture"))
sys.argv = ["run_cv", "--stage", "dev"]
import run_cv as RC                       # noqa: E402  (reuses the v1 comparators and candidate unchanged)
from fuelblend.ml import molpool as MP    # noqa: E402
from fuelblend.ml import dmpnn as D       # noqa: E402

OUT = R / "results/v6/mixture"
STG = R / "data/interim/external_staging"


def mon_head(models_ck, smiles):
    ck = models_ck
    b = D.collate([D.featurize(s) for s in smiles])
    j = ck["tasks"].index("mon")
    return None


def main():
    ck = torch.load(OUT / "encoder.pt", map_location="cpu", weights_only=False)
    enc = MP.load_encoder(OUT / "encoder.pt")
    j = ck["tasks"].index("mon")

    def pred_mon(smiles):
        b = D.collate([D.featurize(s) for s in smiles])
        with torch.no_grad():
            z = torch.stack([m(b) for m in enc]).mean(0).numpy()[:, j]
        return z * ck["sd"][j] + ck["mu"][j]

    pool = pd.read_csv(R / "results/v6/gnn/pool.csv")
    pure = pool[pool.target == "mon"].set_index("inchikey").value.to_dict()
    d = pd.read_csv(STG / "canonical_leenhouts_mon_mix_measurements.csv")
    cmp_ = pd.read_csv(STG / "canonical_leenhouts_mon_compounds.csv")
    smi = dict(zip(cmp_.inchikey, cmp_.smiles))
    missing = sorted({k for s in d.comp_inchikeys for k in json.loads(s) if k not in pure})
    filled = dict(zip(missing, pred_mon([smi[k] for k in missing])))
    full_pure = {**filled, **pure}
    recs = []
    for r in d.itertuples():
        ks = json.loads(r.comp_inchikeys); x = np.array(json.loads(r.x_mole), float); x /= x.sum()
        recs.append((r.system_id, tuple(ks), tuple(np.round(x, 3)), float(r.value)))
    df = pd.DataFrame(recs, columns=["system_id", "ks", "xr", "y"]).groupby(["system_id", "ks", "xr"], as_index=False).y.median()
    rows = []
    for r in df.itertuples():
        ks, x = list(r.ks), np.array(r.xr, float); x /= x.sum()
        cov = float(sum(xi for k, xi in zip(ks, x) if k in pure))
        if cov >= 0.9:
            continue
        rows.append(dict(head="leenhouts", prop="mon", group=r.system_id, system_id=f"leenhouts:{r.system_id}", n_comp=len(ks), T_K=np.nan,
                         y=r.y, base=float(sum(xi * full_pure[k] for k, xi in zip(ks, x))), cov=cov, comps=json.dumps(ks),
                         smiles=json.dumps([smi[k] for k in ks]), x=json.dumps(list(x)), w="[]", role="holdout51"))
    ho = pd.DataFrame(rows)
    assert len(ho) == 51, len(ho)
    tr_df = pd.read_csv(OUT / "rows.csv")
    tr_df = tr_df[tr_df["head"] == "leenhouts"]
    assert not set(tr_df.group) & set(ho.group), "holdout systems overlap training systems"
    full = pd.concat([tr_df, ho], ignore_index=True)
    h = RC.Head(full, enc)
    tr, te = np.arange(len(tr_df)), np.arange(len(tr_df), len(full))
    rec = []
    preds = RC.fit_all(h, tr, te, "octane", full_pure, rec, include_te=True)
    out = ho.copy()
    for k, v in preds.items():
        out[f"pred_{k}"] = v
    out.to_csv(OUT / "holdout51_pred.csv", index=False)
    rng = np.random.default_rng(0)
    g = out.group.to_numpy()
    res = dict(n=len(out), systems=out.group.nunique(), median_cov=float(out["cov"].median()), chosen_cfg=[str(RC.CONFIGS[i]) for i in rec])
    for k in ["lin_molar", "corrected_pure", "desc_gbm", "molpool"]:
        res[f"mae_{k}"] = float(np.abs(out[f"pred_{k}"] - out.y).mean())
    for k in ["lin_molar", "corrected_pure", "desc_gbm"]:
        ec, em = np.abs(out[f"pred_{k}"] - out.y).to_numpy(), np.abs(out.pred_molpool - out.y).to_numpy()
        lo, hi = RC.boot(ec - em, g, rng)
        res[f"gain_vs_{k}"] = [float(lo), float(hi)]
    mc = (out.n_comp >= 3).to_numpy()
    res["n_multi"] = int(mc.sum())
    if mc.any():
        res["mae_multi_desc_gbm"] = float(np.abs(out.pred_desc_gbm - out.y)[mc].mean())
        res["mae_multi_molpool"] = float(np.abs(out.pred_molpool - out.y)[mc].mean())
    res["gate_pass"] = bool(res["gain_vs_desc_gbm"][0] > 0 and (not mc.any() or res["mae_multi_molpool"] <= res["mae_multi_desc_gbm"]))
    json.dump(res, open(OUT / "holdout51_readings.json", "w"), indent=1)
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
