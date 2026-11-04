#!/usr/bin/env python
"""T1: strict leave-family-out CV for the pure flash model.

    python scripts/t1_leakage_audit.py

Re-uses the WP4 pipeline (collect/prepare/train_pure_models helpers) but swaps the
molecule-grouped holdout for scaffold-FAMILY holdout (5 folds, seed 20261030).
Verdict written to reports/T1_LEAKAGE.md.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

BASE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(BASE / "src"))
sys.path.insert(0, str(BATSE := BASE / "scripts"))

import train_pure_models as tpm  # noqa: E402  (helpers collect/prepare/evaluate reuse WP4 exactly)
from fuelblend.features.family import scaffold_family_key, leave_family_folds  # noqa: E402
from fuelblend.models import pure_gbm as pgm  # noqa: E402
from fuelblend.eval import metrics as evm  # noqa: E402

PROP = "flash_point_C"
N_FOLDS = 5
SEED = 20261030
LOG = []


def log(msg) -> None:
    print(msg, flush=True)
    LOG.append(str(msg))


def main() -> None:
    log("T1 leakage audit - strict family holdout re-check of the pure flash claim")
    log(f"property={PROP} n_folds={N_FOLDS} seed={SEED} (WP4 loose holdback claim: MAE 0.5615 RMSE 1.5927, n=142)")
    log("")

    pool = tpm.collect(PROP, sorted(tpm.S1 | tpm.S2))
    for r in pool:
        r["family"] = scaffold_family_key(r["smiles"])
    fams = sorted({r["family"] for r in pool})
    fam_n = {f: sum(1 for r in pool if r["family"] == f) for f in fams}
    log(f"pool rows={len(pool)} molecules={len({r['inchikey'] for r in pool})} families={len(fams)}")
    for f in sorted(fam_n.items(), key=lambda t: -t[1])[:15]:
        log(f"  family {f[0][:40]:42s} rows={f[1]}")

    folds = leave_family_folds(pool, n_folds=N_FOLDS, seed=SEED)
    all_true, all_pred, all_cl = [], [], []
    per_fold = []
    for k, fold in enumerate(folds):
        tr, te = fold["train"], fold["test"]
        tr2, X_tr, y_tr, _, sw_tr = tpm.prepare(PROP, tr)
        te2, X_te, y_te, _, _ = tpm.prepare(PROP, te)
        if not tr2 or not te2:
            log(f"fold {k}: EMPTY after feature prep (train={len(tr2)} test={len(te2)}) - skipped")
            continue
        preds = []
        for seed in tpm.MODEL_SEEDS:
            model = pgm.make_model(random_state=seed)
            pgm.fit_pure_gbm(model, X_tr, y_tr, sample_weight=sw_tr)
            preds.append(pgm.predict_pure_gbm(model, X_te))
        mean_p = np.mean(preds, axis=0)
        cl = np.array([r["system_id"] for r in te2])
        st = tpm.evaluate(PROP, y_te, mean_p, cl)
        per_fold.append(dict(n=st["n"], families=[f.split(":")[-1][:18] for f in fold["families"]],
                             mae=st["mae"], rmse=st["rmse"]))
        all_true.append(np.exp(y_te) if PROP == "viscosity_dyn_mPas" else y_te)
        so_far = np.concatenate(all_true)
        log(f"fold {k}: n={st['n']} mae={st['mae']:.4f} rmse={st['rmse']:.4f} "
            f"test_families={','.join(f.split(':')[-1][:18] for f in fold['families'])}")
        # keep predictions in the same linearized space for a pooled number
        all_true[-1] = y_te
        all_pred.append(mean_p)
        all_cl.append(cl)

    if all_pred:
        yt = np.concatenate(all_true)
        yp = np.concatenate(all_pred)
        cl = np.concatenate(all_cl)
        st = tpm.evaluate(PROP, yt, yp, cl)
        ci = evm.bootstrap_ci(yt, yp, cl, metric_fn=evm.rmse, n_boot=2000, seed=SEED)
        log("")
        log(f"POOLED family-holdout: n={st['n']} mae={st['mae']:.4f} rmse={st['rmse']:.4f} "
            f"bias={st.get('bias', float('nan')):.4f} "
            f"overpred_rate_5={st.get('overpred_rate_5', float('nan')):.4f} "
            f"mean_pos_err={st.get('mean_positive_error', float('nan')):.4f}")
        log(f"bootstrap 95% CI on RMSE: [{ci['lo']:.4f}, {ci['hi']:.4f}]")
        # MAE CI seen as well
        ci_mae = evm.bootstrap_ci(yt, yp, cl, metric_fn=evm.mae, n_boot=2000, seed=SEED)
        log(f"bootstrap 95% CI on MAE: [{ci_mae['lo']:.4f}, {ci_mae['hi']:.4f}]")

    report = BASE / "reports/T1_LEAKAGE.md"
    lines = LOG + ["", "## Verdict", ""]
    lines.append("If pooled family-holdout RMSE is materially above the WP4 heldback RMSE (1.5927), "
                 "the 0.56 C claim reflects molecular interpolation (sibling molecules in train), "
                 "not family generalization; quantify the gap and say so in plain words.")
    report.parent.mkdir(parents=True, exist_ok=True)
    report.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"\nwrote {report}")


if __name__ == "__main__":
    main()
