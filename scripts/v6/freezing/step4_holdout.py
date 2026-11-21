"""Track A step 4: score the locked holdout ONCE.  Models are refit on all 52 dev groups with the settings declared in README.md.
Margins come from dev out-of-fold residuals (dev_oof_predictions.csv).  Refuses to run twice (writes holdout_scored.flag)."""
from __future__ import annotations

import datetime as dt
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common as C  # noqa: E402
from step3_models import fp_const, fp_control, fp_ridge, make_lgbm, stats  # noqa: E402

FLAG = C.OUT / "holdout_scored.flag"


def main() -> None:
    if FLAG.exists():
        raise SystemExit("holdout already scored once; refusing to score it again")
    ds = pd.read_csv(C.OUT / "features_preds.csv", dtype={"fuel_id": str})
    dev, ho = ds[ds.role == "dev"].reset_index(drop=True), ds[ds.role == "holdout"].reset_index(drop=True)
    base = json.load(open(C.OUT / "chosen_base.json"))
    oof = pd.read_csv(C.OUT / "dev_oof_predictions.csv")
    lg_name = "f_lgbm_residual_on_" + base["base"]
    preds = {"a_constant": fp_const(dev, ho), "c2_ridge": fp_ridge(dev, ho), "d_control_curve": fp_control(dev, ho),
             "b_unifac_svl": ho.pred_b.to_numpy(), "c1_weighted_tm": ho.pred_c1.to_numpy(), "e_isomer_tiles": ho.pred_e.to_numpy(),
             lg_name: make_lgbm(base["base_col"])(dev, ho)}
    y = ho.freezing_C.to_numpy()
    rows = []
    for name, p in preds.items():
        s = stats(y, p, preds["a_constant"])
        mg = float(np.quantile(oof.y - oof[name], 0.9))               # dev out-of-fold margin, applied upward
        s.update(model=name, margin90_dev_oof=mg, coverage_with_dev_margin=float((y <= p + mg).mean()),
                 mae_excl_325=float(np.abs(p - y)[ho.fuel_id != "325"].mean()))
        rows.append(s)
    r = pd.DataFrame(rows)
    r.to_csv(C.OUT / "holdout_summary.csv", index=False, float_format="%.4f")
    ho.assign(**{k: v for k, v in preds.items()}).to_csv(C.OUT / "holdout_predictions.csv", index=False, float_format="%.4f")
    FLAG.write_text(dt.datetime.now(dt.timezone.utc).isoformat() + "\n")
    with open(C.ROOT / "reports/LOCKED_ACCESS.log", "a") as f:
        f.write(f"{dt.datetime.now(dt.timezone.utc).isoformat()} v6-freezing holdout scored once (12 groups) by step4_holdout.py\n")
    pd.set_option("display.width", 200)
    print(r[["model", "n", "mae", "bias", "mae_const", "ci_lo", "ci_hi", "margin90_dev_oof", "coverage_with_dev_margin", "mae_excl_325", "verdict"]].round(2).to_string())


if __name__ == "__main__":
    main()
