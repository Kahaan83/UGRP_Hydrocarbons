"""Track A step 2: baselines (a) constant and (b) v5 UNIFAC-SvL, plus the reproduction gate on the 50 v5 headline rows.
The gate passes only if baseline (b) recomputed here equals v5 rows.csv unifac_desig (max abs diff <= 1e-5, the %.8g precision v5 rows.csv was stored at) and
the v5 headline cell (n 50, MAE 8.4246, bias +5.8778) is reproduced."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common as C  # noqa: E402

from fuelblend.eval import composition_route as CR  # noqa: E402


def main() -> None:
    ds = C.dataset()
    ents = C.entities(ds)
    s5 = C.load("s5_v6b", "scripts/v5/s5_freezing_nalkane.py")
    g1 = C.load("g1_v6b", "scripts/v5/eval_ajftd_gcxgc_v5.py")
    from fuelblend.api.predict_route import _tmap_freezing
    from fuelblend.baselines import flash_unifac as fu
    tdes = _tmap_freezing()
    ds["pred_b"] = [C.unifac_designated(ents[e], s5, g1, CR, fu, tdes) for e in ds.entity]
    ds["unmapped"] = [ents[e]["unmapped"] for e in ds.entity]
    ds.to_csv(C.OUT / "baseline_b_preds.csv", index=False, float_format="%.8g")

    v5 = pd.read_csv(C.ROOT / "results/v5/BASE/delta_freeze_t2_unifac/rows.csv").set_index("entity")
    old = ds[ds.v5_role.isin(["headline", "locked"])].set_index("entity")
    diff = (old.pred_b - v5.loc[old.index, "unifac_desig"]).abs()
    head = ds[ds.v5_role == "headline"]
    err = head.pred_b - head.freezing_C
    gate = dict(max_abs_diff=float(diff.max()), n_headline=len(head), mae=float(err.abs().mean()), bias=float(err.mean()))
    ok = gate["max_abs_diff"] <= 1e-5 and gate["n_headline"] == 50 and abs(gate["mae"] - 8.4245951) < 1e-4 and abs(gate["bias"] - 5.8778236) < 1e-4
    print(gate, "GATE", "PASS" if ok else "FAIL")
    pd.Series(gate | dict(passed=ok)).to_json(C.OUT / "repro_gate.json", indent=1)


if __name__ == "__main__":
    main()
