"""v4 Task 19: effect of two suspect anchored-Psat compounds on the lab flash set.

mesitylene (anchored Psat 0.52x the chemicals reference) and n-hexylbenzene (1.20x, status disputed) are flagged in
results/v4/F/pure_inputs_used.csv (column anchored_over_ref). This script recomputes lab-set MAE / bias / worst over-prediction
and the Task 10 style 90 % system-level margin on results/v4/F/flash_lab_loso.csv with and without the records whose
composition contains either compound. No model is rerun and nothing is fitted; it only filters existing rows.
Output: results/v4/F/flash_input_sensitivity.csv
"""
from __future__ import annotations
import json, math
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "v4" / "F"
SUSPECT = {"mesitylene": "Cc1cc(C)cc(C)c1", "n-hexylbenzene": "CCCCCCc1ccccc1"}


def _margin(df: pd.DataFrame, col: str, alpha: float = 0.1) -> tuple[float, int]:
    g = df.groupby("system_id")[col].max().sort_values().to_numpy()
    n = len(g)
    k = math.ceil((n + 1) * (1 - alpha))
    return (float("inf") if k > n else float(g[k - 1])), n


def main() -> None:
    lab = pd.read_csv(OUT / "flash_lab_loso.csv")
    has = {name: lab.comp_smiles.apply(lambda s, smi=smi: smi in json.loads(s)) for name, smi in SUSPECT.items()}
    any_s = has["mesitylene"] | has["n-hexylbenzene"]
    rows = []
    for label, mask in [("all", pd.Series(True, index=lab.index)),
                        ("excl_mesitylene", ~has["mesitylene"]),
                        ("excl_n-hexylbenzene", ~has["n-hexylbenzene"]),
                        ("excl_both", ~any_s),
                        ("only_mesitylene_records", has["mesitylene"]),
                        ("only_n-hexylbenzene_records", has["n-hexylbenzene"])]:
        d = lab[mask]
        for model, col in [("ideal", "ideal_err_C"), ("unifac", "unifac_err_C")]:
            if len(d) == 0:
                continue
            m, ng = _margin(d, col)
            rows.append(dict(subset=label, model=model, n_records=len(d), n_systems=d.system_id.nunique(),
                             MAE_C=d[col].abs().mean(), bias_C=d[col].mean(), worst_over_C=d[col].max(),
                             margin90_C=m, margin_n_groups=ng))
    pd.DataFrame(rows).to_csv(OUT / "flash_input_sensitivity.csv", index=False, float_format="%.6g")
    print(pd.DataFrame(rows).to_string())


if __name__ == "__main__":
    main()
