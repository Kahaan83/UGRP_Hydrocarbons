"""v5 BASE Task 16b: frozen ideal lab flash errors on the Leenhouts benchmark's own split definitions (skip-guarded; zips read in memory).

Outputs (results/v5/BASE/delta_leenhouts_splits/): fold_metrics.csv, delta.csv, readings.json
"""
from __future__ import annotations

import io
import json
import sys
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "src"))
from fuelblend.eval.cellkit import frame, make_row  # noqa: E402

Z = ROOT / "extras/fuelblend_data/data/raw/zenodo/leenhouts_benchmark"
OUT = ROOT / "results" / "v5" / "BASE" / "delta_leenhouts_splits"
SPLITS = ("mixture", "molecule", "pure_to_mixture", "random_rowwise")
RDLogger.DisableLog("rdApp.*")


def key(smiles, x, value) -> tuple:
    return (tuple(sorted((s, round(float(f), 3)) for s, f in zip(smiles, x))), round(float(value), 1))


def main() -> None:
    if not (Z / "splits.zip").exists() or not (Z / "datasets.zip").exists():
        print("SKIPPED: Leenhouts zips not in extras/")
        return
    OUT.mkdir(parents=True, exist_ok=True)
    dz, sz = zipfile.ZipFile(Z / "datasets.zip"), zipfile.ZipFile(Z / "splits.zip")
    mix = pd.read_csv(io.BytesIO(dz.read("datasets/flashpoint/processed_data/flashpoint_mix_exp.csv")))
    can = {}
    keys = []
    for r in mix.itertuples():
        sm, xs = [], []
        for j in range(4):
            inchi = getattr(r, f"component_inchi_{j}")
            if isinstance(inchi, str):
                if inchi not in can:
                    m = Chem.MolFromInchi(inchi)
                    can[inchi] = Chem.MolToSmiles(m) if m is not None else None
                sm.append(can[inchi])
                xs.append(getattr(r, f"fraction_{j}"))
        keys.append(key(sm, xs, r.value) if None not in sm else None)
    lab = pd.read_csv(ROOT / "results/v4/F/flash_lab_loso.csv")
    lab["key"] = [key([Chem.MolToSmiles(Chem.MolFromSmiles(s)) for s in json.loads(a)], json.loads(b), v)
                  for a, b, v in zip(lab.comp_smiles, lab.x_mole, lab.measured_C)]
    kmap = {k: i for i, k in enumerate(keys) if k is not None}
    lab["mix_row"] = [kmap.get(k, -1) for k in lab.key]
    rows, cells = [], []
    for sp in SPLITS:
        for f in range(5):
            fd = pd.read_csv(io.BytesIO(sz.read(f"splits/flashpoint/flashpoint_mix_exp/split_definitions/{sp}/cross_validation_5fold/fold_{f:02d}.csv")))
            test = set(fd[(fd.split == "test") & (fd.source == "mix")].source_row_index.astype(int))
            g = lab[lab.mix_row.isin(test)]
            e = g.ideal_err_C
            rows.append(dict(split=sp, fold=f, n_benchmark_test_rows=len(test), n_ours=len(g), n_systems=g.system_id.nunique(),
                             mae=float(e.abs().mean()) if len(g) else np.nan, rmse=float(np.sqrt((e ** 2).mean())) if len(g) else np.nan,
                             bias=float(e.mean()) if len(g) else np.nan))
            if len(g):
                cells.append(make_row("flash_point_C", "T1_lab", f"leenhouts_{sp}_fold{f:02d}", "ideal_le_chatelier", "designated", "aggregated_v5", "dev",
                                      "results/v4/F/flash_lab_loso.csv; leenhouts splits.zip", err=e, groups=g.system_id, unsafe=lambda x: x > 5.0,
                                      extra=dict(baseline_mae=np.nan, baseline_bias=np.nan, baseline_unsafe_rate=np.nan)))
    fm = pd.DataFrame(rows)
    fm.to_csv(OUT / "fold_metrics.csv", index=False, float_format="%.4f")
    frame(cells).to_csv(OUT / "delta.csv", index=False)
    res = dict(n_lab_rows=len(lab), n_matched=int((lab.mix_row >= 0).sum()),
               by_split={sp: dict(mae_mean_over_folds=float(fm[fm.split == sp].mae.mean()), rmse_mean_over_folds=float(fm[fm.split == sp].rmse.mean()),
                                  rows_covered=int(fm[fm.split == sp].n_ours.sum())) for sp in SPLITS},
               all_rows_mae=float(lab.ideal_err_C.abs().mean()))
    (OUT / "readings.json").write_text(json.dumps(res, indent=2))
    print(fm.to_string())
    print(json.dumps(res, indent=2))


if __name__ == "__main__":
    main()
