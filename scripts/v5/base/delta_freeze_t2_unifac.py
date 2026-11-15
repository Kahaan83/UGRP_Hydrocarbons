"""v5 delta_freeze_t2_unifac: score the designated freezing route (UNIFAC-SvL) on the T2 real-fuel freezing set.

Declared in results/v5/BASE/delta_freeze_t2_unifac/README.md before this script was first run.  Same fuels, rows, comparator
(S5 column C = logo_constant) and verdict rule (cellkit.make_row, exactly as cells_t2_bulk.freezing_cells) as the scoreboard's
ideal-SvL freezing cells.  Reads only; writes rows.csv, summary.csv, readings.json in the delta folder.
"""
from __future__ import annotations

import importlib.util
import json
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "src"))
from fuelblend.eval import composition_route as CR  # noqa: E402
from fuelblend.eval.cellkit import make_row  # noqa: E402

RDLogger.DisableLog("rdApp.*")
OUT = ROOT / "results/v5/BASE/delta_freeze_t2_unifac"
S5P = ROOT / "results/v5/S5/s5_predictions.csv"
RERUN = ROOT / "results/v5/BASE/delta_adopt_rescore/rerun/cells_t2_bulk.csv"
MARGIN_UNIFAC = 13.84274050383641      # configs/uncertainty_v5.yaml freezing_point_C.lab_liquidus.svl_unifac
MARGIN_IDEAL = 16.975450975521966      # configs/uncertainty_v5.yaml freezing_point_C.lab_liquidus.svl_verified
COLS = {"I_s5": ("ideal_svl", "s5_constants_repro"), "unifac_s5": ("svl_unifac", "s5_constants"),
        "ideal_desig": ("ideal_svl", "designated_constants"), "unifac_desig": ("svl_unifac", "designated")}
TIER, EV = "T2_measured", "results/v5/BASE/delta_freeze_t2_unifac/rows.csv"


def load(name: str, rel: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / rel)
    m = importlib.util.module_from_spec(spec)
    sys.modules[name] = m
    spec.loader.exec_module(m)
    return m


def s5_inputs(ent: dict, g1):
    """x, iks, smiles and the S5 constants map, built exactly as s5_freezing_nalkane.ideal_svl_I."""
    mp = ent["mapped"]
    mols = {r.ik: Chem.MolFromSmiles(r.smiles) for r in mp.itertuples()}
    info = pd.DataFrame(dict(smiles=[r.smiles for r in mp.itertuples()],
                             cls=[g1.hc_class(mols[r.ik]) for r in mp.itertuples()],
                             n_carbon=[float(sum(a.GetSymbol() == "C" for a in mols[r.ik].GetAtoms())) for r in mp.itertuples()]),
                        index=list(mp.ik)).loc[lambda d: ~d.index.duplicated()]
    tmap = g1.get_pure_thermal_properties(info)
    w = mp.frac.to_numpy(float)
    x = w / np.array([Descriptors.MolWt(mols[k]) for k in mp.ik])
    return x / x.sum(), list(mp.ik), list(mp.smiles), tmap


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    s5 = load("s5_t2u", "scripts/v5/s5_freezing_nalkane.py")
    g1 = load("g1_t2u", "scripts/v5/eval_ajftd_gcxgc_v5.py")
    from fuelblend.api.predict_route import _tmap_freezing
    from fuelblend.baselines import flash_unifac as fu
    tdes = _tmap_freezing()
    ents = s5.build_entities()
    p = pd.read_csv(S5P)
    rows = []
    for r in p.itertuples():
        x, iks, smi, ts5 = s5_inputs(ents[r.entity], g1)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            gam = fu.make_unifac_gamma_fn_vec(smi)
            v = {}
            v["I_s5"], _ = CR.freezing_C(x, iks, smi, ts5, model="ideal_svl")
            v["unifac_s5"], i1 = CR.freezing_C(x, iks, smi, ts5, model="svl_unifac", gamma=gam)
            v["ideal_desig"], _ = CR.freezing_C(x, iks, smi, tdes, model="ideal_svl")
            v["unifac_desig"], i2 = CR.freezing_C(x, iks, smi, tdes, model="svl_unifac", gamma=gam)
        fb = i2["unifac_fallback"]
        miss = [k for k in iks if not (tdes.get(k) or {}).get("tm_K")]
        nver = sum((tdes.get(k) or {}).get("source") == CR.VERIFIED_SOURCE for k in iks)
        rows.append(dict(entity=r.entity, source=r.source, group=r.group, role=r.role, freezing_C=r.freezing_C, C=r.C, I_scoreboard=r.I,
                         **{k: (float(val) if np.isfinite(val) else np.nan) for k, val in v.items()},
                         n_comp=len(iks), unifac_fallback_n=len(fb), unifac_fallback_x=float(sum(x[i] for i in fb)),
                         n_missing_designated_constants=len(miss), x_missing_designated_constants=float(sum(x[i] for i, k in enumerate(iks) if k in miss)),
                         n_verified_constants=int(nver)))
    d = pd.DataFrame(rows)
    d["I_repro_diff"] = d.I_s5 - d.I_scoreboard
    d.to_csv(OUT / "rows.csv", index=False, float_format="%.8g")

    # sanity check: row level
    both = d.I_s5.notna() & d.I_scoreboard.notna()
    nan_mismatch = int((d.I_s5.isna() != d.I_scoreboard.isna()).sum())
    max_diff = float(d.loc[both, "I_repro_diff"].abs().max())
    repro_ok = nan_mismatch == 0 and max_diff <= 1e-6

    # cells, exactly the freezing_cells() loop
    cells = []
    fin_all = d[list(COLS)].notna().all(axis=1) & d.C.notna() & d.freezing_C.notna()
    for role in ("headline", "locked"):
        rr0 = d[d.role == role]
        for sname, sv in (("pooled", None), ("AV-33-22", "AV-33-22"), ("AJFTD", "AJFTD"), ("FuelLib", "FuelLib")):
            rr = rr0 if sv is None else rr0[rr0.source == sv]
            for col, (model, var) in COLS.items():
                m = rr.dropna(subset=[col, "C", "freezing_C"])
                if m.empty:
                    continue
                row = make_row("freezing_point_C", TIER, sname + ("" if role == "headline" else "_locked"), model, var, "recomputed_v5",
                               "descriptive" if role == "headline" else "already_read", EV,
                               err=m[col] - m.freezing_C, groups=m.group, err_trivial=m.C - m.freezing_C, trivial_model="logo_constant",
                               n_excluded_nonfinite=int(len(rr) - len(m)))
                mg = MARGIN_UNIFAC if model == "svl_unifac" else MARGIN_IDEAL
                row["column"] = col
                row["margin_used"] = mg
                row["coverage_margin"] = float((m.freezing_C <= m[col] + mg).mean())
                row["coverage_unifac_margin_13p84"] = float((m.freezing_C <= m[col] + MARGIN_UNIFAC).mean())
                row["frac_under_predicted"] = float((m[col] < m.freezing_C).mean())
                mi = rr[fin_all.loc[rr.index]]
                e = mi[col] - mi.freezing_C
                row.update(n_intersection=len(mi), bias_intersection=float(e.mean()) if len(mi) else np.nan,
                           mae_intersection=float(e.abs().mean()) if len(mi) else np.nan,
                           mae_trivial_intersection=float((mi.C - mi.freezing_C).abs().mean()) if len(mi) else np.nan)
                cells.append(row)
    s = pd.DataFrame(cells)
    s.to_csv(OUT / "summary.csv", index=False, float_format="%.8g")

    # sanity check: cell level vs the rerun scoreboard cell
    rer = pd.read_csv(RERUN)
    rer = rer[(rer.property == "freezing_point_C") & (rer.model == "ideal_svl")].set_index("set")
    cell_chk = {}
    for sname in ("pooled", "pooled_locked", "AV-33-22", "AJFTD", "FuelLib", "AJFTD_locked"):
        a = s[(s.set == sname) & (s.column == "I_s5")]
        if a.empty or sname not in rer.index:
            continue
        a, b = a.iloc[0], rer.loc[sname]
        diffs = {k: abs(float(a[k]) - float(b[k])) for k in ("n_groups", "bias", "mae", "mae_trivial", "ci_lo", "ci_hi")}
        cell_chk[sname] = dict(max_abs_diff=max(diffs.values()), verdict_match=bool(a.verdict == b.verdict))
    cells_ok = all(v["max_abs_diff"] <= 1e-4 and v["verdict_match"] for v in cell_chk.values()) and "pooled" in cell_chk

    def cell(sname, col):
        r = s[(s.set == sname) & (s.column == col)]
        if r.empty:
            return None
        r = r.iloc[0]
        return {k: (r[k].item() if hasattr(r[k], "item") else r[k]) for k in
                ("n_groups", "n_rows", "n_excluded_nonfinite", "bias", "mae", "rmse", "mae_trivial", "diff_trivial_minus_model", "ci_lo", "ci_hi",
                 "verdict", "coverage_margin", "margin_used", "coverage_unifac_margin_13p84", "frac_under_predicted", "n_intersection",
                 "bias_intersection", "mae_intersection", "mae_trivial_intersection")}

    readings = dict(
        sanity=dict(row_repro_ok=repro_ok, row_max_abs_diff=max_diff, row_nan_mismatch=nan_mismatch, cell_check=cell_chk, cell_repro_ok=cells_ok),
        primary={sn: cell(sn, "unifac_desig") for sn in ("pooled", "pooled_locked", "AV-33-22", "AJFTD", "FuelLib")},
        ideal_designated={sn: cell(sn, "ideal_desig") for sn in ("pooled", "pooled_locked", "AV-33-22", "AJFTD", "FuelLib")},
        unifac_s5={sn: cell(sn, "unifac_s5") for sn in ("pooled", "pooled_locked")},
        ideal_s5={sn: cell(sn, "I_s5") for sn in ("pooled", "pooled_locked")},
        unifac_fallback=dict(rows_with_fallback=int((d.unifac_fallback_n > 0).sum()), max_fallback_x=float(d.unifac_fallback_x.max())),
        designated_constants=dict(rows_with_missing=int((d.n_missing_designated_constants > 0).sum()),
                                  max_missing_x=float(d.x_missing_designated_constants.max()),
                                  rows_all_verified=int((d.n_verified_constants == d.n_comp).sum())),
        margins=dict(svl_unifac=MARGIN_UNIFAC, svl_verified=MARGIN_IDEAL),
    )
    (OUT / "readings.json").write_text(json.dumps(readings, indent=2, default=float))
    pd.set_option("display.width", 250, "display.max_columns", 30)
    print("row repro ok", repro_ok, max_diff, nan_mismatch, "| cell repro ok", cells_ok, cell_chk)
    print(s[["set", "column", "n_groups", "n_excluded_nonfinite", "bias", "mae", "mae_trivial", "ci_lo", "ci_hi", "verdict", "coverage_margin",
             "frac_under_predicted", "n_intersection", "mae_intersection"]].round(4).to_string())
    print(json.dumps({k: readings[k] for k in ("unifac_fallback", "designated_constants")}))


if __name__ == "__main__":
    main()
