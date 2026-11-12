"""v5 BASE Task 12: better representatives (SubRepsV5) against the frozen scoreboard.

Only fuels whose mapped composition changes are re-predicted; all others keep their frozen predictions.
Outputs (results/v5/BASE/delta_representatives/): candidates.csv, candidate_gaps.csv, changed_fuels.csv, delta.csv, diesel_unmapped.csv, readings.json
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
from fuelblend.baselines.flash_unifac import _lookup  # noqa: E402
from fuelblend.data.reconstruction_r3 import SubReps  # noqa: E402
from fuelblend.data.representatives_v5 import SubRepsV5, candidate_table  # noqa: E402
from fuelblend.eval.cellkit import frame, make_row  # noqa: E402
from fuelblend.eval.composition_route import RouteCtx, kind_bins, map_bins, predict_properties  # noqa: E402
from fuelblend.mixture.common import CompoundRegistry, PurePropertyProvider  # noqa: E402

BASE = ROOT / "results" / "v5" / "BASE"
OUT = BASE / "delta_representatives"
V5 = ROOT / "data" / "real_fuels" / "v5"
TARGETS = [("flash_point_C", np.nan), ("viscosity_kin_cSt", 313.15), ("density_kg_m3", 288.15)]
FILL = {"antoine": "ThermoML vapour pressure on disk (refit_antoine), else NIST WebBook", "lfl": "IEC 60079-20-1 / Co-Optima table, Zabetakis 627",
        "viscosity": "ThermoML / Boehm 2022 pure viscosity on disk", "density": "ThermoML pure density on disk",
        "not in registry": "add compound to data/unified/compounds.csv with the above constants", "no SubReps subkind": "new subkind in reconstruction_r3 (code change)"}


def _mod(name, rel):
    s = importlib.util.spec_from_file_location(name, ROOT / rel)
    m = importlib.util.module_from_spec(s)
    s.loader.exec_module(m)
    return m


bulk = _mod("cells_t2_bulk_rep", "scripts/v5/base/cells_t2_bulk.py")
t2f = _mod("cells_t2_flash_rep", "scripts/v5/base/cells_t2_flash.py")


def same(a, b) -> bool:
    return a.idx.shape == b.idx.shape and np.array_equal(a.idx, b.idx) and np.allclose(a.w, b.w)


def score(d: pd.DataFrame, prop: str, T, pred: str) -> list[dict]:
    """Per source (within-source one per group) and pooled (one per group overall), as the frozen scoreboard does."""
    out = []
    d = d[~d.low_coverage.astype(bool) & (d.status == "ok") & np.isfinite(d[pred])]
    err = (lambda x: np.log(x[pred]) - np.log(x.measured)) if prop == "viscosity_kin_cSt" else (lambda x: x[pred] - x.measured)
    unsafe = (lambda e: e > 5.0) if prop == "flash_point_C" else None
    sets = [(s, g) for s, g in d.groupby("source")] + [("pooled", d)]
    for s, g in sets:
        g = g[g.chosen.astype(bool)] if prop == "flash_point_C" else bulk.one_per_group(g)
        out.append(make_row(prop, "T2_measured", s, "", "", "recomputed_v5", "descriptive", "", err=err(g), groups=g.group, unsafe=unsafe,
                            extra=dict(T_K=T)))
    return out


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    warnings.simplefilter("ignore")
    reg = CompoundRegistry(ROOT / "data/unified/compounds.csv")
    prov = PurePropertyProvider(reg, ROOT / "data/cache/pure_grid.npz", verbose=False)
    vp = lambda s: _lookup(s, None)                                                         # noqa: E731
    tab = candidate_table(reg, prov)
    tab.to_csv(OUT / "candidates.csv", index=False)
    gaps = tab[~tab.admitted.astype(bool)].copy()
    gaps["what_would_fill_it"] = [("; ".join(v for k, v in FILL.items() if k in r) or "-") for r in gaps.reason]
    gaps.to_csv(OUT / "candidate_gaps.csv", index=False)
    r0, r1 = SubReps(reg, prov, vp_lookup=vp), SubRepsV5(reg, prov, vp_lookup=vp, table=tab)
    ctxs = {"designated": RouteCtx(reg, prov), "lfl_overlay": RouteCtx(reg, prov, lfl_override=t2f.lfl_overlay_dict())}

    fu = pd.read_csv(BASE / "t2_fuels.csv")
    comps = {"av3322": pd.read_csv(V5 / "crc_av3322_composition.csv"), "fuellib": pd.read_csv(V5 / "fuellib_gcxgc_compositions.csv"),
             "ajftd": pd.read_csv(V5 / "ajftd_gcxgc_composition.csv")}
    changed = []
    for r in fu.itertuples():
        c = comps[r.source]
        fid = int(r.fuel_id) if r.source == "ajftd" else r.fuel_id
        b = kind_bins(r.source, c[c.fuel_id == fid])
        m0, m1 = map_bins(b, r0), map_bins(b, r1)
        if same(m0, m1):
            continue
        tg = pd.DataFrame([dict(property=p, T_K=t) for p, t in TARGETS])
        rec = dict(entity=r.entity)
        for tag, mf in (("old", m0), ("new", m1)):
            for q in predict_properties(ctxs["designated"], mf, tg).itertuples():
                rec[f"{tag}_{q.property}"] = q.pred
            rec[f"{tag}_flash_lfl_overlay"] = predict_properties(ctxs["lfl_overlay"], mf, tg.iloc[:1]).pred.iloc[0]
        changed.append(rec)
    cols = ["entity"] + [f"{t}_{p}" for t in ("old", "new") for p in ("flash_point_C", "viscosity_kin_cSt", "density_kg_m3", "flash_lfl_overlay")]
    ch = pd.DataFrame(changed, columns=cols)
    ch.to_csv(OUT / "changed_fuels.csv", index=False, float_format="%.6g")
    chi = ch.set_index("entity")

    rows = []
    fl = pd.read_csv(BASE / "rows_t2_flash.csv")
    for var, col, ccol in (("designated", "pred_designated", "flash_point_C"), ("lfl_overlay", "pred_lfl_overlay", "flash_lfl_overlay")):
        f = fl.copy()
        f["pred_v5"] = [chi.loc[e, f"new_{ccol}"] if e in chi.index else p for e, p in zip(f.entity, f[col])]
        for b0, b1 in zip(score(f, "flash_point_C", np.nan, col), score(f, "flash_point_C", np.nan, "pred_v5")):
            rows.append({**b1, "model": "ideal_le_chatelier", "variant": f"{var}+subreps_v5", "baseline_mae": b0["mae"], "baseline_bias": b0["bias"],
                         "baseline_unsafe_rate": b0["unsafe_rate"], "evidence_path": "rows_t2_flash.csv + changed_fuels.csv"})
    bk = pd.read_csv(BASE / "rows_t2_bulk.csv")
    for prop, T in TARGETS[1:]:
        g = bk[(bk.property == prop) & (bk.variant == "designated") & (np.round(bk.T_K, 2) == T)].copy()
        g["pred_v5"] = [chi.loc[e, f"new_{prop}"] if e in chi.index else p for e, p in zip(g.entity, g.pred)]
        for b0, b1 in zip(score(g, prop, T, "pred"), score(g, prop, T, "pred_v5")):
            rows.append({**b1, "model": "refutas_raw" if prop.startswith("visc") else "ideal_volume", "variant": "designated+subreps_v5",
                         "baseline_mae": b0["mae"], "baseline_bias": b0["bias"], "baseline_unsafe_rate": b0["unsafe_rate"],
                         "evidence_path": "rows_t2_bulk.csv + changed_fuels.csv"})
    d = frame(rows)
    d.to_csv(OUT / "delta.csv", index=False)

    g1 = _mod("g1_rep", "scripts/v5/eval_ajftd_gcxgc_v5.py")
    dc = pd.read_csv(V5 / "diesel_face_composition.csv")
    du = []
    for fid, c in dc.groupby("fuel_id"):
        b = c.groupby(["hc_class", "carbon_number"], as_index=False).mass_pct.sum().rename(columns={"carbon_number": "n_carbon"})
        du.append(dict(fuel_id=fid, unmapped_subreps=g1.map_fuel(b, r0)[1], unmapped_subreps_v5=g1.map_fuel(b, r1)[1]))
    du = pd.DataFrame(du)
    du.to_csv(OUT / "diesel_unmapped.csv", index=False, float_format="%.5g")

    pooled = d[d.set == "pooled"].set_index("property")
    dv = float(abs(pooled.loc["viscosity_kin_cSt", "bias"]) - abs(pooled.loc["viscosity_kin_cSt", "baseline_bias"]))
    dd = float(abs(pooled.loc["density_kg_m3", "bias"]) - abs(pooled.loc["density_kg_m3", "baseline_bias"]))
    n_below = int((du.unmapped_subreps_v5 < 0.10).sum())
    res = dict(admitted=tab[tab.admitted.astype(bool)].name.tolist(), replaced_bins={f"{k[0]}:{k[1]}": reg.smiles[i] for k, i in r1.replaced.items()},
               n_t2_fuels_changed=len(ch), diesel_below_10pct_v5=n_below, diesel_below_10pct_v4=int((du.unmapped_subreps < 0.10).sum()), n_diesel=len(du),
               reading_diesel_pass=n_below >= 7, visc40_abs_bias_change_ln=dv, density15_abs_bias_change=dd,
               reading_t2_pass=bool(dv <= 0.02 and dd <= 2.0))
    res["reading_pass"] = res["reading_diesel_pass"] and res["reading_t2_pass"]
    (OUT / "readings.json").write_text(json.dumps(res, indent=2))
    print(tab.to_string())
    print(d[["property", "set", "variant", "n_groups", "bias", "baseline_bias", "mae", "baseline_mae", "unsafe_rate", "baseline_unsafe_rate"]].to_string())
    print(du.to_string())
    print(json.dumps(res, indent=2))


if __name__ == "__main__":
    main()
