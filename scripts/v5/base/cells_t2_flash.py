"""v5 BASE Task 5: T2 flash cells (measured composition) with input variants.

Variants (in-memory registry overrides through RouteCtx; nothing on disk changes):
  designated, lfl_overlay (tabulated LFL replaces still-estimated ones), mesitylene_fixed (WebBook Antoine set), lfl_overlay+mesitylene_fixed.
Plus model = api_as_shipped (composition_kind="known"), reported under variant "designated".
Outputs: results/v5/BASE/{cells_t2_flash.csv, rows_t2_flash.csv, rows_t2_flash_outliers.csv, rows_t2_flash_api_differences.csv}
"""
from __future__ import annotations

import json
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "src"))
from fuelblend.api.predict_blend import predict  # noqa: E402
from fuelblend.baselines.flash_unifac import _default_registry, _lookup  # noqa: E402
from fuelblend.data.reconstruction_r3 import SubReps  # noqa: E402
from fuelblend.eval.cellkit import frame, make_row  # noqa: E402
from fuelblend.eval.composition_route import RouteCtx, kind_bins, map_bins, predict_properties  # noqa: E402
from fuelblend.eval.group_stats import logo_constant  # noqa: E402
from fuelblend.mixture.common import CompoundRegistry, PurePropertyProvider  # noqa: E402
from fuelblend.uq.conformal import coverage  # noqa: E402

V5 = ROOT / "data" / "real_fuels" / "v5"
OUT = ROOT / "results" / "v5" / "BASE"
SRC = {"av3322": "av3322", "fuellib": "fuellib", "ajftd": "ajftd_gcxgc"}
TABULATED = ("iec60079_20_1_table", "cooptima_compiled")
MES_SMILES = "Cc1cc(C)cc(C)c1"
OUTLIER_LABELS = ["7695", "9324", "8092", "8438"]
VARIANTS = ["designated", "lfl_overlay", "mesitylene_fixed", "lfl_overlay+mesitylene_fixed"]
LAB_FLASH_MARGIN = 3.698768343688585          # configs/uncertainty_v4.yaml flash_point_C.lab_mixtures.ideal (over side)
OVER5 = lambda e: e > 5.0                                                                   # noqa: E731


def lfl_overlay_dict() -> dict:
    reg = _default_registry().copy()
    ov = pd.read_csv(ROOT / "data/processed/v5/lfl_v5.csv")
    ov = ov[ov.lfl_source_v5.isin(TABULATED)].set_index("inchikey")
    hit = reg.inchikey.isin(ov.index) & reg.lfl_vol_pct.notna()
    ch = reg.loc[hit, ["inchikey", "lfl_source"]]
    ch = ch[ch.lfl_source.astype(str).str.startswith("estimated")]            # never overwrite a measured value
    return {ik: float(ov.lfl_vol_pct[ik]) for ik in ch.inchikey}


def mesitylene_dict() -> dict:
    j = json.load(open(ROOT / "results/v5/S2/mesitylene_best_input.json"))
    fr = _default_registry()
    from rdkit import Chem
    can = Chem.MolToSmiles(Chem.MolFromSmiles(MES_SMILES))
    return {fr.loc[can, "inchikey"]: tuple(j["abc"])}


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    warnings.simplefilter("ignore")
    reg = CompoundRegistry(ROOT / "data/unified/compounds.csv")
    prov = PurePropertyProvider(reg, ROOT / "data/cache/pure_grid.npz", verbose=False)
    reps = SubReps(reg, prov, vp_lookup=lambda s: _lookup(s, None))
    lfl, mes = lfl_overlay_dict(), mesitylene_dict()
    ctxs = {"designated": RouteCtx(reg, prov), "lfl_overlay": RouteCtx(reg, prov, lfl_override=lfl),
            "mesitylene_fixed": RouteCtx(reg, prov, antoine_override=mes),
            "lfl_overlay+mesitylene_fixed": RouteCtx(reg, prov, lfl_override=lfl, antoine_override=mes)}
    t = pd.read_csv(OUT / "t2_targets.csv")
    t = t[t.property == "flash_point_C"]
    fu = pd.read_csv(OUT / "t2_fuels.csv").set_index("entity")
    comps = {"av3322": pd.read_csv(V5 / "crc_av3322_composition.csv"), "fuellib": pd.read_csv(V5 / "fuellib_gcxgc_compositions.csv"),
             "ajftd": pd.read_csv(V5 / "ajftd_gcxgc_composition.csv")}
    rows, api_rows = [], []
    for r in t.itertuples():
        src = fu.loc[r.entity, "source"]
        fid = int(r.fuel_id) if src == "ajftd" else r.fuel_id
        c = comps[src]
        mf = map_bins(kind_bins(src, c[c.fuel_id == fid]), reps)
        row = dict(entity=r.entity, source=SRC[src], fuel_id=str(r.fuel_id), label=r.label, group=r.group, method_family=r.method_family,
                   measured=r.measured, low_coverage=bool(fu.loc[r.entity, "low_coverage"]), chosen=bool(fu.loc[r.entity, "chosen_for_group"]),
                   already_read_locked=bool(r.already_read_locked), unmapped_fraction=mf.unmapped_fraction)
        tg = pd.DataFrame([dict(property="flash_point_C", T_K=np.nan)])
        for v, ctx in ctxs.items():
            o = predict_properties(ctx, mf, tg).iloc[0]
            row["pred_" + v.replace("+", "_and_") if v != "designated" else "pred_designated"] = o.pred
            if v == "designated":
                row["status"] = o.status
        # shipped API on the same mapped composition
        try:
            a = predict([(reg.ik[int(i)], float(w)) for i, w in zip(mf.idx, mf.w)], T_K=298.15, fuel_type="jet", composition_kind="known")
            fl = a["flash_point_C"]
            row.update(api_value=fl["value"], api_lo=fl["lo"], api_margin=fl["margin"], api_in_domain=bool(fl["in_domain"]), api_method=fl["method"])
        except Exception as e:                                                                  # noqa: BLE001
            row.update(api_value=np.nan, api_lo=np.nan, api_margin=np.nan, api_in_domain=False, api_method="error", api_error=str(e)[:80])
        rows.append(row)
    d = pd.DataFrame(rows).rename(columns={"pred_lfl_overlay_and_mesitylene_fixed": "pred_both"})
    d["status"] = d.status.fillna("")
    d.to_csv(OUT / "rows_t2_flash.csv", index=False, float_format="%.8g")
    diff = d[(d.api_value - d.pred_designated).abs() > 0.05]
    diff[["entity", "pred_designated", "api_value", "api_lo", "api_margin"]].to_csv(OUT / "rows_t2_flash_api_differences.csv", index=False)

    # ---- outliers (stay in every statistic)
    ac = comps["ajftd"].assign(label=lambda x: x.fuel_label.astype(str))
    o = ac[ac.label.isin(OUTLIER_LABELS) & (ac.n_carbon < 8)].groupby(["label", "hc_class"]).mass_pct.sum().unstack(fill_value=0)
    o["mass_pct_below_C8"] = o.sum(axis=1)
    o = o.reset_index().merge(d[d.source == "ajftd_gcxgc"][["label", "entity", "measured", "pred_designated"]].assign(label=lambda x: x.label.astype(str)),
                              on="label", how="left")
    o["err_designated"] = o.pred_designated - o.measured
    o.to_csv(OUT / "rows_t2_flash_outliers.csv", index=False, float_format="%.6g")

    # ---- cells
    cells = []
    ev = "results/v5/BASE/rows_t2_flash.csv"
    head = d[d.chosen & ~d.low_coverage & (d.status == "ok")]
    n_low = int((d.chosen & d.low_coverage).sum())
    sets = {"pooled": head, **{s: head[head.source == s] for s in SRC.values()},
            **{f"family_{f}": head[head.method_family == f] for f in sorted(head.method_family.dropna().unique())}}
    predcol = {"designated": "pred_designated", "lfl_overlay": "pred_lfl_overlay", "mesitylene_fixed": "pred_mesitylene_fixed",
               "lfl_overlay+mesitylene_fixed": "pred_both"}
    for sname, s in sets.items():
        if s.empty:
            continue
        const = logo_constant(s.measured, s.group, s.method_family)
        et = const - s.measured.to_numpy()
        for v in VARIANTS:
            e = s[predcol[v]] - s.measured
            k, n = coverage(e, s.group, LAB_FLASH_MARGIN, direction="over")
            cells.append(make_row("flash_point_C", "T2_measured", sname, "ideal_le_chatelier", v, "recomputed_v5", "descriptive", ev,
                                  err=e, groups=s.group, unsafe=OVER5, err_trivial=et, trivial_model="logo_constant_within_method",
                                  margin90=LAB_FLASH_MARGIN, coverage90=k / n, n_low_coverage=n_low))
        if sname in ("pooled", "av3322", "ajftd_gcxgc", "fuellib"):
            ea = s.api_value - s.measured
            cov = float(s.assign(ok=s.measured >= s.api_lo).groupby("group").ok.all().mean())
            mc = s.api_method.value_counts().to_dict()
            cells.append(make_row("flash_point_C", "T2_measured", sname, "api_as_shipped", "designated", "recomputed_v5", "descriptive", ev,
                                  err=ea, groups=s.group, unsafe=OVER5, err_trivial=et, trivial_model="logo_constant_within_method",
                                  margin90=float(s.api_margin.median()) if s.api_margin.notna().any() else np.nan, coverage90=cov,
                                  n_low_coverage=n_low, extra=dict(note="API flash method counts: " + json.dumps(mc)
                                                                    + "; non-Le-Chatelier methods are uncertified fallbacks (no margin)")))
    # low-coverage fuels scored apart (not in the headline)
    lc = d[d.low_coverage & (d.status == "ok")]
    for v in VARIANTS:
        if len(lc):
            cells.append(make_row("flash_point_C", "T2_measured", "low_coverage_apart", "ideal_le_chatelier", v, "recomputed_v5", "descriptive", ev,
                                  err=lc[predcol[v]] - lc.measured, groups=lc.group, unsafe=OVER5, n_low_coverage=len(lc)))
    c = frame(cells)
    c.to_csv(OUT / "cells_t2_flash.csv", index=False, float_format="%.8g")
    print(c[["set", "model", "variant", "n_groups", "n_rows", "bias", "mae", "unsafe_rate", "mae_trivial", "verdict", "margin90", "coverage90"]].to_string())
    print("api differences >0.05 C:", len(diff), "of", len(d))


if __name__ == "__main__":
    main()
