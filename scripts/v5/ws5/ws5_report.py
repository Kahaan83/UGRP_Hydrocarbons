"""WS5 report stage: reads results/v5/BASE/delta_ws5_ghgnn/{rows.csv, ad_pairs.csv} and applies the READMEs declared gate (R-1..R-5).
Writes summary.csv, margins.csv, ad_summary.json, readings.json.  No parameter here is chosen after seeing results."""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import ws5_flash as F  # noqa: E402

OUT = F.OUT
ROOT = F.ROOT


def scoring_frame(d: pd.DataFrame, col: str) -> pd.DataFrame:
    ok = (d.ideal_status == "ok") & np.isfinite(d[col]) & np.isfinite(d.ideal_C)
    return d[ok]


def set_summary(d: pd.DataFrame, name: str, variants: dict[str, str], gate: bool) -> list[dict]:
    out = []
    for vname, col in variants.items():
        s = scoring_frame(d, col)
        if len(s) == 0:
            continue
        e0, e1 = s.ideal_C - s.measured, s[col] - s.measured
        lo, hi = F.cluster_boot_dmae(e0, e1, s.group)
        m0, w0 = F.margin_and_setter(e0, s.group)
        m1, w1 = F.margin_and_setter(e1, s.group)
        d0, d1 = F.describe(e0), F.describe(e1)
        out.append(dict(set=name, variant=vname, n_rows=len(s), n_groups=int(s.group.nunique()), n_dropped=int(len(d) - len(s)),
                        bias_ideal=d0["bias"], bias_gh=d1["bias"], mae_ideal=d0["mae"], mae_gh=d1["mae"], rmse_ideal=d0["rmse"], rmse_gh=d1["rmse"],
                        dMAE=d1["mae"] - d0["mae"], dMAE_ci95_lo=lo, dMAE_ci95_hi=hi, over2_ideal=d0["over2"], over2_gh=d1["over2"],
                        over5_ideal=d0["over5"], over5_gh=d1["over5"], worst_over_ideal=d0["worst_over"], worst_over_gh=d1["worst_over"],
                        margin90_ideal=m0, margin90_gh=m1, margin_setter_ideal=w0, margin_setter_gh=w1,
                        gate_a_ci_below_0=bool(hi < 0), gate_b_over5_not_higher=bool(d1["over5"] <= d0["over5"] + 1e-12),
                        gate_pass=bool(hi < 0 and d1["over5"] <= d0["over5"] + 1e-12), gated=gate and vname == "wilson_all_pairs"))
    return out


def sysname(sid: str, d: pd.DataFrame) -> str:
    return sid


def main() -> None:
    d = pd.read_csv(OUT / "rows.csv")
    ad = pd.read_csv(OUT / "ad_pairs.csv")
    variants = {"wilson_all_pairs": "wilson_C", "wilson_ad_gated": "wilson_ad_C"}
    sm = []
    sm += set_summary(d[d.set == "S-LAB"], "S-LAB", variants, True)
    sm += set_summary(d[d.set == "S-CRC"], "S-CRC", variants, True)
    sm += set_summary(d[d.set == "S-CRCL"], "S-CRCL", variants, False)
    sm += set_summary(d[d.set == "S-T2"], "S-T2", variants, True)
    sm += set_summary(d[d.set == "S-LI"], "S-LI", variants, False)
    sm += set_summary(d[d.set.isin(["S-LAB", "S-LI"])], "S-LAB+S-LI (36 systems)", variants, False)
    # sub-splits (descriptive)
    for src, g in d[d.set == "S-T2"].groupby("method"):
        sm += set_summary(g, f"S-T2/{src}", {"wilson_all_pairs": "wilson_C"}, False)
    crc = d[d.set == "S-CRC"]
    for meth, g in crc.groupby("method"):
        sm += set_summary(g, f"S-CRC/{meth}", {"wilson_all_pairs": "wilson_C"}, False)
    labb = d[(d.set == "S-LAB") & (d.n_comp == 2)]
    sm += set_summary(labb, "S-LAB binaries", {"wilson_all_pairs": "wilson_C", "margules": "margules_C"}, False)
    sm += set_summary(d[(d.set == "S-LAB") & (d.n_comp == 3)], "S-LAB ternaries", {"wilson_all_pairs": "wilson_C"}, False)
    summ = pd.DataFrame(sm)
    summ.to_csv(OUT / "summary.csv", index=False, float_format="%.6g")

    # ---- R-1 overall
    gated = summ[summ.gated]
    per_set = {r.set: dict(dMAE=r.dMAE, ci95=[r.dMAE_ci95_lo, r.dMAE_ci95_hi], over5_ideal=r.over5_ideal, over5_gh=r.over5_gh, gate_pass=r.gate_pass)
               for r in gated.itertuples()}
    overall = bool(len(per_set) == 3 and all(v["gate_pass"] for v in per_set.values()))
    crcl = summ[(summ.set == "S-CRCL") & (summ.variant == "wilson_all_pairs")]
    crcl_unsafe_rise = bool(len(crcl) and crcl.over5_gh.iloc[0] > crcl.over5_ideal.iloc[0] + 1e-12)

    # ---- R-4 margins (lab scopes)
    lab = d[d.set.isin(["S-LAB", "S-LI"])].copy()
    lab["is_old"] = lab.set == "S-LAB"
    names = {}
    for r in lab.itertuples():
        pass
    scopes = {"old_31_systems": lab.is_old, "old_minus_mesitylene": lab.is_old & ~lab.has_mesitylene,
              "old_plus_new": pd.Series(True, index=lab.index), "old_plus_new_minus_mesitylene": ~lab.has_mesitylene, "new_only": ~lab.is_old}
    var_cols = {"ideal": "ideal_C", "wilson": "wilson_C", "ideal_mesitylene_antoine_fixed": "ideal_mesfix_C", "wilson_mesitylene_antoine_fixed": "wilson_mesfix_C",
                "wilson_ad_gated": "wilson_ad_C"}
    mrows = []
    for v, col in var_cols.items():
        for sc, m in scopes.items():
            x = lab[m & np.isfinite(lab[col])]
            e = x[col] - x.measured
            mg, setter = F.margin_and_setter(e, x.group)
            mrows.append(dict(variant=v, scope=sc, margin=mg, n_groups=int(x.group.nunique()), n_points=len(x), setter_system=setter,
                              setter_is_mesitylene=bool(setter is not None and x[x.group == setter].has_mesitylene.any())))
    mt = pd.DataFrame(mrows)
    mt.to_csv(OUT / "margins.csv", index=False, float_format="%.6g")
    mm = mt.set_index(["variant", "scope"]).margin
    r4 = {}
    for v in var_cols:
        w, wo = mm[(v, "old_plus_new")], mm[(v, "old_plus_new_minus_mesitylene")]
        r4[v] = dict(with_mesitylene=float(w), without_mesitylene=float(wo), diff=float(abs(w - wo)), within_0p5=bool(abs(w - wo) < 0.5))
    mes = lab[lab.has_mesitylene]
    mes_err = {k: float((mes[c] - mes.measured).mean()) for k, c in var_cols.items()}
    mes_err_max = {k: float((mes[c] - mes.measured).max()) for k, c in var_cols.items()}

    # ---- AD summaries
    brou = pd.read_csv(F.ROOT / "data/raw/ghgnn_code/repo/data/raw/Brouwer_2021.csv")
    ishc = lambda s: {a.GetSymbol() for a in Chem.MolFromSmiles(s).GetAtoms()} <= {"C"}
    cache = {s: ishc(s) for s in set(brou.Solute_SMILES) | set(brou.Solvent_SMILES)}
    hc = brou[brou.Solute_SMILES.map(cache) & brou.Solvent_SMILES.map(cache)]
    adj = dict(
        n_ordered_pairs=int(len(ad)), share_pairs_in_domain=float(ad.ad_feasible.mean()), n_interpolation_pairs=int(ad.interpolation.sum()),
        n_pairs_out_of_domain=int((~ad.ad_feasible).sum()),
        out_of_domain_reasons=dict(low_tanimoto=int(((ad.tanimoto10 < 0.35) & ~ad.interpolation).sum()), class_count_lt25=int(((ad.n_class_pair < 25) & ~ad.interpolation).sum())),
        class_pair_table=ad.groupby(["class_solvent", "class_solute"]).agg(n_pairs=("ad_feasible", "size"), share_in_domain=("ad_feasible", "mean"),
                                                                             n_class_pair_training=("n_class_pair", "first")).reset_index().to_dict("records"),
        class_source_counts=ad.class_source.value_counts().to_dict(),
        brouwer_hc_rows=int(len(hc)), brouwer_hc_T_C=dict(min=float(hc["T"].min()), p5=float(hc["T"].quantile(0.05)), median=float(hc["T"].median()),
                                                           p95=float(hc["T"].quantile(0.95)), max=float(hc["T"].max())),
        brouwer_hc_unique_pairs=int(hc[["Solute_SMILES", "Solvent_SMILES"]].drop_duplicates().shape[0]),
        per_set={s: dict(pair_weight_in_domain_median=float(g.ad_pair_weight_in_domain.median()), p10=float(g.ad_pair_weight_in_domain.quantile(0.1)),
                         min=float(g.ad_pair_weight_in_domain.min()), n_mixtures=int(len(g)),
                         n_mixtures_all_pairs_in_domain=int((g.ad_pair_weight_in_domain >= 1 - 1e-12).sum()),
                         flash_T_C=dict(min=float(g.flash_T_C.min()), median=float(g.flash_T_C.median()), max=float(g.flash_T_C.max())),
                         pairs_wilson_failed_total=int(g.wilson_fail_pairs.sum()), pairs_total=int(g.n_pairs.sum()),
                         mixtures_with_failed_pairs=int((g.wilson_fail_pairs > 0).sum()),
                         gamma_min=float(g.gamma_min_at_ideal_flash.min()), gamma_max=float(g.gamma_max_at_ideal_flash.max()),
                         wilson_status_not_ok=int((g.wilson_status != "ok").sum()))
                 for s, g in d.groupby("set")})
    (OUT / "ad_summary.json").write_text(json.dumps(adj, indent=1))
    res = dict(gate_per_gated_set=per_set, gate_overall_pass=overall, locked_crc_unsafe_rate_rises=crcl_unsafe_rise,
               R4_lab_margins=r4, R4_mesitylene_system_mean_error=mes_err, R4_mesitylene_system_max_error=mes_err_max,
               expectation_stated="loss or tie likely")
    (OUT / "readings.json").write_text(json.dumps(res, indent=1))
    pd.set_option("display.width", 250, "display.max_columns", 40)
    print(summ[["set", "variant", "n_rows", "n_groups", "bias_ideal", "bias_gh", "mae_ideal", "mae_gh", "dMAE", "dMAE_ci95_lo", "dMAE_ci95_hi",
                "over5_ideal", "over5_gh", "worst_over_ideal", "worst_over_gh", "gate_pass"]].round(3).to_string())
    print(mt.round(3).to_string())
    print(json.dumps(res, indent=1))
    print(json.dumps({k: v for k, v in adj.items() if k != "class_pair_table"}, indent=1))


if __name__ == "__main__":
    main()
