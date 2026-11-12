"""v5 BASE Task 9: assemble every cells_*.csv into scoreboard.csv, data_gaps.csv and SCOREBOARD.md.

Reads only committed cells_*.csv / rows_*.csv in results/v5/BASE.  Nothing is scored here.
"""
from __future__ import annotations

import glob
import re
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
BASE = ROOT / "results" / "v5" / "BASE"
TIERS = ["T1_lab", "T2_measured", "CRC647_r3", "T3_ajftd", "diesel"]
PROPS = ["viscosity_kin_cSt", "density_kg_m3", "flash_point_C", "freezing_point_C", "heat_of_combustion_net_MJ_kg",
         "smoke_point_mm", "distillation_T10_C", "pour_point_C", "dcn"]

# (property, tier) pairs transcribed from the cell table in the plan (every table cell, incl. no_evaluation ones).
# Distillation T10/T50/T90 are one family (distillation_T10_C); cetane/DCN names are one family (dcn).
REQUIRED_CELLS: list[tuple[str, str]] = [(p, t) for p in PROPS for t in TIERS] + [("ron", "T1_lab"), ("mon", "T1_lab"), ("designs_certified", "designs")]


def norm_prop(p: str) -> str:
    if re.fullmatch(r"distillation_T\d+_C", p):
        return "distillation_T10_C"
    if p in ("cetane_number_D613", "dcn_D6890"):
        return "dcn"
    return p


def load_cells() -> pd.DataFrame:
    fs = sorted(glob.glob(str(BASE / "cells_*.csv")))
    d = pd.concat([pd.read_csv(f).assign(cells_file=Path(f).name) for f in fs], ignore_index=True)
    d["property_family"] = d.property.map(norm_prop)
    return d


def what_would_fill_it(r) -> str:
    reason = str(r.get("reason", ""))
    t, p, s = r.tier, r.property_family, str(r.set)
    if r.status == "no_evaluation":
        if "Part C Task 16f" in reason:
            return "Part C Task 16f"
        if "Part C Task 10" in reason:
            return "Part C Task 10"
        if t == "diesel":
            return "Part C Task 12 (representatives) and U4 (Mueller 2016 diesel data)"
        if p in ("smoke_point_mm", "distillation_T10_C") and t == "T2_measured":
            return "Part C Task 12 (representatives that carry YSI and Antoine entries)"
        if p == "pour_point_C":
            return "no scheduled fill: needs pure-component pour points or measured pour points (new data)"
        if p == "dcn":
            return "no scheduled fill: needs a cetane model and measured cetane on known compositions (U4)"
        if t == "T1_lab":
            return "no scheduled fill: lab-mixture literature has no such measurements"
        if t == "CRC647_r3" and p in ("smoke_point_mm", "distillation_T10_C"):
            return "not a target here (input of the reconstruction) or no data; U5 would add measured composition"
        return "no scheduled fill: property is an input of this tier, not a target"
    # n_too_small
    if t == "diesel":
        return "Part C Task 12 and U4"
    if r.label == "already_read" or "locked" in s:
        return "Part C Task 15 (fresh holdout); the locked set is spent"
    if t == "CRC647_r3":
        return "U5 (measured composition for the 57 fuels) or more fuels; fold_b has 12"
    if "fuellib" in s:
        return "Part C Task 15 / U3 (more measured-composition fuels)"
    if t == "T1_lab":
        return "U2 (Li et al. 2014) or more lab systems"
    return "more groups: Part C Task 15, U3"


def gaps(sb: pd.DataFrame) -> pd.DataFrame:
    g = sb[(sb.status == "no_evaluation") | (sb.verdict == "n_too_small")].copy()
    g["reason"] = g.reason.where(g.reason.notna(), "n_groups = " + g.n_groups.astype("Int64").astype(str) + " < 8 (verdict n_too_small)")
    g["what_would_fill_it"] = g.apply(what_would_fill_it, axis=1)
    cols = ["property", "property_family", "tier", "set", "model", "variant", "label", "status", "verdict", "n_groups", "reason", "what_would_fill_it"]
    return g[cols].sort_values(["property_family", "tier", "set", "model", "variant"]).reset_index(drop=True)


# ------------------------------------------------------------------------------------------------ R-BASE readings
def r_base_2() -> pd.DataFrame:
    """recomputed_v5 cells with a committed counterpart: observed vs committed, tolerance per property."""
    rows = []

    def add(cell, observed, committed, tol, source, cause=""):
        diff = observed - committed
        rows.append(dict(cell=cell, observed=observed, committed=committed, diff=diff, tol=tol, within=bool(abs(diff) <= tol), committed_source=source, cause=cause))

    r = pd.read_csv(BASE / "rows_t2_flash.csv")
    a = r[(r.source == "ajftd_gcxgc") & (r.status == "ok")]
    e = a.pred_designated - a.measured
    g1 = pd.read_csv(ROOT / "results/v5/G1/summary.csv")
    g1 = g1[(g1.set == "all eligible") & (g1.property == "flash_point_C") & (g1.model == "ideal")].iloc[0]
    add("flash T2 AJFTD GCxGC, all 37 eligible rows: bias", float(e.mean()), float(g1.bias), 0.05, "results/v5/G1/summary.csv")
    add("flash T2 AJFTD GCxGC, all 37 eligible rows: MAE", float(e.abs().mean()), float(g1.mae), 0.05, "results/v5/G1/summary.csv")
    c = pd.read_csv(BASE / "cells_t2_flash.csv")
    h = c[(c.set == "ajftd_gcxgc") & (c.model == "ideal_le_chatelier") & (c.variant == "designated")].iloc[0]
    add("flash T2 AJFTD GCxGC, headline one-per-group (31): bias", float(h.bias), float(g1.bias), 0.05, "results/v5/G1/summary.csv",
        "BASE scores one fuel per double-count group (31); G1 scored all 37 eligible rows (6 duplicate-group fuels); the 37-row rerun above matches G1")
    vt = pd.read_csv(BASE / "viscosity_error_vs_T.csv")
    av = vt[vt.source == "av3322"].set_index("T_K")
    for T, ref in ((233.15, -0.637), (313.15, 0.066)):
        add(f"ln viscosity AV-33-22 {T} K: bias", float(av.loc[T, "bias_ln"]), ref, 0.005, "results/v5/C1 (R-C1-5, rounded to 3 decimals)",
            "C1 used S3's mapping route; BASE uses the shared composition route" if abs(av.loc[T, "bias_ln"] - ref) > 0.005 else "")
    carried = pd.read_csv(BASE / "cells_carried.csv")
    f = carried[(carried.property == "flash_point_C") & (carried.set == "crc647_r3") & (carried.label == "dev") & (carried.variant == "designated")
                & (carried.model == "ideal_le_chatelier")].iloc[0]
    l1 = pd.read_csv(ROOT / "results/v5/L1/comparison.csv")
    l1 = l1[(l1.scope == "real_fuels_CRC647") & (l1.slice == "D93 dev") & (l1.model == "ideal")].iloc[0]
    add("flash CRC 647 r3 dev (45): bias", float(f.bias), float(l1.bias_C_before), 0.05, "results/v5/L1/comparison.csv (before)")
    z = pd.read_csv(BASE / "cells_t2_bulk.csv")
    fz = z[(z.property == "freezing_point_C") & (z.set == "pooled") & (z.variant == "f1_recalibrated") & (z.label == "descriptive")].iloc[0]
    add("freezing T2 pooled F1: MAE", float(fz.mae), 3.46, 0.05, "results/v5/S5 (README, rounded)")
    return pd.DataFrame(rows)


def r_base_3(sb: pd.DataFrame) -> pd.DataFrame:
    d = sb[(sb.variant == "designated") & sb.label.isin(["dev", "descriptive", "fold_b", "exploratory"]) & sb.T_K.isna() & (sb.status != "no_evaluation")]
    d = d[d.verdict.isin(["beats_trivial", "no_skill_shown", "worse_than_trivial", "n_too_small"])]
    out = []
    for p, g in d.groupby("property_family"):
        def tiers(v):
            return sorted(g[g.verdict == v].tier.unique())
        out.append(dict(property=p, beats=", ".join(tiers("beats_trivial")) or "-", no_skill=", ".join(tiers("no_skill_shown")) or "-",
                        worse=", ".join(tiers("worse_than_trivial")) or "-", too_small_only=", ".join(sorted(set(g.tier) - set(g[g.verdict != "n_too_small"].tier))) or "-"))
    return pd.DataFrame(out)


def r_base_4(sb: pd.DataFrame) -> pd.DataFrame:
    d = sb[sb.unsafe_rate.notna() & (sb.status != "no_evaluation") & (sb.variant == "designated")]
    out = []
    for p, g in d.groupby("property_family"):
        r = g.loc[g.unsafe_rate.idxmax()]
        out.append(dict(property=p, max_unsafe_rate=r.unsafe_rate, tier=r.tier, set=r.set, model=r.model, variant=r.variant, label=r.label, n_groups=r.n_groups))
    return pd.DataFrame(out)


def decisions(sb: pd.DataFrame) -> pd.DataFrame:
    """Designated vs variant, side by side, for the four pending adoptions."""
    groups = {"LFL overlay": ["lfl_overlay", "unifac_lfl_overlay", "lfl_overlay+mesitylene_fixed"],
              "Corrected mesitylene Antoine set": ["mesitylene_fixed", "mesitylene_fix", "lfl_overlay+mesitylene_fixed"],
              "K1 cold offset": ["k1_cold_offset_loso"],
              "SLE2 two-compilations tier": ["sle2_two_compilations", "sle2_single_compilation"]}
    des = sb[sb.variant == "designated"]
    out = []
    for name, vs in groups.items():
        for v in sb[sb.variant.isin(vs) & (sb.status != "no_evaluation")].itertuples():
            tk = -1 if pd.isna(v.T_K) else v.T_K
            m = des[(des.property == v.property) & (des.tier == v.tier) & (des.set == v.set) & (des.T_K.fillna(-1) == tk) & (des.label == v.label)]
            m = m[m.model == ("refutas_raw" if name == "K1 cold offset" else v.model)]
            if len(m) > 1:
                m = m[m.cells_file == v.cells_file]
            if m.empty:
                continue
            d0 = m.iloc[0]
            out.append(dict(decision=name, property=v.property, tier=v.tier, set=v.set, T_K=v.T_K, model=v.model, variant=v.variant, label=v.label,
                            n_groups=v.n_groups, bias_designated=d0.bias, bias_variant=v.bias, mae_designated=d0.mae, mae_variant=v.mae,
                            unsafe_designated=d0.unsafe_rate, unsafe_variant=v.unsafe_rate))
    return pd.DataFrame(out)


# ------------------------------------------------------------------------------------------------ markdown
def fmt(x, nd=3):
    return "" if pd.isna(x) else (f"{x:.{nd}g}" if isinstance(x, (float, np.floating)) else str(x))


def md_table(df: pd.DataFrame) -> str:
    cols = list(df.columns)
    lines = ["| " + " | ".join(cols) + " |", "|" + "|".join("---" for _ in cols) + "|"]
    for r in df.itertuples(index=False):
        lines.append("| " + " | ".join(fmt(v) for v in r) + " |")
    return "\n".join(lines)


def main() -> None:
    sb = load_cells()
    sb = sb.sort_values(["property", "tier", "set", "model", "variant", "label", "T_K"], kind="stable").reset_index(drop=True)
    sb.to_csv(BASE / "scoreboard.csv", index=False, float_format="%.8g")
    gp = gaps(sb)
    gp.to_csv(BASE / "data_gaps.csv", index=False, float_format="%.8g")

    have = set(zip(sb.property_family, sb.tier))
    missing = [c for c in REQUIRED_CELLS if (norm_prop(c[0]), c[1]) not in have]
    r2, r3, r4, dec = r_base_2(), r_base_3(sb), r_base_4(sb), decisions(sb)

    show = ["tier", "set", "model", "variant", "label", "T_K", "n_groups", "bias", "mae", "unsafe_rate", "verdict", "margin90", "coverage90"]
    md = ["# v5 BASE scoreboard", "",
          f"{len(sb)} cells; {int((sb.status == 'no_evaluation').sum())} `no_evaluation`, {int((sb.verdict == 'n_too_small').sum())} `n_too_small` "
          f"(see `data_gaps.csv`).  Every T2 fuel has been read, so T2 cells are `descriptive`, not held-out.  Error is predicted minus measured; "
          "viscosity error is ln(pred/meas).  Bias, MAE and unsafe-side rate are in the property's unit (viscosity: ln).", "",
          "## Readings", "",
          f"**R-BASE-1** ({len(REQUIRED_CELLS)} required (property, tier) cells): " + ("all exist with a status." if not missing else f"MISSING {missing}"), "",
          "**R-BASE-2** (recomputed cells with a committed counterpart):", "", md_table(r2), "",
          "**R-BASE-3** (designated model vs the leave-one-group-out constant, aggregate rows without temperature split, by tier; a tier appears in several columns when its sets or models disagree, see the per-property tables):", "", md_table(r3), "",
          "**R-BASE-4** (largest unsafe-side rate per property, designated rows only; flash = over-prediction > 5 C, cold viscosity = under-prediction beyond the S1 margin):", "", md_table(r4), "",
          f"**R-BASE-5** data gaps: `data_gaps.csv`, {len(gp)} rows ({int((gp.status == 'no_evaluation').sum())} no_evaluation, {int((gp.verdict == 'n_too_small').sum())} n_too_small).", "",
          "## Decisions waiting on the user (U1)", "",
          "Designated (v4 routing) stays unchanged until the user decides.  Each row: the same cell scored with and without the change.", "", md_table(dec), ""]
    for p in PROPS + ["ron", "mon", "designs_certified"]:
        g = sb[(sb.property_family == p) & sb.T_K.isna() & (sb.status != "no_evaluation")]
        if g.empty:
            continue
        md += [f"## {p}", "", "Designated rows first, then variants and alternative models.  Rows with a temperature split are in `scoreboard.csv`.", "",
               md_table(pd.concat([g[g.variant == "designated"], g[g.variant != "designated"]])[show]), ""]
    (BASE / "SCOREBOARD.md").write_text("\n".join(md), encoding="utf-8", newline="\n")
    print(f"scoreboard {len(sb)} rows; gaps {len(gp)}; missing required: {missing}")
    print(r2[["cell", "observed", "committed", "diff", "tol", "within"]].to_string())


if __name__ == "__main__":
    main()
