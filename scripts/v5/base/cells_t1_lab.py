"""v5 BASE Task 4: T1 lab-mixture cells (viscosity bins, flash, SLE tiers, density incl. unseen ThermoML 2020 systems).

Aggregates committed S1 / v4 / L1 / S2 / SLE2 files; recomputes only the density cell on ThermoML 2020 systems that were never scored.
Output: results/v5/BASE/cells_t1_lab.csv, rows_t1_density_thermoml2020.csv
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "src"))
from fuelblend.eval.cellkit import frame, make_row, no_eval  # noqa: E402
from fuelblend.mixture.common import CompoundRegistry, PurePropertyProvider  # noqa: E402
from fuelblend.uq.conformal import coverage  # noqa: E402

R = ROOT / "results"
OUT = R / "v5" / "BASE"
T1 = "T1_lab"
LAB_LABEL = {"dev": "dev", "locked": "already_read"}


def viscosity_rows(margins: dict) -> list[dict]:
    v = pd.read_csv(R / "v5/S1/viscosity_lab_residuals_v5.csv")
    v = v[v.method == "refutas"]
    rows = []
    for (tb, role), g in v.groupby(["t_bin", "split_role"]):
        m = float(margins[f"viscosity_kin_cSt.lab_mixtures.raw_refutas.{tb}"]["applied_margin"])
        k, n = coverage(g.resid, g.system_id, m, two_sided=True)
        unsafe = (lambda e, m=m: e < -m) if tb != "T_ge_283K" else None           # unsafe side only reported below 283 K
        rows.append(make_row("viscosity_kin_cSt", T1, tb, "refutas_raw", "designated", "aggregated_v5", LAB_LABEL[role],
                             "results/v5/S1/viscosity_lab_residuals_v5.csv; results/v5/S1/uncertainty_v5.yaml",
                             err=g.resid, groups=g.system_id, unsafe=unsafe, margin90=m, coverage90=k / n if n else np.nan))
    return rows


def flash_rows(margins: dict) -> list[dict]:
    f = pd.read_csv(R / "v4/F/flash_lab_loso.csv")
    m = float(margins["flash_point_C.lab_mixtures.ideal"]["applied_margin"])
    k, n = coverage(f.ideal_err_C, f.system_id, m, direction="over")
    ev = "results/v4/F/flash_lab_loso.csv"
    over5 = lambda e: e > 5.0                                                                  # noqa: E731
    rows = [make_row("flash_point_C", T1, "lab_mixtures_31_systems", "ideal_le_chatelier", "designated", "aggregated_v5", "dev", ev,
                     err=f.ideal_err_C, groups=f.system_id, unsafe=over5, margin90=m, coverage90=k / n)]
    s2 = pd.read_csv(R / "v5/S2/mesitylene_flash_rerun.csv")
    rows.append(make_row("flash_point_C", T1, "lab_mixtures_31_systems", "ideal_le_chatelier", "mesitylene_fixed", "aggregated_v5", "dev",
                         "results/v5/S2/mesitylene_flash_rerun.csv", err=s2.webbook_ideal_err_C, groups=s2.system_id, unsafe=over5))
    c = pd.read_csv(R / "v5/L1/comparison.csv")
    c = c[(c.scope == "lab_mixtures") & c.slice.str.startswith("all_126") & (c.model == "ideal")].iloc[0]
    rows.append(make_row("flash_point_C", T1, "lab_mixtures_31_systems", "ideal_le_chatelier", "lfl_overlay", "aggregated_v5", "dev",
                         "results/v5/L1/comparison.csv",
                         fixed=dict(n_groups=31, n_rows=int(c.N), bias=c.bias_C_after, mae=c.MAE_C_after, unsafe_rate=c.over_gt5C_rate_after)))
    rows.append(no_eval("flash_point_C", T1, "lab_mixtures_31_systems", "lfl_overlay+mesitylene_fixed not run on lab mixtures",
                        variant="lfl_overlay+mesitylene_fixed", evidence="results/v5/L1/comparison.csv"))
    return rows


def freezing_rows() -> list[dict]:
    s = pd.read_csv(R / "v5/SLE2/summary_by_tier.csv")
    s = s[s.system_class == "ALL"]
    name = {"v4_verified": "designated", "two_compilations": "sle2_two_compilations", "single_compilation": "sle2_single_compilation"}
    rows = []
    for r in s.itertuples():
        rows.append(make_row("freezing_point_C", T1, "sle_alkane_blends", "ideal_svl", name[r.tier], "aggregated_v5", "dev",
                             "results/v5/SLE2/summary_by_tier.csv", margin90=float(r.margin_under_K),
                             fixed=dict(n_groups=int(r.n_systems), n_rows=int(r.n_points), bias=r.bias_K, mae=r.mae_K)))
    return rows


def density_v4_rows() -> list[dict]:
    d = pd.read_csv(R / "v4/L/density_lab_residuals.csv")
    rows = []
    for (tb, role), g in d.groupby(["t_bin", "split_role"]):
        rows.append(make_row("density_kg_m3", T1, tb, "ideal_volume", "designated", "aggregated_v5", LAB_LABEL[role],
                             "results/v4/L/density_lab_residuals.csv", err=g.resid, groups=g.system_id))
    return rows


def density_new_rows(margins: dict) -> tuple[list[dict], pd.DataFrame]:
    t = pd.read_csv(ROOT / "data/real_fuels/v5/thermoml2020_lab.csv")
    d = t[(t.property == "density_kg_m3") & (t.kind != "pure")].copy()
    seen = set(pd.read_csv(ROOT / "data/unified/measurements_dedup.csv", usecols=["system_id"], low_memory=False).system_id)
    old_dois = set(pd.read_csv(ROOT / "data/processed/v2/thermoml_clean_hc.csv", usecols=["file"]).file.str.replace(".xml", "", regex=False))
    drop = {"system_in_unified": int(d.system_id.isin(seen).sum()), "doi_in_v2_extract": int(d.ref.isin(old_dois).sum())}
    d = d[~d.system_id.isin(seen) & ~d.ref.isin(old_dois)].copy()
    reg = CompoundRegistry(ROOT / "data/unified/compounds.csv")
    prov = PurePropertyProvider(reg, ROOT / "data/cache/pure_grid.npz", verbose=False)
    pred, why = [], []
    for r in d.itertuples():
        iks, w = json.loads(r.comp_inchikeys), np.array(json.loads(r.w_mass), float)
        idx = [reg.ik2i.get(k, -1) for k in iks]
        if min(idx) < 0:
            pred.append(np.nan); why.append("component_not_in_registry"); continue
        idx = np.array(idx)
        if not prov.valid[idx].all():
            pred.append(np.nan); why.append("pure_property_invalid"); continue
        if not (prov.T[0] <= r.T_K <= prov.T[-1]):
            pred.append(np.nan); why.append("T_outside_provider_grid"); continue
        _, rho = prov.lookup(idx[None, :], r.T_K)
        pred.append(float(1.0 / (w / rho[0]).sum())); why.append("")
    d["pred"], d["why_no_pred"] = pred, why
    d["resid"] = d.pred - d.value
    d[["system_id", "kind", "n_comp", "T_K", "value", "pred", "resid", "why_no_pred", "ref"]].to_csv(
        OUT / "rows_t1_density_thermoml2020.csv", index=False, float_format="%.8g")
    ok = d[np.isfinite(d.resid)]
    mo = float(margins["density_kg_m3.lab_mixtures.ideal_volume.over"]["applied_margin"])
    mu = float(margins["density_kg_m3.lab_mixtures.ideal_volume.under"]["applied_margin"])
    gs = ok.groupby("system_id").resid.agg(["max", "min"])
    cov = float(((gs["max"] <= mo) & (-gs["min"] <= mu)).mean())
    row = make_row("density_kg_m3", T1, "thermoml2020_new_systems", "ideal_volume", "designated", "recomputed_v5", "dev",
                   "results/v5/BASE/rows_t1_density_thermoml2020.csv", err=d.resid, groups=d.system_id, margin90=mo, coverage90=cov,
                   extra=dict(margin90_under=mu, note="first read; rows dropped before scoring: " + json.dumps(drop)
                              + "; no prediction: " + json.dumps(pd.Series(why).value_counts().drop("", errors="ignore").to_dict())))
    return [row], d


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    margins = yaml.safe_load(open(R / "v5/S1/uncertainty_v5.yaml"))["margins"]
    rows = viscosity_rows(margins) + flash_rows(margins) + freezing_rows() + density_v4_rows()
    new, _ = density_new_rows(margins)
    rows += new + [
        no_eval("heat_of_combustion_net_MJ_kg", T1, "-", "no lab-mixture heat-of-combustion data"),
        no_eval("smoke_point_mm", T1, "-", "no lab-mixture smoke-point data"),
        no_eval("distillation_T50_C", T1, "-", "no lab-mixture distillation data"),
        no_eval("pour_point_C", T1, "-", "no open curated pour-point mixture data; PPBI needs pure-component pour points the registry lacks"),
        no_eval("dcn", T1, "-", "v4: uncertified, no model scored"),
    ]
    df = frame(rows)
    df.to_csv(OUT / "cells_t1_lab.csv", index=False, float_format="%.8g")
    print(df[["property", "set", "variant", "label", "n_groups", "n_rows", "bias", "mae", "margin90", "coverage90"]].to_string())


if __name__ == "__main__":
    main()
