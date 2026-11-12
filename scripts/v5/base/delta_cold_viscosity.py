"""v5 BASE Task 11: cold-viscosity K1 correction with its own cross-fitted margin, and the pure-provider cause check.

Residuals are the leave-source-out K1 residuals (offset fitted on the other source with the scored double-count group removed).
Outputs (results/v5/BASE/delta_cold_viscosity/): delta.csv, proposed_margins.yaml (NOT configs/), pure_cold_extrapolation.csv, readings.json
"""
from __future__ import annotations

import importlib.util
import json
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "src"))
from fuelblend.eval.cellkit import frame, make_row  # noqa: E402
from fuelblend.eval.cold_offset import apply_k1, fit_k1  # noqa: E402
from fuelblend.uq.conformal import coverage, two_sided_margin  # noqa: E402

OUT = ROOT / "results" / "v5" / "BASE" / "delta_cold_viscosity"
V5 = ROOT / "data" / "real_fuels" / "v5"
COLD = {"-40C": 233.15, "-20C": 253.15}
RAW_S1_MARGIN = 0.4153576080066075                    # raw Refutas margin below 263 K (results/v5/S1)

_spec = importlib.util.spec_from_file_location("cells_t2_bulk_mod", ROOT / "scripts/v5/base/cells_t2_bulk.py")
bulk = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(bulk)


def margin(resid, groups) -> float:
    """Two-sided conformal margin on |ln residual|; inf when there are fewer than 9 independent groups."""
    return float(two_sided_margin(np.asarray(resid, float), np.asarray(groups).astype(str), 0.10))


def crossfit(rows: pd.DataFrame) -> pd.DataFrame:
    """Leave-source-out K1 for every scored viscosity row.  Adds err0, err_k1, pred_k1, train_source, train_groups (frozenset)."""
    v = rows[(rows.property == "viscosity_kin_cSt") & (rows.variant == "designated") & (rows.status == "ok") & ~rows.low_coverage.astype(bool)].copy()
    v["resid"] = np.log(v.measured) - np.log(v.pred)
    out = []
    for grp, s in v.groupby("group"):
        for src in s.source.unique():
            tsrc = bulk.TRAIN_FROM[src]
            tr = pd.concat([bulk.one_per_group(g) for _, g in v[(v.source == tsrc) & (v.group != grp)].groupby("T_K")])
            k1 = fit_k1(tr.T_K, tr.resid)
            sc = s[s.source == src]
            off = apply_k1(k1, sc.T_K.to_numpy())
            out.append(sc.assign(pred_k1=np.exp(np.log(sc.pred) + off), train_source=tsrc, train_groups=[frozenset(tr.group)] * len(sc)))
    k = pd.concat(out, ignore_index=True)
    k["err0"] = np.log(k.pred) - np.log(k.measured)
    k["err_k1"] = np.log(k.pred_k1) - np.log(k.measured)
    return k.drop(columns=["resid"])


def scored(cf: pd.DataFrame, T: float) -> pd.DataFrame:
    s = cf[(cf.T_K == T) & np.isfinite(cf.err_k1)]
    return bulk.one_per_group(s) if len(s) else s


def pure_check() -> pd.DataFrame:
    from fuelblend.mixture import common
    from fuelblend.mixture.common import CompoundRegistry, PurePropertyProvider
    reg = CompoundRegistry(ROOT / "data" / "unified" / "compounds.csv")
    prov = PurePropertyProvider(reg, ROOT / "data" / "cache" / "pure_grid.npz", verbose=False)
    ik_index = {k: i for i, k in enumerate(reg.ik)}
    comp = pd.read_csv(V5 / "crc_av3322_composition.csv")
    amap = pd.read_csv(ROOT / "results/v5/S2/av3322_representative_map.csv")
    m = comp.merge(amap[["group", "carbon_number", "rep_inchikey", "mapped"]], on=["group", "carbon_number"], how="left")
    m = m[m.mapped.astype(bool)]
    share = m.groupby(["fuel_id", "rep_inchikey"]).mass_pct.sum().groupby("rep_inchikey").sum() / comp.fuel_id.nunique()
    boehm = pd.read_csv(V5 / "boehm2022_viscosity.csv")
    boehm = boehm[(boehm.kind == "pure") & (boehm.property == "viscosity_kin_cSt")].copy()
    boehm["ik"] = boehm.comp_inchikeys.map(lambda s: json.loads(s)[0])
    fe = common._get_measured_fits()[1]
    rows = []
    for ik, mass in share[share > 2.0].sort_values(ascending=False).items():
        i = ik_index[ik]
        rec = dict(smiles=reg.smiles[i], inchikey=ik, mean_mass_pct=float(mass))
        rec["data_T_min_K"], rec["data_T_max_K"] = (fe[ik][2], fe[ik][3]) if ik in fe else (np.nan, np.nan)
        for lab, T in (("253K", 253.15), ("233K", 233.15)):
            ext = (ik not in fe) or T < fe[ik][2] - 25.0 or T > fe[ik][3] + 25.0
            rec[f"source_at_{lab}"] = "ml_extrapolation" if ext else "fit"
            le, rho = prov.lookup(np.array([[i]]), T)
            nu = float(1000.0 * np.exp(le[0, 0]) / rho[0, 0])
            b = boehm[(boehm.ik == ik) & (np.abs(boehm.T_K - T) < 0.3)]
            rec[f"err_vs_boehm_{lab}"] = float(np.log(nu / b.value.median())) if len(b) else np.nan
        rows.append(rec)
    t = pd.DataFrame(rows)
    t.attrs["total_mapped_mass"] = float(share.sum())
    return t


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    warnings.simplefilter("ignore")
    rows = pd.read_csv(ROOT / "results/v5/BASE/rows_t2_bulk.csv")
    sb = pd.read_csv(ROOT / "results/v5/BASE/scoreboard.csv")
    cf = crossfit(rows)
    drows, props, checks = [], {}, {}
    pairs = (("av3322", "ajftd_gcxgc"), ("ajftd_gcxgc", "av3322"))
    for lab, T in COLD.items():
        s = scored(cf, T)
        by = {src: bulk.one_per_group(g) for src, g in cf[(cf.T_K == T) & np.isfinite(cf.err_k1)].groupby("source")}   # within-source, as the scoreboard does
        m_src = {src: margin(g.err_k1, g.group) for src, g in by.items()}
        props[lab] = dict(T_K=T, pooled=dict(margin=margin(s.err_k1, s.group), n_groups=int(s.group.nunique())),
                          per_source={k: dict(margin=m_src[k], n_groups=int(by[k].group.nunique())) for k in by},
                          raw_s1_margin=RAW_S1_MARGIN)
        cov = {}
        for setter, tested in pairs:
            key = f"margin_from_{setter}_on_{tested}"
            if setter in by and tested in by and np.isfinite(m_src[setter]):
                k, n = coverage(by[tested].err_k1, by[tested].group, m_src[setter], two_sided=True)
                cov[key] = dict(margin=m_src[setter], covered=int(k), n_groups=int(n), coverage=k / n)
            else:
                cov[key] = dict(margin=m_src.get(setter, float("inf")), coverage=None, note="not evaluable: margin inf or source absent")
        checks[lab] = cov
        for src, g in by.items():
            other = bulk.TRAIN_FROM[src]
            mg = m_src.get(other, float("inf"))
            k, n = coverage(g.err_k1, g.group, mg, two_sided=True) if np.isfinite(mg) else (0, 0)
            b = sb[(sb.property == "viscosity_kin_cSt") & (sb.tier == "T2_measured") & (sb.set == src) & (sb.model == "refutas_raw") &
                   (sb.variant == "designated") & (sb.T_K == T)].iloc[0]
            drows.append(make_row("viscosity_kin_cSt", "T2_measured", src, "refutas_raw+k1", "k1_cold_offset_loso", "recomputed_v5", "descriptive",
                                  "results/v5/BASE/rows_t2_bulk.csv (crossfit)", err=g.err_k1, groups=g.group,
                                  unsafe=lambda e, m=RAW_S1_MARGIN: e < -m, margin90=mg, coverage90=(k / n if n else np.nan),
                                  extra=dict(T_K=T, baseline_mae=b.mae, baseline_bias=b.bias, baseline_unsafe_rate=b.unsafe_rate)))
    frame(drows).to_csv(OUT / "delta.csv", index=False)
    yaml.safe_dump(dict(note="PROPOSED, not adopted. Two-sided ln margins on leave-source-out K1 residuals, alpha=0.10; inf = fewer than 9 groups.",
                        margins=props, cross_source_coverage=checks), open(OUT / "proposed_margins.yaml", "w"), sort_keys=False)
    t = pure_check()
    t.to_csv(OUT / "pure_cold_extrapolation.csv", index=False, float_format="%.5g")
    tot = t.attrs["total_mapped_mass"]
    ml = {lab: float(t.loc[t[f"source_at_{lab}"] == "ml_extrapolation", "mean_mass_pct"].sum() / tot) for lab in ("253K", "233K")}
    m20 = props["-20C"]["pooled"]["margin"]
    c20 = [v.get("coverage") for v in checks["-20C"].values()]
    res = dict(reading_margin_below_raw=bool(m20 < RAW_S1_MARGIN), margin_minus20=m20, raw=RAW_S1_MARGIN,
               reading_coverage_other_source_ge_0p90=bool(all(c is not None and c >= 0.90 for c in c20)), coverage_minus20=c20,
               ml_extrapolation_mass_fraction_of_mapped_mass=ml, ml_majority_253K=ml["253K"] > 0.5, ml_majority_233K=ml["233K"] > 0.5)
    (OUT / "readings.json").write_text(json.dumps(res, indent=2))
    print(yaml.safe_dump(dict(margins=props, cross=checks), sort_keys=False))
    print(t.to_string())
    print(json.dumps(res, indent=2))


if __name__ == "__main__":
    main()
