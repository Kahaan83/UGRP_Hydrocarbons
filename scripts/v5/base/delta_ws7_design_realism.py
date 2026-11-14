"""delta_ws7_design_realism runner (README declares every reading; nothing here moves a gate).

    python scripts/v5/base/delta_ws7_design_realism.py blendstocks     # R-W7-1 -> blendstocks.csv
    python scripts/v5/base/delta_ws7_design_realism.py margins         # R-W7-2 / R-W7-3 -> d86_fit_rows.csv, d86_margins.json, tb_sanity.csv
    python scripts/v5/base/delta_ws7_design_realism.py compare         # R-W7-5 -> compare/*.csv
    python scripts/v5/base/delta_ws7_design_realism.py design --pool P1 --spec astm:jet_a   # R-W7-6 -> designs/<pool>/
    python scripts/v5/base/delta_ws7_design_realism.py table           # R-W7-6 table with the counterfactual -> r_designs.csv
"""
from __future__ import annotations

import argparse
import json
import re
import sys
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "src"))
warnings.simplefilter("ignore")

from fuelblend.design import blendstock_v5 as B  # noqa: E402
from fuelblend.design import engine_v4 as E  # noqa: E402

D = B.DELTA
V5 = ROOT / "data" / "real_fuels" / "v5"
MARGINS = ROOT / "configs" / "uncertainty_v5.yaml"
POOL_C = ["n-decane", "n-undecane", "n-dodecane", "butylcyclohexane", "decalin_trans", "1,2,4-trimethylbenzene", "isooctane"]
SPECS = ["astm:jet_a", "bis:jet_a1"]
D86_PROPS = {"distillation_T10_C": "T10", "distillation_T50_C": "T50", "distillation_T90_C": "T90", "distillation_T100_C": "FBP"}
SMOKE_POOL_REF = {"astm:jet_a": 0.0187332302631016, "bis:jet_a1": 0.0155447656676983}     # 6dc66f8 r_pool.csv min_slack_constraints


# ------------------------------------------------------------------------------------------------ R-W7-1
def cmd_blendstocks():
    rows = []
    for name in B.BLENDSTOCKS:
        b = B.load_blendstock(name)
        cf = b.class_fractions()
        rows.append(dict(name=name, label=b.label, source=b.source, fuel_id=b.fuel_id, n_representatives=len(b.smiles),
                         unmapped_fraction=b.unmapped_fraction, synthetic=b.synthetic, status="used",
                         **{f"wt_{k}": round(v, 4) for k, v in sorted(cf.items())}))
    for name, why in B.BLOCKED_BLENDSTOCKS.items():
        rows.append(dict(name=name, label=name, source="", status=f"blocked: {why}"))
    for fid, lab in ((47, "unlabelled aromatic stream"), (139, "unlabelled aromatic stream")):
        rows.append(dict(name=f"ajftd_{fid}", label=lab, source="ajftd", fuel_id=fid, status="not used: no pathway label"))
    out = pd.DataFrame(rows)
    out.to_csv(D / "blendstocks.csv", index=False)
    print(out.to_string())


# ------------------------------------------------------------------------------------------------ R-W7-2 / 3
def _measured() -> pd.DataFrame:
    a = pd.read_csv(V5 / "crc_av3322_properties.csv")
    a = a[a.property.isin(D86_PROPS) & (a.method == "D86")]
    a = a.assign(source="av3322", fuel_id=a.fuel_id.astype(str))
    j = pd.read_csv(V5 / "ajftd_fuels_long.csv", low_memory=False)
    g = pd.read_csv(V5 / "ajftd_gcxgc_fuels.csv")
    j = j[j.property.isin(D86_PROPS) & j.fuel_id.isin(g.fuel_id)]
    n_nomethod = int(j.method.isna().sum())
    j = j[j.method.astype(str).str.replace(" ", "").str.contains("D86")]
    j = j.assign(source="ajftd", fuel_id=j.fuel_id.astype(str))
    m = pd.concat([a[["source", "fuel_id", "property", "value"]], j[["source", "fuel_id", "property", "value"]]])
    m["value"] = pd.to_numeric(m.value, errors="coerce")
    m = m.dropna(subset=["value"]).groupby(["source", "fuel_id", "property"]).value.mean().unstack()
    return m.rename(columns=D86_PROPS).reset_index(), n_nomethod


def cmd_margins():
    from fuelblend.eval.composition_route import kind_bins, map_bins
    reg, prov, reps = B._reg_prov_reps()
    meas, n_nomethod = _measured()
    t2 = pd.read_csv(ROOT / "results/v5/BASE/t2_fuels.csv")
    t2 = t2.assign(fuel_id=t2.fuel_id.astype(str))
    meas = meas.merge(t2[["source", "fuel_id", "group", "chosen_for_group"]], on=["source", "fuel_id"], how="left")
    n_dup = int((meas.chosen_for_group == False).sum())  # noqa: E712
    meas = meas[meas.chosen_for_group.fillna(True).astype(bool)]
    comp = {"av3322": pd.read_csv(V5 / "crc_av3322_composition.csv"), "ajftd": pd.read_csv(V5 / "ajftd_gcxgc_composition.csv")}
    rows, excl = [], dict(low_coverage=0, no_composition=0)
    for r in meas.itertuples():
        c = comp[r.source]
        sub = c[c.fuel_id.astype(str) == r.fuel_id]
        if sub.empty:
            excl["no_composition"] += 1
            continue
        mf = map_bins(kind_bins(r.source, sub), reps)
        tbs = np.array([B.tb_C(reg.smiles[i]) for i in mf.idx])
        ok = np.isfinite(tbs)
        cov = (1.0 - mf.unmapped_fraction) * float(mf.w[ok].sum())
        row = dict(source=r.source, fuel_id=r.fuel_id, group=r.group, coverage=cov)
        for q in ("T10", "T50", "T90", "FBP"):
            row[f"meas_{q}"] = getattr(r, q, np.nan)
        row["meas_T50_T10"] = row["meas_T50"] - row["meas_T10"]
        row["meas_T90_T10"] = row["meas_T90"] - row["meas_T10"]
        if cov < 0.90:
            excl["low_coverage"] += 1
            row["excluded"] = "coverage < 0.90"
            rows.append(row)
            continue
        _, rho = prov.lookup(mf.idx[ok][None, :], B.T15_K)
        vol = mf.w[ok] / rho[0]
        q = B.d86_quantities(tbs[ok], vol / vol.sum())
        row.update({f"pred_{k}": q[k] for k in ("T10", "T50", "T90", "FBP", "T50_T10", "T90_T10")}, tbp_spread=q["tbp_spread"],
                   monotone=q["monotone"], excluded="")
        rows.append(row)
    df = pd.DataFrame(rows)
    df.to_csv(D / "d86_fit_rows.csv", index=False, float_format="%.6g")
    use = df[df.excluded == ""]
    fit = B.fit_d86_margins(use)
    out = dict(margins={q: (v["certifying_margin"] if np.isfinite(v["certifying_margin"]) else "inf") for q, v in fit.items()},
               domain=dict(tbp_spread_min=float(use.tbp_spread.min()), tbp_spread_max=float(use.tbp_spread.max())),
               fit=fit, n_fuels_used=int(len(use)), excluded=excl, n_double_count_dropped=n_dup, ajftd_rows_without_method=n_nomethod,
               alpha=0.10, note="R-W7-2 declared in README; adopting into configs/uncertainty_v5.yaml is a user decision (not done)")
    json.dump(out, open(B.D86_MARGINS, "w", encoding="utf-8"), indent=1, default=float)
    print(json.dumps(out, indent=1, default=float))
    lit = {"n-decane": ("CCCCCCCCCC", 174.1), "n-dodecane": ("CCCCCCCCCCCC", 216.3), "butylcyclohexane": ("CCCCC1CCCCC1", 180.9),
           "isooctane": ("CC(C)CC(C)(C)C", 99.2), "1,2,4-trimethylbenzene": ("Cc1ccc(C)c(C)c1", 169.4),
           "decalin_trans": ("C1CC[C@H]2CCCC[C@@H]2C1", 187.3), "n-undecane": ("CCCCCCCCCCC", 195.9)}
    tb = pd.DataFrame([dict(compound=k, tb_antoine_C=B.tb_C(s), tb_literature_C=v, diff_C=B.tb_C(s) - v) for k, (s, v) in lit.items()])
    tb.to_csv(D / "tb_sanity.csv", index=False, float_format="%.3f")
    print(tb.to_string())


# ------------------------------------------------------------------------------------------------ R-W7-5 / R-W7-6
def make_pool(pool: str) -> B.BlendPool:
    pure = [B.pure_component(n) for n in POOL_C]
    if pool in ("P0", "C"):
        return B.BlendPool(pure)
    if pool == "P1":
        return B.BlendPool([B.load_blendstock("jetA_posf10325")] + pure)
    if pool == "P2":
        return B.BlendPool([B.load_blendstock("jetA_posf10325"), B.load_blendstock("hefa_came"), B.load_blendstock("atj_posf11498"),
                            B.pure_component("1,2,4-trimethylbenzene")])
    raise ValueError(pool)


def _rows(res, top, tag, spec):
    d = dict(design_id=f"{spec.replace(':', '_')}__{tag}", spec=spec, passed=res.passed, category=res.category,
             top_weights=json.dumps({k: round(v, 4) for k, v in top.items()}),
             weights=json.dumps({k: round(v, 5) for k, v in res.weights.items() if v > 0}), **{f"opt_{k}": v for k, v in res.optimiser.items()
                                                                                              if not isinstance(v, (dict, list))})
    c = [dict(design_id=d["design_id"], constraint=x.key, status=x.status, prediction=x.prediction, applied_margin=x.applied_margin,
              slack=x.slack, lo=x.lo, hi=x.hi, model=x.model, spec_verification=x.spec_verification, reason=x.reason) for x in res.constraints]
    return d, c


def run_one(pool, spec, d86, lam, tag, out_dir, n_starts, cap=None, max_iter=120):
    t0 = time.time()
    bp = make_pool(pool)
    models = B.build_models(bp)
    if d86:
        models = B.add_d86_models(models)
    t1 = time.time()
    res, top, info = B.design_pool(spec, bp, models, MARGINS, n_starts=n_starts, seed=42, d86=d86, lam=lam, synthetic_cap=cap, max_iter=max_iter)
    d, c = _rows(res, top, tag, spec)
    d.update(build_s=round(t1 - t0, 1), design_s=round(time.time() - t1, 1), pool=pool, n_compounds=len(bp.comp_ids), d86=d86)
    out_dir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame([d]).to_csv(out_dir / f"{d['design_id']}__designs.csv", index=False)
    pd.DataFrame(c).to_csv(out_dir / f"{d['design_id']}__constraints.csv", index=False)
    print(json.dumps(d, indent=1, default=str))
    for x in c:
        print(f"  {x['constraint']:30s} {x['status']:16s} pred {x['prediction']:.4g} slack {x['slack']:.4g} margin {x['applied_margin']}  {x['reason'][:160]}")
    return d, c


def cmd_compare(args):
    out = D / "compare"
    lam, tag = (B.LAMBDA_UNC, "poolC_lam025") if args.lam_on else (0.0, "poolC_lam0")
    run_one("C", args.spec, False, lam, tag, out, args.n_starts)


def cmd_design(args):
    cap = B.SYNTHETIC_CAP if args.pool == "P2" else None
    run_one(args.pool, args.spec, True, B.LAMBDA_UNC, args.pool, D / "designs" / args.pool, args.n_starts, cap=cap)


# ------------------------------------------------------------------------------------------------ counterfactual (ported from 6dc66f8 delta_smoke_pool.py)
VER_RE = r"spec limit verification: (?:secondary_sources_only|UNVERIFIED) \(an unverified limit is never certified\)"
NAPH = r"\(naphthalenes_vol_pct_max: by_construction \(pred 0\); " + VER_RE + r"\)"


def residual_reason(reason):
    r = "" if pd.isna(reason) else str(reason)
    r = re.sub(NAPH, "(naphthalenes by_construction)", r)
    r = re.sub(VER_RE, "", r)
    if "alt_min 18 met" in r:
        r = r.replace("limit violated after applying the margin", "")
    r = r.replace("alt_min 18 met but not certified (naphthalenes by_construction)", "")
    r = re.sub(r"[;\s]+", " ", r).strip(" ;")
    return r


def counterfactual(cc):
    rem = []
    for k, r in cc.iterrows():
        st = r.status
        if st in ("pass", "by_construction"):
            continue
        res = residual_reason(r.reason)
        if st == "uncertified" and float(r.slack) > 0 and res in ("", "alt_min 18 met but not certified"):
            continue
        rem.append(f"{k}={st}" + (f" [{res}]" if res else "") + (f" slack {float(r.slack):+.3f}" if pd.notna(r.slack) else ""))
    return len(rem) == 0, rem


def cmd_table():
    rows = []
    for f in sorted(list((D / "designs").glob("*/*__designs.csv")) + list((D / "compare").glob("*__designs.csv"))):
        d = pd.read_csv(f).iloc[0]
        c = pd.read_csv(str(f).replace("__designs.csv", "__constraints.csv")).set_index("constraint")
        sl = c.slack.astype(float)
        ok, rem = counterfactual(c)
        row = dict(d)
        row.update(min_slack=float(sl.min()), binding=",".join(k for k, v in sl.items() if v <= sl.min() + 0.003),
                   would_certify_if_specs_verified=ok, counterfactual_blockers=" | ".join(rem),
                   n_pass=int((c.status == "pass").sum()), n_fail=int((c.status == "fail").sum()), n_uncertified=int((c.status == "uncertified").sum()))
        rows.append(row)
    out = pd.DataFrame(rows)
    out.to_csv(D / "r_designs.csv", index=False)
    pd.set_option("display.width", 250)
    pd.set_option("display.max_colwidth", 300)
    print(out[["design_id", "pool", "passed", "min_slack", "binding", "would_certify_if_specs_verified", "n_pass", "n_fail", "n_uncertified",
               "top_weights"]].to_string())
    for _, r in out.iterrows():
        print(r.design_id, "blockers:", r.counterfactual_blockers[:1500])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd")
    ap.add_argument("--pool", default="P1")
    ap.add_argument("--spec", default="astm:jet_a")
    ap.add_argument("--n-starts", type=int, default=10)
    ap.add_argument("--lam-on", action="store_true")
    a = ap.parse_args()
    D.mkdir(parents=True, exist_ok=True)
    {"blendstocks": lambda: cmd_blendstocks(), "margins": lambda: cmd_margins(), "compare": lambda: cmd_compare(a),
     "design": lambda: cmd_design(a), "table": lambda: cmd_table()}[a.cmd]()


if __name__ == "__main__":
    main()
