"""delta_smoke_pool table builder: reads designs/poolC*/{designs,constraints}.csv, writes r_pool.csv.
The column would_certify_if_specs_verified is a COUNTERFACTUAL reading (rule declared in the README), never a certification."""
import re
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
D = ROOT / "results/v5/BASE/delta_smoke_pool"
RUNS = ["poolC_isooctane_tmb"]
VER_RE = r"spec limit verification: (?:secondary_sources_only|UNVERIFIED) \(an unverified limit is never certified\)"
NAPH = r"\(naphthalenes_vol_pct_max: by_construction \(pred 0\); " + VER_RE + r"\)"


def residual_reason(reason):
    """reason text left after removing the clauses that only say 'spec limit is secondary_sources_only'."""
    r = "" if pd.isna(reason) else str(reason)
    r = re.sub(NAPH, "(naphthalenes by_construction)", r)
    r = re.sub(VER_RE, "", r)
    if "alt_min 18 met" in r:   # primary limit violated with margin is the judge encoding of "carried by the D1655 alternative branch"
        r = r.replace("limit violated after applying the margin", "")
    r = r.replace("alt_min 18 met but not certified (naphthalenes by_construction)", "")
    r = re.sub(r"[;\s]+", " ", r).strip(" ;")
    return r


def counterfactual(cc):
    """(bool, list of remaining blockers)"""
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


def main():
    rows = []
    for run in RUNS:
        d = pd.read_csv(D / "designs" / run / "designs.csv")
        c = pd.read_csv(D / "designs" / run / "constraints.csv")
        for _, r in d.iterrows():
            cc = c[c.design_id == r.design_id].set_index("constraint")
            sl = cc["slack"].astype(float)
            ok, rem = counterfactual(cc)
            row = dict(run=run, design=r.design_id.split("__")[0], category=r.category, passed=r.passed,
                       weights=r.weights, min_slack_constraints=float(sl.min()),
                       binding=",".join(k for k, v in sl.items() if v <= sl.min() + 0.003),
                       would_certify_if_specs_verified=ok, counterfactual_blockers=" | ".join(rem))
            for k, v in cc.iterrows():
                row[f"{k}__status"] = v.status
                row[f"{k}__slack"] = v.slack
                row[f"{k}__margin"] = v.applied_margin
                row[f"{k}__pred"] = v.prediction
            row["non_pass_reasons"] = " || ".join(f"{k}: {v.status}: {v.reason}" for k, v in cc.iterrows() if v.status not in ("pass", "by_construction"))
            rows.append(row)
    out = pd.DataFrame(rows)
    out.to_csv(D / "r_pool.csv", index=False)
    pd.set_option("display.width", 250)
    pd.set_option("display.max_colwidth", 400)
    for _, r in out.iterrows():
        print("=====", r.run, r.design, r.category, "passed", r.passed, "min slack", round(r.min_slack_constraints, 4), "binding", r.binding)
        print(" weights", r.weights)
        for k in [c for c in out.columns if c.endswith("__status")]:
            n = k[:-8]
            print(f"  {n:22s} {str(r[k]):16s} slack {r[n+"__slack"]} margin {r[n+'__margin']}  pred {r[n+'__pred']}")
        print(" counterfactual:", r.would_certify_if_specs_verified, r.counterfactual_blockers)
        print(" reasons:", r.non_pass_reasons)


if __name__ == "__main__":
    sys.exit(main())
