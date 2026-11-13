"""R-SC-3 table builder: reads designs/<run>/{designs,constraints}.csv of delta_smoke_cert, writes r_sc_3.csv."""
import json
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "src"))
D = ROOT / "results/v5/BASE/delta_smoke_cert"
RUNS = ["v5pool", "guard4_v5pool", "ext_n30", "ext_n8", "poolA_toluene_isooctane", "poolB_hexylbenzene_isooctane"]


def main():
    from fuelblend.baselines.sooting import smoke_point_from_bulk
    rows = []
    for run in RUNS:
        d = pd.read_csv(D / "designs" / run / "designs.csv")
        c = pd.read_csv(D / "designs" / run / "constraints.csv")
        for _, r in d.iterrows():
            cc = c[c.design_id == r.design_id].set_index("constraint")
            row = dict(run=run, design=r.design_id.split("__")[0], category=r.category, passed=r.passed)
            def g(k, col):
                return cc.loc[k, col] if k in cc.index else None
            ar = g("aromatics_vol_pct", "prediction")
            sm = g("smoke_point_mm", "prediction")
            mg = g("smoke_point_mm", "applied_margin")
            row.update(arom_vol_pct=ar, smoke_pred=sm, smoke_margin=mg,
                       max_margin_certifies=(sm - 18.0) if sm is not None and pd.notna(sm) else None,
                       smoke_status=g("smoke_point_mm", "status"), smoke_reason=g("smoke_point_mm", "reason"),
                       smoke_slack=g("smoke_point_mm", "slack"))
            sr = str(row["smoke_reason"])
            row["class_mix_gate"] = "outside_range" if "class_mix_outside" in sr else ("inside" if pd.notna(mg) else "n/a")
            nk = [k for k in cc.index if "naphth" in k]
            row["naphthalenes"] = (f"{nk[0]}={cc.loc[nk[0], 'prediction']:.2f} {cc.loc[nk[0], 'status']}" if nk else "no constraint row (by construction 0 in the pool)")
            for k, nm in [("freezing_point_C", "freezing"), ("flash_point_C", "flash"), ("density_15C_kg_m3", "density"),
                          ("kin_visc_m20C_cSt", "visc_m20"), ("nhoc_MJ_kg", "nhoc"), ("aromatics_vol_pct", "arom_vol"), ("aromatics_wt_pct", "arom_wt")]:
                row[nm + "_status"] = g(k, "status")
                row[nm + "_slack"] = g(k, "slack")
            bad = [(k, s) for k, s in zip(cc.index, cc.status) if s not in ("pass", "by_construction")]
            row["non_pass"] = "; ".join(f"{k}={s}" for k, s in bad)
            only_ver = [k for k, s, v in zip(cc.index, cc.status, cc.spec_verification.fillna("")) if s not in ("pass", "by_construction")]
            row["blocked_only_by_secondary_sources_only"] = bool(
                len(only_ver) > 0 and all(cc.loc[k, "status"] == "uncertified" and "secondary_sources_only" in str(cc.loc[k, "spec_verification"]) + str(cc.loc[k, "reason"]) for k in only_ver))
            row["blocked_only_by_secondary_sources_only"] = bool(row["blocked_only_by_secondary_sources_only"] and "class_mix_outside" not in str(row["smoke_reason"]) and row["smoke_status"] != "fail")
            row["graded_objective"] = r.get("graded_objective")
            row["best_min_slack"] = r.get("best_min_slack")
            row["weights"] = r.weights
            rows.append(row)
    out = pd.DataFrame(rows)
    out.to_csv(D / "r_sc_3.csv", index=False)
    pd.set_option("display.width", 250)
    pd.set_option("display.max_colwidth", 200)
    print(out[["run", "design", "category", "arom_vol_pct", "smoke_pred", "smoke_margin", "max_margin_certifies", "smoke_status", "class_mix_gate",
               "freezing_status", "freezing_slack", "flash_status", "density_status", "visc_m20_status", "nhoc_status", "blocked_only_by_secondary_sources_only"]].to_string())
    print(out[["run", "design", "non_pass", "smoke_reason", "naphthalenes"]].to_string())


if __name__ == "__main__":
    main()
