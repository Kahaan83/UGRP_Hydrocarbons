"""delta_d86_rayleigh (results/v5/BASE/delta_d86_rayleigh/README.md): Rayleigh D86 model vs the WS7 Engler route.

    python scripts/v5/base/delta_d86_rayleigh.py margins   # R-DR-2 -> d86_rayleigh_fit_rows.csv, d86_rayleigh_margins.json, compare_margins.csv
    python scripts/v5/base/delta_d86_rayleigh.py recipes   # R-DR-3 -> stored_recipes.csv
    python scripts/v5/base/delta_d86_rayleigh.py design --pool P2 --spec astm:jet_a   # R-DR-4
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts" / "v5" / "base"))

from fuelblend.design import blendstock_v5 as B  # noqa: E402

D = ROOT / "results" / "v5" / "BASE" / "delta_d86_rayleigh"
WS7 = ROOT / "results" / "v5" / "BASE" / "delta_ws7_design_realism"
SBO = ROOT / "results" / "v5" / "BASE" / "delta_smoke_branch_opt"
RAY_MARGINS = D / "d86_rayleigh_margins.json"
QS = ("T10", "T50", "T90", "FBP", "T50_T10", "T90_T10")
CONSTRAINED = ("T10", "FBP", "T50_T10", "T90_T10")


def cmd_margins():
    import delta_ws7_design_realism as W
    from fuelblend.eval.composition_route import kind_bins, map_bins
    reg, prov, reps = B._reg_prov_reps()
    eng = pd.read_csv(WS7 / "d86_fit_rows.csv", dtype={"fuel_id": str})
    meas, _ = W._measured()
    comp = {"av3322": pd.read_csv(W.V5 / "crc_av3322_composition.csv"), "ajftd": pd.read_csv(W.V5 / "ajftd_gcxgc_composition.csv")}
    rows = []
    t0 = time.time()
    for r in eng.itertuples():
        row = {k: getattr(r, k) for k in ("source", "fuel_id", "group", "coverage") + tuple(f"meas_{q}" for q in QS)}
        if isinstance(r.excluded, str) and r.excluded:
            row["excluded"] = r.excluded
            rows.append(row)
            continue
        c = comp[r.source]
        sub = c[c.fuel_id.astype(str) == str(r.fuel_id)]
        mf = map_bins(kind_bins(r.source, sub), reps)
        smi = [reg.smiles[i] for i in mf.idx]
        tbs = np.array([B.tb_C(s) for s in smi])
        ok = np.isfinite(tbs)
        _, rho = prov.lookup(mf.idx[ok][None, :], B.T15_K)
        smi_ok = [s for s, k in zip(smi, ok) if k]
        q = B.rayleigh_quantities_smiles(smi_ok, mf.w[ok] / mf.w[ok].sum(), rho[0])
        row.update({f"pred_{k}": q[k] for k in QS}, tbp_spread=q["tbp_spread"], monotone=q["monotone"], excluded="")
        rows.append(row)
    print(f"predicted {len(rows)} fuels in {time.time() - t0:.1f} s")
    df = pd.DataFrame(rows)
    df.to_csv(D / "d86_rayleigh_fit_rows.csv", index=False, float_format="%.6g")
    use = df[df.excluded == ""]
    fit = B.fit_d86_margins(use)
    e_fit = json.load(open(WS7 / "d86_margins.json", encoding="utf-8"))["fit"]
    comp_rows, adopt = [], {}
    for qn in QS:
        f, e = fit[qn], e_fit[qn]
        pref = None
        if qn in CONSTRAINED:
            pref = bool(np.isfinite(f["certifying_margin"]) and f["certifying_margin"] < e["certifying_margin"] and f["mae"] <= e["mae"])
            adopt[qn] = "rayleigh" if pref else "engler"
        comp_rows.append(dict(quantity=qn, direction=f["direction"], n=f["n"], engler_margin=e["certifying_margin"], rayleigh_margin=f["certifying_margin"],
                              engler_mae=e["mae"], rayleigh_mae=f["mae"], engler_bias=e["bias"], rayleigh_bias=f["bias"], mae_constant=f["mae_constant"],
                              rayleigh_per_source=json.dumps(f["per_source"]), rayleigh_preferred=pref))
    cm = pd.DataFrame(comp_rows)
    cm.to_csv(D / "compare_margins.csv", index=False, float_format="%.4f")
    print(cm.drop(columns=["rayleigh_per_source"]).to_string())
    # per-source MAE (reported)
    for src, g in use.groupby("source"):
        print(src, {q: round(float(np.mean(np.abs(g[f"pred_{q}"] - g[f"meas_{q}"]))), 2) for q in QS})
    out = dict(margins={q: (v["certifying_margin"] if np.isfinite(v["certifying_margin"]) else "inf") for q, v in fit.items()},
               domain=dict(tbp_spread_min=float(use.tbp_spread.min()), tbp_spread_max=float(use.tbp_spread.max())),
               fit=fit, n_fuels_used=int(len(use)), adoption_R_DR_2=adopt, source="delta_d86_rayleigh/d86_rayleigh_margins.json",
               alpha=0.10, note="R-DR-2 declared in README; designating the Rayleigh D86 model is a user decision (not done)")
    json.dump(out, open(RAY_MARGINS, "w", encoding="utf-8"), indent=1, default=float)
    print("adoption:", adopt)


def mixed_d86() -> tuple[dict, dict]:
    """(d86 margins dict, per-quantity model choice) under the R-DR-2 adoption rule; domain = the WS7 (Engler) domain (the same TBP rule)."""
    e = B.load_d86_margins()
    r = json.load(open(RAY_MARGINS, encoding="utf-8"))
    choice = r["adoption_R_DR_2"]
    margins = {q: (r["margins"][q] if choice.get(q) == "rayleigh" else e["margins"][q]) for q in B.D86_PROPS}
    return dict(margins=margins, domain=e["domain"], source="delta_d86_rayleigh/d86_rayleigh_margins.json"), choice


def cmd_recipes():
    import delta_smoke_branch_opt as SB
    reg, prov, reps = B._reg_prov_reps()
    tab = B._table_smiles()
    e_m = B.load_d86_margins()["margins"]
    r_m = json.load(open(RAY_MARGINS, encoding="utf-8"))["margins"]

    def rho15(s):
        j = reg.ik2i.get(B._ik(s))
        if j is None or not prov.valid[j]:
            return float("nan")
        _, r = prov.lookup(np.array([[j]]), B.T15_K)
        return float(r[0][0])

    out, pools = [], {}
    for f in (WS7 / "r_designs.csv", SBO / "r_designs.csv"):
        d = pd.read_csv(f)
        for row in d.itertuples():
            pool = "P0" if row.pool == "C" else row.pool
            if pool not in pools:
                pools[pool] = SB.make_pool(pool)
            bp = pools[pool]
            top = json.loads(row.top_weights)
            u = np.array([top.get(m.name, 0.0) for m in bp.members], float)
            wfull = bp.expand(u)
            keep = wfull > 1e-12
            smi = [tab.get(n, n) for n, k in zip(bp.comp_ids, keep) if k]
            w = wfull[keep] / wfull[keep].sum()
            rho = np.array([rho15(s) for s in smi])
            tbs = np.array([B.tb_C(s) for s in smi])
            if not (np.all(np.isfinite(rho)) and np.all(np.isfinite(tbs))):
                print("skip (missing rho/Tb)", row.design_id)
                continue
            vol = w / rho
            qe = B.d86_quantities(tbs, vol / vol.sum())
            qr = B.rayleigh_quantities_smiles(smi, w, rho)
            rec = dict(file=f.parent.name, design_id=row.design_id, pool=pool, spec=row.spec)
            for q in QS:
                rec[f"engler_{q}"] = qe[q]
                rec[f"rayleigh_{q}"] = qr[q]
            for q, (canon, side) in B.D86_PROPS.items():
                lim = {"T10": 205.0, "FBP": 300.0, "T50_T10": 15.0, "T90_T10": 40.0}[q]
                for mdl, qq, mm in (("engler", qe, e_m), ("rayleigh", qr, r_m)):
                    m = float("inf") if isinstance(mm[q], str) else float(mm[q])
                    rec[f"{mdl}_slack_{q}"] = (lim - qq[q] - m) if side == "hi" else (qq[q] - m - lim)
            rec["tbp_spread"] = qe["tbp_spread"]
            out.append(rec)
    df = pd.DataFrame(out)
    df.to_csv(D / "stored_recipes.csv", index=False, float_format="%.3f")
    cols = ["file", "design_id", "engler_T10", "rayleigh_T10", "engler_T50_T10", "rayleigh_T50_T10", "engler_T90_T10", "rayleigh_T90_T10",
            "engler_slack_T90_T10", "rayleigh_slack_T90_T10"]
    pd.set_option("display.width", 250)
    print(df[cols].round(1).to_string())


def cmd_design(args):
    import delta_ws7_design_realism as W
    d86, choice = mixed_d86()
    out_dir = D / "designs" / args.pool
    out_dir.mkdir(parents=True, exist_ok=True)
    bp = W.make_pool(args.pool)
    models = W.LR.build_route_models(bp.comp_ids) if hasattr(W, "LR") else None
    raise SystemExit("design: use run_design() wired in a later step")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["margins", "recipes", "design"])
    ap.add_argument("--pool")
    ap.add_argument("--spec")
    a = ap.parse_args()
    {"margins": lambda: cmd_margins(), "recipes": lambda: cmd_recipes(), "design": lambda: cmd_design(a)}[a.cmd]()


if __name__ == "__main__":
    main()
