"""delta_pool_bicyclohexyl (results/v5/BASE/delta_pool_bicyclohexyl/README.md).

    python scripts/v5/base/delta_pool_bicyclohexyl.py data                          # R-BC-1 -> data/processed/v5/flash_registry_additions.csv, r_bc_1.csv
    python scripts/v5/base/delta_pool_bicyclohexyl.py design --spec astm:jet_a      # R-BC-2 (sets FUELBLEND_FLASH_REGISTRY_ADDITIONS=1 itself)
    python scripts/v5/base/delta_pool_bicyclohexyl.py checks                        # R-BC-2 Rayleigh check + R-BC-3 freezing sensitivity
    python scripts/v5/base/delta_pool_bicyclohexyl.py table
"""
from __future__ import annotations

import argparse
import json
import math
import os
import sys
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "results" / "v5" / "BASE" / "delta_pool_bicyclohexyl"
ADD = ROOT / "data" / "processed" / "v5" / "flash_registry_additions.csv"
BASE_REG = ROOT / "data" / "processed" / "v2" / "compounds_anchored_vp.csv"
NEW = {"bicyclohexyl": ("WVIIMZNLDWSIRH-UHFFFAOYSA-N", "92-51-3", 239.0), "decalin_cis": ("NNBZCPXTIHJBJL-AOOOYVTPSA-N", "493-01-6", 195.8)}
P_ATM = 101.325


def cmd_data():
    from thermo import Chemical
    from scipy.optimize import least_squares
    sys.path.insert(0, str(ROOT / "src"))
    base = pd.read_csv(BASE_REG)
    meas = pd.read_csv(ROOT / "data" / "unified" / "measurements.csv", low_memory=False)
    rows, rep = [], []
    T = np.linspace(273.15, 373.15, 21)
    for name, (ik, cas, tb_lit) in NEW.items():
        b = base[base.inchikey == ik]
        if b.empty:
            raise SystemExit(f"{name}: no base registry row")
        b = b.iloc[0].copy()
        c = Chemical(cas)
        vp = c.VaporPressure
        y = np.log10(np.array([vp.T_dependent_property(t) for t in T]) / 1000.0)
        tb_guess = c.Tb or (tb_lit + 273.15)
        B0 = (tb_guess - 273.15 + 230.0) * 8.0
        x0 = np.array([y[-1] + B0 / (T[-1] - 50.0), B0, -50.0])

        def res(p):
            return p[0] - p[1] / (T + p[2]) - y
        fit = least_squares(res, x0, bounds=([-np.inf, 1.0, -250.0], [np.inf, 1e5, 100.0]), xtol=1e-14, ftol=1e-14)
        A, Bc, Cc = fit.x
        tb = Bc / (A - math.log10(P_ATM)) - Cc - 273.15
        lfl = float(b.lfl_vol_pct)
        # pure Le Chatelier flash: P_sat(T) = LFL/100 * P
        tf = Bc / (A - math.log10(lfl / 100.0 * P_ATM)) - Cc - 273.15
        fm = meas[(meas.kind == "pure") & (meas.property == "flash_point_C") & meas.comp_inchikeys.str.contains(ik, regex=False, na=False)]
        b["antoine_A"], b["antoine_B"], b["antoine_C"] = A, Bc, Cc
        b["antoine_Tlo"], b["antoine_Thi"] = 273.15, 373.15
        b["antoine_source"] = "anchored_thermo_ambient_v5add"
        b["antoine_covers_window"] = True
        rows.append(b)
        gate = abs(tb - tb_lit) <= 10.0
        rep.append(dict(compound=name, inchikey=ik, thermo_method=vp.method, A=A, B=Bc, C=Cc, max_abs_log10_resid=float(np.abs(fit.fun).max()),
                        tb_antoine_C=tb, tb_literature_C=tb_lit, tb_diff_C=tb - tb_lit, gate_tb_within_10C=gate, lfl_vol_pct=lfl,
                        lfl_source=b.get("lfl_source"), flash_pure_pred_C=tf, flash_meas_C=";".join(f"{v:.2f}" for v in fm.value),
                        flash_meas_median_C=float(fm.value.median()) if len(fm) else float("nan")))
    rep = pd.DataFrame(rep)
    OUT.mkdir(parents=True, exist_ok=True)
    rep.to_csv(OUT / "r_bc_1.csv", index=False, float_format="%.6g")
    add = pd.DataFrame(rows)
    add = add[[c for c in base.columns]]
    keep = add[add.inchikey.isin(rep.inchikey[rep.gate_tb_within_10C])]
    ADD.parent.mkdir(parents=True, exist_ok=True)
    keep.to_csv(ADD, index=False)
    pd.set_option("display.width", 250)
    print(rep.to_string())
    print(f"wrote {len(keep)} rows to {ADD}")


def _env():
    os.environ["FUELBLEND_FLASH_REGISTRY_ADDITIONS"] = "1"
    sys.path.insert(0, str(ROOT / "src"))
    sys.path.insert(0, str(ROOT / "scripts" / "v5" / "base"))
    warnings.simplefilter("ignore")


def pool_members(pool="P0y"):
    rep = pd.read_csv(OUT / "r_bc_1.csv")
    import delta_smoke_branch_opt as SB
    if pool == "P0yb":                       # R-BC-4: cis-decalin dropped (MW tier makes flash / NHOC unavailable)
        return SB.POOL_X["P0xb"] + ["bicyclohexyl"]
    return SB.POOL_X["P0xb"] + [c for c in rep.compound[rep.gate_tb_within_10C]]


def cmd_design(a):
    _env()
    from fuelblend.design import blendstock_v5 as B
    import delta_smoke_branch_opt as SB
    W = SB.W
    from fuelblend.baselines.flash_unifac import _lookup
    for n in NEW:
        print(n, "registry entry:", _lookup(B._table_smiles()[n], None))
    t0 = time.time()
    bp = B.BlendPool([B.pure_component(n) for n in pool_members(a.pool)])
    models = B.add_d86_models(B.build_models(bp))
    t1 = time.time()
    res, top, info = B.design_pool(a.spec, bp, models, W.MARGINS, n_starts=a.n_starts, seed=42, d86=True, lam=B.LAMBDA_UNC,
                                   synthetic_cap=None, alt_standin=True)
    d, c = W._rows(res, top, f"{a.pool}_standin", a.spec)
    d.update(build_s=round(t1 - t0, 1), design_s=round(time.time() - t1, 1), pool=a.pool, n_compounds=len(bp.comp_ids), d86=True,
             alt_standin=True, registry_additions=True)
    out = OUT / "designs" / a.pool
    out.mkdir(parents=True, exist_ok=True)
    pd.DataFrame([d]).to_csv(out / f"{d['design_id']}__designs.csv", index=False)
    pd.DataFrame(c).to_csv(out / f"{d['design_id']}__constraints.csv", index=False)
    print(json.dumps(d, indent=1, default=str))
    for x in c:
        print(f"  {x['constraint']:30s} {x['status']:16s} pred {x['prediction']:.5g} slack {x['slack']:+.5f} margin {x['applied_margin']}  "
              f"{str(x['reason'])[:150]}")


def cmd_checks():
    _env()
    from fuelblend.design import blendstock_v5 as B
    from fuelblend.eval.composition_route import build_tmap_verified
    from fuelblend.design import engine_v4 as E
    from fuelblend.baselines.distillation_rayleigh import rayleigh_quantities  # noqa: F401
    rm = json.load(open(ROOT / "results" / "v5" / "BASE" / "delta_d86_rayleigh" / "d86_rayleigh_margins.json", encoding="utf-8"))["margins"]
    reg, prov, _ = B._reg_prov_reps()
    tab = B._table_smiles()
    tmap = build_tmap_verified()
    rows = []
    for f in sorted((OUT / "designs").glob("*/*__designs.csv")):
        d = pd.read_csv(f).iloc[0]
        c = pd.read_csv(str(f).replace("__designs.csv", "__constraints.csv")).set_index("constraint")
        wd = json.loads(d.weights)
        names = [n for n in wd if wd[n] > 0]
        smi = [tab.get(n, n) for n in names]
        w = np.array([wd[n] for n in names], float)
        w /= w.sum()
        rho = []
        for s in smi:
            j = reg.ik2i.get(B._ik(s))
            _, r = prov.lookup(np.array([[j]]), B.T15_K)
            rho.append(float(r[0][0]))
        q = B.rayleigh_quantities_smiles(smi, w, np.array(rho))
        flags = []
        if d.passed and (q["T90_T10"] < 40 + float(rm["T90_T10"]) or q["T50_T10"] < 15 + float(rm["T50_T10"])):
            flags.append("d86_model_dependent")
        # R-BC-3: ideal SvL first-to-freeze with bicyclohexyl dH = 13.94 kJ/mol, vs the governing 6.86
        mw = np.array([B._mw(s) for s in smi])
        x = (w / mw) / (w / mw).sum()
        iks = [B._ik(s) for s in smi]

        def svl(dh_override):
            best = -1e9
            for xi, ik in zip(x, iks):
                if ik not in tmap or xi <= 0:
                    continue
                tm, dh = tmap[ik]["tm_K"], tmap[ik]["dh_fus_J_mol"]
                if ik == NEW["bicyclohexyl"][0] and dh_override is not None:
                    dh = dh_override
                Tsat = 1.0 / (1.0 / tm - 8.314462618 * math.log(xi) / dh)
                best = max(best, Tsat)
            return best - 273.15
        fz_gov, fz_tot = svl(None), svl(13940.0)
        lim = {"astm:jet_a": -40.0, "bis:jet_a1": -47.0}[d.spec]
        fm = 16.975450975521966
        if d.passed and fz_tot + fm > lim:
            flags.append("freezing_dH_choice_dependent")
        rows.append(dict(design_id=d.design_id, spec=d.spec, passed=d.passed, weights=d.weights,
                         judge_T90_T10=c.loc["distillation_T90_minus_T10_C", "prediction"] if "distillation_T90_minus_T10_C" in c.index else np.nan,
                         rayleigh_T10=q["T10"], rayleigh_T50_T10=q["T50_T10"], rayleigh_T90_T10=q["T90_T10"], rayleigh_FBP=q["FBP"],
                         rayleigh_slack_T90_T10=q["T90_T10"] - float(rm["T90_T10"]) - 40.0, rayleigh_slack_T50_T10=q["T50_T10"] - float(rm["T50_T10"]) - 15.0,
                         svl_first_freeze_C_governing=fz_gov, svl_first_freeze_C_bch_total_dH=fz_tot, freezing_limit_C=lim, freezing_margin_C=fm,
                         flags=",".join(flags)))
    out = pd.DataFrame(rows)
    out.to_csv(OUT / "r_bc_checks.csv", index=False, float_format="%.4f")
    pd.set_option("display.width", 250)
    print(out.drop(columns=["weights"]).to_string())


def cmd_flash5(a):
    """R-BC-5: mixture flash of the P0y-b optima with bicyclohexyl's Antoine A shifted so its pure Le Chatelier flash = 74.0 C.
    One process per variant (``--variant overlay|shift74``) because registry lookups are cached."""
    _env()
    from fuelblend.baselines import flash_unifac as fu
    add = pd.read_csv(ADD)
    ik = NEW["bicyclohexyl"][0]
    i = add.index[add.inchikey == ik][0]
    A, Bc, Cc, lfl = (float(add.loc[i, k]) for k in ("antoine_A", "antoine_B", "antoine_C", "lfl_vol_pct"))
    if a.variant == "shift74":
        A = math.log10(lfl / 100.0 * P_ATM) + Bc / (74.0 + 273.15 + Cc)
        add.loc[i, "antoine_A"] = A
        tmp = OUT / "r_bc_5_additions_shift74.csv"
        add.to_csv(tmp, index=False)
        fu.REGISTRY_ADDITIONS_CSV = tmp
    fu._default_registry.cache_clear()
    from fuelblend.design import blendstock_v5 as B
    tb = Bc / (A - math.log10(P_ATM)) - Cc - 273.15
    pure_lc = Bc / (A - math.log10(lfl / 100.0 * P_ATM)) - Cc - 273.15
    bp1 = B.BlendPool([B.pure_component("bicyclohexyl")])
    m1 = B.build_models(bp1)
    pure_model = float(m1["props"]["flash"].fn(list(m1["pool"]), np.array([1.0]), None))
    bp = B.BlendPool([B.pure_component(n) for n in pool_members("P0yb")])
    models = B.build_models(bp)
    rows = []
    for f in sorted((OUT / "designs" / "P0yb").glob("*__designs.csv")):
        d = pd.read_csv(f).iloc[0]
        wd = json.loads(d.weights)
        u = np.array([wd.get(n, 0.0) for n in bp.names], float)
        fl = float(models["props"]["flash"].fn(list(models["pool"]), bp.expand(u), None))
        rows.append(dict(variant=a.variant, design_id=d.design_id, antoine_A=A, tb_antoine_C=tb, pure_lc_flash_C=pure_lc,
                         pure_model_flash_C=pure_model, mixture_flash_C=fl, margin_C=2.206932808388615, limit_C=38.0,
                         bound_C=fl - 2.206932808388615, flash_data_dependent=bool(fl - 2.206932808388615 < 38.0)))
    out = pd.DataFrame(rows)
    fn = OUT / "r_bc_5.csv"
    old = pd.read_csv(fn) if fn.exists() else pd.DataFrame()
    if len(old):
        old = old[old.variant != a.variant]
    pd.concat([old, out], ignore_index=True).to_csv(fn, index=False, float_format="%.5f")
    pd.set_option("display.width", 250)
    print(out.to_string())


def cmd_table():
    rows = []
    for f in sorted((OUT / "designs").glob("*/*__designs.csv")):
        d = pd.read_csv(f).iloc[0]
        c = pd.read_csv(str(f).replace("__designs.csv", "__constraints.csv")).set_index("constraint")
        sl = c.slack.astype(float)
        row = dict(d)
        row.update(min_slack=float(sl.min()), binding=",".join(k for k, v in sl.items() if v <= sl.min() + 0.003),
                   n_pass=int((c.status == "pass").sum()), n_fail=int((c.status == "fail").sum()),
                   n_uncertified=int((c.status == "uncertified").sum()), failing=",".join(c.index[c.status == "fail"]),
                   uncertified=",".join(c.index[c.status == "uncertified"]))
        rows.append(row)
    out = pd.DataFrame(rows)
    out.to_csv(OUT / "r_designs.csv", index=False)
    pd.set_option("display.width", 250)
    print(out[["design_id", "passed", "min_slack", "binding", "n_pass", "n_fail", "n_uncertified", "failing", "uncertified", "top_weights"]].to_string())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["data", "design", "checks", "table", "flash5"])
    ap.add_argument("--variant", default="overlay", choices=["overlay", "shift74"])
    ap.add_argument("--spec", default="astm:jet_a")
    ap.add_argument("--n-starts", type=int, default=10)
    ap.add_argument("--pool", default="P0y", choices=["P0y", "P0yb"])
    a = ap.parse_args()
    {"data": cmd_data, "design": lambda: cmd_design(a), "checks": cmd_checks, "table": cmd_table, "flash5": lambda: cmd_flash5(a)}[a.cmd]()


if __name__ == "__main__":
    main()
