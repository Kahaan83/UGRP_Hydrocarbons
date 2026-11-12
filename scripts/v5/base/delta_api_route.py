"""v5 delta_api_route (Task 5): the re-routed API (predict_fuel) against the BASE scoring route on the T2 fuels, plus a 60-component demo.

Reads (declared in results/v5/BASE/delta_api_route/README.md):
  R-API-1  flash: API designated == rows_t2_flash.pred_lfl_overlay, variants['v4_designated'] == pred_designated, 1e-6 C.
  R-API-2  viscosity (kinematic) and density: API == refutas_nu_safe on the mapped fuel at the 9 BASE nominal temperatures, 1e-9 relative.
  R-API-3  flash, net heat, aromatics (wt, vol): bit-identical across T_K in {233.15, 288.15, 313.15}.
  R-API-4  pending Task 6 (design recipes); recorded as pending, not evaluated here.
Outputs only into results/v5/BASE/delta_api_route/.
"""
from __future__ import annotations

import math
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "src"))
from fuelblend.api import predict_route as pr  # noqa: E402
from fuelblend.eval.composition_route import kind_bins, map_bins, refutas_nu_safe  # noqa: E402

BASE = ROOT / "results" / "v5" / "BASE"
OUT = BASE / "delta_api_route"
V5 = ROOT / "data" / "real_fuels" / "v5"
NOMINAL = [233.15, 253.15, 273.15, 288.15, 293.15, 313.15, 333.15, 353.15, 373.15]      # cells_t2_bulk.NOMINAL
R3_T = [233.15, 288.15, 313.15]
DEMO_T = [233.15, 253.15, 288.15, 313.15]
TOL_FLASH, TOL_REL = 1e-6, 1e-9


def same(a, b) -> bool:
    """Bit-identical, with None == None and NaN == NaN."""
    if a is None or b is None:
        return a is None and b is None
    if isinstance(a, float) and isinstance(b, float) and math.isnan(a) and math.isnan(b):
        return True
    return a == b


def _missing(a) -> bool:
    return a is None or not np.isfinite(a)


def rel(a, b) -> float:
    """Relative difference; 0 when both missing/NaN; inf when only one is."""
    if _missing(a) and _missing(b):
        return 0.0
    if _missing(a) or _missing(b):
        return math.inf
    return abs(a - b) / abs(b) if b != 0 else abs(a - b)


def absd(a, b) -> float:
    if _missing(a) and _missing(b):
        return 0.0
    if _missing(a) or _missing(b):
        return math.inf
    return abs(a - b)


def main() -> None:
    warnings.simplefilter("ignore")
    OUT.mkdir(parents=True, exist_ok=True)
    ctx = pr.get_ctx()
    reg, reps, prov = ctx.reg, ctx.reps, ctx.prov
    t = pd.read_csv(BASE / "t2_targets.csv")
    t = t[t.property.isin(["flash_point_C", "viscosity_kin_cSt", "density_kg_m3"])]
    fu = pd.read_csv(BASE / "t2_fuels.csv").set_index("entity")
    flash = pd.read_csv(BASE / "rows_t2_flash.csv").set_index("entity")
    comps = {"av3322": pd.read_csv(V5 / "crc_av3322_composition.csv"), "fuellib": pd.read_csv(V5 / "fuellib_gcxgc_compositions.csv"),
             "ajftd": pd.read_csv(V5 / "ajftd_gcxgc_composition.csv")}
    bulk = pd.read_csv(BASE / "rows_t2_bulk.csv")
    bulk = bulk[(bulk.variant == "designated") & bulk.property.isin(["viscosity_kin_cSt", "density_kg_m3"])].drop_duplicates(["entity", "property", "T_K"])
    bk = {(r.entity, r.property, round(r.T_K, 2)): r.pred for r in bulk.itertuples()}

    rows = []
    for ent, g in t.groupby("entity", sort=False):
        key = fu.loc[ent, "source"]
        fid = int(g.fuel_id.iloc[0]) if key == "ajftd" else g.fuel_id.iloc[0]
        c = comps[key]
        mf = map_bins(kind_bins(key, c[c.fuel_id == fid]), reps)           # same fuel loop as cells_t2_flash.main
        comp = [(reg.ik[int(i)], float(w)) for i, w in zip(mf.idx, mf.w)]
        row = dict(entity=ent, unmapped_fraction=mf.unmapped_fraction)
        res = {T: pr.predict_fuel(comp, T_K=T, fuel_type="jet") for T in sorted(set(NOMINAL) | set(R3_T))}
        a = res[288.15]
        if ent in flash.index:                                              # R-API-1
            f = flash.loc[ent]
            fl = a["flash_point_C"]
            row.update(in_flash=True, measured=f.measured, chosen=bool(f.chosen), low_coverage=bool(f.low_coverage), route_status=f.status,
                       api_flash=fl["value"], route_overlay=f.pred_lfl_overlay, api_v4=fl["variants"]["v4_designated"],
                       route_designated=f.pred_designated, api_status=fl["status"],
                       d_overlay=absd(fl["value"], f.pred_lfl_overlay), d_v4=absd(fl["variants"]["v4_designated"], f.pred_designated))
        else:
            row.update(in_flash=False)
        dv = dd = db_v = db_d = 0.0                                         # R-API-2: direct route values on the mapped fuel
        for T in NOMINAL:
            nu, rho = refutas_nu_safe(prov, mf.idx, mf.w, T)
            r = res[T]
            dv = max(dv, rel(r["viscosity_kin_cSt"]["value"], nu))
            dd = max(dd, rel(r["density_kg_m3"]["value"], rho))
            b1, b2 = bk.get((ent, "viscosity_kin_cSt", T)), bk.get((ent, "density_kg_m3", T))
            if b1 is not None:
                db_v = max(db_v, rel(r["viscosity_kin_cSt"]["value"], b1))
            if b2 is not None:
                db_d = max(db_d, rel(r["density_kg_m3"]["value"], b2))
        row.update(max_rel_visc=dv, max_rel_dens=dd, max_rel_visc_vs_csv=db_v, max_rel_dens_vs_csv=db_d)
        ok = True                                                           # R-API-3
        for p in ("flash_point_C", "heat_net_MJ_kg", "aromatics_wt_pct", "aromatics_vol_pct"):
            v0 = res[R3_T[0]][p]["value"]
            ok &= all(same(res[T][p]["value"], v0) for T in R3_T[1:])
        row["r3_identical"] = bool(ok)
        rows.append(row)
        print(ent, flush=True)
    d = pd.DataFrame(rows)
    d.to_csv(OUT / "api_vs_route_rows.csv", index=False, float_format="%.10g")

    fl = d[d.in_flash]
    cells = pd.read_csv(BASE / "cells_t2_flash.csv")
    api_cell = cells[(cells.set == "pooled") & (cells.model == "api_as_shipped") & (cells.variant == "designated")].iloc[0]
    head = fl[fl.chosen & ~fl.low_coverage & (fl.route_status == "ok")]
    e = head.api_flash.astype(float) - head.measured
    out = []
    m1 = max(fl.d_overlay.max(), fl.d_v4.max())
    out.append(dict(reading="R-API-1", n_rows=len(fl), max_abs_diff=m1, **{"pass": bool(m1 <= TOL_FLASH)}, before_bias=float(api_cell.bias),
                    after_bias=float(e.mean()), after_mae=float(e.abs().mean())))
    m2 = max(d.max_rel_visc.max(), d.max_rel_dens.max())
    out.append(dict(reading="R-API-2", n_rows=len(d) * len(NOMINAL), max_abs_diff=m2, **{"pass": bool(m2 <= TOL_REL)},
                    before_bias=np.nan, after_bias=np.nan, after_mae=np.nan))
    out.append(dict(reading="R-API-3", n_rows=len(d), max_abs_diff=float((~d.r3_identical).sum()), **{"pass": bool(d.r3_identical.all())},
                    before_bias=np.nan, after_bias=np.nan, after_mae=np.nan))
    out.append(dict(reading="R-API-4", n_rows=0, max_abs_diff=np.nan, **{"pass": "pending_task_6"},
                    before_bias=np.nan, after_bias=np.nan, after_mae=np.nan))
    dd = pd.DataFrame(out)
    dd.to_csv(OUT / "delta.csv", index=False, float_format="%.10g")
    print(dd.to_string())
    print("secondary: max rel vs rows_t2_bulk.csv (8-digit stored): visc %.3g dens %.3g" % (d.max_rel_visc_vs_csv.max(), d.max_rel_dens_vs_csv.max()))
    print("n fuels", len(d), "flash rows", len(fl), "headline", len(head))
    for r in out:
        print(r["reading"], "PASS" if r["pass"] is True else ("FAIL" if r["pass"] is False else r["pass"]))
    demo()


def demo() -> None:
    """60-component blend of Task 2 (tests/test_predict_route.py::test_sixty_component_blend_runs_physics): registry hydrocarbons
    C8-C16 with grid or flash data, seed 0, equal mass weights; every property at four temperatures."""
    ctx = pr.get_ctx()
    P = pr.predicates(ctx)
    cand = [i for i in range(len(ctx.reg.ik)) if (P["HAS_GRID"](i) or P["HAS_FLASH"](i)) and pr.subkind(ctx.reg.smiles[i])
            and 8 <= pr._mol_info(ctx.reg.smiles[i])[2] <= 16 and pr._mol_info(ctx.reg.smiles[i])[3]]
    pick = np.random.default_rng(0).choice(cand, 60, replace=False)
    comp = {ctx.reg.ik[int(i)]: 1.0 for i in pick}
    rows = []
    for T in DEMO_T:
        r = pr.predict_fuel(comp, T_K=T, fuel_type="jet")
        for p, v in r.items():
            if not isinstance(v, dict) or "coverage" not in v:
                continue
            cv = v["coverage"]
            rows.append(dict(T_K=T, property=p, value=v["value"], unit=v["unit"], method=v["method"], status=v["status"], in_domain=v["in_domain"],
                             own_mass=cv["own_mass"], substituted_mass=cv["substituted_mass"], unmapped_mass=cv["unmapped_mass"],
                             n_substitutions=len(cv["substitutions"]), uncertified=bool(v.get("uncertified", False)),
                             warnings=" | ".join(w.split(":")[0] for w in v["warnings"]), variants=",".join(v["variants"])))
    pd.DataFrame(rows).to_csv(OUT / "forward_demo_60.csv", index=False, float_format="%.10g")
    print("demo rows", len(rows))


if __name__ == "__main__":
    main()
