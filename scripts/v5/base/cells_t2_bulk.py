"""v5 BASE Task 6: T2 viscosity, density, freezing, heat of combustion, smoke point, distillation cells.

Viscosity / density: raw Refutas + ideal-volume mixing on the mapped measured composition at every measured temperature (AV-33-22 has 9),
the shipped API in both composition kinds, and the K1 cold offset fitted on the other source(s) with the scored double-count group removed.
Freezing aggregates S5.  Heat / smoke / distillation are new and use the mapped composition.  Nothing is fitted except the K1 offsets.
Outputs: results/v5/BASE/{cells_t2_bulk.csv, rows_t2_bulk.csv, viscosity_error_vs_T.csv, density_error_vs_T.csv}
"""
from __future__ import annotations

import json
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "src"))
from fuelblend.api.predict_blend import predict  # noqa: E402
from fuelblend.baselines.distillation import predict_d86_from_composition  # noqa: E402
from fuelblend.baselines.flash_unifac import _lookup  # noqa: E402
from fuelblend.baselines.heat_of_combustion_v2 import GroupAdditivityGCVModel, mixture_heat_of_combustion_v2  # noqa: E402
from fuelblend.baselines.sooting import smoke_point_from_ysi, ysi_mix, mixture_mw, ysi_table  # noqa: E402
from fuelblend.data.reconstruction_r3 import SubReps  # noqa: E402
from fuelblend.eval.cellkit import frame, make_row, no_eval  # noqa: E402
from fuelblend.eval.cold_offset import apply_k1, fit_k1  # noqa: E402
from fuelblend.eval.composition_route import RouteCtx, kind_bins, map_bins, predict_properties  # noqa: E402
from fuelblend.eval.group_stats import cell_metrics, logo_constant  # noqa: E402
from fuelblend.mixture.common import CompoundRegistry, PurePropertyProvider  # noqa: E402
from fuelblend.uq.conformal import coverage  # noqa: E402

V5 = ROOT / "data" / "real_fuels" / "v5"
R = ROOT / "results"
OUT = R / "v5" / "BASE"
TIER = "T2_measured"
PRIORITY = {"av3322": 0, "fuellib": 1, "ajftd_gcxgc": 2}
SRC = {"av3322": "av3322", "fuellib": "fuellib", "ajftd": "ajftd_gcxgc"}
NOMINAL = np.array([233.15, 253.15, 273.15, 288.15, 293.15, 313.15, 333.15, 353.15, 373.15])
SNAP_TOL = 0.3
S1_MARGIN = {"T_lt_263K": 0.4153576080066075, "T_263_283K": 0.6084381907110359}
TRAIN_FROM = {"av3322": "ajftd_gcxgc", "ajftd_gcxgc": "av3322", "fuellib": "av3322"}
NEW_PROPS = ["heat_of_combustion_net_MJ_kg", "smoke_point_mm", "distillation_T10_C", "distillation_T50_C", "distillation_T90_C"]
COV_MIN = 0.90


def snap(T: float) -> float:
    d = np.abs(NOMINAL - T)
    return float(NOMINAL[d.argmin()]) if d.min() <= SNAP_TOL else round(float(T), 2)


def tbin(T: float) -> str:
    return "T_lt_263K" if T < 263.15 else "T_263_283K" if T < 283.15 else "T_ge_283K"


def one_per_group(df: pd.DataFrame) -> pd.DataFrame:
    d = df.assign(_p=df.source.map(PRIORITY), _c=~df.chosen)
    return d.sort_values(["group", "_c", "_p", "entity"]).groupby("group", sort=False).head(1).drop(columns=["_p", "_c"])


# ----------------------------------------------------------------------------------------------- pure-component data for new properties
class PureData:
    def __init__(self, reg):
        self.reg = reg
        g = pd.read_csv(ROOT / "data/processed/v2/gcv_v2_MJ_kg.csv").set_index("inchikey")
        tier_a = g[g.tier == "A"].reset_index()
        mw = {k: Descriptors.MolWt(Chem.MolFromSmiles(s)) for k, s in zip(tier_a.inchikey, tier_a.smiles)}
        self.model = GroupAdditivityGCVModel()
        self.model.fit_from_tier_a(tier_a.assign(mw=tier_a.inchikey.map(mw)))
        self.gtab = g
        self.comp = pd.read_csv(ROOT / "data/unified/compounds.csv").set_index("inchikey")
        self._ysi = ysi_table()
        self.cache = {}

    def comp_info(self, i: int) -> dict:
        if i not in self.cache:
            ik, smi = self.reg.ik[i], self.reg.smiles[i]
            m = Chem.MolFromSmiles(smi)
            mw = Descriptors.MolWt(m)
            gcv = float(self.gtab.loc[ik, "value"]) if ik in self.gtab.index else self.model.predict_gcv(smi, mw)
            nc = float(sum(a.GetSymbol() == "C" for a in m.GetAtoms()))
            A, B, C = (self.comp.loc[ik, c] if ik in self.comp.index else np.nan for c in ("antoine_A", "antoine_B", "antoine_C"))
            tb = B / (A - np.log10(101.325)) - C if np.all(np.isfinite([A, B, C])) else np.nan   # Antoine inverted at 101.325 kPa (kPa/K form)
            if not (250.0 <= tb <= 800.0):
                tb = np.nan
            can = Chem.MolToSmiles(m)
            self.cache[i] = dict(smiles=smi, mw=mw, gcv=gcv, nc=nc, tb_C=tb - 273.15 if np.isfinite(tb) else np.nan, ysi_ok=can in self._ysi)
        return self.cache[i]


def heat_preds(pd_: PureData, mf) -> dict:
    ci = [pd_.comp_info(int(i)) for i in mf.idx]
    w = mf.w
    kw = dict(w_mass=w, comp_gcv_MJ_kg=[c["gcv"] for c in ci], comp_mw_g_mol=[c["mw"] for c in ci], comp_n_carbon=[c["nc"] for c in ci])
    return dict(ga_comp_H=mixture_heat_of_combustion_v2(use_standard_jet_h=False, **kw)["ncv_MJ_kg"],
                ga_const_H=mixture_heat_of_combustion_v2(use_standard_jet_h=True, **kw)["ncv_MJ_kg"])


def smoke_pred(pd_: PureData, mf):
    ci = [pd_.comp_info(int(i)) for i in mf.idx]
    ok = np.array([c["ysi_ok"] for c in ci])
    cov = float(mf.w[ok].sum())
    if cov < COV_MIN:
        return np.nan, cov
    w = mf.w[ok] / mf.w[ok].sum()
    smi = [c["smiles"] for c, o in zip(ci, ok) if o]
    mw = np.array([c["mw"] for c, o in zip(ci, ok) if o])
    x = (w / mw) / (w / mw).sum()
    return float(smoke_point_from_ysi(ysi_mix(smi, x), mixture_mw(smi, x))), cov


def d86_preds(pd_: PureData, prov, mf) -> tuple[dict, float]:
    ci = [pd_.comp_info(int(i)) for i in mf.idx]
    ok = np.array([np.isfinite(c["tb_C"]) for c in ci])
    cov = float(mf.w[ok].sum())
    if cov < COV_MIN:
        return {}, cov
    idx = mf.idx[ok]
    _, rho = prov.lookup(idx[None, :], 288.15)
    vol = (mf.w[ok] / rho[0])
    vol = vol / vol.sum()
    d = predict_d86_from_composition([c["tb_C"] for c, o in zip(ci, ok) if o], list(vol))
    return {"distillation_T10_C": d["T10_D86_C"], "distillation_T50_C": d["T50_D86_C"], "distillation_T90_C": d["T90_D86_C"]}, cov


# ----------------------------------------------------------------------------------------------- long rows
def build_rows() -> pd.DataFrame:
    warnings.simplefilter("ignore")
    reg = CompoundRegistry(ROOT / "data/unified/compounds.csv")
    prov = PurePropertyProvider(reg, ROOT / "data/cache/pure_grid.npz", verbose=False)
    reps = SubReps(reg, prov, vp_lookup=lambda s: _lookup(s, None))
    ctx = RouteCtx(reg, prov)
    pdat = PureData(reg)
    t = pd.read_csv(OUT / "t2_targets.csv")
    fu = pd.read_csv(OUT / "t2_fuels.csv").set_index("entity")
    comps = {"av3322": pd.read_csv(V5 / "crc_av3322_composition.csv"), "fuellib": pd.read_csv(V5 / "fuellib_gcxgc_compositions.csv"),
             "ajftd": pd.read_csv(V5 / "ajftd_gcxgc_composition.csv")}
    t = t[t.property.isin(["viscosity_kin_cSt", "density_kg_m3"] + NEW_PROPS)].copy()
    rows = []
    api_cache: dict = {}

    def api(ent, mf, T, kind):
        k = (ent, round(T, 2), kind)
        if k not in api_cache:
            try:
                a = predict([(reg.ik[int(i)], float(w)) for i, w in zip(mf.idx, mf.w)], T_K=T, fuel_type="jet", composition_kind=kind)
                v, d = a["viscosity_dyn_mPas"], a["density_kg_m3"]
                rho = d["value"]
                f = lambda x: x / rho * 1000.0 if (x is not None and rho and np.isfinite(x)) else np.nan             # noqa: E731
                api_cache[k] = dict(pred=f(v["value"]), lo=f(v["lo"]), hi=f(v["hi"]), rho=rho, rho_lo=d["lo"], rho_hi=d["hi"],
                                    method=v["method"], dmethod=d["method"], in_domain=bool(v["in_domain"]))
            except Exception as e:                                                                      # noqa: BLE001
                api_cache[k] = dict(pred=np.nan, lo=np.nan, hi=np.nan, rho=np.nan, rho_lo=np.nan, rho_hi=np.nan, method="error:" + str(e)[:40],
                                    dmethod="error", in_domain=False)
        return api_cache[k]

    for ent, g in t.groupby("entity"):
        src = SRC[fu.loc[ent, "source"]]
        key = fu.loc[ent, "source"]
        fid = int(g.fuel_id.iloc[0]) if key == "ajftd" else g.fuel_id.iloc[0]
        c = comps[key]
        mf = map_bins(kind_bins(key, c[c.fuel_id == fid]), reps)
        meta = dict(entity=ent, source=src, group=fu.loc[ent, "group"], chosen=bool(fu.loc[ent, "chosen_for_group"]),
                    low_coverage=bool(fu.loc[ent, "low_coverage"]), label=g.label.iloc[0])
        extra_cache: dict = {}
        for r in g.itertuples():
            base = dict(meta, property=r.property, T_raw=r.T_K, method_family=r.method_family, measured=r.measured)
            if r.property in ("viscosity_kin_cSt", "density_kg_m3"):
                o = predict_properties(ctx, mf, pd.DataFrame([dict(property=r.property, T_K=r.T_K)])).iloc[0]
                Tb = snap(r.T_K)
                rows.append(dict(base, T_K=Tb, model="refutas_raw" if "visc" in r.property else "ideal_volume", variant="designated",
                                 pred=o.pred, status=o.status))
                for kind, var in (("known", "known"), ("reconstructed", "reconstructed_route")):
                    a = api(ent, mf, r.T_K, kind)
                    if "visc" in r.property:
                        rows.append(dict(base, T_K=Tb, model="api_as_shipped", variant=var, pred=a["pred"], status="ok" if np.isfinite(a["pred"]) else "non_finite",
                                         lo=a["lo"], hi=a["hi"], api_method=a["method"], api_in_domain=a["in_domain"]))
                    elif kind == "known":
                        rows.append(dict(base, T_K=Tb, model="api_as_shipped", variant=var, pred=a["rho"], status="ok" if a["rho"] and np.isfinite(a["rho"]) else "non_finite",
                                         lo=a["rho_lo"], hi=a["rho_hi"], api_method=a["dmethod"], api_in_domain=a["in_domain"]))
            else:
                if "heat" not in extra_cache:
                    extra_cache["heat"] = heat_preds(pdat, mf)
                    extra_cache["smoke"] = smoke_pred(pdat, mf)
                    extra_cache["d86"] = d86_preds(pdat, prov, mf)
                if r.property == "heat_of_combustion_net_MJ_kg":
                    h = extra_cache["heat"]
                    for m in ("ga_comp_H", "ga_const_H"):
                        rows.append(dict(base, T_K=np.nan, model=m, variant="designated" if m == "ga_comp_H" else "const_H", pred=h[m],
                                         status="ok" if np.isfinite(h[m]) else "non_finite"))
                elif r.property == "smoke_point_mm":
                    p, cov = extra_cache["smoke"]
                    rows.append(dict(base, T_K=np.nan, model="ysi_route", variant="designated", pred=p, status="ok" if np.isfinite(p) else "ysi_coverage_lt_0.90", cover=cov))
                else:
                    dd, cov = extra_cache["d86"]
                    p = dd.get(r.property, np.nan)
                    rows.append(dict(base, T_K=np.nan, model="antoine_tb_d86", variant="designated", pred=p, status="ok" if np.isfinite(p) else "tb_coverage_lt_0.90", cover=cov))
        extra_cache.clear()
    d = pd.DataFrame(rows)
    d["n_dup"] = d.groupby(["entity", "property", "T_K", "model", "variant"]).pred.transform("size")
    return d


def k1_rows(d: pd.DataFrame) -> pd.DataFrame:
    """K1 per-temperature offset fitted on the other source's one-per-group rows, with the scored double-count group removed."""
    v = d[(d.property == "viscosity_kin_cSt") & (d.variant == "designated") & (d.status == "ok") & ~d.low_coverage].copy()
    v["resid"] = np.log(v.measured) - np.log(v.pred)
    out = []
    for grp, s in v.groupby("group"):
        for src in s.source.unique():
            tr = pd.concat([one_per_group(g) for _, g in v[(v.source == TRAIN_FROM[src]) & (v.group != grp)].groupby("T_K")])  # one per group at each T
            k1 = fit_k1(tr.T_K, tr.resid)
            sc = s[s.source == src]
            off = apply_k1(k1, sc.T_K.to_numpy())
            out.append(sc.assign(variant="k1_cold_offset_loso", pred=np.exp(np.log(sc.pred) + off), status=np.where(np.isfinite(off), "ok", "no_offset_at_T")))
    k = pd.concat(out, ignore_index=True)
    k["model"] = "refutas_raw+k1"
    return k.drop(columns=["resid"])


# ----------------------------------------------------------------------------------------------- cells
def unsafe_visc(T):
    if T >= 283.15:
        return None
    m = S1_MARGIN["T_lt_263K" if T < 263.15 else "T_263_283K"]
    return lambda e: e < -m


def cells(d: pd.DataFrame) -> list[dict]:
    out = []
    ev = "results/v5/BASE/rows_t2_bulk.csv"
    srcs = {"pooled": None, "av3322": "av3322", "ajftd_gcxgc": "ajftd_gcxgc", "fuellib": "fuellib"}
    for (prop, model, var), g in d[d.property.isin(["viscosity_kin_cSt", "density_kg_m3"])].groupby(["property", "model", "variant"]):
        visc = prop == "viscosity_kin_cSt"
        for sname, sv in srcs.items():
            gs = g if sv is None else g[g.source == sv]
            keys = [(T, f"{sname}") for T in sorted(gs.T_K.unique())] + [(b, f"{sname}_{b}") for b in ("T_lt_263K", "T_263_283K", "T_ge_283K")]
            for tk, setname in keys:
                cell = gs[gs.T_K == tk] if isinstance(tk, float) else gs[gs.T_K.map(tbin) == tk]
                if cell.empty:
                    continue
                nlow = int(cell.low_coverage.sum())
                hd = one_per_group(cell[~cell.low_coverage])
                nonf = int((hd.status != "ok").sum())
                hd = hd[hd.status == "ok"]
                if hd.empty:
                    continue
                if visc:
                    e = np.log(hd.pred) - np.log(hd.measured)
                    meas = np.log(hd.measured)
                else:
                    e = hd.pred - hd.measured
                    meas = hd.measured
                # constant: leave-one-group-out per temperature (per bin when pooling temperatures)
                strata = hd.T_K.to_numpy()
                et = logo_constant(meas, hd.group, strata) - meas.to_numpy()
                T0 = tk if isinstance(tk, float) else (253.15 if tk == "T_lt_263K" else 273.15 if tk == "T_263_283K" else 300.0)
                unsafe = unsafe_visc(T0) if visc else None
                m90, c90, m90s1, c90s1, extra = np.nan, np.nan, np.nan, np.nan, {}
                if visc and model == "refutas_raw":
                    bn = tbin(T0)
                    m90 = np.inf if bn != "T_ge_283K" else 0.1855823            # configs/uncertainty_v4.yaml lab raw_refutas: inf below 283 K
                    m90s1 = S1_MARGIN.get(bn, 0.1855823)
                    if np.isfinite(m90):
                        k, n = coverage(e, hd.group, m90, two_sided=True); c90 = k / n
                    k, n = coverage(e, hd.group, m90s1, two_sided=True); c90s1 = k / n
                    extra = dict(margin90_s1=m90s1, coverage90_s1=c90s1)
                elif visc and model == "api_as_shipped" and hd.lo.notna().any():
                    inside = ((hd.measured >= hd.lo) & (hd.measured <= hd.hi))[hd.lo.notna()]
                    c90 = float(inside.groupby(hd.group[hd.lo.notna()]).all().mean())
                    extra = dict(note="API methods: " + json.dumps(hd.api_method.value_counts().to_dict()))
                elif (not visc) and model == "api_as_shipped" and hd.lo.notna().any():
                    inside = ((hd.measured >= hd.lo) & (hd.measured <= hd.hi))[hd.lo.notna()]
                    c90 = float(inside.groupby(hd.group[hd.lo.notna()]).all().mean())
                    extra = dict(note="API methods: " + json.dumps(hd.api_method.value_counts().to_dict()))
                row = make_row(prop, TIER, setname, model, var, "recomputed_v5", "descriptive", ev, err=e, groups=hd.group, unsafe=unsafe,
                               err_trivial=et, trivial_model="logo_constant_per_T", margin90=m90, coverage90=c90,
                               n_excluded_nonfinite=nonf, n_low_coverage=nlow, extra=dict(T_K=tk if isinstance(tk, float) else np.nan, **extra))
                out.append(row)
    # heat, smoke, distillation
    for (prop, model, var), g in d[d.property.isin(NEW_PROPS)].groupby(["property", "model", "variant"]):
        label = "exploratory" if prop != "heat_of_combustion_net_MJ_kg" else "descriptive"
        if prop == "heat_of_combustion_net_MJ_kg":
            sets = {"measured_pooled": g[g.method_family == "measured"], "calculated_pooled": g[g.method_family == "calculated"],
                    "av3322": g[g.source == "av3322"], "ajftd_gcxgc_measured": g[(g.source == "ajftd_gcxgc") & (g.method_family == "measured")]}
        else:
            sets = {"pooled": g, "av3322": g[g.source == "av3322"], "fuellib": g[g.source == "fuellib"], "ajftd_gcxgc": g[g.source == "ajftd_gcxgc"]}
        for sname, s in sets.items():
            if s.empty:
                continue
            nlow = int(s.low_coverage.sum())
            hd = one_per_group(s[~s.low_coverage])
            nonf = int((hd.status != "ok").sum())
            excl = hd[hd.status != "ok"].status.value_counts().to_dict()
            hd = hd[hd.status == "ok"]
            if hd.group.nunique() < 2:  # nothing scorable: say so explicitly (Review Focus 4), never drop the cell
                r = no_eval(prop, TIER, sname, f"{hd.group.nunique()} scorable group(s); excluded by status: {excl}; low_coverage (mapping) {nlow}",
                            model=model, variant=var, evidence=ev, label=label)
                r["n_excluded_nonfinite"], r["n_low_coverage"] = nonf, nlow
                out.append(r)
                continue
            e = hd.pred - hd.measured
            et = logo_constant(hd.measured, hd.group, hd.method_family.fillna("-")) - hd.measured.to_numpy()
            out.append(make_row(prop, TIER, sname, model, var, "recomputed_v5", label, ev, err=e, groups=hd.group, err_trivial=et,
                                trivial_model="logo_constant_within_method" if prop.startswith("heat") else "logo_constant",
                                n_excluded_nonfinite=nonf, n_low_coverage=nlow))
    return out


def freezing_cells() -> list[dict]:
    s5 = pd.read_csv(R / "v5/S5/s5_predictions.csv")
    out = []
    ev = "results/v5/S5/s5_predictions.csv"
    names = {"I": ("ideal_svl", "designated"), "F1": ("coutinho_f1", "f1_recalibrated"), "P": ("ridge_p", "p_ridge")}
    for role, label in (("headline", "descriptive"), ("locked", "already_read")):
        r = s5[s5.role == role]
        for sname, sv in (("pooled", None), ("AV-33-22", "AV-33-22"), ("AJFTD", "AJFTD"), ("FuelLib", "FuelLib")):
            rr = r if sv is None else r[r.source == sv]
            for col, (model, var) in names.items():
                m = rr.dropna(subset=[col, "C", "freezing_C"])
                if m.empty:
                    continue
                out.append(make_row("freezing_point_C", TIER, sname + ("" if role == "headline" else "_locked"), model, var, "aggregated_v5", label, ev,
                                    err=m[col] - m.freezing_C, groups=m.group, err_trivial=m.C - m.freezing_C, trivial_model="logo_constant",
                                    n_excluded_nonfinite=int(len(rr) - len(m))))
    return out


def error_vs_T(d: pd.DataFrame) -> None:
    for prop, fn in (("viscosity_kin_cSt", "viscosity_error_vs_T.csv"), ("density_kg_m3", "density_error_vs_T.csv")):
        x = d[(d.property == prop) & (d.variant == "designated") & (d.status == "ok") & ~d.low_coverage]
        rows = []
        for sv, gs in list(x.groupby("source")) + [("pooled", x)]:
            for T, g in gs.groupby("T_K"):
                h = one_per_group(g)
                e = np.log(h.pred / h.measured) if prop.startswith("visc") else h.pred - h.measured
                rows.append(dict(source=sv, T_K=T, n_groups=h.group.nunique(), bias_ln=e.mean() if prop.startswith("visc") else np.nan,
                                 bias=e.mean(), mae_ln=e.abs().mean(), sd_ln=e.std(ddof=1) if len(e) > 1 else np.nan))
        pd.DataFrame(rows).to_csv(OUT / fn, index=False, float_format="%.6g")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    if "--from-rows" in sys.argv:  # reuse saved predictions (build_rows is slow); K1 rows are refitted
        d = pd.read_csv(OUT / "rows_t2_bulk.csv")
        d = d[d.variant != "k1_cold_offset_loso"].reset_index(drop=True)
        d["n_dup"] = d.groupby(["entity", "property", "T_K", "model", "variant"]).pred.transform("size")
    else:
        d = build_rows()
    d = pd.concat([d, k1_rows(d)], ignore_index=True)
    keep = ["entity", "source", "group", "chosen", "low_coverage", "property", "T_K", "T_raw", "method_family", "measured", "model", "variant", "pred", "status",
            "lo", "hi", "api_method", "api_in_domain", "cover"]
    d[[k for k in keep if k in d.columns]].to_csv(OUT / "rows_t2_bulk.csv", index=False, float_format="%.8g")
    error_vs_T(d)
    c = cells(d) + freezing_cells()
    c += [no_eval("viscosity_kin_cSt", TIER, "pooled", f"{m}: trained inside build_viscosity_leaderboard_v4.py::main on CRC 647 fuels; no saved model to apply to T2 (Part C Task 16f)",
                  model=m, evidence="scripts/v4/build_viscosity_leaderboard_v4.py") for m in ("structured_ridge", "bilinear_pairwise", "structured_gp")]
    c += [no_eval("pour_point_C", TIER, "-", "no pour-point measurements for T2 fuels; PPBI needs pure-component pour points the registry lacks")]
    c += [no_eval("dcn", TIER, "-", "no cetane model scored (v4: uncertified)")]
    df = frame(c)
    df.to_csv(OUT / "cells_t2_bulk.csv", index=False, float_format="%.8g")
    print(df[df.status != "no_evaluation"].groupby(["property", "model", "variant"]).size().to_string())


if __name__ == "__main__":
    main()
