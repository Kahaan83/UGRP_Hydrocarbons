"""v5 S3: composition error versus physics error.  Readings R-S3-1..5 are declared in results/v5/S3/README.md (before scoring).

Every fuel with a measured GCxGC composition is scored twice with the SAME property call (`predict_comp`):
  M  measured composition, mapped to representative compounds (G1 map_fuel for AJFTD; S2 maps for AV-33-22 and FuelLib)
  R  composition reconstructed by the r2 route (reconstruction_r2.reconstruct) from bulk descriptors of the same fuel only:
     D86 (IBP, T10, T20, T50, T90, FBP), aromatics vol%, naphthalenes vol% (diaromatics), olefins vol% (absent -> 0).
Table G.1 hydrocarbon types exist only for CRC 647 fuels, so none of these fuels can use r3; all R compositions are r2.
Nothing is fitted.  Outputs: results/v5/S3/*.csv.  The README is not touched by this script.
"""
from __future__ import annotations

import importlib.util
import json
import re
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors

RDLogger.DisableLog("rdApp.*")
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from fuelblend.baselines.flash_unifac import UnifacFallbackWarning, _lookup, flash_point_unifac_detail  # noqa: E402
from fuelblend.data.reconstruction_r2 import RealisticReps, nominal_sat_split_detail, reconstruct  # noqa: E402
from fuelblend.data.reconstruction_r3 import SubReps  # noqa: E402
from fuelblend.mixture.common import CompoundRegistry, PurePropertyProvider  # noqa: E402
from fuelblend.thermo.coutinho_wax import get_pure_thermal_properties  # noqa: E402

V5 = ROOT / "data" / "real_fuels" / "v5"
OUT = ROOT / "results" / "v5" / "S3"
S2 = ROOT / "results" / "v5" / "S2"
JETCLASS = ROOT / "data" / "real_fuels" / "v4" / "jet_class_compositions.csv"
SPLITS = ROOT / "results" / "v5" / "S4" / "splits_ajftd_v5.json"

MAX_UNMAPPED = 0.10
N_BOOT, SEED, MIN_GROUPS = 5000, 0, 8
SRC_PRIORITY = {"AV-33-22": 0, "FuelLib": 1, "AJFTD": 2}       # tie-break when two group members have equally many scored rows
# targets: label -> (property, T_K or None, tolerance K)
TARGETS = {"density_15C": ("density_kg_m3", 288.65, 0.6), "visc_-40C": ("viscosity_kin_cSt", 233.15, 0.3),
           "visc_-20C": ("viscosity_kin_cSt", 253.15, 0.3), "visc_40C": ("viscosity_kin_cSt", 313.15, 0.3),
           "flash": ("flash_point_C", None, 0), "freezing": ("freezing_point_C", None, 0)}
# fixed tolerances for the near-identical-bulk double-count check (set before any error was looked at)
DUP_TOL = {"dens15": 1.5, "aromatics": 0.5, "t10": 2.0, "t90": 2.0, "flash": 2.0, "freezing": 2.0}
DUP_MIN_FEATURES = 3

# ---- R route: whitelisted descriptors only ------------------------------------------------------------------------
R_KEYS = ("d86_ibp", "d86_t10", "d86_t20", "d86_t50", "d86_t90", "d86_fbp", "aromatics_vol_pct", "olefins_vol_pct",
          "diaromatics_vol_pct")
D86_PROP = {"distillation_T0_C": "d86_ibp", "distillation_T10_C": "d86_t10", "distillation_T20_C": "d86_t20",
            "distillation_T50_C": "d86_t50", "distillation_T90_C": "d86_t90", "distillation_T100_C": "d86_fbp"}
R_PROPS = set(D86_PROP) | {"aromatics_vol_pct", "naphthalenes_vol_pct", "olefins_vol_pct"}


def r_descriptors(bulk_long: pd.DataFrame) -> dict:
    """{R_KEYS: value} from long rows (property, value, method) of ONE fuel.  Only exact whitelisted property names are used
    (D86 points, aromatics, naphthalenes, olefins); anything else in the rows is ignored."""
    b = bulk_long[bulk_long.property.isin(R_PROPS) & bulk_long.value.notna()]
    out = {}
    for prop, key in D86_PROP.items():
        v = b.loc[b.property == prop, "value"]
        if len(v):
            out[key] = float(v.median())
    a = b[b.property == "aromatics_vol_pct"]
    if len(a):
        pref = a[a.method.astype(str).str.contains("1319")]
        out["aromatics_vol_pct"] = float((pref if len(pref) else a).value.median())
    for prop, key in (("naphthalenes_vol_pct", "diaromatics_vol_pct"), ("olefins_vol_pct", "olefins_vol_pct")):
        v = b.loc[b.property == prop, "value"]
        if len(v):
            out[key] = float(v.median())
    assert set(out) <= set(R_KEYS)
    return out


def prior_split(exclude_posf: set[str] | None = None) -> tuple[float, float, float]:
    """r2 fixed saturate split (mean n/iso/cyclo share of the GC jet fuels); optionally leave out named POSF fuels."""
    (n, i, c), used, _ = nominal_sat_split_detail(JETCLASS)
    if not exclude_posf:
        return n, i, c
    d = pd.read_csv(JETCLASS)
    sh = []
    for fid in used:
        if fid in exclude_posf:
            continue
        cf = json.loads(d.loc[d.fuel_id == fid, "class_fractions"].iloc[0])
        sat = np.array([cf.get(k, np.nan) for k in ("n_paraffin", "iso_paraffin", "cycloalkane")])
        sh.append(sat / sat.sum())
    m = np.mean(sh, axis=0)
    return float(m[0]), float(m[1]), float(m[2])


# ---- one property call for both M and R -----------------------------------------------------------------------------
class Ctx:
    def __init__(self, reg, prov, refutas_nu, ideal_sle_C, tmap=None):
        self.reg, self.prov, self.refutas_nu, self.ideal_sle_C, self.tmap = reg, prov, refutas_nu, ideal_sle_C, tmap or {}
        self._mw: dict[str, float] = {}

    def mw(self, smi: str) -> float:
        if smi not in self._mw:
            self._mw[smi] = Descriptors.MolWt(Chem.MolFromSmiles(smi))
        return self._mw[smi]


def predict_comp(ctx: Ctx, idx, w, targets: pd.DataFrame) -> pd.DataFrame:
    """Existing property models on a composition (registry indices, mass fractions).  Returns targets + pred, status.
    Flash: ideal Le Chatelier (as G1/S2); viscosity: Refutas; density: ideal volume mixing; freezing: ideal SvL."""
    idx = np.asarray(idx, np.int64)
    w = np.asarray(w, float)
    w = w / w.sum()
    iks = [ctx.reg.ik[i] for i in idx]
    smi = [ctx.reg.smiles[i] for i in idx]
    x = w / np.array([ctx.mw(s) for s in smi])
    x = x / x.sum()
    rows = []
    for r in targets.itertuples():
        pred, status = np.nan, "ok"
        if r.property == "flash_point_C":
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", UnifacFallbackWarning)
                d = flash_point_unifac_detail(smi, x, use_unifac=False)
            pred, status = d["flash_C"], d["status"]
        elif r.property in ("viscosity_kin_cSt", "density_kg_m3"):
            nu, rho, _, _ = ctx.refutas_nu(ctx.prov, idx, w, float(r.T_K))
            pred = nu if r.property == "viscosity_kin_cSt" else rho
        else:
            pred = ctx.ideal_sle_C(x, iks, ctx.tmap)
        rows.append(dict(pred=float(pred) if pred is not None else np.nan, status=status))
    out = targets.reset_index(drop=True).copy()
    out["pred"] = [r["pred"] for r in rows]
    out["status"] = [r["status"] for r in rows]
    return out


def err_of(prop: str, pred, meas):
    pred, meas = np.asarray(pred, float), np.asarray(meas, float)
    if prop == "viscosity_kin_cSt":
        with np.errstate(all="ignore"):
            return np.log(pred / meas)
    return pred - meas


def score_entity(ctx: Ctx, m_comp, r_comp, targets: pd.DataFrame, r_comp_lo=None) -> pd.DataFrame:
    """Score ONE fuel's targets with M and R through the same predict_comp.  Compositions are (idx, w) or None."""
    t = targets.reset_index(drop=True).copy()
    for tag, comp in (("M", m_comp), ("R", r_comp), ("Rlo", r_comp_lo if r_comp_lo is not None else r_comp)):
        if comp is None:
            t[f"pred_{tag}"], t[f"status_{tag}"] = np.nan, "no_composition"
        else:
            p = predict_comp(ctx, comp[0], comp[1], t)
            t[f"pred_{tag}"], t[f"status_{tag}"] = p.pred.to_numpy(), p.status.to_numpy()
        e = np.full(len(t), np.nan)
        for prop in t.property.unique():
            m = (t.property == prop).to_numpy()
            e[m] = err_of(prop, t.loc[m, f"pred_{tag}"], t.loc[m, "measured"])
        e[(t[f"status_{tag}"] != "ok").to_numpy()] = np.nan
        t[f"err_{tag}"] = e
    return t


# ---- groups ------------------------------------------------------------------------------------------------------------
class UF:
    def __init__(self, items):
        self.p = {i: i for i in items}

    def find(self, a):
        while self.p[a] != a:
            self.p[a] = self.p[self.p[a]]
            a = self.p[a]
        return a

    def union(self, a, b):
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            self.p[max(ra, rb)] = min(ra, rb)


def select_one_per_group(ent: pd.DataFrame) -> pd.DataFrame:
    """ent: columns entity, group, source, n_pair_rows, pairable, group_locked.  One entity per group: pairable members only,
    most scored (M and R finite) rows, ties by source priority then entity id; groups containing a locked fuel are dropped."""
    e = ent[ent.pairable & ~ent.group_locked].copy()
    e["prio"] = e.source.map(SRC_PRIORITY)
    e = e.sort_values(["group", "n_pair_rows", "prio", "entity"], ascending=[True, False, True, True])
    return e.drop_duplicates("group", keep="first").drop(columns="prio")


# ---- statistics --------------------------------------------------------------------------------------------------------
def paired_bootstrap(gid, eM, eR, n=N_BOOT, seed=SEED):
    """R-minus-M MAE difference, rows resampled by group (same indices for M and R).  -> (diff, lo, hi)."""
    gid, eM, eR = np.asarray(gid), np.abs(np.asarray(eM, float)), np.abs(np.asarray(eR, float))
    ug, inv = np.unique(gid, return_inverse=True)
    G = len(ug)
    sM, sR, cnt = (np.bincount(inv, weights=v, minlength=G) for v in (eM, eR, np.ones(len(eM))))
    rng = np.random.default_rng(seed)
    ix = rng.integers(0, G, size=(n, G))
    d = sR[ix].sum(1) / cnt[ix].sum(1) - sM[ix].sum(1) / cnt[ix].sum(1)
    return float(eR.mean() - eM.mean()), float(np.percentile(d, 2.5)), float(np.percentile(d, 97.5))


def verdict_s31(n_groups, diff, lo, hi, biasM, biasR):
    if n_groups < MIN_GROUPS:
        return "n too small"
    if lo > 0:
        return "composition error matters"
    if hi < 0:
        return "R better than M (measured composition worse)"
    return "physics-limited" if abs(biasM) >= abs(biasR) else "inconclusive (CI includes 0, M bias smaller)"


def central_rows(d: pd.DataFrame, setname: str, which: str = "R") -> list[dict]:
    """d: long predictions restricted to one set; paired rows only (M and `which` finite)."""
    rows = []
    for tl, (prop, _, _) in TARGETS.items():
        s = d[(d.target == tl) & np.isfinite(d.err_M) & np.isfinite(d[f"err_{which}"])]
        if s.empty:
            rows.append(dict(set=setname, target=tl, n_groups=0, n_rows=0, verdict="not evaluable (no paired rows)"))
            continue
        eM, eR = s.err_M.to_numpy(), s[f"err_{which}"].to_numpy()
        ng = s.group.nunique()
        diff, lo, hi = paired_bootstrap(s.group, eM, eR)
        rec = dict(set=setname, target=tl, unit="ln(pred/meas)" if prop == "viscosity_kin_cSt" else ("kg/m3" if prop.startswith("dens") else "C"),
                   n_groups=ng, n_rows=len(s), bias_M=eM.mean(), bias_R=eR.mean(), mae_M=np.abs(eM).mean(), mae_R=np.abs(eR).mean(),
                   mae_diff_R_minus_M=diff, ci_lo=lo, ci_hi=hi, verdict=verdict_s31(ng, diff, lo, hi, eM.mean(), eR.mean()))
        if prop == "viscosity_kin_cSt":
            rec.update(bias_M_pct=np.expm1(eM.mean()) * 100, bias_R_pct=np.expm1(eR.mean()) * 100)
        if prop == "flash_point_C":
            rec.update(frac_over5_M=float((eM > 5).mean()), frac_over5_R=float((eR > 5).mean()))
        rows.append(rec)
    return rows


# ---- data assembly -----------------------------------------------------------------------------------------------------
def load_module(name, rel):
    spec = importlib.util.spec_from_file_location(name, ROOT / rel)
    m = importlib.util.module_from_spec(spec)
    sys.modules[name] = m
    spec.loader.exec_module(m)
    return m


def tag_targets(df: pd.DataFrame) -> pd.DataFrame:
    """df: property, T_K (nan for flash/freezing), measured -> adds `target` label, keeps only declared S3 targets."""
    out = []
    for r in df.itertuples():
        lab = None
        for tl, (prop, T, tol) in TARGETS.items():
            if r.property != prop:
                continue
            if T is None or (np.isfinite(r.T_K) and abs(r.T_K - T) <= tol):
                lab = tl
        out.append(lab)
    df = df.copy()
    df["target"] = out
    return df[df.target.notna()].reset_index(drop=True)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    g1 = load_module("g1_eval", "scripts/v5/eval_ajftd_gcxgc_v5.py")
    vlb = load_module("vlb_s3", "scripts/v4/build_viscosity_leaderboard_v4.py")
    reg = CompoundRegistry(ROOT / "data" / "unified" / "compounds.csv")
    prov = PurePropertyProvider(reg, ROOT / "data" / "cache" / "pure_grid.npz", verbose=False)
    vp = lambda s: _lookup(s, None)                                                    # noqa: E731
    reps_m = SubReps(reg, prov, vp_lookup=vp)
    reps_r = RealisticReps(reg, prov, vp_lookup=vp)
    ik_index = {k: i for i, k in enumerate(reg.ik)}
    roles = json.load(open(SPLITS))
    role, dup_map = roles["fuel_role"], roles["dup_map"]

    ents: dict[str, dict] = {}                                                         # entity id -> record

    # ---------------- AJFTD GCxGC
    mapped, summ = g1.build_mapping(reps_m)                                            # in memory only (G1 files not rewritten)
    tg_aj = g1.targets()
    long_aj = pd.read_csv(V5 / "ajftd_fuels_long.csv", low_memory=False)
    gfu = pd.read_csv(V5 / "ajftd_gcxgc_fuels.csv")
    labels = dict(zip(gfu.fuel_id, gfu.fuel_label.astype(str)))
    for s in summ.itertuples():
        fid = int(s.fuel_id)
        g = mapped[mapped.fuel_id == fid]
        t = tg_aj[tg_aj.fuel_id == fid].rename(columns={"T_key": "T_K"})
        t["T_K"] = t.T_K.where(t.T_K > 0, np.nan)
        canon = dup_map.get(str(fid), fid)
        ents[f"AJFTD:{fid}"] = dict(
            source="AJFTD", fid=fid, label=labels[fid], n_bins=int(s.n_bins), unmapped_M=float(s.unmapped_fraction),
            M=(g.reg_index.to_numpy(np.int64), g.fraction.to_numpy(float)) if len(g) else None, comp_hash=s.comp_hash,
            locked=role.get(str(canon)) == "locked", targets=tag_targets(t[["property", "T_K", "measured", "method"]]),
            bulk=long_aj[long_aj.fuel_id == fid][["property", "value", "method"]], canon=f"AJFTD:{canon}")
    # ---------------- CRC AV-33-22
    comp = pd.read_csv(V5 / "crc_av3322_composition.csv")
    amap = pd.read_csv(S2 / "av3322_representative_map.csv")
    props = pd.read_csv(V5 / "crc_av3322_properties.csv")
    bulk = pd.read_csv(V5 / "crc_av3322_bulk.csv")
    for fid, gdf in comp.groupby("fuel_id"):
        m = gdf.merge(amap[["group", "kind", "carbon_number", "rep_inchikey", "mapped"]], on=["group", "carbon_number"], how="left")
        tot = m.mass_pct.sum()
        mm = m[m.mapped.astype(bool)]
        agg = mm.groupby("rep_inchikey").mass_pct.sum()
        idx = np.array([ik_index[k] for k in agg.index], np.int64)
        assert prov.valid[idx].all(), "AV-33-22 representative without valid pure eta/rho"
        p = props[props.fuel_id == fid]
        b = bulk[bulk.fuel_id == fid].iloc[0]
        tt = [dict(property="flash_point_C", T_K=np.nan, measured=float(b.coa_flash_point_C), method="CoA")]
        for prop in ("density_kg_m3", "viscosity_kin_cSt"):
            q = p[(p.property == prop) & p.value.notna() & p.T_K.notna()]
            for T, qq in q.groupby(q.T_K.round(2)):
                tt.append(dict(property=prop, T_K=float(T), measured=float(qq.value.median()), method="+".join(sorted(set(qq.method.astype(str))))))
        q = p[(p.property == "freezing_point_C") & p.value.notna() & p.qualifier.isna()]
        if len(q):
            tt.append(dict(property="freezing_point_C", T_K=np.nan, measured=float(q.value.median()), method="+".join(sorted(set(q.method.astype(str))))))
        ents[f"AV-33-22:{fid}"] = dict(
            source="AV-33-22", fid=fid, label=fid, n_bins=len(gdf), unmapped_M=float(1 - agg.sum() / tot), M=(idx, agg.to_numpy(float)), comp_hash="",
            locked=False, targets=tag_targets(pd.DataFrame(tt)), bulk=p[["property", "value", "method"]], canon=f"AV-33-22:{fid}")
    # ---------------- FuelLib (4 fuels with a measured composition and properties)
    fl = pd.read_csv(V5 / "fuellib_gcxgc_compositions.csv")
    fmap = pd.read_csv(S2 / "fuellib_representative_map.csv")
    flp = pd.read_csv(V5 / "fuellib_properties.csv")
    for fid in sorted(f for f in flp.fuel_id.unique() if f in set(fl.fuel_id) and (flp[flp.fuel_id == f].property == "flash_point_C").any()):
        m = fl[fl.fuel_id == fid].merge(fmap[["gc_bin", "n_carbon", "rep_inchikey", "mapped"]], on=["gc_bin", "n_carbon"], how="left")
        tot = m.fraction.sum()
        agg = m[m.mapped.astype(bool)].groupby("rep_inchikey").fraction.sum()
        idx = np.array([ik_index[k] for k in agg.index], np.int64)
        assert prov.valid[idx].all()
        p = flp[flp.fuel_id == fid].copy()
        p["T_K"] = p.T_C + 273.15
        t = p[p.property.isin(["flash_point_C", "freezing_point_C", "density_kg_m3", "viscosity_kin_cSt"])].rename(columns={"value": "measured"})
        bl = p[p.property.isin(R_PROPS)][["property", "value", "method"]]
        label = fid.replace("FUELLIB_posf", "")
        ents[f"FuelLib:{fid}"] = dict(
            source="FuelLib", fid=fid, label=label, n_bins=len(m), unmapped_M=float(1 - agg.sum() / tot), M=(idx, agg.to_numpy(float)), comp_hash="",
            locked=False, targets=tag_targets(t[["property", "T_K", "measured", "method"]]), bulk=bl, canon=f"FuelLib:{fid}")

    # ---------------- R compositions (descriptor route)
    std_prior = prior_split()
    prior_used = set(nominal_sat_split_detail(JETCLASS)[1])
    prior_ids = {"POSF-" + lab for lab in {e["label"] for e in ents.values()}}
    r_rows, fuel_rows = [], []
    for eid, e in ents.items():
        row = r_descriptors(e["bulk"])
        e["desc"] = row
        n_d86 = sum(k in row for k in R_KEYS if k.startswith("d86"))
        why = ""
        if "aromatics_vol_pct" not in row:
            why = "no aromatics descriptor (non-GCxGC) in the fuel's own bulk data"
        elif n_d86 < 4:
            why = f"fewer than 4 D86 points ({n_d86})"
        e["R"] = e["R_lo"] = None
        e["fallback_R"], e["prior_variant"] = np.nan, "standard"
        if not why:
            def build(prior):
                fr, fb = reconstruct(row, prior, reps_r)
                if fr.empty:
                    return None, 1.0, fr
                return (np.array([ik_index[k] for k in fr.inchikey], np.int64), fr.fraction.to_numpy(float)), float(fb), fr
            res = build(std_prior)
            if res[0] is None:
                why = "reconstruct returned no composition"
            else:
                e["R"], e["fallback_R"], fr = res
                posf = "POSF-" + e["label"]
                if posf in prior_used:
                    e["prior_variant"] = "label in prior set; Rlo uses leave-out prior"
                    lo = build(prior_split({posf}))
                    e["R_lo"] = lo[0]
                else:
                    e["R_lo"] = e["R"]
                for r in fr.itertuples():
                    r_rows.append(dict(entity=eid, source=e["source"], label=e["label"], inchikey=r.inchikey, smiles=r.smiles, kind=r.kind,
                                       fraction_of_total_mass=r.fraction, fallback_mass=e["fallback_R"]))
        e["R_missing_reason"] = why
        e["n_d86"] = n_d86
    # eligibility of the two compositions
    for eid, e in ents.items():
        e["M_ok"] = e["M"] is not None and e["unmapped_M"] <= MAX_UNMAPPED
        e["R_ok"] = e["R"] is not None and e["fallback_R"] <= MAX_UNMAPPED
        if e["R"] is not None and not e["R_ok"]:
            e["R_missing_reason"] = f"r2 fallback mass {e['fallback_R']:.3f} > {MAX_UNMAPPED}"

    # ---------------- thermal-properties table over the UNION of M and R compounds (freezing)
    all_idx = set()
    for e in ents.values():
        for c in (e["M"], e["R"], e["R_lo"]):
            if c is not None:
                all_idx |= set(int(i) for i in c[0])
    iks = [reg.ik[i] for i in sorted(all_idx)]
    info = pd.DataFrame(index=iks)
    info["smiles"] = [reg.smiles[i] for i in sorted(all_idx)]
    mols = [Chem.MolFromSmiles(s) for s in info.smiles]
    info["cls"] = [g1.hc_class(m) for m in mols]
    info["n_carbon"] = [float(sum(a.GetSymbol() == "C" for a in m.GetAtoms())) for m in mols]
    tmap = get_pure_thermal_properties(info)
    ctx = Ctx(reg, prov, vlb.refutas_nu, g1.ideal_sle_C, tmap)

    # ---------------- scoring
    parts = []
    for eid, e in ents.items():
        if e["targets"].empty:
            continue
        t = e["targets"].copy()
        sc = score_entity(ctx, e["M"] if e["M_ok"] else None, e["R"] if e["R_ok"] else None, t,
                          r_comp_lo=e["R_lo"] if e["R_ok"] else None)
        sc.insert(0, "entity", eid)
        parts.append(sc)
    d = pd.concat(parts, ignore_index=True)
    d["source"] = d.entity.str.split(":").str[0]
    d["label"] = d.entity.map({k: v["label"] for k, v in ents.items()})
    d["pair_ok"] = np.isfinite(d.err_M) & np.isfinite(d.err_R)

    # ---------------- double-count groups (key built before pooling; double_count_group columns are constant per source and unused)
    uf = UF(list(ents))
    links = []
    for a, e in ents.items():
        if e["canon"] != a and e["canon"] in ents:
            uf.union(a, e["canon"]); links.append((a, e["canon"], "S4 dup_map"))
    byl: dict[str, list] = {}
    for a, e in ents.items():
        byl.setdefault(e["label"], []).append(a)
    for lab, mem in byl.items():
        for b in mem[1:]:
            uf.union(mem[0], b); links.append((mem[0], b, f"same POSF label {lab}"))
    byh: dict[str, list] = {}
    for a, e in ents.items():
        if e["comp_hash"]:
            byh.setdefault(e["comp_hash"], []).append(a)
    for h, mem in byh.items():
        for b in mem[1:]:
            uf.union(mem[0], b); links.append((mem[0], b, f"identical mapped composition {h}"))

    def feats(eid):
        e, t = ents[eid], ents[eid]["targets"]
        f = {}
        for k, tl in (("dens15", "density_15C"), ("flash", "flash"), ("freezing", "freezing")):
            v = t.loc[t.target == tl, "measured"]
            if len(v):
                f[k] = float(v.iloc[0])
        r = e["desc"]
        for k, rk in (("aromatics", "aromatics_vol_pct"), ("t10", "d86_t10"), ("t90", "d86_t90")):
            if rk in r:
                f[k] = r[rk]
        return f
    F = {a: feats(a) for a in ents}
    pair_rows = []
    keys = list(ents)
    for i, a in enumerate(keys):
        for b in keys[i + 1:]:
            common = [k for k in DUP_TOL if k in F[a] and k in F[b]]
            if len(common) < DUP_MIN_FEATURES:
                continue
            within = all(abs(F[a][k] - F[b][k]) <= DUP_TOL[k] for k in common)
            if within:
                uf.union(a, b); links.append((a, b, "near-identical bulk (" + ",".join(common) + ")"))
            pair_rows.append(dict(a=a, b=b, n_features=len(common), near_identical=within,
                                  max_scaled_diff=max(abs(F[a][k] - F[b][k]) / DUP_TOL[k] for k in common)))
    # composition distance AV-33-22 x AJFTD (corroboration only)
    nm = {"n-alkane": "n_paraffin", "iso-alkane": "iso_paraffin", "monocycloalkane": "monocycloparaffin", "dicycloalkane": "dicycloparaffin",
          "tricycloalkane": "tricycloparaffin", "alkylbenzene": "alkylbenzene", "cycloaromatic": "cycloaromatic", "diaromatic": "diaromatic", "alkene": "alkene"}
    cajf = pd.read_csv(V5 / "ajftd_gcxgc_composition.csv").groupby(["fuel_id", "hc_class", "n_carbon"]).mass_pct.sum()
    vec = {}
    for fid in comp.fuel_id.unique():
        vec[f"AV-33-22:{fid}"] = {k: v for k, v in comp[comp.fuel_id == fid].assign(cls=lambda x: x.group.map(nm)).groupby(["cls", "carbon_number"]).mass_pct.sum().items()}
    for fid in gfu.fuel_id:
        vec[f"AJFTD:{fid}"] = {k: v for k, v in cajf.xs(fid, level=0).items()}
    comp_dist = []
    for a in [k for k in vec if k.startswith("AV-33-22")]:
        for b in [k for k in vec if k.startswith("AJFTD")]:
            ks = set(vec[a]) | set(vec[b])
            tv = 0.5 * sum(abs(vec[a].get(k, 0) - vec[b].get(k, 0)) for k in ks) / 100.0
            comp_dist.append(dict(av=a, ajftd=b, ajftd_label=ents[b]["label"] if b in ents else labels[int(b.split(":")[1])], tv_distance=tv))
    cd = pd.DataFrame(comp_dist).sort_values("tv_distance")
    cd.head(25).to_csv(OUT / "s3_av3322_vs_ajftd_composition_distance_top25.csv", index=False, float_format="%.5g")
    pd.DataFrame(pair_rows).sort_values("max_scaled_diff").head(60).to_csv(OUT / "s3_bulk_overlap_pairs_top60.csv", index=False, float_format="%.5g")
    pd.DataFrame(links, columns=["entity_a", "entity_b", "reason"]).to_csv(OUT / "s3_group_links.csv", index=False)

    grp = {a: uf.find(a) for a in ents}
    uf2 = UF(list(ents))
    for a, b, why in links:
        if not why.startswith("near-identical"):
            uf2.union(a, b)
    grp2 = {a: uf2.find(a) for a in ents}
    for a, e in ents.items():
        e["group"] = grp[a]
    locked_groups = {grp[a] for a, e in ents.items() if e["locked"]}
    npair = d[d.pair_ok].groupby("entity").size()
    erows = []
    for a, e in ents.items():
        erows.append(dict(entity=a, source=e["source"], fuel_id=e["fid"], label=e["label"], group=e["group"], locked_ajftd=e["locked"],
                          group_locked=e["group"] in locked_groups, n_bins=e["n_bins"], M_unmapped=e["unmapped_M"], M_ok=e["M_ok"],
                          n_components_M=0 if e["M"] is None else len(e["M"][0]),
                          R_available=e["R"] is not None, R_fallback=e["fallback_R"], R_ok=e["R_ok"], R_missing_reason=e["R_missing_reason"],
                          n_d86_points=e["n_d86"], aromatics_vol_pct=e["desc"].get("aromatics_vol_pct", np.nan),
                          diaromatics_vol_pct=e["desc"].get("diaromatics_vol_pct", np.nan), prior_variant=e["prior_variant"],
                          n_target_rows=len(e["targets"]), n_pair_rows=int(npair.get(a, 0)), pairable=bool(e["M_ok"] and e["R_ok"] and npair.get(a, 0) > 0)))
    ent = pd.DataFrame(erows)
    ent.to_csv(OUT / "s3_fuels.csv", index=False, float_format="%.6g")
    pd.DataFrame(r_rows).to_csv(OUT / "s3_r_compositions.csv", index=False, float_format="%.6g")
    sel = select_one_per_group(ent)
    ent2 = ent.assign(group=ent.entity.map(grp2))
    ent2["group_locked"] = ent2.group.isin(set(ent2[ent2.locked_ajftd].group))
    sel2 = select_one_per_group(ent2)
    # groups table
    gt = ent.groupby("group").agg(members=("entity", lambda s: ";".join(s)), n_members=("entity", "size"), sources=("source", lambda s: "+".join(sorted(set(s)))),
                                  any_locked=("group_locked", "max"), n_pairable=("pairable", "sum")).reset_index()
    gt["chosen"] = gt.group.map(dict(zip(sel.group, sel.entity))).fillna("")
    gt.to_csv(OUT / "s3_groups.csv", index=False)

    d = d.merge(ent[["entity", "group", "pairable", "group_locked", "locked_ajftd"]], on="entity", how="left")
    d["headline"] = d.entity.isin(set(sel.entity))
    d["err_Rlo_pair"] = d.err_Rlo
    d.drop(columns=["target_dummy"], errors="ignore").to_csv(OUT / "s3_predictions.csv", index=False, float_format="%.6g")

    # ---------------- central table
    P = d[d.pairable]
    sets = {"headline (one per group, locked groups excluded)": d[d.headline],
            "AJFTD (all pairable)": P[P.source == "AJFTD"], "AV-33-22 (all pairable)": P[P.source == "AV-33-22"], "FuelLib (all pairable)": P[P.source == "FuelLib"],
            "AJFTD S4-locked (separate)": P[P.locked_ajftd]}
    rows = []
    for nm_, s in sets.items():
        rows += central_rows(s, nm_)
    hl = d[d.headline].copy()
    hl_lo = hl.assign(err_R=hl.err_Rlo)
    sens = pd.DataFrame(central_rows(hl_lo, "headline, leave-out prior sensitivity"))
    # drop-group sensitivity: remove groups whose label is in the prior set
    inprior = set(ent[ent.prior_variant != "standard"].group)
    sens2 = pd.DataFrame(central_rows(hl[~hl.group.isin(inprior)], "headline minus groups in prior set"))
    d2 = d[d.entity.isin(set(sel2.entity))].assign(group=lambda x: x.entity.map(grp2))
    sens3 = pd.DataFrame(central_rows(d2, "headline with strict groups (label/dup_map/identical composition only; no bulk-match links)"))
    cen = pd.DataFrame(rows + sens.to_dict("records") + sens2.to_dict("records") + sens3.to_dict("records"))
    cen.to_csv(OUT / "s3_central_table.csv", index=False, float_format="%.5g")

    # ---------------- readings
    rd = []
    H = cen[cen.set.str.startswith("headline (")].set_index("target")
    for tag, Hc in (("headline", H), ("sens_leaveout", sens.set_index("target")), ("sens_dropgroups", sens2.set_index("target"))):
        if "density_15C" in Hc.index and Hc.loc["density_15C"].get("n_groups", 0) > 0:
            r = Hc.loc["density_15C"]
            ok = abs(r.bias_M) <= 7 and abs(r.bias_R) >= 10
            rd.append(dict(reading="R-S3-2", variant=tag, n_groups=int(r.n_groups), bias_M=r.bias_M, bias_R=r.bias_R,
                           met=bool(ok), note=f"M side {'ok' if abs(r.bias_M) <= 7 else 'FAILED'} (|bias_M|<=7); R side {'ok' if abs(r.bias_R) >= 10 else 'FAILED'} (|bias_R|>=10)"))
        else:
            rd.append(dict(reading="R-S3-2", variant=tag, met="not evaluable", note="no paired density rows"))
        if "flash" in Hc.index and Hc.loc["flash"].get("n_groups", 0) > 0:
            r = Hc.loc["flash"]
            ok = r.frac_over5_M <= 0.5 * r.frac_over5_R and r.bias_M <= 2.0
            rd.append(dict(reading="R-S3-3", variant=tag, n_groups=int(r.n_groups), bias_M=r.bias_M, bias_R=r.bias_R, frac_over5_M=r.frac_over5_M,
                           frac_over5_R=r.frac_over5_R, met=bool(ok)))
        else:
            rd.append(dict(reading="R-S3-3", variant=tag, met="not evaluable", note="no paired flash rows"))
        for tl in ("visc_-20C", "visc_-40C"):
            if tl in Hc.index and Hc.loc[tl].get("n_groups", 0) > 0:
                r = Hc.loc[tl]
                rd.append(dict(reading="R-S3-5", variant=tag, target=tl, n_groups=int(r.n_groups), bias_M=r.bias_M, bias_R=r.bias_R,
                               bias_M_pct=r.bias_M_pct, bias_R_pct=r.bias_R_pct, mae_M=r.mae_M, mae_R=r.mae_R,
                               met=("R worse than M" if r.mae_R > r.mae_M else "R not worse than M") if r.n_groups >= (1 if tl == "visc_-20C" else MIN_GROUPS) else "n too small for -40C"))
        for tl in H.index:
            if tag == "headline":
                r = H.loc[tl]
                rd.append(dict(reading="R-S3-1", variant=tag, target=tl, n_groups=int(r.n_groups), met=r.verdict,
                               bias_M=r.get("bias_M"), bias_R=r.get("bias_R"), mae_M=r.get("mae_M"), mae_R=r.get("mae_R"),
                               diff=r.get("mae_diff_R_minus_M"), ci_lo=r.get("ci_lo"), ci_hi=r.get("ci_hi")))
            else:
                r = Hc.loc[tl]
                rd.append(dict(reading="R-S3-1", variant=tag, target=tl, n_groups=int(r.n_groups), met=r.verdict))
    # R-S3-4
    for src, g in ent.groupby("source"):
        rd.append(dict(reading="R-S3-4", variant="M unmapped", source=src, n_fuels=len(g), mean=g.M_unmapped.mean(), max=g.M_unmapped.max(),
                       excluded=";".join(g[~g.M_ok].label.astype(str)) or "none"))
        gg = g[g.R_available]
        rd.append(dict(reading="R-S3-4", variant="R fallback (r2)", source=src, n_fuels=len(gg), mean=gg.R_fallback.mean() if len(gg) else np.nan,
                       max=gg.R_fallback.max() if len(gg) else np.nan, excluded=";".join(gg[gg.R_fallback > MAX_UNMAPPED].label.astype(str)) or "none"))
        rd.append(dict(reading="R-S3-4", variant="R missing", source=src, n_fuels=int((~g.R_available).sum()),
                       excluded=";".join(f"{r.label}({r.R_missing_reason})" for r in g[~g.R_available].itertuples())[:600] or "none"))
    pd.DataFrame(rd).to_csv(OUT / "s3_readings.csv", index=False, float_format="%.5g")

    pd.set_option("display.width", 250); pd.set_option("display.max_columns", 40); pd.set_option("display.max_colwidth", 70)
    print("strict groups:", len(set(grp2.values())), "chosen:", len(sel2))
    print(ent.groupby("source")[["M_ok", "R_available", "R_ok", "pairable"]].sum())
    print("groups:", len(gt), "multi-member groups:", int((gt.n_members > 1).sum()), "locked groups:", int(gt.any_locked.sum()), "chosen:", len(sel))
    print(gt[gt.n_members > 1].to_string())
    print(cen.to_string())
    print(pd.DataFrame(rd).to_string())


if __name__ == "__main__":
    main()
