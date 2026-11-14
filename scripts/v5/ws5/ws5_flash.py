"""WS5: GH-GNN-derived Wilson activity coefficients inside Le Chatelier flash.  Readings are declared in
results/v5/BASE/delta_ws5_ghgnn/README.md (committed before this file produced any score).

Stages (each writes into results/v5/BASE/delta_ws5_ghgnn/):
  mixtures   build the four mixture sets and the (cached) T2 mapped compositions           -> mixtures.json, t2_mixtures.json
  pairs      ln gamma_inf = A + B/T for every ordered pair + AD table                      -> pair_cache.csv, ad_pairs.csv
  baseline   gamma = 1 through the new wrapper must equal the stored values (stops on mismatch > 0.01 C)
  run        flash with Wilson (primary, all pairs), Wilson AD-gated (secondary), Margules (lab binaries), mesitylene-fixed Antoine
  report     gate, margins, AD summaries -> rows.csv, summary.csv, margins.csv, readings.json
Usage: python scripts/v5/ws5/ws5_flash.py [--main-checkout PATH] [--stage all]
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import math
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors

ROOT = Path(__file__).resolve().parents[3]
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(HERE))
from fuelblend.baselines.flash import predict_flash_le_chatelier  # noqa: E402
from fuelblend.baselines.flash_unifac import _prepare  # noqa: E402
from fuelblend.uq.conformal import conformal_rank, group_scores, one_sided_margin  # noqa: E402
import wilson as W  # noqa: E402

OUT = ROOT / "results" / "v5" / "BASE" / "delta_ws5_ghgnn"
MAIN_DEFAULT = Path(r"C:\Users\kahaa\OneDrive\Desktop\UGRP_hydrocarbons")
B_BOOT, SEED = 5000, 20261109
T_PAIR_C = (25.0, 75.0)
OVER = (2.0, 5.0)
SETS_GATED = ("S-LAB", "S-CRC", "S-T2")


def canon(s: str) -> str:
    return Chem.MolToSmiles(Chem.MolFromSmiles(s))


def _load_script(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(spec)
    sys.modules[name] = m
    spec.loader.exec_module(m)
    return m


# --------------------------------------------------------------------------- flash with an arbitrary gamma
def flash_with_gamma(smiles, x, registry, gamma_fn=None):
    """Same call chain as ``flash_point_unifac_detail`` (P = 101.325 kPa so lfl_eff == lfl) but with a caller-supplied gamma_fn."""
    xn, ant, lfl, missing = _prepare(list(smiles), x, registry)
    if missing:
        return float("nan"), "missing_antoine_or_lfl"
    return predict_flash_le_chatelier(list(xn), ant, list(lfl), lfl_temp=True, gamma_fn=gamma_fn)


# --------------------------------------------------------------------------- mixture sets
def mass_to_mole(smi, w):
    mw = np.array([Descriptors.MolWt(Chem.MolFromSmiles(s)) for s in smi])
    n = np.asarray(w, float) / mw
    return n / n.sum()


def build_lab(t2f) -> list[dict]:
    from fuelblend.data.loaders_li2014 import load_li2014
    old = pd.read_csv(ROOT / "results/v4/F/flash_lab_loso.csv")
    new = load_li2014()
    pts = pd.read_csv(ROOT / "results/v5/BASE/delta_lab_flash_li2014/points.csv")
    assert len(pts) == len(old) + len(new) == 159
    mes_sys = set(pd.read_csv(ROOT / "results/v5/S2/mesitylene_flash_rerun.csv").query("has_mesitylene").system_id)
    rows = []
    for k, r in enumerate(old.itertuples()):
        rows.append(dict(set="S-LAB", row_id=f"lab{k}", group=r.system_id, measured=float(r.measured_C), smiles=json.loads(r.comp_smiles),
                         x=json.loads(r.x_mole), has_mesitylene=r.system_id in mes_sys, label="dev_old", stored_overlay=float(pts.lfl_overlay[k])))
    for k, r in enumerate(new.itertuples()):
        j = len(old) + k
        assert pts.system_id[j] == r.system_id and abs(pts.measured_C[j] - r.measured_C) < 1e-9
        rows.append(dict(set="S-LI", row_id=f"li{k}", group=r.system_id, measured=float(r.measured_C), smiles=json.loads(r.comp_smiles),
                         x=json.loads(r.x_mole), has_mesitylene=False, label="dev", stored_overlay=float(pts.lfl_overlay[j])))
    return rows


def build_crc() -> list[dict]:
    comp = pd.read_csv(ROOT / "data/real_fuels/v4/compositions_crc647_r3.csv")
    f = pd.read_csv(ROOT / "results/v5/L1/F_r3/flash_real_fuels.csv")
    f = f[f.statistics_allowed.astype(bool) & f.fuel_id.str.startswith("CRC647")]
    rows = []
    for r in f.itertuples():
        d = comp[(comp.system_id == f"fuel:{r.fuel_id}") & ~comp.is_lumped.astype(bool) & (comp.inchikey != "LUMP_NO_REALISTIC_REP")]
        smi, w = list(d.smiles), d.fraction.to_numpy(float)
        rows.append(dict(set="S-CRCL" if r.split_role == "locked" else "S-CRC", row_id=r.fuel_id, group=r.fuel_id, measured=float(r.measured_C),
                         smiles=smi, x=list(mass_to_mole(smi, w)), method=r.method, split_role=r.split_role, stored_overlay=float(r.ideal_C)))
    return rows


def build_t2(t2f, main: Path) -> list[dict]:
    cache = OUT / "t2_mixtures.json"
    if cache.exists():
        return json.loads(cache.read_text())
    from fuelblend.baselines.flash_unifac import _lookup
    from fuelblend.data.reconstruction_r3 import SubReps
    from fuelblend.eval.composition_route import RouteCtx, kind_bins, map_bins
    from fuelblend.mixture.common import CompoundRegistry, PurePropertyProvider
    warnings.simplefilter("ignore")
    reg = CompoundRegistry(main / "data/unified/compounds.csv")
    prov = PurePropertyProvider(reg, main / "data/cache/pure_grid.npz", verbose=False)
    reps = SubReps(reg, prov, vp_lookup=lambda s: _lookup(s, None))
    ctx = RouteCtx(reg, prov)
    V5 = ROOT / "data/real_fuels/v5"
    comps = {"av3322": pd.read_csv(V5 / "crc_av3322_composition.csv"), "fuellib": pd.read_csv(V5 / "fuellib_gcxgc_compositions.csv"),
             "ajftd": pd.read_csv(V5 / "ajftd_gcxgc_composition.csv")}
    inv = {"av3322": "av3322", "fuellib": "fuellib", "ajftd_gcxgc": "ajftd"}
    d = pd.read_csv(ROOT / "results/v5/BASE/rows_t2_flash.csv")
    d["status"] = d.status.fillna("")
    d = d[d.chosen & ~d.low_coverage & (d.status == "ok")]
    rows = []
    for r in d.itertuples():
        src = inv[r.source]
        fid = int(r.fuel_id) if src == "ajftd" else r.fuel_id
        c = comps[src]
        mf = map_bins(kind_bins(src, c[c.fuel_id == fid]), reps)
        keep = mf.w > 0
        idx, w = mf.idx[keep], mf.w[keep]
        w = w / w.sum()
        smi = [reg.smiles[i] for i in idx]
        x = w / np.array([ctx.mw(s) for s in smi])
        x = x / x.sum()
        rows.append(dict(set="S-T2", row_id=r.entity, group=r.group, measured=float(r.measured), smiles=smi, x=list(map(float, x)),
                         source=r.source, method_family=r.method_family, stored_overlay=float(r.pred_lfl_overlay)))
    cache.write_text(json.dumps(rows))
    return rows


# --------------------------------------------------------------------------- pair cache and AD table
def build_pairs(mixtures: list[dict], lib) -> tuple[pd.DataFrame, pd.DataFrame]:
    mols = sorted({canon(s) for m in mixtures for s in m["smiles"]})
    print(f"unique molecules {len(mols)}; ordered pairs {len(mols) * (len(mols) - 1)}")
    sol, slv = [], []
    for a in mols:
        for b in mols:
            if a != b:
                sol.append(a), slv.append(b)
    l1 = lib.predict_ln_gamma_inf(sol, slv, [T_PAIR_C[0]] * len(sol))
    l2 = lib.predict_ln_gamma_inf(sol, slv, [T_PAIR_C[1]] * len(sol))
    T1, T2 = T_PAIR_C[0] + 273.15, T_PAIR_C[1] + 273.15
    Bc = (l1 - l2) / (1 / T1 - 1 / T2)
    Ac = l1 - Bc / T1
    pairs = pd.DataFrame(dict(solute=sol, solvent=slv, A=Ac, B=Bc, ln_g_25C=l1, ln_g_75C=l2))
    ad = []
    for a in mols:
        for b in mols:
            if a != b:
                ad.append(dict(solute=a, solvent=b, **lib.ad_pair(a, b)))
    return pairs, pd.DataFrame(ad)


# --------------------------------------------------------------------------- gamma functions
class Gammas:
    def __init__(self, pairs: pd.DataFrame, ad: pd.DataFrame):
        self.table = W.PairTable({(r.solute, r.solvent): r.A for r in pairs.itertuples()}, {(r.solute, r.solvent): r.B for r in pairs.itertuples()})
        feas = {(r.solute, r.solvent): bool(r.ad_feasible) for r in ad.itertuples()}
        self.feas = feas

    def ad_unordered(self, smiles):
        s = [canon(x) for x in smiles]
        return {frozenset((s[i], s[j])): self.feas[(s[i], s[j])] and self.feas[(s[j], s[i])] for i in range(len(s)) for j in range(i + 1, len(s))}

    def wilson_fn(self, smiles, gated=False, stats=None):
        s = [canon(x) for x in smiles]
        adok = self.ad_unordered(smiles) if gated else None

        def fn(T_K, x):
            lnL, nf, npair = self.table.lnL_matrix(s, T_K, adok)
            if stats is not None:
                stats["fail"], stats["pairs"] = nf, npair
            return W.wilson_gamma(np.asarray(x, float), lnL)
        return fn

    def margules_fn(self, smiles):
        s = [canon(x) for x in smiles]
        assert len(s) == 2

        def fn(T_K, x):
            return W.margules_gamma_binary(x, self.table.lng_inf(s[0], s[1], T_K), self.table.lng_inf(s[1], s[0], T_K))
        return fn


# --------------------------------------------------------------------------- statistics
def cluster_boot_dmae(e_a, e_b, groups, B=B_BOOT, seed=SEED):
    """Paired cluster bootstrap of MAE(b) - MAE(a)."""
    e_a, e_b = np.abs(np.asarray(e_a, float)), np.abs(np.asarray(e_b, float))
    g = pd.factorize(np.asarray(groups).astype(str))[0]
    G = g.max() + 1
    n_g = np.bincount(g, minlength=G).astype(float)
    sa, sb = np.bincount(g, weights=e_a, minlength=G), np.bincount(g, weights=e_b, minlength=G)
    rng = np.random.default_rng(seed)
    cnt = rng.multinomial(G, np.full(G, 1.0 / G), size=B).astype(float)
    d = (cnt @ sb - cnt @ sa) / (cnt @ n_g)
    return float(np.percentile(d, 2.5)), float(np.percentile(d, 97.5))


def margin_and_setter(err, groups, alpha=0.10):
    err, groups = np.asarray(err, float), np.asarray(groups).astype(str)
    gs = pd.Series(err).groupby(groups).max()
    n = len(gs)
    k = conformal_rank(n, alpha)
    if k > n:
        return float("inf"), None
    srt = gs.sort_values(kind="stable")
    return float(srt.iloc[k - 1]), srt.index[k - 1]


def describe(err):
    e = np.asarray(err, float)
    return dict(n=int(len(e)), bias=float(e.mean()), mae=float(np.abs(e).mean()), rmse=float(np.sqrt((e ** 2).mean())), worst_over=float(e.max()),
                over2=float((e > 2.0).mean()), over5=float((e > 5.0).mean()))


# --------------------------------------------------------------------------- driver
def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--main-checkout", type=Path, default=MAIN_DEFAULT)
    ap.add_argument("--stage", default="all", choices=["all", "mixtures", "pairs", "baseline", "run", "report"])
    a = ap.parse_args(argv)
    OUT.mkdir(parents=True, exist_ok=True)
    warnings.simplefilter("ignore")
    t2f = _load_script(ROOT / "scripts/v5/base/cells_t2_flash.py", "cells_t2_flash_mod")
    from fuelblend.eval.composition_route import RouteCtx
    lfl, mes = t2f.lfl_overlay_dict(), t2f.mesitylene_dict()
    reg_ov = RouteCtx(None, None, lfl_override=lfl).flash_reg
    reg_mes = RouteCtx(None, None, lfl_override=lfl, antoine_override=mes).flash_reg

    mixtures = build_lab(t2f) + build_crc() + build_t2(t2f, a.main_checkout)
    (OUT / "mixtures.json").write_text(json.dumps(mixtures))
    print(pd.Series([m["set"] for m in mixtures]).value_counts().to_dict())

    # ---- baseline reproduction (gamma = 1), before any GH-GNN number
    base_rows = []
    for m in mixtures:
        t, st = flash_with_gamma(m["smiles"], m["x"], reg_ov, None)
        base_rows.append((t, st))
    diffs = np.array([abs(t - m["stored_overlay"]) if st == "ok" else np.nan for (t, st), m in zip(base_rows, mixtures)])
    worst = np.nanmax(diffs)
    print("baseline reproduction: max |new - stored| =", worst, " non-ok:", sum(st != "ok" for _, st in base_rows))
    bad = [(m["row_id"], m["set"], d) for m, d in zip(mixtures, diffs) if np.isfinite(d) and d > 0.01]
    (OUT / "baseline_reproduction.json").write_text(json.dumps(dict(max_abs_diff_C=float(worst), n_rows=len(mixtures), n_mismatch_gt_0p01=len(bad),
                                                                    mismatches=[list(map(str, b)) for b in bad[:30]]), indent=1))
    if bad:
        print("MISMATCH", bad[:10])
        sys.exit("baseline reproduction failed; stop (declared rule)")
    if a.stage in ("baseline",):
        return

    import ghgnn_lib as lib
    pc, adc = OUT / "pair_cache.csv", OUT / "ad_pairs.csv"
    if pc.exists() and adc.exists():
        pairs, ad = pd.read_csv(pc), pd.read_csv(adc)
    else:
        pairs, ad = build_pairs(mixtures, lib)
        pairs.to_csv(pc, index=False, float_format="%.8g")
        ad.to_csv(adc, index=False, float_format="%.6g")
    if a.stage == "pairs":
        return
    G = Gammas(pairs, ad)

    # ---- run
    rows = []
    for m, (t0, st0) in zip(mixtures, base_rows):
        smi, x = m["smiles"], m["x"]
        r = dict(set=m["set"], row_id=m["row_id"], group=m["group"], measured=m["measured"], n_comp=len(smi), ideal_C=t0, ideal_status=st0,
                 has_mesitylene=m.get("has_mesitylene", False), method=m.get("method", m.get("method_family", m.get("source", ""))))
        stats = {}
        r["wilson_C"], r["wilson_status"] = flash_with_gamma(smi, x, reg_ov, G.wilson_fn(smi, False, stats))
        r["wilson_fail_pairs"], r["n_pairs"] = stats.get("fail", 0), stats.get("pairs", 0)
        r["wilson_ad_C"], r["wilson_ad_status"] = flash_with_gamma(smi, x, reg_ov, G.wilson_fn(smi, True))
        if m["set"] in ("S-LAB", "S-LI"):
            r["ideal_mesfix_C"], _ = flash_with_gamma(smi, x, reg_mes, None)
            r["wilson_mesfix_C"], _ = flash_with_gamma(smi, x, reg_mes, G.wilson_fn(smi, False))
            if len(smi) == 2:
                r["margules_C"], r["margules_status"] = flash_with_gamma(smi, x, reg_ov, G.margules_fn(smi))
        # pair-weight in domain, and gamma at the ideal flash temperature
        s = [canon(q) for q in smi]
        xs = np.asarray(x, float) / np.sum(x)
        adm = G.ad_unordered(smi)
        wt = [(xs[i] * xs[j], adm[frozenset((s[i], s[j]))]) for i in range(len(s)) for j in range(i + 1, len(s))]
        tot = sum(w for w, _ in wt)
        r["ad_pair_weight_in_domain"] = sum(w for w, ok in wt if ok) / tot if tot > 0 else 1.0
        r["n_pairs_in_domain"] = int(sum(ok for _, ok in wt))
        if np.isfinite(t0):
            gam = G.wilson_fn(smi, False)(t0 + 273.15, xs)
            r["gamma_min_at_ideal_flash"], r["gamma_max_at_ideal_flash"] = float(np.min(gam)), float(np.max(gam))
            r["flash_T_C"] = float(t0)
        rows.append(r)
        if len(rows) % 20 == 0:
            print(len(rows), "/", len(mixtures), flush=True)
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "rows.csv", index=False, float_format="%.8g")
    print("rows written")


if __name__ == "__main__":
    main()
