"""
fuel_hindcast - test the mixture models on REAL many-component fuels (density + kinematic viscosity)
====================================================================================================
Install at fuelblend_data/src/fuelblend/mixture/fuel_hindcast.py

Two ways to build a composition for each fuel
  bulk route  aromatics / olefins / diaromatics vol% (D1319/D6379) + D86 curve -> carbon-number distribution
              via n-alkane boiling points; saturates split into n / iso / cyclo by --sat-split; every
              (class, carbon number) bin is represented by one registry compound.  Up to ~40 components.
              --mc N redraws the saturate split from a Dirichlet and reports the prediction spread.
  GC route    speciated compositions (smiles + wt%): compounds in the registry are used directly, the rest are
              lumped onto the nearest (class, carbon number) representative; matched / lumped / dropped mass
              is printed.
Predictions: physics baseline (ideal volume + Refutas) vs every --bundle (DeepSets / pairwise, from T8).
Flash is NOT evaluated here (class lumps cannot resolve light ends - your T4 audit); use the D86 flash model.

INPUT  --fuels  fuels.csv  (YOU must export this from AJFTD / CRC; I cannot see their layout), one row per fuel:
    fuel_id, aromatics_vol_pct, [olefins_vol_pct], [diaromatics_vol_pct],
    d86_ibp, d86_t5, d86_t10, d86_t20, ... d86_t90, d86_t95, d86_fbp         (degC; any subset, >= 4 points)
    density_15C_kg_m3, density_20C_kg_m3, ...        (measured targets: density_<T>C_kg_m3)
    kv_m40C_cSt, kv_m20C_cSt, kv_20C_cSt, kv_40C_cSt  (measured: kv_<T>C_cSt, 'm' prefix = minus)
  --gc-csv  columns (auto-detected): fuel id, smiles, wt% / mass fraction; fuel ids must match fuels.csv.

    python -m fuelblend.mixture.fuel_hindcast --root . --fuels fuels.csv --bundle results/pairwise/extrap_deepsets_bundle.pt \
        results/pairwise/extrap_pairwise_bundle.pt --mc 16
Also writes T9 inputs:  --write-pool pseudo_pool.csv  and  --write-reference ref.csv --reference-fuel ID [--route gc]
"""
from __future__ import annotations

import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd
import torch

from fuelblend.mixture.cold_flow import CLASSES, compound_meta
from fuelblend.mixture.common import (CompoundRegistry, PurePropertyProvider, RegistryTensors, assemble_batch,
                                      ctx_np, physics_baselines)
from fuelblend.mixture.t8_pairwise import load_any_bundle

NALKANE_TB_C = {5: 36.1, 6: 68.7, 7: 98.4, 8: 125.7, 9: 150.8, 10: 174.1, 11: 195.9, 12: 216.3, 13: 235.4,
                14: 253.6, 15: 270.6, 16: 287.0, 17: 302.0, 18: 316.0, 19: 329.7, 20: 343.0}
D86 = list(zip(["ibp", "t5", "t10", "t20", "t30", "t40", "t50", "t60", "t70", "t80", "t90", "t95", "fbp"],
               [0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 95, 100]))
RHO_CLASS = {"aromatic": 0.87, "olefin": 0.78, "sat": 0.79}      # g/cm3, for vol% -> mass% (assumption)


# --------------------------------------------------------------------------- chemistry helpers
def classify_smiles(smi):
    from rdkit import Chem
    m = Chem.MolFromSmiles(smi)
    if m is None or any(a.GetSymbol() not in ("C", "H") for a in m.GetAtoms()):
        return None
    nc = sum(a.GetSymbol() == "C" for a in m.GetAtoms())
    if any(a.GetIsAromatic() for a in m.GetAtoms()):
        return "aromatic", nc
    unsat = any(b.GetBondType() != Chem.BondType.SINGLE for b in m.GetBonds())
    if m.GetRingInfo().NumRings() > 0:
        return ("alkene_alkyne" if unsat else "cycloalkane"), nc
    if unsat:
        return "alkene_alkyne", nc
    return ("n_paraffin" if all(a.GetDegree() <= 2 for a in m.GetAtoms()) else "iso_paraffin"), nc


class Reps:
    """One registry compound (with valid pure properties) per (class, carbon number) bin."""

    def __init__(self, reg, cls, nc, provider, counts=None):
        cand = [i for i in range(len(reg.ik)) if reg.usable[i] and reg.has_flash[i] and cls[i] in CLASSES
                and np.isfinite(nc[i])]
        provider.ensure(cand)
        self.by = {}
        for i in (i for i in cand if provider.valid[i]):
            self.by.setdefault((cls[i], int(nc[i])), []).append(i)
        if counts is not None:
            for k in self.by:
                self.by[k].sort(key=lambda i: -counts[i])
        self.ns = {c: sorted(n for (cc, n) in self.by if cc == c) for c in CLASSES}

    def get(self, c, n, max_gap=3):
        ns = self.ns.get(c, [])
        if not ns:
            return None
        best = min(ns, key=lambda k: (abs(k - n), k))
        return self.by[(c, best)][0] if abs(best - n) <= max_gap else None


# --------------------------------------------------------------------------- compositions
def carbon_distribution(row):
    pts = [(float(row[f"d86_{k}"]), p) for k, p in D86 if f"d86_{k}" in row and np.isfinite(row[f"d86_{k}"])]
    if len(pts) < 4:
        return None
    pts.sort()
    T = np.array([q[0] for q in pts])
    F = np.maximum.accumulate(np.array([q[1] for q in pts]) / 100.0)
    ns = np.array(sorted(NALKANE_TB_C))
    tb = np.array([NALKANE_TB_C[n] for n in ns])
    cdf = np.concatenate([[0.0], np.interp((tb[:-1] + tb[1:]) / 2, T, F, left=0.0, right=1.0), [1.0]])
    mass = np.diff(cdf)
    return dict(zip(ns.tolist(), (mass / mass.sum()).tolist())) if mass.sum() > 0 else None


def bulk_items(row, sat_split):
    cd = carbon_distribution(row)
    A = float(row.get("aromatics_vol_pct", np.nan))
    if cd is None or not np.isfinite(A):
        return None
    get = lambda k: float(row[k]) if k in row and np.isfinite(row[k]) else 0.0
    O, D = get("olefins_vol_pct"), min(get("diaromatics_vol_pct"), A)
    S = max(100.0 - A - O, 0.0)
    m = {"aromatic": A * RHO_CLASS["aromatic"], "alkene_alkyne": O * RHO_CLASS["olefin"], "sat": S * RHO_CLASS["sat"]}
    tot = sum(m.values())
    m = {k: v / tot for k, v in m.items()}
    items = {}
    for cls_, share in (("n_paraffin", sat_split[0]), ("iso_paraffin", sat_split[1]), ("cycloalkane", sat_split[2])):
        for n, f in cd.items():
            items[(cls_, n)] = items.get((cls_, n), 0.0) + m["sat"] * share * f
    for n, f in cd.items():
        items[("alkene_alkyne", n)] = items.get(("alkene_alkyne", n), 0.0) + m["alkene_alkyne"] * f
        items[("aromatic", n)] = items.get(("aromatic", n), 0.0) + m["aromatic"] * (1 - D / A if A > 0 else 1) * f
    if D > 0:                                                     # diaromatics -> naphthalene bin (C11)
        items[("aromatic", 11)] = items.get(("aromatic", 11), 0.0) + m["aromatic"] * D / A
    return items


def to_vector(items, reps):
    agg, dropped = {}, 0.0
    for (c, n), mass in items.items():
        i = reps.get(c, n)
        if i is None:
            dropped += mass
        else:
            agg[i] = agg.get(i, 0.0) + mass
    s = sum(agg.values())
    idx = np.array(list(agg), dtype=np.int64)
    return idx, np.array([agg[i] for i in idx]) / s, dropped / max(dropped + s, 1e-12)


def gc_items(df_f, reg, reps, valid, scol, wcol):
    from rdkit import Chem
    direct, lumped, lost, items = 0.0, 0.0, 0.0, {}
    ws = df_f[wcol].to_numpy(float)
    ws = ws / ws.sum()
    for smi, w in zip(df_f[scol], ws):
        m = Chem.MolFromSmiles(smi)
        i = reg.ik2i.get(Chem.MolToInchiKey(m)) if m is not None else None
        if i is not None and valid[i]:
            items[("__direct__", i)] = items.get(("__direct__", i), 0.0) + w
            direct += w
            continue
        c = classify_smiles(smi)
        if c is None or reps.get(*c) is None:
            lost += w
        else:
            items[c] = items.get(c, 0.0) + w
            lumped += w
    return items, direct, lumped, lost


def gc_vector(items, reps):
    agg = {}
    for (c, n), mass in items.items():
        i = n if c == "__direct__" else reps.get(c, n)
        agg[i] = agg.get(i, 0.0) + mass
    idx = np.array(list(agg), dtype=np.int64)
    w = np.array([agg[i] for i in idx])
    return idx, w / w.sum()


# --------------------------------------------------------------------------- prediction
@torch.no_grad()
def predict_batch(bundles, rt, provider, comps, T_K, dev):
    B, N = len(comps), max(len(i) for i, _ in comps)
    idx, w = -np.ones((B, N), np.int64), np.zeros((B, N))
    for r, (i, ww) in enumerate(comps):
        idx[r, :len(i)], w[r, :len(i)] = i, ww
    ln_i, rho_i = provider.lookup(idx, T_K)
    tt = lambda a: torch.as_tensor(a, dtype=torch.float32, device=dev)
    it = torch.as_tensor(idx, device=dev)
    base, valid, _ = physics_baselines(rt, it, tt(w), tt(ln_i), tt(rho_i))
    ctx = torch.as_tensor(np.repeat(ctx_np([T_K], [False]), B, 0), device=dev)
    batch = assemble_batch(rt, it, tt(w), tt(np.nan_to_num(ln_i)), tt(np.nan_to_num(rho_i, nan=750.0)), ctx, base, valid)
    out = {"physics": (base.cpu().numpy(), np.zeros((B, 3)))}
    for name, b in bundles.items():
        P = torch.stack([m(batch) for m in b["models"]])
        out[name] = (P.mean(0).cpu().numpy(), P.std(0).cpu().numpy() if P.shape[0] > 1 else np.zeros((B, 3)))
    return out


def parse_targets(cols):
    t = []
    for c in cols:
        m = re.fullmatch(r"kv_(m?)(\d+(?:\.\d+)?)C_cSt", c)
        if m:
            t.append((c, "kv", (-1 if m.group(1) else 1) * float(m.group(2))))
        m = re.fullmatch(r"density_(m?)(\d+(?:\.\d+)?)C_kg_m3", c)
        if m:
            t.append((c, "density", (-1 if m.group(1) else 1) * float(m.group(2))))
    return t


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=".")
    ap.add_argument("--fuels", required=True)
    ap.add_argument("--bundle", nargs="*", default=[])
    ap.add_argument("--gc-csv", default=None)
    ap.add_argument("--sat-split", nargs=3, type=float, default=[0.30, 0.40, 0.30],
                    help="saturates split n-paraffin / iso-paraffin / cycloalkane (ASSUMPTION - report sensitivity)")
    ap.add_argument("--mc", type=int, default=0, help="Dirichlet redraws of the saturate split for a spread estimate")
    ap.add_argument("--write-pool", default=None)
    ap.add_argument("--write-reference", default=None)
    ap.add_argument("--reference-fuel", default=None)
    ap.add_argument("--route", default="bulk", choices=["bulk", "gc"])
    ap.add_argument("--out", default="results/hindcast")
    ap.add_argument("--device", default="cuda" if torch.cuda.is_available() else "cpu")
    a = ap.parse_args()

    root, dev = Path(a.root), torch.device(a.device)
    reg = CompoundRegistry(root / "data/unified/compounds.csv")
    rt = RegistryTensors(reg, dev)
    provider = PurePropertyProvider(reg, root / "data/cache/pure_grid.npz")
    cdf, cls, nc = compound_meta(root / "data/unified/compounds.csv")
    bundles = {Path(p).stem.replace("_bundle", ""): load_any_bundle(p, rt, dev) for p in a.bundle}
    counts = None
    if a.bundle:
        counts = np.asarray(torch.load(a.bundle[0], map_location="cpu", weights_only=False)["comp_count"])
    reps = Reps(reg, cls, nc, provider, counts)
    print(f"[reps] {sum(len(v) for v in reps.by.values())} candidate compounds in {len(reps.by)} (class, C#) bins")
    if a.write_pool:
        rows = [{"inchikey": reg.ik[v[0]], "smiles": reg.smiles[v[0]], "cls": c, "n_carbon": n}
                for (c, n), v in sorted(reps.by.items())]
        pd.DataFrame(rows).to_csv(a.write_pool, index=False)
        print(f"[pool] wrote {len(rows)} representatives -> {a.write_pool}  (use with T9 --pool {a.write_pool})")

    fuels = pd.read_csv(a.fuels)
    targets = parse_targets(fuels.columns)
    gc = pd.read_csv(a.gc_csv) if a.gc_csv else None
    if gc is not None:
        pick = lambda keys: next(c for c in gc.columns if c.lower() in keys)
        gid, gs = pick({"fuel_id", "fuel", "sample", "sample_id", "name"}), pick({"smiles"})
        gw = pick({"wt_pct", "wt%", "wt", "mass_frac", "mass_fraction", "weight_pct", "wt_percent"})
    sat = tuple(a.sat_split)
    rng = np.random.default_rng(0)
    recs, ref_vec = [], {}
    for _, row in fuels.iterrows():
        fid = row["fuel_id"]
        routes = {}
        items = bulk_items(row, sat)
        if items is not None:
            idx, w, drop = to_vector(items, reps)
            draws = []
            for _ in range(a.mc):
                s = rng.dirichlet(np.array(sat) / 0.08)                    # ~ +-0.05 spread around the default
                it2 = bulk_items(row, tuple(s))
                if it2 is not None:
                    draws.append(to_vector(it2, reps)[:2])
            routes["bulk"] = ([(idx, w)] + draws, f"dropped {drop:.1%}")
        if gc is not None and (gc[gid] == fid).any():
            valid = provider.valid
            it, d, l, x = gc_items(gc[gc[gid] == fid], reg, reps, valid, gs, gw)
            idx, w = gc_vector(it, reps)
            routes["gc"] = ([(idx, w)], f"direct {d:.1%}, lumped {l:.1%}, dropped {x:.1%}")
        for rname, (comps, note) in routes.items():
            ref_vec[(fid, rname)] = comps[0]
            for col, kind, T_C in targets:
                if not np.isfinite(row[col]):
                    continue
                pr = predict_batch(bundles, rt, provider, comps, T_C + 273.15, dev)
                rec = {"fuel_id": fid, "route": rname, "n_comp": len(comps[0][0]), "target": col, "kind": kind,
                       "T_C": T_C, "measured": float(row[col]), "composition_note": note}
                for name, (mu, sd) in pr.items():
                    if kind == "density":
                        rec[name] = float(mu[0, 1])
                        spread = float(mu[1:, 1].std()) if len(mu) > 1 else np.nan
                    else:
                        nu = 1000.0 * np.exp(mu[:, 0]) / mu[:, 1]
                        rec[name] = float(nu[0])
                        spread = float(np.log(nu[1:]).std()) if len(nu) > 1 else np.nan
                    rec[f"{name}_mc_spread"] = spread
                    rec[f"{name}_ens_sd"] = float(sd[0, 1] if kind == "density" else sd[0, 0])
                recs.append(rec)
    res = pd.DataFrame(recs)
    Path(a.out).mkdir(parents=True, exist_ok=True)
    res.to_csv(Path(a.out) / "hindcast_predictions.csv", index=False)

    models = ["physics"] + list(bundles)
    rows = []
    for (route, kind, T), g in res.groupby(["route", "kind", "T_C"]):
        row = {"route": route, "kind": kind, "T_C": T, "n_fuels": len(g), "mean_n_comp": g["n_comp"].mean()}
        for mdl in models:
            e = (g[mdl] - g["measured"]) if kind == "density" else np.log(g[mdl] / g["measured"])
            row[f"{mdl}_{'rmse' if kind == 'density' else 'rmse_ln'}"] = float(np.sqrt((e ** 2).mean()))
            row[f"{mdl}_bias"] = float(e.mean())
        rows.append(row)
    pd.set_option("display.width", 250, "display.max_columns", 40, "display.float_format", "{:.4f}".format)
    print("\n" + pd.DataFrame(rows).to_string(index=False))
    print("\n(density in kg/m3; viscosity as RMSE of ln(pred/measured), 0.05 ~ 5%; bias = mean signed error)")
    if a.write_reference:
        key = next((k for k in ref_vec if k == (a.reference_fuel, a.route)), None)
        if key is None:
            raise SystemExit("--reference-fuel / --route not found among processed fuels")
        idx, w = ref_vec[key]
        pd.DataFrame({"inchikey": [reg.ik[i] for i in idx], "mass_frac": w}).to_csv(a.write_reference, index=False)
        print(f"[reference] wrote {len(idx)}-component reference -> {a.write_reference}")


if __name__ == "__main__":
    main()
