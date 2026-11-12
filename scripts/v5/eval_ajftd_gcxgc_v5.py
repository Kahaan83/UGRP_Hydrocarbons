"""v5 G1: unfitted physics baselines on AJFTD fuels with a measured GCxGC class x carbon-number composition.

Readings, data and models are declared in results/v5/G1/README.md (written before this script was first run).
Step 1 maps every bin to one representative compound and writes the mapping; step 2 scores.  Nothing is fitted.

Outputs: data/real_fuels/v5/ajftd_gcxgc_mapped.csv, results/v5/G1/{mapping_summary.csv, predictions.csv, summary.csv, RESULTS.md}
"""
from __future__ import annotations

import hashlib
import importlib.util
import json
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from fuelblend.baselines.flash_unifac import UnifacFallbackWarning, _lookup, flash_point_unifac_detail  # noqa: E402
from fuelblend.data.reconstruction_r3 import SubReps  # noqa: E402
from fuelblend.mixture.common import CompoundRegistry, PurePropertyProvider  # noqa: E402
from fuelblend.thermo.coutinho_wax import R_GAS, get_pure_thermal_properties  # noqa: E402

V5 = ROOT / "data" / "real_fuels" / "v5"
OUT = ROOT / "results" / "v5" / "G1"
SPLITS = ROOT / "results" / "v5" / "S4" / "splits_ajftd_v5.json"
MAX_UNMAPPED = 0.10
SHARED_LABELS = {"10325", "10289", "4658", "11769"}
# AJFTD class -> (r3 sub-kind, max carbon-number gap); None = no representative
KIND = {"n_paraffin": ("n_paraffin", 3), "iso_paraffin": ("iso_paraffin", 3), "monocycloparaffin": ("cyclo_mono", 3),
        "alkylbenzene": ("alkylbenzene", 3), "alkene": ("alkene_alkyne", 3), "dicycloparaffin": ("cyclo_di", 99),
        "cycloaromatic": ("indan_tetralin", 99), "diaromatic": ("aromatic_di", 99), "tricycloparaffin": None, "diene": None}


def hc_class(m) -> str:
    if any(a.GetIsAromatic() for a in m.GetAtoms()):
        return "aromatic"
    if any(b.GetBondType() != Chem.BondType.SINGLE for b in m.GetBonds()):
        return "alkene_alkyne"
    if m.GetRingInfo().NumRings() > 0:
        return "cycloalkane"
    return "iso_paraffin" if any(a.GetSymbol() == "C" and a.GetDegree() > 2 for a in m.GetAtoms()) else "n_paraffin"


def map_fuel(bins: pd.DataFrame, reps) -> tuple[dict[int, float], float, dict[str, float]]:
    """bins: rows of one fuel (hc_class, n_carbon, mass_pct) -> ({registry index: mass}, unmapped fraction, {reason: mass})."""
    agg: dict[int, float] = {}
    lost: dict[str, float] = {}
    tot = float(bins.mass_pct.sum())
    for r in bins.itertuples():
        k = KIND.get(r.hc_class)
        i = reps.get(k[0], int(r.n_carbon), max_gap=k[1]) if k else None
        if i is None:
            key = f"no_rep_{r.hc_class}"
            lost[key] = lost.get(key, 0.0) + r.mass_pct / tot
        else:
            agg[i] = agg.get(i, 0.0) + r.mass_pct
    mapped = sum(agg.values())
    return agg, 1.0 - mapped / tot, lost


def build_mapping(reps) -> tuple[pd.DataFrame, pd.DataFrame]:
    comp = pd.read_csv(V5 / "ajftd_gcxgc_composition.csv")
    fuels = pd.read_csv(V5 / "ajftd_gcxgc_fuels.csv")
    ok = fuels[fuels.composition_complete]
    rows, summ = [], []
    for f in ok.itertuples():
        agg, unmapped, lost = map_fuel(comp[comp.fuel_id == f.fuel_id], reps)
        tot = sum(agg.values())
        idx = sorted(agg)
        h = hashlib.md5("|".join(f"{reps.ik[i]}:{agg[i] / tot:.4f}" for i in idx).encode()).hexdigest()[:12] if idx else ""
        summ.append(dict(fuel_id=f.fuel_id, fuel_label=str(f.fuel_label), n_bins=int(f.n_bins), n_components=len(idx),
                         unmapped_fraction=unmapped, eligible=bool(idx) and unmapped <= MAX_UNMAPPED, comp_hash=h,
                         shared_with_other_source=str(f.fuel_label) in SHARED_LABELS, unmapped_detail=json.dumps(lost)))
        for i in idx:
            rows.append(dict(fuel_id=f.fuel_id, fuel_label=str(f.fuel_label), reg_index=i, inchikey=reps.ik[i], smiles=reps.smiles[i],
                             fraction=agg[i] / tot, basis="mass", unmapped_fraction=unmapped, provenance="real_fuel_gcxgc_bins_mapped_r3reps"))
    return pd.DataFrame(rows), pd.DataFrame(summ)


def targets(exclude_censored: bool = False) -> pd.DataFrame:
    """Per fuel/property/T median target.  exclude_censored=False (default) keeps G1's committed behaviour; True drops rows with a
    qualifier (e.g. '<' freezing limits), which are bounds, not measurements."""
    t = pd.read_csv(V5 / "ajftd_fuels_long.csv", low_memory=False)
    t = t[t.property.isin(["flash_point_C", "viscosity_kin_cSt", "density_kg_m3", "freezing_point_C"]) & t.value.notna()]
    if exclude_censored:
        t = t[t.qualifier.isna()]
    t = t[~((t.property == "viscosity_kin_cSt") & (t.T_source.astype(str) != "reported"))]
    t = t[~(t.property.isin(["viscosity_kin_cSt", "density_kg_m3"]) & t.T_K.isna())]
    t["T_key"] = t.T_K.round(1).fillna(-1.0)
    t["method_n"] = [str(m).replace(" ", "") if isinstance(m, str) else "unstated" for m in t.method]
    g = t.groupby(["fuel_id", "property", "T_key"])
    out = g.agg(measured=("value", "median"), n_tests=("value", "size"), spread=("value", lambda v: float(v.max() - v.min())),
                method=("method_n", lambda v: "+".join(sorted(set(v))))).reset_index()
    return out


def ideal_sle_C(x, iks, tmap) -> float:
    best = np.nan
    for xi, ik in zip(x, iks):
        t = tmap.get(ik)
        if t and np.isfinite(t["tm_K"]) and np.isfinite(t["dh_fus_J_mol"]) and xi > 1e-4:
            inv = 1.0 / t["tm_K"] - R_GAS / t["dh_fus_J_mol"] * np.log(xi)
            if inv > 0:
                T = 1.0 / inv - 273.15
                best = T if not np.isfinite(best) else max(best, T)
    return best


def score(mapped: pd.DataFrame, summ: pd.DataFrame, prov, refutas_nu) -> pd.DataFrame:
    tg = targets()
    role = json.load(open(SPLITS))["fuel_role"]
    info = mapped.drop_duplicates("inchikey").set_index("inchikey")[["smiles"]].copy()
    mols = {k: Chem.MolFromSmiles(s) for k, s in info.smiles.items()}
    info["cls"] = [hc_class(mols[k]) for k in info.index]
    info["n_carbon"] = [float(sum(a.GetSymbol() == "C" for a in mols[k].GetAtoms())) for k in info.index]
    tmap = get_pure_thermal_properties(info)
    mw = {k: Descriptors.MolWt(m) for k, m in mols.items()}
    rows = []
    for s in summ[summ.eligible].itertuples():
        g = mapped[mapped.fuel_id == s.fuel_id]
        t = tg[tg.fuel_id == s.fuel_id]
        if t.empty:
            continue
        w = g.fraction.to_numpy(float)
        iks, smi, idx = list(g.inchikey), list(g.smiles), g.reg_index.to_numpy(np.int64)
        x = w / np.array([mw[k] for k in iks])
        x = x / x.sum()
        base = dict(fuel_id=s.fuel_id, fuel_label=s.fuel_label, split_role=role.get(str(s.fuel_id), "not_in_split"),
                    shared_with_other_source=s.shared_with_other_source, comp_hash=s.comp_hash,
                    unmapped_fraction=s.unmapped_fraction, n_components=len(iks))
        for r in t.itertuples():
            rec = dict(base, property=r.property, T_K=np.nan if r.T_key < 0 else r.T_key, measured=r.measured,
                       n_tests=r.n_tests, test_spread=r.spread, method=r.method)
            if r.property == "flash_point_C":
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore", UnifacFallbackWarning)
                    ide = flash_point_unifac_detail(smi, x, use_unifac=False)
                    uni = flash_point_unifac_detail(smi, x, use_unifac=True)
                rows.append(dict(rec, model="ideal", pred=ide["flash_C"], status=ide["status"]))
                rows.append(dict(rec, model="unifac", pred=uni["flash_C"], status=uni["status"]))
            elif r.property in ("viscosity_kin_cSt", "density_kg_m3"):
                nu, rho_mix, _, _ = refutas_nu(prov, idx, w, float(r.T_key))
                val = nu if r.property == "viscosity_kin_cSt" else rho_mix
                rows.append(dict(rec, model="refutas" if r.property == "viscosity_kin_cSt" else "ideal_volume", pred=val, status="ok"))
            else:
                rows.append(dict(rec, model="ideal_svl", pred=ideal_sle_C(x, iks, tmap), status="ok"))
    d = pd.DataFrame(rows)
    visc = d.property == "viscosity_kin_cSt"
    d["err"] = d.pred - d.measured
    d.loc[visc, "err"] = np.log(d.loc[visc, "pred"] / d.loc[visc, "measured"])
    return d


def stat(e: pd.Series, flash=False) -> dict:
    e = e[np.isfinite(e)]
    out = dict(n=int(len(e)), bias=float(e.mean()) if len(e) else np.nan, mae=float(e.abs().mean()) if len(e) else np.nan)
    if flash and len(e):
        out.update(over_gt5_rate=float((e > 5).mean()), worst_over=float(e.max()), worst_under=float(e.min()))
    return out


def summarise(d: pd.DataFrame) -> pd.DataFrame:
    rows = []
    sets = {"all eligible": np.ones(len(d), bool), "not shared with another source": ~d.shared_with_other_source.to_numpy(bool)}
    for role in sorted(d.split_role.unique()):
        sets[f"S4 role: {role}"] = (d.split_role == role).to_numpy()
    d56 = d.method.str.contains("D56")
    for name, m in sets.items():
        for (prop, model), g in d[m].groupby(["property", "model"]):
            per_T = prop in ("viscosity_kin_cSt", "density_kg_m3")
            keys = [(f"{t - 273.15:.0f} C", gg) for t, gg in g.groupby("T_K")] if per_T else [("all T", g)]
            for tk, gg in keys:
                rows.append(dict(set=name, property=prop, model=model, T=tk, **stat(gg.err, flash=prop == "flash_point_C")))
            if prop == "flash_point_C":
                parts = {k: gg for k, gg in g.groupby(d56.loc[g.index].map({True: "D56", False: "D93"}))}
                if len(parts) == 2 and all(p.fuel_id.nunique() >= 5 for p in parts.values()):
                    for meth, gg in parts.items():
                        rows.append(dict(set=name, property=prop, model=model, T=f"method {meth}", **stat(gg.err, flash=True)))
    return pd.DataFrame(rows)


def loo_constant(meas: np.ndarray) -> float:
    n = len(meas)
    return float(np.mean(np.abs((meas.sum() - meas) / (n - 1) - meas))) if n > 2 else np.nan


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    spec = importlib.util.spec_from_file_location("vlb", ROOT / "scripts" / "v4" / "build_viscosity_leaderboard_v4.py")
    vlb = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(vlb)
    reg = CompoundRegistry(ROOT / "data" / "unified" / "compounds.csv")
    prov = PurePropertyProvider(reg, ROOT / "data" / "cache" / "pure_grid.npz", verbose=False)
    reps = SubReps(reg, prov, vp_lookup=lambda s: _lookup(s, None))
    mapped, summ = build_mapping(reps)
    mapped.to_csv(V5 / "ajftd_gcxgc_mapped.csv", index=False)          # written before any score is computed
    summ.to_csv(OUT / "mapping_summary.csv", index=False)
    d = score(mapped, summ, prov, vlb.refutas_nu)
    d.to_csv(OUT / "predictions.csv", index=False)
    s = summarise(d)
    s.to_csv(OUT / "summary.csv", index=False)

    el = summ[summ.eligible]
    L = ["# G1 results (generated by scripts/v5/eval_ajftd_gcxgc_v5.py; readings are in README.md)", "",
         f"Compositions: {len(summ)} complete, {len(el)} eligible (unmapped mass at most {MAX_UNMAPPED:.0%}); median unmapped "
         f"{el.unmapped_fraction.median():.3f}, largest {el.unmapped_fraction.max():.3f}. "
         f"{d.fuel_id.nunique()} eligible fuels have at least one scored target; "
         f"{d[['fuel_id', 'comp_hash']].drop_duplicates().comp_hash.nunique()} distinct mapped compositions among them.",
         f"Not eligible: {', '.join(f'{r.fuel_label} ({r.unmapped_fraction:.2f})' for r in summ[~summ.eligible].itertuples()) or 'none'}.", "",
         "Error = predicted - measured (C, kg/m3); viscosity error = ln(pred/meas).", ""]
    for name in s.set.unique():
        L += [f"## {name}", "", "| property | model | T / method | n | bias | MAE | over by >5 C | worst over | worst under |", "|---|---|---|---|---|---|---|---|---|"]
        for r in s[s.set == name].itertuples():
            ex = (f"{r.over_gt5_rate:.2f} | {r.worst_over:+.1f} | {r.worst_under:+.1f}" if r.property == "flash_point_C" and r.n else " | | ")
            L.append(f"| {r.property} | {r.model} | {r.T} | {r.n} | {r.bias:+.3f} | {r.mae:.3f} | {ex} |")
        L.append("")
    a = s[(s.set == "all eligible") & (s["T"] == "all T")].set_index(["property", "model"])
    L += ["## Readings", ""]
    if ("flash_point_C", "ideal") in a.index:
        fi, fu = a.loc[("flash_point_C", "ideal")], a.loc[("flash_point_C", "unifac")]
        r2 = fi.bias <= 2.0 and fi.over_gt5_rate <= 0.20
        L.append(f"* R2 (ideal flash, all eligible): bias {fi.bias:+.2f} C, over by >5 C {fi.over_gt5_rate:.2f} (n={int(fi.n)}): "
                 f"\"measured composition removes the unsafe bias\" {'HOLDS' if r2 else 'DOES NOT HOLD'}.")
        r3 = fu.mae < fi.mae and fu.over_gt5_rate <= fi.over_gt5_rate
        L.append(f"* R3 (UNIFAC): MAE {fu.mae:.2f} vs ideal {fi.mae:.2f} C; over by >5 C {fu.over_gt5_rate:.2f} vs {fi.over_gt5_rate:.2f}: "
                 f"\"UNIFAC helps\" {'HOLDS' if r3 else 'DOES NOT HOLD'}.")
    fz = d[(d.property == "freezing_point_C") & np.isfinite(d.err)]
    if len(fz) > 2:
        cm = loo_constant(fz.measured.to_numpy(float))
        L.append(f"* R4 (freezing): ideal SvL MAE {fz.err.abs().mean():.2f} C vs leave-one-out constant mean {cm:.2f} C (n={len(fz)}): "
                 f"\"the model has skill\" {'HOLDS' if fz.err.abs().mean() < cm else 'DOES NOT HOLD'}.")
    nz = d[(d.property == "freezing_point_C") & ~np.isfinite(d.err)]
    if len(nz):
        L.append(f"* Freezing not evaluable for {len(nz)} fuels (no component with library Tm and dHfus above the 1e-4 mole-fraction floor).")
    bad = d[(d.property == "flash_point_C") & (d.status != "ok")]
    if len(bad):
        L.append(f"* Flash solver status not ok on {len(bad)} rows: {bad.status.value_counts().to_dict()}.")
    L.append("")
    (OUT / "RESULTS.md").write_text("\n".join(L), encoding="utf-8")
    print("\n".join(L))


if __name__ == "__main__":
    main()
