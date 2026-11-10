"""v4 Task 9: ideal vs UNIFAC flash point on (a) the 126 lab mixtures (31 systems, leave-one-system-out) and
(b) v4 real fuels that have a raw_file flash point.

Nothing here is fitted: ideal and UNIFAC are zero-parameter physics (anchored Antoine + LFL + original UNIFAC), so
"leave-one-system-out" reduces to scoring each system with no parameter that has seen it; no threshold, margin or
constant was chosen on these rows.  The locked real-fuel set is scored only as a separate, never-tuned slice.

Outputs: results/v4/F/flash_lab_loso.csv, flash_real_fuels.csv, flash_summary.csv, pure_inputs_used.csv
"""
from __future__ import annotations

import argparse
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
from fuelblend.baselines.flash_unifac import (  # noqa: E402
    UnifacFallbackWarning, _default_registry, flash_point_unifac_detail, unifac_subgroups)

OUT = ROOT / "results" / "v4" / "F"
LAB = ROOT / "results" / "eval" / "flash_loso_predictions.csv"
REG = ROOT / "data" / "processed" / "v2" / "compounds_anchored_vp.csv"
COMP_V4 = ROOT / "data" / "real_fuels" / "v4" / "compositions_v4.csv"
JETCLS = ROOT / "data" / "real_fuels" / "v4" / "jet_class_compositions.csv"
PROPS = ROOT / "data" / "real_fuels" / "v4" / "jet_properties_v4.csv"
SPLITS = ROOT / "splits" / "v4" / "real_fuels_v4_splits.json"
PURE_V4 = ROOT / "data" / "processed" / "v4" / "pure_props_v4.csv"
RNG_SEED = 20261109


def stats(err: np.ndarray) -> dict:
    err = np.asarray(err, float)
    err = err[np.isfinite(err)]
    if err.size == 0:
        return dict(N=0)
    return dict(N=int(err.size), MAE_C=float(np.mean(np.abs(err))), RMSE_C=float(np.sqrt(np.mean(err ** 2))),
                bias_C=float(np.mean(err)), over_gt5C_rate=float(np.mean(err > 5.0)),
                n_over_gt5C=int(np.sum(err > 5.0)), worst_over_C=float(err.max()), worst_under_C=float(err.min()))


def run_detail(smiles, x, use_unifac):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UnifacFallbackWarning)   # recorded via fallback_components, not lost
        return flash_point_unifac_detail(smiles, x, use_unifac=use_unifac)


# --------------------------------------------------------------------------- (a) lab mixtures
def lab_loso() -> pd.DataFrame:
    df = pd.read_csv(LAB)
    assert len(df) == 126 and df.system_id.nunique() == 31
    reg = pd.read_csv(REG).set_index("inchikey")
    rows = []
    for r in df.itertuples():
        iks, x = json.loads(r.comp_inchikeys), np.array(json.loads(r.x_mole), float)
        smi = [reg.loc[k, "smiles"] for k in iks]
        ide, uni = run_detail(smi, x, False), run_detail(smi, x, True)
        rows.append(dict(system_id=r.system_id, n_comp=len(iks), measured_C=r.value, method=r.method,
                         source=r.sources, comp_smiles=json.dumps(smi), x_mole=json.dumps(list(x)),
                         ideal_C=ide["flash_C"], unifac_C=uni["flash_C"], ideal_status=ide["status"],
                         unifac_status=uni["status"], n_unifac_fallback=len(uni["fallback_components"]),
                         gamma_min_at_flash=float(np.nanmin(uni["gamma_at_flash"])) if uni["status"] == "ok" else np.nan,
                         gamma_max_at_flash=float(np.nanmax(uni["gamma_at_flash"])) if uni["status"] == "ok" else np.nan))
    out = pd.DataFrame(rows)
    out["ideal_err_C"] = out.ideal_C - out.measured_C
    out["unifac_err_C"] = out.unifac_C - out.measured_C
    out["unifac_minus_ideal_C"] = out.unifac_C - out.ideal_C
    return out


def system_bootstrap_delta(lab: pd.DataFrame, B=5000) -> dict:
    """Paired system-cluster bootstrap of MAE(unifac) - MAE(ideal) and of the >5 C over-prediction rates."""
    rng = np.random.default_rng(RNG_SEED)
    g = {s: d for s, d in lab.groupby("system_id")}
    ids = np.array(list(g))
    dm, do = [], []
    for _ in range(B):
        d = pd.concat([g[s] for s in rng.choice(ids, len(ids))])
        dm.append(np.mean(np.abs(d.unifac_err_C)) - np.mean(np.abs(d.ideal_err_C)))
        do.append(np.mean(d.unifac_err_C > 5) - np.mean(d.ideal_err_C > 5))
    q = lambda v: [float(np.percentile(v, 2.5)), float(np.percentile(v, 97.5))]
    return dict(dMAE_ci95=q(dm), dOver_ci95=q(do), B=B)


# --------------------------------------------------------------------------- (b) real fuels
def mass_to_mole(smiles, w):
    mw = np.array([Descriptors.MolWt(Chem.MolFromSmiles(s)) for s in smiles])
    n = np.asarray(w, float) / mw
    return n / n.sum()


def representative_registry() -> pd.DataFrame:
    reg = _default_registry().reset_index()
    return reg[reg.antoine_A.notna() & reg.lfl_vol_pct.notna()].sort_values("inchikey")


def posf_items(rng=None):
    """POSF-10325 class x carbon-number pseudo-components -> one registry molecule each (nearest cn within class,
    |dcn| <= 1).  rng=None picks the first registry entry by InChIKey; an rng draws a random one (sensitivity)."""
    j = pd.read_csv(JETCLS)
    comps = json.loads(j[j.fuel_id == "POSF-10325"].iloc[0].components)
    reg = representative_registry()
    smiles, w, dropped = [], [], 0.0
    for c in comps:
        cand = reg[(reg.cls == c["class"]) & ((reg.n_carbon - c["cn"]).abs() <= 1)]
        if cand.empty:
            dropped += c["frac"]
            continue
        best = cand[(cand.n_carbon - c["cn"]).abs() == (cand.n_carbon - c["cn"]).abs().min()]
        pick = best.iloc[0] if rng is None else best.iloc[int(rng.integers(len(best)))]
        smiles.append(pick.can)
        w.append(c["frac"])
    return smiles, np.array(w), dropped


def real_fuels() -> pd.DataFrame:
    props = pd.read_csv(PROPS)
    fl = props[(props.property == "flash_point_C") & (props.provenance == "raw_file")]
    sp = json.load(open(SPLITS))
    role, family = sp["fuel_role"], sp["fuel_source_family"]
    comp = pd.read_csv(COMP_V4)
    rows = []
    for r in fl.itertuples():
        fid = r.fuel_id
        meth = "D93" if "D93" in str(r.method) else ("D56" if "D56" in str(r.method) else str(r.method))
        base = dict(fuel_id=fid, method=meth, measured_C=r.value, source_file=r.source_file, flash_provenance=r.provenance,
                    split_role=role.get(fid, "not_in_split"), family=family.get(fid, "jet_gc" if fid == "POSF-10325" else "?"))
        if fid.startswith("CRC647"):
            d = comp[(comp.system_id == f"fuel:{fid}") & ~comp.is_lumped.astype(bool) & (comp.inchikey != "LUMP_NO_REALISTIC_REP")]
            if d.empty:
                continue
            smi, w = list(d.smiles), d.fraction.values
            x = mass_to_mole(smi, w)
            ide, uni = run_detail(smi, x, False), run_detail(smi, x, True)
            fb = uni["fallback_components"]
            base.update(composition_provenance=",".join(sorted(set(d.provenance))), composition_basis="mass->mole (RDKit MW)",
                        n_components=len(smi), unexplained_fraction_dropped=float(d.unexplained_fraction.iloc[0]),
                        mass_frac_unifac_fallback=float(w[fb].sum() / w.sum()) if fb else 0.0,
                        ideal_C=ide["flash_C"], unifac_C=uni["flash_C"], ideal_status=ide["status"], unifac_status=uni["status"],
                        statistics_allowed=True)
        elif fid == "POSF-10325":
            smi, w, dropped = posf_items()
            x = mass_to_mole(smi, w)
            ide, uni = run_detail(smi, x, False), run_detail(smi, x, True)
            rng = np.random.default_rng(RNG_SEED)
            draws = []
            for _ in range(20):
                s2, w2, _ = posf_items(rng)
                x2 = mass_to_mole(s2, w2)
                draws.append((run_detail(s2, x2, False)["flash_C"], run_detail(s2, x2, True)["flash_C"]))
            dr = np.array(draws)
            fb = uni["fallback_components"]
            base.update(composition_provenance="class x carbon-number, formula_derived (D2887 same-fuel, R18); molecules = registry representatives",
                        composition_basis="mass->mole", n_components=len(smi), unexplained_fraction_dropped=dropped,
                        mass_frac_unifac_fallback=float(w[fb].sum() / w.sum()) if fb else 0.0,
                        ideal_C=ide["flash_C"], unifac_C=uni["flash_C"], ideal_status=ide["status"], unifac_status=uni["status"],
                        ideal_C_rep_draw_min=float(dr[:, 0].min()), ideal_C_rep_draw_max=float(dr[:, 0].max()),
                        unifac_C_rep_draw_min=float(dr[:, 1].min()), unifac_C_rep_draw_max=float(dr[:, 1].max()),
                        statistics_allowed=False)
        else:
            continue
        rows.append(base)
    out = pd.DataFrame(rows)
    out["ideal_err_C"] = out.ideal_C - out.measured_C
    out["unifac_err_C"] = out.unifac_C - out.measured_C
    return out


# --------------------------------------------------------------------------- pure-input provenance
def pure_inputs(smiles_set) -> pd.DataFrame:
    reg = _default_registry()
    pv = pd.read_csv(PURE_V4)
    pv["can"] = pv.smiles.map(lambda s: Chem.MolToSmiles(Chem.MolFromSmiles(s)) if isinstance(s, str) else None)
    rows = []
    for s in sorted(smiles_set):
        can = Chem.MolToSmiles(Chem.MolFromSmiles(s))
        r = reg.loc[can]
        A, B, C = r.antoine_A, r.antoine_B, r.antoine_C
        p378 = 10 ** (A - B / (310.95 + C))
        ps = pv[(pv.can == can) & (pv.prop == "Psat37_8")]
        fl = pv[(pv.can == can) & (pv.prop == "flash")]
        rows.append(dict(smiles=can, antoine_source=r.antoine_source, lfl_vol_pct=r.lfl_vol_pct, lfl_source=r.lfl_source,
                         psat37_8_anchored_kPa=p378,
                         pure_v4_psat37_8_status=ps.status.iloc[0] if len(ps) else "no_v4_row",
                         pure_v4_psat37_8_ref_kPa=ps.ref_value.iloc[0] if len(ps) else np.nan,
                         anchored_over_ref=p378 / ps.ref_value.iloc[0] if len(ps) and pd.notna(ps.ref_value.iloc[0]) else np.nan,
                         pure_v4_flash_status=fl.status.iloc[0] if len(fl) else "no_v4_row",
                         unifac_groups=json.dumps(unifac_subgroups(can))))
    return pd.DataFrame(rows)


def main(argv=None):
    global COMP_V4, OUT
    ap = argparse.ArgumentParser()
    ap.add_argument("--compositions", type=Path, default=COMP_V4, help="composition table (default: the v4 file)")
    ap.add_argument("--out", type=Path, default=OUT, help="output directory (default: results/v4/F)")
    a = ap.parse_args(argv)
    COMP_V4, OUT = a.compositions, a.out
    OUT.mkdir(parents=True, exist_ok=True)
    lab = lab_loso()
    lab.to_csv(OUT / "flash_lab_loso.csv", index=False)
    real = real_fuels()
    real.to_csv(OUT / "flash_real_fuels.csv", index=False)

    summ = []

    def add(scope, slice_, d, model, err):
        summ.append(dict(scope=scope, slice=slice_, model=model, **stats(err)))

    ok = lab[(lab.ideal_status == "ok") & (lab.unifac_status == "ok")]
    for name, d in [("all_126", ok), ("binary", ok[ok.n_comp == 2]), ("ternary", ok[ok.n_comp == 3]),
                    ("any_unifac_fallback", ok[ok.n_unifac_fallback > 0])]:
        for m in ("ideal", "unifac"):
            add("lab_mixtures", f"{name} (n_systems={d.system_id.nunique()})", d, m, d[f"{m}_err_C"].values)
    # per-system unit of analysis (mean error within system)
    ps = ok.groupby("system_id")[["ideal_err_C", "unifac_err_C"]].mean()
    for m in ("ideal", "unifac"):
        add("lab_mixtures", "per_system_mean_error (31 systems)", ps, m, ps[f"{m}_err_C"].values)
    boot = system_bootstrap_delta(ok)

    rf = real[(real.ideal_status == "ok") & (real.unifac_status == "ok") & real.statistics_allowed]
    rf = rf.assign(role=np.where(rf.split_role == "locked", "locked", "dev"))
    for (meth, role), d in rf.groupby(["method", "role"]):
        for m in ("ideal", "unifac"):
            add("real_fuels_CRC647", f"{meth} {role}", d, m, d[f"{m}_err_C"].values)
    for meth, d in rf.groupby("method"):
        for m in ("ideal", "unifac"):
            add("real_fuels_CRC647", f"{meth} all (dev+locked, reported separately above)", d, m, d[f"{m}_err_C"].values)
    s = pd.DataFrame(summ)
    s.to_csv(OUT / "flash_summary.csv", index=False)
    json.dump(boot, open(OUT / "flash_lab_bootstrap.json", "w"), indent=1)

    used = set()
    for col in lab.comp_smiles:
        used.update(json.loads(col))
    comp = pd.read_csv(COMP_V4)
    used.update(comp[comp.system_id.str.startswith("fuel:CRC647") & ~comp.is_lumped.astype(bool) & (comp.inchikey != "LUMP_NO_REALISTIC_REP")].smiles)
    used.update(posf_items()[0])
    pure_inputs(used).to_csv(OUT / "pure_inputs_used.csv", index=False)
    pd.set_option("display.width", 250)
    pd.set_option("display.max_columns", 30)
    print(s.round(3).to_string())
    print(boot)
    print(real[["fuel_id", "method", "split_role", "measured_C", "ideal_C", "unifac_C", "mass_frac_unifac_fallback"]].head(8))
    print(real[real.fuel_id == "POSF-10325"].T)


if __name__ == "__main__":
    main()
