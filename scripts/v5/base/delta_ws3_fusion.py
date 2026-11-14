"""v5 delta_ws3_pure_gaps Part B: Tm / dHfus gap table, leave-one-family-out estimator (gate R-W3-3), filled values, propagation.

Declared readings: results/v5/BASE/delta_ws3_pure_gaps/README.md (committed before this script existed).
Needs data/processed/v5/thermoml_fusion_pure_hc_by_compound.csv (python -m fuelblend.data.thermoml_fusion).
Run: PYTHONPATH=src .venv/Scripts/python.exe scripts/v5/base/delta_ws3_fusion.py
Writes gap_table.csv, fusion_oof.csv, fusion_family_summary.csv, propagation.csv, readings_fusion.json, and
data/processed/v5/fusion_ml_estimated.csv.  Nothing in the API or design path reads these files.
"""
from __future__ import annotations

import itertools
import json
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors
from sklearn.model_selection import KFold

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "src"))
from fuelblend.eval.composition_route import ideal_sle_C  # noqa: E402
from fuelblend.thermo import fusion_estimate as FE  # noqa: E402

OUT = ROOT / "results" / "v5" / "BASE" / "delta_ws3_pure_gaps"
VERIFIED = ROOT / "data" / "raw" / "sle_constants" / "tm_dhfus_verified.csv"
TML = ROOT / "data" / "processed" / "v5" / "thermoml_fusion_pure_hc_by_compound.csv"
ML_OUT = ROOT / "data" / "processed" / "v5" / "fusion_ml_estimated.csv"
DESIGNS = ROOT / "results" / "v5" / "BASE" / "delta_smoke_cert" / "designs"
POOLS = {"v5pool": "v5pool", "ext": "ext_n30", "poolA": "poolA_toluene_isooctane", "poolB": "poolB_hexylbenzene_isooctane"}
GATE_K = 5.0


def ik(smi):
    return Chem.MolToInchiKey(Chem.MolFromSmiles(smi))


def measured_table() -> pd.DataFrame:
    """One row per InChIKey with tier per quantity: verified first, else thermoml_measured."""
    v = pd.read_csv(VERIFIED)
    t = pd.read_csv(TML)
    rows = {}
    for r in t.itertuples():
        rows[r.inchikey] = dict(inchikey=r.inchikey, smiles=r.smiles, name=r.name,
                                tm_K=r.tm_K, tm_source="thermoml_measured" if np.isfinite(r.tm_K) else None,
                                dh_kJ=r.dhfus_kJ_mol, dh_source="thermoml_measured" if np.isfinite(r.dhfus_kJ_mol) else None)
    for r in v.itertuples():
        e = rows.setdefault(r.inchikey, dict(inchikey=r.inchikey, smiles=r.smiles, name=r.name, tm_K=np.nan, tm_source=None,
                                             dh_kJ=np.nan, dh_source=None))
        e["name"] = r.name
        if np.isfinite(r.tm_K):
            e.update(tm_K=float(r.tm_K), tm_source="verified")
        if np.isfinite(r.svl_dhfus_kJ_mol):
            e.update(dh_kJ=float(r.svl_dhfus_kJ_mol), dh_source="verified")
    return pd.DataFrame(rows.values())


def gap_compounds() -> pd.DataFrame:
    """Pool members (by name, via pure_props_v5) and the API representatives used for the AJFTD GCxGC fuels."""
    pp = pd.read_csv(ROOT / "data/processed/v5/pure_props_v5.csv").drop_duplicates("compound").set_index("compound").smiles
    recs = {}
    for pool, run in POOLS.items():
        names = json.loads(pd.read_csv(DESIGNS / run / "designs.csv").pool.iloc[0])
        for n in names:
            s = pp[n]
            k = ik(s)
            recs.setdefault(k, dict(inchikey=k, smiles=Chem.MolToSmiles(Chem.MolFromSmiles(s)), label=n, used_in=set()))["used_in"].add(pool)
    from fuelblend.api.predict_route import get_ctx
    from fuelblend.data.representatives_v5 import CANDIDATES
    from fuelblend.eval.composition_route import kind_bins, map_bins
    ctx = get_ctx()
    c = pd.read_csv(ROOT / "data/real_fuels/v5/ajftd_gcxgc_composition.csv")
    for fid, g in c.groupby("fuel_id"):
        mf = map_bins(kind_bins("ajftd", g), ctx.reps)
        for i in mf.idx:
            k = ctx.reg.ik[int(i)]
            recs.setdefault(k, dict(inchikey=k, smiles=ctx.reg.smiles[int(i)], label=None, used_in=set()))["used_in"].add("ajftd_gcxgc_rep")
    for n, s in CANDIDATES.items():
        k = ik(s)
        recs.setdefault(k, dict(inchikey=k, smiles=Chem.MolToSmiles(Chem.MolFromSmiles(s)), label=n, used_in=set()))["used_in"].add("subreps_v5_candidate")
    df = pd.DataFrame(recs.values())
    df["used_in"] = df.used_in.map(lambda s: "+".join(sorted(s)))
    return df


def _tier(m, k, col):
    v = m.loc[k, col] if k in m.index else None
    return v if isinstance(v, str) and v else "needs_estimate"


def oof_tables(train: pd.DataFrame):
    X = FE.features(train.smiles)
    fam = train.smiles.map(FE.family_of).to_numpy()
    y_tm = train.tm_K.to_numpy(float)
    p_tm = FE.oof_predict(X, y_tm, fam)
    # random 5-fold by compound (descriptive)
    p_rand = np.full(len(y_tm), np.nan)
    for tr, te in KFold(5, shuffle=True, random_state=0).split(X):
        p_rand[te] = FE.fit_predict(X.iloc[tr], y_tm[tr], X.iloc[te])
    has_dh = np.isfinite(train.dh_kJ.to_numpy(float))
    ds = train.dh_kJ.to_numpy(float) * 1000.0 / y_tm
    p_ds = np.full(len(y_tm), np.nan)
    p_ds[has_dh] = FE.oof_predict(X[has_dh], ds[has_dh], fam[has_dh])
    oof = train[["inchikey", "smiles", "name", "tm_K", "tm_source", "dh_kJ", "dh_source"]].copy()
    oof["family"] = fam
    oof["n_carbon"] = X.n_carbon.to_numpy()
    oof["tm_oof"] = p_tm
    oof["tm_random5"] = p_rand
    oof["tm_joback"] = X.joback_tm.to_numpy()
    oof["dsfus_oof"] = p_ds
    oof["dh_oof_kJ"] = p_ds * p_tm / 1000.0                      # predicted dS x predicted Tm (what a fill would give)
    oof["dh_joback_kJ"] = X.joback_dhfus.to_numpy()
    return oof, X


def summarize(oof: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for scope, d in [("all_hc", oof), ("hc_C8_C20", oof[(oof.n_carbon >= 8) & (oof.n_carbon <= 20)])] + \
            [(f"family:{f}", g) for f, g in oof.groupby("family")]:
        e = d.tm_oof - d.tm_K
        ej = d.tm_joback - d.tm_K
        er = d.tm_random5 - d.tm_K
        dd = d[np.isfinite(d.dh_kJ) & np.isfinite(d.dh_oof_kJ)]
        eh = dd.dh_oof_kJ - dd.dh_kJ
        ehj = (dd.dh_joback_kJ - dd.dh_kJ).dropna()
        rows.append(dict(scope=scope, n=len(d), mae_tm_lofo=e.abs().mean(), rmse_tm_lofo=float(np.sqrt((e ** 2).mean())), bias_tm_lofo=e.mean(),
                         mae_tm_joback=ej.abs().mean(), mae_tm_random5=er.abs().mean(),
                         n_dh=len(dd), mae_dh_lofo_kJ=eh.abs().mean(), rmse_dh_lofo_kJ=float(np.sqrt((eh ** 2).mean())) if len(dd) else np.nan,
                         mape_dh_lofo=float((eh.abs() / dd.dh_kJ).mean() * 100) if len(dd) else np.nan,
                         mae_dh_joback_kJ=ehj.abs().mean() if len(ehj) else np.nan))
    return pd.DataFrame(rows)


def propagation(const: dict, rmse_tm: float, rmse_dh: float) -> pd.DataFrame:
    pp = pd.read_csv(ROOT / "data/processed/v5/pure_props_v5.csv").drop_duplicates("compound").set_index("compound").smiles
    rows = []
    for pool in ("poolA", "poolB"):
        d = pd.read_csv(DESIGNS / POOLS[pool] / "designs.csv")
        for r in d.itertuples():
            w = json.loads(r.weights)
            names = [n for n, v in w.items() if v > 0]
            smi = [pp[n] for n in names]
            iks = [ik(s) for s in smi]
            x = np.array([w[n] / Descriptors.MolWt(Chem.MolFromSmiles(s)) for n, s in zip(names, smi)])
            x = x / x.sum()
            tmap = {k: {"tm_K": const[k]["tm_K"], "dh_fus_J_mol": const[k]["dh_kJ"] * 1000.0} for k in iks if k in const}
            base = ideal_sle_C(x, iks, tmap)
            est = [k for k in iks if k in const and const[k]["est"]]
            worst = 0.0
            for st, sh in itertools.product((-1, 1), (-1, 1)):
                tm2 = {k: dict(v) for k, v in tmap.items()}
                for k in est:
                    tm2[k]["tm_K"] += st * rmse_tm
                    tm2[k]["dh_fus_J_mol"] = max(tm2[k]["dh_fus_J_mol"] + sh * rmse_dh * 1000.0, 100.0)
                v = ideal_sle_C(x, iks, tm2)
                worst = max(worst, abs(v - base)) if np.isfinite(v) and np.isfinite(base) else np.nan
            rows.append(dict(pool=pool, design_id=r.design_id, spec=r.spec_name, components=json.dumps(names),
                             estimated_components=json.dumps([names[iks.index(k)] for k in est]),
                             mole_frac_estimated=float(sum(x[iks.index(k)] for k in est)),
                             freezing_C_tiered=base, max_abs_change_pm1rmse_C=worst,
                             missing_constants=json.dumps([n for n, k in zip(names, iks) if k not in const])))
    return pd.DataFrame(rows)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    meas = measured_table()
    train = meas[np.isfinite(meas.tm_K)].reset_index(drop=True)
    oof, X = oof_tables(train)
    summ = summarize(oof)
    oof.to_csv(OUT / "fusion_oof.csv", index=False)
    summ.to_csv(OUT / "fusion_family_summary.csv", index=False)
    g = summ[summ.scope == "hc_C8_C20"].iloc[0]
    a = summ[summ.scope == "all_hc"].iloc[0]
    gate = bool(g.mae_tm_lofo <= GATE_K)
    rmse_tm, rmse_dh = float(g.rmse_tm_lofo), float(g.rmse_dh_lofo_kJ)

    gaps = gap_compounds()
    m = meas.set_index("inchikey")
    gaps["tm_tier"] = [_tier(m, k, "tm_source") for k in gaps.inchikey]
    gaps["dh_tier"] = [_tier(m, k, "dh_source") for k in gaps.inchikey]
    gaps["family"] = gaps.smiles.map(FE.family_of)
    gaps["n_carbon"] = gaps.smiles.map(FE.n_carbon)
    need = gaps[(gaps.tm_tier == "needs_estimate") | (gaps.dh_tier == "needs_estimate")].reset_index(drop=True)

    # final models on all measured compounds -> estimates for the gap compounds only
    est_rows = []
    if len(need):
        Xn = FE.features(need.smiles)
        tm_hat = FE.fit_predict(X, train.tm_K.to_numpy(float), Xn)
        has_dh = np.isfinite(train.dh_kJ.to_numpy(float))
        ds = train.dh_kJ.to_numpy(float)[has_dh] * 1000.0 / train.tm_K.to_numpy(float)[has_dh]
        ds_hat = FE.fit_predict(X[has_dh], ds, Xn)
        for i, r in need.iterrows():
            tm_meas = m.loc[r.inchikey, "tm_K"] if r.tm_tier != "needs_estimate" else np.nan
            tm_use = tm_meas if np.isfinite(tm_meas) else tm_hat[i]
            est_rows.append(dict(inchikey=r.inchikey, smiles=r.smiles, label=r.label, used_in=r.used_in, family=r.family, n_carbon=r.n_carbon,
                                 tm_K=float(tm_hat[i]) if r.tm_tier == "needs_estimate" else np.nan,
                                 tm_measured_K=tm_meas, dhfus_kJ_mol=float(ds_hat[i] * tm_use / 1000.0) if r.dh_tier == "needs_estimate" else np.nan,
                                 source=FE.SOURCE, oof_rmse_Tm_K=rmse_tm, oof_rmse_dHfus_kJ_mol=rmse_dh,
                                 gate_R_W3_3="passed" if gate else "failed",
                                 note="estimate; never replaces a measured value; not read by the API or designs" +
                                      ("" if gate else "; gate failed: not for use in any margin")))
    est = pd.DataFrame(est_rows, columns=["inchikey", "smiles", "label", "used_in", "family", "n_carbon", "tm_K", "tm_measured_K", "dhfus_kJ_mol",
                                          "source", "oof_rmse_Tm_K", "oof_rmse_dHfus_kJ_mol", "gate_R_W3_3", "note"])
    ML_OUT.parent.mkdir(parents=True, exist_ok=True)
    est.to_csv(ML_OUT, index=False)
    gaps.to_csv(OUT / "gap_table.csv", index=False)

    # propagation constants: verified > thermoml_measured > ml_estimated
    const = {}
    for r in meas.itertuples():
        if np.isfinite(r.tm_K) and np.isfinite(r.dh_kJ):
            const[r.inchikey] = dict(tm_K=r.tm_K, dh_kJ=r.dh_kJ, est=False)
    for r in est.itertuples():
        tm = r.tm_K if np.isfinite(r.tm_K) else r.tm_measured_K
        dh = r.dhfus_kJ_mol if np.isfinite(r.dhfus_kJ_mol) else (m.loc[r.inchikey, "dh_kJ"] if r.inchikey in m.index else np.nan)
        if r.inchikey not in const and np.isfinite(tm) and np.isfinite(dh):
            const[r.inchikey] = dict(tm_K=float(tm), dh_kJ=float(dh), est=True)
    prop = propagation(const, rmse_tm, rmse_dh)
    prop.to_csv(OUT / "propagation.csv", index=False)

    readings = {
        "training": dict(n_compounds=int(len(train)), n_verified_tm=int((train.tm_source == "verified").sum()),
                         n_thermoml_tm=int((train.tm_source == "thermoml_measured").sum()), n_with_dh=int(np.isfinite(train.dh_kJ).sum()),
                         families=oof.family.value_counts().to_dict()),
        "R-W3-3": dict(n_C8_C20=int(g.n), mae_tm_lofo_K=float(g.mae_tm_lofo), mae_tm_joback_K=float(g.mae_tm_joback), threshold_K=GATE_K,
                       passed=gate, expected="fail"),
        "R-W3-4": dict(n=int(g.n_dh), mae_dh_lofo_kJ=float(g.mae_dh_lofo_kJ), mape_dh_lofo_pct=float(g.mape_dh_lofo),
                       mae_dh_joback_kJ=float(g.mae_dh_joback_kJ), rmse_dh_lofo_kJ=rmse_dh),
        "R-W3-5": dict(mae_tm_random5_C8_C20_K=float(g.mae_tm_random5), mae_tm_lofo_all_K=float(a.mae_tm_lofo)),
        "R-W3-6": dict(n_gap_compounds=int(len(gaps)), tm_tiers=gaps.tm_tier.value_counts().to_dict(), dh_tiers=gaps.dh_tier.value_counts().to_dict(),
                       needs_estimate=need[["label", "smiles", "used_in"]].fillna("").to_dict("records")),
        "R-W3-7": prop[["pool", "spec", "estimated_components", "freezing_C_tiered", "max_abs_change_pm1rmse_C"]].to_dict("records"),
    }
    (OUT / "readings_fusion.json").write_text(json.dumps(readings, indent=2, default=float), encoding="utf-8")
    pd.set_option("display.width", 250, "display.max_columns", 30)
    print(summ.round(2).to_string())
    print(gaps.to_string())
    print(est.drop(columns=["note"]).round(2).to_string())
    print(prop.round(3).to_string())
    print(json.dumps({k: readings[k] for k in ("training", "R-W3-3", "R-W3-4", "R-W3-5")}, indent=2, default=float))


if __name__ == "__main__":
    warnings.simplefilter("ignore")
    main()
