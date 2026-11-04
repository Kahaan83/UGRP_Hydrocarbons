#!/usr/bin/env python
"""T4: score the composition->property pipeline on REAL fuels (CRC 647), never trained on.

    python scripts/t4_real_fuels.py

Real fuels are L0 (excluded by R1 from every training table), so they are a legitimate
out-of-distribution test set. Surrogates: 4 representative molecules, one per measured
class (saturates split 0.6/0.25/0.15 n/iso/cyclo by default, olefins -> iso-paraffin,
aromatics incl. naphthalenes -> ethylbenzene; vol% -> mass% via equal-density approx).
Predictions: physics-only paths (pure GBM component props -> id-volume density /
Refutas viscosity; Le Chatelier flash with compounds.csv Antoine+LFL).
Outputs: data/runs/t4_real_fuel_scores.csv, reports/T4_REAL_FUELS.md.
"""
from __future__ import annotations

import sys
from pathlib import Path

import joblib
import numpy as np
import pandas as pd

BASE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(BASE / "src"))
from fuelblend.features import pure_features as pf  # noqa: E402
from fuelblend.eval import metrics as evm  # noqa: E402

DEDUP = BASE / "data/unified/measurements_dedup.csv"
COMPOUNDS = BASE / "data/unified/compounds.csv"
MODELS = BASE / "runs/pure_models"
SCORES = BASE / "data/runs/t4_real_fuel_scores.csv"
REPORT = BASE / "reports/T4_REAL_FUELS.md"

LOG: list[str] = []
SEEDS = [11, 23, 37, 41, 59]


def log(msg) -> None:
    print(msg, flush=True)
    LOG.append(str(msg))


# ----------------------------------------------------------------- surrogate setup
REPS = {  # class -> (smiles)
    "n_paraffin": "CCCCCCCCCCC",      # n-undecane, Antoine 313.2-458.6 K
    "iso_paraffin": "CC(C)CC(C)(C)C", # TMP, 283.6-413.2 K
    "cycloalkane": "CC1CCCC1",        # methylcyclohexane, 300.4-345.1 K
    "aromatic": "CCc1ccccc1",         # ethylbenzene, 305.5-447.2 K
}
SAT_DIV = {"n_paraffin": 0.60, "iso_paraffin": 0.25, "cycloalkane": 0.15}  # default split of saturates

_dd = pd.read_csv(DEDUP, low_memory=False)
_comp = pd.read_csv(COMPOUNDS, low_memory=False).drop_duplicates("smiles").set_index("smiles")
MODELS_LOADED: dict = {}


def pure_pred(prop: str, smiles: str, T_K: float):
    """Ensemble mean+std of the trained pure GBM in ORIGINAL units."""
    key = (prop, smiles)
    if key not in MODELS_LOADED:
        cls = _comp.loc[smiles, "cls"] if smiles in _comp.index else "other"
        fer = pf.pure_feature_frame([dict(inchikey=str(_comp.loc[smiles, "inchikey"]), smiles=smiles, T_K=T_K, cls=cls)])
        preds = [joblib.load(MODELS / f"{prop}__seed_{s}.joblib").predict(fer)[0] for s in SEEDS]
        MODELS_LOADED[key] = (float(np.mean(preds)), float(np.std(preds)), fer)
    mean, std, fer = MODELS_LOADED[key]
    return mean, std, fer


def pure_value(prop: str, smiles: str, T_K: float):
    mean, std, fer = pure_pred(prop, smiles, T_K)
    if key_out := None:
        pass
    if prop == "viscosity_dyn_mPas":
        return float(np.exp(mean)), fer
    return float(mean), fer


def chat_flash(w_mass: list[float], smiles_list: list[str]) -> float:
    """Le Chatelier + Raoult with compounds.csv Antoine (kPa) + LFL. Bisection."""
    comp = _comp
    xs = np.asarray(w_mass, float)
    xs = xs / xs.sum()
    mw = np.asarray([float(comp.loc[s, "mw"]) for s in smiles_list])
    x = (xs / mw) / (xs / mw).sum()
    A = np.asarray([float(comp.loc[s, "antoine_A"]) for s in smiles_list])
    B = np.asarray([float(comp.loc[s, "antoine_B"]) for s in smiles_list])
    C = np.asarray([float(comp.loc[s, "antoine_C"]) for s in smiles_list])
    lfl = np.asarray([float(comp.loc[s, "lfl_vol_pct"]) / 100.0 for s in smiles_list])

    def H(T):
        Psat = 10 ** (A - B / (T + C))
        return float(np.sum(x * Psat / 101.325 / lfl))

    lo, hi = 200.0, 500.0
    if H(lo) >= 1:
        return lo - 273.15
    for _ in range(80):
        mid = (lo + hi) / 2
        if H(mid) < 1:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2 - 273.15


def surrogate_w(vol_fracs: dict) -> tuple[list[str], list[float]]:
    """Measured class vol% -> mass fractions on the 4 representatives (equal-density approx)."""
    sats = vol_fracs.get("saturates_vol_pct", 0.0)
    aro = vol_fracs.get("aromatics_vol_pct", 0.0)
    ole = vol_fracs.get("olefins_vol_pct", 0.0)
    v = {"n_paraffin": sats * SAT_DIV["n_paraffin"],
         "iso_paraffin": sats * SAT_DIV["iso_paraffin"] + ole,
         "cycloalkane": sats * SAT_DIV["cycloalkane"],
         "aromatic": aro}
    smiles = [REPS[c] for c, val in v.items() if val > 0]
    w = [val for val in v.values() if val > 0]
    if not smiles:
        return [], []
    tot = sum(w)
    return smiles, [x / tot for x in w]


def predict_fuel(vol_fracs: dict, T_K: float) -> dict:
    sm, w = surrogate_w(vol_fracs)
    if not sm:
        return {}
    out = {}
    # viscosity: Refutas over GBM component nu; density: ideal volume
    vbn = 0.0
    rho_ids = []
    for s, wi in zip(sm, w):
        eta, fer = pure_value("viscosity_dyn_mPas", s, T_K)
        rho, _ = pure_value("density_kg_m3", s, T_K)
        nu = eta / (rho / 1000.0)
        vbn += wi * (14.534 * np.log(np.log(nu + 0.8)) + 10.975)
        rho_ids.append((wi, rho))
    nu_mix = np.exp(np.exp((vbn - 10.975) / 14.534)) - 0.8
    rho_mix = 1.0 / sum(wi / ri for wi, ri in rho_ids)
    out["viscosity_kin_cSt"] = nu_mix * rho_mix / 1000.0
    out["density_kg_m3"] = rho_mix
    out["flash_point_C"] = chat_flash(w, sm)
    out["_smiles"] = ','.join(sm)
    out["_w"] = ','.join(f"{x:.3f}" for x in w)
    return out


def main() -> None:
    log("T4 real-fuel scoring (CRC 647 surrogates, physics-only paths)")
    log(f"representatives: {REPS}")
    log("")
    rf = _dd[_dd["kind"] == "real_fuel"]
    crc = rf[rf["sources"].astype(str).str.contains("crc647", na=False)]

    # measured class composition per fuel (median across rows)
    cls_props = ["saturates_vol_pct", "olefins_vol_pct", "aromatics_vol_pct", "naphthalenes_vol_pct"]
    comp_rows = []
    for sid, g in crc[crc["property"].isin(cls_props)].groupby("system_id"):
        comp_rows.append(dict(system_id=sid, **{p: float(g.loc[g["property"] == p, "value"].median()) for p in cls_props}))
    comp_df = pd.DataFrame(comp_rows)
    log(f"fuels with measured class composition: {len(comp_df)}")

    # sanity self-checks before trusting the pipeline
    dec_eta, _ = pure_value("viscosity_dyn_mPas", REPS["n_paraffin"], 293.15)
    dec_rho, _ = pure_value("density_kg_m3", REPS["n_paraffin"], 293.15)
    assert 600 < dec_rho < 800, f"GBM density sanity {dec_rho}"
    assert 0.5 < dec_eta < 5, f"GBM eta sanity {dec_eta}"
    ch_flash = chat_flash([1.0], [REPS["cycloalkane"]])
    assert -35 < ch_flash < 10, f"chat cyclohexane sanity {ch_flash}"
    log(f"sanity: GBM n-decane@293.15K eta={dec_eta:.3f} rho={dec_rho:.1f}; chat cyclohexane flash={ch_flash:.1f} C (lit -18 C)")

    # records: measured fuel rows for the three scored properties
    scored = crc[crc["property"].isin(["density_kg_m3", "flash_point_C", "viscosity_kin_cSt"])].copy()
    out_rows = []
    for r in scored.itertuples():
        c = comp_df[comp_df["system_id"] == r.system_id]
        if c.empty:
            continue
        T = float(r.T_K) if not pd.isna(r.T_K) else float("nan")
        # CRC density rows carry T=15.000 -> that is deg C in a K column (data bug); fix for scoring only
        t_bug = bool(r.property == "density_kg_m3" and T == 15.0)
        T_eff = T + 273.15 if t_bug else T
        p = predict_fuel({k: float(c.iloc[0][k]) for k in cls_props}, T_eff)
        if not p or r.property not in p:
            continue
        out_rows.append(dict(
            fuel=str(r.fuel_id) if hasattr(r, "fuel_id") else "", system_id=r.system_id, property=r.property,
            T_K=T, T_bug_celsius=t_bug, measured=float(r.value), predicted=p[r.property],
            smiles=p["_smiles"], w=p["_w"],
        ))
        if r.property == "flash_point_C":
            # also evaluate the surrogate Le Chatelier for the *co-surrogate* path & a GBM mix fallback
            gbm_mix = sum(wi * pure_value("flash_point_C", s, T_eff)[0] for s, wi in zip(p["_smiles"].split(','), [float(x) for x in p["_w"].split(',')]))
            out_rows[-1]["predicted_gbm_mix"] = gbm_mix

    od = pd.DataFrame(out_rows)
    SCORES.parent.mkdir(parents=True, exist_ok=True)
    od.to_csv(SCORES, index=False)
    log(f"\nwrote {SCORES} rows={len(od)}")

    lines = ["# T4 — Real-fuel scoring (CRC 647, surrogates, physics-only paths)", "",
             "## surrogate scheme"] + LOG
    REPORT.write_text("\n".join(lines) + "\n", encoding="utf-8")

    # stats
    for prop, mode in [("density_kg_m3", "generic"), ("viscosity_kin_cSt", "generic"), ("flash_point_C", "flash")]:
        g = od[od["property"] == prop].dropna(subset=["measured", "predicted"])
        if not len(g):
            log(f"\n{prop}: no rows")
            continue
        cl = np.array(g["system_id"])
        st = evm.report(g["measured"].values, g["predicted"].values, cl, mode=mode)
        mape = float(np.mean(np.abs(g["measured"] - g["predicted"]) / np.clip(np.abs(g["measured"]), 1e-9, None)) * 100)
        within10 = float(np.mean(np.abs(g["measured"] - g["predicted"]) <= 0.10 * np.abs(g["measured"])))
        log(f"\n== {prop} (n={len(g)} fuels={g['system_id'].nunique()}) "
            f"MAE={st['mae']:.3f} RMSE={st['rmse']:.3f} bias={st.get('bias', float('nan')):.3f} MAPE={mape:.1f}% within10%={within10:.1%}")
        if prop == "flash_point_C" and "predicted_gbm_mix" in g.columns:
            st2 = evm.report(g["measured"].values, g["predicted_gbm_mix"].values, cl, mode=mode)
            log(f"   vs GBM-mix fallback path: MAE={st2['mae']:.3f} RMSE={st2['rmse']:.3f} "
                f"overpred5={st2.get('overpred_rate_5', float('nan')):.1%} mean_pos_err={st2.get('mean_positive_error', float('nan')):.3f}")

    log("\nDONE")


if __name__ == "__main__":
    main()
