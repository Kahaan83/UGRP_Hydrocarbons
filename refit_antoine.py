#!/usr/bin/env python
"""Refit Antoine constants from union of WebBook, ThermoML, and T01 Cox data.

    python build_unified.py      # needs thermoml_hc.csv  -> data/unified/measurements.csv
    python refit_antoine.py      # -> data/processed/antoine_refit.csv
    python build_unified.py      # again: compounds.csv now carries the refit constants

Why: WebBook Antoine sets for C10+ alkanes start at ~95-190 C, but flash points need 40-100 C.
Units of THIS file's output:  log10( P / kPa ) = A - B / ( T/K + C )     (WebBook is bar: A_kPa = A_bar + 2)

A compound is refit only if it has >= 6 points spanning >= 30 K within 230-520 K and P <= 300 kPa.
If a WebBook set exists, it is used as a weak anchor (weight 0.3) and cross-checked via `webbook_maxdiff`.
"""
from __future__ import annotations
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import least_squares
from rdkit import Chem

UNI, PROC, RAW = Path("data/unified"), Path("data/processed"), Path("data/raw")
WINDOW = (313.15, 373.15)


def antoine_logp(T, A, B, C):
    return A - B / (np.asarray(T, float) + C)


LO = np.array([-5.0, 50.0, -250.0])
HI = np.array([20.0, 12000.0, 150.0])


def fit_antoine(T, logP, w=None):
    """Robust 3-parameter fit. Starting points are clipped into the bounds (bad data - falling or flat pressure -
    used to give a negative B and crash scipy with 'Initial guess is outside of provided bounds'), and a
    physically impossible trend raises ValueError so the caller can mark just that compound as failed."""
    T, logP = np.asarray(T, float), np.asarray(logP, float)
    w = np.ones_like(T) if w is None else np.asarray(w, float)
    ok = np.isfinite(T) & np.isfinite(logP) & np.isfinite(w)
    T, logP, w = T[ok], logP[ok], w[ok]
    if len(T) < 4 or np.ptp(T) <= 0:
        raise ValueError("not enough distinct temperatures")
    if np.corrcoef(T, logP)[0, 1] <= 0.5:
        raise ValueError("vapour pressure does not increase with temperature (bad/mixed data)")
    best = None
    for C0 in (-120.0, -80.0, -50.0, -20.0, 0.0):
        slope, icpt = np.polyfit(1.0 / (T + C0), logP, 1)
        p0 = np.clip([icpt, -slope, C0], LO + 1e-6, HI - 1e-6)
        try:
            res = least_squares(lambda p: w * (antoine_logp(T, *p) - logP), p0, bounds=(LO, HI), loss="soft_l1", f_scale=0.02)
        except ValueError:
            continue
        if best is None or res.cost < best.cost:
            best = res
    if best is None:
        raise ValueError("no starting point converged")
    return best.x


def load_t01_cox() -> dict[str, tuple[str, float, float, float, float, float]]:
    """Parse T01 Cox parameters (JPCRD 23, 1994) for C5-C20 n-alkanes.
    Returns: InChIKey -> (smiles, A0, A1, A2, To, Po)
    """
    t01_path = RAW / "ruzicka_majer" / "t01_transcription.md"
    if not t01_path.exists():
        return {}
    alkane_smiles = {
        "pentane": "CCCCC", "hexane": "CCCCCC", "heptane": "CCCCCCC",
        "octane": "CCCCCCCC", "nonane": "CCCCCCCCC", "decane": "CCCCCCCCCC",
        "undecane": "CCCCCCCCCCC", "dodecane": "CCCCCCCCCCCC", "tridecane": "CCCCCCCCCCCCC",
        "tetradecane": "CCCCCCCCCCCCCC", "pentadecane": "CCCCCCCCCCCCCCC",
        "hexadecane": "CCCCCCCCCCCCCCCC", "heptadecane": "CCCCCCCCCCCCCCCCC",
        "octadecane": "CCCCCCCCCCCCCCCCCC", "nonadecane": "CCCCCCCCCCCCCCCCCCC",
        "eicosane": "CCCCCCCCCCCCCCCCCCCC",
    }
    out = {}
    lines = t01_path.read_text(encoding="utf-8").splitlines()
    in_table = False
    for line in lines:
        if "## Table 10" in line:
            in_table = True
            continue
        if in_table:
            if line.startswith("## "):
                break
            parts = [p.strip() for p in line.split(",")]
            if len(parts) == 6 and parts[0] in alkane_smiles:
                name, a0, a1, a2, to, po = parts
                smi = alkane_smiles[name]
                m = Chem.MolFromSmiles(smi)
                if m is not None:
                    ik = Chem.MolToInchiKey(m)
                    out[ik] = (smi, float(a0), float(a1), float(a2), float(to), float(po))
    return out


def webbook_anchors() -> dict[str, tuple[float, float, float, float, float]]:
    """InChIKey -> (A_kPa, B, C, T_lo, T_hi) from antoine_best.csv or antoine_all_sets.csv (in kPa)."""
    out = {}
    p_best = PROC / "antoine_best.csv"
    if p_best.exists():
        b = pd.read_csv(p_best)
        if "smiles" not in b and Path("compounds.csv").exists():
            b = b.merge(pd.read_csv("compounds.csv")[["name", "smiles"]], on="name", how="left")
        for r in b.dropna(subset=["smiles"]).itertuples():
            m = Chem.MolFromSmiles(r.smiles)
            if m is not None:
                out[Chem.MolToInchiKey(m)] = (r.A + 2.0, r.B, r.C, r.T_lo, r.T_hi)
    p_all = PROC / "antoine_all_sets.csv"
    if p_all.exists():
        b_all = pd.read_csv(p_all)
        if "smiles" in b_all.columns:
            for r in b_all.dropna(subset=["smiles"]).itertuples():
                m = Chem.MolFromSmiles(r.smiles)
                if m is not None:
                    ik = Chem.MolToInchiKey(m)
                    if ik not in out:
                        out[ik] = (r.A + 2.0, r.B, r.C, r.T_lo, r.T_hi)
    return out


def main():
    meas_p = UNI / "measurements.csv"
    if not meas_p.exists():
        raise SystemExit("run build_unified.py first (it needs data/processed/thermoml_hc.csv)")
    m = pd.read_csv(meas_p, low_memory=False)
    vp = m[(m.property == "vapor_pressure_kPa") & (m.n_comp == 1)].copy()
    if vp.empty:
        raise SystemExit("no pure-component vapour-pressure rows in measurements.csv - did thermoml_extract.py produce any?")
    vp["key"] = vp.comp_inchikeys.map(lambda s: json.loads(s)[0])

    comp_p = UNI / "compounds.csv"
    comp = pd.read_csv(comp_p).set_index("inchikey") if comp_p.exists() else pd.DataFrame()

    anchors = webbook_anchors()
    t01_cox = load_t01_cox()

    all_keys = sorted(set(vp["key"]).union(t01_cox.keys()).union(anchors.keys()))
    vp_groups = {k: g for k, g in vp.groupby("key")}
    out = []

    for key in all_keys:
        name = comp.loc[key, "smiles"] if key in comp.index and "smiles" in comp.columns else key
        if name == key and key in t01_cox:
            name = t01_cox[key][0]

        # 1. ThermoML measured points
        T_meas, P_meas = np.array([]), np.array([])
        if key in vp_groups:
            g = vp_groups[key]
            g = g[(g.T_K.between(230, 520)) & (g.value > 0) & (g.value <= 300)]
            if not g.empty:
                g = g.assign(Tr=g.T_K.round(1)).groupby("Tr").agg(P=("value", "median"), T=("T_K", "mean")).reset_index(drop=True)
                T_meas = g["T"].to_numpy(dtype=float).flatten()
                P_meas = g["P"].to_numpy(dtype=float).flatten()

        # 2. T01 Cox equation grid points for C5-C20 n-alkanes across 313.15 to 373.15 K
        T_t01, P_t01 = np.array([]), np.array([])
        if key in t01_cox:
            smi, a0, a1, a2, to, po = t01_cox[key]
            name = smi
            T_t01 = np.linspace(WINDOW[0], WINDOW[1], 13)
            # Ruzicka & Majer 1994 Eq. (19): ln(P/Po) = (1 - To/T) * exp(A0 + A1*T + A2*T^2)
            ln_p_po = (1.0 - to / T_t01) * np.exp(a0 + a1 * T_t01 + a2 * T_t01**2)
            P_t01 = po * np.exp(ln_p_po)

        T_data = np.r_[T_meas, T_t01]
        P_data = np.r_[P_meas, P_t01]
        w_data = np.r_[np.ones_like(T_meas), np.ones_like(T_t01)]

        wb = anchors.get(key)

        # 3. If no ThermoML and no T01, but WebBook exists (e.g. n-propylbenzene, n-butylbenzene, tetralin)
        if len(T_data) == 0:
            if wb:
                T_data = np.linspace(wb[3], wb[4], 15)
                P_data = 10 ** antoine_logp(T_data, *wb[:3])
                w_data = np.ones_like(T_data)
            else:
                continue

        row = dict(inchikey=key, smiles=name, n_points=len(T_data))
        if len(T_data) < 6 or (T_data.max() - T_data.min()) < 30:
            out.append({**row, "status": "skipped: <6 points or span <30 K"})
            continue

        # Anchoring with WebBook if both exist
        Ts, ls, ws = T_data, np.log10(P_data), w_data
        if wb and (len(T_meas) > 0 or len(T_t01) > 0):
            Ta = np.linspace(wb[3], wb[4], 8)
            Ts = np.r_[T_data, Ta]
            ls = np.r_[np.log10(P_data), antoine_logp(Ta, *wb[:3])]
            ws = np.r_[w_data, np.full(8, 0.3)]

        try:
            A, B, C = fit_antoine(Ts, ls, ws)
        except ValueError as e:
            out.append({**row, "status": f"failed: {e}"})
            continue

        resid = antoine_logp(T_data, A, B, C) - np.log10(P_data)
        rms = float(np.sqrt(np.mean(resid ** 2)))
        max_abs = float(np.max(np.abs(resid)))
        maxdiff = np.nan
        if wb:
            Tc = np.linspace(wb[3], wb[4], 30)
            maxdiff = float(np.max(np.abs(antoine_logp(Tc, A, B, C) - antoine_logp(Tc, *wb[:3]))))

        covers = bool(T_data.min() <= WINDOW[0] and T_data.max() >= WINDOW[1])
        status = "ok" if rms < 0.03 else "check: rms > 0.03 log10 units (~7 %)"

        out.append({
            **row,
            "A": A, "B": B, "C": C,
            "T_min": float(T_data.min()), "T_max": float(T_data.max()),
            "rms_log10P": rms,
            "max_abs_log10P": max_abs,
            "covers_window": covers,
            "anchored": bool(wb),
            "webbook_maxdiff": maxdiff,
            "status": status,
        })

    df = pd.DataFrame(out)
    PROC.mkdir(parents=True, exist_ok=True)
    df.to_csv(PROC / "antoine_refit.csv", index=False)
    ok = df[df.status == "ok"]
    print(f"compounds with vapour-pressure data: {len(df)}   refit ok: {len(ok)}   "
          f"covering 40-100 C: {int(ok.get('covers_window', pd.Series(dtype=bool)).sum())}")
    if len(df) - len(ok):
        print(df[df.status != "ok"][["smiles", "n_points", "status"]].head(15).to_string(index=False))
    print(f"wrote {PROC/'antoine_refit.csv'}")


if __name__ == "__main__":
    main()
