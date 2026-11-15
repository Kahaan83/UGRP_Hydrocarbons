"""Transcribed primary-paper calorimetry data -> data/raw/calorimetry_primary/*.csv  (data only; nothing adopted).

Sources (all read from page images, not from the OCR text layer):
  F54  Finke, Gross, Waddington, Huffman, J. Am. Chem. Soc. 76, 333-341 (1954)   Tables I, IV, V
  DH46 Douslin, Huffman, J. Am. Chem. Soc. 68, 1704-1708 (1946)                  Tables I, IV, V, VI
  PS43 Pitzer, Scott, J. Am. Chem. Soc. 65, 803-829 (1943)                       Tables I, VI, VII, X + text (heats of vaporization, eqs 1-5)
  S62  Scott, Guthrie, Messerly, Todd, Berg, Hossenlopp, McCullough, J. Phys. Chem. 66, 911-914 (1962)  Tables V, VII, eq 2

Run:  .venv/Scripts/python.exe scripts/v5/base/ingest_calorimetry_primary.py [--integrate]
(--integrate appends the data to data/unified/measurements*.csv and data/processed/v5/thermoml_sle_hc.csv; idempotent.)
Unit rules: kJ = cal * 4.184 / 1000 (thermochemical calorie; reproduces the repo's existing values to the last printed digit);
mmHg -> kPa x 0.133322368.  Temperatures are kept AS PUBLISHED (0 degC = 273.16 K in all four papers); nothing is shifted.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pandas as pd
from rdkit import Chem

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "data" / "raw" / "calorimetry_primary"
VERIFIED = ROOT / "data" / "raw" / "sle_constants" / "tm_dhfus_verified.csv"
CAL = 4.184
MMHG_KPA = 0.133322368

CIT = {
    "F54": "Finke, Gross, Waddington, Huffman, J. Am. Chem. Soc. 1954, 76, 333-341",
    "DH46": "Douslin, Huffman, J. Am. Chem. Soc. 1946, 68, 1704-1708",
    "PS43": "Pitzer, Scott, J. Am. Chem. Soc. 1943, 65, 803-829",
    "S62": "Scott, Guthrie, Messerly, Todd, Berg, Hossenlopp, McCullough, J. Phys. Chem. 1962, 66, 911-914",
}

SMILES = {
    "n-octane": "CCCCCCCC", "n-nonane": "CCCCCCCCC", "n-decane": "CCCCCCCCCC", "n-undecane": "CCCCCCCCCCC",
    "n-dodecane": "CCCCCCCCCCCC", "n-tridecane": "CCCCCCCCCCCCC", "n-tetradecane": "CCCCCCCCCCCCCC",
    "n-pentadecane": "CCCCCCCCCCCCCCC", "n-hexadecane": "CCCCCCCCCCCCCCCC",
    "n-hexane": "CCCCCC", "2-methylpentane": "CC(C)CCC", "3-methylpentane": "CCC(C)CC",
    "2,3-dimethylbutane": "CC(C)C(C)C", "2,2-dimethylbutane": "CC(C)(C)CC",
    "toluene": "Cc1ccccc1", "o-xylene": "Cc1ccccc1C", "m-xylene": "Cc1cccc(C)c1", "p-xylene": "Cc1ccc(C)cc1",
}


def ik(name: str) -> str:
    return Chem.MolToInchiKey(Chem.MolFromSmiles(SMILES[name]))


def kj(cal):
    return None if cal is None else round(cal * CAL / 1000.0, 3)


# ---------------------------------------------------------------- phase-change constants
# (paper, table, name, tm_K, tm_note, impurity_mol_pct, impurity_unc, dh_cal, dh_unc, [(T_trans, dH_trans_cal, unc), ...], dh_kind, note)
PC = [
    ("F54", "I (Pure row); IV", "n-octane", 216.375, "triple point (extrapolated to pure)", 0.06, 0.03, 4957, 1, [], "fusion only", ""),
    ("F54", "I; IV; V", "n-nonane", 219.659, "triple point (extrapolated to pure)", 0.12, 0.05, 3697, None,
     [(217.2, 1501, None)], "fusion only",
     "Table IV footnote c: fusion heat obtained from enthalpy measurements of transition + fusion combined (not a separate calorimetric run); transition T 217.2 (Table V) / 217.19 (Table VI)"),
    ("F54", "I; IV", "n-decane", 243.512, "triple point (extrapolated to pure)", 0.09, 0.04, 6863, 2, [], "fusion only", ""),
    ("F54", "I; IV; V", "n-undecane", 247.592, "triple point (extrapolated to pure)", 0.02, 0.01, 5301, None,
     [(236.6, 1639, None)], "fusion only", "Table IV footnote c: obtained from transition + fusion heats combined; Table V transition 1639 cal/mol at 236.6 K"),
    ("F54", "I; IV", "n-dodecane", 263.589, "triple point (extrapolated to pure)", 0.07, 0.05, 8804, 5, [], "fusion only", ""),
    ("F54", "I; IV; V", "n-tridecane", 267.792, "triple point (extrapolated to pure)", 0.05, 0.04, 6812, 1,
     [(255.0, 1831, 2)], "fusion only", "Table V: expts I 1832, II 1829, mean 1831 +/- 2"),
    ("F54", "I; IV", "n-tetradecane", 279.027, "triple point (extrapolated to pure)", 0.07, 0.04, 10772, 1, [], "fusion only", ""),
    ("F54", "I; IV; V", "n-pentadecane", 283.108, "triple point (extrapolated to pure)", 0.05, 0.04, 8268, 4,
     [(270.9, 2191, 5)], "fusion only", "Table V: expts I 2186, II 2195, mean 2191 +/- 5"),
    ("F54", "I; IV", "n-hexadecane", 291.340, "triple point (extrapolated to pure)", 0.12, 0.06, 12753, 3, [], "fusion only", ""),

    ("DH46", "I; VI", "n-hexane", 177.84, "triple point +/- 0.05 K", 0.0075, None, 3126.1, 2.2, [], "fusion only", "expts 3127.0 / 3127.4 / 3123.9"),
    ("DH46", "I; VI", "2-methylpentane", 119.55, "triple point +/- 0.05 K", 0.029, None, 1497.6, 1.6, [], "fusion only",
     "expts 1496.4 / 1499.2 / 1497.2. Crystal form; the paper also measured a glass"),
    ("DH46", "I; IV; V; VI", "2,3-dimethylbutane", 145.19, "triple point +/- 0.05 K (calculated; Raoult's law not obeyed)", 0.015, None, 191.42, 0.2,
     [(136.07, 1552.0, 1.0)], "fusion only",
     "two crystalline forms; transition 136.07 +/- 0.05 K, 1552.0 +/- 1.0 cal/mol (expts 1551.9 / 1552.2); paper: T_tp and impurity 'more or less arbitrary'"),
    ("DH46", "I; IV; V; VI", "2,2-dimethylbutane", 174.28, "triple point +/- 0.05 K (calculated; Raoult's law not obeyed)", 0.015, None, 138.40, 0.2,
     [(126.81, 1292.7, 1.0), (140.79, 68.20, 0.12)], "fusion only",
     "three crystalline forms; lower transition 126.81 K 1292.7 +/- 1.0 cal/mol; upper transition 140.79 K 68.20 +/- 0.12 cal/mol; fusion expts 138.24 / 138.57"),
    ("DH46", "Table I text", "3-methylpentane", None, "not crystallisable (glass); no triple point", None, None, None, None, [], "none",
     "paper: impossible to crystallise 3-methylpentane, purity could not be checked by freezing point; Tm and dHfus intentionally blank"),

    ("PS43", "I", "o-xylene", 247.82, "true melting point +/- 0.1 K", None, None, 3250, 10, [], "fusion only",
     "Table VI expts 3250/3245/3255 (heating times 19/24/27 min); other labs: PHD 3116, Rossini 3330, Nakatsuchi 3100"),
    ("PS43", "I; VI", "m-xylene", 225.27, "true melting point +/- 0.15 K (impure sample melted at 224.38 K)", 2.34, None, 2765, 20, [], "fusion only (corrected for p-xylene impurity)",
     "sample: total impurity 2.34 mol% (0.22% non-para); apparent dHfus of the impure sample 2784 +/- 15 cal/mol; 2765 is the pure-m-xylene value after correcting for p-xylene impurity using a DCp assumption (p-xylene apparent Cp = 1.004 x liquid m-xylene); a constant DCp of 5.4 would change it by 2 cal/mol"),
    ("PS43", "I; VI", "p-xylene", 286.39, "true melting point +/- 0.05 K", None, None, 4090, 20, [], "fusion only", "Table VI expts 4089/4085/4097"),

    ("S62", "V; text", "toluene", 178.15, "triple point 178.15 +/- 0.05 K (calculated pure 178.1504)", 0.0012, 0.001, 1586, 1, [], "fusion only",
     "impurity N2* = 0.0000012 +/- 0.00001 mole fraction (printed 0.00001_2); cryoscopic A = 0.02515 deg^-1, B = 0.00216 deg^-1; sample purity 99.999 mol%; DHm from three determinations. Pitzer & Scott 1943 did not measure toluene fusion, so the repo's toluene row has one primary source (S62), not two"),
]


def phase_change_rows():
    rows = []
    for paper, table, name, tm, tm_note, imp, imp_u, dh, dh_u, trans, kind, note in PC:
        r = {"inchikey": ik(name), "name": name, "smiles": SMILES[name], "paper": paper, "citation": CIT[paper], "table_ref": table,
             "tm_K": tm, "tm_note": tm_note, "impurity_mol_pct": imp, "impurity_unc_mol_pct": imp_u,
             "dhfus_cal_mol": dh, "dhfus_unc_cal_mol": dh_u, "dhfus_kJ_mol": kj(dh), "dhfus_kind": kind}
        for i in (1, 2):
            if len(trans) >= i:
                t, d, u = trans[i - 1]
                r.update({f"trans{i}_T_K": t, f"trans{i}_dH_cal_mol": d, f"trans{i}_dH_unc_cal_mol": u, f"trans{i}_dH_kJ_mol": kj(d)})
            else:
                r.update({f"trans{i}_T_K": None, f"trans{i}_dH_cal_mol": None, f"trans{i}_dH_unc_cal_mol": None, f"trans{i}_dH_kJ_mol": None})
        tot = None if dh is None else round((dh + sum(d for _, d, _ in trans)) * CAL / 1000, 3)
        r["dhfus_total_incl_trans_kJ_mol"] = tot if trans else None
        r["selected"] = True
        r["note"] = note
        rows.append(r)
    return pd.DataFrame(rows)


# ---------------------------------------------------------------- vapour pressure (PS43 Table VII), T in degC as published
VP = {  # name: [(T_C, p_obs_mmHg, p_calc_mmHg, phase)]
    "toluene": [(0.0, 6.75, 6.75), (12.5, 14.40, 14.38), (25.0, 28.40, 28.45), (37.5, 52.70, 52.65), (50.0, 92.00, 91.98)],
    "o-xylene": [(0.0, 1.30, 1.27), (12.5, 3.00, 3.04), (25.0, 6.60, 6.64), (37.5, 13.50, 13.45), (50.0, 25.50, 25.48), (60.0, 40.60, 40.63)],
    "m-xylene": [(0.0, 1.75, 1.71), (12.5, 3.90, 3.94), (25.0, 8.35, 8.39), (37.5, 16.65, 16.66), (50.0, 31.20, 31.17), (60.0, 49.50, 49.49)],
    "p-xylene": [(0.0, 1.35, 1.29, "solid"), (12.5, 4.05, 4.10, "solid"), (25.0, 8.90, 8.87), (37.5, 17.50, 17.55), (50.0, 32.60, 32.61), (60.0, 51.35, 51.32)],
}
# log10 p(mmHg) = a/T + b*log10 T + c ; T in K (= t_C + 273.16 in the paper); b=0 for the solid p-xylene equation
VP_EQ = [
    ("toluene", "(1)", "liquid", 0, 50, -2866.53, -6.7, 27.6470),
    ("o-xylene", "(2)", "liquid", 0, 60, -3327.16, -8.0, 31.7771),
    ("m-xylene", "(3)", "liquid", 0, 60, -2871.66, -5.0, 22.9270),
    ("p-xylene", "(4)", "liquid", 13.23, 60, -3080.31, -6.7, 27.8581),
    ("p-xylene", "(5)", "solid", 0, 13.23, -3141.33, 0.0, 11.6092),
]


def vp_rows():
    rows = []
    for name, pts in VP.items():
        for p in pts:
            t, po, pc = p[:3]
            phase = p[3] if len(p) > 3 else "liquid"
            rows.append({"inchikey": ik(name), "name": name, "paper": "PS43", "table_ref": "VII", "T_C_published": t,
                         "T_K_published_basis": round(t + 273.16, 2), "phase": phase,
                         "p_obs_mmHg": po, "p_calc_mmHg": pc, "p_obs_kPa": round(po * MMHG_KPA, 4),
                         "p_calc_kPa": round(pc * MMHG_KPA, 4),
                         "note": "manometer readings reduced to 0 degC and standard gravity; 0.05 mm cathetometer sensitivity"})
    return pd.DataFrame(rows)


def vp_eq_rows():
    return pd.DataFrame([{"inchikey": ik(n), "name": n, "paper": "PS43", "equation": eq, "phase": ph, "T_C_lo": lo, "T_C_hi": hi,
                          "form": "log10(p/mmHg) = a/T + b*log10(T) + c, T in K (t_C + 273.16)", "a": a, "b": b, "c": c,
                          "note": "paper's own fit to Table VII; eq (5) is derived from eq (4) and the measured heat of fusion"} for n, eq, ph, lo, hi, a, b, c in VP_EQ])


# ---------------------------------------------------------------- heats of vaporization
DHVAP = [  # name, paper, table, T_K, P_atm, cal, unc, note
    ("toluene", "PS43", "text after Table VII", 298.16, None, 9115, 50, "Clapeyron on Table VII, Berthelot gas-imperfection correction 26 cal/mol"),
    ("o-xylene", "PS43", "text after Table VII", 298.16, None, 10470, 50, "gas-imperfection correction 10 cal/mol"),
    ("m-xylene", "PS43", "text after Table VII", 298.16, None, 10160, 50, "gas-imperfection correction 12 cal/mol"),
    ("p-xylene", "PS43", "text after Table VII", 298.16, None, 10110, 50, "gas-imperfection correction 11 cal/mol"),
    ("toluene", "S62", "VII", 341.27, 0.250, 8521, 4, "vapor-flow calorimetry"),
    ("toluene", "S62", "VII", 361.06, 0.500, 8254, 3, "vapor-flow calorimetry"),
    ("toluene", "S62", "VII", 383.77, 1.000, 7933, 3, "vapor-flow calorimetry (normal boiling point region)"),
    ("toluene", "S62", "VII", 410.11, 2.000, 7538, 2, "vapor-flow calorimetry"),
]


def dhvap_rows():
    rows = [{"inchikey": ik(n), "name": n, "paper": p, "citation": CIT[p], "table_ref": t, "T_K": T, "P_atm": P,
             "dHvap_cal_mol": c, "dHvap_unc_cal_mol": u, "dHvap_kJ_mol": kj(c), "selected": True, "note": note}
            for n, p, t, T, P, c, u, note in DHVAP]
    rows.append({"inchikey": ik("toluene"), "name": "toluene", "paper": "S62", "citation": CIT["S62"], "table_ref": "eq (2)", "T_K": None, "P_atm": None,
                 "dHvap_cal_mol": None, "dHvap_unc_cal_mol": None, "dHvap_kJ_mol": None, "selected": True,
                 "note": "empirical: dHv = 11637 - 4.823 T - 1.260e-2 T^2 cal/mol, T in K, valid 298-410 K (fits Osborne-Ginnings at 298.15 within +/-3 cal/mol)"})
    return pd.DataFrame(rows)


# ---------------------------------------------------------------- xylene eutectics (PS43 Table X)
X_NAMES = ["o-xylene", "m-xylene", "p-xylene"]
# (label, x_o, x_m, x_p (percent), solid phases, calc_K, exp_K, exp_unc, nakatsuchi_K, nak_unc, exp_composition_note)
EUT = [
    ("m-p binary eutectic", 0.0, 87.0, 13.0, ["m-xylene", "p-xylene"], 220.3, 220.2, 0.2, 216.2, 3.0,
     "experimental composition (this research) 12.9% para, 87.1% meta (deduced from the intersection of the melting curve and the heats of fusion)"),
    ("o-p binary eutectic", 75.7, 0.0, 24.3, ["o-xylene", "p-xylene"], 237.6, None, None, 237.7, 0.4, "no experiment in this work for this eutectic"),
    ("o-m binary eutectic", 66.6, 33.4, 0.0, ["o-xylene", "m-xylene"], 211.2, 211.4, 1.0, 210.2, 1.0, ""),
    ("o-m-p ternary eutectic", 30.5, 61.4, 8.1, ["o-xylene", "m-xylene", "p-xylene"], 208.5, 209.0, 1.0, None, None, ""),
]


def eutectic_rows():
    rows = []
    for label, xo, xm, xp, solids, calc, exp, eu, nak, nu, cnote in EUT:
        pct = {"o-xylene": xo, "m-xylene": xm, "p-xylene": xp}
        comps = [n for n in X_NAMES if pct[n] > 0]
        iks = sorted(ik(n) for n in comps)
        by_ik = {ik(n): pct[n] / 100.0 for n in comps}
        base = {"label": label, "n_comp": len(comps), "inchikeys": json.dumps(iks),
                "names": json.dumps([next(n for n in comps if ik(n) == k) for k in iks]),
                "x_liquid": json.dumps([by_ik[k] for k in iks]),
                "comp_scale": "mole fraction (ASSUMED: paper builds the diagram from Raoult's law; the table does not say 'mole')",
                "crystal_inchikeys": json.dumps(sorted(ik(n) for n in solids)), "P_kPa": 101.325,
                "paper": "PS43", "citation": CIT["PS43"], "table_ref": "X", "T_scale_note": "0 degC = 273.16 K as published"}
        rows.append({**base, "point_kind": "eutectic_calculated_ideal", "T_K": calc, "T_unc_K": None,
                     "note": "ideal-solution (Raoult) calculation with DCp of fusion = 3R, from the PS43 pure-component Tm and dHfus; NOT a measurement"})
        if exp is not None:
            rows.append({**base, "point_kind": "eutectic_measured", "T_K": exp, "T_unc_K": eu, "note": cnote})
        if nak is not None:
            rows.append({**base, "point_kind": "eutectic_literature_Nakatsuchi1929", "T_K": nak, "T_unc_K": nu,
                         "note": "older literature value quoted in Table X; secondary"})
    df = pd.DataFrame(rows)
    df.insert(0, "system_id", [hashlib.md5("|".join(sorted(json.loads(s))).encode()).hexdigest()[:12] for s in df.inchikeys])  # repo rule: md5 of sorted InChIKeys joined by '|'
    return df


# ---------------------------------------------------------------- verification vs existing registry
def verification(pc: pd.DataFrame) -> pd.DataFrame:
    v = pd.read_csv(VERIFIED).set_index("inchikey")
    rows = []
    for _, r in pc.iterrows():
        rec = {"inchikey": r.inchikey, "name": r["name"], "paper": r.paper, "in_tm_dhfus_verified": r.inchikey in v.index,
               "paper_tm_K": r.tm_K, "paper_dhfus_kJ_mol": r.dhfus_kJ_mol, "paper_trans1_dH_kJ_mol": r.trans1_dH_kJ_mol,
               "paper_trans1_T_K": r.trans1_T_K}
        if r.inchikey in v.index:
            x = v.loc[r.inchikey]
            rec.update({"verified_tm_K": x.tm_K, "verified_dhfus_kJ_mol": x.dhfus_kJ_mol, "verified_ss_dH_kJ_mol": x.dhfus_ss_transition_kJ_mol,
                        "verified_ss_T_K": x.ss_transition_T_K})
            d_tm = None if pd.isna(r.tm_K) or pd.isna(x.tm_K) else round(x.tm_K - r.tm_K, 4)
            d_dh = None if pd.isna(r.dhfus_kJ_mol) or pd.isna(x.dhfus_kJ_mol) else round(x.dhfus_kJ_mol - r.dhfus_kJ_mol, 4)
            d_ss = None if pd.isna(r.trans1_dH_kJ_mol) or pd.isna(x.dhfus_ss_transition_kJ_mol) else round(x.dhfus_ss_transition_kJ_mol - r.trans1_dH_kJ_mol, 4)
            rec.update({"d_tm_K": d_tm, "d_dhfus_kJ_mol": d_dh, "d_ss_dH_kJ_mol": d_ss})
            if r["name"] == "3-methylpentane":
                st = "paper gives no Tm/dHfus (glass); repo row comes from another source, not verifiable here"
            elif r["name"] == "2,2-dimethylbutane" and pd.isna(x.dhfus_ss_transition_kJ_mol):
                st = "Tm/dHfus agree; repo has no solid-solid transitions (paper: 126.81 K 5.409 kJ/mol and 140.79 K 0.285 kJ/mol), flagged only"
            else:
                bad = []
                if d_tm is not None and abs(d_tm) > 0.011:
                    bad.append("Tm")
                if d_dh is not None and abs(d_dh) > 0.0015:
                    bad.append("dHfus")
                if d_ss is not None and abs(d_ss) > 0.0015:
                    bad.append("ss_dH")
                if not bad:
                    st = "agrees (|dTm|<=0.011 K, |dH|<=0.0015 kJ/mol = last printed digit)"
                elif bad == ["dHfus"] and not pd.isna(r.dhfus_unc_cal_mol) and abs(d_dh) <= r.dhfus_unc_cal_mol * CAL / 1000:
                    st = (f"dHfus differs by {d_dh:+.4f} kJ/mol (> last digit) but inside the paper's stated uncertainty "
                          f"(+/-{r.dhfus_unc_cal_mol * CAL / 1000:.4f} kJ/mol); repo value is not reproduced exactly, values unchanged")
                else:
                    st = "DISCREPANCY in " + ",".join(bad)
            rec["status"] = st
        else:
            rec["status"] = "new compound: not in tm_dhfus_verified.csv" if not pd.isna(r.tm_K) else "no data in paper"
        rows.append(rec)
    return pd.DataFrame(rows)


SRC_TAG = "calorimetry_primary"


def integrate():
    """Idempotently add the new data to the pipeline files (all gitignored/generated; rerun this after regenerating them):
      * data/unified/measurements.csv + measurements_dedup.csv : 21 liquid vapour pressures (PS43 Table VII)
      * data/processed/v5/thermoml_sle_hc.csv                  : 3 measured xylene eutectics (point_kind 'eutectic'; liquidus-only consumers are unaffected)
    """
    import numpy as np
    U = ROOT / "data" / "unified"
    vp = pd.read_csv(OUT / "vapor_pressure_pitzer_scott1943.csv")
    vp = vp[vp.phase == "liquid"]                      # measurements are liquid-phase only (schema contract)
    raw = pd.read_csv(U / "measurements.csv", low_memory=False)
    raw = raw[raw.source != SRC_TAG]
    new = pd.DataFrame({
        "source": SRC_TAG, "property": "vapor_pressure_kPa", "value": vp.p_obs_kPa, "T_K": vp.T_K_published_basis, "P_kPa": np.nan,
        "kind": "pure", "n_comp": 1, "is_hc": True,
        "system_id": [hashlib.md5(k.encode()).hexdigest()[:12] for k in vp.inchikey],
        "comp_inchikeys": [json.dumps([k]) for k in vp.inchikey], "x_mole": "[1.0]", "w_mass": "[1.0]", "fuel_id": np.nan,
        "ref": "Pitzer&Scott1943 Table VII", "raw_value": vp.p_obs_mmHg, "raw_unit": "mmHg", "composition_level": "L3",
        "comp_coverage": 1.0, "comp_basis": "mole"})
    new = new.reindex(columns=raw.columns)
    allraw = pd.concat([raw, new], ignore_index=True)
    allraw.to_csv(U / "measurements.csv", index=False)

    dd = pd.read_csv(U / "measurements_dedup.csv", low_memory=False)
    dd = dd[~dd.sources.astype(str).eq(SRC_TAG)]
    new = new.assign(T_key=lambda d: (d.T_K / 0.5).round() * 0.5, P_key=101.0, x_key=1.0)
    keys = new[["property", "system_id", "x_key", "T_key", "P_key"]].drop_duplicates()
    fixed = []
    for _, k in keys.iterrows():
        m = ((dd.property == k.property) & (dd.system_id == k.system_id) & (dd.x_key == k.x_key) & (dd["T_key"] == k.T_key) & (dd.P_key == k.P_key))
        old = dd[m]
        grp = new[(new.system_id == k.system_id) & (new.T_key == k.T_key)]
        if len(old):  # merge with an existing condition: recompute the group statistics (existing medians/spread are not recoverable
            # per-observation from the dedup row, so use its n_obs copies of its median plus the new observation(s))
            o = old.iloc[0]
            vals = list(np.repeat(o.value, int(o.n_obs))) + list(grp.value)
            nref = int(o.n_refs) + grp.ref.nunique()
            srcs = sorted(set(str(o.sources).split("+")) | {SRC_TAG})
            base = o.to_dict()
            dd = dd[~m]
        else:
            vals, nref, srcs = list(grp.value), grp.ref.nunique(), [SRC_TAG]
            g0 = grp.iloc[0]
            base = {c: g0.get(c, np.nan) for c in dd.columns}
            base.update(property=k.property, system_id=k.system_id, x_key=k.x_key, T_key=k.T_key, P_key=k.P_key)
        v = np.array(vals, float)
        med = float(np.median(v))
        iqr = np.percentile(v, 75) - np.percentile(v, 25) if len(v) >= 3 else v.max() - v.min()
        base.update(value=med, vmin=v.min(), vmax=v.max(), q25=np.percentile(v, 25), q75=np.percentile(v, 75), n_obs=len(v), n_refs=nref,
                    n_sources=len(srcs), sources="+".join(srcs), n_outlier=int((np.abs(v - med) > 0.25 * abs(med)).sum()) if len(v) >= 3 else 0,
                    conflict=bool(iqr > 0.05 * abs(med)))
        fixed.append(base)
    dd = pd.concat([dd, pd.DataFrame(fixed)[dd.columns]], ignore_index=True)
    dd.to_csv(U / "measurements_dedup.csv", index=False)

    sle_p = ROOT / "data" / "processed" / "v5" / "thermoml_sle_hc.csv"
    sle = pd.read_csv(sle_p, low_memory=False)
    sle = sle[~sle.file.astype(str).str.startswith(SRC_TAG)]
    e = eutectic_rows()
    e = e[e.point_kind == "eutectic_measured"]
    rows = []
    for _, r in e.iterrows():
        iks = json.loads(r.inchikeys)
        nm = json.loads(r.names)
        smi = [SMILES[n] for n in nm]
        rows.append({"file": f"{SRC_TAG}:PS43_TableX", "doi": "", "n_comp": r.n_comp,
                     "inchis": json.dumps([Chem.MolToInchi(Chem.MolFromSmiles(s)) for s in smi]), "names": r.names, "smiles": json.dumps(smi),
                     "x_liquid": r.x_liquid, "T_K": r["T_K"], "P_kPa": np.nan, "solid_phase_desc": "eutectic (Crystal: " + "; ".join(nm) + ")",
                     "point_kind": "eutectic", "is_hc": True, "system_id": r.system_id, "inchikeys": r.inchikeys,
                     "crystal_inchikeys": r.crystal_inchikeys, "comp_scale": "mole fraction", "prop_name": "Eutectic temperature, K", "prop_phase": "Crystal"})
    pd.concat([sle, pd.DataFrame(rows)[sle.columns]], ignore_index=True).to_csv(sle_p, index=False)
    print("integrated:", len(vp), "vapour-pressure obs,", len(rows), "eutectic rows")


def main():
    import sys
    OUT.mkdir(parents=True, exist_ok=True)
    pc = phase_change_rows()
    pc.to_csv(OUT / "phase_change_constants.csv", index=False)
    vp_rows().to_csv(OUT / "vapor_pressure_pitzer_scott1943.csv", index=False)
    vp_eq_rows().to_csv(OUT / "vapor_pressure_equations_pitzer_scott1943.csv", index=False)
    dhvap_rows().to_csv(OUT / "dhvap.csv", index=False)
    eutectic_rows().to_csv(OUT / "xylene_eutectics_pitzer_scott1943.csv", index=False)
    ver = verification(pc)
    ver.to_csv(OUT / "verification_vs_tm_dhfus_verified.csv", index=False)
    print(ver[["name", "paper", "in_tm_dhfus_verified", "d_tm_K", "d_dhfus_kJ_mol", "d_ss_dH_kJ_mol", "status"]].to_string())
    if "--integrate" in sys.argv:
        integrate()


if __name__ == "__main__":
    main()
