"""v5 S2 Question A: mesitylene vapour pressure at the flash temperature (readings R-A1..R-A3 in results/v5/S2/README.md).

Reads only; writes only under results/v5/S2/. Nothing is fitted and no v4 output is overwritten.

Outputs
  mesitylene_psat_compare.csv      Psat(T) of the stored WebBook set, the registry ("anchored") set, and the `thermo` library
                                   correlation (a REFERENCE, not a measurement), plus the measured point(s)
  mesitylene_measured_psat.csv     every mesitylene pure vapour-pressure row found in the ThermoML-derived tables, kept/dropped + reason
  thermoml_bad_vp_rows_summary.csv rows in v2/thermoml_clean_hc.csv absent from the fixed v5 re-extraction (the known-bad rows)
  mesitylene_refit_decision.csv    R-A2 decision (count of usable measured points <= 345 K)
  mesitylene_flash_rerun.csv       the 126 lab records re-scored with each candidate mesitylene input (ideal, UNIFAC)
  mesitylene_flash_sensitivity.csv the v4 flash_input_sensitivity subsets, for each candidate input
  mesitylene_readings.csv          R-A1 and R-A3 values
  mesitylene_crc647_usage.csv      share of mesitylene in the r3 CRC 647 reconstructions (context only; CRC flash NOT rerun)
  mesitylene_best_input.json       the input chosen as "best available" (read by s2_t2_flash.py)
"""
from __future__ import annotations

import json
import math
import re
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from fuelblend.baselines.flash_unifac import UnifacFallbackWarning, _default_registry, flash_point_unifac_detail  # noqa: E402

OUT = ROOT / "results" / "v5" / "S2"
MES_SMILES = "Cc1cc(C)cc(C)c1"
MES_IK = "AUHZEENZYGFFBQ-UHFFFAOYSA-N"
MES_INCHI = "InChI=1S/C9H12/c1-7-4-8(2)6-9(3)5-7/h4-6H,1-3H3"
CAS = "108-67-8"
SRC = {"thermoml_hc.csv (old parser)": ROOT / "data/processed/thermoml_hc.csv",
       "v2/thermoml_clean_hc.csv": ROOT / "data/processed/v2/thermoml_clean_hc.csv",
       "v5 re-extract thermoml_hc_v5.csv": ROOT / "results/v5/thermoml_reextract/thermoml_hc_v5.csv"}
T_LOW_FOR_REFIT_K = 345.0
MIN_POINTS_REFIT = 6


def psat(T, abc):
    A, B, C = abc
    return 10.0 ** (A - B / (np.asarray(T, float) + C))


def margin90(df: pd.DataFrame, col: str, alpha: float = 0.1):
    """Same as scripts/v4/flash_input_sensitivity.py: split-conformal order statistic over per-system max error."""
    g = df.groupby("system_id")[col].max().sort_values().to_numpy()
    n = len(g)
    k = math.ceil((n + 1) * (1 - alpha))
    return (float("inf") if k > n else float(g[k - 1])), n


def stored_sets():
    src = pd.read_csv(ROOT / "data/antoine_sources.csv")
    r = src[src.inchikey == MES_IK].iloc[0]
    m = re.search(r"WebBook parameters \(.*?A=([\d.]+), B=([\d.]+), C=(-?[\d.]+)\)", r["flags/notes"])
    wb = tuple(float(m.group(i)) for i in (1, 2, 3))
    reg = _default_registry()
    can = Chem.MolToSmiles(Chem.MolFromSmiles(MES_SMILES))
    rr = reg.loc[can]
    anch = (float(rr.antoine_A), float(rr.antoine_B), float(rr.antoine_C))
    return dict(webbook=wb, anchored=anch, lo=float(r.range_T_lo_K), hi=float(r.range_T_hi_K),
                anchored_source=str(rr.antoine_source), can=can)


def measured_rows() -> pd.DataFrame:
    rows = []
    bad_summary = []
    v5 = pd.read_csv(SRC["v5 re-extract thermoml_hc_v5.csv"], low_memory=False)
    key = ["file", "prop", "value", "T_K", "inchis", "comp_values", "n_comp"]
    v5k = set(map(tuple, v5[key].astype(str).values))
    for name, p in SRC.items():
        d = pd.read_csv(p, low_memory=False)
        vp = d[d.prop.str.contains("Vapor", na=False)]
        in_v5 = np.array([tuple(r) in v5k for r in vp[key].astype(str).values])
        n_bad = int((~in_v5).sum())
        bad_summary.append(dict(table=name, n_rows=len(d), n_vp_rows=len(vp), n_vp_rows_absent_from_v5_reextract=n_bad,
                                n_all_rows_absent_from_v5_reextract=int(sum(tuple(r) not in v5k for r in d[key].astype(str).values))))
        sel = vp[vp.inchis.str.contains(MES_INCHI, regex=False, na=False)].copy()
        sel["in_v5_reextract"] = [tuple(r) in v5k for r in sel[key].astype(str).values]
        sel["table"] = name
        rows.append(sel)
    um = pd.read_csv(ROOT / "data/unified/measurements.csv", low_memory=False)
    um = um[(um.comp_inchikeys == f'["{MES_IK}"]') & (um.property == "vapor_pressure_kPa")].copy()
    um = um.rename(columns={"ref": "file"})
    um["table"] = "unified/measurements.csv"
    um["prop_phase"] = "Liquid (loader keeps liquid only)"
    um["n_comp"] = 1
    um["in_v5_reextract"] = [any((v5.file == r.file) & (v5.T_K == r.T_K) & (v5.value == r.value)) for r in um.itertuples()]
    out = pd.concat(rows + [um], ignore_index=True)[["table", "file", "prop", "prop_phase", "value", "T_K", "P_kPa", "n_comp", "in_v5_reextract"]]
    reasons = []
    for r in out.itertuples():
        why = []
        if r.n_comp != 1:
            why.append("not pure")
        if not str(r.prop_phase).startswith("Liquid"):
            why.append(f"phase={r.prop_phase}")
        if not (pd.notna(r.value) and r.value > 0):
            why.append("non-positive value")
        if pd.isna(r.T_K):
            why.append("no temperature")
        if pd.notna(r.P_kPa) and r.P_kPa > 150:
            why.append("pressure flag")
        if not r.in_v5_reextract:
            why.append("row absent from fixed v5 re-extraction (known-bad class)")
        reasons.append("; ".join(why))
    out["drop_reason"] = reasons
    out["kept"] = out.drop_reason == ""
    return out, pd.DataFrame(bad_summary)


def lab_rerun(lab: pd.DataFrame, reg_override: dict[str, pd.DataFrame]):
    rows = []
    for r in lab.itertuples():
        smi, x = json.loads(r.comp_smiles), json.loads(r.x_mole)
        rec = dict(system_id=r.system_id, measured_C=r.measured_C, comp_smiles=r.comp_smiles, x_mole=r.x_mole,
                   has_mesitylene=MES_SMILES in smi, v4_ideal_C=r.ideal_C, v4_unifac_C=r.unifac_C)
        for tag, reg in reg_override.items():
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", UnifacFallbackWarning)
                i = flash_point_unifac_detail(smi, x, registry=reg, use_unifac=False)
                u = flash_point_unifac_detail(smi, x, registry=reg, use_unifac=True)
            rec[f"{tag}_ideal_C"], rec[f"{tag}_unifac_C"] = i["flash_C"], u["flash_C"]
            rec[f"{tag}_ideal_err_C"], rec[f"{tag}_unifac_err_C"] = i["flash_C"] - r.measured_C, u["flash_C"] - r.measured_C
        rows.append(rec)
    return pd.DataFrame(rows)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    st = stored_sets()
    # ---- measured points and the drop list
    meas, bad = measured_rows()
    meas.to_csv(OUT / "mesitylene_measured_psat.csv", index=False)
    bad.to_csv(OUT / "thermoml_bad_vp_rows_summary.csv", index=False)
    kept = meas[meas.kept].drop_duplicates(subset=["file", "T_K", "value"])
    kept_low = kept[kept.T_K <= T_LOW_FOR_REFIT_K]

    # ---- Psat comparison (kPa)
    import thermo
    chem = thermo.Chemical(CAS)
    grid = [310.95, 317.0, 320.0, 323.0, 345.0, st["lo"], 400.0, 437.83, st["hi"]] + list(kept.T_K)
    rows = []
    for T in sorted(set(grid)):
        chem.calculate(T=T)
        mp = kept[np.isclose(kept.T_K, T)]
        rows.append(dict(T_K=T, stored_webbook_kPa=float(psat(T, st["webbook"])), anchored_registry_kPa=float(psat(T, st["anchored"])),
                         thermo_reference_kPa=chem.Psat / 1000.0, measured_kPa=float(mp.value.iloc[0]) if len(mp) else np.nan,
                         inside_stored_webbook_range=bool(st["lo"] <= T <= st["hi"])))
    cmp_ = pd.DataFrame(rows)
    cmp_["anchored_over_stored"] = cmp_.anchored_registry_kPa / cmp_.stored_webbook_kPa
    cmp_["stored_over_thermo_ref"] = cmp_.stored_webbook_kPa / cmp_.thermo_reference_kPa
    cmp_["anchored_over_thermo_ref"] = cmp_.anchored_registry_kPa / cmp_.thermo_reference_kPa
    cmp_.to_csv(OUT / "mesitylene_psat_compare.csv", index=False, float_format="%.6g")

    # ---- R-A2 decision
    dec = dict(n_rows_found_all_tables=len(meas), n_rows_kept=int(meas.kept.sum()), n_distinct_points_kept=len(kept),
               n_distinct_points_kept_le_345K=len(kept_low), min_points_required=MIN_POINTS_REFIT,
               refit_made=bool(len(kept_low) >= MIN_POINTS_REFIT),
               statement=("refit made" if len(kept_low) >= MIN_POINTS_REFIT else
                          "no independent low-T measurement on disk (D7 stays open)"))
    pd.DataFrame([dec]).to_csv(OUT / "mesitylene_refit_decision.csv", index=False)

    # ---- choose the best available input: smallest log-error at the measured point(s), fixed rule
    score = {}
    for tag in ("webbook", "anchored"):
        score[tag] = float(np.mean([abs(math.log(psat(r.T_K, st[tag]) / r.value)) for r in kept.itertuples()])) if len(kept) else np.nan
    best = min(score, key=score.get)
    best_abc = st[best]
    json.dump(dict(best=best, abc=best_abc, abc_convention="log10(P/kPa)=A-B/(T/K+C)", mean_abs_ln_error_at_measured=score,
                   rule="smallest mean |ln(Pcalc/Pmeas)| over kept measured points (only 437.83 K available)",
                   valid_range_K=[st["lo"], st["hi"]]), open(OUT / "mesitylene_best_input.json", "w"), indent=1)

    # ---- flash rerun
    lab = pd.read_csv(ROOT / "results/v4/F/flash_lab_loso.csv")
    base = _default_registry()
    regs = {"registry": base}
    for tag in ("webbook",):
        r2 = base.copy()
        r2.loc[st["can"], ["antoine_A", "antoine_B", "antoine_C"]] = list(st[tag])
        regs[tag] = r2
    rr = lab_rerun(lab, regs)
    chk = float(np.nanmax(np.abs(rr.registry_ideal_C - rr.v4_ideal_C)))
    chk_u = float(np.nanmax(np.abs(rr.registry_unifac_C - rr.v4_unifac_C)))
    rr.to_csv(OUT / "mesitylene_flash_rerun.csv", index=False, float_format="%.6g")
    assert chk < 1e-3 and chk_u < 1e-3, f"registry rerun does not reproduce v4 predictions ({chk}, {chk_u})"

    # ---- v4 filter logic, per input
    HEX = "CCCCCCc1ccccc1"
    has_h = rr.comp_smiles.apply(lambda s: HEX in json.loads(s))
    has_m = rr.has_mesitylene
    srows = []
    for inp in ("registry", "webbook"):
        d0 = rr.assign(ideal_err_C=rr[f"{inp}_ideal_err_C"], unifac_err_C=rr[f"{inp}_unifac_err_C"])
        for label, mask in [("all", pd.Series(True, index=rr.index)), ("excl_mesitylene", ~has_m), ("excl_n-hexylbenzene", ~has_h),
                            ("excl_both", ~(has_m | has_h)), ("only_mesitylene_records", has_m)]:
            d = d0[mask]
            for model in ("ideal", "unifac"):
                col = f"{model}_err_C"
                m, ng = margin90(d, col)
                srows.append(dict(mesitylene_input=inp, subset=label, model=model, n_records=len(d), n_systems=d.system_id.nunique(),
                                  MAE_C=d[col].abs().mean(), bias_C=d[col].mean(), worst_over_C=d[col].max(), margin90_C=m, margin_n_groups=ng))
    sens = pd.DataFrame(srows)
    sens.to_csv(OUT / "mesitylene_flash_sensitivity.csv", index=False, float_format="%.6g")

    # ---- readings
    mrec = rr[has_m]
    ra = []
    for inp in ("registry", "webbook"):
        for model in ("ideal", "unifac"):
            T = mrec[f"{inp}_{model}_C"] + 273.15
            ra.append(dict(reading="R-A1", input=inp, model=model, n_mesitylene_records=len(mrec), T_flash_min_K=T.min(), T_flash_max_K=T.max(),
                           set_valid_K=f"{st['lo']}-{st['hi']}", n_outside=int(((T < st["lo"]) | (T > st["hi"])).sum()),
                           extrapolated=bool(((T < st["lo"]) | (T > st["hi"])).any())))
    ra1 = pd.DataFrame(ra)
    ra3 = []
    for model in ("ideal", "unifac"):
        a = sens[(sens.mesitylene_input == best) & (sens.subset == "all") & (sens.model == model)].margin90_C.iloc[0]
        e = sens[(sens.mesitylene_input == best) & (sens.subset == "excl_mesitylene") & (sens.model == model)].margin90_C.iloc[0]
        v4a = sens[(sens.mesitylene_input == "registry") & (sens.subset == "all") & (sens.model == model)].margin90_C.iloc[0]
        ok = bool(np.isfinite(a) and np.isfinite(e) and abs(a - e) <= 1.0)
        ra3.append(dict(reading="R-A3", model=model, best_input=best, margin90_all_best_input_C=a, margin90_excl_mesitylene_C=e,
                        abs_diff_C=abs(a - e), margin90_all_registry_as_v4_C=v4a, holds=ok))
    ra3 = pd.DataFrame(ra3)
    ra3["R_A3_overall_holds"] = bool(ra3.holds.all())
    pd.concat([ra1, ra3], ignore_index=True).to_csv(OUT / "mesitylene_readings.csv", index=False, float_format="%.6g")

    # ---- context: CRC 647 r3 reconstruction usage of mesitylene
    comp = pd.read_csv(ROOT / "data/real_fuels/v4/compositions_crc647_r3.csv")
    comp = comp[~comp.is_lumped.astype(bool)]
    tot = comp.groupby("system_id").fraction.sum()
    mm = comp[comp.smiles == MES_SMILES].groupby("system_id").fraction.sum()
    share = (mm / tot.reindex(mm.index))
    pd.DataFrame([dict(n_fuels=comp.system_id.nunique(), n_fuels_with_mesitylene=len(mm), mass_share_min=share.min(),
                       mass_share_mean=share.mean(), mass_share_max=share.max(),
                       note="CRC 647 r3 flash NOT rerun with the corrected mesitylene input (only the lab filter was)")]).to_csv(
        OUT / "mesitylene_crc647_usage.csv", index=False, float_format="%.6g")

    pd.set_option("display.width", 250)
    pd.set_option("display.max_columns", 30)
    print(cmp_.to_string())
    print(meas.to_string())
    print(bad.to_string())
    print(pd.DataFrame([dec]).T)
    print("best input:", best, score, "registry reproduction max |diff| ideal/unifac:", chk, chk_u)
    print(sens.to_string())
    print(ra1.to_string())
    print(ra3.to_string())


if __name__ == "__main__":
    main()
