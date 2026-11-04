#!/usr/bin/env python3
"""Reproduce the Leenhouts et al. (arXiv 2605.29698) mixture-benchmark baselines.

Runs from the repo root (fuelblend_data/) using ONLY the released benchmark files
(read-only under data/interim/leenhouts_benchmark/) plus this repo's
data/unified/compounds.csv for the Le Chatelier flash baseline.

  1. Loads their processed mixture/pure CSVs and their 5-fold split definitions.
  2. Reproduces their "ideal-mixture" baseline by replicating
     MixturePureMatcher.compute_ideal_value (mole-fraction-weighted average of
     matched pure-compound values; temperature tolerance +/-1 K for viscosity):
       - viscosity: values are ln(mPa*s), so the ideal is exactly the
         Grunberg-Nissan log-linear rule with k=0 (their viscosity baseline).
       - flash point: plain linear average (their flash baseline; NOT Le
         Chatelier).
     Validated against their precomputed *_with_ideal.csv `ideal_value` column.
  3. Adds the two physics baselines requested for this repo:
       - flash: Le Chatelier + Raoult + Antoine (compounds.csv, log10(P/kPa)
         convention), LFL from compounds.csv with the Zabetakis T-correction.
         A blend is skipped (counted, never guessed) if any positive-fr
         component lacks Antoine constants or LFL.
       - viscosity: Refutas VBN averaging using mass fractions converted from
         mole fractions via rdkit MW. shortcut: pure dynamic viscosity in
         mPa*s is used directly as cSt (no densities in their pure file).
  4. Scores every baseline on the held-out test rows of every fold of every
     released split type; reports RMSE/MAE (mean +/- std over completed folds)
     and total evaluated rows, mirroring the paper's tables.

Determinism: no RNG anywhere; np.random.seed(0) set as a documented formality.
Output: stdout only (READ-ONLY on data/).
"""

from __future__ import annotations

import glob
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import kendalltau

np.random.seed(0)  # no RNG is used anywhere; formality only

REPO = Path(__file__).resolve().parents[1]
ROOT = REPO / "data" / "interim" / "leenhouts_benchmark"
SPLITS = ROOT / "splits" / "splits"
FLASH_DIR = ROOT / "datasets" / "datasets" / "flashpoint" / "processed_data"
VISC_DIR = ROOT / "datasets" / "datasets" / "viscosity" / "processed_data"
REGISTRY = REPO / "data" / "unified" / "compounds.csv"

P_ATM_KPA = 101.325
FLASH_LO_K, FLASH_HI_K = 233.15, 573.15   # closed-cup search window
TEMP_TOL_K = 1.0                          # paper's default temperature tolerance
DECANE_INCHI = "InChI=1S/C10H22/c1-3-5-7-9-10-8-6-4-2/h3-10H2,1-2H3"


# ---------------------------------------------------------------------------
# data loading + ideal-value machinery (replica of their matcher semantics)
# ---------------------------------------------------------------------------
def load_pair(folder: Path, name: str) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    mix = pd.read_csv(folder / f"{name}.csv")
    ideal = pd.read_csv(folder / f"{name}_with_ideal.csv")
    pure = pd.read_csv(folder / f"{name.replace('_mix_', '_pure_')}.csv")
    return mix, ideal, pure


def pure_lookup(pure: pd.DataFrame, with_temperature: bool) -> dict[tuple, float]:
    """Group-mean pure values, keyed as the paper does.
    flash:  (unit, inchi);   viscosity: (unit, temperature, temperature_unit, inchi)."""
    key_cols = ["unit"] + (["temperature", "temperature_unit"] if with_temperature else [])
    meaned = pure.groupby(key_cols + ["component_inchi_0"], as_index=False, dropna=False)["value"].mean()
    out: dict[tuple, float] = {}
    for d in meaned.to_dict(orient="records"):
        if pd.isna(d["value"]):
            continue
        t = round(float(d["temperature"]), 8) if with_temperature else None
        tu = d["temperature_unit"] if with_temperature else None
        out[(d["unit"], t, tu, d["component_inchi_0"])] = float(d["value"])
    return out


def positive_components(row) -> list[int]:
    out = []
    for i in range(4):
        comp = row.get(f"component_inchi_{i}")
        frac = pd.to_numeric(row.get(f"fraction_{i}"), errors="coerce")
        if isinstance(comp, str) and comp.startswith("InChI=") and pd.notna(frac) and float(frac) > 0.0:
            out.append(i)
    return out


def compute_ideal(row, lookup: dict[tuple, float], with_temperature: bool) -> float | None:
    """Weighted avg of matched pure values (their compute_ideal_value).
    Uses +frac_sum denominator and +/-1 K temperature candidates, exactly as theirs."""
    comps = positive_components(row)
    if not comps:
        return None
    unit = row.get("unit")
    if with_temperature:
        t_raw = pd.to_numeric(row.get("temperature"), errors="coerce")
        if pd.isna(t_raw):
            return None
        t_center = float(t_raw)
        tunit = row.get("temperature_unit")
        candidates = sorted({round(t_center + d, 8) for d in (-TEMP_TOL_K, 0.0, TEMP_TOL_K)})
    else:
        candidates = [None]
        tunit = None
    for t_c in candidates:
        keys, ok = [], True
        for i in comps:
            key = (unit, t_c, tunit, str(row[f"component_inchi_{i}"]))
            if key not in lookup:
                ok = False
                break
            keys.append((float(row[f"fraction_{i}"]), lookup[key]))
        if ok:
            den = sum(f for f, _ in keys)
            if den <= 0.0:
                return None
            return sum(f * v for f, v in keys) / den
    return None


# ---------------------------------------------------------------------------
# baseline functions
# ---------------------------------------------------------------------------
def refutas_ln_eta_from_mass(w_mass: list[float], eta_cSt: list[float]) -> float | None:
    """Refutas VBN: VBN = 14.534*ln(ln(v+0.8)) + 10.975, volume... no, MASS-fraction
    average VBN (standard oil-blending form). Needs v > 0.2 (VBN undefined below)."""
    if any(v <= 0.2 for v in eta_cSt):
        return None
    vbn = [14.534 * np.log(np.log(v + 0.8)) + 10.975 for v in eta_cSt]
    wtot = sum(w_mass)
    if wtot <= 0:
        return None
    vbn_mix = sum(w * b for w, b in zip(w_mass, vbn)) / wtot
    eta = np.exp(np.exp((vbn_mix - 10.975) / 14.534)) - 0.8
    if not np.isfinite(eta) or eta <= 0.0:
        return None
    return float(np.log(eta))


def lechatelier_flash_C(x_mole, A, B, C, lfl_vol):
    """Closed-cup flash: solve sum_i x_i*Psat_i(T)/(P_atm*LFL_i(T)) = 1 for T [ -> °C].
    Antoine (kPa): log10(P/kPa) = A - B/(T_K + C). LFL: vol%/100 mole fraction,
    Zabetakis LFL(T)=LFL_25*(1-7.21e-4*(T_C-25)). Bisection; H increases in T."""
    def f(T):
        h = 0.0
        for xi, Ai, Bi, Ci, Li in zip(x_mole, A, B, C, lfl_vol):
            psat = 10.0 ** (Ai - Bi / (T + Ci))
            lf_t = max((Li / 100.0) * (1.0 - 0.000721 * (T - 273.15 - 25.0)), 1e-9)
            h += xi * psat / (P_ATM_KPA * lf_t)
        return h - 1.0

    f_lo, f_hi = f(FLASH_LO_K), f(FLASH_HI_K)
    if f_lo > 0.0:
        return None
    if f_hi < 0.0:
        return None
    lo, hi = FLASH_LO_K, FLASH_HI_K
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if f(mid) > 0.0:
            hi = mid
        else:
            lo = mid
        if hi - lo < 1e-7:
            break
    return 0.5 * (lo + hi) - 273.15


# ---------------------------------------------------------------------------
# metrics
# ---------------------------------------------------------------------------
def fold_metrics(y: np.ndarray, yhat: np.ndarray) -> dict[str, float]:
    err = y - yhat
    tau = kendalltau(y, yhat).correlation
    return {"rmse": float(np.sqrt(np.mean(err**2))), "mae": float(np.mean(np.abs(err))),
            "kendall": float(tau) if tau == tau else np.nan, "n": int(len(y))}


def agg(per_fold: list[dict]) -> dict:
    got = [d for d in per_fold if d["n"] > 0]
    if not got:
        return {"rmse_mean": np.nan, "rmse_std": np.nan, "mae_mean": np.nan, "mae_std": np.nan,
                "kendall_mean": np.nan, "n_total": 0, "n_folds": 0}
    rm = np.array([d["rmse"] for d in got]); ma = np.array([d["mae"] for d in got])
    ta = np.array([d["kendall"] for d in got if d["kendall"] == d["kendall"]])
    return {"rmse_mean": rm.mean(), "rmse_std": rm.std(ddof=0), "mae_mean": ma.mean(),
            "mae_std": ma.std(ddof=0),
            "kendall_mean": ta.mean() if len(ta) else np.nan,
            "n_total": sum(d["n"] for d in got), "n_folds": len(got)}


# ---------------------------------------------------------------------------
# evaluation driver
# ---------------------------------------------------------------------------
def evaluate(tag: str, mix: pd.DataFrame, split_dirs: list[str], baselines: dict[str, tuple],
             unit_note: str, mix_basename: str, ideal_series: pd.Series) -> dict[str, dict]:
    """Score baselines on each fold's held-out test rows.
    Per baseline:
      - if require_ideal: scored on test rows with pure coverage (same subset as
        the paper's ideal baseline -> apples-to-apples),
      - else: scored on ALL test rows where the baseline itself can predict.
    Skips (inapplicable rows: no coverage / missing constants / VBN undefined)
    are counted, never guessed."""
    print(f"\n{'=' * 96}\n{tag}  (unit: {unit_note}; test rows of the released splits)\n{'=' * 96}")
    results = {}
    for sd in split_dirs:
        skips = {b: 0 for b in baselines}
        folds = sorted(glob.glob(str(SPLITS / sd / "cross_validation_5fold" / "fold_*.csv")))
        print(f"\n--- split: {Path(sd).parts[-1]}  ({len(folds)} folds) ---")
        print(f"{'baseline':34} {'RMSE mean±std':>22}   {'MAE mean±std':>20} {'tau':>5} {'n':>8} {'skip':>6}")
        per_fold = {b: [] for b in baselines}
        for ff in folds:
            spl = pd.read_csv(ff)
            spl = spl[spl["source"].eq("mix") & spl["source_file"].str.endswith(mix_basename)]
            test_rows = spl.loc[spl["split"].eq("test"), "source_row_index"].astype(int).to_numpy()
            for bname, (fn, req) in baselines.items():
                y, yhat = [], []
                for ridx in test_rows:
                    row = mix.iloc[ridx]
                    if req and not np.isfinite(ideal_series.iloc[ridx]):
                        skips[bname] += 1  # invisible to the coverage-based comparison set
                        continue
                    pred = fn(row)
                    if pred is None:
                        skips[bname] += 1  # baseline-inapplicable (missing constants etc.)
                        continue
                    y.append(float(row["value"])); yhat.append(float(pred))
                if len(y) < 2:
                    per_fold[bname].append({"rmse": 0.0, "mae": 0.0, "kendall": np.nan, "n": 0})
                else:
                    per_fold[bname].append(fold_metrics(np.array(y), np.array(yhat)))
        for bname in baselines:
            a = agg(per_fold[bname])
            results[(Path(sd).parts[-1], bname)] = a
            ks = f"{a['kendall_mean']:.2f}" if a["kendall_mean"] == a["kendall_mean"] else "-"
            print(f"{bname:34} {a['rmse_mean']:>11.4f} ±{a['rmse_std']:<9.4f}  "
                  f"{a['mae_mean']:>9.4f} ±{a['mae_std']:<9.4f} {ks:>5} {a['n_total']:>8d} {skips[bname]:>6d}")
    return results


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------
def main() -> int:
    from rdkit.Chem.inchi import MolFromInchi, MolToInchiKey
    from rdkit import Chem
    from rdkit.Chem.Descriptors import MolWt

    print("Leenhouts et al. mixture-benchmark baseline reproduction")
    print("(no RNG used; seed 0 set as determinism statement)\n")

    # row-count sanity against the task description
    mix, ideal, pure = load_pair(FLASH_DIR, "flashpoint_mix_exp")
    print(f"[loaded] flashpoint_pure_exp.csv={len(pure)} rows (expect 585), "
          f"flashpoint_mix_exp.csv={len(mix)} rows (expect 1006)")
    mix2, ideal2, pure2 = load_pair(VISC_DIR, "viscosity_mix_exp")
    print(f"[loaded] viscosity_pure_exp.csv={len(pure2)} rows (expect 3027), "
          f"viscosity_mix_exp.csv={len(mix2)} rows (expect 35938)")

    # unit sanity checks with verifiable numbers
    fp_dec = pure.loc[pure["component_inchi_0"].eq(DECANE_INCHI), "value"]
    print(f"[unit check] n-decane pure flash point = {float(fp_dec.iloc[0]) if len(fp_dec) else 'N/A'} °C (lit ~46)")
    dec_v = pure2.loc[pure2["component_inchi_0"].eq(DECANE_INCHI)]
    if len(dec_v):
        print(f"[unit check] n-decane pure eta = {np.exp(dec_v['value'].iloc[0]):.3f} mPa·s at {dec_v['temperature'].iloc[0]} K (lit ~0.92 @298 K)")

    # ------------------ flash ------------------
    print("\n" + "#" * 96 + "\n# FLASH POINT MIXTURES  (their baseline: mole-weighted linear avg of pure FPs)\n" + "#" * 96)
    lookup_f = pure_lookup(pure, with_temperature=False)
    ideal_f = pd.Series([compute_ideal(mix.loc[i], lookup_f, False) for i in mix.index], index=mix.index)
    LCHAT_SKIP = {"not_in_registry": 0, "missing_antoin_or_lfl": 0}
    n_cov = int(ideal_f.notna().sum())
    print(f"[coverage] rows with full pure-compound coverage: {n_cov}/{len(mix)}")

    # validate vs precomputed ideal_value
    theirs = ideal["ideal_value"].to_numpy(dtype=float)
    mask_pre = np.isfinite(theirs)
    mine = ideal_f.to_numpy(dtype=float)[mask_pre]
    agree = np.isfinite(mine) & (np.abs(mine - theirs[mask_pre]) <= 1e-9)
    print(f"[self-check] recomputed ideal matches their *_with_ideal.csv: {int(agree.sum())}/{int(mask_pre.sum())} "
          f"exact to 1e-9; max|diff| on agreed rows = 0 by construction"
          f" (disagreements where mine invalid: {int(np.sum(~np.isfinite(mine)))})")

    reg = pd.read_csv(REGISTRY)
    by_inchi = {r.inchi: r for r in reg.itertuples(index=False) if isinstance(r.inchi, str)
                and isinstance(r.inchi, str)}
    by_key = {r.inchikey: r for r in reg.itertuples(index=False) if isinstance(r.inchikey, str)}
    comp_rec: dict[str, object | None] = {}

    def comp_record(inchi: str):
        if inchi not in comp_rec:
            r = by_inchi.get(inchi)
            if r is None:
                mol = MolFromInchi(inchi)
                k = MolToInchiKey(mol) if mol is not None else None
                r = by_key.get(k) if k else None
            comp_rec[inchi] = r
        return comp_rec[inchi]

    def flash_lechatelier(row) -> float | None:
        xs, A, B, C, L = [], [], [], [], []
        for i in positive_components(row):
            r = comp_record(str(row[f"component_inchi_{i}"]))
            if r is None:
                LCHAT_SKIP["not_in_registry"] += 1
                return None  # hydrocarbon registry gap: do not guess
            if not (pd.notna(r.antoine_A) and pd.notna(r.antoine_B) and pd.notna(r.antoine_C)
                    and pd.notna(r.lfl_vol_pct)):
                LCHAT_SKIP["missing_antoin_or_lfl"] += 1
                return None
            xs.append(float(row[f"fraction_{i}"]))
            A.append(float(r.antoine_A)); B.append(float(r.antoine_B)); C.append(float(r.antoine_C))
            L.append(float(r.lfl_vol_pct))
        s = sum(xs)
        if s <= 0:
            return None
        return lechatelier_flash_C([x / s for x in xs], A, B, C, L)

    def flash_ideal(row) -> float | None:
        return compute_ideal(row, lookup_f, False)

    res_f = evaluate("FLASHPOINT", mix,
                     ["flashpoint/flashpoint_mix_exp/split_definitions/random_rowwise",
                      "flashpoint/flashpoint_mix_exp/split_definitions/mixture",
                      "flashpoint/flashpoint_mix_exp/split_definitions/molecule",
                      "flashpoint/flashpoint_mix_exp/split_definitions/pure_to_mixture"],
                     {"ideal (their baseline)": (flash_ideal, True),
                      "Le Chatelier (Antoine+LFL)": (flash_lechatelier, False)},
                     "°C (mix value & prediction)", "flashpoint_mix_exp.csv", ideal_f)
    print(f"[Le Chatelier skips: by reason (components, cumulative over evaluated rows)] {LCHAT_SKIP}")

    # ------------------ viscosity ------------------
    print("\n" + "#" * 96 + "\n# VISCOSITY MIXTURES  (their baseline: Grunberg-Nissan k=0 log-linear in ln(mPa·s))\n" + "#" * 96)
    lookup_v = pure_lookup(pure2, with_temperature=True)
    ideal_v = pd.Series([compute_ideal(mix2.loc[i], lookup_v, True) for i in mix2.index], index=mix2.index)
    n_cov2 = int(ideal_v.notna().sum())
    print(f"[coverage] rows with full pure-compound coverage (T±1 K): {n_cov2}/{len(mix2)}")
    theirs2 = ideal2["ideal_value"].to_numpy(dtype=float)
    mask2 = np.isfinite(theirs2)
    mine2 = ideal_v.to_numpy(dtype=float)[mask2]
    agree2 = np.isfinite(mine2) & (np.abs(mine2 - theirs2[mask2]) <= 1e-9)
    print(f"[self-check] recomputed ideal matches their *_with_ideal.csv: {int(agree2.sum())}/{int(mask2.sum())} exact to 1e-9")

    mw_cache: dict[str, float] = {}
    def mw_of(inchi: str) -> float:
        if inchi not in mw_cache:
            mol = MolFromInchi(inchi)
            mw_cache[inchi] = float(MolWt(mol)) if mol is not None else float("nan")
        return mw_cache[inchi]

    def visc_loglinear(row) -> float | None:
        return compute_ideal(row, lookup_v, True)

    def visc_refutas(row) -> float | None:
        comps = positive_components(row)
        if not comps:
            return None
        t_raw = pd.to_numeric(row.get("temperature"), errors="coerce")
        if pd.isna(t_raw):
            return None
        t0, unit, tunit = float(t_raw), row.get("unit"), row.get("temperature_unit")
        x, eta_cSt, mw = [], [], []
        for i in comps:
            comp = str(row[f"component_inchi_{i}"])
            val = None
            for tc in sorted({round(t0 + d, 8) for d in (-TEMP_TOL_K, 0.0, TEMP_TOL_K)}):
                k = (unit, tc, tunit, comp)
                if k in lookup_v:
                    val = lookup_v[k]; break
            if val is None:
                return None  # pure eta at matched T unavailable -> Ideal/Refutas both inapplicable
            m = mw_of(comp)
            if not np.isfinite(m):
                return None
            x.append(float(row[f"fraction_{i}"])); eta_cSt.append(float(np.exp(val))); mw.append(m)
        tot = sum(xi * mi for xi, mi in zip(x, mw))
        if tot <= 0 or any(mi <= 0 for mi in mw):
            return None
        w = [xi * mi / tot for xi, mi in zip(x, mw)]
        return refutas_ln_eta_from_mass(w, eta_cSt)

    res_v = evaluate("VISCOSITY", mix2,
                     ["viscosity/viscosity_mix_exp/split_definitions/random_rowwise",
                      "viscosity/viscosity_mix_exp/split_definitions/mixture",
                      "viscosity/viscosity_mix_exp/split_definitions/molecule",
                      "viscosity/viscosity_mix_exp/split_definitions/pure_to_mixture",
                      "viscosity/viscosity_mix_exp/split_definitions/mixture_temperature"],
                     {"log-linear (their baseline)": (visc_loglinear, True),
                      "Refutas VBN (mass frac)": (visc_refutas, False)},
                     "ln(mPa·s) (mix value & prediction)", "viscosity_mix_exp.csv", ideal_v)

    # ------------------ summary table (compare vs paper) ------------------
    print("\n" + "#" * 96 + "\n# SUMMARY vs PAPER (paper: mean ± std over folds; T_flash decimals clipped in their PDF)\n" + "#" * 96)
    paper = {
        ("random_rowwise", "ideal"): (0.30, 0.015, "ln(mPa·s) excess/abs RMSE, Tables 3 & 11"),
        ("mixture", "ideal"): (0.30, 0.035, "ln(mPa·s) excess/abs RMSE, Tables 3 & 11"),
        ("molecule", "ideal"): (0.21, 0.062, "ln(mPa·s) excess/abs RMSE, Tables 3 & 11"),
        ("random_rowwise", "flash"): ("14.x", 0.78, "°C — PDF renders '14.±0.78' (mean decimal clipped in their table column)"),
        ("mixture", "flash"): ("13.x", 5.4, "°C — PDF renders '13.±5.4'"),
        ("molecule", "flash"): ("13.x", 8.3, "°C — PDF renders '13.±8.3'"),
    }
    print("\n[viscosity ideal-baseline RMSE reproduction]")
    for split in ["random_rowwise", "mixture", "molecule"]:
        a = res_v.get((split, "log-linear (their baseline)"))
        if a is None: continue
        p = paper[(split, "ideal")]
        print(f"  {split:>16}: reproduced {a['rmse_mean']:.3f}±{a['rmse_std']:.3f}   paper {p[0]}±{p[1]}")
    print("\n[flashpoint ideal-baseline RMSE reproduction]")
    for split in ["random_rowwise", "mixture", "molecule"]:
        a = res_f.get((split, "ideal (their baseline)"))
        if a is None: continue
        p = paper[(split, "flash")]
        print(f"  {split:>16}: reproduced {a['rmse_mean']:.3f}±{a['rmse_std']:.3f}   paper {p[0]}±{p[1]} ({p[2]})")

    print("\nDone.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
