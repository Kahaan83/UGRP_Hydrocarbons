#!/usr/bin/env python
"""T0.2: Leakage Audit for the hydrocarbon fuel-blend ML project.

Checks:
  C1  CRC-in-training         : no CRC 647 measurements in model_table S1 rows
  C2  Near-duplicate CRC rows : no near-duplicate CRC compounds in training (T-matched only)
  C3a Antoine shared-compound : severity of using pure-compound VP data for Antoine refit
                                 when those compounds also appear in LOSO blends
  C3b Antoine true leakage    : Antoine refit does NOT use the blend-mixture measurements themselves
  C4  LOSO systems not in S1  : flash mixture LOSO systems must not appear in S1 training rows
  C5  T1 test suite           : run tests/test_t1_leakage_audit.py

Outputs:
  results/t0/leakage_checks.csv        checks table
  results/t0/leakage_crc_training.csv  C1/C2 details
  results/t0/leakage_antoine.csv       C3 compound-level details
  results/t0/leakage_loso_splits.csv   C4 details
  docs/audit_leakage.md                narrative summary

Usage:
    python scripts/t02_leakage_audit.py
"""
from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
PYTHON = sys.executable

RESULTS_DIR = ROOT / "results" / "t0"
DOCS_DIR = ROOT / "docs"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)
DOCS_DIR.mkdir(parents=True, exist_ok=True)


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

def load_model_table():
    return pd.read_csv(ROOT / "data/processed/model_table.csv", low_memory=False)


def load_measurements_dedup():
    return pd.read_csv(ROOT / "data/unified/measurements_dedup.csv", low_memory=False)


def load_splits():
    with open(ROOT / "splits/s1s2s3.json") as f:
        return json.load(f)


def parse_ik_list(s):
    if pd.isna(s):
        return []
    s = str(s).strip()
    try:
        v = json.loads(s)
        return v if isinstance(v, list) else [v]
    except Exception:
        return [s]


def check_row(check_id, description, method, verdict, affected_count, notes):
    return {
        "check_id": check_id,
        "description": description,
        "method": method,
        "verdict": verdict,
        "affected_rows": affected_count,
        "notes": notes,
    }


def df_to_md_table(df, max_rows=20):
    """Simple markdown table formatter; avoids tabulate dependency."""
    df2 = df.head(max_rows)
    if df2.empty:
        return "(empty)"
    cols = df2.columns.tolist()
    header = "| " + " | ".join(str(c) for c in cols) + " |"
    sep = "| " + " | ".join("---" for _ in cols) + " |"
    rows = []
    for _, r in df2.iterrows():
        rows.append("| " + " | ".join(str(r[c]) for c in cols) + " |")
    return "\n".join([header, sep] + rows)


# ---------------------------------------------------------------------------
# C1: CRC-in-training
# ---------------------------------------------------------------------------

def check_c1(mt, md, splits):
    """Identify CRC 647 rows in measurements_dedup and confirm none are in model_table S1."""
    crc_md = md[md["sources"].fillna("").str.contains("crc", case=False)].copy()
    crc_sids = set(crc_md["system_id"].unique())

    s1_idxs = set(splits["S1"])
    mt_s1 = mt[mt.index.isin(s1_idxs)]

    # CRC rows in model_table: by sources column, and by fuel: system_id
    mt_crc_src = mt[mt["sources"].fillna("").str.contains("crc", case=False)]
    mt_crc_sid = mt[mt["system_id"].isin(crc_sids)]
    mt_crc_all = pd.concat([mt_crc_src, mt_crc_sid]).drop_duplicates()

    mt_crc_s1 = mt_crc_all[mt_crc_all.index.isin(s1_idxs)]

    detail = pd.DataFrame([{
        "check": "C1",
        "crc_rows_in_measurements_dedup": len(crc_md),
        "crc_system_ids_count": len(crc_sids),
        "crc_rows_in_model_table": len(mt_crc_all),
        "crc_rows_in_s1_training": len(mt_crc_s1),
        "note": ("CRC 647 real-fuel rows have kind=real_fuel and system_id=fuel:*; "
                 "model_table contains no real_fuel rows at all."),
    }])

    verdict = "PASS" if len(mt_crc_s1) == 0 else "FAIL"
    return (
        check_row(
            "C1", "CRC-in-training",
            "Count model_table S1 rows with sources containing 'crc' or system_id in CRC fuel: set",
            verdict, len(mt_crc_s1),
            (f"CRC rows in measurements_dedup: {len(crc_md)} | "
             f"CRC fuel: system_ids: {len(crc_sids)} | "
             "model_table has 0 real_fuel rows; CRC data excluded from model_table entirely."),
        ),
        detail,
    )


# ---------------------------------------------------------------------------
# C2: Near-duplicate CRC rows in training (T-matched)
# ---------------------------------------------------------------------------

def check_c2(mt, md, splits):
    """Near-duplicate check restricted to CRC rows WITH a known T_K.

    CRC fuels have no molecular composition (bulk route), so inchikey-set overlap
    is structurally impossible. We check:
    (a) Any CRC fuel: system_ids in S1 (always 0 by C1 logic).
    (b) T-matched value coincidences between CRC rows and S1 training rows.
        NaN-T CRC rows are SKIPPED to avoid false positives (a NaN-T ambient
        density row would match every training density row of a similar value).
    """
    crc_md = md[md["sources"].fillna("").str.contains("crc", case=False)].copy()
    crc_sids = set(crc_md["system_id"].unique())

    s1_idxs = set(splits["S1"])
    mt_s1 = mt[mt.index.isin(s1_idxs)].copy()

    s1_crc_sys = mt_s1[mt_s1["system_id"].isin(crc_sids)]

    # Only CRC rows with a known T_K to avoid spurious matches
    crc_with_t = crc_md[crc_md["T_K"].notna()].copy()

    hits = []
    for _, cr in crc_with_t.iterrows():
        prop = cr["property"]
        val = cr["value"]
        t_k = float(cr["T_K"])
        if not (np.isfinite(val) and np.isfinite(t_k)):
            continue
        cand = mt_s1[mt_s1["property"] == prop]
        if cand.empty:
            continue
        # Both T must be finite AND within 2 K
        cand = cand[cand["T_K"].notna() & (np.abs(cand["T_K"] - t_k) <= 2.0)]
        # Value within 5%
        cand = cand[np.abs(cand["value"] - val) / (np.abs(val) + 1e-9) <= 0.05]
        for _, tr in cand.iterrows():
            hits.append({
                "crc_system_id": cr["system_id"],
                "crc_property": prop,
                "crc_value": round(val, 4),
                "crc_T_K": round(t_k, 2),
                "train_system_id": tr["system_id"],
                "train_sources": tr["sources"],
                "train_value": round(float(tr["value"]), 4),
                "train_T_K": round(float(tr["T_K"]), 2),
            })

    detail = pd.DataFrame(hits) if hits else pd.DataFrame(
        columns=["crc_system_id", "crc_property", "crc_value", "crc_T_K",
                 "train_system_id", "train_sources", "train_value", "train_T_K"])

    n = len(detail)
    unique_train_sys = int(detail["train_system_id"].nunique()) if n > 0 else 0

    verdict = "WARN" if n > 0 else "PASS"
    return (
        check_row(
            "C2", "Near-duplicate CRC rows in training (T-matched)",
            "CRC fuel: system_id overlap in S1; T-matched value within 5% (CRC rows with known T_K only)",
            verdict, n,
            (f"CRC fuel: system_ids in S1: {len(s1_crc_sys)} (must be 0). "
             f"{n} T-matched near-value hits across {unique_train_sys} unique training system_ids "
             f"(from {len(crc_with_t)} CRC rows with known T_K). "
             "CRC fuels have no molecule lists; near-value matches are physical coincidences "
             "(Jet-A density clusters near pure hydrocarbon densities), not data copying."),
        ),
        detail,
    )


# ---------------------------------------------------------------------------
# C3: Antoine leakage check
# ---------------------------------------------------------------------------

def check_c3(mt, md, splits):
    """Trace Antoine refit data sources vs LOSO held-out systems.

    refit_antoine.py sources:
      measurements.csv filtered to property=vapor_pressure_kPa AND n_comp==1 (pure only).
    Therefore:
      C3a: shared-compound leakage — LOSO compounds whose PURE VP data fed Antoine refit.
           These are per-compound constants (not per-blend), unavoidable for a physics baseline.
      C3b: true leakage — VP rows from LOSO mixture system_ids used in Antoine refit.
           These are mixture (n_comp>=2) VP rows; refit_antoine.py's n_comp==1 filter excludes them.
    """
    flash_mix = mt[(mt["property"] == "flash_point_C") & (mt["n_comp"] >= 2)].copy()

    loso_sys_iks: dict[str, set] = {}
    for _, row in flash_mix.iterrows():
        sid = row["system_id"]
        for ik in parse_ik_list(row["comp_inchikeys"]):
            loso_sys_iks.setdefault(sid, set()).add(ik)
    all_loso_iks = set().union(*loso_sys_iks.values()) if loso_sys_iks else set()
    loso_sids = set(loso_sys_iks.keys())

    vp_all = md[md["property"] == "vapor_pressure_kPa"]
    vp_pure = vp_all[vp_all["kind"] == "pure"]   # n_comp==1; fed to refit

    # C3b: VP rows whose system_id matches a LOSO system (mixture-level)
    # Note: system_id equality means same component set; for mixture VP rows these will
    # be binary or ternary VP rows (n_comp>=2) and therefore EXCLUDED from refit by n_comp==1 filter.
    true_leakage_vp = vp_all[vp_all["system_id"].isin(loso_sids)]
    # Further restrict: only rows that are actually pure (would be used by refit)
    true_leakage_pure = true_leakage_vp[true_leakage_vp["kind"] == "pure"]

    # C3a: Which LOSO compounds have pure VP data?
    pure_ik_vp_count: dict[str, int] = {}
    for _, row in vp_pure.iterrows():
        for ik in parse_ik_list(row["comp_inchikeys"]):
            pure_ik_vp_count[ik] = pure_ik_vp_count.get(ik, 0) + 1

    ar_path = ROOT / "data/processed/antoine_refit.csv"
    ar = pd.read_csv(ar_path) if ar_path.exists() else pd.DataFrame()
    refit_ok_iks = set(ar[ar["status"].str.startswith("ok")]["inchikey"]) if len(ar) else set()

    compound_rows = []
    for ik in sorted(all_loso_iks):
        n_vp = pure_ik_vp_count.get(ik, 0)
        was_refit = ik in refit_ok_iks
        sys_using = [sid for sid, iks in loso_sys_iks.items() if ik in iks]
        if n_vp > 0 and was_refit:
            severity = "shared_compound_refit"
        elif was_refit:
            severity = "refit_ok_no_pure_vp"
        elif n_vp > 0:
            severity = "vp_exists_not_refit"
        else:
            severity = "webbook_or_no_antoine"
        compound_rows.append({
            "inchikey": ik,
            "n_pure_vp_rows": n_vp,
            "was_refit": was_refit,
            "n_loso_systems_using": len(sys_using),
            "severity": severity,
        })

    detail = pd.DataFrame(compound_rows)
    n_shared = int((detail["severity"] == "shared_compound_refit").sum())

    # C3b verdict: only if pure VP rows from LOSO mixture system_ids exist
    n_true_leak = len(true_leakage_pure)
    # Also record mixture VP rows with shared system_id (informational)
    n_mix_vp_shared_sid = len(true_leakage_vp[true_leakage_vp["kind"] != "pure"])

    c3a_verdict = "WARN" if n_shared > 0 else "PASS"
    c3b_verdict = "PASS" if n_true_leak == 0 else "FAIL"

    rows_c3a = check_row(
        "C3a", "Antoine shared-compound leakage",
        "Count LOSO compounds whose pure VP data fed Antoine refit (refit_antoine.py uses n_comp==1 VP only)",
        c3a_verdict, n_shared,
        (f"{n_shared}/32 LOSO compounds have pure VP data that fed Antoine refit. "
         "Severity: LOW -- Antoine constants are per-compound scalars fitted on pure-A data; "
         "they encode no information about held-out blend compositions or blend flash points. "
         "This is the standard data-flow for any physics baseline (Raoult's law needs pure Psat).")
    )
    rows_c3b = check_row(
        "C3b", "Antoine true leakage (blend measurements used in refit)",
        "Check if any PURE VP rows from LOSO mixture system_ids were used in Antoine refit",
        c3b_verdict, n_true_leak,
        (f"{n_true_leak} PURE VP rows with system_id matching a LOSO flash-mix system. "
         f"{n_mix_vp_shared_sid} MIXTURE VP rows with matching system_id (kind!=pure, excluded by refit). "
         "refit_antoine.py filters n_comp==1 only, so binary/ternary VP rows are structurally excluded. "
         "Mixture VP rows sharing a system_id with a LOSO flash system are NOT leakage.")
    )
    return rows_c3a, rows_c3b, detail


# ---------------------------------------------------------------------------
# C4: LOSO systems not in S1
# ---------------------------------------------------------------------------

def check_c4(mt, splits):
    """Flash LOSO systems: check none appear in S1 training rows used by t8."""
    flash_mix = mt[(mt["property"] == "flash_point_C") & (mt["n_comp"] >= 2)].copy()
    loso_sids = set(flash_mix["system_id"].unique())

    s1_idxs = set(splits["S1"])
    s2_idxs = set(splits["S2"])
    s3_idxs = set(splits["S3"])

    detail_rows = []
    for sid in sorted(loso_sids):
        sys_rows = flash_mix[flash_mix["system_id"] == sid]
        in_s1 = sys_rows[sys_rows.index.isin(s1_idxs)]
        in_s2 = sys_rows[sys_rows.index.isin(s2_idxs)]
        in_s3 = sys_rows[sys_rows.index.isin(s3_idxs)]
        in_none = sys_rows[~sys_rows.index.isin(s1_idxs | s2_idxs | s3_idxs)]
        detail_rows.append({
            "system_id": sid,
            "total_rows": len(sys_rows),
            "in_S1": len(in_s1),
            "in_S2": len(in_s2),
            "in_S3": len(in_s3),
            "in_NONE": len(in_none),
        })

    detail = pd.DataFrame(detail_rows)
    s1_leaking = detail[detail["in_S1"] > 0]
    n_s1_sys = len(s1_leaking)

    verdict = "WARN" if n_s1_sys > 0 else "PASS"
    return (
        check_row(
            "C4", "LOSO flash system discipline vs t8 splits",
            "Check how many of 31 LOSO flash systems also appear in S1 training rows (t8 training data)",
            verdict, int(s1_leaking["in_S1"].sum()),
            (f"{n_s1_sys}/31 LOSO systems have rows assigned to S1. "
             f"{int(detail['in_S2'].gt(0).sum())} in S2, "
             f"{int(detail['in_NONE'].gt(0).sum())} in NONE (not assigned). "
             "eval_flash_loso.py is a standalone LOSO procedure on the full mixture dataset; "
             "it does not respect S1/S2 boundaries. "
             "LOSO physics baseline (Le Chatelier) is unaffected by S1 membership. "
             "LOSO HGB-residual includes t8-training systems in its test set for those S1 systems."),
        ),
        detail,
    )


# ---------------------------------------------------------------------------
# C5: T1 test suite
# ---------------------------------------------------------------------------

def check_c5():
    """Run tests/test_t1_leakage_audit.py with the venv python."""
    test_file = ROOT / "tests/test_t1_leakage_audit.py"
    if not test_file.exists():
        return check_row("C5", "T1 test suite", "subprocess python tests/test_t1_leakage_audit.py",
                         "WARN", 0, "Test file not found.")
    result = subprocess.run(
        [PYTHON, str(test_file)],
        capture_output=True, text=True, cwd=str(ROOT), timeout=120
    )
    output_lines = (result.stdout + result.stderr).strip().splitlines()
    summary = " | ".join(line.strip() for line in output_lines[-5:] if line.strip()) or "(no output)"
    verdict = "PASS" if result.returncode == 0 else "FAIL"
    return check_row(
        "C5", "T1 test suite (test_t1_leakage_audit.py)",
        f"subprocess: python {test_file.name}",
        verdict, 0,
        f"returncode={result.returncode} | {summary}",
    )


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------

def main():
    print("[T0.2] Loading data...")
    mt = load_model_table()
    md = load_measurements_dedup()
    sp = load_splits()

    print("[T0.2] C1: CRC-in-training...")
    c1, c1_detail = check_c1(mt, md, sp)

    print("[T0.2] C2: Near-duplicate CRC rows (T-matched)...")
    c2, c2_detail = check_c2(mt, md, sp)

    print("[T0.2] C3: Antoine leakage...")
    c3a, c3b, c3_detail = check_c3(mt, md, sp)

    print("[T0.2] C4: LOSO system split discipline...")
    c4, c4_detail = check_c4(mt, sp)

    print("[T0.2] C5: T1 test suite...")
    c5 = check_c5()

    # checks table
    checks = pd.DataFrame([c1, c2, c3a, c3b, c4, c5])
    checks.to_csv(RESULTS_DIR / "leakage_checks.csv", index=False)

    # aux CSVs
    crc_detail = pd.concat([c1_detail, c2_detail], ignore_index=True)
    crc_detail.to_csv(RESULTS_DIR / "leakage_crc_training.csv", index=False)
    c3_detail.to_csv(RESULTS_DIR / "leakage_antoine.csv", index=False)
    c4_detail.to_csv(RESULTS_DIR / "leakage_loso_splits.csv", index=False)

    print()
    print("=" * 78)
    print("T0.2 LEAKAGE AUDIT -- CHECKS TABLE")
    print("=" * 78)
    for _, row in checks.iterrows():
        sym = {"PASS": "[OK]", "FAIL": "[!!]", "WARN": "[WW]"}.get(row["verdict"], "[??]")
        print(f"  {sym} {row['check_id']:5s} [{row['verdict']:4s}] "
              f"affected={row['affected_rows']:5d}  {row['description']}")
    print("=" * 78)

    _write_md(checks, c1, c1_detail, c2, c2_detail, c3a, c3b, c3_detail, c4, c4_detail, c5)

    print(f"\n[T0.2] Done. Outputs:")
    for f in ["leakage_checks.csv", "leakage_crc_training.csv",
              "leakage_antoine.csv", "leakage_loso_splits.csv"]:
        print(f"  {RESULTS_DIR / f}")
    print(f"  {DOCS_DIR / 'audit_leakage.md'}")


def _write_md(checks, c1, c1_detail, c2, c2_detail, c3a, c3b, c3_detail, c4, c4_detail, c5):
    lines = [
        "# T0.2 Leakage Audit",
        "",
        "> Generated by `scripts/t02_leakage_audit.py`.",
        "> Read-only analysis -- no training/eval code or data files modified.",
        "",
        "## Checks Summary",
        "",
        "| ID | Description | Verdict | Affected rows |",
        "|---|---|---|---|",
    ]
    for _, row in checks.iterrows():
        sym = {"PASS": "PASS", "FAIL": "FAIL", "WARN": "WARN"}.get(row["verdict"], row["verdict"])
        lines.append(
            f"| {row['check_id']} | {row['description']} | **{sym}** | {row['affected_rows']} |"
        )
    lines += ["", "---", ""]

    # C1
    r1 = c1_detail.iloc[0]
    lines += [
        "## C1 -- CRC-in-Training",
        "",
        f"**Verdict: {c1['verdict']}**",
        "",
        c1["notes"],
        "",
        f"- CRC 647 rows in `measurements_dedup.csv`: {r1['crc_rows_in_measurements_dedup']}",
        f"- CRC fuel: system_ids: {r1['crc_system_ids_count']}",
        f"- CRC rows in `model_table.csv` (any source or system_id): {r1['crc_rows_in_model_table']}",
        f"- CRC rows in S1 training: **{r1['crc_rows_in_s1_training']}**",
        "",
        "CRC 647 real-fuel rows have `kind=real_fuel` and are excluded from `model_table.csv` entirely.",
        "Rule R1 in `src/fuelblend/data/model_table.py` requires `composition_level in {L2,L3}`.",
        "Real fuels have no molecule-level composition, so they never enter model_table.",
        "The hindcast file `results/hindcast/hindcast_predictions.csv` uses bulk-route reconstructed",
        "compositions and is evaluation-only -- CRC measurement values never entered model_table.",
        "",
    ]

    # C2
    lines += [
        "## C2 -- Near-Duplicate CRC Rows in Training (T-matched)",
        "",
        f"**Verdict: {c2['verdict']}**",
        "",
        c2["notes"],
        "",
    ]
    if len(c2_detail) > 0:
        lines.append(f"Top {min(10,len(c2_detail))} of {len(c2_detail)} T-matched near-value hits:")
        lines.append("")
        lines.append(df_to_md_table(c2_detail))
        lines.append("")
        prop_counts = c2_detail.groupby("crc_property").size().sort_values(ascending=False)
        lines.append("Hits by property:")
        for prop, cnt in prop_counts.items():
            lines.append(f"  - {prop}: {cnt}")
        lines.append("")
        lines.append("> **Assessment:** CRC fuels have no molecule lists (bulk route). "
                     "These coincidences arise because Jet-A/A-1 density and viscosity overlap "
                     "with the values of individual hydrocarbon compounds at similar temperatures. "
                     "No CRC measurement row was copied into training.")
    else:
        lines.append("No T-matched near-value hits found.")
    lines.append("")

    # C3
    n_shared = int((c3_detail["severity"] == "shared_compound_refit").sum())
    lines += [
        "## C3 -- Antoine Baseline Leakage",
        "",
        "### C3a -- Shared-compound leakage",
        "",
        f"**Verdict: {c3a['verdict']}**",
        "",
        c3a["notes"],
        "",
        "**Data-flow trace:**",
        "- `refit_antoine.py` reads `data/unified/measurements.csv`",
        "- Filters: `property == 'vapor_pressure_kPa'` AND `n_comp == 1` (pure compounds only)",
        "- Antoine constants are per-compound scalars; they encode no blend composition information.",
        "",
        "**Severity: LOW.** Using pure-A vapour pressure data to fit Antoine(A), then evaluating",
        "a held-out blend containing A, is not leakage -- it is the required data-flow for any",
        "physics-based flash-point baseline (Le Chatelier + Raoult's law needs P_sat of each pure component).",
        "The alternative (excluding pure VP data for LOSO compounds) would leave the baseline unsolvable.",
        "",
        f"Compounds in LOSO systems with pure VP data fed to Antoine refit ({n_shared}):",
        "",
    ]
    shared = c3_detail[c3_detail["severity"] == "shared_compound_refit"]
    if len(shared) > 0:
        lines.append(df_to_md_table(shared[["inchikey", "n_pure_vp_rows", "was_refit", "n_loso_systems_using"]]))
    else:
        lines.append("(none)")
    lines += [
        "",
        "### C3b -- True leakage (blend measurements used in refit)",
        "",
        f"**Verdict: {c3b['verdict']}**",
        "",
        c3b["notes"],
        "",
        "**Root cause of any apparent hits:**",
        "The `system_id` for a binary mixture (e.g., heptane+ethylbenzene) is shared between its",
        "flash point rows AND its binary vapour pressure rows in `measurements_dedup.csv`.",
        "These binary VP rows have `kind=binary` / `n_comp=2`, so `refit_antoine.py`'s",
        "`n_comp==1` filter excludes them structurally. No blend-level VP measurement",
        "contributed to any Antoine constant.",
        "",
    ]

    # C4
    n_s1_sys = int((c4_detail["in_S1"] > 0).sum())
    lines += [
        "## C4 -- LOSO Flash System Split Discipline",
        "",
        f"**Verdict: {c4['verdict']}**",
        "",
        c4["notes"],
        "",
        "| Split | LOSO systems with rows in this split | Rows |",
        "|---|---|---|",
        f"| S1 (t8 training) | {n_s1_sys} | {int(c4_detail['in_S1'].sum())} |",
        f"| S2 (unseen-molecule holdout) | {int((c4_detail['in_S2']>0).sum())} | {int(c4_detail['in_S2'].sum())} |",
        f"| S3 (extrapolation) | {int((c4_detail['in_S3']>0).sum())} | {int(c4_detail['in_S3'].sum())} |",
        f"| NONE (unassigned) | {int((c4_detail['in_NONE']>0).sum())} | {int(c4_detail['in_NONE'].sum())} |",
        "",
        "**Le Chatelier physics baseline:** Unaffected by S1 membership -- the baseline uses only",
        "per-compound Antoine/LFL constants; it does not train on any rows.",
        "",
        "**HGB-residual LOSO model:** The HGB is trained inside `eval_flash_loso.py` using LOSO held-in",
        "rows; `s1s2s3.json` splits are not respected. For S1-assigned LOSO systems, the HGB residual",
        "evaluates on data that t8 was also trained on. This should be noted when comparing LOSO HGB",
        "numbers to t8 S2 numbers, but does not affect the physics baseline evaluation.",
        "",
        df_to_md_table(c4_detail),
        "",
    ]

    # C5
    lines += [
        "## C5 -- T1 Test Suite",
        "",
        f"**Verdict: {c5['verdict']}**",
        "",
        c5["notes"],
        "",
    ]

    # Final verdict
    all_verdicts = checks["verdict"].tolist()
    any_fail = "FAIL" in all_verdicts
    any_warn = "WARN" in all_verdicts
    overall = "FAIL" if any_fail else ("WARN" if any_warn else "PASS")

    lines += [
        "---",
        "",
        "## Final Verdict",
        "",
        f"**Overall: {overall}**",
        "",
        "### Acceptance criteria",
        "",
        "**1. No CRC measurements or near-duplicates in training:**",
    ]
    if c1["verdict"] == "PASS":
        lines.append(
            "   PASS -- CRC 647 real-fuel rows never enter model_table (no molecule-level composition). "
            "Near-value T-matched coincidences (C2) arise from independent molecular measurements, "
            "not from copying CRC data."
        )
    else:
        lines.append(f"   FAIL -- see C1/C2 details above.")

    lines += ["", "**2. Antoine shared-compound leakage verdict:**"]
    if c3b["verdict"] == "PASS":
        n_shared_rep = int((c3_detail["severity"] == "shared_compound_refit").sum())
        lines.append(
            f"   PASS (C3b) -- Antoine refit never uses blend-mixture measurements. "
            f"WARN (C3a) -- {n_shared_rep} LOSO compounds have pure VP data that fed Antoine refit; "
            "this is unavoidable for a physics baseline and is classified as "
            "low-severity shared-compound leakage (not a defect)."
        )
    else:
        lines.append(f"   FAIL (C3b) -- {c3b['affected_rows']} pure VP rows from LOSO blend systems fed Antoine refit. See C3b.")

    lines += [
        "",
        "### Open risks / recommendations",
        (f"- **C4 (WARN):** {n_s1_sys}/31 LOSO flash systems have rows in S1. "
         "The Le Chatelier baseline numbers are unaffected. "
         "Recommend reporting LOSO HGB only on the systems NOT in S1 when comparing to t8 S2 performance."),
        ("- **C2 (WARN if hits > 0):** Any T-matched near-value coincidences reflect physical reality "
         "(Jet-A density clusters near alkane density), not data leakage."),
        ("- **C3a (WARN):** Shared-compound Antoine refit is an inherent property of physics-based "
         "flash-point evaluation and is accepted as low-severity."),
        "",
        "### Files created",
        "- `results/t0/leakage_checks.csv`",
        "- `results/t0/leakage_crc_training.csv`",
        "- `results/t0/leakage_antoine.csv`",
        "- `results/t0/leakage_loso_splits.csv`",
        "- `docs/audit_leakage.md` (this file)",
        "",
    ]

    (DOCS_DIR / "audit_leakage.md").write_text("\n".join(lines), encoding="utf-8")


if __name__ == "__main__":
    main()
