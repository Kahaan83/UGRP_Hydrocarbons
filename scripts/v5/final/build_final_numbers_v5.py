#!/usr/bin/env python
"""Build the v5 final numbers pack (tables, figures, master document) from COMMITTED outputs only.

Every source file is read with ``git show <commit>:<path>`` at a pinned commit (default ``a5afab1``), so the
working tree (including files that a running pytest may rewrite) is never read.  Nothing is recomputed: the
script selects, joins and plots values that already exist in committed CSV / JSON / YAML / Markdown files.

Outputs (all new): ``<out>/table_*.csv``, ``<out>/README.md``, ``<out>/figures/*.png`` and the master document
(``--doc``, default ``docs/final_numbers_v5.md``).  Writing anywhere under ``results/v5/BASE`` is refused.

No ``fuelblend`` import, no network.  Usage::

    .venv/Scripts/python.exe scripts/v5/final/build_final_numbers_v5.py
    .venv/Scripts/python.exe scripts/v5/final/build_final_numbers_v5.py --out <dir> --doc <file> --commit a5afab1
"""
from __future__ import annotations

import argparse
import io
import json
import re
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

REPO = Path(__file__).resolve().parents[3]
DEFAULT_COMMIT = "a3c5967"

RERUN = "results/v5/BASE/delta_adopt_rescore/rerun"
FROZEN = "results/v5/BASE"
SRC = {
    "sb_rerun": f"{RERUN}/scoreboard.csv",
    "sb_frozen": f"{FROZEN}/scoreboard.csv",
    "drift": "results/v5/BASE/delta_adopt_rescore/drift.csv",
    "rescore_readme": "results/v5/BASE/delta_adopt_rescore/README.md",
    "flash_rows_rerun": f"{RERUN}/rows_t2_flash.csv",
    "yaml": "configs/uncertainty_v5.yaml",
    "spec_astm": "configs/specs_astm.yaml",
    "spec_bis": "configs/specs_bis.yaml",
    "freeze_summary": "results/v5/BASE/delta_freeze_adopt/summary.csv",
    "freeze_systems": "results/v5/BASE/delta_freeze_adopt/per_system_S_F.csv",
    "freeze_readme": "results/v5/BASE/delta_freeze_adopt/README.md",
    "freeze_rows": "results/v5/BASE/delta_freeze_adopt/rows.csv",
    "freeze_t2": "results/v5/BASE/delta_freeze_t2_unifac/summary.csv",
    "freeze_t2_rows": "results/v5/BASE/delta_freeze_t2_unifac/rows.csv",
    "flash_lab_margin": "results/v5/BASE/delta_lab_flash_li2014/lab_margin.csv",
    "mlv": "results/v5/BASE/delta_ml_viscosity/readings.json",
    "mlv5": "results/v5/BASE/delta_ml_viscosity/readings_R-MLV-5.json",
    "mlv5_csv": "results/v5/BASE/delta_ml_viscosity/realfuel_transfer_R-MLV-5.csv",
    "mlv_rows": "results/v5/BASE/delta_ml_viscosity/rows_scored.csv",
    "mld": "results/v5/BASE/delta_ml_density/readings.json",
    "mld_rows": "results/v5/BASE/delta_ml_density/rows_scored.csv",
    "mlysi": "results/v5/BASE/delta_ml_ysi/readings.json",
    "designs_table": "results/v5/BASE/delta_adopt_designs/results_table.csv",
    "designs_csv": "results/v5/BASE/delta_adopt_designs/designs.csv",
    "designs_readme": "results/v5/BASE/delta_adopt_designs/README.md",
    "negative": "docs/v5_negative_results.md",
    "model_decision": "docs/model_decision.md",
    "claims": "docs/claims_register.md",
    "decisions": "docs/decisions_2026-11-15.md",
}

# Okabe-Ito colour-blind-safe palette
OI = {"black": "#000000", "orange": "#E69F00", "sky": "#56B4E9", "green": "#009E73", "yellow": "#F0E442",
      "blue": "#0072B2", "verm": "#D55E00", "purple": "#CC79A7", "grey": "#999999"}

HEADLINE_PROPS = ["viscosity_kin_cSt", "density_kg_m3", "flash_point_C", "freezing_point_C", "pour_point_C",
                  "smoke_point_mm", "distillation_T10_C", "distillation_T50_C", "distillation_T90_C",
                  "heat_of_combustion_net_MJ_kg"]
JOIN_KEY = ["property", "tier", "set", "model", "variant", "label", "T_K", "cells_file", "status", "evidence_path"]


# ----------------------------------------------------------------------------- git readers
class Src:
    def __init__(self, commit: str):
        self.commit = commit
        self.sha = subprocess.run(["git", "-C", str(REPO), "rev-parse", commit], capture_output=True, text=True,
                                  check=True).stdout.strip()

    def raw(self, path: str) -> bytes:
        return subprocess.run(["git", "-C", str(REPO), "show", f"{self.commit}:{path}"], capture_output=True,
                              check=True).stdout

    def text(self, path: str) -> str:
        return self.raw(path).decode("utf-8")

    def csv(self, path: str, **kw) -> pd.DataFrame:
        return pd.read_csv(io.BytesIO(self.raw(path)), low_memory=False, **kw)

    def json(self, path: str):
        return json.loads(self.text(path))

    def yaml(self, path: str):
        return yaml.safe_load(self.text(path))

    def subject(self, sha: str) -> str:
        r = subprocess.run(["git", "-C", str(REPO), "log", "-1", "--format=%h|%ad|%s", "--date=short", sha],
                           capture_output=True, text=True)
        return r.stdout.strip()


# ----------------------------------------------------------------------------- helpers
def fmt(v, nd=4):
    if v is None or (isinstance(v, float) and np.isnan(v)):
        return ""
    if isinstance(v, (bool, np.bool_)):
        return str(bool(v))
    if isinstance(v, (int, np.integer)):
        return str(int(v))
    if isinstance(v, (float, np.floating)):
        if np.isinf(v):
            return "inf" if v > 0 else "-inf"
        if v == 0:
            return "0"
        a = abs(v)
        if a >= 1000:
            return f"{v:.1f}"
        return f"{v:.{nd}g}" if a < 1e-3 else f"{float(f'{v:.{nd}g}'):g}"
    return str(v).replace("|", "/").replace("\n", " ")


def md_table(df: pd.DataFrame, cols=None, nd=4) -> str:
    cols = cols or list(df.columns)
    lines = ["| " + " | ".join(cols) + " |", "|" + "---|" * len(cols)]
    for _, r in df[cols].iterrows():
        lines.append("| " + " | ".join(fmt(r[c], nd) for c in cols) + " |")
    return "\n".join(lines)


def key_frame(df: pd.DataFrame) -> pd.DataFrame:
    k = df.copy()
    for c in JOIN_KEY:
        k[c] = k[c].astype(str) if c != "T_K" else k[c].fillna(-1.0).round(2).astype(str)
    return k


def style():
    plt.rcParams.update({"figure.facecolor": "white", "axes.facecolor": "white", "savefig.facecolor": "white",
                         "font.size": 10, "axes.titlesize": 11, "axes.labelsize": 10, "axes.grid": True,
                         "grid.color": "#dddddd", "grid.linewidth": 0.6, "axes.spines.top": False,
                         "axes.spines.right": False, "legend.frameon": False})


def save(fig, path: Path):
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)


# ----------------------------------------------------------------------------- table A
def headline_rule(sb: pd.DataFrame) -> pd.Series:
    """Mechanical selection of headline rows (declared in the README):
    designated rows of the headline properties with an aggregate temperature (no T_K split) and
      * T2_measured: sets starting with 'pooled' or equal to 'measured_pooled';
      * CRC647_r3: labels dev and fold_b;
      * T1_lab: label dev;
      * T3_ajftd: label dev;
    excluding rows with status 'no_evaluation' and the trivial comparator models ('const', 'const_mean', 'C_all').
    """
    d = (sb.variant == "designated") & sb.property.isin(HEADLINE_PROPS) & sb.T_K.isna()
    d &= (sb.status != "no_evaluation") & ~sb.model.isin(["const", "const_mean", "C_all", "constant_mean"])
    t2 = (sb.tier == "T2_measured") & (sb.set.astype(str).str.startswith("pooled") | (sb.set == "measured_pooled"))
    crc = (sb.tier == "CRC647_r3") & sb.label.isin(["dev", "fold_b"])
    t1 = (sb.tier == "T1_lab") & (sb.label == "dev")
    t3 = (sb.tier == "T3_ajftd") & (sb.label == "dev")
    return d & (t2 | crc | t1 | t3)


def build_A(S: Src):
    cur, frz = S.csv(SRC["sb_rerun"]), S.csv(SRC["sb_frozen"])
    sel = cur[headline_rule(cur)].copy()
    kc, kf = key_frame(sel), key_frame(frz)
    assert not kc.duplicated(JOIN_KEY).any(), "headline join key not unique (re-run)"
    assert not kf.duplicated(JOIN_KEY).any(), "join key not unique (frozen)"
    keep = ["n_groups", "bias", "mae", "trivial_model", "mae_trivial", "verdict", "margin90", "unsafe_rate"]
    m = kc.merge(kf[JOIN_KEY + keep], on=JOIN_KEY, how="left", suffixes=("_current", "_frozen"), validate="1:1")
    m["T_K"] = np.nan
    for c in ["bias", "mae", "verdict"]:
        a, b = m[f"{c}_current"], m[f"{c}_frozen"]
        if c == "verdict":
            m[f"{c}_changed"] = a.astype(str) != b.astype(str)
        else:
            m[f"{c}_changed"] = ~np.isclose(a.astype(float), b.astype(float), rtol=1e-6, atol=1e-9, equal_nan=True)
    m["changed"] = m[["bias_changed", "mae_changed", "verdict_changed"]].any(axis=1)
    m["source_current"] = SRC["sb_rerun"]
    m["source_frozen"] = SRC["sb_frozen"]
    order = {p: i for i, p in enumerate(HEADLINE_PROPS)}
    tier_order = {"T2_measured": 0, "CRC647_r3": 1, "T1_lab": 2, "T3_ajftd": 3}
    m["_o"] = m.property.map(order) * 10 + m.tier.map(tier_order)
    m = m.sort_values(["_o", "set", "model", "label"]).drop(columns="_o")
    cols = ["property", "tier", "set", "model", "variant", "label", "status", "n_groups_current", "bias_current",
            "mae_current", "trivial_model_current", "mae_trivial_current", "verdict_current", "margin90_current",
            "unsafe_rate_current", "n_groups_frozen", "bias_frozen", "mae_frozen", "verdict_frozen", "changed",
            "cells_file", "evidence_path", "source_current", "source_frozen"]
    A = m[cols].reset_index(drop=True)
    A.insert(A.columns.get_loc("variant") + 1, "scoreboard_variant", A["variant"])

    # superseded flash margin90 (the scoreboard carries the pre-Q2 yaml `value`, not the applied margin)
    yf = S.yaml(SRC["yaml"])["margins"]["flash_point_C.lab_mixtures.ideal"]
    A["margin90_note"] = ""
    sup = (A.property == "flash_point_C") & np.isclose(A.margin90_current.astype(float), yf["value"], rtol=1e-6)
    A.loc[sup, "margin90_note"] = (f"superseded (pre-Q2) {yf['value']:.4f}; applied {yf['applied_margin']:.4f} "
                                   f"({yf['n_groups_applied']} groups)")

    # freezing: scoreboard cells are ideal SvL, not the adopted designated route (UNIFAC SvL, CLM-73)
    isv = (A.property == "freezing_point_C") & (A.model == "ideal_svl")
    A.loc[isv, "variant"] = "ideal_svl, scoreboard cell, frozen/re-scored (not the adopted route)"
    t2 = S.csv(SRC["freeze_t2"])
    u = t2[(t2.model == "svl_unifac") & (t2.variant == "designated") & t2.set.isin(["pooled", "pooled_locked"])]
    urows = pd.DataFrame({
        "property": u.property, "tier": u.tier, "set": u.set, "model": u.model,
        "variant": "designated (CLM-73; delta_freeze_t2_unifac)", "scoreboard_variant": "(not a scoreboard cell)",
        "label": u.label, "status": u.status, "n_groups_current": u.n_groups, "bias_current": u.bias,
        "mae_current": u.mae, "trivial_model_current": u.trivial_model, "mae_trivial_current": u.mae_trivial,
        "verdict_current": u.verdict, "margin90_current": np.nan, "unsafe_rate_current": np.nan,
        "n_groups_frozen": np.nan, "bias_frozen": np.nan, "mae_frozen": np.nan, "verdict_frozen": "",
        "changed": False, "cells_file": "", "evidence_path": u.evidence_path, "source_current": SRC["freeze_t2"],
        "source_frozen": "", "margin90_note": [f"coverage at applied {m:.4f} C margin: {c:g} (uninformative, bias "
                                               f"{b:+.2f} C)" for m, c, b in zip(u.margin_used, u.coverage_margin, u.bias)]})
    A = pd.concat([A, urows], ignore_index=True)
    A["_o"] = A.property.map(order) * 10 + A.tier.map(tier_order)
    A = A.sort_values("_o", kind="mergesort").drop(columns="_o").reset_index(drop=True)

    # UNIFAC SvL freezing (designated route since CLM-73) is scored only on the lab liquidus set S_F.
    fs = S.csv(SRC["freeze_summary"])
    sf = fs[fs.scope == "S_F"].set_index("model")
    A_frz = pd.DataFrame([{
        "property": "freezing_point_C", "tier": "T1_lab", "set": "lab_liquidus_S_F", "model": mdl,
        "variant": "designated (CLM-73)" if mdl == "svl_unifac" else "variant (ideal SvL)",
        "label": "dev", "n_systems": int(sf.loc[mdl, "n_systems"]), "n_rows": int(sf.loc[mdl, "n_scored"]),
        "mae": sf.loc[mdl, "mae"], "bias": sf.loc[mdl, "bias"], "mae_constant": sf.loc[mdl, "mae_const"],
        "margin_alpha10": sf.loc[mdl, "margin_alpha10"], "frac_unsafe": sf.loc[mdl, "frac_unsafe"],
        "source": SRC["freeze_summary"]} for mdl in ["svl_unifac", "svl_verified"]])

    # Per-source T2 designated rows (supplement)
    t2s = cur[(cur.variant == "designated") & cur.property.isin(HEADLINE_PROPS) & cur.T_K.isna()
              & (cur.tier == "T2_measured") & ~headline_rule(cur)].copy()
    kt = key_frame(t2s).merge(kf[JOIN_KEY + keep], on=JOIN_KEY, how="left", suffixes=("_current", "_frozen"),
                              validate="1:1")
    A2 = kt[["property", "set", "model", "n_groups_current", "bias_current", "mae_current", "mae_trivial_current",
             "verdict_current", "bias_frozen", "mae_frozen", "verdict_frozen"]].sort_values(
        ["property", "set", "model"]).reset_index(drop=True)
    A2["source_current"], A2["source_frozen"] = SRC["sb_rerun"], SRC["sb_frozen"]
    return A, A_frz, A2, t2


# ----------------------------------------------------------------------------- table B
def build_B(S: Src):
    dr = S.csv(SRC["drift"])
    s = dr[(dr.path == SRC["sb_frozen"]) & (dr.kind == "cell")]
    B = s.groupby("cause").size().rename("n_cells").reset_index().sort_values("n_cells", ascending=False)
    B["n_designated_rows"] = B.cause.map(s.groupby("cause").designated.apply(lambda x: int(x.astype(str).eq("True").sum())))
    B["source"] = SRC["drift"]
    B = B.reset_index(drop=True)
    total = int(B.n_cells.sum())
    files = dr[dr.kind != "cell"].groupby(["path", "cause"]).size().rename("n_rows").reset_index()
    return B, total, files


# ----------------------------------------------------------------------------- table C
def spec_limits(S: Src):
    a = S.yaml(SRC["spec_astm"])["families"]["jet_a"]["constraints"]
    b = S.yaml(SRC["spec_bis"])["families"]["jet_a1"]["constraints"]

    def lim(c, k):
        d = c.get(k, {})
        out = []
        if "min" in d:
            out.append(f"min {d['min']:g}")
        if "max" in d:
            out.append(f"max {d['max']:g}")
        return " / ".join(out)

    keys = {"flash_point_C": "flash_point_C", "freezing_point_C": "freezing_point_C",
            "viscosity_kin_cSt": "kin_visc_m20C_cSt", "density_kg_m3": "density_15C_kg_m3",
            "heat_of_combustion_net_MJ_kg": "nhoc_MJ_kg", "smoke_point_mm": "smoke_point_mm"}
    names = {"viscosity_kin_cSt": " cSt at -20 C", "density_kg_m3": " kg/m3 at 15 C", "flash_point_C": " C",
             "freezing_point_C": " C", "heat_of_combustion_net_MJ_kg": " MJ/kg", "smoke_point_mm": " mm"}
    return {p: (lim(a, k) + names[p], lim(b, k) + names[p]) for p, k in keys.items()}


MARGIN_ROWS = [
    # (yaml key, role)
    ("flash_point_C.lab_mixtures.ideal", "designated (design judge)"),
    ("flash_point_C.lab_mixtures.unifac", "variant"),
    ("flash_point_C.real_fuels_D93_reconstructed.ideal", "check (floored at 0)"),
    ("flash_point_C.real_fuels_D56.ideal", "no valid margin (n=1)"),
    ("freezing_point_C.lab_liquidus.svl_unifac", "designated (CLM-73)"),
    ("freezing_point_C.lab_liquidus.svl_verified", "variant (ideal SvL)"),
    ("freezing_point_C.crc647.wat_ideal_minus_loo_offset", "CRC 647 reading"),
    ("freezing_point_C.crc647.constant_mean_loo", "comparator (no skill)"),
    ("viscosity_kin_cSt.lab_mixtures.raw_refutas.T_lt_263K", "Refutas (out of ML domain)"),
    ("viscosity_kin_cSt.lab_mixtures.raw_refutas.T_263_283K", "Refutas (out of ML domain)"),
    ("viscosity_kin_cSt.lab_mixtures.raw_refutas.T_ge_283K", "Refutas (out of ML domain)"),
    ("viscosity_kin_cSt.crc647.recalibrated_refutas.all_T", "CRC 647 reading"),
    ("viscosity_kin_cSt.crc647.structured_ridge.all_T", "CRC 647 reading"),
    ("density_kg_m3.lab_mixtures.ideal_volume.over", "designated (ideal volume)"),
    ("density_kg_m3.lab_mixtures.ideal_volume.under", "designated (ideal volume)"),
    ("density_kg_m3.crc647_ideal_volume.all_T", "CRC 647 reading"),
    ("heat_of_combustion_net_MJ_kg.crc647.group_additivity_composition_H", "designated"),
    ("smoke_point_mm.ajftd_v5_in_domain.bulk_regression", "designated in GCxGC domain (R-SC-2b)"),
    ("smoke_point_mm.lab_mixtures.bulk_regression", "lab regression domain"),
]


def build_C(S: Src, A: pd.DataFrame):
    y = S.yaml(SRC["yaml"])["margins"]
    lim = spec_limits(S)
    lm = S.csv(SRC["flash_lab_margin"])
    fsum = S.csv(SRC["freeze_summary"])
    sf = fsum[fsum.scope == "S_F"].set_index("model")
    mlv = S.json(SRC["mlv"])["R-MLV-6"]
    rows = []
    for key, role in MARGIN_ROWS:
        v = y[key]
        prop = key.split(".")[0]
        n = v.get("n_groups_applied", v.get("n_groups"))
        note, src = "", SRC["yaml"]
        if "n_groups_applied" in v:
            note = (f"applied margin derived on {v['n_groups_applied']} groups; the yaml `value` "
                    f"{v['value']:.4f} (superseded) was derived on {v['n_groups']} groups")
        if key == "flash_point_C.lab_mixtures.ideal":
            r = lm[(lm.variant == "designated") & np.isclose(lm.margin, v["applied_margin"])].iloc[0]
            assert int(r.n_groups) == int(n)
            note += f"; {int(r.n_groups)} systems / {int(r.n_points)} points, scope {r.scope}"
            src += "; " + SRC["flash_lab_margin"]
        if key.startswith("freezing_point_C.lab_liquidus."):
            mdl = key.split(".")[-1]
            r = sf.loc[mdl]
            same = np.isclose(r.margin_alpha10, v["applied_margin"])
            note = (f"stored config value derived on {v['n_groups']} groups; the {int(r.n_systems)}-system S_F reading "
                    f"(delta_freeze_adopt) gives {'the same' if same else 'a different'} {r.margin_alpha10:.4f}")
            src += "; " + SRC["freeze_summary"]
        if key == "viscosity_kin_cSt.lab_mixtures.raw_refutas.T_ge_283K":
            r6 = mlv["refutas"]["T_ge_283K"]
            note = (f"yaml config entry (lab raw Refutas, {v['n_groups']} groups); the ML delta recomputes Refutas on "
                    f"its own {r6['n_dev_systems']} dev systems (R-MLV-6): {r6['margin']:.4f}")
        rows.append({"property": prop, "margin_key": key, "role": role, "model": v.get("model"),
                     "applied_margin": v.get("applied_margin", np.nan), "value": v.get("value"),
                     "unit": v.get("unit"), "direction": v.get("direction"), "n_groups": n, "n_note": note,
                     "spec_jet_a": lim[prop][0], "spec_jet_a1": lim[prop][1], "source": src})
    for b in ["T_lt_263K", "T_263_283K", "T_ge_283K"]:
        rows.append({"property": "viscosity_kin_cSt", "margin_key": f"ml_pair_residual.{b}",
                     "role": "designated in domain (ml_pair_residual)", "model": "LightGBM pair residual on Refutas",
                     "applied_margin": mlv["ml"][b]["margin"], "value": mlv["ml"][b]["margin"], "unit": "ln-ratio",
                     "direction": "two_sided", "n_groups": mlv["ml"][b]["n_dev_systems"],
                     "n_note": (f"R-MLV-6 on the ML delta's dev systems; Refutas on the same systems "
                                f"{mlv['refutas'][b]['margin']:.4f}"),
                     "spec_jet_a": lim["viscosity_kin_cSt"][0], "spec_jet_a1": lim["viscosity_kin_cSt"][1],
                     "source": SRC["mlv"] + " R-MLV-6"})
    mld = S.json(SRC["mld"])["R-MLD-4"]["ml"]
    rows.append({"property": "density_kg_m3", "margin_key": "ml_pair_excess_volume.locked",
                 "role": "variant, known binaries only (computed on locked systems)",
                 "model": "LightGBM pair excess volume", "applied_margin": mld["locked_margin_kg_m3"],
                 "value": mld["locked_margin_kg_m3"], "unit": "kg/m3", "direction": "two_sided",
                 "n_groups": S.json(SRC["mld"])["locked_systems"],
                 "n_note": "computed on the locked systems; locked coverage is therefore not an independent check",
                 "spec_jet_a": lim["density_kg_m3"][0],
                 "spec_jet_a1": lim["density_kg_m3"][1], "source": SRC["mld"] + " R-MLD-4"})
    rows.append({"property": "density_kg_m3", "margin_key": "ml_pair_excess_volume.dev_oof",
                 "role": "variant, dev out-of-fold reading", "model": "LightGBM pair excess volume",
                 "applied_margin": np.nan, "value": mld["dev_oof_margin_kg_m3"], "unit": "kg/m3",
                 "direction": "two_sided", "n_groups": S.json(SRC["mld"])["dev_systems"],
                 "spec_jet_a": lim["density_kg_m3"][0], "spec_jet_a1": lim["density_kg_m3"][1],
                 "source": SRC["mld"] + " R-MLD-4"})
    # T2 real-fuel cold viscosity: scoreboard margin90 (inf)
    t2c = A[(A.property == "viscosity_kin_cSt") & (A.set == "pooled_T_lt_263K")].iloc[0]
    rows.append({"property": "viscosity_kin_cSt", "margin_key": "scoreboard T2 pooled_T_lt_263K refutas_raw margin90",
                 "role": "T2 real fuels, cold bin (no valid margin)", "model": "refutas_raw",
                 "applied_margin": np.nan, "value": t2c.margin90_current, "unit": "ln-ratio", "direction": "two_sided",
                 "n_groups": t2c.n_groups_current, "spec_jet_a": lim["viscosity_kin_cSt"][0],
                 "spec_jet_a1": lim["viscosity_kin_cSt"][1], "source": SRC["sb_rerun"]})
    return pd.DataFrame(rows)


# ----------------------------------------------------------------------------- table D
def build_D(S: Src):
    v, d = S.json(SRC["mlv"]), S.json(SRC["mld"])
    v5 = S.json(SRC["mlv5"])
    rows = []

    def gate(comp, g, set_, r, base, pass_, n_rows=None, note=""):
        rows.append({"component": comp, "gate": g, "set": set_, "n_systems": r.get("n_systems"), "n_rows": n_rows,
                     "baseline": base, "mae_baseline": r["mae_a"], "mae_ml": r["mae_b"], "gain": r["gain"],
                     "gain_ci95_lo": r["lo"], "gain_ci95_hi": r["hi"],
                     "ci_basis": "per-system group bootstrap, same basis as gain", "pass": pass_, "note": note})

    gate("ml_pair_residual (viscosity, ln)", "R-MLV-1", "dev (system-grouped OOF)", v["R-MLV-1"]["ml"], "Refutas",
         v["R-MLV-1"]["pass"], note=f"constant-offset comparator MAE {v['R-MLV-1']['const']['mae_b']:.4f}")
    gate("ml_pair_residual (viscosity, ln)", "R-MLV-2", "locked systems", v["R-MLV-2"]["ml"], "Refutas",
         v["R-MLV-2"]["pass"], note=f"constant-offset comparator MAE {v['R-MLV-2']['const']['mae_b']:.4f}")
    t = v["R-MLV-3"]["all"]
    rows.append({"component": "ml_pair_residual (viscosity, ln)", "gate": "R-MLV-3", "set": "ternaries+ (row-level)",
                 "n_systems": t["systems"], "n_rows": t["rows"], "baseline": "Refutas", "mae_baseline": t["mae_refutas"],
                 "mae_ml": t["mae_ml"], "gain": t["mae_refutas"] - t["mae_ml"], "gain_ci95_lo": np.nan,
                 "gain_ci95_hi": np.nan, "ci_basis": "no CI (row-level reading)", "pass": v["R-MLV-3"]["pass"],
                 "note": "row-level MAE over all rows"})
    gate("ml_pair_residual (viscosity, ln)", "R-MLV-3", "ternaries+ (per-system macro)", t["macro"], "Refutas",
         v["R-MLV-3"]["pass"], n_rows=t["rows"], note="system-macro MAE with bootstrap CI")
    r4 = v["R-MLV-4"]
    gate("ml_pair_residual (viscosity, ln)", "R-MLV-4", "cold (< 283 K)", r4, "Refutas", r4["pass"],
         n_rows=r4["rows"], note=f"under-prediction rate {r4['under_pred_rate_refutas']:.3f} -> {r4['under_pred_rate_ml']:.3f}")
    rows.append({"component": "ml_pair_residual (viscosity, ln)", "gate": "R-MLV-5 (descriptive)",
                 "set": "T2 real fuels (frozen rows_t2_bulk)", "n_systems": np.nan, "n_rows": v5["n_rows"],
                 "baseline": "Refutas", "mae_baseline": v5["mae_ln_refutas"], "mae_ml": v5["mae_ln_ml"],
                 "gain": v5["mae_ln_refutas"] - v5["mae_ln_ml"], "gain_ci95_lo": np.nan, "gain_ci95_hi": np.nan,
                 "ci_basis": "no CI (descriptive)", "pass": "descriptive (MIXED)",
                 "note": f"ML not worse in {v5['temps_ml_not_worse']} of {v5['temps_total']} source x T cells"})
    gate("ml_pair_excess_volume (density, kg/m3)", "R-MLD-1", "dev (system-grouped OOF)", d["R-MLD-1"]["ml"],
         "ideal volume", d["R-MLD-1"]["pass"])
    gate("ml_pair_excess_volume (density, kg/m3)", "R-MLD-2", "locked systems", d["R-MLD-2"]["ml"], "ideal volume",
         d["R-MLD-2"]["pass"])
    t = d["R-MLD-3"]["all"]
    rows.append({"component": "ml_pair_excess_volume (density, kg/m3)", "gate": "R-MLD-3", "set": "ternaries+ (row-level)",
                 "n_systems": t["systems"], "n_rows": t["rows"], "baseline": "ideal volume", "mae_baseline": t["mae_ideal"],
                 "mae_ml": t["mae_ml"], "gain": t["mae_ideal"] - t["mae_ml"], "gain_ci95_lo": t["macro"]["lo"],
                 "gain_ci95_hi": t["macro"]["hi"],
                 "ci_basis": (f"MIXED BASIS: gain is row-level; the CI is the per-system macro CI (macro gain "
                              f"{t['macro']['gain']:.4f})"),
                 "pass": d["R-MLD-3"]["pass"], "note": "FAILED"})
    D = pd.DataFrame(rows)
    D["source"] = np.where(D.component.str.startswith("ml_pair_residual"), SRC["mlv"], SRC["mld"])
    D.loc[D.gate.str.startswith("R-MLV-5"), "source"] = SRC["mlv5"]
    D2 = []
    ym = S.yaml(SRC["yaml"])["margins"]
    for b in ["T_lt_263K", "T_263_283K", "T_ge_283K"]:
        m = v["R-MLV-6"]
        yv = ym[f"viscosity_kin_cSt.lab_mixtures.raw_refutas.{b}"]
        D2.append({"component": "viscosity", "T_bin": b, "margin_ml": m["ml"][b]["margin"],
                   "margin_physics_baseline": m["refutas"][b]["margin"],
                   "physics_baseline": f"Refutas recomputed on the ML delta's {m['refutas'][b]['n_dev_systems']} dev systems (R-MLV-6)",
                   "applied_physics_margin_yaml": f"{yv['applied_margin']:.4f} ({yv['n_groups']} groups, yaml raw_refutas.{b})",
                   "n_dev_systems": m["ml"][b]["n_dev_systems"],
                   "locked_systems": m["locked_coverage"]["ml"][b]["locked_systems"],
                   "locked_covered_ml": m["locked_coverage"]["ml"][b]["covered"],
                   "locked_covered_physics": m["locked_coverage"]["refutas"][b]["covered"],
                   "source": SRC["mlv"] + "; " + SRC["yaml"]})
    m = d["R-MLD-4"]
    yo, yu = ym["density_kg_m3.lab_mixtures.ideal_volume.over"], ym["density_kg_m3.lab_mixtures.ideal_volume.under"]
    D2.append({"component": "density", "T_bin": "all", "margin_ml": m["ml"]["locked_margin_kg_m3"],
               "margin_physics_baseline": m["ideal"]["locked_margin_kg_m3"],
               "physics_baseline": f"ideal volume, computed on the {d['locked_systems']} locked systems (R-MLD-4)",
               "applied_physics_margin_yaml": (f"{yo['applied_margin']:.4f} over / {yu['applied_margin']:.4f} under "
                                               f"({yo['n_groups']} systems, yaml ideal_volume)"),
               "n_dev_systems": d["dev_systems"],
               "locked_systems": m["ml"]["locked_coverage_by_locked_margin"][1],
               "locked_covered_ml": m["ml"]["locked_coverage_by_locked_margin"][0],
               "locked_covered_physics": m["ideal"]["locked_coverage_by_locked_margin"][0],
               "source": SRC["mld"] + "; " + SRC["yaml"]})
    return D, pd.DataFrame(D2), v, d


# ----------------------------------------------------------------------------- table E
def build_E(S: Src):
    r = S.csv(SRC["designs_table"])
    g = S.csv(SRC["designs_csv"])
    cols = ["pool", "spec", "certified", "n_pass", "n_fail", "n_uncertified", "min_slack", "flash_pred_C",
            "flash_margin_C", "flash_bound_C", "flash_limit_C", "freezing_unifac_C", "freezing_unifac_margin_C",
            "freezing_unifac_bound_C", "freezing_limit_C", "freezing_ideal_svl_bound_C",
            "ideal_svl_bound_minus_limit_C", "engler_T50_T10_C", "engler_T90_T10_C", "rayleigh_T50_T10_C",
            "rayleigh_T90_T10_C", "dicyclo_wt_pct", "labels", "failing"]
    E = r[cols].copy()
    E["source"] = SRC["designs_table"]
    label_counts = g["freezing_model_y"].value_counts().to_dict() if "freezing_model_y" in g else {}
    route_counts = g["freezing_route"].value_counts().to_dict() if "freezing_route" in g else {}
    return E, label_counts, route_counts


# ----------------------------------------------------------------------------- table F / G
def parse_md_table(text: str, header_start: str):
    lines = text.splitlines()
    i = next(k for k, l in enumerate(lines) if l.startswith(header_start))
    hdr = [c.strip() for c in lines[i].strip("|").split("|")]
    rows = []
    for l in lines[i + 2:]:
        if not l.startswith("|"):
            break
        cells = [c.strip() for c in l.strip().strip("|").split("|")]
        rows.append(dict(zip(hdr, cells)))
    return pd.DataFrame(rows)


def build_F(S: Src, A: pd.DataFrame, D: pd.DataFrame):
    neg = parse_md_table(S.text(SRC["negative"]), "| Workstream")
    F = pd.DataFrame({"item": neg["Workstream"], "what": neg["What was tried"],
                      "outcome": neg["Outcome (declared gate, not moved)"], "where": neg["Where"],
                      "source": SRC["negative"]})
    md = S.text(SRC["model_decision"])
    extra = []
    for pat, item in [(r"^- ML cetane number.*$", "ML DCN"), (r"^- ML YSI.*$", "ML YSI")]:
        m = re.search(pat, md, flags=re.M)
        extra.append({"item": item, "what": "", "outcome": m.group(0)[2:] if m else "(line not found)",
                      "where": SRC["model_decision"] + " (v5 ML components)", "source": SRC["model_decision"]})
    ysi = S.json(SRC["mlysi"])["R-MLY-1"]
    extra.append({"item": "ML YSI (numbers)", "what": "ridge on molecular descriptors, outer CV",
                  "outcome": (f"median relative error ridge {ysi['median_rel_err']['ridge']:.3f} (gate <= 0.20); "
                              f"MAE ridge {ysi['mae']['ridge']:.1f} vs class-mean {ysi['mae']['classmean']:.1f}, "
                              f"gain CI [{ysi['gain_vs']['classmean']['lo95']:.1f}, {ysi['gain_vs']['classmean']['hi95']:.1f}]"),
                  "where": SRC["mlysi"], "source": SRC["mlysi"]})
    for _, r in A[(A.property == "freezing_point_C") & (A.model == "ideal_svl")
                  & A.tier.isin(["T2_measured", "CRC647_r3"])].iterrows():
        extra.append({"item": f"freezing, ideal SvL on real fuels ({r.tier} {r.set}, {r.label})",
                      "what": f"{r.tier} {r.set} ({r.label})",
                      "outcome": f"MAE {fmt(r.mae_current)} C vs trivial {fmt(r.mae_trivial_current)}; verdict {r.verdict_current}",
                      "where": SRC["sb_rerun"], "source": SRC["sb_rerun"]})
    for _, r in A[(A.property == "freezing_point_C") & (A.model == "svl_unifac")].iterrows():
        extra.append({"item": f"freezing, UNIFAC SvL (designated) on real fuels ({r.tier} {r.set}, {r.label})",
                      "what": "designated route on T2, delta_freeze_t2_unifac",
                      "outcome": (f"MAE {fmt(r.mae_current)} C vs trivial {fmt(r.mae_trivial_current)} "
                                  f"(n {int(r.n_groups_current)}), bias {fmt(r.bias_current)} C; verdict {r.verdict_current}"),
                      "where": SRC["freeze_t2"], "source": SRC["freeze_t2"]})
    t2r = S.csv(SRC["freeze_t2_rows"]).set_index("entity")
    a325 = t2r.loc["AJFTD:325"]
    sbn = A[(A.property == "freezing_point_C") & (A.model == "ideal_svl") & (A.tier == "T2_measured")
            & (A.set == "pooled")].n_groups_current
    un = A[(A.property == "freezing_point_C") & (A.model == "svl_unifac") & (A.set == "pooled")].n_groups_current
    extra.append({"item": "freezing caveat: AJFTD:325 under the designated constants",
                  "what": "UNIFAC SvL designated constants on AJFTD:325 (12 components)",
                  "outcome": (f"predicted {a325.unifac_desig:.1f} C (UNIFAC; ideal {a325.ideal_desig:.1f} C) vs measured "
                              f"{a325.freezing_C:.0f} C, error {a325.unifac_desig - a325.freezing_C:+.1f} C; the fuel has no "
                              f"scoreboard prediction (I_scoreboard empty), hence the designated pooled n = "
                              f"{int(un.iloc[0]) if len(un) else '?'} vs the scoreboard's "
                              f"{int(sbn.iloc[0]) if len(sbn) else '?'}"),
                  "where": SRC["freeze_t2_rows"], "source": SRC["freeze_t2_rows"]})
    for _, r in A[(A.property == "viscosity_kin_cSt") & (A.set == "pooled_T_lt_263K")].iterrows():
        extra.append({"item": "cold viscosity on real fuels (< 263 K)", "what": f"{r.model}, T2 pooled",
                      "outcome": f"MAE {fmt(r.mae_current)} ln, unsafe rate {fmt(r.unsafe_rate_current)}, verdict "
                                 f"{r.verdict_current}, margin90 {fmt(r.margin90_current)} (no valid real-fuel margin; "
                                 "the lab cold-bin margins are finite)",
                      "where": SRC["sb_rerun"], "source": SRC["sb_rerun"]})
    r3 = D[D.gate == "R-MLD-3"].iloc[0]
    extra.append({"item": "ML density on ternaries", "what": "ml_pair_excess_volume R-MLD-3",
                  "outcome": f"FAILED: MAE {fmt(r3.mae_baseline)} (ideal) -> {fmt(r3.mae_ml)} kg/m3 (ML)",
                  "where": SRC["mld"], "source": SRC["mld"]})
    return pd.concat([F, pd.DataFrame(extra)], ignore_index=True)


def build_G(S: Src, E: pd.DataFrame, D2: pd.DataFrame, label_counts):
    y = S.yaml(SRC["yaml"])["margins"]
    fsys = S.csv(SRC["freeze_systems"]).set_index("group")
    frd = S.text(SRC["freeze_readme"])
    drd = S.text(SRC["designs_readme"])
    m_rise = re.search(r"Its predictions rise by\s+([\d.]+)-([\d.]+) C", frd)
    m_swap = re.search(r"\(C\) MAE ([\d.]+) \(ideal\) / ([\d.]+) \(UNIFAC\)", frd)
    m_max = re.search(r"max ([\d.]+) %", drd)
    bch = fsys.loc["ffd8e7cf5721"]
    xyl = fsys.loc["74aede3c6579"]
    cert = E[E.certified == "Y"]
    sb = S.csv(SRC["sb_rerun"])
    crc_fb = sb[(sb.tier == "CRC647_r3") & (sb.label == "fold_b") & (sb.variant == "designated")].n_groups.dropna().max()
    fr = S.csv(SRC["freeze_rows"])
    fr = fr[(fr.group == "ffd8e7cf5721") & (fr.in_SF == True)]  # noqa: E712
    res_i, res_u = fr.pred_verified - fr.measured, fr.pred_unifac - fr.measured
    i_max = res_u.idxmax()
    G = [
        {"item": "single GC fuel for D56 flash", "value": y["flash_point_C.real_fuels_D56.ideal"]["n_groups"],
         "unit": "fuels", "detail": y["flash_point_C.real_fuels_D56.ideal"]["reason"], "source": SRC["yaml"]},
        {"item": "thin locked sets: ML viscosity cold bins", "value": int(D2.loc[0, "locked_systems"]),
         "unit": "locked systems (< 263 K)",
         "detail": f"263-283 K: {int(D2.loc[1, 'locked_systems'])} locked systems; ML margin in that bin equals Refutas",
         "source": SRC["mlv"]},
        {"item": "thin locked sets: CRC 647 fold_b", "value": int(crc_fb), "unit": "fuels",
         "detail": "n_groups of the CRC647_r3 fold_b scoreboard rows", "source": SRC["sb_rerun"]},
        {"item": "bicyclohexyl dHfus 13.94 kJ/mol over-predicts freezing (safe side)", "value": bch["mae_unifac"],
         "unit": "C system MAE (UNIFAC)",
         "detail": (f"all {len(fr)} S_F rows over-predict; largest over-prediction +{res_i.max():.1f} C (ideal) / "
                    f"+{res_u.max():.1f} C (UNIFAC) at T = {fr.loc[i_max, 'T_K']:.2f} K; smallest over-prediction "
                    f"+{bch['worst_unifac']:.2f} C (UNIFAC; the `worst_unifac` column of per_system_S_F.csv); predictions rise by "
                    f"{m_rise.group(1) if m_rise else '?'}-{m_rise.group(2) if m_rise else '?'} C vs the old constant"),
         "source": SRC["freeze_systems"] + "; " + SRC["freeze_rows"] + "; " + SRC["freeze_readme"]},
        {"item": "o-/m-xylene eutectic composition likely transposed (ingest_calorimetry_primary.py EUT), NOT fixed",
         "value": xyl["worst_unifac"], "unit": "C residual on that row",
         "detail": (f"with the swap MAE would be {m_swap.group(1) if m_swap else '?'} (ideal) / "
                    f"{m_swap.group(2) if m_swap else '?'} (UNIFAC); margins unchanged"),
         "source": SRC["freeze_systems"] + "; " + SRC["freeze_readme"]},
        {"item": "certified Jet A-1 design outside AJFTD dicycloparaffin range",
         "value": float(cert[cert.spec == "bis:jet_a1"].dicyclo_wt_pct.iloc[0]), "unit": "wt % dicyclo",
         "detail": f"AJFTD GCxGC maximum {m_max.group(1) if m_max else '?'} %", "source": SRC["designs_table"]},
        {"item": "engine freezing label", "value": int(sum(label_counts.values())), "unit": "design rows",
         "detail": "designs.csv freezing_model_y = " + ", ".join(f"{k} ({v})" for k, v in label_counts.items())
                   + "; value and margin are UNIFAC SvL (freezing_route column)", "source": SRC["designs_csv"]},
        {"item": "no lab validation", "value": 0, "unit": "designs measured",
         "detail": "lab test ruled out by the user (2026-11-14); certification is model-level only",
         "source": SRC["negative"]},
        {"item": "re-score pending adoption (U4)", "value": np.nan, "unit": "",
         "detail": "current-default numbers live in delta_adopt_rescore/rerun; frozen BASE and SHA256SUMS_v5_base.txt unchanged",
         "source": SRC["rescore_readme"]},
    ]
    return pd.DataFrame(G)


def build_H(S: Src):
    shas = ["a5afab1", "28dfe8d", "adb5972", "fc618b5", "457f236", "1b196e3", "9fbff24", "5986a6b", "5a62e20"]
    rows = []
    for s in shas:
        sub = S.subject(s)
        if sub:
            h, dte, msg = sub.split("|", 2)
            rows.append({"commit": h, "date": dte, "subject": msg})
    return pd.DataFrame(rows)


def p1p3_wording(E: pd.DataFrame) -> str:
    q = E[E.pool.isin(["P1", "P2", "P3"])]
    sp = lambda x: x.replace("astm:", "").replace("bis:", "")
    fail = q[q.n_fail > 0]
    rest = q[q.n_fail == 0]
    parts = [f"{r.pool} {sp(r.spec)} fails {int(r.n_fail)} constraint{'s' if r.n_fail != 1 else ''} "
             f"({str(r.failing).replace(',', ', ')})" for r in fail.itertuples()]
    un = sorted(set(int(x) for x in rest.n_uncertified))
    parts.append(f"the other {len(rest)} P1-P3 rows fail none but leave "
                 f"{' / '.join(str(x) for x in un)} constraints uncertified (no declared margins)")
    return "P1-P3 are not certified: " + "; ".join(parts)


def dicyclo_max(S: Src):
    m = re.search(r"max ([\d.]+) %", S.text(SRC["designs_readme"]))
    return float(m.group(1)) if m else float("nan")


# ----------------------------------------------------------------------------- figures
def fig_flash(S: Src, A: pd.DataFrame, margin: float, out: Path):
    t = S.csv(SRC["flash_rows_rerun"])
    t = t[(t.chosen == True) & (t.status == "ok")].copy()  # noqa: E712
    t["err"] = t.pred_designated - t.measured
    # consistency: per-source bias and MAE must reproduce the re-run scoreboard (designated ideal_le_chatelier)
    sb = S.csv(SRC["sb_rerun"])
    checks = {}
    for src in ["ajftd_gcxgc", "av3322", "fuellib"]:
        r = sb[(sb.property == "flash_point_C") & (sb.tier == "T2_measured") & (sb.set == src)
               & (sb.model == "ideal_le_chatelier") & (sb.variant == "designated")].iloc[0]
        e = t[t.source == src].err
        assert len(e) == int(r.n_groups) and np.isclose(e.mean(), r.bias, atol=1e-5) and \
            np.isclose(e.abs().mean(), r.mae, atol=1e-5), f"flash rows do not reproduce scoreboard for {src}"
        checks[src] = (len(e), e.mean(), e.abs().mean())
    style()
    fig, ax = plt.subplots(figsize=(5.6, 5.2))
    lo, hi = min(t.measured.min(), t.pred_designated.min()) - 3, max(t.measured.max(), t.pred_designated.max()) + 3
    xs = np.array([lo, hi])
    ax.fill_between(xs, xs, xs + margin, color=OI["sky"], alpha=0.25, lw=0,
                    label=f"over-prediction within applied margin ({margin:.2f} C)")
    ax.plot(xs, xs, color=OI["black"], lw=1, label="parity")
    ax.plot(xs, xs + margin, color=OI["blue"], lw=1, ls="--", label="unsafe beyond this line")
    names = {"ajftd_gcxgc": ("AJFTD GCxGC", OI["orange"], "o"), "av3322": ("CRC AV-33-22", OI["green"], "s"),
             "fuellib": ("FuelLib", OI["purple"], "^")}
    for src, (nm, c, mk) in names.items():
        s = t[t.source == src]
        n, b, mae = checks[src]
        ax.scatter(s.measured, s.pred_designated, s=26, c=c, marker=mk, edgecolor="white", lw=0.5,
                   label=f"{nm} (n={n}, MAE {mae:.2f}, bias {b:+.2f} C)", zorder=3)
    pooled = A[(A.property == "flash_point_C") & (A.set == "pooled") & (A.model == "ideal_le_chatelier")].iloc[0]
    ax.set_xlim(lo, hi); ax.set_ylim(lo, hi); ax.set_aspect("equal")
    ax.set_xlabel("Measured flash point (C)"); ax.set_ylabel("Predicted flash point (C), designated route")
    ax.set_title(f"Flash point, T2 real fuels (pooled n={int(pooled.n_groups_current)}, "
                 f"MAE {pooled.mae_current:.2f} C)")
    ax.legend(fontsize=7.5, loc="upper left")
    save(fig, out / "a_flash_parity_t2.png")
    return checks


def fig_visc(S: Src, out: Path):
    r = S.csv(SRC["mlv_rows"])
    t5 = S.csv(SRC["mlv5_csv"])
    style()
    fig, axs = plt.subplots(1, 2, figsize=(10.5, 4.4), gridspec_kw={"width_ratios": [1.25, 1]})
    ax = axs[0]
    bins = ["T_lt_263K", "T_263_283K", "T_ge_283K"]
    labels = {"T_lt_263K": "< 263 K", "T_263_283K": "263-283 K", "T_ge_283K": ">= 283 K"}
    pos, data, cols = [], [], []
    for i, b in enumerate(bins):
        s = r[r.t_bin == b]
        for j, (c, col) in enumerate([("resid_refutas", OI["grey"]), ("resid_ml", OI["blue"])]):
            pos.append(i * 3 + j); data.append(s[c].abs().dropna().values); cols.append(col)
    bp = ax.boxplot(data, positions=pos, widths=0.8, showfliers=False, patch_artist=True)
    for p, c in zip(bp["boxes"], cols):
        p.set_facecolor(c); p.set_alpha(0.75)
    for med in bp["medians"]:
        med.set_color("black")
    ax.set_xticks([i * 3 + 0.5 for i in range(3)])
    ax.set_xticklabels([f"{labels[b]}\n(n rows={int((r.t_bin == b).sum())})" for b in bins])
    ax.set_ylabel("|ln(pred / measured)|")
    ax.set_title("Lab mixtures: Refutas vs ML pair residual")
    ax.legend([bp["boxes"][0], bp["boxes"][1]], ["Refutas (physics)", "ML pair residual"], loc="upper right")
    ax = axs[1]
    mk = {"ajftd_gcxgc": ("o", OI["orange"], "AJFTD GCxGC"), "av3322": ("s", OI["green"], "CRC AV-33-22"),
          "fuellib": ("^", OI["purple"], "FuelLib")}
    mx = max(t5.mae_ln_refutas.max(), t5.mae_ln_ml.max()) * 1.1
    ax.plot([0, mx], [0, mx], color="black", lw=1)
    for s, (m, c, nm) in mk.items():
        q = t5[t5.source == s]
        ax.scatter(q.mae_ln_refutas, q.mae_ln_ml, marker=m, c=c, s=18 + 2 * q.n_groups, edgecolor="white", lw=0.5,
                   label=nm, zorder=3)
    ax.set_xlim(0, mx); ax.set_ylim(0, mx)
    ax.set_xlabel("MAE ln, Refutas"); ax.set_ylabel("MAE ln, ML pair residual")
    ax.set_title("Real fuels (R-MLV-5, descriptive): per source x T")
    ax.text(0.97 * mx, 0.05 * mx, "below the line = ML better", ha="right", fontsize=8)
    ax.legend(fontsize=8, loc="upper left")
    save(fig, out / "b_viscosity_refutas_vs_ml.png")


def fig_density(S: Src, out: Path):
    r = S.csv(SRC["mld_rows"])
    r["grp"] = np.where(r.n_comp >= 3, "ternary+", "binary") + " / " + r.role
    style()
    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    groups = ["binary / dev", "binary / locked", "ternary+ / dev", "ternary+ / locked"]
    pos, data, cols = [], [], []
    for i, g in enumerate(groups):
        s = r[r.grp == g]
        for j, (c, col) in enumerate([("err_ideal", OI["grey"]), ("err_ml", OI["verm"])]):
            pos.append(i * 3 + j); data.append(s[c].abs().dropna().values); cols.append(col)
    bp = ax.boxplot(data, positions=pos, widths=0.8, showfliers=False, patch_artist=True)
    for p, c in zip(bp["boxes"], cols):
        p.set_facecolor(c); p.set_alpha(0.75)
    for med in bp["medians"]:
        med.set_color("black")
    ax.set_xticks([i * 3 + 0.5 for i in range(len(groups))])
    ax.set_xticklabels([f"{g}\n(n rows={int((r.grp == g).sum())})" for g in groups], fontsize=8.5)
    ax.set_ylabel("|predicted - measured| (kg/m3)")
    ax.set_title("Lab mixture density: ideal volume vs ML excess volume (ternaries FAIL)")
    ax.legend([bp["boxes"][0], bp["boxes"][1]], ["ideal volume (physics)", "ML pair excess volume"], loc="upper left")
    save(fig, out / "c_density_ideal_vs_ml.png")


def fig_margins(C: pd.DataFrame, out: Path):
    style()
    panels = [("flash_point_C", "Flash point (C, over-side)"), ("freezing_point_C", "Freezing point (C, under-side)"),
              ("viscosity_kin_cSt", "Kinematic viscosity (ln ratio)"), ("density_kg_m3", "Density (kg/m3)")]
    fig, axs = plt.subplots(2, 2, figsize=(11, 7.5))
    for ax, (p, ttl) in zip(axs.ravel(), panels):
        q = C[C.property == p].copy()
        q["shown"] = q.applied_margin.where(q.applied_margin.notna(), q.value).astype(float)
        fin = q[np.isfinite(q.shown)]
        inf = q[~np.isfinite(q.shown)]
        y = np.arange(len(fin))
        col = [OI["blue"] if "designated" in r else OI["grey"] for r in fin.role]
        ax.barh(y, fin.shown, color=col, alpha=0.85)
        for yi, v in zip(y, fin.shown):
            ax.text(v, yi, f" {v:.3g}", va="center", fontsize=7.5)
        lab = [k.replace(p + ".", "").replace("lab_mixtures.", "lab.") for k in fin.margin_key]
        ax.set_yticks(y); ax.set_yticklabels(lab, fontsize=7.2); ax.invert_yaxis()
        spec = q.spec_jet_a.iloc[0]; spec1 = q.spec_jet_a1.iloc[0]
        extra = "" if inf.empty else f"\nno valid margin (inf): {len(inf)} entr{'y' if len(inf) == 1 else 'ies'}"
        sp = f"Jet A and Jet A-1: {spec}" if spec == spec1 else f"Jet A {spec}; Jet A-1 {spec1}"
        ax.set_title(f"{ttl}\n{sp}{extra}", fontsize=9)
        ax.set_xlim(0, fin.shown.max() * 1.3 if len(fin) else 1)
    fig.suptitle("90 % one-/two-sided margins (blue = designated / applied route)", fontsize=11)
    save(fig, out / "d_margins_vs_spec.png")


def fig_ml_bars(D: pd.DataFrame, out: Path):
    style()
    q = D.copy()
    fig, axs = plt.subplots(1, 2, figsize=(11, 4.4), gridspec_kw={"width_ratios": [1.6, 1]})
    for ax, comp, unit in [(axs[0], "ml_pair_residual", "MAE (ln)"), (axs[1], "ml_pair_excess_volume", "MAE (kg/m3)")]:
        s = q[q.component.str.startswith(comp)].reset_index(drop=True)
        x = np.arange(len(s))
        ax.bar(x - 0.2, s.mae_baseline, 0.4, color=OI["grey"], label="physics baseline")
        ax.bar(x + 0.2, s.mae_ml, 0.4, color=OI["blue"] if comp.endswith("residual") else OI["verm"], label="ML")
        for xi, (a, b, p) in enumerate(zip(s.mae_baseline, s.mae_ml, s["pass"])):
            tag = "PASS" if p is True or p == "True" else ("FAIL" if p is False or p == "False" else "descr.")
            ax.text(xi, max(a, b) * 1.04, tag, ha="center", fontsize=7.5,
                    color=OI["green"] if tag == "PASS" else (OI["verm"] if tag == "FAIL" else OI["black"]))
        ax.set_xticks(x)
        short = {"dev (system-grouped OOF)": "dev", "locked systems": "locked", "ternaries+ (row-level)": "tern. rows",
                 "ternaries+ (per-system macro)": "tern. macro", "cold (< 283 K)": "cold",
                 "T2 real fuels (frozen rows_t2_bulk)": "real fuels"}
        ax.set_xticklabels([f"{g.replace(' (descriptive)', '')}\n{short.get(st, st)}" for g, st in zip(s.gate, s.set)],
                           fontsize=8)
        ax.set_ylabel(unit)
        ax.set_ylim(0, max(s.mae_baseline.max(), s.mae_ml.max()) * 1.18)
        ax.set_title("Viscosity: ml_pair_residual" if comp.endswith("residual") else "Density: ml_pair_excess_volume")
        ax.legend(fontsize=8, loc="upper left")
    save(fig, out / "e_ml_before_after.png")


def fig_designs(E: pd.DataFrame, C: pd.DataFrame, out: Path):
    q = E[E.flash_bound_C.notna()].copy()
    fmg = float(q.flash_margin_C.iloc[0])
    flim = float(q.flash_limit_C.iloc[0])
    umg = float(q.freezing_unifac_margin_C.iloc[0])
    img = float(C.loc[C.margin_key == "freezing_point_C.lab_liquidus.svl_verified", "applied_margin"].iloc[0])
    q["name"] = q.pool + " " + q.spec.str.replace("astm:", "").str.replace("bis:", "")
    q["flash_head"] = q.flash_bound_C - q.flash_limit_C
    q["frz_unifac_head"] = q.freezing_limit_C - q.freezing_unifac_bound_C
    q["frz_ideal_head"] = -q.ideal_svl_bound_minus_limit_C
    style()
    fig, axs = plt.subplots(1, 2, figsize=(11, 5.4), sharey=True)
    y = np.arange(len(q))
    cert = (q.certified == "Y").values
    for ax, cols, ttl in [(axs[0], [("flash_head", OI["blue"], f"flash: bound - {flim:g} C limit")],
                           f"Flash (bound = prediction - {fmg:.2f} C margin)"),
                          (axs[1], [("frz_unifac_head", OI["blue"], f"UNIFAC SvL (designated, +{umg:.2f} C)"),
                                    ("frz_ideal_head", OI["grey"], f"ideal SvL variant (+{img:.2f} C)")],
                           "Freezing (limit - safe bound)")]:
        w = 0.8 / len(cols)
        for k, (c, col, lab) in enumerate(cols):
            ax.barh(y + (k - (len(cols) - 1) / 2) * w, q[c], w, color=col, alpha=0.85, label=lab)
        ax.axvline(0, color="black", lw=1)
        ax.set_title(ttl, fontsize=10)
        ax.set_xlabel("headroom (C); > 0 passes this constraint")
        ax.legend(fontsize=7.5, loc="upper center", bbox_to_anchor=(0.5, -0.16), ncol=2)
    axs[0].set_yticks(y)
    axs[0].set_yticklabels([("[certified] " if c else "") + n for n, c in zip(q.name, cert)], fontsize=8)
    for t, c in zip(axs[0].get_yticklabels(), cert):
        if c:
            t.set_fontweight("bold")
    axs[0].invert_yaxis()
    fig.suptitle("Designs with declared margins (P0 family; P1-P3 have no declared margins and are not shown)",
                 fontsize=10)
    save(fig, out / "f_design_headroom.png")


def fig_drift(B: pd.DataFrame, total: int, out: Path):
    style()
    fig, ax = plt.subplots(figsize=(7.5, 3.4))
    q = B.sort_values("n_cells")
    ax.barh(q.cause, q.n_cells, color=[OI["orange"] if c.startswith("D3") else OI["sky"] for c in q.cause])
    for i, v in enumerate(q.n_cells):
        ax.text(v, i, f" {v}", va="center", fontsize=8)
    ax.set_xlabel("changed scoreboard cells")
    ax.set_title(f"Re-score drift by cause ({total} cells; D1/D4/D5/D7 and unattributed: 0)", fontsize=10)
    ax.set_xlim(0, q.n_cells.max() * 1.15)
    ax.tick_params(axis="y", labelsize=8)
    save(fig, out / "g_rescore_drift_by_cause.png")


# ----------------------------------------------------------------------------- document
def write_doc(path: Path, S: Src, A, A_frz, A2, t2, B, total, C, D, D2, E, F, G, H, label_counts, flash_checks, out):
    rel = lambda p: Path(p).resolve().relative_to(REPO).as_posix() if Path(p).resolve().is_relative_to(REPO) else str(p)
    o = rel(out)
    ncert = int((E.certified == "Y").sum())
    certs = ", ".join(f"{p} {s}" for p, s in E[E.certified == "Y"][["pool", "spec"]].itertuples(index=False))
    fl = A[(A.property == "flash_point_C") & (A.set == "pooled") & (A.model == "ideal_le_chatelier")].iloc[0]
    aj = A2[(A2.property == "flash_point_C") & (A2.set == "ajftd_gcxgc") & (A2.model == "ideal_le_chatelier")].iloc[0]
    fm = C[C.margin_key == "flash_point_C.lab_mixtures.ideal"].iloc[0]
    fz_u = C[C.margin_key == "freezing_point_C.lab_liquidus.svl_unifac"].iloc[0]
    fz_i = C[C.margin_key == "freezing_point_C.lab_liquidus.svl_verified"].iloc[0]
    ta = ["property", "tier", "set", "model", "variant", "label", "status", "n_groups_current", "bias_current",
          "mae_current", "mae_trivial_current", "verdict_current", "margin90_current", "margin90_note", "bias_frozen",
          "mae_frozen", "verdict_frozen"]
    yy = S.yaml(SRC["yaml"])["margins"]
    yfl = yy["flash_point_C.lab_mixtures.ideal"]
    fu = A[(A.property == "freezing_point_C") & (A.model == "svl_unifac")].set_index("set")
    tu = t2[(t2.model == "svl_unifac") & (t2.variant == "designated")].set_index("set")
    fsf = A_frz.set_index("model").loc["svl_unifac"]
    bias_lo, bias_hi = fu.bias_current.min(), fu.bias_current.max()
    n_sb_frz = int(A[(A.property == "freezing_point_C") & (A.model == "ideal_svl") & (A.tier == "T2_measured")
                     & (A.set == "pooled")].n_groups_current.iloc[0])
    v1 = D2[D2.T_bin == "T_ge_283K"].iloc[0]
    cert = E[E.certified == "Y"]
    umg = float(cert.freezing_unifac_margin_C.iloc[0])
    all_ideal_exceed = bool((cert.ideal_svl_bound_minus_limit_C > 0).all())
    pyb = E[(E.pool == "P0y-b") & (E.spec == "bis:jet_a1")].iloc[0]
    dmax = dicyclo_max(S)
    s = []
    s.append("# v5 final numbers (report and presentation pack)\n")
    s.append(f"Generated by `scripts/v5/final/build_final_numbers_v5.py` from files committed at `{S.sha[:7]}` "
             f"(read with `git show`, never from the working tree). Every number below is copied from a committed file; "
             f"the CSV behind each table is in `{o}/` and `{o}/README.md` maps each table to its sources. No gate, margin, "
             "spec limit or verdict was changed; verdict strings are the scoreboard's own.\n")
    s.append("**Status of the numbers.** \"Current defaults\" = the BASE re-score under the decisions of 2026-11-15 "
             "(`results/v5/BASE/delta_adopt_rescore/rerun/`). Adopting that re-score as the frozen BASE is user question "
             "**U4** and has NOT happened: the frozen BASE (`results/v5/BASE/`, commit `5a62e20`) and "
             "`SHA256SUMS_v5_base.txt` are unchanged, and are shown next to every current value.\n")
    s.append("**Certification** means: passes the design judge under the model and the declared margins. No design has "
             "been measured and no lab test is possible (user decision 2026-11-14).\n")

    s.append("## A. Headline results per property\n")
    s.append("Selection rule (mechanical, in the script): designated scoreboard rows of the headline properties with no "
             "temperature split, for T2 real fuels (pooled sets), CRC 647 r3 (dev and fold_b), T1 lab (dev) and T3 AJFTD "
             "bulk models (dev); rows with status `no_evaluation` and the trivial comparator models are excluded "
             "(the trivial MAE is a column). Error = predicted - measured; viscosity in ln(pred/meas). "
             "`mae_trivial` is the scoreboard's trivial comparator (leave-one-group-out constant). Tier: T1 = lab "
             "pure/mixture data, T2 = real fuels with measured composition, T3 = AJFTD bulk-property regressions, "
             "CRC647_r3 = real fuels with reconstructed composition.\n")
    s.append(f"**Freezing on real fuels, in plain words.** Neither freezing route beats the trivial constant on real fuels. "
             f"The designated route (UNIFAC SvL, CLM-73) reads `{fu.loc['pooled', 'verdict_current']}` on the pooled T2 set "
             f"(n {int(fu.loc['pooled', 'n_groups_current'])}, MAE {fmt(fu.loc['pooled', 'mae_current'])} C vs constant "
             f"{fmt(fu.loc['pooled', 'mae_trivial_current'])} C) and `{fu.loc['pooled_locked', 'verdict_current']}` on the "
             f"locked set (n {int(fu.loc['pooled_locked', 'n_groups_current'])}, MAE "
             f"{fmt(fu.loc['pooled_locked', 'mae_current'])} C vs {fmt(fu.loc['pooled_locked', 'mae_trivial_current'])} C). "
             f"UNIFAC SvL is validated only on the {int(fsf.n_systems)}-system lab liquidus set S_F. The real-fuel coverage "
             f"of its {fmt(tu.loc['pooled', 'margin_used'])} C margin ({tu.loc['pooled', 'coverage_margin']:g} pooled / "
             f"{tu.loc['pooled_locked', 'coverage_margin']:g} locked) is uninformative, because the bias is "
             f"{bias_lo:+.1f} to {bias_hi:+.1f} C (over-prediction), so a one-sided margin covers it almost "
             f"automatically. The UNIFAC rows are read from `{SRC['freeze_t2']}` (not scoreboard cells); the "
             f"`ideal_svl` rows are the scoreboard cells. The pooled UNIFAC n is "
             f"{int(fu.loc['pooled', 'n_groups_current'])} vs the scoreboard's {n_sb_frz}, because AJFTD:325 is scored "
             "with the designated constants (see section F).\n")
    s.append(md_table(A, ta) + "\n")
    s.append(f"Flash `margin90_current` is the scoreboard's pre-Q2 figure ({fmt(yfl['value'])} C, superseded); the applied "
             f"flash margin is {fmt(yfl['applied_margin'])} C on {yfl['n_groups_applied']} groups (section C).\n")
    s.append(f"Table: `{o}/table_A_headline.csv` (full columns incl. margin90, unsafe rate, cells_file, evidence_path). "
             f"Per-source T2 rows: `{o}/table_A2_t2_per_source.csv`.\n")
    s.append("**Two route notes (read before quoting).**\n")
    s.append("* **Freezing.** The adopted designated freezing route is UNIFAC Schroeder-van Laar (CLM-73), but every BASE "
             "freezing cell aggregates committed CSVs, so the scoreboard freezing rows above are still the ideal SvL "
             "(`ideal_svl`) route. UNIFAC SvL is scored on real fuels (T2) only in `delta_freeze_t2_unifac` (the "
             "`svl_unifac` rows of Table A, no skill shown); its gated evidence is the lab liquidus set S_F:\n")
    s.append(md_table(A_frz, ["set", "model", "variant", "n_systems", "n_rows", "mae", "bias", "mae_constant",
                              "margin_alpha10", "frac_unsafe"]) + "\n")
    s.append("* **Viscosity.** The scoreboard's designated T2 viscosity route is `refutas_raw` (v4 routing), and its "
             "cold bins are `worse_than_trivial`. `ml_pair_residual`, the designated in-domain viscosity of "
             "`predict_fuel`, has gate evidence on lab mixtures only (section D); its real-fuel transfer is descriptive "
             "and mixed.\n")
    s.append("* **D86.** All distillation cells are `no_evaluation` (exploratory) in the scoreboard; no D86 accuracy is "
             "claimed here.\n")

    s.append("## B. What changed versus the frozen BASE, and why\n")
    s.append(f"{total} scoreboard cells changed between the frozen BASE and the re-score. Cause per cell from "
             "`drift.csv` (no cell traces to decisions D1, D4, D5 or D7, and none is unattributed):\n")
    s.append(md_table(B, ["cause", "n_cells", "n_designated_rows"]) + "\n")
    s.append("`D3_overlay` = decision 3, flash-registry additions (bicyclohexyl, cis-decalin Antoine + LFL) on by default. "
             "`PRE_api_route` = code that changed before the decisions (`predict()` moved onto the BASE route in "
             "delta_api_route); it moves the `api_as_shipped` rows. `PRE_flash_margin_q2` = the "
             f"{fmt(fm.applied_margin)} C flash margin now reported for the API rows.\n")
    s.append(f"Headline example: flash T2 AJFTD GCxGC designated ({int(aj.n_groups_current)} fuels) bias {fmt(aj.bias_frozen)} -> "
             f"{fmt(aj.bias_current)} C, MAE {fmt(aj.mae_frozen)} -> {fmt(aj.mae_current)} C, verdict "
             f"{aj.verdict_frozen} -> {aj.verdict_current} (D3). Changed headline rows:\n")
    ch = A[A.changed]
    s.append(md_table(ch, ["property", "set", "model", "label", "bias_frozen", "bias_current", "mae_frozen",
                           "mae_current", "verdict_frozen", "verdict_current"]) + "\n")
    s.append(f"Table: `{o}/table_B_drift_by_cause.csv`; figure `{o}/figures/g_rescore_drift_by_cause.png`.\n")

    s.append("## C. Margins and the spec limits they certify against\n")
    s.append("`applied_margin` is the value the design judge / API applies (`configs/uncertainty_v5.yaml`, read only; "
             "`value` is the raw conformal quantile). ML margins are not in the yaml; they are read from the ML delta "
             "`readings.json`. The `uncertainty_margin` fields inside the spec yaml files are NOT applied by the judge "
             "and are not shown.\n")
    s.append(md_table(C, ["property", "margin_key", "role", "applied_margin", "value", "unit", "direction", "n_groups",
                          "n_note", "spec_jet_a", "spec_jet_a1"]) + "\n")
    s.append(f"Notes: the flash margin applied is {fmt(fm.applied_margin)} C (the yaml `value` "
             f"{fmt(fm.value)} C and the scoreboard `margin90` {fmt(fl.margin90_current)} C are the pre-Q2 figure). Freezing: UNIFAC "
             f"{fmt(fz_u.applied_margin)} C vs ideal {fmt(fz_i.applied_margin)} C. The -20 C viscosity spec maps to the "
             "< 263 K bin. The 263-283 K ML margin equals the Refutas margin (no gain in that bin). The ML density margin "
             "was computed on the locked systems (dev out-of-fold reading shown alongside), so its locked coverage is not "
             f"an independent check. Table `{o}/table_C_margins.csv`; figure `{o}/figures/d_margins_vs_spec.png`.\n")

    s.append("## D. Machine learning: components whose declared gates passed\n")
    s.append("Two ML components passed their pre-declared gates and ship: `ml_pair_residual` (LightGBM residual on "
             "Refutas, pairwise, system-grouped) for viscosity, designated in domain by user decision; "
             "`ml_pair_excess_volume` for density, a labelled variant (dev and locked passed, ternaries FAILED, "
             "`promotion_offered` false). Gains are baseline MAE - ML MAE with a group-bootstrap 95 % CI.\n")
    s.append(md_table(D, ["component", "gate", "set", "n_systems", "n_rows", "baseline", "mae_baseline", "mae_ml",
                          "gain", "gain_ci95_lo", "gain_ci95_hi", "ci_basis", "pass", "note"]) + "\n")
    s.append("Margins and locked coverage (90 %, two-sided, ln for viscosity, kg/m3 for density):\n")
    s.append(md_table(D2, ["component", "T_bin", "margin_ml", "margin_physics_baseline", "physics_baseline",
                           "applied_physics_margin_yaml", "n_dev_systems", "locked_systems", "locked_covered_ml",
                           "locked_covered_physics"]) + "\n")
    s.append("`margin_physics_baseline` is computed on the same systems as `margin_ml` (Refutas on the ML delta's dev "
             "systems; ideal volume on the density delta's locked systems). `applied_physics_margin_yaml` is the margin "
             "the yaml applies to the physics route, derived on its own system set; each carries its own n.\n")
    s.append("Reading: R-MLV-3 passed; the gate is the per-system macro row (bootstrap CI). The row-level figure is shown "
             "too. R-MLD-3 mixes bases: its gain is row-level while its CI is the per-system macro CI (`ci_basis`). "
             f"At >= 283 K the ML margin is tighter ({v1.margin_ml:.3f} vs {v1.margin_physics_baseline:.3f} ln) but covers "
             f"{int(v1.locked_covered_ml)} of {int(v1.locked_systems)} locked systems against "
             f"{int(v1.locked_covered_physics)} of {int(v1.locked_systems)} for Refutas. R-MLV-5 (real fuels, descriptive) is mixed: a large gain at 233/253 K, worse at several "
             "warm cells; no claim beyond lab mixtures is made. Sources: "
             "`results/v5/BASE/delta_ml_viscosity/readings.json`, `readings_R-MLV-5.json`, "
             "`results/v5/BASE/delta_ml_density/readings.json`; `docs/model_decision.md` (v5 ML components). Figures "
             f"`{o}/figures/e_ml_before_after.png`, `b_viscosity_refutas_vs_ml.png`, `c_density_ideal_vs_ml.png` "
             "(residual distributions; the row files store residuals only, so no parity plot).\n")

    s.append("## E. Inverse-design results (16 designs: 8 pools x 2 specs)\n")
    s.append(f"Judge-certified under the model: **{ncert} of {len(E)}** ({certs}). None is lab-certified. All {ncert} "
             "carry `d86_model_dependent`, `freezing_route_dependent` and `overlay_dependent`, so none is certified under "
             f"every reported route. The freezing bound uses UNIFAC SvL + {umg:.2f} C; the ideal SvL bound "
             f"(+{fz_i.applied_margin:.2f} C) exceeds the limit for "
             f"{'all ' + str(ncert) if all_ideal_exceed else 'some'} of them. P0y-b Jet A-1 contains "
             f"{pyb.dicyclo_wt_pct:.1f} wt % dicycloparaffin, above the AJFTD GCxGC maximum of {dmax:g} %. The declared "
             "expectation that P0y-b Jet A-1 fails freezing is refuted for the designated route (true only on the ideal "
             f"SvL variant). {p1p3_wording(E)}.\n")
    s.append(md_table(E, ["pool", "spec", "certified", "n_pass", "n_fail", "n_uncertified", "min_slack",
                          "flash_bound_C", "freezing_unifac_bound_C", "freezing_limit_C",
                          "ideal_svl_bound_minus_limit_C", "rayleigh_T50_T10_C", "rayleigh_T90_T10_C",
                          "dicyclo_wt_pct", "labels"]) + "\n")
    lab = ", ".join(f"`{k}` ({v} rows)" for k, v in label_counts.items())
    s.append(f"Label caveat: the engine's freezing constraint label still reads {lab} in `designs.csv` "
             "(`freezing_model_y`), although the value and margin are UNIFAC SvL (`freezing_route`). A label issue only. "
             f"Table `{o}/table_E_designs.csv`; figure `{o}/figures/f_design_headroom.png`.\n")

    s.append("## F. Negative results (kept, not hidden)\n")
    s.append(md_table(F, ["item", "outcome", "where"]) + "\n")
    s.append(f"Table `{o}/table_F_negative_results.csv`.\n")

    s.append("## G. Limitations and honesty box\n")
    s.append(md_table(G, ["item", "value", "unit", "detail", "source"]) + "\n")
    s.append("Further: compositions of CRC 647 fuels are reconstructed, not measured; the T2 cells are descriptive "
             "(every T2 fuel has been read); lab flash and freezing margins rest on mostly binary systems.\n")

    s.append("## H. Reproduction\n")
    s.append("* Re-score (current defaults): `.venv/Scripts/python.exe scripts/v5/base/reproduce_base.py --full` "
             "(writes into `results/v5/BASE/`; see `results/v5/BASE/delta_adopt_rescore/README.md` for the copy/restore "
             "protocol). Frozen check: `reproduce_base.py` (default).\n")
    s.append("* Designs: `results/v5/BASE/delta_adopt_designs/README.md` (runner committed in `fc618b5`).\n")
    s.append(f"* This pack: `.venv/Scripts/python.exe scripts/v5/final/build_final_numbers_v5.py --commit {S.sha[:7]}`; "
             "test `tests/test_final_numbers_v5.py`.\n")
    s.append(md_table(H, ["commit", "date", "subject"]) + "\n")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(s), encoding="utf-8", newline="\n")


def write_readme(out: Path, S: Src, figs_skipped: list[str], fm: float = float('nan')):
    fm_txt = f"{fm:.2f}"
    t = f"""# v5 final numbers pack: tables and figures

Built by `scripts/v5/final/build_final_numbers_v5.py` from files committed at `{S.sha[:7]}` (`git show`, never the working tree).
Master document: `docs/final_numbers_v5.md`. Nothing here is recomputed; values are selected, joined and plotted.

| file | content | sources |
|---|---|---|
| `table_A_headline.csv` | headline designated rows (current defaults = re-score, plus frozen BASE side by side; join on the full key incl. cells_file/status/evidence_path); plus the UNIFAC SvL designated freezing rows on T2 pooled / pooled_locked (not scoreboard cells; `scoreboard_variant` = `(not a scoreboard cell)`) | `{SRC['sb_rerun']}`, `{SRC['sb_frozen']}`, `{SRC['freeze_t2']}` |
| `table_A_freezing_unifac_lab.csv` | UNIFAC vs ideal SvL on lab liquidus S_F | `{SRC['freeze_summary']}` |
| `table_A2_t2_per_source.csv` | T2 per-source designated rows | same as table A |
| `table_B_drift_by_cause.csv` | changed scoreboard cells by cause | `{SRC['drift']}` |
| `table_B2_drift_other_files.csv` | non-cell drift rows (rows_* summaries) by file and cause | `{SRC['drift']}` |
| `table_C_margins.csv` | applied margins + spec limits | `{SRC['yaml']}`, `{SRC['spec_astm']}`, `{SRC['spec_bis']}`, `{SRC['mlv']}`, `{SRC['mld']}`, `{SRC['sb_rerun']}` |
| `table_D_ml_gates.csv` | ML gate readings (R-MLV-1..5, R-MLD-1..3) | `{SRC['mlv']}`, `{SRC['mlv5']}`, `{SRC['mld']}` |
| `table_D2_ml_margins.csv` | ML margins and locked coverage | `{SRC['mlv']}` (R-MLV-6), `{SRC['mld']}` (R-MLD-4) |
| `table_E_designs.csv` | 16 designs | `{SRC['designs_table']}` |
| `table_F_negative_results.csv` | negative results | `{SRC['negative']}`, `{SRC['model_decision']}`, `{SRC['mlysi']}`, `{SRC['mld']}`, `{SRC['sb_rerun']}`, `{SRC['freeze_t2']}`, `{SRC['freeze_t2_rows']}` |
| `table_G_limitations.csv` | limitations | `{SRC['yaml']}`, `{SRC['freeze_systems']}`, `{SRC['freeze_rows']}`, `{SRC['freeze_readme']}`, `{SRC['designs_table']}`, `{SRC['designs_csv']}` |
| `table_H_commits.csv` | commits for reproduction | git log |

Headline selection rule: designated rows of the headline properties, no temperature split, T2 sets starting with `pooled`
(or `measured_pooled`), CRC647_r3 labels dev and fold_b, T1_lab label dev, T3_ajftd label dev excluding the trivial models.

## Figures (`figures/`, 200 dpi, Okabe-Ito palette)

| figure | content | source |
|---|---|---|
| `a_flash_parity_t2.png` | flash parity, T2 real fuels (70 chosen rows), designated `pred_designated`; one-sided over-prediction band of the applied margin ({fm_txt} C). The script asserts the rows reproduce the scoreboard per source. | `{SRC['flash_rows_rerun']}` |
| `b_viscosity_refutas_vs_ml.png` | lab |ln resid| by T bin, Refutas vs ML; real-fuel R-MLV-5 MAE per source x T | `{SRC['mlv_rows']}`, `{SRC['mlv5_csv']}` |
| `c_density_ideal_vs_ml.png` | lab |error| ideal vs ML by binary/ternary and dev/locked | `{SRC['mld_rows']}` |
| `d_margins_vs_spec.png` | margins per property with the Jet A / Jet A-1 limits in the panel title | table C |
| `e_ml_before_after.png` | ML gate MAEs, baseline vs ML, with PASS/FAIL | table D |
| `f_design_headroom.png` | flash and freezing headroom of the 10 designs with declared margins | table E |
| `g_rescore_drift_by_cause.png` | changed scoreboard cells by cause | table B |

Deviations from the brief:
""" + "\n".join(f"* {x}" for x in figs_skipped) + "\n"
    (out / "README.md").write_text(t, encoding="utf-8", newline="\n")


# ----------------------------------------------------------------------------- main
def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", default=str(REPO / "results/v5/final"))
    ap.add_argument("--doc", default=str(REPO / "docs/final_numbers_v5.md"))
    ap.add_argument("--commit", default=DEFAULT_COMMIT)
    a = ap.parse_args(argv)
    out = Path(a.out).resolve()
    base = (REPO / "results/v5/BASE").resolve()
    if out == base or base in out.parents or Path(a.doc).resolve().is_relative_to(base):
        raise SystemExit("refusing to write under results/v5/BASE")
    (out / "figures").mkdir(parents=True, exist_ok=True)
    S = Src(a.commit)

    A, A_frz, A2, t2 = build_A(S)
    B, total, B2 = build_B(S)
    C = build_C(S, A)
    D, D2, _, _ = build_D(S)
    E, label_counts, route_counts = build_E(S)
    F = build_F(S, A, D)
    G = build_G(S, E, D2, label_counts)
    H = build_H(S)
    for name, df in [("table_A_headline", A), ("table_A_freezing_unifac_lab", A_frz), ("table_A2_t2_per_source", A2),
                     ("table_B_drift_by_cause", B), ("table_B2_drift_other_files", B2), ("table_C_margins", C),
                     ("table_D_ml_gates", D), ("table_D2_ml_margins", D2), ("table_E_designs", E),
                     ("table_F_negative_results", F), ("table_G_limitations", G), ("table_H_commits", H)]:
        df.to_csv(out / f"{name}.csv", index=False, lineterminator="\n")

    fm = float(C.loc[C.margin_key == "flash_point_C.lab_mixtures.ideal", "applied_margin"].iloc[0])
    flash_checks = fig_flash(S, A, fm, out / "figures")
    fig_visc(S, out / "figures")
    fig_density(S, out / "figures")
    fig_margins(C, out / "figures")
    fig_ml_bars(D, out / "figures")
    fig_designs(E, C, out / "figures")
    fig_drift(B, total, out / "figures")
    notes = [
        "(a) the brief asked for a +-margin band; the flash margin is one-sided (direction `over`), so a one-sided band is drawn.",
        "(b), (c) the brief asked for parity plots; `rows_scored.csv` of both ML deltas stores residuals only (no measured or "
        "predicted values), so residual distributions are plotted instead. No values were rebuilt by joining other files.",
        "(f) drawn as headroom bars (bound vs limit) for the 10 designs that have declared margins, instead of a heat map.",
    ]
    write_readme(out, S, notes, fm)
    write_doc(Path(a.doc), S, A, A_frz, A2, t2, B, total, C, D, D2, E, F, G, H, label_counts, flash_checks, out)
    print(f"wrote {len(list(out.glob('*.csv')))} CSVs, {len(list((out / 'figures').glob('*.png')))} figures to {out}; doc {a.doc}")


if __name__ == "__main__":
    main()
