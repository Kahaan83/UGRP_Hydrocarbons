"""S4b: heat of combustion on measured (calorimeter) rows only, with a within-method constant (v5).

Scores the readings declared in results/v5/S4b/README.md. Nothing in the declared readings is changed here.
Reuses S4's table, features, split file, models and leave-source-group-out (cell) procedure from
scripts/v5/eval_bulk_descriptor.py. Step 0 reproduces S4's published heat numbers with C_all; the script
stops if it cannot.

    .venv/Scripts/python.exe scripts/v5/s4b_heat_measured.py
"""
from __future__ import annotations

import importlib.util
import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "v5" / "S4b"
S4DIR = ROOT / "results" / "v5" / "S4"

_spec = importlib.util.spec_from_file_location("s4_eval", ROOT / "scripts" / "v5" / "eval_bulk_descriptor.py")
s4 = importlib.util.module_from_spec(_spec)
sys.modules["s4_eval"] = s4
_spec.loader.exec_module(s4)

HEAT = "heat_of_combustion_net_MJ_kg"
BOOT_B, BOOT_SEED, MIN_GROUPS = 5000, 0, 8
sys.path.insert(0, str(ROOT / "src"))
from fuelblend.eval.group_stats import CALC, CLASS_EVIDENCE, MEAS, UNC, classify_method, norm_code  # noqa: E402,F401  (moved, bodies unchanged)


# --------------------------------------------------------------------------- labels
def build_heat_labels(lg: pd.DataFrame, dup_map: dict) -> pd.DataFrame:
    """Per de-duplicated fuel: y_meas / y_calc (median of the per-raw-fuel class medians), n rows, class set."""
    h = lg[lg["property"] == HEAT].copy()
    h["rep"] = h["fuel_id"].map(lambda i: dup_map.get(int(i), int(i)))
    h["cls"] = h["method"].map(classify_method)
    per = h.groupby(["rep", "fuel_id", "cls"])["value"].median().reset_index()
    agg = per.groupby(["rep", "cls"])["value"].median().unstack("cls")
    for c in (MEAS, CALC, UNC):
        if c not in agg:
            agg[c] = np.nan
    out = pd.DataFrame({"fuel_id": agg.index, "y_meas": agg[MEAS].to_numpy(), "y_calc": agg[CALC].to_numpy()})
    out["fuel_class"] = out["fuel_id"].map(h.groupby("rep")["cls"].agg(lambda s: "+".join(sorted(set(s)))))
    return out


# --------------------------------------------------------------------------- statistics
def boot_pos(err_ref, err_model, groups, B: int = BOOT_B, seed: int = BOOT_SEED):
    """Paired bootstrap over groups of pooled MAE(ref) - MAE(model); positive = model better. 95% percentile."""
    d = np.asarray(err_ref, float) - np.asarray(err_model, float)
    g = pd.Series(groups).astype(str).to_numpy()
    ug, inv = np.unique(g, return_inverse=True)
    s = np.bincount(inv, weights=d)
    n = np.bincount(inv).astype(float)
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, len(ug), size=(B, len(ug)))
    stat = s[idx].sum(1) / n[idx].sum(1)
    return float(d.mean()), float(np.percentile(stat, 2.5)), float(np.percentile(stat, 97.5)), len(ug)


def within_const_logo(df: pd.DataFrame, ycol: str, class_col: str, group_col: str = "cell"):
    """C_within per row: mean ycol over training rows (scored group excluded) of the same class.
    Returns (pred array, n_fallback) where fallback = training mean when no same-class training row."""
    pred = np.full(len(df), np.nan)
    g, c, y = df[group_col].to_numpy(), df[class_col].to_numpy(), df[ycol].to_numpy(float)
    nfb = 0
    for grp in np.unique(g):
        te = g == grp
        tr = ~te
        if not tr.any():
            continue
        for k in np.unique(c[te]):
            m = te & (c == k)
            same = tr & (c == k)
            if same.any():
                pred[m] = y[same].mean()
            else:
                pred[m] = y[tr].mean()
                nfb += int(m.sum())
    return pred, nfb


def run_logo(df: pd.DataFrame, ycol: str, class_col: str, all_df: pd.DataFrame, reading: str) -> pd.DataFrame:
    """Leave-cell-out on df (S4's procedure). Models fitted on df rows of the training cells. C_all = mean of
    S4's y over ALL-method training rows of `all_df` (scored cell excluded)."""
    d = df.copy()
    d["y"] = d[ycol]
    cw, _ = within_const_logo(d, "y", class_col)
    ridge, hgb, call = (np.full(len(d), np.nan) for _ in range(3))
    g = d["cell"].to_numpy()
    for grp in np.unique(g):
        te = g == grp
        if (~te).sum() < 1:
            continue
        ridge[te] = s4.fitpred("ridge", "heat", d[~te], d[te])
        hgb[te] = s4.fitpred("hgb", "heat", d[~te], d[te])
        call[te] = all_df.loc[all_df["cell"] != grp, "y"].mean()
    rows = []
    for name, p in (("C_all", call), ("C_within", cw), ("ridge", ridge), ("hgb", hgb)):
        r = d[["fuel_id", "cell", "source_group", "y"]].copy()
        r["pred"], r["model"], r["reading"], r["split"] = p, name, reading, "dev_logo"
        rows.append(r)
    return pd.concat(rows, ignore_index=True)


def score(preds: pd.DataFrame, reading: str, note: str = "") -> list[dict]:
    P = preds[preds["pred"].notna()].copy()
    P["ae"] = (P["pred"] - P["y"]).abs()
    ref = P[P["model"] == "C_within"].set_index("fuel_id")
    call = P[P["model"] == "C_all"].set_index("fuel_id")
    rows = []
    for model in ("C_all", "ridge", "hgb"):
        M = P[P["model"] == model].set_index("fuel_id")
        idx = M.index.intersection(ref.index)
        M, R = M.loc[idx], ref.loc[idx]
        if len(idx) == 0:
            continue
        diff, lo, hi, ng = boot_pos(R["ae"].to_numpy(), M["ae"].to_numpy(), M["cell"].to_numpy())
        ok = ng >= MIN_GROUPS
        rows.append({"reading": reading, "model": model, "n_fuels": len(idx), "n_groups": ng,
                     "mae_model": M["ae"].mean(), "mae_C_within": R["ae"].mean(),
                     "mae_C_all": call.loc[call.index.intersection(idx), "ae"].mean() if len(call) else np.nan,
                     "diff_Cwithin_minus_model": diff, "ci_lo": lo, "ci_hi": hi,
                     "evaluable": ok,
                     "skill": bool(ok and lo > 0) if model != "C_all" else np.nan,
                     "note": note if ok else f"not evaluable: {ng} groups (< {MIN_GROUPS}). {note}"})
    return rows


# --------------------------------------------------------------------------- stages
def reproduce(t, roles, out: Path):
    """Step 0: reproduce S4's published heat numbers (dev LOGO and locked) with C_all."""
    dv = pd.read_csv(S4DIR / "dev_scores.csv")
    lk = pd.read_csv(S4DIR / "locked_scores.csv")
    d = s4.eligible(t, "heat", roles, "dev")
    logo = s4.add_flags(s4.logo_preds(d, "heat"), t)
    sc = s4.summarize(logo, {}, "all")
    te = s4.eligible(t, "heat", roles, "locked")
    lp = s4.add_flags(s4.collect(te, "heat", "locked", {m: s4.fitpred(m, "heat", d, te) for m in s4.MODELS["heat"]}), t)
    sl = s4.summarize(lp, {}, "all")
    rows = []
    for scheme, mine, pub in (("dev_logo", sc, dv[(dv.scheme == "logo") & (dv["subset"] == "all")]),
                              ("locked", sl, lk[lk["subset"] == "all"])):
        for model in ("const", "ridge", "hgb"):
            a = mine[mine.model == model].iloc[0]
            b = pub[(pub.target == "heat") & (pub.model == model)].iloc[0]
            for q in ("n_fuels", "n_groups", "mae", "ci_lo", "ci_hi"):
                rows.append({"scheme": scheme, "model": model if model != "const" else "const (=C_all)", "quantity": q,
                             "published": b[q], "reproduced": a[q], "abs_diff": abs(float(a[q]) - float(b[q])),
                             "ok": bool(abs(float(a[q]) - float(b[q])) < 5e-4)})
    rc = pd.DataFrame(rows)
    rc.to_csv(out / "s4b_reproduction_check.csv", index=False)
    return rc, logo, d


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    sp, t, roles = s4.load_all(S4DIR)
    lg = pd.read_csv(s4.LONG)
    dup = {int(k): int(v) for k, v in sp["dup_map"].items()}

    # ---- step 0: reproduction
    rc, logo, dev_all = reproduce(t, roles, OUT)
    print("reproduction check: all ok =", bool(rc["ok"].all()))
    if not rc["ok"].all():
        print(rc[~rc["ok"]].to_string(index=False))
        print("CANNOT REPRODUCE S4 heat numbers; stopping. See s4b_reproduction_check.csv")
        return 1

    # ---- method classes (evidence per code)
    h = lg[lg["property"] == HEAT].copy()
    h["code"] = h["method"].map(norm_code)
    h["rep"] = h["fuel_id"].map(lambda i: dup.get(int(i), int(i)))
    mc = h.groupby("code").agg(n_rows=("value", "size"), n_fuels=("rep", "nunique"),
                               raw_variants=("method", lambda s: "|".join(sorted(set(map(str, s)))))).reset_index()
    mc["class"] = mc["code"].map(lambda c: CLASS_EVIDENCE.get(c, (UNC, ""))[0])
    mc["evidence"] = mc["code"].map(lambda c: CLASS_EVIDENCE.get(c, (UNC, "no documentation found"))[1])
    mc = mc[["code", "raw_variants", "class", "n_rows", "n_fuels", "evidence"]].sort_values("n_rows", ascending=False)
    mc.to_csv(OUT / "method_classes.csv", index=False)

    # ---- labels
    lab = build_heat_labels(lg, dup)
    t2 = t.merge(lab, on="fuel_id", how="left")
    assert t2[HEAT].notna().sum() == t[HEAT].notna().sum()
    dev = s4.eligible(t2, "heat", roles, "dev")        # S4's dev heat rows (y = S4 target)
    lock = s4.eligible(t2, "heat", roles, "locked")
    parts, rows = [], []

    # R-S4b-1 / 2: dev LOGO on one class
    for reading, ycol in (("R-S4b-1_measured", "y_meas"), ("R-S4b-2_calculated", "y_calc")):
        sub = dev[dev[ycol].notna()].copy()
        sub["cls"] = reading
        pr = run_logo(sub, ycol, "cls", dev, reading)
        parts.append(pr)
        rows += score(pr, reading)

    # R-S4b-3: S4 as published (mixed), C_within by fuel class set vs C_all
    cw, nfb = within_const_logo(dev, "y", "fuel_class")
    pr3 = []
    for name in ("const", "ridge", "hgb"):
        r = logo[logo["model"] == name][["fuel_id", "cell", "source_group", "y", "pred"]].copy()
        r["model"] = "C_all" if name == "const" else name
        pr3.append(r)
    r = dev[["fuel_id", "cell", "source_group", "y"]].copy()
    r["pred"], r["model"] = cw, "C_within"
    pr3.append(r)
    pr3 = pd.concat(pr3, ignore_index=True)
    pr3["reading"], pr3["split"] = "R-S4b-3_mixed", "dev_logo"
    parts.append(pr3)
    r3 = score(pr3, "R-S4b-3_mixed",
               f"C_within = mean of training fuels with the same set of row classes; {nfb} fuel-scorings fell back to the training mean")
    rows += r3

    # R-S4b-4: locked, measured only (second read, not a clean test)
    trm = dev[dev["y_meas"].notna()].copy()
    trm["y"] = trm["y_meas"]
    tem = lock[lock["y_meas"].notna()].copy()
    tem["y"] = tem["y_meas"]
    lk = []
    if len(tem) and len(trm):
        cwv = trm["y"].mean()
        for name, p in (("C_all", np.full(len(tem), dev["y"].mean())), ("C_within", np.full(len(tem), cwv)),
                        ("ridge", s4.fitpred("ridge", "heat", trm, tem)), ("hgb", s4.fitpred("hgb", "heat", trm, tem))):
            q = tem[["fuel_id", "cell", "source_group", "y"]].copy()
            q["pred"], q["model"], q["reading"], q["split"] = p, name, "R-S4b-4_locked_measured", "locked"
            lk.append(q)
        lk = pd.concat(lk, ignore_index=True)
        parts.append(lk)
        rows += score(lk, "R-S4b-4_locked_measured", "NOT A CLEAN TEST: locked rows were already scored once in S4")
        for row in rows:
            if row["reading"] == "R-S4b-4_locked_measured":
                row["evaluable"] = row["n_groups"] >= 2
                row["skill"] = np.nan      # no group threshold or skill claim declared for the locked second read
                row["note"] = ("NOT A CLEAN TEST: locked rows were already scored once in S4; "
                               f"{row['n_groups']} locked groups, no skill claim rests on it")
    pred = pd.concat(parts, ignore_index=True)
    pred.to_csv(OUT / "s4b_predictions.csv", index=False)
    mt = pd.DataFrame(rows)
    mt.to_csv(OUT / "s4b_metrics.csv", index=False)
    mt[["reading", "model", "n_fuels", "n_groups", "diff_Cwithin_minus_model", "ci_lo", "ci_hi", "evaluable", "skill"]] \
        .to_csv(OUT / "s4b_bootstrap.csv", index=False)

    # ---- readings
    rd = []

    def get(reading, model):
        s = mt[(mt.reading == reading) & (mt.model == model)]
        return s.iloc[0] if len(s) else None

    for reading, label in (("R-S4b-1_measured", "R-S4b-1"), ("R-S4b-2_calculated", "R-S4b-2")):
        for m in ("ridge", "hgb"):
            x = get(reading, m)
            if x is None:
                rd.append({"reading": label, "model": m, "outcome": "not evaluable: no rows"})
            elif not x["evaluable"]:
                rd.append({"reading": label, "model": m, "n_fuels": x["n_fuels"], "n_groups": x["n_groups"],
                           "outcome": f"not evaluable: {x['n_groups']} groups (< {MIN_GROUPS})"})
            else:
                rd.append({"reading": label, "model": m, "n_fuels": x["n_fuels"], "n_groups": x["n_groups"],
                           "mae_model": x["mae_model"], "mae_C_within": x["mae_C_within"],
                           "diff": x["diff_Cwithin_minus_model"], "ci_lo": x["ci_lo"], "ci_hi": x["ci_hi"],
                           "outcome": "skill" if x["skill"] else "no skill shown"})
    for m in ("C_all", "ridge", "hgb"):
        x = get("R-S4b-3_mixed", m)
        rd.append({"reading": "R-S4b-3", "model": m, "n_fuels": x["n_fuels"], "n_groups": x["n_groups"],
                   "mae_model": x["mae_model"], "mae_C_within": x["mae_C_within"], "diff": x["diff_Cwithin_minus_model"],
                   "ci_lo": x["ci_lo"], "ci_hi": x["ci_hi"],
                   "outcome": ("C_all vs C_within: positive diff = C_within is the better constant" if m == "C_all"
                               else "model vs C_within")})
    for m in ("C_all", "ridge", "hgb"):
        x = get("R-S4b-4_locked_measured", m)
        if x is None:
            rd.append({"reading": "R-S4b-4", "model": m, "outcome": "not evaluable: no locked measured rows"})
        else:
            rd.append({"reading": "R-S4b-4", "model": m, "n_fuels": x["n_fuels"], "n_groups": x["n_groups"],
                       "mae_model": x["mae_model"], "mae_C_within": x["mae_C_within"],
                       "diff": x["diff_Cwithin_minus_model"], "ci_lo": x["ci_lo"], "ci_hi": x["ci_hi"],
                       "outcome": "not a clean test (already scored once in S4)"})
    ev = [get("R-S4b-1_measured", m) for m in ("ridge", "hgb")]
    ev = [x for x in ev if x is not None]
    if not ev or not any(x["evaluable"] for x in ev):
        verdict = "not evaluable"
    elif any(x["skill"] for x in ev if x["evaluable"]):
        verdict = "heat skill survives on measured values"
    else:
        verdict = "not supported"
    rd.append({"reading": "R-S4b-5", "outcome": verdict})
    rd.append({"reading": "R-S4b-6", "outcome": "nothing written to configs/; S4 flash correlations remain unverified "
                                                "(coefficients quoted from memory), not re-run"})
    pd.DataFrame(rd).to_csv(OUT / "s4b_readings.csv", index=False)

    pd.set_option("display.width", 250)
    print(mc[["code", "class", "n_rows", "n_fuels"]].to_string(index=False))
    print(mt.drop(columns=["note"]).round(4).to_string(index=False))
    print("R-S4b-5:", verdict)
    return 0


if __name__ == "__main__":
    sys.exit(main())
