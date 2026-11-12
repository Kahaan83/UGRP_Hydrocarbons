"""S4: large held-out evaluation of bulk-descriptor models on AJFTD (v5).

Sub-commands (binding order): ``splits`` -> (write README.md by hand) -> ``dev`` -> ``locked`` (exactly once) -> ``report``.

    .venv/Scripts/python.exe scripts/v5/eval_bulk_descriptor.py splits
    .venv/Scripts/python.exe scripts/v5/eval_bulk_descriptor.py dev
    .venv/Scripts/python.exe scripts/v5/eval_bulk_descriptor.py locked
    .venv/Scripts/python.exe scripts/v5/eval_bulk_descriptor.py report

Hyper-parameters are fixed in this file (no tuning anywhere); the locked stage refuses to run twice.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.impute import SimpleImputer
from sklearn.linear_model import Ridge
from sklearn.model_selection import KFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from fuelblend.uq.conformal import coverage as conf_coverage  # noqa: E402
from fuelblend.uq.conformal import two_sided_margin  # noqa: E402

WIDE = ROOT / "data" / "real_fuels" / "v5" / "ajftd_fuels_wide.csv"
LONG = ROOT / "data" / "real_fuels" / "v5" / "ajftd_fuels_long.csv"
OUT = ROOT / "results" / "v5" / "S4"

SEED = "S4-v5-split-seed-1"
LOCK_FRACTION = 0.20
TIME_TEST_SHARE = 0.25
ALPHA = 0.10
BOOT_B = 5000
BOOT_SEED = 0
MIN_GROUPS = 5
MIN_FEATURES = 3

FEATURES = ["distillation_T0_C", "distillation_T10_C", "distillation_T50_C", "distillation_T90_C",
            "distillation_T100_C", "aromatics_vol_pct", "density_kg_m3_15C"]
TARGET_COL = {"flash": "flash_point_C", "freezing": "freezing_point_C", "visc_m20": "visc_m20",
              "heat": "heat_of_combustion_net_MJ_kg"}
TARGETS = list(TARGET_COL)
MODELS = {"flash": ["const", "api_2b7_t10", "butler_t10", "ridge", "hgb"],
          "freezing": ["const", "ridge", "hgb"], "visc_m20": ["const", "ridge", "hgb"],
          "heat": ["const", "ridge", "hgb"]}
CORRELATIONS = ("api_2b7_t10", "butler_t10")


# --------------------------------------------------------------------------- data
def sha256_file(p: Path) -> str:
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


from fuelblend.eval.group_stats import flash_family as _flash_family  # noqa: E402  (moved, body unchanged)



def build_table(wide_path: Path = WIDE, long_path: Path = LONG):
    """One row per de-duplicated fuel.  Returns (table, dup_map {raw_id: rep_id})."""
    w = pd.read_csv(wide_path)
    lg = pd.read_csv(long_path)

    v = lg[(lg["property"] == "viscosity_kin_cSt") & np.isclose(lg["T_K"], 253.15, atol=0.01)]
    meas = v[v["T_source"].isin(["reported", "reported_parsed_text"])].groupby("fuel_id")["value"].median()
    ass = v[v["T_source"] == "assumed_spec_-20C"].groupby("fuel_id")["value"].median()
    ass = ass[~ass.index.isin(meas.index)]          # a fuel with any measured -20 C row is "measured"
    w["visc_m20"] = w["fuel_id"].map(meas)
    w["visc_m20_assumed"] = w["fuel_id"].map(ass)

    fm = (lg[lg["property"] == "flash_point_C"].groupby("fuel_id")["method"]
          .agg(lambda s: s.mode().iloc[0] if s.notna().any() else "").map(_flash_family))
    w["flash_family"] = w["fuel_id"].map(fm)

    num_cols = [c for c in w.columns if c not in ("fuel_id", "fuel_label", "source_group", "test_year",
                                                  "provenance_class", "density_is_derived", "flash_family")
                and pd.api.types.is_numeric_dtype(w[c])]
    w = w.sort_values("fuel_id").reset_index(drop=True)

    # --- de-duplication (documented in README) ---
    rep = {int(i): int(i) for i in w["fuel_id"]}
    key = w[num_cols].round(6).apply(lambda r: tuple(None if pd.isna(x) else float(x) for x in r), axis=1)
    nn = w[num_cols].notna().sum(axis=1)
    first = {}
    for i, k, n, sg in zip(w["fuel_id"], key, nn, w["source_group"]):
        if n < 6:
            continue
        kk = (sg, k)
        if kk in first:
            rep[int(i)] = first[kk]          # sorted by fuel_id -> lowest id is kept
        else:
            first[kk] = int(i)
    labels = {(sg, lab): int(i) for i, sg, lab in zip(w["fuel_id"], w["source_group"], w["fuel_label"])}
    for i, sg, lab in zip(w["fuel_id"], w["source_group"], w["fuel_label"]):
        m = re.match(r"^(.*)_(\d)$", str(lab))
        if m and (sg, m.group(1)) in labels:
            rep[int(i)] = rep[labels[(sg, m.group(1))]]
    w["rep_id"] = w["fuel_id"].map(lambda i: rep[int(i)])
    dup_map = {int(i): rep[int(i)] for i in w["fuel_id"] if rep[int(i)] != int(i)}

    agg = {c: "median" for c in sorted(set(num_cols) | {"visc_m20", "visc_m20_assumed"})}
    med = w.groupby("rep_id").agg(agg)
    base = w[w["fuel_id"] == w["rep_id"]].set_index("rep_id")[["fuel_label", "source_group", "test_year",
                                                                "flash_family"]]
    n_members = w.groupby("rep_id").size().rename("n_members")
    t = base.join(med, rsuffix="_med").join(n_members)
    t = t.reset_index().rename(columns={"rep_id": "fuel_id"})
    yr = t["test_year"].map(lambda y: "NA" if pd.isna(y) else str(int(y)))
    t["cell"] = t["source_group"] + "|" + yr
    tcols = [TARGET_COL[k] for k in TARGETS] + ["visc_m20_assumed"]
    t = t[t[tcols].notna().any(axis=1)].reset_index(drop=True)
    t["n_feat"] = t[FEATURES].notna().sum(axis=1)
    return t, dup_map


def _vec_key(row, cols):
    return tuple(None if pd.isna(row[c]) else round(float(row[c]), 6) for c in cols)


def make_splits(t: pd.DataFrame, dup_map: dict) -> dict:
    cells = sorted(t["cell"].unique())
    order = sorted(cells, key=lambda c: hashlib.sha256(f"{SEED}|{c}".encode()).hexdigest())
    n_lock = int(math.ceil(LOCK_FRACTION * len(cells)))
    locked_cells = set(order[:n_lock])
    role = {int(i): ("locked" if c in locked_cells else "dev") for i, c in zip(t["fuel_id"], t["cell"])}

    # identical numeric vectors spanning dev/locked are leakage: keep them on the dev side
    cols = FEATURES + [TARGET_COL[k] for k in TARGETS] + ["visc_m20_assumed"]
    clusters: dict = {}
    for _, r in t.iterrows():
        if r[cols].notna().sum() < 6:
            continue
        clusters.setdefault(_vec_key(r, cols), []).append(int(r["fuel_id"]))
    excluded = []
    for ids in clusters.values():
        roles = {role[i] for i in ids}
        if roles == {"dev", "locked"}:
            for i in ids:
                if role[i] == "locked":
                    role[i] = "excluded_leak"
                    excluded.append(i)

    # time split (dev only): latest years until >=25 % of dev fuels with a known year
    dev = t[[role[int(i)] == "dev" for i in t["fuel_id"]]]
    yrs = dev.dropna(subset=["test_year"]).groupby("test_year").size().sort_index(ascending=False)
    tot, acc, test_years = int(yrs.sum()), 0, []
    for y, n in yrs.items():
        test_years.append(int(y))
        acc += int(n)
        if acc >= TIME_TEST_SHARE * tot:
            break
    cutoff = min(test_years)

    cell_info = {}
    for c in cells:
        s = t[t["cell"] == c]
        cell_info[c] = {"role": "locked" if c in locked_cells else "dev", "n_fuels": int(len(s)),
                        "source_group": str(s["source_group"].iloc[0])}
    n_locked_fuels = sum(1 for v in role.values() if v == "locked")
    return {
        "task": "S4", "seed": SEED, "hash": "sha256(seed|cell) ascending",
        "rule": ("group unit = source_group x test_year ('cell'); cells sorted by sha256(seed|cell); the first "
                 f"ceil({LOCK_FRACTION} * n_cells) cells are LOCKED; no re-rolling"),
        "lock_fraction": LOCK_FRACTION, "features": FEATURES,
        "input_sha256": {"wide": sha256_file(WIDE), "long": sha256_file(LONG)},
        "n_raw_fuels_in_wide": int(len(dup_map) + len(t)), "n_fuels_dedup": int(len(t)),
        "n_cells": len(cells), "n_locked_cells": n_lock,
        "locked_fuel_fraction": round(n_locked_fuels / len(t), 4),
        "n_dev_fuels": sum(1 for v in role.values() if v == "dev"),
        "n_locked_fuels": n_locked_fuels, "n_excluded_leak": len(excluded), "excluded_leak_ids": sorted(excluded),
        "cells": cell_info,
        "fuel_role": {str(k): v for k, v in sorted(role.items())},
        "dup_map": {str(k): v for k, v in sorted(dup_map.items())},
        "dedup_rule": ("(1) rows of the same source_group with identical numeric vectors (>=6 non-null values) "
                       "collapse to the lowest fuel_id; (2) 'X_<digit>' labels whose base label X exists in the "
                       "same source_group collapse onto X; remaining members are merged by per-column median."),
        "time_split": {"rule": f"dev fuels with known test_year; test = latest years until >={TIME_TEST_SHARE:.0%} "
                               "of them; train = earlier years; fuels without test_year excluded",
                       "test_years": sorted(test_years), "cutoff_year": cutoff},
    }


def write_splits_guarded(obj: dict, path: Path) -> str:
    text = json.dumps(obj, indent=1, sort_keys=True)
    path = Path(path)
    if path.exists():
        if path.read_text() == text:
            return "unchanged"
        raise FileExistsError(f"{path} exists with different content; refusing to overwrite")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)
    return "written"


def cmd_splits(out: Path = OUT):
    t, dup = build_table()
    sp = make_splits(t, dup)
    status = write_splits_guarded(sp, out / "splits_ajftd_v5.json")
    big = sorted(((v["n_fuels"], c, v["role"]) for c, v in sp["cells"].items()), reverse=True)[:6]
    print(status, "| fuels", sp["n_fuels_dedup"], "cells", sp["n_cells"], "locked cells", sp["n_locked_cells"],
          "locked fuel frac", sp["locked_fuel_fraction"], "excluded_leak", sp["n_excluded_leak"])
    print("largest cells:", big)
    print("time split:", sp["time_split"]["test_years"])


# --------------------------------------------------------------------------- models
class Model:
    def __init__(self, name: str, target: str):
        self.name, self.target = name, target

    def fit(self, X: pd.DataFrame, y):
        y = np.asarray(y, float)
        self.mean_ = float(y.mean())
        yy = np.log(y) if self.target == "visc_m20" else y
        if self.name == "ridge":
            self.est_ = make_pipeline(SimpleImputer(strategy="median", add_indicator=True, keep_empty_features=True),
                                      StandardScaler(), Ridge(alpha=1.0)).fit(X.values, yy)
        elif self.name == "hgb":
            self.est_ = HistGradientBoostingRegressor(max_iter=200, learning_rate=0.05, max_depth=3,
                                                      min_samples_leaf=20, l2_regularization=1.0,
                                                      random_state=0).fit(X.values, yy)
        return self

    def predict(self, X: pd.DataFrame):
        n = len(X)
        if self.name == "const":
            return np.full(n, self.mean_)
        if self.name in CORRELATIONS:
            t10 = X["distillation_T10_C"].to_numpy(float)
            if self.name == "api_2b7_t10":     # API TDB 2B7.1 form, T in K
                T = t10 + 273.15
                inv = -0.024209 + 2.84947 / T + 0.0034254 * np.log(T)
                p = 1.0 / inv - 273.15
            else:                               # Butler (1956): FP[F] = 0.69*T10[F] - 118.2
                p = ((0.69 * (t10 * 1.8 + 32.0) - 118.2) - 32.0) / 1.8
            return np.where(np.isfinite(p), p, self.mean_)
        raw = self.est_.predict(X.values)
        return np.exp(raw) if self.target == "visc_m20" else raw


def eligible(t: pd.DataFrame, target: str, roles: dict, role: str) -> pd.DataFrame:
    col = TARGET_COL[target]
    m = t[col].notna() & (t["n_feat"] >= MIN_FEATURES) & np.array([roles[int(i)] == role for i in t["fuel_id"]])
    d = t[m].copy()
    d["y"] = d[col]
    return d


def fitpred(name, target, train, test):
    return Model(name, target).fit(train[FEATURES], train["y"]).predict(test[FEATURES])


def collect(df, target, scheme, preds_by_model, extra=None):
    rows = []
    for mname, p in preds_by_model.items():
        r = df[["fuel_id", "cell", "source_group", "y", "flash_family"]].copy()
        r["pred"], r["model"], r["target"], r["scheme"] = p, mname, target, scheme
        rows.append(r)
    return pd.concat(rows)


def logo_preds(df, target, group_col="cell", scheme="logo"):
    out = {m: np.full(len(df), np.nan) for m in MODELS[target]}
    g = df[group_col].to_numpy()
    for grp in np.unique(g):
        te = g == grp
        for m in MODELS[target]:
            out[m][te] = fitpred(m, target, df[~te], df[te])
    return collect(df, target, scheme, out)


def kfold_preds(df, target):
    out = {m: np.full(len(df), np.nan) for m in MODELS[target]}
    for tr, te in KFold(5, shuffle=True, random_state=0).split(df):
        for m in MODELS[target]:
            out[m][te] = fitpred(m, target, df.iloc[tr], df.iloc[te])
    return collect(df, target, "within_group_5fold", out)


def time_preds(df, target, cutoff):
    d = df.dropna(subset=["test_year"])
    tr, te = d[d["test_year"] < cutoff], d[d["test_year"] >= cutoff]
    if len(tr) < 30 or len(te) < 10:
        return None
    out = {m: fitpred(m, target, tr, te) for m in MODELS[target]}
    return collect(te, target, "time_split", out)


# --------------------------------------------------------------------------- statistics
def boot_diff(err_m, err_c, groups, seed=BOOT_SEED):
    """Paired bootstrap over groups of the pooled MAE difference (model - const)."""
    d = np.asarray(err_m) - np.asarray(err_c)
    g = pd.Series(groups).astype(str).to_numpy()
    ug, inv = np.unique(g, return_inverse=True)
    s = np.bincount(inv, weights=d)
    n = np.bincount(inv).astype(float)
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, len(ug), size=(BOOT_B, len(ug)))
    stat = s[idx].sum(1) / n[idx].sum(1)
    return float(d.mean()), float(np.percentile(stat, 2.5)), float(np.percentile(stat, 97.5)), len(ug)


def verdict(model, lo, hi, n_groups):
    if model == "const":
        return "baseline"
    if n_groups < MIN_GROUPS:
        return f"not decidable (<{MIN_GROUPS} groups)"
    if hi < 0:
        return "beats const"
    if lo > 0:
        return "worse than const"
    return "no skill shown"


def summarize(preds, margins=None, subset="all"):
    """One row per (scheme, target, model).  margins: {target: {model: margin}} -> conformal coverage."""
    rows = []
    for (scheme, target), P in preds.groupby(["scheme", "target"]):
        P = P.copy()
        P["ae"] = (P["pred"] - P["y"]).abs()
        C = P[P["model"] == "const"].set_index("fuel_id")
        for model, M in P.groupby("model"):
            M = M.set_index("fuel_id")
            Cc = C.loc[M.index]
            diff, lo, hi, ng = boot_diff(M["ae"].to_numpy(), Cc["ae"].to_numpy(), M["cell"].to_numpy())
            mg = (margins or {}).get(target, {}).get(model)
            cov_n = cov_tot = np.nan
            if mg is not None and np.isfinite(mg):
                cov_n, cov_tot = conf_coverage((M["pred"] - M["y"]).to_numpy(), M.index.to_numpy(), mg,
                                               two_sided=True)
            rows.append({"scheme": scheme, "subset": subset, "target": target, "model": model,
                         "n_fuels": len(M), "n_groups": ng, "mae": M["ae"].mean(), "mae_const": Cc["ae"].mean(),
                         "mae_diff_vs_const": diff, "ci_lo": lo, "ci_hi": hi, "verdict": verdict(model, lo, hi, ng),
                         "bias_pred_minus_meas": (M["pred"] - M["y"]).mean(),
                         "n_input_fallback": (int(M["t10_missing"].sum()) if model in CORRELATIONS and
                                              "t10_missing" in M else 0),
                         "margin90": mg if mg is not None else np.nan,
                         "coverage90": (cov_n / cov_tot) if cov_tot and cov_tot == cov_tot else np.nan,
                         "cov_n": cov_n, "cov_tot": cov_tot})
    return pd.DataFrame(rows)


def compute_margins(logo: pd.DataFrame) -> dict:
    """Split-conformal two-sided margins from dev LOGO out-of-fold residuals; conformal group = fuel."""
    out: dict = {}
    for (target, model), M in logo.groupby(["target", "model"]):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            mg = two_sided_margin((M["pred"] - M["y"]).to_numpy(), M["fuel_id"].to_numpy(), ALPHA)
        out.setdefault(target, {})[model] = mg
    return out


# --------------------------------------------------------------------------- stages
def load_all(out: Path):
    sp = json.loads((out / "splits_ajftd_v5.json").read_text())
    t, dup = build_table()
    roles = {int(k): v for k, v in sp["fuel_role"].items()}
    missing = set(int(i) for i in t["fuel_id"]) - set(roles)
    if missing or sp["input_sha256"]["wide"] != sha256_file(WIDE) or sp["input_sha256"]["long"] != sha256_file(LONG):
        raise RuntimeError("splits file does not match current data; refusing to continue")
    return sp, t, roles


def add_flags(preds, t):
    f = t.set_index("fuel_id")["distillation_T10_C"].isna()
    preds["t10_missing"] = preds["fuel_id"].map(f)
    preds["year"] = preds["fuel_id"].map(t.set_index("fuel_id")["test_year"])
    return preds


def assumed_rows(t, roles, role):
    m = t["visc_m20_assumed"].notna() & (t["n_feat"] >= MIN_FEATURES) & \
        np.array([roles[int(i)] == role for i in t["fuel_id"]])
    d = t[m].copy()
    d["y"] = d["visc_m20_assumed"]
    return d


def assumed_table(train, test, scheme):
    rows = []
    if len(test) == 0:
        return pd.DataFrame(columns=["scheme", "target", "model", "n_fuels", "n_groups", "mae",
                                     "bias_pred_minus_meas"])
    for m in MODELS["visc_m20"]:
        p = fitpred(m, "visc_m20", train, test)
        e = np.abs(p - test["y"].to_numpy())
        rows.append({"scheme": scheme, "target": "visc_m20_ASSUMED_T", "model": m, "n_fuels": len(test),
                     "n_groups": test["cell"].nunique(), "mae": e.mean(),
                     "bias_pred_minus_meas": float((p - test["y"].to_numpy()).mean())})
    return pd.DataFrame(rows)


def cmd_dev(out: Path = OUT):
    readme = out / "README.md"
    if not readme.exists() or "DECLARED BEFORE SCORING" not in readme.read_text():
        raise RuntimeError("README.md with declared readings must exist before any scoring")
    sp, t, roles = load_all(out)
    cutoff = sp["time_split"]["cutoff_year"]
    logo_all, within_all, time_all, coarse_all = [], [], [], []
    for tg in TARGETS:
        d = eligible(t, tg, roles, "dev")
        print(tg, "dev fuels", len(d), "cells", d["cell"].nunique(), flush=True)
        logo_all.append(logo_preds(d, tg))
        within_all.append(kfold_preds(d, tg))
        tp = time_preds(d, tg, cutoff)
        if tp is not None:
            time_all.append(tp)
        coarse_all.append(logo_preds(d, tg, group_col="source_group", scheme="logo_coarse_source_group"))
    logo = add_flags(pd.concat(logo_all), t)
    margins = compute_margins(logo)
    allp = add_flags(pd.concat(logo_all + within_all + time_all + coarse_all), t)
    scores = [summarize(allp, margins, "all")]
    nj = allp[allp["source_group"] != "JFK"]
    scores.append(summarize(nj, margins, "excl_JFK"))
    sc = pd.concat(scores, ignore_index=True)
    sc.to_csv(out / "dev_scores.csv", index=False)
    allp.drop(columns=["year"]).to_csv(out / "dev_predictions.csv", index=False)
    (out / "margins_dev.json").write_text(json.dumps(
        {t_: {m: (v if np.isfinite(v) else "inf") for m, v in d_.items()} for t_, d_ in margins.items()}, indent=1))

    fam = logo[(logo["target"] == "flash")].copy()
    fam["ae"] = (fam["pred"] - fam["y"]).abs()
    fam.groupby(["flash_family", "model"]).agg(n=("ae", "size"), mae=("ae", "mean"),
                                               bias=("pred", lambda s: 0)).drop(columns="bias") \
        .to_csv(out / "dev_flash_by_method_family.csv")

    # assumed-temperature viscosity rows, scored separately (models trained on measured dev rows only)
    tr = eligible(t, "visc_m20", roles, "dev")
    te = assumed_rows(t, roles, "dev")
    assumed_table(tr, te, "dev_fit_all_dev_measured").to_csv(out / "dev_assumed_visc.csv", index=False)

    meta = {"splits_sha256": sha256_file(out / "splits_ajftd_v5.json"), "readme_sha256": sha256_file(readme),
            "script_sha256": sha256_file(Path(__file__)), "n_dev_fuels": sp["n_dev_fuels"]}
    (out / "dev_meta.json").write_text(json.dumps(meta, indent=1))
    show = sc[(sc["subset"] == "all") & (sc["scheme"] == "logo")][
        ["target", "model", "n_fuels", "n_groups", "mae", "mae_const", "ci_lo", "ci_hi", "verdict", "margin90"]]
    print(show.to_string(index=False))


def cmd_locked(out: Path = OUT):
    ls = out / "locked_scores.csv"
    if ls.exists() or (out / "locked_run.json").exists():
        raise FileExistsError("locked set already scored; it may be scored exactly once")
    meta_p = out / "dev_meta.json"
    if not meta_p.exists():
        raise RuntimeError("dev stage must be run first")
    meta = json.loads(meta_p.read_text())
    if meta["splits_sha256"] != sha256_file(out / "splits_ajftd_v5.json") or \
            meta["readme_sha256"] != sha256_file(out / "README.md"):
        raise RuntimeError("splits/README changed since the dev stage")
    sp, t, roles = load_all(out)
    margins = {t_: {m: (float("inf") if v == "inf" else float(v)) for m, v in d_.items()}
               for t_, d_ in json.loads((out / "margins_dev.json").read_text()).items()}
    (out / "locked_run.json").write_text(json.dumps({"note": "locked set scored once", "script_sha256":
                                                     sha256_file(Path(__file__))}))
    parts = []
    for tg in TARGETS:
        tr, te = eligible(t, tg, roles, "dev"), eligible(t, tg, roles, "locked")
        out_p = {m: fitpred(m, tg, tr, te) for m in MODELS[tg]}
        parts.append(collect(te, tg, "locked", out_p))
        print(tg, "locked fuels", len(te), "cells", te["cell"].nunique(), flush=True)
    P = add_flags(pd.concat(parts), t)
    sc = pd.concat([summarize(P, margins, "all"), summarize(P[P["source_group"] != "JFK"], margins, "excl_JFK")],
                   ignore_index=True)
    sc.to_csv(ls, index=False)
    P.drop(columns=["year"]).to_csv(out / "locked_predictions.csv", index=False)
    assumed_table(eligible(t, "visc_m20", roles, "dev"), assumed_rows(t, roles, "locked"),
                  "locked_fit_all_dev_measured").to_csv(out / "locked_assumed_visc.csv", index=False)
    print(sc[sc["subset"] == "all"][["target", "model", "n_fuels", "n_groups", "mae", "mae_const", "ci_lo", "ci_hi",
                                      "verdict", "coverage90"]].to_string(index=False))


def cmd_report(out: Path = OUT):
    dv, lk = pd.read_csv(out / "dev_scores.csv"), pd.read_csv(out / "locked_scores.csv")
    sp = json.loads((out / "splits_ajftd_v5.json").read_text())
    cols = ["target", "model", "n_fuels", "n_groups", "mae", "mae_const", "ci_lo", "ci_hi", "verdict"]
    L = ["# S4 results (generated by `eval_bulk_descriptor.py report`)", "",
         f"dev fuels {sp['n_dev_fuels']}, locked fuels {sp['n_locked_fuels']} "
         f"({sp['locked_fuel_fraction']:.1%} of {sp['n_fuels_dedup']}), cells {sp['n_cells']} "
         f"({sp['n_locked_cells']} locked). Verdicts follow the README reading.", ""]
    def sect(title, df):
        L.append(f"## {title}\n")
        L.append("```" + chr(10) + df.round(3).to_string(index=False) + chr(10) + "```")
        L.append("")
    for scheme in ["logo", "within_group_5fold", "time_split", "logo_coarse_source_group"]:
        sect(f"dev / {scheme}", dv[(dv.scheme == scheme) & (dv["subset"] == "all")][cols + ["margin90", "coverage90"]])
    sect("dev / logo, excluding JFK", dv[(dv.scheme == "logo") & (dv["subset"] == "excl_JFK")][cols])
    sect("LOCKED (scored once)", lk[lk["subset"] == "all"][cols + ["margin90", "coverage90"]])
    sect("LOCKED, excluding JFK", lk[lk["subset"] == "excl_JFK"][cols])
    for f in ["dev_assumed_visc.csv", "locked_assumed_visc.csv", "dev_flash_by_method_family.csv"]:
        sect(f, pd.read_csv(out / f))
    (out / "RESULTS.md").write_text("\n".join(L))
    print("wrote RESULTS.md")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("stage", choices=["splits", "dev", "locked", "report"])
    a = ap.parse_args()
    {"splits": cmd_splits, "dev": cmd_dev, "locked": cmd_locked, "report": cmd_report}[a.stage]()


if __name__ == "__main__":
    main()
