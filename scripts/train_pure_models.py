#!/usr/bin/env python
"""WP4: pure-compound property models - train (HGB, multi-seed), evaluate, write artifacts.

    python scripts/train_pure_models.py

Data: data/processed/model_table.csv (HC-only) + data/unified/compounds.csv.
Splits: splits/s1s2s3.json (S1/S2 train pool, S3 all-novel-molecule holdout).

For every trained pure property:
  - train pool = pure rows of S1 + S2 (property-scoped), minus a held-back 20 % molecule set
    (GroupShuffleSplit by inchikey, seed 20261030) -> held-back = eval (interpolation)
  - in-training validation = GroupShuffleSplit by system_id (seed 20261030, test_size=0.15);
    HGB early stopping (validation_fraction=0.15) on that split
  - final eval = S3 (all-novel molecules, extrapolation) + held-back molecule set

Outputs:
  data/runs/pure_preds.csv            per-row predictions (pred_mean / pred_std over seeds)
  data/runs/pure_feature_cache.csv   keyed by (inchikey, T_K)
  data/runs/PURE_EVAL_LOG.md          verbatim run log
  runs/pure_models/<prop>__seed_<s>.joblib   one fitted model per property+seed (full config)
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd

BASE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(BASE / "src"))

from fuelblend.features import pure_features as pf  # noqa: E402
from fuelblend.models import pure_gbm as pgm  # noqa: E402
from fuelblend.eval import metrics as evm  # noqa: E402

from sklearn.inspection import permutation_importance  # noqa: E402
from sklearn.model_selection import GroupShuffleSplit  # noqa: E402
import joblib  # noqa: E402

MODEL_TABLE = BASE / "data/processed/model_table.csv"
COMPOUNDS = BASE / "data/unified/compounds.csv"
SPLIT_FILE = BASE / "splits/s1s2s3.json"
PRED_OUT = BASE / "data/runs/pure_preds.csv"
LOG_OUT = BASE / "data/runs/PURE_EVAL_LOG.md"
CACHE_OUT = BASE / "data/runs/pure_feature_cache.csv"
MODEL_OUT = BASE / "runs/pure_models"

LOG: list[str] = []


def log(msg) -> None:
    s = str(msg)
    LOG.append(s)
    print(s, flush=True)


def flush_log() -> None:
    LOG_OUT.parent.mkdir(parents=True, exist_ok=True)
    LOG_OUT.write_text("\n".join(LOG) + "\n", encoding="utf-8")


SPLITS = json.loads(SPLIT_FILE.read_text(encoding="utf-8"))
S1 = set(int(x) for x in SPLITS["S1"])
S2 = set(int(x) for x in SPLITS["S2"])
S3 = set(int(x) for x in SPLITS["S3"])

HELDOUT_SEED = 20261030
HELDOUT_FRAC = 0.2
MIN_PURE_TRAIN_ROWS = 200
MODEL_SEEDS = [11, 23, 37, 41, 59]
TRAINED_PROPS = ["viscosity_dyn_mPas", "density_kg_m3", "flash_point_C", "vapor_pressure_kPa"]
SKIPPED_PROPS = ["freezing_point_C", "pour_point_C"]
CONFIGS = [
    ("full", dict(sw=True, with_classes=True)),
    ("no_sw", dict(sw=False, with_classes=True)),
    ("no_class", dict(sw=True, with_classes=False)),
]

_mt = pd.read_csv(MODEL_TABLE, low_memory=False)
_comp = pd.read_csv(COMPOUNDS, low_memory=False)
_COMP_LOOKUP = {str(r.inchikey): dict(smiles=r.smiles, cls=r.cls, mw=r.mw, n_carbon=r.n_carbon) for r in _comp.itertuples()}


def _jload(x):
    if isinstance(x, list):
        return x
    if x is None or (isinstance(x, float) and math.isnan(x)):
        return []
    try:
        return json.loads(str(x))
    except Exception:
        return []


def _eff_T(r) -> float:
    t = float(r["T_K"]) if not pd.isna(r["T_K"]) else float("nan")
    if math.isnan(t):
        v = float(r["value"]) if not pd.isna(r["value"]) else float(pf.DEFAULT_T_K)
        t = v + 273.15
    return t


def collect(prop: str, positions) -> list[dict]:
    recs = []
    for i in sorted(positions):
        r = _mt.iloc[i]
        if r["kind"] != "pure" or r["property"] != prop:
            continue
        ik = _jload(r["comp_inchikeys"])
        if len(ik) != 1:
            continue
        c = _COMP_LOOKUP.get(str(ik[0]), {})
        recs.append(dict(
            inchikey=str(ik[0]), smiles=c.get("smiles"), T_K=_eff_T(r), P_kPa=float(r["P_kPa"]) if not pd.isna(r["P_kPa"]) else float("nan"),
            value=float(r["value"]), quality_weight=float(r["quality_weight"]), cls=c.get("cls"),
            system_id=str(r["system_id"]),
        ))
    return recs


def _target_y(prop: str, recs):
    vals = np.array([r["value"] for r in recs], dtype=float)
    if prop == "viscosity_dyn_mPas":
        return np.log(vals), float
    if prop == "vapor_pressure_kPa":
        return np.log10(vals), lambda p: float(10 ** p)
    return vals, lambda p: float(p)


def prepare(prop: str, recs):
    rows = [dict(inchikey=r["inchikey"], smiles=r["smiles"], T_K=r["T_K"], cls=r["cls"]) for r in recs]
    X = pf.pure_feature_frame(rows)
    good = np.isfinite(X).all(axis=1) if len(X) else np.zeros((0, X.shape[1]), dtype=bool)
    recs = [r for r, g in zip(recs, good) if bool(g)]
    X = X[list(good)] if len(X) else X
    y, back = _target_y(prop, recs)
    sw = np.array([r["quality_weight"] for r in recs], dtype=float)
    return recs, X, np.asarray(y, dtype=float), back, sw


def pred_value(prop: str, p) -> float:
    if prop == "viscosity_dyn_mPas":
        return float(np.exp(p))
    if prop == "vapor_pressure_kPa":
        return float(10 ** p)
    return float(p)


def evaluate(prop: str, y_true: np.ndarray, y_pred: np.ndarray, clusters) -> dict:
    mode = "viscosity" if "viscosity" in prop else ("flash" if prop == "flash_point_C" else "generic")
    y_true, y_pred = np.asarray(y_true, float), np.asarray(y_pred, float)
    rep = evm.report(y_true, y_pred, np.asarray(clusters), mode=mode)
    out = dict(rep)
    out["n"] = int(len(y_true))
    out["bootstrap_ci"] = evm.bootstrap_ci(y_true, y_pred, np.asarray(clusters), metric_fn=evm.rmse, n_boot=2000, seed=11)
    if prop == "viscosity_dyn_mPas":
        yt, yp = np.exp(y_true), np.exp(y_pred)
        out["rmse_log"] = float(np.sqrt(np.mean((y_true - y_pred) ** 2)))
        out["mape_back_pct"] = float(np.mean(np.abs(yt - yp) / np.clip(np.abs(yt), 1e-12, None)) * 100)
    if prop == "vapor_pressure_kPa":
        yt, yp = 10 ** y_true, 10 ** y_pred
        out["rmse_log10"] = float(np.sqrt(np.mean((y_true - y_pred) ** 2)))
        out["mape_back_pct"] = float(np.mean(np.abs(yt - yp) / np.clip(np.abs(yt), 1e-12, None)) * 100)
    return out


def fmt(x, nd=4) -> str:
    return f"{float(x):.{nd}f}"


def main() -> None:
    log("WP4 pure-compound property models - train + evaluate + write artifacts")
    log("command: python scripts/train_pure_models.py")
    log("model seeds: [11, 23, 37, 41, 59]; split seed 20261030; held-out molecule fraction 0.20")
    log("")

    cache = pf.load_cache(CACHE_OUT)
    pred_rows: list[dict] = []
    results: dict = {}

    for prop in TRAINED_PROPS:
        log(f"=== property {prop}")
        pool = collect(prop, S1 | S2)
        mols = sorted({r["inchikey"] for r in pool})
        rng = np.random.default_rng(HELDOUT_SEED)
        k = max(1, int(round(HELDOUT_FRAC * len(mols))))
        hold_mols = set(mols[i] for i in rng.choice(len(mols), size=k, replace=False).tolist())
        heldback_recs = [r for r in pool if r["inchikey"] in hold_mols]
        train_recs = [r for r in pool if r["inchikey"] not in hold_mols]
        s3_recs = collect(prop, S3)

        train2, X_tr, y_tr, back_tr, sw_tr = prepare(prop, train_recs)
        groups = np.array([r["system_id"] for r in train2])
        trn, val = list(GroupShuffleSplit(n_splits=5, test_size=0.15, random_state=HELDOUT_SEED).split(X_tr, y_tr, groups=groups))[0]
        train_gss_recs, X_tr_gss, y_tr_gss, back_tr_gss, sw_tr_gss = prepare(prop, [train2[i] for i in trn])
        val_recs, X_val, y_val, back_val, sw_val = prepare(prop, [train2[i] for i in val])
        hold_recs, X_hold, y_hold, back_hold, sw_hold = prepare(prop, heldback_recs)
        s3_r, X_s3, y_s3, back_s3, sw_s3 = prepare(prop, s3_recs)

        need: dict = {}
        for r in heldback_recs + s3_recs:
            key = pf.cache_key(r["inchikey"], r["T_K"])
            if key not in need:
                need[key] = (r["inchikey"], r["T_K"], r["smiles"], r["cls"])
        cache_rows: dict = {}
        cached = pf.load_cache(CACHE_OUT)
        feat_names = [n for n in pf.FEATURE_NAMES if n != "T_K"]
        if len(cached):
            for r in cached.to_dict("records"):
                key = (str(r["inchikey"]), float(r["T_K"]))
                cache_rows[key] = [key[0], key[1]] + [float(r[n]) for n in feat_names]
        for key, (ik, Tk, smi, cl) in sorted(need.items()):
            if key in cache_rows:
                continue
            row = pf.pure_feature_row(ik, smi, Tk, cl)
            cache_rows[key] = [ik, Tk] + [float(row[pf.FEATURE_NAMES.index(n)]) for n in feat_names]
        cache_df = pd.DataFrame(
            [[k[0], k[1]] + v[2:] for k, v in sorted(cache_rows.items())], columns=pf.CACHE_COLUMNS)
        pf.save_cache(cache_df, CACHE_OUT)
        cache = cache_df

        per_cfg: dict = {}
        imp_model = None
        for cfg_name, cfg in CONFIGS:
            per_seed: dict = {"heldback": {}, "S3": {}}
            for seed in MODEL_SEEDS:
                model = pgm.make_model(random_state=seed)
                Xc_fit = X_tr_gss if cfg["with_classes"] else pf.drop_class_columns(X_tr_gss)
                swc = sw_tr_gss if cfg["sw"] else None
                pgm.fit_pure_gbm(model, Xc_fit, y_tr_gss, sample_weight=swc)
                if cfg_name == "full":
                    MODEL_OUT.mkdir(parents=True, exist_ok=True)
                    joblib.dump(model, MODEL_OUT / f"{prop}__seed_{seed}.joblib")
                    if seed == MODEL_SEEDS[0]:
                        imp_model = model
                p_hold = pgm.predict_pure_gbm(model, X_hold if cfg["with_classes"] else pf.drop_class_columns(X_hold))
                p_s3 = pgm.predict_pure_gbm(model, X_s3 if cfg["with_classes"] else pf.drop_class_columns(X_s3))
                per_seed["heldback"][seed] = p_hold
                per_seed["S3"][seed] = p_s3

            split_stats: dict = {}
            for split_name, X_e, y_e, recs_e, per_seed_map in (
                ("heldback", X_hold, y_hold, hold_recs, per_seed["heldback"]),
                ("S3", X_s3, y_s3, s3_r, per_seed["S3"]),
            ):
                P = np.array([per_seed_map[s] for s in MODEL_SEEDS], dtype=float)
                mean_p = P.mean(axis=0)
                std_p = P.std(axis=0)
                stats = {}
                for s in MODEL_SEEDS:
                    stats[s] = evaluate(prop, y_e, per_seed_map[s], np.array([r["system_id"] for r in recs_e]))
                stats["mean"] = evaluate(prop, y_e, mean_p, np.array([r["system_id"] for r in recs_e]))
                stats["pred_mean"] = mean_p
                stats["pred_std"] = std_p
                split_stats[split_name] = stats

            log(f"RESULT {prop} cfg={cfg_name} splits=" + str({k: float(v['mean']['rmse']) for k, v in split_stats.items()}))

            if cfg_name == "full":
                for split_name, X_e, y_e, recs_e, stats in (
                    ("heldback", X_hold, y_hold, hold_recs, split_stats["heldback"]),
                    ("S3", X_s3, y_s3, s3_r, split_stats["S3"]),
                ):
                    pred_mean, pred_std = stats["pred_mean"], stats["pred_std"]
                    for i, r in enumerate(recs_e):
                        pred_rows.append(dict(
                            property=prop, inchikey=r["inchikey"], T_K=float(r["T_K"]),
                            P_kPa=float(r["P_kPa"]) if not pd.isna(r["P_kPa"]) else float("nan"),
                            value=float(r["value"]), pred_mean=float(pred_mean[i]), pred_std=float(pred_std[i]),
                            split=split_name, quality_weight=float(r["quality_weight"]),
                        ))
                for split_name in ("heldback", "S3"):
                    st = split_stats[split_name]
                    m = st["mean"]
                    log(f"RESULT {prop} cfg={cfg_name} split={split_name} n={m['n']} mae={fmt(m['mae'], 4)} "
                        f"rmse={fmt(m['rmse'], 4)} median_ape={fmt(m.get('median_ape', float('nan')), 4)} "
                        f"frac_within_10={fmt(m.get('frac_within_10', float('nan')), 4)} bias={fmt(m.get('bias', float('nan')), 4)} "
                        f"overpred_rate_5={fmt(m.get('overpred_rate_5', float('nan')), 4)} mean_pos_err={fmt(m.get('mean_positive_error', float('nan')), 4)} "
                        f"rmse_log={fmt(m.get('rmse_log', float('nan')), 6)} mape_back_pct={fmt(m.get('mape_back_pct', float('nan')), 4)} "
                        f"ci_lo={fmt(m['bootstrap_ci']['lo'], 4)} ci_hi={fmt(m['bootstrap_ci']['hi'], 4)} "
                        f"pred_mean_spread={float(np.max([st[s]['rmse'] for s in MODEL_SEEDS]) - np.min([st[s]['rmse'] for s in MODEL_SEEDS])):.4f}")

        # permutation importance on in-training validation (full config, first seed model)
        imp = permutation_importance(imp_model, X_val, y_val, sample_weight=sw_val, n_repeats=3, random_state=11)
        imps = np.asarray(imp.importances, dtype=float).ravel().tolist()
        top = sorted(zip(pf.FEATURE_NAMES, imps), key=lambda t: -abs(t[1]))[:5]
        results[prop] = dict(
            pool_n=len(pool), heldback_n=len(hold_recs), s3_n=len(s3_r),
            train_gss_n=len(train_gss_recs), val_n=len(val_recs),
            top_features=top,
        )
        log(f"RESULT {prop} pool_n={len(pool)} heldback_n={len(hold_recs)} s3_n={len(s3_r)} "
            f"train_gss_n={len(train_gss_recs)} val_n={len(val_recs)} top_features={top}")

    PRED_OUT.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(pred_rows, columns=["property", "inchikey", "T_K", "P_kPa", "value", "pred_mean", "pred_std", "split", "quality_weight"]) \
        .to_csv(PRED_OUT, index=False)
    log(f"wrote {PRED_OUT} rows={len(pred_rows)}")

    for prop, res in results.items():
        log(f"{prop}: pool={res['pool_n']} heldback={res['heldback_n']} s3={res['s3_n']} "
            f"train_gss={res['train_gss_n']} val={res['val_n']} top_features={res['top_features']}")

    log("DONE")


if __name__ == "__main__":
    main()
    flush_log()
