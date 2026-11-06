"""Metrics with cluster bootstrap by system_id."""
import numpy as np


def mae(y_true, y_pred):
    return float(np.mean(np.abs(y_true - y_pred)))


def rmse(y_true, y_pred):
    return float(np.sqrt(np.mean((y_true - y_pred) ** 2)))


def median_ape(y_true, y_pred):
    return float(np.median(np.abs((y_true - y_pred) / np.clip(np.abs(y_true), 1e-9, None))) * 100)


def frac_within(y_true, y_pred, tol=0.1):
    rel = np.abs((y_true - y_pred) / np.clip(np.abs(y_true), 1e-9, None))
    return float(np.mean(rel <= tol))


def flash_overpred_rate(y_true, y_pred, thr=5.0):
    return float(np.mean((y_pred - y_true) > thr))


def mean_positive_error(y_true, y_pred):
    d = y_pred - y_true
    return float(np.mean(d[d > 0])) if len(d[d > 0]) else 0.0


def bootstrap_ci(y_true, y_pred, clusters, metric_fn, n_boot=2000, seed=11, alpha=0.05):
    y_true = np.asarray(y_true, float)
    y_pred = np.asarray(y_pred, float)
    clusters = np.asarray(clusters)
    center = float(metric_fn(y_true, y_pred))
    rng = np.random.default_rng(seed)
    cl = np.unique(clusters)
    idx_by_cl = {c: np.where(clusters == c)[0] for c in cl}
    boots = []
    for _ in range(n_boot):
        sampled = rng.choice(cl, size=len(cl), replace=True)
        idx = np.concatenate([idx_by_cl[c] for c in sampled])
        try:
            boots.append(metric_fn(y_true[idx], y_pred[idx]))
        except Exception:
            boots.append(np.nan)
    lo = float(np.nanpercentile(boots, 100 * alpha / 2))
    hi = float(np.nanpercentile(boots, 100 * (1 - alpha / 2)))
    return {"point": center, "lo": lo, "hi": hi}


def report(y_true, y_pred, clusters, mode="generic"):
    out = {"n": len(y_true), "mae": mae(y_true, y_pred), "rmse": rmse(y_true, y_pred)}
    if mode == "viscosity":
        out["median_ape"] = median_ape(y_true, y_pred)
        out["frac_within_10"] = frac_within(y_true, y_pred, 0.1)
    elif mode == "flash":
        out["bias"] = float(np.mean(y_pred - y_true))
        out["overpred_rate_5"] = flash_overpred_rate(y_true, y_pred, 5.0)
        out["mean_positive_error"] = mean_positive_error(y_true, y_pred)
    return out
