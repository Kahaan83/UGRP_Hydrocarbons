"""Nested Optuna harness: inner GroupKFold(3) on the training groups of each outer fold, resumable SQLite studies.

Hyperparameters are chosen in inner folds only; each outer fold is scored once with its chosen config (an existing
predictions file for a fold means that fold is skipped)."""
from __future__ import annotations

import inspect
import json
import os
from pathlib import Path
from typing import Callable

import numpy as np
import optuna
import pandas as pd
import yaml
from sklearn.model_selection import GroupKFold

from fuelblend.v7.guard import ROOT

SPACES_YAML = ROOT / "configs/v7/search_spaces.yaml"


def default_storage() -> str:
    runs = os.environ.get("FUELBLEND_RUNS", r"C:\fuelblend-gpu\runs")
    Path(runs).mkdir(parents=True, exist_ok=True)
    return "sqlite:///" + (Path(runs) / "optuna.db").as_posix()


def load_space(name: str) -> dict:
    return yaml.safe_load(SPACES_YAML.read_text(encoding="utf-8"))[name]


def suggest(trial: optuna.Trial, space: dict) -> dict:
    """Sample the 'params' of a space and merge its 'fixed' entries."""
    p = dict(space.get("fixed", {}))
    for k, s in space["params"].items():
        t = s["type"]
        if t == "float":
            p[k] = trial.suggest_float(k, s["low"], s["high"], log=s.get("log", False))
        elif t == "int":
            p[k] = trial.suggest_int(k, s["low"], s["high"], log=s.get("log", False))
        elif t == "cat":
            p[k] = trial.suggest_categorical(k, s["choices"])
        else:
            raise ValueError(f"bad space type {t!r} for {k}")
    return p


def _call_fit_score(fit_score, params, tr, va, logger, on_epoch):
    if on_epoch is not None and "on_epoch" in inspect.signature(fit_score).parameters:
        return fit_score(params, tr, va, logger, on_epoch=on_epoch)
    return fit_score(params, tr, va, logger)


def nested_cv(df, groups, outer_folds, *, fit_score: Callable, space: dict, n_trials: int, timeout_s: int,
              study_prefix: str, storage: str | None = None, final_fit: Callable | None = None,
              out_dir=None, logger=None, inner_k: int = 3, seed: int = 0) -> pd.DataFrame:
    """One row per outer fold: fold, best_params (json), inner_score, outer_pred_path.

    fit_score(params, idx_tr, idx_va, logger[, on_epoch]) -> float, lower is better (indices are positions in df).
    final_fit(params, idx_tr, idx_te) -> array of outer predictions, written to out_dir/pred_fold{k}.npy once.
    Pruning: inner fold 0 reports (epoch, val_loss) through on_epoch to a MedianPruner(n_warmup_steps=20)."""
    storage = storage or default_storage()
    groups = np.asarray(groups)
    outer = np.asarray(outer_folds)
    out = Path(out_dir) if out_dir is not None else None
    if out is not None:
        out.mkdir(parents=True, exist_ok=True)
    rows = []
    for k in sorted(np.unique(outer)):
        te = np.where(outer == k)[0]
        tr = np.where(outer != k)[0]
        pred_path = None if out is None else out / f"pred_fold{k}.npy"
        done = out is not None and pred_path.exists()
        study = optuna.create_study(
            study_name=f"{study_prefix}_fold{k}", storage=storage, load_if_exists=True, direction="minimize",
            sampler=optuna.samplers.TPESampler(seed=seed + int(k)),
            pruner=optuna.pruners.MedianPruner(n_warmup_steps=20))
        inner = list(GroupKFold(inner_k).split(tr, groups=groups[tr]))

        def objective(trial, tr=tr, inner=inner, k=k):
            params = suggest(trial, space)
            scores = []
            for j, (a, b) in enumerate(inner):
                cb = None
                if j == 0:
                    def cb(ep, vl, trial=trial):
                        trial.report(float(vl), ep)
                        if trial.should_prune():
                            raise optuna.TrialPruned()
                lg = None if logger is None else logger.child(fold=k, trial=trial.number, inner=j)
                scores.append(float(_call_fit_score(fit_score, params, tr[a], tr[b], lg, cb)))
            return float(np.mean(scores))

        if not done:
            finished = [t for t in study.trials if t.state.is_finished()]
            remaining = max(0, n_trials - len(finished))
            if remaining:
                study.optimize(objective, n_trials=remaining, timeout=timeout_s)
        bp = {**space.get("fixed", {}), **study.best_params}
        row = dict(fold=int(k), best_params=json.dumps(bp, sort_keys=True), inner_score=study.best_value,
                   outer_pred_path="" if pred_path is None else str(pred_path))
        if final_fit is not None and out is not None and not done:
            np.save(pred_path, np.asarray(final_fit(bp, tr, te)))
        rows.append(row)
    return pd.DataFrame(rows)
