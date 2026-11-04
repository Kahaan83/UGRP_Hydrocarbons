"""Pure-compound GBM model wrapper (WP4).

HistGradientBoostingRegressor with in-training early stopping (15 % validation split).
min_samples_leaf is kept small so small pure-compound training sets still fit.
"""
from __future__ import annotations

import numpy as np
from sklearn.ensemble import HistGradientBoostingRegressor

MODEL_SEEDS = [11, 23, 37, 41, 59]
HGB_MAX_ITER = 1000
HGB_LEARNING_RATE = 0.05
HGB_VALIDATION_FRACTION = 0.15
HGB_N_ITER_NO_CHANGE = 10
HGB_MIN_SAMPLES_LEAF = 5


def make_model(random_state: int = 11) -> HistGradientBoostingRegressor:
    return HistGradientBoostingRegressor(
        loss="squared_error",
        learning_rate=HGB_LEARNING_RATE,
        max_iter=HGB_MAX_ITER,
        min_samples_leaf=HGB_MIN_SAMPLES_LEAF,
        early_stopping=True,
        validation_fraction=HGB_VALIDATION_FRACTION,
        n_iter_no_change=HGB_N_ITER_NO_CHANGE,
        random_state=random_state,
    )


def fit_pure_gbm(model, X, y, sample_weight=None) -> HistGradientBoostingRegressor:
    model.fit(X, y, sample_weight=sample_weight)
    return model


def predict_pure_gbm(model, X) -> np.ndarray:
    return np.asarray(model.predict(X), dtype=float)
