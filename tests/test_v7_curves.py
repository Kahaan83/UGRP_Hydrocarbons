import numpy as np
import pandas as pd
import pytest

pytest.importorskip("torch")
from fuelblend.ml.curves import CurveLogger, lgb_curve_callback  # noqa: E402
from fuelblend.v7 import golden as G  # noqa: E402
from fuelblend.v7.guard import ROOT  # noqa: E402


def _g(n):
    return np.load(ROOT / f"tests/data/v7_golden_{n}.npy")


class Boom(Exception):
    pass


def test_fit_ensemble_logger_none_bitwise():
    assert np.array_equal(G.dmpnn_run(), _g("dmpnn"))


def test_molpool_defaults_bitwise():
    assert np.array_equal(G.molpool_run(), _g("molpool"))


def test_set_models_defaults_bitwise():
    assert np.array_equal(G.set_models_run(), _g("set_models"))


def test_on_epoch_prune_propagates():
    def cb(ep, vl):
        if ep == 3:
            raise Boom()
    with pytest.raises(Boom):
        G.dmpnn_run(on_epoch=cb)
    with pytest.raises(Boom):
        G.molpool_run(on_epoch=cb)
    with pytest.raises(Boom):
        G.set_models_run(on_epoch=cb)


def test_logger_rows_per_epoch(tmp_path):
    lg = CurveLogger(tmp_path / "c.csv", run={"model": "dmpnn", "fold": 0})
    G.dmpnn_run(logger=lg)
    lg.close()
    d = pd.read_csv(tmp_path / "c.csv")
    assert {"train_loss", "val_loss", "val_0", "val_1"} <= set(d.columns)
    for s, g in d.groupby("seed"):
        assert g.epoch.is_monotonic_increasing and g.epoch.iloc[0] == 1
    assert set(d.seed) == {0, 1} and (d.groupby("seed").size() <= 6).all()


def test_molpool_early_stop_restores_and_logs(tmp_path):
    lg = CurveLogger(tmp_path / "m.csv", run={"model": "molpool"})
    G.molpool_run(logger=lg, val_idx=np.arange(5), patience=2)
    lg.close()
    d = pd.read_csv(tmp_path / "m.csv")
    assert d.val_loss.notna().all() and len(d) >= 2


def test_curve_csv_partial_read(tmp_path):
    lg = CurveLogger(tmp_path / "p.csv", run={"model": "x"})
    for e in range(1, 4):
        lg.log(e, train_loss=1.0 / e, val_loss=2.0 / e)
    assert len(pd.read_csv(tmp_path / "p.csv")) == 3
    lg.close()


def test_lgb_callback_records(tmp_path):
    import lightgbm as lgb
    g = np.random.default_rng(0)
    X = g.normal(size=(50, 4)); y = X[:, 0] + g.normal(size=50) * .1
    lg = CurveLogger(tmp_path / "l.csv", run={"model": "lgb"})
    dt, dv = lgb.Dataset(X[:35], y[:35]), lgb.Dataset(X[35:], y[35:])
    lgb.train({"objective": "l2", "verbose": -1, "min_data_in_leaf": 3}, dt, num_boost_round=20, valid_sets=[dv],
              callbacks=[lgb_curve_callback(lg)])
    lg.close()
    d = pd.read_csv(tmp_path / "l.csv")
    assert len(d) == 20 and d.val_loss.notna().all()
