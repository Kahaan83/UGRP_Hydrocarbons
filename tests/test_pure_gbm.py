"""Pure GBM fit test (WP4): synthetic exact-linear truth validates the fitting code only. Run: python tests/test_pure_gbm.py"""
import sys
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(HERE / "src"))
from fuelblend.models import pure_gbm as pgm

def test_pure_gbm_recovers_linear_truth():
    rng = np.random.default_rng(7)
    n = 60
    mwt = rng.uniform(50.0, 200.0, n)
    logT = np.log(rng.uniform(250.0, 400.0, n))
    y = 2.0 * mwt + 0.05 * logT + 10.0
    X = np.column_stack([mwt, logT, np.ones(n)])
    model = pgm.make_model(random_state=11)
    pgm.fit_pure_gbm(model, X, y)
    pred = pgm.predict_pure_gbm(model, X)
    rmse = float(np.sqrt(np.mean((y - pred) ** 2)))
    # fixture truth is linear in (mWt, logT): GBM should recover it to well under 5 % relative error
    assert rmse < 0.05 * float(np.mean(np.abs(y))), rmse

if __name__ == "__main__":
    test_pure_gbm_recovers_linear_truth()
    print("ALL PURE GBM TESTS PASSED")
