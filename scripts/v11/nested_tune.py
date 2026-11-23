"""v11b: nested-CV Optuna tuning of the viscosity / density pair models. Protocol: docs/v11b_predeclaration.md.

Usage: python scripts/v11/nested_tune.py visc|dens   then   python scripts/v11/nested_tune.py combine
"""
import sys, json, warnings, importlib.util
from pathlib import Path
import numpy as np, pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
warnings.simplefilter("ignore")
import lightgbm as lgb  # noqa: E402
import optuna  # noqa: E402
from sklearn.model_selection import GroupKFold  # noqa: E402

optuna.logging.set_verbosity(optuna.logging.WARNING)
OUT = ROOT / "results/v11/nested"
OUT.mkdir(parents=True, exist_ok=True)
N_TRIALS, SEED, NJ = 60, 11, 4
THRESH = {"visc": 0.002, "dens": 0.04}


def load(name, rel):
    spec = importlib.util.spec_from_file_location(name, ROOT / rel)
    m = importlib.util.module_from_spec(spec)
    sys.modules[name] = m
    spec.loader.exec_module(m)
    return m


class Setup:
    def __init__(self, which):
        self.which = which
        if which == "visc":
            from fuelblend.ml import viscosity_pair as mod
            self.D = load("dmv", "scripts/v5/base/delta_ml_viscosity.py")
            ctx = mod.context()
            self.R = self.D.regenerate_rows(ctx)
            self.P = self.D.build_pairs(self.R, ctx)
            self.y = (self.R.ln_meas - self.R.ln_pred).to_numpy()
            self.base_err = (self.R.ln_pred - self.R.ln_meas).to_numpy()
            self.mask_primary = (self.R.role.to_numpy() == "dev") & (self.R.T_K.to_numpy() >= self.D.T_WARM)
        else:
            from fuelblend.ml import density_pair as mod
            self.D = load("dmd", "scripts/v5/base/delta_ml_density.py")
            self.R, desc = self.D.regenerate()
            self.P = self.D.build_pairs(self.R, desc)
            self.y = (self.R.v_meas - self.R.v_ideal).to_numpy()
            self.base_err = (self.R.pred - self.R.meas).to_numpy()
            self.mask_primary = self.R.role.to_numpy() == "dev"
        self.mod = mod
        self.wins, self.xmin = self.D.WINS, self.D.XMIN
        self.sysarr = self.R.system_id.to_numpy()
        self.dev_rows = np.where(self.R.role.to_numpy() == "dev")[0]
        self.default = dict(mod.LGB_PARAMS)

    def err_from(self, corr):
        if self.which == "visc":
            return self.base_err + corr
        R = self.R
        return R.mm.to_numpy() / (R.v_ideal.to_numpy() + corr) * 1000.0 - R.meas.to_numpy()

    def fit_predict(self, params, p, train_sys, test_rows):
        """Corrections for test_rows (array of row indices), model trained on binary pair rows of train_sys."""
        Ptr = self.D.train_set(self.R, self.P, train_sys)
        ridx = Ptr.row.to_numpy()
        xx = Ptr.xx.to_numpy()
        g = np.clip(self.y[ridx], -self.wins, self.wins) / xx
        w = xx ** p
        tot = pd.Series(w).groupby(self.sysarr[ridx]).transform("sum").to_numpy()
        m = lgb.LGBMRegressor(**params)
        m.fit(Ptr[self.mod.FEATURES], g, sample_weight=w / tot, categorical_feature=["class_pair"])
        Pte = self.P[self.P.row.isin(test_rows)]
        corr = np.bincount(Pte.row.to_numpy(), weights=Pte.xx.to_numpy() * m.booster_.predict(Pte[self.mod.FEATURES]),
                           minlength=len(self.R))
        return corr

    def macro(self, err, rows):
        mk = np.zeros(len(self.R), bool)
        mk[rows] = True
        mk &= self.mask_primary
        return self.D.macro(self.R, err, mk)


def params_from(trial):
    loss = trial.suggest_categorical("loss", ["regression", "huber", "l1"])
    prm = dict(learning_rate=0.03, n_estimators=trial.suggest_int("n_estimators", 100, 1200, log=True),
               num_leaves=trial.suggest_int("num_leaves", 4, 48), min_child_samples=trial.suggest_int("min_child_samples", 10, 100),
               colsample_bytree=trial.suggest_float("colsample_bytree", 0.4, 1.0), subsample=trial.suggest_float("subsample", 0.6, 1.0),
               subsample_freq=1, reg_lambda=trial.suggest_float("reg_lambda", 0.1, 50, log=True), objective=loss,
               random_state=42, deterministic=True, force_row_wise=True, n_jobs=NJ, verbose=-1)
    if loss == "huber":
        prm["alpha"] = trial.suggest_float("huber_delta", 0.3, 3.0, log=True)
    return prm, trial.suggest_categorical("p", [0, 1, 2])


def tune_fold(S, tr_rows):
    inner = list(GroupKFold(n_splits=3).split(tr_rows, groups=S.sysarr[tr_rows]))

    def obj(trial):
        prm, p = params_from(trial)
        errs = []
        for a, b in inner:
            te = tr_rows[b]
            corr = S.fit_predict(prm, p, set(S.sysarr[tr_rows[a]]), te)
            errs.append(S.macro(S.err_from(corr), te))
        return float(pd.concat(errs).mean())

    st = optuna.create_study(sampler=optuna.samplers.TPESampler(seed=SEED), direction="minimize")
    st.optimize(obj, n_trials=N_TRIALS)
    bt = st.best_trial
    class _T:  # rebuild the best params
        def __init__(s, pr): s.pr = pr
        def suggest_categorical(s, n, c): return s.pr[n]
        def suggest_int(s, n, a, b, log=False): return s.pr[n]
        def suggest_float(s, n, a, b, log=False): return s.pr[n]
    return params_from(_T(bt.params)), bt.params, bt.value


def margins(S, err):
    out = {}
    R = S.R
    dev = R.role.to_numpy() == "dev"
    for b in ("T_lt_263K", "T_263_283K", "T_ge_283K"):
        msk = dev & (R.t_bin.to_numpy() == b)
        d = pd.DataFrame({"s": S.sysarr[msk], "e": np.abs(err[msk])}).groupby("s").e.max()
        out[b] = S.D.conformal(d.to_numpy())[0] if len(d) >= 9 else float("inf")
    return out


def run(which):
    S = Setup(which)
    n = len(S.R)
    cd, ct = np.zeros(n), np.zeros(n)
    folds = []
    for f, (a, b) in enumerate(GroupKFold(n_splits=5).split(S.dev_rows, groups=S.sysarr[S.dev_rows])):
        tr, te = S.dev_rows[a], S.dev_rows[b]
        trs = set(S.sysarr[tr])
        (prm, p), raw, inner_val = tune_fold(S, tr)
        cd[te] = S.fit_predict({**S.default, "n_jobs": NJ}, 2, trs, te)[te]
        ct[te] = S.fit_predict(prm, p, trs, te)[te]
        folds.append(dict(fold=f, inner_macro=inner_val, params=raw))
        print("outer fold", f, "inner", round(inner_val, 5), raw, flush=True)
    ed, et = S.err_from(cd), S.err_from(ct)
    a, b = S.macro(ed, S.dev_rows), S.macro(et, S.dev_rows)
    gain = (a - b.reindex(a.index)).to_numpy()
    rng = np.random.default_rng(42)
    bs = np.array([gain[rng.integers(0, len(gain), len(gain))].mean() for _ in range(2000)])
    res = dict(model=which, n_systems=int(len(gain)), mae_v5=float(a.mean()), mae_tuned=float(b.mean()), mean_gain=float(gain.mean()),
               ci95=[float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5))], p_one_sided=float((np.sum(bs <= 0) + 1) / 2001),
               threshold=THRESH[which], folds=folds)
    if which == "visc":
        res["margins_v5"], res["margins_tuned"] = margins(S, ed), margins(S, et)
        res["margin_guard_ok"] = all(res["margins_tuned"][k] <= res["margins_v5"][k] + 1e-12 for k in res["margins_v5"])
    else:
        tern = (S.R.role.to_numpy() == "dev") & (S.R.n_comp.to_numpy() >= 3)
        res["dev_ternary_rows"] = int(tern.sum())
        res["dev_ternary_mae"] = dict(ideal=float(np.abs(S.base_err[tern]).mean()), v5=float(np.abs(ed[tern]).mean()), tuned=float(np.abs(et[tern]).mean()))
    (OUT / f"{which}_summary.json").write_text(json.dumps(res, indent=1, default=float))
    print(json.dumps({k: v for k, v in res.items() if k != "folds"}, indent=1, default=float))


def combine():
    r = {k: json.load(open(OUT / f"{k}_summary.json")) for k in ("visc", "dens")}
    order = sorted(r, key=lambda k: r[k]["p_one_sided"])
    out = {}
    for i, k in enumerate(order):
        padj = min(1.0, r[k]["p_one_sided"] * (2 - i))
        ok = r[k]["mean_gain"] >= r[k]["threshold"] and padj < 0.05
        if k == "visc":
            ok = ok and r[k]["margin_guard_ok"]
        out[k] = dict(mean_gain=r[k]["mean_gain"], p_holm=padj, adopt=bool(ok))
    (OUT / "combined.json").write_text(json.dumps(out, indent=1))
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    {"visc": lambda: run("visc"), "dens": lambda: run("dens"), "combine": combine}[sys.argv[1]]()
