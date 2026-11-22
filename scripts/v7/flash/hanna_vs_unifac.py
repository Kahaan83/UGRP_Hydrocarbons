"""Task 17: ideal vs UNIFAC vs HANNA gamma in the v5 Le Chatelier flash route, Leenhouts hydrocarbon binary dev rows."""
import sys, json, warnings
from pathlib import Path
import numpy as np, pandas as pd

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "src"))
warnings.filterwarnings("ignore")
from fuelblend.baselines import flash_unifac as FU  # noqa: E402
from fuelblend.baselines.flash import predict_flash_le_chatelier  # noqa: E402
from fuelblend.thermo.gamma_hanna import make_hanna_gamma_fn  # noqa: E402
from fuelblend.v7 import guard  # noqa: E402

OUT = ROOT / "results/v7/flash/hanna"
rng = np.random.default_rng(7)


def flash(smiles, x, gamma_fn):
    x, ant, lfl, missing = FU._prepare(smiles, x, None)
    if missing:
        return np.nan
    T, st = predict_flash_le_chatelier(list(x), ant, list(lfl), lfl_temp=True, gamma_fn=gamma_fn)
    return T if st == "ok" else np.nan


def boot(d, g, n=4000):
    ug = np.unique(g)
    idx = {k: np.where(g == k)[0] for k in ug}
    pm = np.array([d[idx[k]].mean() for k in ug])
    b = np.array([pm[rng.integers(0, len(ug), len(ug))].mean() for _ in range(n)])
    return float(pm.mean()), float(np.percentile(b, 5)), float(np.percentile(b, 95))


def main():
    d = pd.read_csv(ROOT / "data/unified/measurements_dedup.csv")
    b = d[(d.property == "flash_point_C") & (d.kind == "binary") & d.is_hc.astype(bool) & ~d.conflict.astype(bool)].copy()
    loc = set()
    for n, s in guard.LOCKS.items():
        if s.id_kind == "system_id":
            loc |= {guard.canon(i) for i in guard.locked_ids(n)}
    b = b[~b.system_id.isin(loc)].reset_index(drop=True)
    smi = pd.read_csv(ROOT / "data/unified/compounds.csv").set_index("inchikey").smiles
    rows = []
    for r in b.itertuples():
        iks = json.loads(r.comp_inchikeys)
        x = np.array(json.loads(r.x_mole), float)
        s = [smi[k] for k in iks]
        try:
            ideal = flash(s, x, None)
            ug, _ = FU.make_unifac_gamma_fn(s)
            unf = flash(s, x, ug)
            han = flash(s, x, make_hanna_gamma_fn(s))
        except Exception as e:  # noqa: BLE001
            ideal = unf = han = np.nan
        rows.append(dict(system_id=r.system_id, measured=r.value, ideal=ideal, unifac=unf, hanna=han))
    R = pd.DataFrame(rows)
    ok = R[["ideal", "unifac", "hanna"]].notna().all(axis=1)
    n_drop = int((~ok).sum())
    R = R[ok].reset_index(drop=True)
    R.to_csv(OUT / "predictions.csv", index=False, float_format="%.4f")
    g = R.system_id.to_numpy()
    e = {k: np.abs(R[k] - R.measured).to_numpy() for k in ("ideal", "unifac", "hanna")}
    o = {k: np.clip(R[k] - R.measured, 0, None).to_numpy() for k in e}
    summ = dict(n_rows=len(R), n_systems=int(R.system_id.nunique()), n_dropped=n_drop,
                mae={k: float(v.mean()) for k, v in e.items()}, mean_over={k: float(v.mean()) for k, v in o.items()})
    for comp in ("ideal", "unifac"):
        gm = boot(e[comp] - e["hanna"], g)
        ov = boot(o["hanna"] - o[comp], g)
        summ[f"hanna_vs_{comp}"] = dict(mae_gain=gm[0], mae_gain_lo=gm[1], mae_gain_hi=gm[2], over_increase=ov[0], over_increase_lo=ov[1], over_increase_hi=ov[2],
                                        gate_pass=bool(gm[1] > 0 and ov[2] < 0.5))
    (OUT / "summary.json").write_text(json.dumps(summ, indent=1))
    print(json.dumps(summ, indent=1))


if __name__ == "__main__":
    main()
