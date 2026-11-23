"""Score the v7 fresh flash lock once: HANNA vs UNIFAC gamma in the v5 Le Chatelier flash route.

Pre-registered in results/v7/lock/LOCK.md and lock.json. Refuses to run if an output already exists.
Gate: MAE gain (UNIFAC minus HANNA) 90 % bootstrap CI lower bound > 0 and upper bound of the
mean over-prediction increase (HANNA minus UNIFAC) < 0.5 degC. Nothing is fitted.
"""
import sys, json, glob, datetime, warnings
from pathlib import Path
import numpy as np, pandas as pd

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
warnings.filterwarnings("ignore")
from fuelblend.baselines import flash_unifac as FU  # noqa: E402
from fuelblend.thermo.gamma_hanna import make_hanna_gamma_fn  # noqa: E402
from hanna_vs_unifac import flash, boot  # noqa: E402

LOCK = ROOT / "results/v7/lock"
OUT = LOCK / "scored"
STAGE = ROOT / "data/interim/external_staging"


def main():
    if OUT.exists():
        sys.exit("lock already scored; refusing to rescore")
    OUT.mkdir()
    ids = set(json.load(open(LOCK / "lock.json"))["locked_system_ids"])
    d = pd.read_csv(STAGE / "canonical_leenhouts_flash_mix_measurements.csv")
    d = d[d.system_id.isin(ids)].reset_index(drop=True)
    smi = {}
    for f in [ROOT / "data/unified/compounds.csv"] + [Path(p) for p in glob.glob(str(STAGE / "*compounds*.csv"))]:
        c = pd.read_csv(f)
        for k, v in zip(c.inchikey, c.smiles):
            if isinstance(v, str):
                smi.setdefault(k, v)
    rows = []
    for r in d.itertuples():
        s = [smi[k] for k in json.loads(r.comp_inchikeys)]
        x = np.array(json.loads(r.x_mole), float)
        try:
            ug, _ = FU.make_unifac_gamma_fn(s)
            unf = flash(s, x, ug)
            han = flash(s, x, make_hanna_gamma_fn(s))
            ide = flash(s, x, None)
        except Exception:  # noqa: BLE001
            ide = unf = han = np.nan
        rows.append(dict(system_id=r.system_id, measured=r.value, ideal=ide, unifac=unf, hanna=han))
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
    gm = boot(e["unifac"] - e["hanna"], g)
    ov = boot(o["hanna"] - o["unifac"], g)
    summ["primary_hanna_vs_unifac"] = dict(mae_gain=gm[0], mae_gain_lo=gm[1], mae_gain_hi=gm[2],
                                           over_increase=ov[0], over_increase_lo=ov[1], over_increase_hi=ov[2],
                                           gate_pass=bool(gm[1] > 0 and ov[2] < 0.5))
    (OUT / "summary.json").write_text(json.dumps(summ, indent=1))
    with open(ROOT / "reports/LOCKED_ACCESS.log", "a") as fh:
        fh.write(f"{datetime.datetime.now(datetime.timezone.utc).isoformat()} v7-flash-fresh lock scored once "
                 f"({summ['n_systems']} systems) by score_fresh_lock.py\n")
    print(json.dumps(summ, indent=1))


if __name__ == "__main__":
    main()
