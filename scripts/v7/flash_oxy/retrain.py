"""v7 Task 11: flash oxygenate arm with MolPool heads, long-cap encoder and nested Optuna (GPU). See PREDECLARATION.md.
Usage: FUELBLEND_DEVICE=cuda python scripts/v7/flash_oxy/retrain.py [--quick]   (needs results/v7/mixture/encoder.pt)"""
import sys, json, time, os, warnings
warnings.filterwarnings("ignore")
sys.path.insert(0, "src"); sys.path.insert(0, "scripts/v7/mixture")
from pathlib import Path
import numpy as np, pandas as pd
from sklearn.model_selection import GroupKFold
from fuelblend.ml import molpool as MP
from fuelblend.v7 import guard, runner
import _mp

R = Path(__file__).resolve().parents[3]
QUICK = "--quick" in sys.argv
OUT = R / "results/v7/flash_oxy" / ("retrain_quick" if QUICK else "retrain")
OUT.mkdir(parents=True, exist_ok=True)
RUNS = Path(os.environ.get("FUELBLEND_RUNS", r"C:\fuelblend-gpu\runs"))
NT, TO = (2, 120) if QUICK else (60, 1800)
V6 = R / "results/v6/flash_oxy"
THETAS = [10, 20, 40, 80, 160, 1e4]   # as v6; 1e4 is linear molar
df = pd.read_csv(V6 / "rows.csv")
dev = df[df.role == "dev"].reset_index(drop=True)
guard.assert_dev_only(dev.system_id.unique(), locks=("flash_oxy_lock",))
ref = pd.read_csv(V6 / "dev_oof.csv")
assert len(ref) == len(dev) and np.allclose(ref.y, dev.y), "v6 OOF rows do not match"
kinds = json.load(open(V6 / "chosen_dev.json"))["base_kind"]
gate_comp = json.load(open(V6 / "readings.json"))["gate_comparator"]
enc = MP.load_encoder(R / "results/v7/mixture" / ("encoder_quick.pt" if QUICK else "encoder.pt"))
sm = [json.loads(s) for s in dev.smiles]
x = [np.array(json.loads(s), float) for s in dev.x]
pf = [np.array([np.nan if v is None else v for v in json.loads(s)], float) for s in dev.pure_flash]
allsmi = sorted({s for r in sm for s in r})
E = dict(zip(allsmi, MP.embed(enc, allsmi)))
H, X = MP.pad_sets([np.stack([E[s] for s in r]) for r in sm], x)
y, g, cls = dev.y.to_numpy(float), dev.system_id.to_numpy(), dev.chem_class.to_numpy()


def blend(theta):
    out = np.empty(len(y))
    for i, (xx, t) in enumerate(zip(x, pf)):
        h = ~np.isnan(t)
        w = xx[h] / xx[h].sum()
        out[i] = (w * t[h]).sum() if theta >= 1e4 else -theta * np.log((w * np.exp(-t[h] / theta)).sum())
    return out


BL = {th: blend(th) for th in THETAS}
lin = BL[1e4]
folds = np.full(len(y), -1)
for f, (_, te) in enumerate(GroupKFold(5).split(np.arange(len(y)), groups=g)):
    folds[te] = f
t0 = time.time()
pred = {"B": np.full(len(y), np.nan), "A": np.full(len(y), np.nan)}
thetas = {}
# arm B: all dev rows; the baseline is refit per outer fold (blend law) from the outer training rows only
r_b = np.zeros(len(y))
used = {}


def on_fold(k, tr):
    th = min(THETAS, key=lambda t: np.abs(BL[t][tr] - y[tr]).mean()) if kinds["B"] == "blend_law" else 1e4
    thetas[k] = th
    r_b[:] = y - BL[th]
    used[k] = BL[th].copy()


oof_b, res_b = _mp.tune(H, X, None, r_b, g, folds, prefix="v7_flash_B" + ("_quick" if QUICK else ""), out_dir=OUT / "armB", nt=NT, to=TO,
                        tb=RUNS / "flash_oxy" / "retrain" / "B", quick=QUICK, on_fold=on_fold)
for k in np.unique(folds):
    m = folds == k
    pred["B"][m] = used[int(k)][m] + oof_b[m]
print("arm B", f"{time.time() - t0:.0f}s", flush=True)
# arm A: hydrocarbon rows only, baseline linear molar
hc = np.where(cls == "hc")[0]
r_a = (y - lin)[hc]
oof_a, res_a = _mp.tune(H[hc], X[hc], None, r_a, g[hc], folds[hc], prefix="v7_flash_A" + ("_quick" if QUICK else ""), out_dir=OUT / "armA",
                        nt=NT, to=TO, tb=RUNS / "flash_oxy" / "retrain" / "A", quick=QUICK)
pred["A"][hc] = lin[hc] + oof_a
print("arm A", f"{time.time() - t0:.0f}s", flush=True)

o = dev[["system_id", "chem_class", "n_comp", "y"]].copy()
o["predB_v7"], o["predA_v7"] = pred["B"], pred["A"]
for c in ("predB_molpool", "predA_molpool", f"predB_{gate_comp}", "predB_lin_molar", "predB_blend_law"):
    o[c] = ref[c].to_numpy()
ox = cls == "oxy"
m3 = ox & (dev.n_comp.to_numpy() >= 3)
e7, ec, e1 = pred["B"] - y, ref[f"predB_{gate_comp}"].to_numpy() - y, ref.predB_molpool.to_numpy() - y


def safety(e):
    return dict(mean_over=float(np.clip(e, 0, None).mean()), p90_signed=float(np.percentile(e, 90)))


gc, lc, hc_, pc = _mp.boot_p(np.abs(ec[ox]) - np.abs(e7[ox]), g[ox])
g1, l1, h1, p1 = _mp.boot_p(np.abs(e1[ox]) - np.abs(e7[ox]), g[ox])
ph = _mp.holm(np.array([pc, p1]))
eb, ea = pred["B"][hc] - y[hc], pred["A"][hc] - y[hc]
dm, dlo, dhi, _ = _mp.boot_p(np.abs(eb) - np.abs(ea), g[hc])
sv7, sc = safety(e7[ox]), safety(ec[ox])
S = pd.DataFrame([dict(
    n_oxy=int(ox.sum()), systems_oxy=int(len(np.unique(g[ox]))), gate_comparator=gate_comp, mae_v7=np.abs(e7[ox]).mean(), mae_v1=np.abs(e1[ox]).mean(),
    mae_comparator=np.abs(ec[ox]).mean(), gain_vs_comp=gc, comp_lo=lc, comp_hi=hc_, p_holm_comp=ph[0], gain_vs_v1=g1, v1_lo=l1, v1_hi=h1, p_holm_v1=ph[1],
    beats_comparator=bool(ph[0] < 0.10 and gc > 0), beats_v1=bool(ph[1] < 0.10 and g1 > 0),
    n_multi=int(m3.sum()), mae_v7_multi=np.abs(e7[m3]).mean() if m3.any() else np.nan, mae_comp_multi=np.abs(ec[m3]).mean() if m3.any() else np.nan,
    safety_v7=json.dumps(sv7), safety_comp=json.dumps(sc), safety_v1=json.dumps(safety(e1[ox])),
    unsafe_tail=bool(sv7["p90_signed"] - sc["p90_signed"] > 1.0),
    hc_n=len(hc), hc_mae_armA=np.abs(ea).mean(), hc_mae_armB=np.abs(eb).mean(), hc_B_minus_A=dm, hc_lo=dlo, hc_hi=dhi, hc_noninferior=bool(dhi <= 1.0),
    thetas=json.dumps({int(k): v for k, v in thetas.items()}))])
runner.write_results(OUT, summary=S, oof=o, meta=dict(task="11_flash_oxy", n_trials=NT, timeout_s=TO, total_s=time.time() - t0,
                                                     best_params_B=res_b.best_params.tolist(), best_params_A=res_a.best_params.tolist()))
print(S.T.to_string())
