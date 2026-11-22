"""v7 Task 13: CPU tree models with reproduction gate, curves and nested Optuna. See PREDECLARATION.md.
Usage: python scripts/v7/cpu/retrain_trees.py <model> [--quick]   model in viscosity_pair, density_pair, fusion_estimate, ysi_lgbm, freezing_step3, mixture_gbr"""
import sys, json, time, os, importlib, warnings
warnings.filterwarnings("ignore")
sys.path.insert(0, "src"); sys.path.insert(0, "scripts/v7/cpu")
from pathlib import Path
import pandas as pd
import _trees

R = Path(__file__).resolve().parents[3]
name = sys.argv[1]
QUICK = "--quick" in sys.argv
NT, TO = (2, 120) if QUICK else (60, 1800)
OUT = R / "results/v7/retrain_v5" / (name + ("_quick" if QUICK else ""))
OUT.mkdir(parents=True, exist_ok=True)
RUNS = Path(os.environ.get("FUELBLEND_RUNS", r"C:\fuelblend-gpu\runs"))
GUARD = ("src", "results/v5/BASE")
h0 = _trees.tree_hash(*GUARD)
t0 = time.time()
mod = importlib.import_module("m_" + name)
gate, summ = mod.run(NT, TO, QUICK, OUT, RUNS / "cpu" / name, time.time)
gate["total_s"] = time.time() - t0
(OUT / "repro_gate.json").write_text(json.dumps(gate, indent=2, default=float))
print("gate", json.dumps(gate, default=float), flush=True)
if summ is not None:
    bp = summ.pop("best_params")
    (OUT / "best_params.json").write_text(json.dumps(bp, indent=2))
    pd.DataFrame([summ]).to_csv(OUT / "summary.csv", index=False)
    cv = [pd.read_csv(p) for p in OUT.glob("curves_*.csv")]
    if cv:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        d = pd.concat(cv)
        fig, ax = plt.subplots(figsize=(6, 4))
        for f_, g_ in d.groupby("fold"):
            ax.plot(g_.groupby("epoch").val_loss.median(), label=f"fold {f_}")
        ax.set_yscale("log"); ax.set_xlabel("boosting round"); ax.set_ylabel("val l2"); ax.legend(fontsize=7)
        fig.tight_layout(); fig.savefig(OUT / "curves.png", dpi=120)
        d.to_csv(OUT / "curves.csv", index=False)
    print(pd.DataFrame([summ]).round(4).T.to_string())
assert _trees.tree_hash(*GUARD) == h0, "files under src/ or results/v5/BASE changed during the run"
