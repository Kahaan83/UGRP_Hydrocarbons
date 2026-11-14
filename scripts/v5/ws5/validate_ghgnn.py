"""WS5 step 0: does the wrapper reproduce the GH-GNN repo's own stored predictions (hydrocarbon-only rows)?
Writes results/v5/BASE/delta_ws5_ghgnn/validation_repo_predictions.json.  Uses only files shipped in the cloned repo."""
import json
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
warnings.simplefilter("ignore")
import ghgnn_lib as L  # noqa: E402

OUT = L.ROOT / "results" / "v5" / "BASE" / "delta_ws5_ghgnn"
OUT.mkdir(parents=True, exist_ok=True)
base = L.REPO / "models" / "temperature_dependency" / "continuous_extrapolation"
d = pd.concat([pd.read_csv(base / f).assign(src=f) for f in ("interpolation_predictions.csv", "lower_extrapolation_predictions.csv", "upper_extrapolation_predictions.csv")])


def ishc(s):
    return {a.GetSymbol() for a in L.Chem.MolFromSmiles(s).GetAtoms()} <= {"C"}


d = d[d.Solute_SMILES.map(ishc) & d.Solvent_SMILES.map(ishc)]
n_all = len(d)
d = d.head(400)
p = L.predict_ln_gamma_inf(d.Solute_SMILES.tolist(), d.Solvent_SMILES.tolist(), d["T"].tolist())
diff = np.abs(p - d["GH-GNN"].to_numpy())
res = dict(n_hc_rows_available=int(n_all), n_checked=int(len(d)), max_abs_diff=float(diff.max()), mean_abs_diff=float(diff.mean()),
           note="stored column GH-GNN in the repo's *_predictions.csv (ln gamma_inf); wrapper uses the repo model/features/loader unchanged")
print(res)
(OUT / "validation_repo_predictions.json").write_text(json.dumps(res, indent=2))
for s, v in [("CCCCCCCCCC", "CCCCCCCCCCCC"), ("c1ccccc1", "CCCCCCCCCCCC"), ("CCCCCCCCCCCC", "c1ccccc1"), ("Cc1cc(C)cc(C)c1", "CCCCCCCCCC")]:
    print(s, v, L.predict_ln_gamma_inf([s], [v], [25.0])[0], L.ad_pair(s, v))
