"""G1 / G4: ML-constant substitution on lab mixtures. Usage: substitution_test.py --set dev|locked [--props a,b]."""
import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
import pandas as pd  # noqa: E402

from fuelblend.v8.substitution import TOL, score_substitution, summarize  # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument("--set", default="dev", choices=["dev", "locked"])
ap.add_argument("--props", default=",".join(TOL))
a = ap.parse_args()
out = ROOT / "results/v8/substitution"
out.mkdir(parents=True, exist_ok=True)
rows, summ = [], []
for p in a.props.split(","):
    df = score_substitution(p, a.set)
    rows.append(df.assign(prop=p, set=a.set))
    summ.append(dict(summarize(df, p), set=a.set))
    print(p, summ[-1], flush=True)
suffix = "" if a.set == "dev" else "_locked"
pd.concat(rows).to_csv(out / f"rows{suffix}.csv", index=False)
pd.DataFrame(summ).to_csv(out / f"summary{suffix}.csv", index=False)
