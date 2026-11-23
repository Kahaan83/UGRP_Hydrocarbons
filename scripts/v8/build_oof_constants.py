"""Write results/v8/oof_constants.csv (long) from the committed v7 Track B OOF predictions (no refit needed)."""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
import pandas as pd  # noqa: E402

from fuelblend.v8.ml_constants import PROPS, oof_constants  # noqa: E402

out = ROOT / "results/v8"
out.mkdir(parents=True, exist_ok=True)
df = pd.concat([oof_constants(p).assign(prop=p) for p in PROPS], ignore_index=True)
df.to_csv(out / "oof_constants.csv", index=False)
print(df.groupby("prop").size().to_string())
