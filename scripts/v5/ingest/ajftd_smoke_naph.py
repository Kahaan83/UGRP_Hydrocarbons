"""Write data/real_fuels/v5/ajftd_smoke_naph.csv and the R-SC-1 counts."""
import json
from pathlib import Path

from fuelblend.data.loaders_ajftd_v5 import load_ajftd_smoke_naph

REPO = Path(__file__).resolve().parents[3]
df = load_ajftd_smoke_naph()
df.to_csv(REPO / "data/real_fuels/v5/ajftd_smoke_naph.csv", index=False)
sm = df["smoke_point_mm"].notna()
ar = df["aromatics_vol_pct"].notna()
na = df["naphthalenes_vol_pct"].notna()
counts = {
    "fuels_in_table": int(len(df)),
    "fuels_with_smoke": int(sm.sum()),
    "fuels_censored_only": int((~sm & (df["smoke_qualifier"] != "")).sum()),
    "smoke_and_aromatics": int((sm & ar).sum()),
    "R-SC-1_smoke_aromatics_naphthalenes": int((sm & ar & na).sum()),
    "in_domain_aromatics_8_23_smoke_uncensored": int((sm & df["aromatics_vol_pct"].between(8, 23)).sum()),
    "source_group_counts_fuels_with_aromatics": {k: int(v) for k, v in df.loc[ar, "source_group"].value_counts().items()},
    "source_group_counts_fuels_with_smoke_and_aromatics": {k: int(v) for k, v in df.loc[sm & ar, "source_group"].value_counts().items()},
}
(REPO / "results/v5/BASE/delta_smoke_cert/r_sc_1_counts.json").write_text(json.dumps(counts, indent=2))
print(json.dumps(counts, indent=2))
