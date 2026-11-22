"""Task 20: how many Leenhouts flash-mixture systems in the external staging file are unseen by every dataset and lock? Writes results/v7/lock/overlap.json."""
import sys, json
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "src"))
from fuelblend.v7 import guard  # noqa: E402

st = pd.read_csv(ROOT / "data/interim/external_staging/canonical_leenhouts_flash_mix_measurements.csv", usecols=["system_id", "n_comp", "kind", "property"])
sys_ids = set(st.system_id.astype(str))
dd = pd.read_csv(ROOT / "data/unified/measurements_dedup.csv", usecols=["system_id", "property", "kind"])
seen_all = set(dd.system_id.astype(str))
seen_flash = set(dd[dd.property == "flash_point_C"].system_id.astype(str))
out = dict(staging_rows=int(len(st)), staging_systems=len(sys_ids), by_n_comp={str(k): int(v) for k, v in st.groupby("n_comp").system_id.nunique().items()},
           in_dedup_any_property=len(sys_ids & seen_all), in_dedup_flash=len(sys_ids & seen_flash), locks={})
for name, spec in guard.LOCKS.items():
    if spec.id_kind == "system_id":
        out["locks"][name] = len(sys_ids & {guard.canon(i) for i in guard.locked_ids(name)})
    else:
        out["locks"][name] = f"not comparable (id kind {spec.id_kind})"
unseen = sys_ids - seen_all - {guard.canon(i) for n, s in guard.LOCKS.items() if s.id_kind == "system_id" for i in guard.locked_ids(n)}
out["unseen_systems"] = len(unseen)
out["unseen_rows"] = int(st.system_id.astype(str).isin(unseen).sum())
out["unseen_by_n_comp"] = {str(k): int(v) for k, v in st[st.system_id.astype(str).isin(unseen)].groupby("n_comp").system_id.nunique().items()}
(ROOT / "results/v7/lock/overlap.json").write_text(json.dumps(out, indent=1))
print(json.dumps(out, indent=1))
