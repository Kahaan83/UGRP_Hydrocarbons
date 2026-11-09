"""v4 Task 2: run the InChIKey-vs-SMILES check over every composition table and registry.

Writes results/v4/P1/identity_mismatches.csv (one row per flagged row; `table`, `row_id`
(0-based row position in the file), `system_id`/`fuel_id`, smiles, stored and recomputed key,
`block1_match`, `kind`) and prints a per-table summary. Usage: python scripts/v4/check_identities.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from fuelblend.data.identities import inchikey_mismatches  # noqa: E402

COMPOSITIONS = [
    "data/real_fuels/compositions.csv",
    "data/real_fuels/compositions_measured.csv",
    "data/real_fuels/compositions_reconstructed.csv",
    "data/real_fuels/v3/compositions_v3.csv",
]
# (path, smiles column, key column). Only registries that carry both columns can be checked.
REGISTRIES = [
    ("data/unified/compounds.csv", "smiles", "inchikey"),
    ("data/compounds_canonical.csv", "smiles_canonical", "inchikey"),
    ("data/processed/v2/compounds_anchored_vp.csv", "smiles", "inchikey"),
    ("extras/fuelblend_data/data/compounds_canonical.csv", "smiles_canonical", "inchikey"),
    ("compounds.csv", "smiles", "inchikey"),                      # name,cas,smiles -> no key
    ("extras/fuelblend_data/compounds.csv", "smiles", "inchikey"),  # name,cas,smiles -> no key
]
OUT = ROOT / "results/v4/P1/identity_mismatches.csv"


def _run(rel: str, smi: str, key: str, id_cols: list[str]):
    p = ROOT / rel
    if not p.exists():
        return None, "missing"
    df = pd.read_csv(p)
    if smi not in df.columns or key not in df.columns:
        return None, f"no {key} column (cannot check)"
    ids = [c for c in id_cols if c in df.columns]
    out = inchikey_mismatches(df, smi, key, id_cols=ids)
    out.insert(0, "table", rel)
    out = out.rename(columns={"row_index": "row_id"})
    return (out, len(df)), "ok"


def main() -> int:
    frames, summary = [], []
    for rel in COMPOSITIONS:
        res, st = _run(rel, "smiles", "inchikey", ["system_id", "fuel_id"])
        if res is not None:
            frames.append(res[0]); summary.append((rel, res[1], res[0]))
        else:
            summary.append((rel, 0, st))
    for rel, smi, key in REGISTRIES:
        res, st = _run(rel, smi, key, ["name", "cas"])
        if res is not None:
            frames.append(res[0]); summary.append((rel, res[1], res[0]))
        else:
            summary.append((rel, 0, st))
    allr = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
    OUT.parent.mkdir(parents=True, exist_ok=True)
    allr.to_csv(OUT, index=False)
    print(f"{'table':58s} {'rows':>6s} {'flagged':>8s}  by kind")
    for rel, n, r in summary:
        if isinstance(r, str):
            print(f"{rel:58s} {'-':>6s} {'-':>8s}  {r}")
        else:
            print(f"{rel:58s} {n:6d} {len(r):8d}  {r['kind'].value_counts().to_dict()}")
    print(f"\nwrote {OUT.relative_to(ROOT)}: {len(allr)} rows")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
