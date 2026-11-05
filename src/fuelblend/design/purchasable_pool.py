"""Purchasable Single-Compound Pool Module (Phase 3 Task U3).

Loads curated purchasable compounds from `data/processed/v3/purchasable_compounds.csv`
and provides verified property maps, commercial isomer flags, and single-isomer filtering.
"""
from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Optional
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
POOL_CSV = ROOT / "data" / "processed" / "v3" / "purchasable_compounds.csv"


def load_purchasable_pool(only_single_isomers: bool = False) -> pd.DataFrame:
    """Load purchasable compounds catalog."""
    df = pd.read_csv(POOL_CSV)
    if only_single_isomers:
        df = df[df["is_single_isomer"] == True].copy()
    return df


def get_compound_registry(only_single_isomers: bool = False) -> Dict[str, Dict[str, Any]]:
    """Return dictionary of compound_name -> property dict."""
    df = load_purchasable_pool(only_single_isomers=only_single_isomers)
    registry = {}
    for _, row in df.iterrows():
        name = str(row["compound_name"])
        registry[name] = row.to_dict()
    return registry


def is_blend_purchasable_single_isomers(components: List[str]) -> bool:
    """Check if all components in a blend are verified single chemical isomers."""
    registry = get_compound_registry(only_single_isomers=True)
    return all(c in registry for c in components)
