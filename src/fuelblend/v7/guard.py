"""Lock guard: v7 code must never score a row already scored in v6.

Every call names its locks explicitly; there is no global default.
Ids are compared in canonical form (known source prefixes stripped).
"""
from __future__ import annotations

import json
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
_STRIP = ("leenhouts:", "boehm:", "kaust:", "fuel:")


class LockedRowError(ValueError):
    """A v7 row overlaps a v6 lock."""


@dataclass(frozen=True)
class LockSpec:
    name: str
    path: str
    id_kind: str  # inchikey | system_id | group | fuel_id
    scope: str


LOCKS: dict[str, LockSpec] = {
    s.name: s
    for s in (
        LockSpec("ysi_holdout", "results/v6/gnn/ysi_holdout/lock.json", "inchikey", "ysi"),
        LockSpec("freezing_holdout", "results/v6/freezing/holdout_lock.json", "group", "freezing"),
        LockSpec("mixture_lock", "results/v6/mixture/lock.json", "system_id", "mixture heads"),
        LockSpec("mixture_holdout51", "results/v6/mixture/holdout51_pred.csv", "system_id", "mixture heads"),
        LockSpec("flash_oxy_lock", "results/v6/flash_oxy/lock.json", "system_id", "flash"),
        LockSpec("splits_locked", "splits/locked.json", "fuel_id", "all properties"),
    )
}


def canon(x: object) -> str:
    s = str(x)
    for p in _STRIP:
        if s.startswith(p):
            return s[len(p):]
    return s


def _json(path: str) -> dict:
    return json.loads((ROOT / path).read_text(encoding="utf-8"))


@lru_cache(maxsize=None)
def locked_ids(name: str) -> frozenset[str]:
    if name not in LOCKS:
        raise KeyError(f"unknown lock {name!r}; known: {sorted(LOCKS)}")
    path = LOCKS[name].path
    if name == "ysi_holdout":
        ids = _json(path)["locked_inchikeys"]
    elif name == "freezing_holdout":
        ids = _json(path)["holdout_groups"]
    elif name == "mixture_lock":
        ids = [i for h in _json(path)["heads"].values() for i in h["locked"]]
    elif name == "mixture_holdout51":
        ids = pd.read_csv(ROOT / path, usecols=["system_id"])["system_id"].unique()
    elif name == "flash_oxy_lock":
        ids = [i for c in _json(path)["classes"].values() for i in c["locked"]]
    else:  # splits_locked
        d = _json(path)
        ids = list(d["real_fuel_ids_locked"]) + list(d["multi_system_ids_locked"])
    return frozenset(canon(i) for i in ids)


def _overlap(ids, locks: tuple[str, ...]) -> dict[str, set[str]]:
    if not locks:
        raise ValueError("name the locks explicitly")
    cid = {canon(i) for i in ids}
    return {n: cid & locked_ids(n) for n in locks if cid & locked_ids(n)}


def assert_dev_only(ids, *, locks: tuple[str, ...]) -> None:
    hit = _overlap(ids, locks)
    if hit:
        n = next(iter(hit))
        raise LockedRowError(
            f"{sum(len(v) for v in hit.values())} locked id(s); lock {n!r} "
            f"({LOCKS[n].path}), e.g. {sorted(hit[n])[:3]}"
        )


def mask_locked(df: pd.DataFrame, id_col: str, *, locks: tuple[str, ...],
                label_cols: list[str] | None = None) -> tuple[pd.DataFrame, int]:
    """Drop locked rows, or with label_cols set only those labels to NaN.

    Returns (frame, number of locked rows)."""
    if not locks:
        raise ValueError("name the locks explicitly")
    cid = df[id_col].map(canon)
    bad = pd.Series(False, index=df.index)
    for n in locks:
        bad |= cid.isin(locked_ids(n))
    out = df.copy()
    if label_cols:
        for c in label_cols:
            out[c] = out[c].astype(float).where(~bad)
    else:
        out = out.loc[~bad]
    return out, int(bad.sum())
