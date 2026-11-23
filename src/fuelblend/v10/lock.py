"""v10 lock rule: a deterministic quarter of the v10-new systems and fuels."""
from __future__ import annotations

import hashlib
from collections.abc import Iterable

RULE = "md5 mod 4 == 0"


def is_v10_locked(group_id: str) -> bool:
    return int(hashlib.md5(group_id.encode("utf-8")).hexdigest(), 16) % 4 == 0


def v10_lock_ids(new_ids: Iterable[str], existing_ids: Iterable[str]) -> frozenset[str]:
    """Locked subset of new_ids, never including an id already in existing_ids."""
    existing = set(existing_ids)
    return frozenset(i for i in new_ids if i not in existing and is_v10_locked(i))
