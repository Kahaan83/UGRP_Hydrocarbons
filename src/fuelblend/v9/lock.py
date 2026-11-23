"""v9 lock rule: a deterministic quarter of the NEW systems and fuels."""
from __future__ import annotations

import hashlib
from collections.abc import Iterable

RULE = "md5 mod 4 == 0"


def is_v9_locked(group_id: str) -> bool:
    return int(hashlib.md5(group_id.encode("utf-8")).hexdigest(), 16) % 4 == 0


def v9_lock_ids(new_ids: Iterable[str]) -> frozenset[str]:
    return frozenset(i for i in new_ids if is_v9_locked(i))
