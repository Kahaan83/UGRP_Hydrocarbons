from fuelblend.v9.lock import is_v9_locked, v9_lock_ids


def test_lock_rule_is_deterministic():
    assert is_v9_locked("fuel:ajftd:2") == is_v9_locked("fuel:ajftd:2")


def test_lock_share_near_quarter():
    ids = [f"fuel:x:{i}" for i in range(4000)]
    assert 0.22 < sum(map(is_v9_locked, ids)) / 4000 < 0.28


def test_lock_ids_subset():
    ids = [f"s{i}" for i in range(100)]
    assert v9_lock_ids(ids) == {i for i in ids if is_v9_locked(i)}
