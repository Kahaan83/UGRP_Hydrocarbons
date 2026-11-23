import hashlib

from fuelblend.v10.lock import is_v10_locked, v10_lock_ids


def test_lock_rule_matches_v9_form():
    for i in range(20):
        g = f"fuel:x:{i}"
        assert is_v10_locked(g) == (int(hashlib.md5(g.encode()).hexdigest(), 16) % 4 == 0)


def test_lock_only_new():
    ids = [f"s{i}" for i in range(200)]
    existing = ids[:100]
    out = v10_lock_ids(ids, existing)
    assert not (out & set(existing))
    assert out == {i for i in ids[100:] if is_v10_locked(i)}
