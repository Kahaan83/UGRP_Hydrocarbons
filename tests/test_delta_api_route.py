"""Pin the frozen provider grid used by the delta_api_route work."""
import hashlib
from pathlib import Path

import pytest

GRID_SHA256 = "b931e969926698dff52d69ff51dd99a990c24ff12eca545aa98e970afd3e7daf"
GRID = Path(__file__).resolve().parents[1] / "data" / "cache" / "pure_grid.npz"


def test_grid_pinned():
    if not GRID.exists():
        pytest.skip(f"{GRID} absent (gitignored binary)")
    h = hashlib.sha256()
    with open(GRID, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    assert h.hexdigest() == GRID_SHA256
