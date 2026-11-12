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


def test_delta_readings_pass():
    """R-API-1..3 pass in delta.csv; pooled after_bias equals the frozen lfl_overlay T2 cell bias (read from cells_t2_flash.csv)."""
    import pandas as pd
    base = Path(__file__).resolve().parents[1] / "results" / "v5" / "BASE"
    d = pd.read_csv(base / "delta_api_route" / "delta.csv").set_index("reading")
    for r in ("R-API-1", "R-API-2", "R-API-3"):
        assert str(d.loc[r, "pass"]) == "True", (r, d.loc[r].to_dict())
    c = pd.read_csv(base / "cells_t2_flash.csv")
    cell = c[(c.set == "pooled") & (c.variant == "lfl_overlay") & (c.model == "ideal_le_chatelier")]
    assert len(cell) == 1
    assert abs(float(d.loc["R-API-1", "after_bias"]) - float(cell.bias.iloc[0])) < 1e-6
