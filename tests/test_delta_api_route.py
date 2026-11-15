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
    # delta_api_route was measured against the pre-adoption BASE.  After U4 (2026-11-15) the BASE cell moved -0.69151175 -> 0.096130515
    # (decision 3, D3_overlay); the committed after_bias is checked against the pre-adoption value recorded in drift.csv, and the adopted
    # cell against drift.csv `new` and against the adopted api_as_shipped pooled cell (the API still equals the lfl_overlay route).
    dr = pd.read_csv(base / "delta_adopt_rescore" / "drift.csv")
    key = ("property=flash_point_C;tier=T2_measured;set=pooled;model=ideal_le_chatelier;variant=lfl_overlay;label=descriptive|bias")
    row = dr[(dr.path == "results/v5/BASE/cells_t2_flash.csv") & (dr.cell == key)]
    assert len(row) == 1
    assert abs(float(d.loc["R-API-1", "after_bias"]) - float(row.old.iloc[0])) < 1e-6
    assert abs(float(cell.bias.iloc[0]) - float(row.new.iloc[0])) < 1e-6
    api = c[(c.set == "pooled") & (c.model == "api_as_shipped")]
    assert len(api) == 1 and abs(float(api.bias.iloc[0]) - float(cell.bias.iloc[0])) < 1e-6
