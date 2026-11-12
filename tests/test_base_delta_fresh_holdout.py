"""v5 BASE Task 15: fresh holdout is scored once."""
import importlib.util
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
_s = importlib.util.spec_from_file_location("dfh", ROOT / "scripts/v5/base/delta_fresh_holdout.py")
dfh = importlib.util.module_from_spec(_s)
_s.loader.exec_module(dfh)


def test_refuses_second_run(tmp_path):
    dfh.guard(tmp_path)                                   # first run allowed
    (tmp_path / "scored_once.json").write_text("{}")
    with pytest.raises(SystemExit):
        dfh.guard(tmp_path)


def test_committed_holdout_is_marked_scored():
    assert (dfh.OUT / "scored_once.json").exists()
    with pytest.raises(SystemExit):
        dfh.main()
