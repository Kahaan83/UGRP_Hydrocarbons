"""Manifest drift check in scripts/v4/reproduce_all.py: read-only by default, --update-manifest is the only writer."""
import importlib.util
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
_spec = importlib.util.spec_from_file_location("reproduce_all_under_test", ROOT / "scripts/v4/reproduce_all.py")
ra = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(ra)


def _setup(tmp_path):
    (tmp_path / "a.csv").write_text("x\n1\n")
    (tmp_path / "b.csv").write_text("y\n2\n")
    paths = [tmp_path / "a.csv", tmp_path / "b.csv"]
    man = tmp_path / "SUMS.txt"
    ra.write_manifest(man, root=tmp_path, paths=paths)
    return paths, man


def test_untouched_passes(tmp_path):
    paths, man = _setup(tmp_path)
    assert ra.check_manifest(man, root=tmp_path, paths=paths) is False


def test_tampered_missing_new_are_drift_and_manifest_untouched(tmp_path, capsys):
    paths, man = _setup(tmp_path)
    before = man.read_text()
    paths[0].write_text("x\n999\n")
    assert ra.check_manifest(man, root=tmp_path, paths=paths) is True
    assert "CHANGED: a.csv" in capsys.readouterr().out
    paths[1].unlink()
    assert ra.check_manifest(man, root=tmp_path, paths=paths) is True
    assert "MISSING: b.csv" in capsys.readouterr().out
    (tmp_path / "c.csv").write_text("z\n")
    assert ra.check_manifest(man, root=tmp_path, paths=paths + [tmp_path / "c.csv"]) is True
    assert "NEW: c.csv" in capsys.readouterr().out
    assert man.read_text() == before  # check never writes


def test_missing_manifest_is_drift_and_not_created(tmp_path):
    paths, man = _setup(tmp_path)
    man.unlink()
    assert ra.check_manifest(man, root=tmp_path, paths=paths) is True
    assert not man.exists()


def _run(monkeypatch, *argv):
    monkeypatch.setattr(sys, "argv", ["x", *argv])
    try:
        ra.main()
    except SystemExit as e:
        return e.code
    return 0


def test_cli_only_update_manifest_writes(tmp_path, monkeypatch):
    man = tmp_path / "SUMS.txt"
    f = tmp_path / "a.csv"
    f.write_text("x\n1\n")
    monkeypatch.setattr(ra, "manifest_paths", lambda: [f])
    monkeypatch.setattr(ra, "ROOT", tmp_path)
    monkeypatch.setattr(ra, "verify_deliverables", lambda: True)
    assert _run(monkeypatch, "--verify-only", "--manifest", str(man)) == 1
    assert not man.exists()
    assert _run(monkeypatch, "--verify-only", "--update-manifest", "--manifest", str(man)) == 0
    assert man.exists()
    assert _run(monkeypatch, "--verify-only", "--manifest", str(man)) == 0
    f.write_text("tampered\n")
    assert _run(monkeypatch, "--verify-only", "--manifest", str(man)) == 1
