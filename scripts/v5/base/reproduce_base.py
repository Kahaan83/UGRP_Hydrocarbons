"""v5 BASE: rebuild the scoreboard, run the BASE tests, and diff results/v5/BASE/** against SHA256SUMS_v5_base.txt.

    python scripts/v5/base/reproduce_base.py                  # fast stages, tests, drift check; exit 1 on drift
    python scripts/v5/base/reproduce_base.py --full           # also re-run the slow prediction stages (hours)
    python scripts/v5/base/reproduce_base.py --update-manifest   # the ONLY manifest writer

Fast stages (default) rebuild from the committed rows_*.csv: cells_t2_bulk.py --from-rows (K1 offsets refitted, no new predictions),
cells_t3_bulk.py, cells_carried.py, build_scoreboard.py.  --full adds build_t2_targets, cells_t1_lab, cells_t2_flash and the full
cells_t2_bulk prediction pass (the last alone took over 80 minutes on the development machine).
"""
from __future__ import annotations

import importlib.util
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
PY = sys.executable
MANIFEST = ROOT / "SHA256SUMS_v5_base.txt"
BASE = ROOT / "results" / "v5" / "BASE"

_spec = importlib.util.spec_from_file_location("reproduce_all_v4", ROOT / "scripts/v4/reproduce_all.py")
ra = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(ra)

SLOW = [["scripts/v5/base/build_t2_targets.py"], ["scripts/v5/base/cells_t1_lab.py"], ["scripts/v5/base/cells_t2_flash.py"], ["scripts/v5/base/cells_t2_bulk.py"]]
FAST = [["scripts/v5/base/cells_t2_bulk.py", "--from-rows"], ["scripts/v5/base/cells_t3_bulk.py"], ["scripts/v5/base/cells_carried.py"],
        ["scripts/v5/base/build_scoreboard.py"]]


def manifest_paths() -> list[Path]:
    return sorted(p for p in BASE.rglob("*") if p.is_file())


def run(cmd: list[str]) -> None:
    print("[RUN]", " ".join(cmd), flush=True)
    r = subprocess.run([PY, *cmd], cwd=ROOT)
    if r.returncode:
        print(f"[FAIL] {' '.join(cmd)} exited {r.returncode}")
        sys.exit(1)


def main() -> int:
    update = "--update-manifest" in sys.argv
    stages = SLOW + FAST[1:] if "--full" in sys.argv else FAST   # full cells_t2_bulk replaces the --from-rows pass
    for s in stages:
        run(s)
    run(["-m", "pytest", *sorted(str(p.relative_to(ROOT)) for p in (ROOT / "tests").glob("test_base_*.py")), "-q", "-p", "no:cacheprovider"])
    if update:
        ra.write_manifest(MANIFEST, root=ROOT, paths=manifest_paths())
        return 0
    return 1 if ra.check_manifest(MANIFEST, root=ROOT, paths=manifest_paths()) else 0


if __name__ == "__main__":
    sys.exit(main())
