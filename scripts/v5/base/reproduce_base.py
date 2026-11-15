"""v5 BASE: rebuild the scoreboard, run the BASE tests, and diff results/v5/BASE/** against SHA256SUMS_v5_base.txt.

    python scripts/v5/base/reproduce_base.py                  # fast stages, tests, drift check; exit 1 on drift
    python scripts/v5/base/reproduce_base.py --full           # also re-run the slow prediction stages (hours)
    python scripts/v5/base/reproduce_base.py --update-manifest   # the ONLY manifest writer

Fast stages (default) rebuild from the committed rows_*.csv: cells_t2_bulk.py --from-rows (K1 offsets refitted, no new predictions),
cells_t3_bulk.py, cells_carried.py, build_scoreboard.py.  --full adds build_t2_targets, cells_t1_lab, cells_t2_flash and the full
cells_t2_bulk prediction pass (the last alone took over 80 minutes on the development machine).
"""
from __future__ import annotations

import hashlib
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
PY = sys.executable
MANIFEST = ROOT / "SHA256SUMS_v5_base.txt"
BASE = ROOT / "results" / "v5" / "BASE"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def compute_hashes(paths, root: Path) -> dict[str, str]:
    """{posix relative path: sha256} for every existing path."""
    return {p.relative_to(root).as_posix(): sha256(p) for p in paths if p.exists()}


def read_manifest(path: Path) -> dict[str, str]:
    out = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.strip():
            digest, rel = line.split(None, 1)
            out[rel.strip()] = digest
    return out


def check_manifest(manifest_path: Path, root: Path, paths) -> bool:
    """Read-only. Print drift vs the committed manifest; return True iff drift (or no manifest)."""
    if not manifest_path.exists():
        print(f"[MANIFEST] DRIFT: {manifest_path} does not exist (create it deliberately with --update-manifest)")
        return True
    manifest, current = read_manifest(manifest_path), compute_hashes(paths, root)
    changed = sorted(r for r in manifest if r in current and manifest[r] != current[r])
    missing = sorted(r for r in manifest if r not in current)
    new = sorted(r for r in current if r not in manifest)
    for tag, rels in (("CHANGED", changed), ("MISSING", missing), ("NEW", new)):
        for r in rels:
            print(f"[MANIFEST] {tag}: {r}")
    drift = bool(changed or missing or new)
    print(f"[MANIFEST] {'DRIFT' if drift else 'OK'}: {len(changed)} changed, {len(missing)} missing, {len(new)} new "
          f"({len(current)} files hashed)")
    return drift


def write_manifest(out_path: Path, root: Path, paths) -> int:
    """The ONLY manifest writer; called only for --update-manifest."""
    hashes = compute_hashes(paths, root)
    out_path.write_text("".join(f"{h}  {r}" + chr(10) for r, h in sorted(hashes.items())), encoding="utf-8", newline="\n")
    print(f"[MANIFEST] wrote {len(hashes)} hashes to {out_path.name}")
    return len(hashes)


SLOW = [["scripts/v5/base/build_t2_targets.py"], ["scripts/v5/base/cells_t1_lab.py"], ["scripts/v5/base/cells_t2_flash.py"], ["scripts/v5/base/cells_t2_bulk.py"]]
FAST = [["scripts/v5/base/cells_t2_bulk.py", "--from-rows"], ["scripts/v5/base/cells_t3_bulk.py"], ["scripts/v5/base/cells_carried.py"],
        ["scripts/v5/base/build_scoreboard.py"]]


def manifest_paths() -> list[Path]:
    return sorted(p for p in BASE.rglob("*") if p.is_file() and not p.relative_to(BASE).parts[0].startswith("delta_"))   # deltas are measured against the frozen files, not part of them


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
        write_manifest(MANIFEST, root=ROOT, paths=manifest_paths())
        return 0
    return 1 if check_manifest(MANIFEST, root=ROOT, paths=manifest_paths()) else 0


if __name__ == "__main__":
    sys.exit(main())
