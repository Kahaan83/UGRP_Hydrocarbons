"""make_v2_manifest.py - Generates archive/v2_manifest.json with SHA-256 hashes and file sizes.
"""
from __future__ import annotations

import hashlib
import json
from datetime import datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def hash_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while chunk := f.read(65536):
            h.update(chunk)
    return h.hexdigest()


def build_manifest() -> None:
    manifest_path = ROOT / "archive" / "v2_manifest.json"
    manifest_path.parent.mkdir(parents=True, exist_ok=True)

    files_to_index: list[Path] = []

    # 1. results/v2/
    for p in sorted((ROOT / "results" / "v2").rglob("*")):
        if p.is_file() and not p.name.startswith("."):
            files_to_index.append(p)

    # 2. configs/
    for p in sorted((ROOT / "configs").glob("*.yaml")):
        if p.is_file():
            files_to_index.append(p)

    # 3. scripts/v2/
    for p in sorted((ROOT / "scripts" / "v2").glob("*.py")):
        if p.is_file():
            files_to_index.append(p)

    # 4. Root execution files
    for p in [ROOT / "scripts" / "reproduce_all_v2.sh", ROOT / "Makefile"]:
        if p.is_file():
            files_to_index.append(p)

    # 5. Key real fuels & processed data
    for p in sorted((ROOT / "data" / "real_fuels").glob("*.csv")):
        if p.is_file():
            files_to_index.append(p)
    if (ROOT / "splits" / "real_fuels_nested_v2.json").is_file():
        files_to_index.append(ROOT / "splits" / "real_fuels_nested_v2.json")
    for p in sorted((ROOT / "data" / "processed" / "v2").glob("*.csv")):
        if p.is_file():
            files_to_index.append(p)

    # 6. Core v2 documentation
    for p in [
        ROOT / "docs" / "FULL_REPORT_v2.md",
        ROOT / "docs" / "experiment_plan.md",
        ROOT / "docs" / "claims_register.md",
    ]:
        if p.is_file():
            files_to_index.append(p)

    entries = []
    for p in files_to_index:
        rel_path = p.relative_to(ROOT).as_posix()
        sha = hash_file(p)
        size = p.stat().st_size
        entries.append({
            "path": rel_path,
            "sha256": sha,
            "bytes": size,
        })

    manifest = {
        "tag": "v2-2026-11-05",
        "created": datetime.now().isoformat(),
        "description": "Checksum manifest for UGRP Hydrocarbons v2 verification freeze, deliverables, configs, scripts, and docs.",
        "n_files": len(entries),
        "files": entries,
    }

    with open(manifest_path, "w", encoding="utf-8") as f:
        json.dump(manifest, f, indent=2)

    print(f"Wrote v2 manifest to {manifest_path} ({len(entries)} files indexed).")


if __name__ == "__main__":
    build_manifest()
