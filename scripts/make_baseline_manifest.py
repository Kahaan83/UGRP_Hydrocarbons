"""T0.1 - baseline freeze manifest + code/splits/results archive.

Writes:
  archive/baseline_manifest.json  (path, sha256, bytes for whitelisted files)
  archive/baseline_2026-11-03_code.zip (code, tests, scripts, configs, reports, splits, CSV/JSON results)
Run: .venv python scripts/make_baseline_manifest.py  (from repo root)
"""
import hashlib
import json
import zipfile
from datetime import datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "archive"
OUT.mkdir(exist_ok=True)

ZIP_DIRS = ["src", "tests", "scripts", "configs", "reports", "splits"]
ZIP_RESULTS_EXT = {".csv", ".json"}          # results/: CSVs + JSONs only (skip *.pt)
DATA_FILES = [                                # key data files: checksum only, not zipped
    "data/processed/model_table.csv",
    "data/processed/excluded_by_R1.csv",
    "data/processed/gcv_MJ_kg.csv",
    "data/processed/dcn_hc.csv",
    "data/processed/octane_ron_mon_hc.csv",
    "data/processed/freezing_point_pure_hc.csv",
    "data/processed/heat_of_combustion_real_fuels.csv",
    "data/unified/measurements_dedup.csv",
    "data/unified/measurements.csv",
    "splits/s1s2s3.json",
    "splits/locked.json",
]


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


files = []
zip_paths = []
for d in ZIP_DIRS:
    base = ROOT / d
    if not base.exists():
        continue
    for p in sorted(base.rglob("*")):
        if p.is_file() and "__pycache__" not in p.parts:
            zip_paths.append(p)
for p in sorted((ROOT / "results").rglob("*")):
    if p.is_file() and p.suffix in ZIP_RESULTS_EXT and "__pycache__" not in p.parts:
        zip_paths.append(p)
tracked = {p.relative_to(ROOT).as_posix() for p in zip_paths}

for rel in DATA_FILES:
    p = ROOT / rel
    if not p.exists():
        files.append({"path": rel, "status": "missing"})
        continue
    files.append({"path": rel, "sha256": sha256(p), "bytes": p.stat().st_size,
                  "in_zip": rel in tracked})

zip_path = OUT / "baseline_2026-11-03_code.zip"
with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as z:
    for p in zip_paths:
        z.write(p, p.relative_to(ROOT).as_posix())

manifest = {
    "tag": "baseline-2026-11-03",
    "created": datetime.now().isoformat(timespec="seconds"),
    "note": ("Full raw data (data/, debug/) is intentionally not copied: 4.9 GB + 2.1 GB. "
             "Key data files are checksummed above so state can be verified; code+configs+"
             "splits+results(CSV/JSON) are in the zip."),
    "zip": {"path": "archive/baseline_2026-11-03_code.zip",
            "sha256": sha256(zip_path), "bytes": zip_path.stat().st_size,
            "n_files": len(zip_paths)},
    "files": files,
}
(OUT / "baseline_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
n_missing = sum(f.get("status") == "missing" for f in files)
n_neg = sum(1 for f in files if f.get("bytes", 0) == 0)
print(f"zip: {len(zip_paths)} files, {zip_path.stat().st_size/1e6:.1f} MB")
print(f"data checksums: {len(files)} entries, {n_missing} missing, {n_neg} zero-byte")
