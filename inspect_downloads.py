#!/usr/bin/env python
"""Unpack the Zenodo zips and write a compact inventory of every table inside.

    python inspect_downloads.py                 # scans data/raw/zenodo/**/*.zip
    python inspect_downloads.py some_folder     # any folder of zips / csv / xlsx
    python inspect_downloads.py --splits        # also list every split file (default: summarised)

Output: data/interim/<label>/<zip name>/...   data/interim/INVENTORY.md   data/interim/read_errors.txt

What it ignores on purpose
  * macOS junk inside zips: __MACOSX/, ._* (AppleDouble binary shadows), .DS_Store
    (these have .csv/.json extensions but are not tables -> caused hundreds of parse errors)
  * model weights / logs (.pt, tfevents, .log) - counted, not parsed
  * `splits/` folders: summarised per folder (they only hold row-index/train-test labels)
"""
from __future__ import annotations
import collections, json, re, shutil, sys, zipfile
from pathlib import Path

import pandas as pd

args = [a for a in sys.argv[1:] if not a.startswith("--")]
SHOW_SPLITS = "--splits" in sys.argv
ROOT = Path(args[0]) if args else Path("data/raw/zenodo")
INTERIM = Path("data/interim")
KEY = re.compile(r"flash|visc|dcn|cetane|mon|ron|pour|cloud|density|vapou?r", re.I)
TABLE = {".csv", ".tsv", ".xlsx", ".xls", ".parquet", ".json"}
errors: collections.Counter = collections.Counter()
error_examples: dict[str, list[str]] = collections.defaultdict(list)


def is_junk(p: str | Path) -> bool:
    parts = Path(str(p).replace("\\", "/")).parts
    name = parts[-1] if parts else ""
    return "__MACOSX" in parts or name.startswith("._") or name.startswith("~$") or name == ".DS_Store" or name == "Thumbs.db"


def log_error(p: Path, e: Exception):
    kind = type(e).__name__
    errors[kind] += 1
    if len(error_examples[kind]) < 15:
        error_examples[kind].append(f"{p}: {str(e)[:150]}")


def purge_junk():
    n = 0
    for p in list(INTERIM.rglob("*")):
        if p.exists() and is_junk(p.relative_to(INTERIM)):
            shutil.rmtree(p, ignore_errors=True) if p.is_dir() else p.unlink(missing_ok=True)
            n += 1
    if n:
        print(f"removed {n} macOS junk entries from earlier extraction")


def unzip_all():
    for z in sorted(ROOT.rglob("*.zip")):
        dest = INTERIM / z.parent.name / z.stem
        marker = dest / ".extracted"
        if marker.exists():
            continue
        dest.mkdir(parents=True, exist_ok=True)
        print(f"unzipping {z.name} -> {dest}")
        skipped = failed = 0
        with zipfile.ZipFile(z) as f:
            for m in f.infolist():
                if is_junk(m.filename):
                    skipped += 1
                    continue
                try:
                    f.extract(m, dest)
                except Exception as e:  # noqa: BLE001  (e.g. Windows long-path)
                    failed += 1
                    log_error(Path(m.filename), e)
        marker.write_text("ok")
        print(f"   skipped {skipped} junk entries, {failed} extraction failures")


def count_lines(p: Path) -> int:
    n = 0
    with p.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 22), b""):
            n += chunk.count(b"\n")
    return max(n - 1, 0)


def describe(p: Path) -> str | None:
    """Return markdown for one table file, or None if unreadable (error logged)."""
    try:
        suf = p.suffix.lower()
        if suf in (".csv", ".tsv"):
            head = pd.read_csv(p, sep="\t" if suf == ".tsv" else ",", nrows=3, low_memory=False,
                               encoding_errors="replace")
            rows = count_lines(p)
            sheets = ""
        elif suf in (".xlsx", ".xls"):
            xl = pd.ExcelFile(p)
            head = xl.parse(xl.sheet_names[0], nrows=3)
            rows = len(xl.parse(xl.sheet_names[0]))
            sheets = f"\n- sheets: {xl.sheet_names}"
        elif suf == ".parquet":
            df = pd.read_parquet(p)
            head, rows, sheets = df.head(3), len(df), ""
        else:  # json: describe structure instead of forcing a table
            obj = json.loads(p.read_text(encoding="utf-8", errors="replace"))
            if isinstance(obj, dict):
                return f"- json object, keys: {list(obj)[:25]}"
            if isinstance(obj, list) and obj and isinstance(obj[0], dict):
                head = pd.DataFrame(obj[:3])
                rows, sheets = len(obj), ""
            else:
                return f"- json {type(obj).__name__}, length {len(obj) if hasattr(obj, '__len__') else 1}"
        return (f"- rows: {rows}, cols: {len(head.columns)}{sheets}\n- columns: {list(head.columns)[:40]}\n"
                f"```\n{head.head(2).to_string(max_colwidth=36)}\n```")
    except Exception as e:  # noqa: BLE001
        log_error(p, e)
        return None


def priority(p: Path):
    s = str(p).lower().replace("\\", "/")
    stage = 0 if "processed_data" in s else 1 if "raw_data" in s else 3 if "rejected" in s else 2
    return (not KEY.search(p.name), stage, s)


def main():
    INTERIM.mkdir(parents=True, exist_ok=True)
    purge_junk()
    unzip_all()
    all_files = [p for p in INTERIM.rglob("*") if p.is_file() and not is_junk(p.relative_to(INTERIM))
                 and p.name != ".extracted"]
    tables = [p for p in all_files if p.suffix.lower() in TABLE]
    other = [p for p in all_files if p.suffix.lower() not in TABLE]
    split_tables = [p for p in tables if "splits" in p.relative_to(INTERIM).parts]
    main_tables = sorted((p for p in tables if p not in set(split_tables)), key=priority)

    lines = ["# Inventory\n", f"- table files: {len(tables)} ({len(split_tables)} in splits/)",
             f"- other files (weights, logs, ...): {len(other)}\n"]
    n_ok = 0
    for p in main_tables:
        md = describe(p)
        if md:
            n_ok += 1
            lines.append(f"\n## {p.relative_to(INTERIM)}\n{md}")

    # splits: one schema example + counts per folder (they are only row-index / train-test labels)
    by_dir = collections.defaultdict(list)
    for p in split_tables:
        by_dir[p.parent].append(p)
    lines.append("\n\n# Splits (summary)\n")
    for d in sorted(by_dir):
        files = by_dir[d]
        lines.append(f"- `{d.relative_to(INTERIM)}`: {len(files)} files, e.g. {files[0].name}")
    if by_dir:
        md = describe(next(iter(by_dir.values()))[0])
        lines.append(f"\nExample split file schema:\n{md}")
    if SHOW_SPLITS:
        for p in split_tables:
            md = describe(p)
            if md:
                lines.append(f"\n## {p.relative_to(INTERIM)}\n{md}")

    ext = collections.Counter(p.suffix.lower() or "(none)" for p in other)
    lines.append(f"\n\n# Other files by extension\n{dict(ext)}")
    (INTERIM / "INVENTORY.md").write_text("\n".join(lines), encoding="utf-8")

    err_txt = ["Read/extract errors by type\n"]
    for k, v in errors.most_common():
        err_txt.append(f"{k}: {v}")
        err_txt += [f"    {x}" for x in error_examples[k]]
    (INTERIM / "read_errors.txt").write_text("\n".join(err_txt), encoding="utf-8")

    print(f"\ntables described: {n_ok}/{len(main_tables)} main, {len(split_tables)} split files summarised")
    print(f"errors: {sum(errors.values())} {dict(errors)}  (details: data/interim/read_errors.txt)")
    print("wrote data/interim/INVENTORY.md")


if __name__ == "__main__":
    main()
