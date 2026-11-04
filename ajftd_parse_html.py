#!/usr/bin/env python
"""Re-read the AJFTD detail pages your scraper already saved (debug/detail_<id>.html) - NO site access needed.

    python ajftd_parse_html.py                         # reads debug/detail_*.html
    python ajftd_parse_html.py path/to/html_folder
    python ajftd_parse_html.py --compare ajftd_data.csv    # also cross-check against your scraper's CSV

Why: the scraper only exports metrics whose label is in its alias lists, so any other metric on the page
(freezing point, heat of combustion, GC composition ...) was silently skipped, and the fuel name was taken
from the first heading ("Alternative Jet Fuels" for every record). The HTML is already on disk, so this
re-extracts EVERYTHING without touching the slow website again.

Outputs (written next to the html folder's parent, in data/raw/ajftd/):
  ajftd_all_metrics.csv       same columns as ajftd_data.csv, but EVERY metric button/table + real fuel name/label
  ajftd_metrics_inventory.csv metric label x (fuels, rows, rows with a real number): shows what exists on the site
  ajftd_fuel_details.csv      any other tables on the page (fuel description, producer, blend info ...)
Then:  python ingest_ajftd.py   (it reads data/raw/ajftd/*.csv; delete or move the old ajftd_data.csv to avoid duplicates)

Page structure assumed (taken from your scraper, which was written against the real HTML):
  <button onclick="myFunction('tbl_15')">metric_label</button>  +  <table id="tbl_15"> header row, then data rows
  headings: "Fuel Details for (02239-070416)T108", "Composition Metrics:", "Property Metrics:"
"""
from __future__ import annotations
import argparse, re, sys
from pathlib import Path

import pandas as pd
from bs4 import BeautifulSoup

OUT_DIR = Path("data/raw/ajftd")
NOT_A_VALUE = {"", "no test", "no search", "not found", "nan", "none", "n/a"}


def clean(v) -> str:
    return re.sub(r"\s+", " ", v or "").strip()


def norm_metric(v: str) -> str:
    v = re.sub(r"\([^)]*\)", " ", clean(v).lower())
    return re.sub(r"[^a-z0-9]+", "_", v).strip("_")


def table_rows(table) -> list[dict]:
    rows = table.select("tr")
    if not rows:
        return []
    headers = [clean(c.get_text(" ", strip=True)) for c in rows[0].select("th, td")]
    out = []
    for tr in rows[1:]:
        cells = [clean(c.get_text(" ", strip=True)) for c in tr.select("th, td")]
        if not any(cells):
            continue
        row = {}
        for i, v in enumerate(cells):
            k = headers[i] if i < len(headers) and headers[i] else f"column_{i + 1}"
            if k in row:
                k = f"{k}_{i + 1}"
            row[k] = v or None
        out.append(row)
    return out


def parse_page(html: str, fuel_id: str):
    soup = BeautifulSoup(html, "lxml")
    heads = [clean(h.get_text(" ", strip=True)) for h in soup.select("h1, h2, h3, h4")]
    heads = [h for h in heads if h]
    detail = next((h for h in heads if h.lower().startswith("fuel details for")), None)
    label = re.sub(r"(?i)^fuel details for\s*", "", detail).strip() if detail else None
    title = clean(soup.title.get_text(" ", strip=True)) if soup.title else None

    metrics, used_tables = [], set()
    for b in soup.select("button[onclick*='myFunction']"):
        m = re.search(r"""myFunction\(\s*['"]([^'"]+)['"]\s*\)""", b.get("onclick", ""))
        if not m:
            continue
        tid = m.group(1)
        table = soup.find(id=tid)
        if table is None:
            continue
        used_tables.add(tid)
        sec = b.find_previous(["h1", "h2", "h3", "h4"])
        section = clean(sec.get_text(" ", strip=True)).rstrip(":") if sec else ""
        metrics.append(dict(metric=norm_metric(b.get_text(" ", strip=True)), label=clean(b.get_text(" ", strip=True)),
                            table_id=tid, section=section, rows=table_rows(table)))

    other = []                                    # tables not tied to a metric button -> fuel description etc.
    for i, t in enumerate(soup.find_all("table")):
        if t.get("id") in used_tables or t.find_parent("table") is not None:
            continue
        for ri, r in enumerate(table_rows(t), 1):
            for k, v in r.items():
                if v:
                    other.append(dict(fuel_id=fuel_id, fuel_label=label, table_index=i, row=ri, column=k, value=v))
    return label, title, metrics, other


def _process_one_file(f_path):
    fid = re.sub(r"^detail_|\.html$", "", f_path.name)
    html = f_path.read_text(encoding="utf-8", errors="replace")
    label, title, metrics, other = parse_page(html, fid)
    rows = []
    for mt in metrics:
        base = dict(fuel_id=fid, fuel_name=label or "", fuel_label=label or fid, record_url="", page_title=title,
                    category=mt["section"], metric=mt["metric"], table_id=mt["table_id"])
        if not mt["rows"]:
            rows.append({**base, "measurement_index": None, "Value": None, "Unit": None})
        for i, r in enumerate(mt["rows"], 1):
            rows.append({**base, "measurement_index": i, **r})
    is_empty = not metrics
    return rows, other, is_empty, f_path.name

def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("folder", nargs="?", default="debug")
    ap.add_argument("--compare", help="scraper CSV (ajftd_data.csv) to cross-check overlapping metrics against")
    a = ap.parse_args(argv)
    files = sorted(Path(a.folder).glob("detail_*.html"))
    if not files:
        sys.exit(f"no detail_*.html files in {a.folder}/ - run from the folder that contains the scraper's debug/ directory")

    from concurrent.futures import ProcessPoolExecutor
    import os

    csv_rows, details, empty_pages = [], [], []
    workers = min(os.cpu_count() or 4, 8)
    print(f"Parsing {len(files)} detail HTML files using {workers} parallel processes...")

    with ProcessPoolExecutor(max_workers=workers) as executor:
        for rows, other, is_empty, fname in executor.map(_process_one_file, files, chunksize=25):
            csv_rows.extend(rows)
            details.extend(other)
            if is_empty:
                empty_pages.append(fname)

    df = pd.DataFrame(csv_rows)
    # same column naming as the scraper so ingest_ajftd.py needs no changes
    if "Value" in df:
        df["_num"] = pd.to_numeric(df["Value"].astype(str).str.replace(r"[<>]", "", regex=True), errors="coerce")
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    df.drop(columns=["_num"], errors="ignore").to_csv(OUT_DIR / "ajftd_all_metrics.csv", index=False, encoding="utf-8-sig")
    if details:
        pd.DataFrame(details).to_csv(OUT_DIR / "ajftd_fuel_details.csv", index=False, encoding="utf-8-sig")

    inv = (df.groupby(["category", "metric"]).agg(fuels=("fuel_id", "nunique"), rows=("fuel_id", "size"),
                                                  numeric_rows=("_num", lambda s: int(s.notna().sum())),
                                                  fuels_with_number=("fuel_id", lambda s: s[df.loc[s.index, "_num"].notna()].nunique()))
           .reset_index().sort_values(["category", "metric"]))
    inv.to_csv(OUT_DIR / "ajftd_metrics_inventory.csv", index=False)

    print(f"pages parsed: {len(files)}   pages with no metric tables: {len(empty_pages)}   metric rows: {len(df)}")
    print(f"fuel labels found: {df.fuel_label.nunique()}  e.g. {list(df.fuel_label.drop_duplicates().head(4))}")
    print("\nWHAT EXISTS ON THE SITE (metric: fuels with a real number / rows):")
    for r in inv.itertuples():
        print(f"  [{r.category[:22]:22}] {r.metric:50} fuels_with_number={r.fuels_with_number:4}  rows={r.rows}")
    if empty_pages:
        print(f"\nWARNING: no metric tables found in: {empty_pages[:5]}{' ...' if len(empty_pages) > 5 else ''} - page layout may differ; send me one.")

    if a.compare and Path(a.compare).exists():
        old = pd.read_csv(a.compare, encoding="utf-8-sig", dtype=str)
        key = ["fuel_id", "metric", "measurement_index"]
        new = df.astype({"fuel_id": str, "measurement_index": "Int64"}).astype({"measurement_index": str})
        mg = old.merge(new[key + ["Value"]], on=key, how="inner", suffixes=("_scraper", "_html"))
        bad = mg[mg["Value_scraper"].fillna("") != mg["Value_html"].fillna("")]
        print(f"\ncross-check vs {a.compare}: {len(mg)} overlapping rows, {len(bad)} differ")
        if len(bad):
            print(bad.head(5).to_string())
    print(f"\nwrote {OUT_DIR/'ajftd_all_metrics.csv'}, ajftd_metrics_inventory.csv" + (", ajftd_fuel_details.csv" if details else ""))


if __name__ == "__main__":
    main()
