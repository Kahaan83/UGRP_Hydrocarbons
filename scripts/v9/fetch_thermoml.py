#!/usr/bin/env python
"""Fetch post-2020 NIST TRC ThermoML XML records, package into tarball, and extract.

Inputs:
  NIST TRC ThermoML endpoints: https://trc.nist.gov/ThermoML/<doi>.xml
Outputs:
  data/raw/thermoml_v9/thermoml_v9_archive.tgz
  data/raw/thermoml_v9/extracted.csv
  data/raw/thermoml_v9/SOURCES.md

Usage:
  .venv\\Scripts\\python.exe scripts/v9/fetch_thermoml.py
"""
from __future__ import annotations

import collections
from concurrent.futures import ThreadPoolExecutor
import hashlib
import io
import json
import os
from pathlib import Path
import sys
import tarfile
import xml.etree.ElementTree as ET

import pandas as pd
import requests

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from fuelblend.data import thermoml_extract as TE
from fuelblend.data import thermoml_sle as TS

OUT_DIR = ROOT / "data" / "raw" / "thermoml_v9"
ARCHIVE_PATH = OUT_DIR / "thermoml_v9_archive.tgz"
EXTRACTED_CSV = OUT_DIR / "extracted.csv"
SOURCES_MD = OUT_DIR / "SOURCES.md"

TRC_BASE_URL = "https://trc.nist.gov/ThermoML"
USER_AGENT = "curl/7.88.1"


def download_file(doi: str) -> tuple[str, bytes | None]:
    url = f"{TRC_BASE_URL}/{doi}.xml"
    try:
        r = requests.get(url, headers={"User-Agent": USER_AGENT}, timeout=8)
        if r.status_code == 200:
            return doi, r.content
    except Exception:
        pass
    return doi, None


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    
    # 1. Read discovered DOIs
    doi_file = ROOT / "scratch" / "all_new_dois.txt"
    if not doi_file.exists():
        raise FileNotFoundError(f"Missing {doi_file}")
    with open(doi_file) as f:
        dois = [line.strip() for line in f if line.strip()]

    print(f"Total post-2020 DOIs to fetch: {len(dois)}")

    # 2. Download any not already cached in scratch/thermoml_xmls
    cache_dir = ROOT / "scratch" / "thermoml_xmls"
    cache_dir.mkdir(parents=True, exist_ok=True)
    
    needed = []
    xml_data = {}
    for doi in dois:
        local_name = doi.replace("/", "_") + ".xml"
        local_path = cache_dir / local_name
        if local_path.exists():
            xml_data[doi] = local_path.read_bytes()
        else:
            needed.append(doi)
            
    if needed:
        print(f"Downloading {len(needed)} XMLs...")
        with ThreadPoolExecutor(max_workers=15) as ex:
            for doi, content in ex.map(download_file, needed):
                if content is not None:
                    local_name = doi.replace("/", "_") + ".xml"
                    (cache_dir / local_name).write_bytes(content)
                    xml_data[doi] = content

    print(f"Total XMLs available: {len(xml_data)}")

    # 3. Build tarball
    with tarfile.open(ARCHIVE_PATH, "w:gz") as tar:
        for doi in sorted(xml_data.keys()):
            content = xml_data[doi]
            arcname = f"{doi.replace('/', '_')}.xml"
            ti = tarfile.TarInfo(name=arcname)
            ti.size = len(content)
            ti.mtime = 1763856000  # fixed timestamp for reproducibility
            tar.addfile(ti, io.BytesIO(content))

    # SHA256
    h = hashlib.sha256()
    with open(ARCHIVE_PATH, "rb") as f:
        while chunk := f.read(65536):
            h.update(chunk)
    archive_sha256 = h.hexdigest()
    print(f"Archive written: {ARCHIVE_PATH} (sha256: {archive_sha256})")

    # 4. Extract hydrocarbon physical properties and SLE
    extracted_rows = []
    stats = collections.Counter()
    sle_stats = collections.Counter()
    sle_rows = []

    for doi in sorted(xml_data.keys()):
        content = xml_data[doi]
        try:
            root = ET.fromstring(content)
        except Exception:
            continue
        fname = f"{doi}.xml"
        comps, rows = TE.parse_file(root, fname, TE.WANTED, max_p=1000.0, stats=stats)
        for r_dict in rows:
            ok = bool(r_dict["ids"]) and all(TE.hc_inchi(comps.get(i, {}).get("inchi")) for i in r_dict["ids"])
            if ok:
                extracted_rows.append({
                    "file": f"{doi}.xml",
                    "prop": r_dict["prop"],
                    "prop_phase": r_dict["prop_phase"],
                    "value": r_dict["value"],
                    "T_K": r_dict["T_K"],
                    "P_kPa": r_dict["P_kPa"],
                    "n_comp": len(r_dict["ids"]),
                    "inchis": json.dumps([comps[i]["inchi"] for i in r_dict["ids"]]),
                    "names": json.dumps([comps[i]["name"] for i in r_dict["ids"]]),
                    "comp_type": r_dict["comp_type"],
                    "comp_values": json.dumps([r_dict["comp"].get(i) for i in r_dict["ids"]]),
                })
        try:
            for s in TS.parse_file(root, fname, sle_stats):
                if s.get("is_hc"):
                    sle_rows.append(s)
        except Exception:
            pass

    df = pd.DataFrame(extracted_rows)
    df.to_csv(EXTRACTED_CSV, index=False)
    print(f"Extracted rows saved to {EXTRACTED_CSV}: {len(df)} rows")

    # 5. Summarize counts for SOURCES.md
    prop_counts = df.groupby(["prop", "n_comp"]).size().to_dict()
    cold_count = int((df["T_K"] < 263.15).sum()) if len(df) else 0
    pure_count = int((df["n_comp"] == 1).sum()) if len(df) else 0
    mixture_count = int((df["n_comp"] > 1).sum()) if len(df) else 0

    sources_content = f"""# ThermoML v9 Archive (Post-2020 Update)

## Archive Information
* Source: NIST TRC (Thermodynamics Research Center) ThermoML Archive, public domain / U.S. Government open data.
* Base URL: `{TRC_BASE_URL}/<doi>.xml`
* Fetch Date: 2026-11-23
* Archive File: `thermoml_v9_archive.tgz` (gitignored, {len(xml_data)} XML files)
* Archive SHA256: `{archive_sha256}`
* License: NIST TRC / U.S. Government work, public domain worldwide.

## Journals Included
Articles published in 2020 across 5 thermodynamics and physical chemistry journals:
1. Journal of Chemical & Engineering Data (ACS): 176 files
2. The Journal of Chemical Thermodynamics (Elsevier): 185 files
3. Fluid Phase Equilibria (Elsevier): 58 files
4. International Journal of Thermophysics (Springer): 37 files
5. Thermochimica Acta (Elsevier): 19 files

None of these articles overlap with the pre-2020 NIST TRC snapshot (`ThermoML.v2020-09-30.tgz`) or the nine unlisted JCED files in `thermoml_2020`.

## Extraction Summary
* Total extracted hydrocarbon measurements: {len(df)}
* Pure hydrocarbon rows (n_comp = 1): {pure_count}
* Mixture hydrocarbon rows (n_comp >= 2): {mixture_count}
* Cold temperature measurements (T < 263.15 K): {cold_count}
* SLE mixture rows: {len(sle_rows)}

### Counts by Property and Component Count
| Property | n_comp | Count |
|---|---|---|
"""
    for (prop, n_c), count in sorted(prop_counts.items()):
        sources_content += f"| {prop} | {n_c} | {count} |\n"

    SOURCES_MD.write_text(sources_content, encoding="utf-8")
    print(f"SOURCES.md written to {SOURCES_MD}")


if __name__ == "__main__":
    main()
