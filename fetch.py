#!/usr/bin/env python
"""Download every source that can be fetched automatically.

    python fetch.py all                    # everything below
    python fetch.py pdfs thermoml          # only some groups
    python fetch.py all --skip zenodo      # everything except some groups
    python fetch.py zenodo --list          # just list Zenodo files + sizes, download nothing
    python fetch.py zenodo --include "flash|visc" --max-gb 2   # filter by filename regex / size cap

Ctrl-C skips the CURRENT group and continues with the next (press Ctrl-C twice
within 3 s to stop everything). The summary is always printed at the end.

Groups: git | zenodo | pdfs | thermoml
Files land in data/raw/<group>/ ; data/raw/MANIFEST.json records url, size,
sha256 and time for every file so the dataset version is reproducible.
Anything that fails is listed at the end with what to do by hand.
Run it on a machine with open internet access (university network is fine).
"""
from __future__ import annotations
import argparse, hashlib, json, os, re, subprocess, sys, time
from datetime import datetime, timezone
from pathlib import Path

import requests

RAW = Path("data/raw")
MANIFEST = RAW / "MANIFEST.json"
UA = {"User-Agent": "fuelblend-data/1.0 (academic research)"}
BROWSER = {"User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/124.0 Safari/537.36", "Accept": "application/pdf,*/*;q=0.8", "Accept-Language": "en-US,en;q=0.9"}
ZENODO_API = os.environ.get("ZENODO_API", "https://zenodo.org/api/records")
CONNECT_TO, READ_TO = 15, 30          # seconds; READ_TO applies per chunk
FILE_DEADLINE = 1800                  # hard cap per file (s) so a slow trickle can't hang forever
BACKOFF_MAX = 120

# ----------------------------------------------------------------- sources
GIT = {
    # contains nist-logV (239k viscosity rows) and MON datasets
    "chemixhub": "https://github.com/chemcognition-lab/chemixhub.git",
    # benchmark code (flash point / viscosity / DCN / MON) - may need VPN/uni network
    "mixture_evaluation": "https://gitlab.kuleuven.be/creas/vermeiregroup/mixture_evaluation.git",
}

ZENODO = {
    # Leenhouts et al. curated pure+mixture datasets (flash point, viscosity, DCN, MON)
    "leenhouts_benchmark": "19914760",
    # Bilodeau et al.: pretrained_models.zip (95 MB) CONTAINS the raw NIST viscosity table
    # nist_dippr_source/NIST_Visc_Data.csv (239,200 rows; Visc in Pa.s) plus bp/mp data. Keep it.
    "bilodeau_nist_viscosity": "8042966",
}

PDFS = {
    "NISTIR_6646_RP1.pdf": "https://nvlpubs.nist.gov/nistpubs/Legacy/IR/nistir6646.pdf",
    "CRC_647_world_fuel_sampling.pdf": "https://crcao.org/wp-content/uploads/2025/07/CRC-647.pdf",
    "CRC_AV-27-18_fuel_properties.pdf": "https://crcao.org/wp-content/uploads/2022/03/CRC-AV-27-18_Final-Report_UDR-TR-2021-90.pdf",
    "CRC_AV-11-09_freezing_point.pdf": "https://crcao.org/wp-content/uploads/2019/05/FINAL-REPORT-Freezing-point-CRC-project-AV-11-09-v3a.pdf",
    "AFRL_TR-2020-0017_jet_fuel_properties.pdf": "https://apps.dtic.mil/sti/pdfs/AD1093317.pdf",
    "CRC_handbook_aviation_fuel_properties.pdf": "https://apps.dtic.mil/sti/pdfs/ADA132106.pdf",
    "Zabetakis_Bulletin_627_flammability.pdf": "https://digital.library.unt.edu/ark:/67531/metadc1444558/m2/1/high_res_d/7328370.pdf",
    "Li_2020_flash_point_GNN_arXiv2002.11315.pdf": "https://arxiv.org/pdf/2002.11315",
    "Leenhouts_2026_mixture_eval_arXiv2605.29698.pdf": "https://arxiv.org/pdf/2605.29698",
    "Chew_2024_viscosity_SI.pdf": "https://media.springernature.com/original/springer-static/esm/art%3A10.1186%2Fs13321-024-00820-5/MediaObjects/13321_2024_820_MOESM1_ESM.pdf",
}

# ThermoML archive (NIST TRC, doi 10.18434/mds2-2422). File name taken from the
# catalogue entry; the direct URL pattern is a best guess -> manual fallback below.
THERMOML_URLS = [
    # (url, size_bytes, sha256) - taken from the NIST record metadata (data.nist.gov/od/id/mds2-2422)
    ("https://data.nist.gov/od/ds/mds2-2422/ThermoML.v2020-09-30.tgz", 189433115,
     "231161b5e443dc1ae0e5da8429d86a88474cb722016e5b790817bb31c58d7ec2"),
]

MANUAL = [
    "ThermoML (if auto-download failed): open https://data.nist.gov/od/id/mds2-2422 , "
    "download ThermoML.v2020-09-30.tgz (189 MB) into data/raw/thermoml/ then run: python thermoml_extract.py",
    "DTIC PDFs (if 403 persists): open the two apps.dtic.mil links from the FAILED list in a normal browser "
    "(may ask you to accept a notice) and save them in data/raw/pdfs/ with the exact filenames shown",
    "AJFTD jet-fuel database: register at https://altjetfuels.illinois.edu (admin approval needed), "
    "then export fuel records to data/raw/ajftd/",
    "NIST RP-1 blend data (11 blends): Energy & Fuels 2012, doi 10.1021/ef3004009 (publisher/SI download)",
    "Pour/cloud point literature (Lira thesis d.lib.msu.edu/etd/3010, Daridon 2002, Coutinho, Joshi/IIP): extract tables by hand",
    "Li et al. flash-point dataset (10,575 molecules): supporting information of arXiv 2002.11315 / J. Chem. Inf. Model.",
    "Chew et al. viscosity dataset: check the 'Availability of data' section of doi 10.1186/s13321-024-00820-5",
    "CRC AV-33-22 World Jet Fuel Survey 2023-24 (published by CRC/Volpe, check crcao.org)",
    "DIPPR is proprietary - only use if your library licenses it",
]

# ----------------------------------------------------------------- helpers
failed: list[str] = []


def sha256(p: Path) -> str:
    h = hashlib.sha256()
    with p.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def load_manifest() -> dict:
    return json.loads(MANIFEST.read_text()) if MANIFEST.exists() else {}


def record(path: Path, url: str):
    m = load_manifest()
    m[str(path)] = {
        "url": url,
        "bytes": path.stat().st_size,
        "sha256": sha256(path),
        "fetched_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
    }
    RAW.mkdir(parents=True, exist_ok=True)
    MANIFEST.write_text(json.dumps(m, indent=2))


def valid_pdf(p: Path) -> bool:
    try:
        with p.open("rb") as f:
            return f.read(5) == b"%PDF-"
    except OSError:
        return False


def _retry_after(r) -> float | None:
    try:
        return float(r.headers.get("Retry-After"))
    except (TypeError, ValueError):
        return None


def download(url: str, dest: Path, retries: int = 3, expect_size: int | None = None,
             expect_sha256: str | None = None) -> bool:
    """Resumable, rate-limit aware, bounded. Returns True on success."""
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.suffix.lower() == ".pdf" and dest.exists() and not valid_pdf(dest):
        print(f"  [bad]  {dest.name} is not a real PDF (earlier run saved an HTML page) - deleting and refetching")
        dest.unlink()
    if dest.exists() and dest.stat().st_size > 0 and (expect_size is None or dest.stat().st_size == expect_size):
        print(f"  [skip] {dest.name} already present")
        return True
    tmp = dest.with_suffix(dest.suffix + ".part")
    last_err = ""
    use_browser = False
    for attempt in range(1, retries + 1):
        try:
            have = tmp.stat().st_size if tmp.exists() else 0
            hdr = dict(BROWSER if use_browser else UA, **({"Range": f"bytes={have}-"} if have else {}))
            t0 = time.time()
            with requests.get(url, headers=hdr, stream=True, timeout=(CONNECT_TO, READ_TO),
                              allow_redirects=True) as r:
                if r.status_code == 429 or r.status_code >= 500:
                    wait = min(_retry_after(r) or 10 * 2 ** attempt, BACKOFF_MAX)
                    last_err = f"HTTP {r.status_code}"
                    print(f"  [try {attempt}/{retries}] {dest.name}: {last_err}, waiting {wait:.0f}s")
                    time.sleep(wait)
                    continue
                if r.status_code == 403 and not use_browser:
                    use_browser = True
                    last_err = "HTTP 403 (retrying with browser-like headers)"
                    print(f"  [try {attempt}/{retries}] {dest.name}: {last_err}")
                    continue
                if r.status_code == 416:               # range not satisfiable -> .part already complete/bad
                    tmp.unlink(missing_ok=True)
                    last_err = "HTTP 416 (stale partial file removed)"
                    continue
                r.raise_for_status()
                resumed = r.status_code == 206
                total = int(r.headers.get("Content-Length", 0)) + (have if resumed else 0)
                done = have if resumed else 0
                shown = 0
                with tmp.open("ab" if resumed else "wb") as f:
                    for chunk in r.iter_content(1 << 20):
                        f.write(chunk)
                        done += len(chunk)
                        if time.time() - t0 > FILE_DEADLINE:
                            raise TimeoutError(f"exceeded {FILE_DEADLINE}s for this file")
                        if done - shown >= 25 << 20:    # progress every ~25 MB
                            shown = done
                            pct = f" ({100*done/total:.0f}%)" if total else ""
                            print(f"         {dest.name}: {done/1e6:.0f} MB{pct}", flush=True)
                if dest.suffix.lower() == ".pdf" and not valid_pdf(tmp):
                    tmp.unlink(missing_ok=True)
                    last_err = ("server returned a web page (consent/login/bot check), not a PDF - "
                                "open the URL in a normal browser and save it into data/raw/pdfs/ with this file name")
                    print(f"  [html] {dest.name}: {last_err}")
                    break
                tmp.rename(dest)
            if expect_sha256 and sha256(dest) != expect_sha256:
                dest.unlink()
                raise ValueError("sha256 mismatch - file deleted, will retry")
            record(dest, url)
            print(f"  [ok]   {dest.name} ({dest.stat().st_size/1e6:.1f} MB)")
            return True
        except KeyboardInterrupt:
            raise
        except Exception as e:  # noqa: BLE001
            last_err = f"{type(e).__name__}: {str(e)[:120]}"
            print(f"  [try {attempt}/{retries}] {dest.name}: {last_err}")
            time.sleep(min(5 * attempt, BACKOFF_MAX))
    failed.append(f"{dest.name}  <-  {url}   ({last_err})")
    return False


# ----------------------------------------------------------------- groups
def fetch_git():
    print("== git repositories")
    for name, url in GIT.items():
        dest = RAW / "git" / name
        if (dest / ".git").exists():
            subprocess.run(["git", "-C", str(dest), "pull", "--ff-only"], check=False)
            continue
        dest.parent.mkdir(parents=True, exist_ok=True)
        r = subprocess.run(["git", "clone", "--depth", "1", url, str(dest)])
        if r.returncode != 0:
            failed.append(f"git clone {url}")


def zenodo_files(record_json: dict) -> list[tuple[str, str, int | None]]:
    """Return (filename, url, size) from a Zenodo API record (old and new schema)."""
    out = []
    for f in record_json.get("files", []):
        name = f.get("key") or f.get("filename")
        links = f.get("links", {})
        url = links.get("self") or links.get("download")
        if not url:
            url = f"https://zenodo.org/records/{record_json.get('id')}/files/{name}?download=1"
        out.append((name, url, f.get("size")))
    return out


def fetch_zenodo(list_only=False, include=None, max_gb=10.0):
    print("== zenodo")
    pat = re.compile(include, re.I) if include else None
    for label, rec in ZENODO.items():
        js = None
        for attempt in range(1, 4):
            try:
                r = requests.get(f"{ZENODO_API}/{rec}", headers=UA, timeout=(CONNECT_TO, READ_TO))
                if r.status_code == 429:
                    time.sleep(min(_retry_after(r) or 15 * attempt, BACKOFF_MAX))
                    continue
                r.raise_for_status()
                js = r.json()
                break
            except KeyboardInterrupt:
                raise
            except Exception as e:  # noqa: BLE001
                print(f"  record {rec} attempt {attempt}/3: {type(e).__name__}: {str(e)[:100]}")
                time.sleep(5 * attempt)
        if js is None:
            failed.append(f"zenodo record {rec} ({label}) - browse https://zenodo.org/records/{rec} and download manually")
            continue
        files = zenodo_files(js)
        if pat:
            files = [x for x in files if pat.search(x[0])]
        tot = sum(x[2] or 0 for x in files)
        print(f"  record {rec} ({label}): {len(files)} files, {tot/1e9:.2f} GB")
        for name, url, size in files[:50]:
            print(f"     {name}  {(size or 0)/1e6:.1f} MB")
        if len(files) > 50:
            print(f"     ... and {len(files)-50} more (use --include REGEX to select)")
        if list_only:
            continue
        if tot / 1e9 > max_gb:
            failed.append(f"zenodo {rec}: {tot/1e9:.1f} GB exceeds --max-gb {max_gb}; rerun with --include or a higher --max-gb")
            continue
        for name, url, size in files:
            download(url, RAW / "zenodo" / label / name, expect_size=size)
            time.sleep(1.0)              # be polite; Zenodo rate-limits anonymous clients


def fetch_pdfs():
    print("== pdfs / reports")
    for name, url in PDFS.items():
        download(url, RAW / "pdfs" / name)


def fetch_thermoml():
    print("== thermoml archive")
    for url, size, digest in THERMOML_URLS:
        if download(url, RAW / "thermoml" / url.rsplit("/", 1)[-1], expect_size=size, expect_sha256=digest):
            return


GROUPS = {"git": fetch_git, "zenodo": fetch_zenodo, "pdfs": fetch_pdfs, "thermoml": fetch_thermoml}


def summary():
    print("\n=== FAILED / SKIPPED (retry later or download manually) ===")
    print("\n".join(failed) or "none")
    RAW.mkdir(parents=True, exist_ok=True)
    if failed:
        (RAW / "fetch_failures.txt").write_text("\n".join(failed))
    else:
        (RAW / "fetch_failures.txt").unlink(missing_ok=True)   # don't leave a stale failure list
    print("\n=== ALWAYS MANUAL / GATED ===")
    print("\n".join(f"- {m}" for m in MANUAL))


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("groups", nargs="*", default=["all"], choices=["all", *GROUPS])
    ap.add_argument("--skip", nargs="*", default=[], choices=list(GROUPS))
    ap.add_argument("--list", action="store_true", help="zenodo: list files only")
    ap.add_argument("--include", help="zenodo: regex on filenames to download")
    ap.add_argument("--max-gb", type=float, default=10.0, help="zenodo: refuse records bigger than this")
    a = ap.parse_args(argv)
    todo = [g for g in (GROUPS if "all" in a.groups else a.groups) if g not in a.skip]
    try:
        for g in todo:
            try:
                if g == "zenodo":
                    fetch_zenodo(a.list, a.include, a.max_gb)
                else:
                    GROUPS[g]()
            except KeyboardInterrupt:
                failed.append(f"group '{g}' interrupted by Ctrl-C (rerun: python fetch.py {g})")
                print(f"\n  ^C -> skipping rest of '{g}'. Press Ctrl-C again within 3 s to stop everything.")
                time.sleep(3)
    except KeyboardInterrupt:
        print("\n  stopped by user")
    summary()


if __name__ == "__main__":
    main()
