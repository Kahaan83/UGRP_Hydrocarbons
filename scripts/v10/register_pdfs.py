"""Register user-supplied PDFs: header check, hash, page count, text layer.

Reads only the %PDF- header, page count and the first 3 pages' text length.
"""
from __future__ import annotations

import hashlib
import sys
from pathlib import Path

import pandas as pd

EXPECTED = ["AFRL-TR-2020-0017", "AIAA 2017-0146"]
COLS = ["filename", "sha256", "is_pdf", "n_pages", "has_text_layer",
        "title_from_sources_txt", "status"]


def _sources(pdf_dir: Path) -> dict[str, str]:
    p = pdf_dir / "sources.txt"
    out: dict[str, str] = {}
    if p.exists():
        for line in p.read_text(encoding="utf-8").splitlines():
            parts = [s.strip() for s in line.split("|")]
            if len(parts) >= 2 and parts[0]:
                out[parts[0]] = parts[1]
    return out


def register(pdf_dir: Path, expected: list[str] | None = None) -> pd.DataFrame:
    pdf_dir = Path(pdf_dir)
    src = _sources(pdf_dir)
    rows = []
    for f in sorted(pdf_dir.iterdir()) if pdf_dir.exists() else []:
        if not f.is_file() or f.name == "sources.txt":
            continue
        data = f.read_bytes()
        is_pdf = data[:5] == b"%PDF-"
        n_pages, has_text = 0, False
        if is_pdf:
            try:
                import pdfplumber
                with pdfplumber.open(f) as pdf:
                    n_pages = len(pdf.pages)
                    txt = "".join((pg.extract_text() or "") for pg in pdf.pages[:3])
                    has_text = len(txt) >= 200
            except Exception:
                is_pdf = False
        rows.append(dict(filename=f.name, sha256=hashlib.sha256(data).hexdigest(),
                         is_pdf=is_pdf, n_pages=n_pages, has_text_layer=has_text,
                         title_from_sources_txt=src.get(f.name, ""),
                         status="ok" if is_pdf else "not_pdf"))
    names = " ".join(r["filename"] for r in rows).lower()
    for key in (EXPECTED if expected is None else expected):
        if key.lower().replace(" ", "") not in names.replace(" ", ""):
            rows.append(dict(filename=key, sha256="", is_pdf=False, n_pages=0,
                             has_text_layer=False, title_from_sources_txt="",
                             status="not provided"))
    return pd.DataFrame(rows, columns=COLS)


if __name__ == "__main__":
    d = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("data/raw/v10_pdfs")
    df = register(d)
    Path("results/v10").mkdir(parents=True, exist_ok=True)
    df.to_csv("results/v10/pdf_register.csv", index=False)
    print(df.drop(columns="sha256").to_string())
