"""Parse CRC Report 647 Appendix G, Table G.1 (hydrocarbon species types, PDF pages 125-127) for the 57 fuels.

Writes data/raw/pour_point/crc647_hydrocarbon_types_g1.csv in the same long contract as crc647_fuel_properties.csv:
    fuel_code,fuel_name,property,value_raw,unit_raw,value_std,unit_std,T_C,method,source_pdf,page

* Page 125: ASTM D6379 mono/di/total aromatics and total saturates (vol%), GC/MS benzene and naphthalene (ug/mL).
* Pages 126-127: ASTM D2425 hydrocarbon types (12 columns).  The table prints no unit for these; the saturate types
  sum to the D6379 saturates and the aromatic types to the D6379 aromatics (checked in ``consistency``), so they are
  D2425 type shares scaled to the two HPLC fractions.  unit_raw is therefore "% (unit not printed)".
* Cells below the detection limit ("<0.2", "<0.1", "<0.5") keep the qualifier in value_raw; value_std holds the
  printed limit (same convention as the other CRC tables).  Consumers decide how to treat them.
* Nothing is fabricated: a row without exactly 6 numeric cells stops the script.
"""
from __future__ import annotations

import csv
import re
import sys
from pathlib import Path

import pdfplumber

ROOT = Path(__file__).resolve().parents[2]
PDF_NAME = "CRC_647_world_fuel_sampling.pdf"
PDF_CANDIDATES = [ROOT / "data/raw/pdfs" / PDF_NAME, ROOT / "extras/fuelblend_data/data/raw/pdfs" / PDF_NAME]
OUT = ROOT / "data/raw/pour_point/crc647_hydrocarbon_types_g1.csv"
COLUMNS = ["fuel_code", "fuel_name", "property", "value_raw", "unit_raw", "value_std", "unit_std", "T_C", "method",
           "source_pdf", "page"]
N_FUELS = 57
HW = 4.0            # half-width of a column bucket (px); columns are aligned to < 1 px
GAP = 1.0           # max gap (px) between fragments of one split cell (observed: <= 0.1)
U2425 = "% (unit not printed)"

# page -> (x-centre of the sample-number column, [(x-centre, property, unit_raw, unit_std, method)])
PAGES = {
    125: (228.6, [(265.1, "monoaromatics_vol_pct", "vol%", "vol%", "ASTM D6379"),
                  (312.3, "diaromatics_vol_pct", "vol%", "vol%", "ASTM D6379"),
                  (359.6, "aromatics_vol_pct", "vol%", "vol%", "ASTM D6379"),
                  (406.8, "saturates_vol_pct", "vol%", "vol%", "ASTM D6379"),
                  (454.0, "benzene_ug_ml", "ug/mL", "ug/mL", "GC/MS"),
                  (501.2, "naphthalene_ug_ml", "ug/mL", "ug/mL", "GC/MS")]),
    126: (206.9, [(244.7, "hc_type_paraffins_pct", U2425, "%", "ASTM D2425"),
                  (294.6, "hc_type_cycloparaffins_pct", U2425, "%", "ASTM D2425"),
                  (344.4, "hc_type_dicycloparaffins_pct", U2425, "%", "ASTM D2425"),
                  (394.2, "hc_type_tricycloparaffins_pct", U2425, "%", "ASTM D2425"),
                  (443.9, "hc_type_alkylbenzenes_pct", U2425, "%", "ASTM D2425"),
                  (493.8, "hc_type_indans_tetralins_pct", U2425, "%", "ASTM D2425")]),
    127: (200.9, [(239.5, "hc_type_indenes_pct", U2425, "%", "ASTM D2425"),
                  (290.4, "hc_type_naphthalene_pct", U2425, "%", "ASTM D2425"),
                  (341.4, "hc_type_naphthalenes_pct", U2425, "%", "ASTM D2425"),
                  (392.3, "hc_type_acenaphthenes_pct", U2425, "%", "ASTM D2425"),
                  (443.2, "hc_type_acenaphthylenes_pct", U2425, "%", "ASTM D2425"),
                  (494.1, "hc_type_tricyclic_aromatics_pct", U2425, "%", "ASTM D2425")]),
}
SAT_TYPES = ("paraffins", "cycloparaffins", "dicycloparaffins", "tricycloparaffins")
ARO_TYPES = ("alkylbenzenes", "indans_tetralins", "indenes", "naphthalene", "naphthalenes", "acenaphthenes",
             "acenaphthylenes", "tricyclic_aromatics")
CELL = re.compile(r"^(<?)(\d+(?:\.\d+)?)$")


def _merge_fragments(toks, gap=GAP):
    """pdfplumber sometimes splits one printed cell into touching tokens ('<','0.1' or '2','.70'); join tokens of one
    line (sorted by x0) whose horizontal gap is < ``gap`` px.  Columns are ~47 px apart, so nothing joins across them."""
    out = []
    for t, x0, x1 in toks:
        if out and x0 - out[-1][2] < gap:
            out[-1] = (out[-1][0] + t, out[-1][1], x1)
        else:
            out.append((t, x0, x1))
    return out


def parse_page(page, pno):
    code_x, cols = PAGES[pno]
    words = [w for w in page.extract_words(keep_blank_chars=False, extra_attrs=["upright"]) if w["upright"]]
    codes = [w for w in words if re.fullmatch(r"\d{3}", w["text"]) and abs((w["x0"] + w["x1"]) / 2 - code_x) <= HW]
    rows = []
    for cw in sorted(codes, key=lambda w: w["top"]):
        line = sorted((w for w in words if abs(w["top"] - cw["top"]) <= 4.0), key=lambda w: w["x0"])
        name = " ".join(w["text"] for w in line if w["x1"] < cw["x0"])
        cells = _merge_fragments([(w["text"], w["x0"], w["x1"]) for w in line if w["x0"] > cw["x1"]])
        got = {}
        for t, x0, x1 in cells:
            xc = (x0 + x1) / 2
            j = min(range(len(cols)), key=lambda k: abs(cols[k][0] - xc))
            if abs(cols[j][0] - xc) > HW + (3.0 if t.startswith("<") else 0.0) or not CELL.match(t) or j in got:
                raise ValueError(f"p{pno} fuel {cw['text']}: cannot place cell {t!r} at x={xc:.1f}")
            got[j] = t
        if len(got) != len(cols):
            raise ValueError(f"p{pno} fuel {cw['text']}: {len(got)} of {len(cols)} cells")
        for j, (_x, prop, u_raw, u_std, method) in enumerate(cols):
            rows.append(dict(fuel_code=cw["text"], fuel_name=name, property=prop, value_raw=got[j], unit_raw=u_raw,
                             value_std=f"{float(CELL.match(got[j]).group(2)):.10g}", unit_std=u_std, T_C="",
                             method=method, source_pdf=PDF_NAME, page=pno))
    n = len({r["fuel_code"] for r in rows})
    if n != N_FUELS:
        raise ValueError(f"p{pno}: {n} fuels parsed, expected {N_FUELS}")
    return rows


def wide(rows, censored_as=0.0):
    """{fuel_code: {property: value}}; cells printed with '<' take ``censored_as``."""
    out = {}
    for r in rows:
        v = censored_as if str(r["value_raw"]).startswith("<") else float(r["value_std"])
        out.setdefault(r["fuel_code"], {})[r["property"]] = v
    return out


def consistency(rows):
    """Per-fuel closure checks (censored cells counted as 0).  Returns a list of (fuel, check, value)."""
    out = []
    for code, d in sorted(wide(rows).items()):
        sat = sum(d[f"hc_type_{k}_pct"] for k in SAT_TYPES)
        aro = sum(d[f"hc_type_{k}_pct"] for k in ARO_TYPES)
        out.append((code, "mono+di - total aromatics (D6379)", d["monoaromatics_vol_pct"] + d["diaromatics_vol_pct"] - d["aromatics_vol_pct"]))
        out.append((code, "aromatics + saturates - 100 (D6379)", d["aromatics_vol_pct"] + d["saturates_vol_pct"] - 100.0))
        out.append((code, "D2425 saturate types - D6379 saturates", sat - d["saturates_vol_pct"]))
        out.append((code, "D2425 aromatic types - D6379 aromatics", aro - d["aromatics_vol_pct"]))
        out.append((code, "D2425 all types - 100", sat + aro - 100.0))
    return out


def main():
    pdf_path = next((p for p in PDF_CANDIDATES if p.exists()), None)
    if pdf_path is None:
        sys.exit("CRC 647 PDF not found in: " + ", ".join(str(p) for p in PDF_CANDIDATES))
    rows = []
    with pdfplumber.open(pdf_path) as pdf:
        for pno in sorted(PAGES):
            rows += parse_page(pdf.pages[pno - 1], pno)
    codes = {pno: sorted({r["fuel_code"] for r in rows if r["page"] == pno}) for pno in PAGES}
    assert codes[125] == codes[126] == codes[127], "fuel codes differ between the three pages"
    OUT.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=COLUMNS)
        w.writeheader()
        w.writerows(rows)
    n_cens = sum(str(r["value_raw"]).startswith("<") for r in rows)
    print(f"wrote {OUT.relative_to(ROOT).as_posix()}: {len(rows)} rows, {len(codes[125])} fuels, {n_cens} censored cells")
    chk = consistency(rows)
    for name in dict.fromkeys(c[1] for c in chk):
        v = [c[2] for c in chk if c[1] == name]
        worst = max((c for c in chk if c[1] == name), key=lambda c: abs(c[2]))
        print(f"  {name:42s} mean {sum(v) / len(v):+.3f}  min {min(v):+.3f}  max {max(v):+.3f}  (worst fuel {worst[0]})")


if __name__ == "__main__":
    main()
