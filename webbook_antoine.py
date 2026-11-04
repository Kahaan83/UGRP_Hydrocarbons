#!/usr/bin/env python
"""Scrape Antoine constants from the NIST Chemistry WebBook (phase-change page).

    python webbook_antoine.py                      # uses compounds.csv
    python webbook_antoine.py --csv my_list.csv --out data/processed/antoine.csv

WebBook form (units!):   log10( P / bar ) = A - B / ( T/K + C )
Your methods document uses kPa, so we also emit  P_kPa = 100 * 10**(A - B/(T+C)).

For every compound ALL WebBook parameter sets are saved (long file), and one
'best' set is chosen = the one whose validity range covers most of the flash-point
window 313.15-373.15 K. `covers_window` tells you whether Antoine is being used
inside its range; if False, re-fit from vapour-pressure data instead (e.g. ThermoML).
`psat_40C_kPa` is the 40 degC cross-check value called for in the methods document.
`page_name` is the compound name WebBook returned - CHECK it matches `name`.
"""
from __future__ import annotations
import argparse, io, re, time
from pathlib import Path

import pandas as pd
import requests

URL = "https://webbook.nist.gov/cgi/cbook.cgi?ID=C{cas}&Mask=4"
UA = {"User-Agent": "fuelblend-data/1.0 (academic research)"}
WINDOW = (313.15, 373.15)


def psat_kpa(A: float, B: float, C: float, T_K: float) -> float:
    return 100.0 * 10 ** (A - B / (T_K + C))


def parse_antoine(html: str) -> list[dict]:
    """Return one dict per Antoine row: T_lo, T_hi, A, B, C, reference, comment."""
    out = []
    try:
        tables = pd.read_html(io.StringIO(html))
    except ValueError:
        return out
    for t in tables:
        cols = [str(c).strip() for c in t.columns]
        if {"A", "B", "C"} <= set(cols) and any(c.startswith("Temperature") for c in cols):
            tcol = next(c for c in cols if c.startswith("Temperature"))
            t.columns = cols
            for _, r in t.iterrows():
                m = re.match(r"\s*([\d.]+)\s*(?:to|-|–)\s*([\d.]+)", str(r[tcol]))
                if not m:
                    continue
                try:
                    A, B, C = float(r["A"]), float(r["B"]), float(r["C"])
                except (TypeError, ValueError):
                    continue
                out.append(dict(T_lo=float(m.group(1)), T_hi=float(m.group(2)), A=A, B=B, C=C,
                                reference=str(r.get("Reference", "")), comment=str(r.get("Comment", ""))))
    return out


def page_name(html: str) -> str:
    """Compound name. The first <h1> on a real WebBook page is the site banner, so prefer
    <h1 id="Top">, then <title>, and never return the banner."""
    for pat in (r'<h1[^>]*id="Top"[^>]*>(.*?)</h1>', r"<title>(.*?)</title>"):
        m = re.search(pat, html, re.S | re.I)
        if m:
            t = re.sub(r"\s+", " ", re.sub(r"<.*?>", "", m.group(1))).strip()
            if t and "WebBook" not in t:
                return t
    return ""


def coverage(row: dict) -> float:
    lo, hi = max(row["T_lo"], WINDOW[0]), min(row["T_hi"], WINDOW[1])
    return max(0.0, hi - lo) / (WINDOW[1] - WINDOW[0])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", default="compounds.csv")
    ap.add_argument("--out", default="data/processed/antoine_best.csv")
    ap.add_argument("--delay", type=float, default=1.0)
    a = ap.parse_args()
    comps = pd.read_csv(a.csv)
    best, allrows = [], []
    for _, c in comps.iterrows():
        cas = str(c["cas"]).replace("-", "")
        try:
            r = requests.get(URL.format(cas=cas), headers=UA, timeout=30)
            r.raise_for_status()
        except Exception as e:  # noqa: BLE001
            print(f"[fail] {c['name']} ({c['cas']}): {e}")
            continue
        pname = page_name(r.text)
        rows = parse_antoine(r.text)
        for x in rows:
            x.update(name=c["name"], cas=c["cas"], smiles=c.get("smiles", ""), page_name=pname, coverage=coverage(x))
        allrows += rows
        if not rows:
            print(f"[none] {c['name']}: no Antoine table on page ({pname})")
        else:
            b = max(rows, key=coverage)
            b = dict(b, covers_window=b["coverage"] >= 0.999,
                     psat_40C_kPa=psat_kpa(b["A"], b["B"], b["C"], 313.15))
            best.append(b)
            print(f"[ok]   {c['name']:<24} {pname[:28]:<28} range {b['T_lo']:.0f}-{b['T_hi']:.0f} K "
                  f"covers_window={b['covers_window']}  Psat(40C)={b['psat_40C_kPa']:.4g} kPa")
        time.sleep(a.delay)
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(best).to_csv(a.out, index=False)
    pd.DataFrame(allrows).to_csv(Path(a.out).with_name("antoine_all_sets.csv"), index=False)
    print(f"\nsaved {len(best)} compounds -> {a.out}")


if __name__ == "__main__":
    main()
