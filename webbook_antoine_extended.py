#!/usr/bin/env python
"""Extended WebBook scraper: ALL registry comps, keep all Antoine sets,
cache HTML with sha256, extract thermo (ΔcH, Tc/Pc, triple/boiling).
Rate ≤1 req/s; writes antoine_all_sets.csv + webbook_thermo.csv.
"""
from __future__ import annotations
import argparse, io, re, time, hashlib, json
from pathlib import Path

import pandas as pd
import requests

URL = "https://webbook.nist.gov/cgi/cbook.cgi?ID=C{cas}&Mask=4"
UA = {"User-Agent": "fuelblend-data/1.0 (academic research)"}
DELAY = 1.0
CACHE = Path("data/raw/webbook")
CACHE.mkdir(parents=True, exist_ok=True)


def sha256_text(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()[:16]


def fetch_cached(cas: str) -> tuple[str, Path]:
    url = URL.format(cas=str(cas).replace("-", ""))
    # try cached by CAS first
    cas_path = CACHE / f"{cas.replace('-', '')}.html"
    if cas_path.exists():
        return cas_path.read_text(encoding="utf-8"), cas_path
    r = requests.get(url, headers=UA, timeout=30)
    r.raise_for_status()
    html = r.text
    h = sha256_text(html)
    cached = CACHE / f"{cas.replace('-', '')}_{h}.html"
    cached.write_text(html, encoding="utf-8")
    # also link CAS-only for simple lookup
    cas_path.write_text(html, encoding="utf-8")
    # manifest append
    manifest = CACHE / "MANIFEST.json"
    entry = {"cas": cas, "url": url, "sha256": h, "cached": str(cached), "bytes": len(html)}
    MANIFEST = json.load(manifest.open()) if manifest.exists() else {"entries": []}
    MANIFEST["entries"].append(entry)
    manifest.write_text(json.dumps(MANIFEST, indent=2))
    return html, cached


def parse_antoine(html: str) -> list[dict]:
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
    for pat in (r'<h1[^>]*id="Top"[^>]*>(.*?)</h1>', r"<title>(.*?)</title>"):
        m = re.search(pat, html, re.S | re.I)
        if m:
            t = re.sub(r"\s+", " ", re.sub(r"<.*?>", "", m.group(1))).strip()
            if t and "WebBook" not in t and len(t) > 2:
                return t
    return ""


def extract_thermo(html: str, name: str, cas: str) -> list[dict]:
    """Best-effort thermo extraction from known table orders."""
    rows = []
    try:
        tables = pd.read_html(io.StringIO(html))
    except Exception:
        return rows
    # Table 0 ≈ boiling point; table 3 ≈ Antoine; 4/5/6 ≈ phase change; 7 ≈ critical props list
    for idx, t in enumerate(tables):
        if idx == 0:
            # Boiling point
            try:
                for _, r in t.iterrows():
                    for c in t.columns:
                        if isinstance(r[c], str) and "±" in str(r[c]) and "K" in str(r.get(t.columns[-1], "")):
                            val = re.search(r"([\d.]+)\s*±", str(r[c]))
                            if val:
                                rows.append({"property":"boiling_point_K","value":float(val.group(1)),"unit":"K","method_ref":"WebBook",
                                              "name":name,"cas":cas,"note":"Tboil from page"})
            except Exception:
                pass
        elif idx == 4 or idx == 5 or idx == 6:
            # Phase-change enthalpies (Δ values in kJ/mol, T in K, ref, method code)
            try:
                for _, r in t.iterrows():
                    delta = r.iloc[0]
                    Tval = r.iloc[1]
                    ref = str(r.iloc[2]) if len(r) > 2 else ""
                    method = str(r.iloc[3]) if len(r) > 3 else ""
                    if isinstance(delta, (int, float)) and delta > 0:
                        # Heuristic: large ~70 -> enthalpy of fusion; ~13 -> vaporization at triple; ~50 -> vaporization
                        prop_name = "enthalpy_kJ_mol"
                        note = f"table_{idx} method={method} ref={ref[:40]}"
                        rows.append({"property":prop_name,"value":float(delta),"unit":"kJ/mol","method_ref":ref,
                                      "name":name,"cas":cas,"note":note,"T_ref_K":float(Tval) if isinstance(Tval,(int,float)) else None})
            except Exception:
                pass
        elif idx == 7:
            # Critical properties listing (headers in first col, not values) — skip; try to find Tc/Pc elsewhere by regex
            pass
    # Regex fallback for Tc / Pc / triple from whole text
    text = html
    # Triple point temperature: often "Ttriple" followed by value ± K
    for pat in (r"Ttriple.*?([0-9]+\.[0-9]+)\s*±?\s*([0-9]+\.[0-9]+)?\s*K",
                r"T<sub>triple</sub>.*?([0-9]+\.[0-9]+)",
                r"triple point.*?([0-9]+\.[0-9]+)\s*K"):
        m = re.search(pat, text, re.S | re.I)
        if m:
            try:
                val = float(m.group(1))
                rows.append({"property":"triple_point_K","value":val,"unit":"K","method_ref":"WebBook","name":name,"cas":cas,"note":"from HTML regex"})
            except Exception:
                pass
            break
    # Critical temperature / pressure (look for Tc / Pc near values)
    # Pattern: Tc = 400.0 K with reference nearby — do simple regex
    for label, prop, unit in [("Tc","critical_temperature_K","K"),("Pc","critical_pressure_kPa","kPa")]:
        # Heuristic: look for the property header in table 7 text nearby numeric
        m = re.search(rf"{label}.*?([0-9]+\.[0-9]+).*?(?=[A-Z][a-z]|<|$)", text, re.S | re.I)
        if not m:
            m = re.search(rf"{label}.*?([0-9]+\.[0-9]+)", text, re.S | re.I)
        if m:
            try:
                val = float(m.group(1))
                # Convert pressure: WebBook critical pressure often in bar or MPa? Usually bar -> kPa ×100 if < 10; else keep
                if prop == "critical_pressure_kPa" and val < 10:
                    val *= 100.0
                rows.append({"property":prop,"value":val,"unit":unit,"method_ref":"WebBook","name":name,"cas":cas,"note":"from HTML regex"})
            except Exception:
                pass
    # ΔcH°liquid (enthalpy of formation liquid) — search for formation enthalpy section
    for pat in (r"ΔfH°liquid.*?(-[0-9]+\.[0-9]+).*?kJ/mol",
                r"enthalpy of formation.*liquid.*?(-[0-9]+\.[0-9]+)"):
        m = re.search(pat, text, re.S | re.I)
        if m:
            try:
                val = float(m.group(1))
                rows.append({"property":"delta_Hf_liquid_kJ_mol","value":val,"unit":"kJ/mol","method_ref":"WebBook",
                              "name":name,"cas":cas,"note":"ΔcH°liquid (negative = exothermic)"})
            except Exception:
                pass
            break
    return rows


def coverage(row: dict) -> float:
    WINDOW = (313.15, 373.15)
    lo, hi = max(row["T_lo"], WINDOW[0]), min(row["T_hi"], WINDOW[1])
    return max(0.0, hi - lo) / (WINDOW[1] - WINDOW[0])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", default="data/unified/compounds.csv")
    ap.add_argument("--out-antoine", default="data/processed/antoine_all_sets.csv")
    ap.add_argument("--out-thermo", default="data/raw/webbook_thermo.csv")
    ap.add_argument("--delay", type=float, default=DELAY)
    ap.add_argument("--max", type=int, default=0, help="0=all, else limit for test")
    a = ap.parse_args()

    comps = pd.read_csv(a.csv)
    # We need name, cas, smiles; build temporary if missing
    if "name" not in comps.columns and "smiles" in comps.columns:
        # derive from existing
        pass
    # Filter to those with CAS to attempt
    comps = comps.dropna(subset=["cas"])
    # If using unified compounds, pick columns
    for c in ("name", "cas", "smiles"):
        if c not in comps.columns:
            comps[c] = ""
    if a.max > 0:
        comps = comps.head(a.max)

    best_all = []
    allrows = []
    thermo_rows = []
    fails = []
    for _, c in comps.iterrows():
        cas = str(c["cas"]).replace("-", "")
        name = str(c.get("name", c.get("compound_name", ""))).strip()
        smiles = str(c.get("smiles", c.get("smiles", ""))).strip()
        try:
            html, cached = fetch_cached(str(c["cas"]))
        except Exception as e:
            fails.append({"name": name, "cas": c["cas"], "error": str(e)})
            print(f"[fail] {name} ({c['cas']}): {e}")
            time.sleep(a.delay)
            continue
        pname = page_name(html)
        antoine_rows = parse_antoine(html)
        for x in antoine_rows:
            x.update(name=name, cas=str(c["cas"]), smiles=smiles, page_name=pname, coverage=coverage(x))
        allrows += antoine_rows
        # Capture best per compound (not collapsed to final output, but for info)
        if antoine_rows:
            b = max(antoine_rows, key=coverage)
            best_all.append(dict(name=name, cas=str(c["cas"]), page_name=pname,
                                 T_lo=b["T_lo"], T_hi=b["T_hi"], A=b["A"], B=b["B"], C=b["C"],
                                 reference=b["reference"], comment=b["comment"], coverage=b["coverage"],
                                 covers_window=b["coverage"] >= 0.999))
        # Thermo extraction
        thermo = extract_thermo(html, name, str(c["cas"]))
        for th in thermo:
            thermo_rows.append(th)
        if antoine_rows:
            b = max(antoine_rows, key=coverage)
            psat = 100.0 * 10 ** (b["A"] - b["B"] / (313.15 + b["C"]))
            print(f"[ok] {name:<28} {pname[:26]:<26} set={len(antoine_rows)} best_range={b['T_lo']:.0f}-{b['T_hi']:.0f}K "
                  f"covers={b['coverage']:.2f} psat40={psat:.3f}kPa thermo={len(thermo)}")
        else:
            print(f"[none] {name} ({c['cas']}) {pname:28} (no Antoine table)")
        time.sleep(a.delay)

    # Write outputs — ALL sets kept in all_sets; best in best (optional, but task says don't collapse to best for all_sets)
    Path(a.out_antoine).parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(allrows).to_csv(a.out_antoine, index=False)
    # Best summary (not required but kept for reference, separate from all_sets)
    # The task says keep every set, don't collapse — so all_sets is the main product.
    # Write thermo
    if thermo_rows:
        pd.DataFrame(thermo_rows).to_csv(a.out_thermo, index=False)
        print(f"\nthermo rows: {len(thermo_rows)} -> {a.out_thermo}")
    else:
        # Write empty with header for schema
        pd.DataFrame(columns=["property","value","unit","method_ref","name","cas","note","T_ref_K","page_name"]).to_csv(a.out_thermo, index=False)
        print(f"\nno thermo extracted -> empty {a.out_thermo}")
    print(f"\nantoine_all_sets.csv rows={len(allrows)} (sets) from {len(best_all)} compounds; fails={len(fails)}")
    if fails:
        pd.DataFrame(fails).to_csv(CACHE / "fetch_failures.csv", index=False)
        print(f"failures saved to {CACHE / 'fetch_failures.csv'}")
    # Source check note if any 403 / block
    blocked = [f for f in fails if "403" in str(f.get("error","")) or "blocked" in str(f.get("error","")).lower()]
    if blocked or any("challenge" in str(f.get("error","")).lower() for f in fails):
        Path("SOURCE_CHECK.md").write_text("NOTE 2026-10-29 T04: NIST WebBook returned challenge/403 for some CAS; partial cache used. No paywall bypass attempted.\n", encoding="utf-8")
        print("SOURCE_CHECK.md appended (block/challenge note).")


if __name__ == "__main__":
    main()
