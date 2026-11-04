"""
fetch_pubchem.py - T08: PubChem experimental properties fetch (PUG-View API)

Queries PubChem PUG-View API for experimental properties of pure hydrocarbons:
- Flash Point (converted to °C)
- Vapor Pressure (converted to kPa)
- Heat of Combustion (converted to MJ/kg)
- Boiling Point (converted to °C)

Output:
- data/raw/pubchem/pubchem_experimental.csv (Tier B)
- data/raw/pubchem/pubchem_experimental_rejects.csv
- Cached raw JSON in data/raw/pubchem/json/<CID>.json
"""

import os
import sys
import re
import json
import time
import urllib.parse
import requests
import pandas as pd
from pathlib import Path

BASE_DIR = Path(__file__).resolve().parent
COMPOUNDS_CSV = BASE_DIR / "data" / "unified" / "compounds.csv"
DEDUP_CSV = BASE_DIR / "data" / "unified" / "measurements_dedup.csv"
PUBCHEM_DIR = BASE_DIR / "data" / "raw" / "pubchem"
JSON_CACHE_DIR = PUBCHEM_DIR / "json"
OUT_CSV = PUBCHEM_DIR / "pubchem_experimental.csv"
REJECTS_CSV = PUBCHEM_DIR / "pubchem_experimental_rejects.csv"
CID_MAP_FILE = PUBCHEM_DIR / "inchikey_cid_map.json"

HEADERS = {
    "User-Agent": "UGRP-Hydrocarbons-Research/1.0 (academic research; contact: ugrp-hydrocarbons@gmail.com)"
}
RATE_LIMIT_DELAY = 0.25  # seconds (max 4-5 req/s)


SESSION = requests.Session()
SESSION.headers.update(HEADERS)


def clean_text(text):
    if not isinstance(text, str):
        return ""
    # Normalize unicode degree signs and apostrophes
    s = text.replace('\u00b0', ' deg ').replace('\u00ba', ' deg ').replace('&deg;', ' deg ')
    # Remove thousands separators in numbers (e.g. 4,040 -> 4040, 16,720 -> 16720)
    s = re.sub(r'(\d),(\d)', r'\1\2', s)
    # Normalize scientific notation like 5.21X10+4 or 388.8X10+5 to 5.21e4
    s = re.sub(r'([0-9.]+)\s*[xX]\s*10\+?([0-9+-]+)', r'\1e\2', s)
    return s


def parse_temperature_condition(text):
    """Extract temperature condition in Celsius if reported in text."""
    text = clean_text(text)
    # at 25 C / @ 20 deg C / (25 C)
    m_tc = re.search(r'(?:at|@|\()\s*([-+]?\d+(?:\.\d+)?)\s*(?:deg)?\s*C\b', text, re.IGNORECASE)
    if m_tc:
        return float(m_tc.group(1))
    # at 68 F / @ 77 deg F / (70 deg F)
    m_tf = re.search(r'(?:at|@|\()\s*([-+]?\d+(?:\.\d+)?)\s*(?:deg)?\s*F\b', text, re.IGNORECASE)
    if m_tf:
        return (float(m_tf.group(1)) - 32.0) * 5.0 / 9.0
    return None


def parse_flash_point(text):
    """
    Parse flash point string to Celsius and method.
    Returns: (val_c, unit, method, raw_val, raw_unit) or None
    """
    text = clean_text(text)
    lower = text.lower()
    method = 'unspecified'
    if any(k in lower for k in ['closed cup', 'c.c.', 'closed-cup', ' c.c', '(c.c.)', 'closed', 'tag closed', 'pensky']):
        method = 'closed_cup'
    elif any(k in lower for k in ['open cup', 'o.c.', 'open-cup', ' o.c', '(o.c.)', 'open', 'cleveland']):
        method = 'open_cup'

    # Check for direct Celsius: e.g. 46.0 deg C, 46 C, -18 deg C, 4.4 dec C, 134 deg. C
    m_c = re.search(r'([-+]?\d+(?:\.\d+)?)\s*(?:-|to)?\s*([-+]?\d+(?:\.\d+)?)?\s*(?:deg\.?|dec\.?)?\s*C\b', text, re.IGNORECASE)
    if m_c:
        val1 = float(m_c.group(1))
        val2 = float(m_c.group(2)) if m_c.group(2) else val1
        val_c = (val1 + val2) / 2.0
        return val_c, 'C', method, val1, 'deg C'

    # Check for Fahrenheit: e.g. 115 deg F, 134 deg. F, -4 F
    m_f = re.search(r'([-+]?\d+(?:\.\d+)?)\s*(?:-|to)?\s*([-+]?\d+(?:\.\d+)?)?\s*(?:deg\.?|dec\.?)?\s*F\b', text, re.IGNORECASE)
    if m_f:
        val1 = float(m_f.group(1))
        val2 = float(m_f.group(2)) if m_f.group(2) else val1
        val_f = (val1 + val2) / 2.0
        val_c = (val_f - 32.0) * 5.0 / 9.0
        return val_c, 'C', method, val_f, 'deg F'

    # Check for Kelvin: e.g. 319.15 K
    m_k = re.search(r'(\d+(?:\.\d+)?)\s*K\b', text)
    if m_k:
        val_k = float(m_k.group(1))
        val_c = val_k - 273.15
        return val_c, 'C', method, val_k, 'K'

    return None


def parse_boiling_point(text):
    """
    Parse boiling point string to Celsius and pressure condition in kPa.
    Returns: (val_c, unit, p_kpa, raw_val, raw_unit) or None
    """
    text = clean_text(text)
    p_kpa = None
    m_p = re.search(r'(?:at|@)\s*(\d+(?:\.\d+)?)\s*\[?\s*(?:mm\s*Hg|mmHg|torr)', text, re.IGNORECASE)
    if m_p:
        p_kpa = float(m_p.group(1)) * 0.133322
    else:
        m_pkpa = re.search(r'(?:at|@)\s*(\d+(?:\.\d+)?)\s*kPa', text, re.IGNORECASE)
        if m_pkpa:
            p_kpa = float(m_pkpa.group(1))

    m_c = re.search(r'([-+]?\d+(?:\.\d+)?)\s*(?:-|to)?\s*([-+]?\d+(?:\.\d+)?)?\s*(?:deg\.?|dec\.?)?\s*C\b', text, re.IGNORECASE)
    if m_c:
        val1 = float(m_c.group(1))
        val2 = float(m_c.group(2)) if m_c.group(2) else val1
        val_c = (val1 + val2) / 2.0
        return val_c, 'C', p_kpa, val1, 'deg C'

    m_f = re.search(r'([-+]?\d+(?:\.\d+)?)\s*(?:-|to)?\s*([-+]?\d+(?:\.\d+)?)?\s*(?:deg\.?|dec\.?)?\s*F\b', text, re.IGNORECASE)
    if m_f:
        val1 = float(m_f.group(1))
        val2 = float(m_f.group(2)) if m_f.group(2) else val1
        val_f = (val1 + val2) / 2.0
        val_c = (val_f - 32.0) * 5.0 / 9.0
        return val_c, 'C', p_kpa, val_f, 'deg F'

    m_k = re.search(r'(\d+(?:\.\d+)?)\s*K\b', text)
    if m_k:
        val_k = float(m_k.group(1))
        val_c = val_k - 273.15
        return val_c, 'C', p_kpa, val_k, 'K'

    # Pattern like "146.3 deg" or "135 DEG @ 10 MM HG" (defaults to Celsius)
    m_deg = re.search(r'([-+]?\d+(?:\.\d+)?)\s*(?:-|to)?\s*([-+]?\d+(?:\.\d+)?)?\s*(?:deg\.?|deg)\b', text, re.IGNORECASE)
    if m_deg:
        val1 = float(m_deg.group(1))
        val2 = float(m_deg.group(2)) if m_deg.group(2) else val1
        val_c = (val1 + val2) / 2.0
        return val_c, 'C', p_kpa, val1, 'deg'

    return None


def parse_vapor_pressure(text):
    """
    Parse vapor pressure string into list of (val_kpa, T_C, raw_val, raw_unit).
    """
    text = clean_text(text)
    parts = text.split(';')
    results = []

    for part in parts:
        part = part.strip()
        if not part:
            continue
        t_c = parse_temperature_condition(part)

        # 1. atm
        m_atm = re.search(r'(\d+(?:\.\d+)?(?:[eE][-+]?\d+)?)\s*atm\b', part, re.IGNORECASE)
        if m_atm:
            val = float(m_atm.group(1))
            results.append((val * 101.325, t_c, val, 'atm'))
            continue

        # 2. kPa
        m_kpa = re.search(r'(?:kPa\s*(?:at|@)\s*[-+]?\d+(?:\.\d+)?\s*(?:deg)?\s*[CF]:\s*|[:=]\s*|\b)(\d+(?:\.\d+)?(?:[eE][-+]?\d+)?)\s*kPa\b', part, re.IGNORECASE)
        if not m_kpa:
            m_kpa = re.search(r'Vapor\s+pressure\s*,\s*kPa[^:]*:\s*(\d+(?:\.\d+)?(?:[eE][-+]?\d+)?)', part, re.IGNORECASE)
        if m_kpa:
            val = float(m_kpa.group(1))
            results.append((val, t_c, val, 'kPa'))
            continue

        # 3. Pa
        m_pa = re.search(r'Vapor\s+pressure\s*,\s*Pa[^:]*:\s*(\d+(?:\.\d+)?(?:[eE][-+]?\d+)?)', part, re.IGNORECASE)
        if not m_pa:
            m_pa = re.search(r'(\d+(?:\.\d+)?(?:[eE][-+]?\d+)?)\s*Pa\b', part)
        if m_pa:
            val = float(m_pa.group(1))
            results.append((val / 1000.0, t_c, val, 'Pa'))
            continue

        # 4. mmHg / mm Hg / torr (including bracketed [mmHg] / [mm Hg])
        m_mm = re.search(r'(\d+(?:\.\d+)?(?:[eE][-+]?\d+)?)\s*(?:\[\s*mm\s*Hg\s*\]|\[\s*mmHg\s*\]|(?:mm\s*Hg|mmHg|torr)\b)', part, re.IGNORECASE)
        if m_mm:
            val = float(m_mm.group(1))
            results.append((val * 0.133322, t_c, val, 'mmHg'))
            continue

        # 5. bar / mbar
        m_mbar = re.search(r'(\d+(?:\.\d+)?(?:[eE][-+]?\d+)?)\s*mbar\b', part, re.IGNORECASE)
        if m_mbar:
            val = float(m_mbar.group(1))
            results.append((val * 0.1, t_c, val, 'mbar'))
            continue
        m_bar = re.search(r'(\d+(?:\.\d+)?(?:[eE][-+]?\d+)?)\s*bar\b', part, re.IGNORECASE)
        if m_bar:
            val = float(m_bar.group(1))
            results.append((val * 100.0, t_c, val, 'bar'))
            continue

        # 6. psi
        m_psi = re.search(r'(\d+(?:\.\d+)?(?:[eE][-+]?\d+)?)\s*psi\b', part, re.IGNORECASE)
        if m_psi:
            val = float(m_psi.group(1))
            results.append((val * 6.89476, t_c, val, 'psi'))
            continue

    return results


def parse_heat_of_combustion(text, mw=None):
    """
    Parse heat of combustion string to MJ/kg.
    Returns: (val_mj_kg, T_C, raw_val, raw_unit) or None
    """
    text = clean_text(text).replace(',', '')
    t_c = parse_temperature_condition(text)

    # 1. Joules/kg or J/kg (e.g. -388.8X10+5 Joules/kg, -408.41 X 105 J/kg)
    m_jkg = re.search(r'([-+]?\d+(?:\.\d+)?)\s*(?:[xX]\s*10\+?(\d+)|[eE]([-+]?\d+))\s*(?:Joules/kg|J/kg)\b', text, re.IGNORECASE)
    if m_jkg:
        base = float(m_jkg.group(1))
        exp = int(m_jkg.group(2) if m_jkg.group(2) is not None else m_jkg.group(3))
        val_j_kg = abs(base * (10 ** exp))
        return val_j_kg / 1e6, t_c, base, 'J/kg'

    m_jkg_simple = re.search(r'([-+]?\d+(?:\.\d+)?)\s*(?:Joules/kg|J/kg)\b', text, re.IGNORECASE)
    if m_jkg_simple:
        val_j_kg = abs(float(m_jkg_simple.group(1)))
        return val_j_kg / 1e6, t_c, float(m_jkg_simple.group(1)), 'J/kg'

    # 2. MJ/kg
    m_mj = re.search(r'([-+]?\d+(?:\.\d+)?)\s*MJ/kg\b', text, re.IGNORECASE)
    if m_mj:
        val = abs(float(m_mj.group(1)))
        return val, t_c, float(m_mj.group(1)), 'MJ/kg'

    # 3. kJ/kg (e.g. 41.2X10+3 kJ/kg, 40110 kJ/kg)
    m_kjkg_sci = re.search(r'([-+]?\d+(?:\.\d+)?)\s*[xX]\s*10\+?(\d+)\s*kJ/kg\b', text, re.IGNORECASE)
    if m_kjkg_sci:
        base = float(m_kjkg_sci.group(1))
        exp = int(m_kjkg_sci.group(2))
        val_kj = abs(base * (10 ** exp))
        return val_kj / 1000.0, t_c, base, 'kJ/kg'

    m_kjkg = re.search(r'([-+]?\d+(?:\.\d+)?)\s*kJ/kg\b', text, re.IGNORECASE)
    if m_kjkg:
        val = abs(float(m_kjkg.group(1)))
        return val / 1000.0, t_c, float(m_kjkg.group(1)), 'kJ/kg'

    # 4. J/kmol (e.g. -4.8229X10+09 J/kmol)
    m_jkmol = re.search(r'([-+]?\d+(?:\.\d+)?)\s*(?:[xX]\s*10\+?(\d+)|[eE]([-+]?\d+))\s*J/kmol\b', text, re.IGNORECASE)
    if m_jkmol and mw and mw > 0:
        base = float(m_jkmol.group(1))
        exp = int(m_jkmol.group(2) if m_jkmol.group(2) is not None else m_jkmol.group(3))
        val_j_kmol = abs(base * (10 ** exp))
        val_kj_mol = val_j_kmol / 1e6
        return val_kj_mol / mw, t_c, base, 'J/kmol'

    # 5. BTU/lb (1 BTU/lb = 0.002326 MJ/kg)
    m_btu = re.search(r'([-+]?\d+(?:\.\d+)?)\s*(?:BTU/lb|Btu/lb)\b', text, re.IGNORECASE)
    if m_btu:
        val = abs(float(m_btu.group(1)))
        val_mj_kg = val * 0.002326
        return val_mj_kg, t_c, float(m_btu.group(1)), 'BTU/lb'

    # 6. cal/g (1 cal/g = 0.004184 MJ/kg)
    m_calg = re.search(r'([-+]?\d+(?:\.\d+)?)\s*(?:cal/g)\b', text, re.IGNORECASE)
    if m_calg:
        val = abs(float(m_calg.group(1)))
        val_mj_kg = val * 0.004184
        return val_mj_kg, t_c, float(m_calg.group(1)), 'cal/g'

    # 7. kJ/mol (MJ/kg = kJ/mol / MW)
    m_kjmol = re.search(r'([-+]?\d+(?:\.\d+)?)\s*(?:kJ/mol|kJ/mole|KJ/mol)\b', text, re.IGNORECASE)
    if m_kjmol:
        raw_val = abs(float(m_kjmol.group(1)))
        if mw and mw > 0:
            val_mj_kg = raw_val / mw
            return val_mj_kg, t_c, raw_val, 'kJ/mol'

    # 8. kcal/mol (1 kcal = 4.184 kJ)
    m_kcalmol = re.search(r'([-+]?\d+(?:\.\d+)?)\s*(?:kcal/mol|kcal/mole)\b', text, re.IGNORECASE)
    if m_kcalmol:
        raw_val = abs(float(m_kcalmol.group(1)))
        if mw and mw > 0:
            val_mj_kg = (raw_val * 4.184) / mw
            return val_mj_kg, t_c, raw_val, 'kcal/mol'

    return None


def get_cid_map(df_compounds):
    """
    Get or load InChIKey -> CID map.
    """
    cid_map = {}
    if CID_MAP_FILE.exists():
        try:
            with open(CID_MAP_FILE, 'r', encoding='utf-8') as f:
                cid_map = json.load(f)
            print(f"Loaded {len(cid_map)} cached InChIKey -> CID mappings from {CID_MAP_FILE.name}")
        except Exception as e:
            print(f"Warning: Failed loading {CID_MAP_FILE}: {e}")

    # Check unmapped
    all_iks = df_compounds['inchikey'].dropna().unique().tolist()
    missing_iks = [ik for ik in all_iks if ik not in cid_map]
    if missing_iks:
        print(f"Resolving {len(missing_iks)} InChIKeys via PubChem PUG REST...")
        batch_size = 100
        for i in range(0, len(missing_iks), batch_size):
            batch = missing_iks[i:i + batch_size]
            url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/property/InChIKey/JSON"
            try:
                r = requests.post(url, data={"inchikey": ",".join(batch)}, headers=HEADERS, timeout=20)
                if r.status_code == 200:
                    for item in r.json().get("PropertyTable", {}).get("Properties", []):
                        cid_map[item["InChIKey"]] = item["CID"]
            except Exception as ex:
                print(f"Batch {i} error: {ex}")
            time.sleep(RATE_LIMIT_DELAY)

        # For remaining unmapped, try SMILES
        still_missing = [ik for ik in missing_iks if ik not in cid_map]
        if still_missing:
            print(f"Trying SMILES lookup for {len(still_missing)} remaining compounds...")
            for ik in still_missing:
                sm = df_compounds[df_compounds['inchikey'] == ik]['smiles'].iloc[0]
                url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/{urllib.parse.quote(sm)}/cids/JSON"
                try:
                    r = requests.get(url, headers=HEADERS, timeout=10)
                    if r.status_code == 200:
                        cids = r.json().get('IdentifierList', {}).get('CID', [])
                        if cids and cids[0] > 0:
                            cid_map[ik] = cids[0]
                except Exception:
                    pass
                time.sleep(RATE_LIMIT_DELAY)

        PUBCHEM_DIR.mkdir(parents=True, exist_ok=True)
        with open(CID_MAP_FILE, 'w', encoding='utf-8') as f:
            json.dump(cid_map, f, indent=2)
        print(f"Saved {len(cid_map)} InChIKey -> CID mappings to {CID_MAP_FILE}")

    return cid_map


def fetch_pug_view(cid):
    """
    Fetch PUG-View JSON for a CID with rate limiting and local caching.
    """
    cache_file = JSON_CACHE_DIR / f"{cid}.json"
    if cache_file.exists():
        try:
            with open(cache_file, 'r', encoding='utf-8') as f:
                return json.load(f)
        except Exception:
            pass

    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/{cid}/JSON?heading=Experimental+Properties"
    max_retries = 3
    for attempt in range(max_retries):
        try:
            r = SESSION.get(url, timeout=15)
            if r.status_code == 200:
                data = r.json()
                with open(cache_file, 'w', encoding='utf-8') as f:
                    json.dump(data, f)
                time.sleep(RATE_LIMIT_DELAY)
                return data
            elif r.status_code == 404:
                # No experimental properties for this compound
                data = {"status": 404, "message": "No data found"}
                with open(cache_file, 'w', encoding='utf-8') as f:
                    json.dump(data, f)
                time.sleep(RATE_LIMIT_DELAY)
                return data
            elif r.status_code in (429, 503):
                wait = (attempt + 1) * 3
                print(f"Rate limited or server busy (status {r.status_code}) on CID {cid}. Waiting {wait}s...")
                time.sleep(wait)
            else:
                time.sleep(RATE_LIMIT_DELAY)
                break
        except Exception as e:
            if attempt < max_retries - 1:
                time.sleep(2)
            else:
                print(f"Error fetching CID {cid}: {e}")

    return None


def extract_properties_from_json(data, cid, ik, sm, mw, compound_name):
    """
    Extract target experimental properties from PUG-View JSON.
    Returns: (list_of_parsed_records, list_of_rejects)
    """
    records = []
    rejects = []

    if not data or data.get("status") == 404:
        return records, rejects

    record = data.get("Record", {})
    title = record.get("RecordTitle", compound_name)
    top_refs = {r.get("ReferenceNumber"): r for r in record.get("Reference", [])}

    # Find Experimental Properties section
    exp_sections = []
    for s0 in record.get("Section", []):
        if s0.get("TOCHeading") == "Chemical and Physical Properties":
            for s1 in s0.get("Section", []):
                if s1.get("TOCHeading") == "Experimental Properties":
                    exp_sections = s1.get("Section", [])
                    break
        elif s0.get("TOCHeading") == "Experimental Properties":
            exp_sections = s0.get("Section", [])
            break

    target_headings = {
        "Flash Point": "flash_point_C",
        "Vapor Pressure": "vapor_pressure_kPa",
        "Heat of Combustion": "heat_of_combustion_MJ_kg",
        "Boiling Point": "boiling_point_C"
    }

    for sec in exp_sections:
        heading = sec.get("TOCHeading")
        if heading not in target_headings:
            continue
        canon_prop = target_headings[heading]

        for info in sec.get("Information", []):
            ref_num = info.get("ReferenceNumber")
            ref_text = ""
            if info.get("Reference"):
                ref_text = "; ".join(info["Reference"])
            elif ref_num in top_refs:
                top_ref = top_refs[ref_num]
                ref_text = top_ref.get("Citation") or top_ref.get("SourceName") or top_ref.get("Name") or ""

            val_dict = info.get("Value", {})
            raw_strings = [item.get("String") for item in val_dict.get("StringWithMarkup", []) if item.get("String")]
            
            # If string is empty but Number is present:
            if not raw_strings and "Number" in val_dict:
                for num in val_dict["Number"]:
                    unit = val_dict.get("Unit", "")
                    raw_strings.append(f"{num} {unit}")

            if not raw_strings:
                continue

            for raw_str in raw_strings:
                raw_str_clean = raw_str.strip()
                if not raw_str_clean:
                    continue

                if heading == "Flash Point":
                    res = parse_flash_point(raw_str_clean)
                    if res:
                        val_c, unit, method, raw_val, raw_unit = res
                        # Sanity filter
                        if -150 <= val_c <= 400:
                            records.append({
                                "cid": cid,
                                "inchikey": ik,
                                "smiles": sm,
                                "compound_name": title,
                                "property": canon_prop,
                                "value": round(val_c, 3),
                                "unit": unit,
                                "T_C": "",
                                "P_kPa": "",
                                "method": method,
                                "source": ref_text,
                                "raw_value": raw_val,
                                "raw_unit": raw_unit,
                                "raw_text": raw_str_clean,
                                "tier": "B"
                            })
                        else:
                            rejects.append({
                                "cid": cid, "inchikey": ik, "smiles": sm,
                                "property_heading": heading, "raw_text": raw_str_clean,
                                "reason": f"Value out of physical range: {val_c} C"
                            })
                    else:
                        rejects.append({
                            "cid": cid, "inchikey": ik, "smiles": sm,
                            "property_heading": heading, "raw_text": raw_str_clean,
                            "reason": "Failed to parse flash point"
                        })

                elif heading == "Boiling Point":
                    res = parse_boiling_point(raw_str_clean)
                    if res:
                        val_c, unit, p_kpa, raw_val, raw_unit = res
                        if -170 <= val_c <= 600:
                            records.append({
                                "cid": cid,
                                "inchikey": ik,
                                "smiles": sm,
                                "compound_name": title,
                                "property": canon_prop,
                                "value": round(val_c, 3),
                                "unit": unit,
                                "T_C": "",
                                "P_kPa": round(p_kpa, 3) if p_kpa is not None else "",
                                "method": "",
                                "source": ref_text,
                                "raw_value": raw_val,
                                "raw_unit": raw_unit,
                                "raw_text": raw_str_clean,
                                "tier": "B"
                            })
                        else:
                            rejects.append({
                                "cid": cid, "inchikey": ik, "smiles": sm,
                                "property_heading": heading, "raw_text": raw_str_clean,
                                "reason": f"Value out of physical range: {val_c} C"
                            })
                    else:
                        rejects.append({
                            "cid": cid, "inchikey": ik, "smiles": sm,
                            "property_heading": heading, "raw_text": raw_str_clean,
                            "reason": "Failed to parse boiling point"
                        })

                elif heading == "Vapor Pressure":
                    res_list = parse_vapor_pressure(raw_str_clean)
                    if res_list:
                        for val_kpa, t_c, raw_val, raw_unit in res_list:
                            if 0 < val_kpa < 10000:
                                records.append({
                                    "cid": cid,
                                    "inchikey": ik,
                                    "smiles": sm,
                                    "compound_name": title,
                                    "property": canon_prop,
                                    "value": round(val_kpa, 6),
                                    "unit": "kPa",
                                    "T_C": round(t_c, 2) if t_c is not None else "",
                                    "P_kPa": "",
                                    "method": "",
                                    "source": ref_text,
                                    "raw_value": raw_val,
                                    "raw_unit": raw_unit,
                                    "raw_text": raw_str_clean,
                                    "tier": "B"
                                })
                            else:
                                rejects.append({
                                    "cid": cid, "inchikey": ik, "smiles": sm,
                                    "property_heading": heading, "raw_text": raw_str_clean,
                                    "reason": f"Vapor pressure out of physical range: {val_kpa} kPa"
                                })
                    else:
                        rejects.append({
                            "cid": cid, "inchikey": ik, "smiles": sm,
                            "property_heading": heading, "raw_text": raw_str_clean,
                            "reason": "Failed to parse vapor pressure"
                        })

                elif heading == "Heat of Combustion":
                    res = parse_heat_of_combustion(raw_str_clean, mw=mw)
                    if res:
                        val_mj_kg, t_c, raw_val, raw_unit = res
                        # Typical hydrocarbon heat of combustion ~35 to 55 MJ/kg
                        if 20 <= val_mj_kg <= 70:
                            records.append({
                                "cid": cid,
                                "inchikey": ik,
                                "smiles": sm,
                                "compound_name": title,
                                "property": canon_prop,
                                "value": round(val_mj_kg, 4),
                                "unit": "MJ/kg",
                                "T_C": round(t_c, 2) if t_c is not None else "",
                                "P_kPa": "",
                                "method": "",
                                "source": ref_text,
                                "raw_value": raw_val,
                                "raw_unit": raw_unit,
                                "raw_text": raw_str_clean,
                                "tier": "B"
                            })
                        else:
                            rejects.append({
                                "cid": cid, "inchikey": ik, "smiles": sm,
                                "property_heading": heading, "raw_text": raw_str_clean,
                                "reason": f"Heat of combustion out of hydrocarbon range: {val_mj_kg} MJ/kg"
                            })
                    else:
                        rejects.append({
                            "cid": cid, "inchikey": ik, "smiles": sm,
                            "property_heading": heading, "raw_text": raw_str_clean,
                            "reason": f"Failed to parse heat of combustion (MW={mw})"
                        })

    return records, rejects


def main():
    print("=== T08: PubChem Experimental Properties Fetch ===")
    if not COMPOUNDS_CSV.exists():
        print(f"Error: {COMPOUNDS_CSV} not found!")
        sys.exit(1)

    df_compounds = pd.read_csv(COMPOUNDS_CSV)
    print(f"Total compounds in {COMPOUNDS_CSV.name}: {len(df_compounds)}")

    JSON_CACHE_DIR.mkdir(parents=True, exist_ok=True)

    # 1. InChIKey -> CID mapping
    cid_map = get_cid_map(df_compounds)

    # 2. Prioritize:
    # Check existing flash points in dedup
    missing_fp_iks = set()
    if DEDUP_CSV.exists():
        try:
            dedup = pd.read_csv(DEDUP_CSV, low_memory=False)
            pure_fp = dedup[(dedup['kind'] == 'pure') & (dedup['property'] == 'flash_point_C')]
            fp_iks = set()
            for s in pure_fp['comp_inchikeys'].dropna():
                try:
                    fp_iks.update(json.loads(s))
                except Exception:
                    pass
            missing_fp_iks = set(df_compounds['inchikey']) - fp_iks
            print(f"Compounds missing flash point in dedup: {len(missing_fp_iks)}")
        except Exception as e:
            print(f"Warning reading dedup: {e}")

    # Build ordered list of (ik, row)
    priority_list = []
    normal_list = []
    for _, row in df_compounds.iterrows():
        ik = row['inchikey']
        if ik in missing_fp_iks:
            priority_list.append((ik, row))
        else:
            normal_list.append((ik, row))

    ordered_compounds = priority_list + normal_list
    print(f"Ordered queue: {len(priority_list)} priority (missing FP) + {len(normal_list)} other = {len(ordered_compounds)}")

    # 3. Fetch PUG-View JSONs
    unique_cids_needed = set()
    for ik, _ in ordered_compounds:
        cid = cid_map.get(ik)
        if cid:
            unique_cids_needed.add(cid)

    cached_count = sum(1 for cid in unique_cids_needed if (JSON_CACHE_DIR / f"{cid}.json").exists())
    print(f"Unique CIDs needed: {len(unique_cids_needed)}. Already cached: {cached_count}. To fetch: {len(unique_cids_needed) - cached_count}")

    # Helper to parse and save
    def parse_and_save():
        all_records = []
        all_rejects = []
        for ik, row in ordered_compounds:
            cid = cid_map.get(ik)
            if not cid:
                continue
            cache_file = JSON_CACHE_DIR / f"{cid}.json"
            if not cache_file.exists():
                continue
            try:
                with open(cache_file, 'r', encoding='utf-8') as f:
                    data = json.load(f)
            except Exception:
                continue

            mw = row.get('mw')
            if pd.isna(mw) or mw <= 0:
                mw = None
            comp_name = row.get('formula', '')

            recs, rejs = extract_properties_from_json(data, cid, ik, row['smiles'], mw, comp_name)
            all_records.extend(recs)
            all_rejects.extend(rejs)

        df_records = pd.DataFrame(all_records)
        if not df_records.empty:
            df_records = df_records.drop_duplicates(subset=['cid', 'inchikey', 'property', 'value', 'T_C', 'P_kPa', 'raw_text'])
            df_records.to_csv(OUT_CSV, index=False)
            print(f"  [Checkpoint] Saved {len(df_records)} records to {OUT_CSV}")
            print(df_records['property'].value_counts().to_dict())

        df_rejects = pd.DataFrame(all_rejects)
        if not df_rejects.empty:
            df_rejects = df_rejects.drop_duplicates(subset=['cid', 'property_heading', 'raw_text'])
            df_rejects.to_csv(REJECTS_CSV, index=False)

        return df_records, df_rejects

    fetch_start = time.time()
    n_fetched = 0
    
    # Priority pass (compounds missing flash point)
    print("\n--- Phase 1: Priority compounds (missing flash points) ---")
    for idx, (ik, row) in enumerate(priority_list):
        cid = cid_map.get(ik)
        if not cid:
            continue
        cache_file = JSON_CACHE_DIR / f"{cid}.json"
        if not cache_file.exists():
            data = fetch_pug_view(cid)
            n_fetched += 1
            if n_fetched % 20 == 0 or n_fetched == 1:
                elapsed = time.time() - fetch_start
                rate = n_fetched / elapsed if elapsed > 0 else 0
                print(f"  Priority queue: fetched {n_fetched} CIDs ({rate:.1f} req/s)... CID {cid}")

    print(f"Phase 1 complete. Fetched {n_fetched} CIDs. Parsing priority checkpoint...")
    parse_and_save()

    # Normal pass (all remaining compounds)
    print("\n--- Phase 2: Systematic fetch of all remaining registry compounds ---")
    phase2_fetched = 0
    for idx, (ik, row) in enumerate(normal_list):
        cid = cid_map.get(ik)
        if not cid:
            continue
        cache_file = JSON_CACHE_DIR / f"{cid}.json"
        if not cache_file.exists():
            data = fetch_pug_view(cid)
            phase2_fetched += 1
            n_fetched += 1
            if phase2_fetched % 25 == 0:
                elapsed = time.time() - fetch_start
                rate = n_fetched / elapsed if elapsed > 0 else 0
                print(f"  Phase 2: fetched {phase2_fetched} additional CIDs (total {n_fetched}, {rate:.1f} req/s)... CID {cid}")
                parse_and_save()

    print(f"\nAll fetching complete. Total newly fetched: {n_fetched} in {time.time() - fetch_start:.1f}s")

    # Final parse and save
    print("\nFinal parsing and summary generation...")
    df_records, df_rejects = parse_and_save()


if __name__ == '__main__':
    main()
