"""v5 S6: diesel evidence from two published TABLES (CRC FACE-1 and NIST TN 1949).  Readings: results/v5/S6/README.md.

Subcommands (run with the venv python; the PDFs are gitignored under data/raw/diesel/):
  build  parse the PDF tables -> data/real_fuels/v5/diesel_*.csv      (needs pdfplumber; run with python -I)
  gate   apply R-S6-1 to the CSVs  -> results/v5/S6/gate_decision.csv
  score  R-S6-2 / R-S6-3 / R-S6-4 -> results/v5/S6/*.csv

There is NO reconstruction path: the loader below reads only rows whose composition came from a PDF table
(source_kind == "table") and the gate counts only such fuels.  Nothing in this file estimates a composition from
distillation, bulk descriptors or properties.
"""
from __future__ import annotations

import importlib.util
import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
V5 = ROOT / "data" / "real_fuels" / "v5"
OUT = ROOT / "results" / "v5" / "S6"
PDF_FACE = ROOT / "data" / "raw" / "diesel" / "CRC_FACE-1.pdf"
PDF_TN = ROOT / "data" / "raw" / "diesel" / "NIST_TN_1949.pdf"

FUELS = [f"FD{i}A" for i in range(1, 10)]
RHO_WATER_15_56 = 999.016           # repo convention (AJFTD ingest): density = SG(60/60 F) * 999.016 kg/m3 at 15.56 C
T_DENS_K = 288.71
T_VISC_K = 313.15
MAX_UNMAPPED = 0.10                 # same cut-off as v5 G1; declared before scoring (README note)
DCN_MIN_COVERAGE = 0.90             # DCN scored only if >= 90 % of mass has a verified pure DCN in pure_props_v4.csv

# B-table (GC-FIMS) leaf rows -> (hc_class used by map_fuel in eval_ajftd_gcxgc_v5.py, note)
FIMS_LEAF = {
    "isoparaffins": "iso_paraffin", "n-paraffins": "n_paraffin",
    "monocycloparaffins": "monocycloparaffin", "dicycloparaffins": "dicycloparaffin", "polycycloparaffins": "tricycloparaffin",
    "alkylbenzenes": "alkylbenzene", "benzocycloalkanes": "cycloaromatic", "benzodicycloalkanes": "benzodicycloalkane",
    "naphthalenes": "diaromatic", "biphenyls": "diaromatic", "naphthocycloalkanes": "naphthocycloalkane",
    "fluorenes": "fluorene", "phenanthrenes": "triaromatic", "phenanthrocyclolalkanes": "triaromatic",
    "pyrenes": "tetraaromatic", "crysenes": "tetraaromatic",
    "benzothiophenes": "aromatic_sulfur", "dibenzothiophenes": "aromatic_sulfur", "benzonaphthothiophenes": "aromatic_sulfur",
}
FIMS_ALIASES = {"monaaromatics": "monoaromatics"}       # typo in the PDF
FIMS_PARENT = {"saturates", "iso- + n-paraffins", "cycloparaffins", "aromatics", "monoaromatics", "diaromatics",
               "triaromatics", "tetraaromatics", "aromatic sulfur"}
CARBONS = list(range(8, 22))                            # C8 .. C21 (the first column is headed "C1 to C8")
# PIONA (Appendix C) columns after the carbon number: Monocyclo, Iso, Normal, Cyclic olefin, Branched olefin, Straight olefin, Monoaromatic, Totals
PIONA_COLS = [("monocycloparaffin", "monocyclo paraffins"), ("iso_paraffin", "iso paraffins"), ("n_paraffin", "normal paraffins"),
              ("alkene", "cyclic olefins"), ("alkene", "branched olefins"), ("alkene", "straight chain olefins"),
              ("alkylbenzene", "monoaromatics")]

NUM = re.compile(r"^-?\d+(\.\d*)?$")


# ----------------------------------------------------------------------------------------------- parsing
def _num(tok: str):
    return float(tok) if NUM.match(tok) else None


def parse_fims_row(line: str):
    """One GC-FIMS table row -> (label, z_number or None, [14 carbon values], tail tokens) or None if not a data row."""
    toks = line.split()
    lab = []
    while toks and not any(ch.isdigit() for ch in toks[0]):
        lab.append(toks.pop(0))
    if not lab:
        return None
    z = None
    if toks and re.fullmatch(r"-?\d+S?", toks[0]):       # integer Z number (no decimal point), e.g. 2, -6, -10S
        z = toks.pop(0)
    vals = [_num(t) for t in toks[:14]]
    if len(toks) < 14 or any(v is None for v in vals):
        return None
    return " ".join(lab), z, vals, toks[14:]


def parse_fims_page(text: str) -> dict:
    """-> {'label' (lower-case): [14 values]} for every full row of one Appendix B page, plus the printed page label."""
    lines = text.splitlines()
    out, total = {}, None
    for ln in lines:
        if ln.startswith("Total "):
            total = [float(t) for t in ln.split()[1:4] if _num(t) is not None]
            continue
        r = parse_fims_row(ln)
        if r:
            lab = r[0].lower()
            out[FIMS_ALIASES.get(lab, lab)] = dict(vals=r[2], z=r[1], tail=r[3])
    return dict(rows=out, total=total, printed=lines[0].strip())


def parse_piona_tables(text: str) -> dict:
    """Appendix C page -> {'C.n': {carbon: [8 values]}, with key 'Totals' per table}."""
    tables, cur = {}, None
    for ln in text.splitlines():
        m = re.match(r"Table (C\.\d+)\. PIONA Data for (FD\dA)", ln)
        if m:
            cur = m.group(2)
            tables[cur] = dict(id=m.group(1), rows={}, totals=None)
            continue
        if cur is None:
            continue
        toks = ln.split()
        if len(toks) == 9 and toks[0].isdigit() and all(NUM.match(t) for t in toks[1:]):
            tables[cur]["rows"][int(toks[0])] = [float(t) for t in toks[1:]]
        elif toks and toks[0] == "Totals" and len(toks) == 9:
            tables[cur]["totals"] = [float(t) for t in toks[1:]]
    return tables


def parse_table_52(text: str) -> pd.DataFrame:
    rows = []
    for ln in text.splitlines():
        t = ln.split()
        if len(t) == 9 and re.fullmatch(r"(FD\dA|Face5)", t[0]) and all(NUM.match(x) for x in t[1:]):
            lab = "FD5A" if t[0] == "Face5" else t[0]          # "Face5" is a typo for FD5A in the PDF
            v = [float(x) for x in t[1:]]
            rows.append(dict(fuel_id=lab, n_paraffins=v[0], paraffins=v[1], isoparaffins=v[2], cycloalkanes=v[3],
                             known_aromatics=v[4], unknown_aromatics=v[5], total_aromatics=v[6], total=v[7]))
    return pd.DataFrame(rows)


def row_values(text: str, start_pat: str, n: int = 9, after: str | None = None):
    """First line matching start_pat (after the anchor line 'after') whose remainder is n numbers."""
    lines = text.splitlines()
    i0 = 0
    if after:
        i0 = next(i for i, ln in enumerate(lines) if after in ln)
    for ln in lines[i0:]:
        if re.match(start_pat, ln):
            t = [x for x in ln.split() if NUM.match(x.replace(",", ""))]
            v = [float(x.replace(",", "")) for x in t]
            if len(v) >= n:
                return v[-n:]
    raise ValueError(f"row not found: {start_pat}")


# ----------------------------------------------------------------------------------------------- build
def _pages(pdf_path: Path) -> list[str]:
    import pdfplumber
    with pdfplumber.open(pdf_path) as pdf:
        return [(p.extract_text() or "") for p in pdf.pages]


def build() -> None:
    face, tn = _pages(PDF_FACE), _pages(PDF_TN)
    P = lambda n: face[n - 1]                              # 1-based PDF page
    comp_rows, checks = [], []
    # ---- GC-FIMS (>200 C) by type x carbon number: Appendix B, tables B.1..B.9 on PDF pages 185..193
    for k, fid in enumerate(FUELS):
        pg = 185 + k
        d = parse_fims_page(P(pg))
        assert f"Table B.{k + 1}. GC-FIMS Data for {fid}" in P(pg), (fid, pg)
        leaf_total = np.zeros(14)
        for lab, cls in FIMS_LEAF.items():
            if lab not in d["rows"]:
                raise ValueError(f"{fid}: row {lab} missing on page {pg}")
            vals = d["rows"][lab]["vals"]
            leaf_total += np.array(vals)
            for c, v in zip(CARBONS, vals):
                if v > 0:
                    comp_rows.append(dict(fuel_id=fid, boiling_part=">200C", method="GC-FIMS (D2786/D3239 types, Appendix B)",
                                          type_raw=lab, hc_class=cls, carbon_number=c, open_ended_low=(c == 8), mass_pct=v,
                                          source_table=f"B.{k + 1}", source_page=f"pdf p{pg} (printed {d['printed']})"))
        # internal consistency: parent rows vs sum of leaves, printed >200 total vs sum of carbon columns
        tot200, pionat, tsum = d["total"]
        checks.append(dict(fuel_id=fid, check="leaf_sum_vs_printed_total_200C", printed=tot200, computed=float(leaf_total.sum())))
        for par, leaves in (("saturates", ["isoparaffins", "n-paraffins", "monocycloparaffins", "dicycloparaffins", "polycycloparaffins"]),
                            ("aromatics", ["alkylbenzenes", "benzocycloalkanes", "benzodicycloalkanes", "naphthalenes", "biphenyls",
                                           "naphthocycloalkanes", "fluorenes", "phenanthrenes", "phenanthrocyclolalkanes", "pyrenes", "crysenes",
                                           "benzothiophenes", "dibenzothiophenes", "benzonaphthothiophenes"])):
            pv = np.array(d["rows"][par]["vals"])
            lv = sum(np.array(d["rows"][x]["vals"]) for x in leaves)
            for c, a, b in zip(CARBONS, pv, lv):
                if abs(a - b) > 0.025:
                    checks.append(dict(fuel_id=fid, check=f"parent_{par}_vs_leaves_C{c}", printed=float(a), computed=float(b)))
    # ---- PIONA (<200 C) by carbon number: Appendix C, PDF pages 197..201
    piona = {}
    for pg in range(197, 202):
        for fid, t in parse_piona_tables(P(pg)).items():
            piona[fid] = dict(t, page=f"pdf p{pg} (printed {P(pg).splitlines()[0].strip()})")
    for fid in FUELS:
        t = piona[fid]
        colsum = np.zeros(8)
        for n, vals in t["rows"].items():
            colsum += np.array(vals)
            for (cls, raw), v in zip(PIONA_COLS, vals[:7]):
                if v > 0:
                    comp_rows.append(dict(fuel_id=fid, boiling_part="<200C", method="PIONA (Appendix C)", type_raw=raw, hc_class=cls,
                                          carbon_number=n, open_ended_low=False, mass_pct=v, source_table=t["id"], source_page=t["page"]))
            checks.append(dict(fuel_id=fid, check=f"piona_row_sum_C{n}", printed=vals[7], computed=float(sum(vals[:7]))))
        checks.append(dict(fuel_id=fid, check="piona_column_totals_vs_printed_Totals_row", printed=t["totals"][7], computed=float(colsum[7])))
    comp = pd.DataFrame(comp_rows)
    comp.insert(1, "fuel_type", "FACE diesel")
    comp["source_kind"] = "table"
    comp["provenance"] = "crc_face1_tables_B_and_C"
    comp["data_tier"] = "T2-diesel"
    comp["double_count_group"] = "CRC_FACE1"
    comp.to_csv(V5 / "diesel_face_composition.csv", index=False)
    pd.DataFrame(checks).to_csv(OUT / "composition_internal_checks.csv", index=False)

    # ---- class-level Table 5.2 (kept as a cross-check of the type x carbon table, not used for mapping)
    t52 = parse_table_52(P(80))
    assert len(t52) == 9
    t52["source_table"], t52["source_page"] = "5.2", "pdf p80 (printed 58)"
    t52.to_csv(V5 / "diesel_face_class_table52.csv", index=False)

    # ---- properties (Table 2.2 pages 39-40, Table 3.5 page 53)
    prop = []

    def add(fid, name, val, unit, T_K, lab, method, raw, raw_unit, tab, page, note=""):
        prop.append(dict(fuel_id=fid, property=name, value=val, unit_std=unit, T_K=T_K, lab=lab, method=method, raw_value=raw,
                         raw_unit=raw_unit, source_table=tab, source_page=page, note=note, provenance="crc_face1_tables",
                         data_tier="T2-diesel", double_count_group="CRC_FACE1"))
    p39, p40, p53 = P(39), P(40), P(53)
    sgp = p39.split("Specific Gravity")[1]
    sg = {"CPChem": row_values(sgp, r"-CPChem C of A"), "ORNL/SwRI": row_values(sgp, r"-ORNL/SwRI"), "NCUT": row_values(sgp, r"-NCUT")}
    for lab, v in sg.items():
        for fid, x in zip(FUELS, v):
            add(fid, "density_kg_m3", x * RHO_WATER_15_56, "kg/m3", T_DENS_K, lab, "D4052 specific gravity x 999.016",
                x, "specific gravity (60/60 F by API-gravity consistency)", "2.2", "pdf p39 (printed 17)",
                "derived from specific gravity; density at 15.56 C")
    visc = {"ORNL/SwRI": row_values(p40.split("Kinematic Viscosity")[1], r"-ORNL/SwRI"),
            "NCUT": row_values(p40.split("Kinematic Viscosity")[1], r"-NCUT")}
    for lab, v in visc.items():
        for fid, x in zip(FUELS, v):
            add(fid, "viscosity_kin_cSt", x, "cSt", T_VISC_K, lab, "D445", x, "cSt", "2.2", "pdf p40 (printed 18)")
    # CPChem C of A viscosity: only FD6A and FD7A are filled in (read from the page image; columns are blank for the rest)
    for fid, x in (("FD6A", 3.09), ("FD7A", 1.55)):
        add(fid, "viscosity_kin_cSt", x, "cSt", T_VISC_K, "CPChem", "D445", x, "cSt", "2.2", "pdf p40 (printed 18)")
    fl = row_values(p40.split("Flash Point by Closed Cup")[1], r"-ORNL/SwRI")
    for fid, x in zip(FUELS, fl):
        add(fid, "flash_point_C", (x - 32) / 1.8, "C", np.nan, "ORNL/SwRI", "D93 closed cup", x, "F", "2.2", "pdf p40 (printed 18)")
    for fid, x in (("FD1A", 140.4), ("FD2A", 134.4), ("FD4A", 130.4), ("FD5A", 138.1), ("FD7A", 150.3)):    # blank for the other four
        add(fid, "flash_point_C", (x - 32) / 1.8, "C", np.nan, "CPChem", "D93 closed cup", x, "F", "2.2", "pdf p40 (printed 18)")
    cl = p40.split("Cloud Point")[1].split("Kinematic")[0]
    for lab, pat in (("ORNL/SwRI", r"-ORNL/SwRI"), ("NCUT", r"-NCUT")):
        for fid, x in zip(FUELS, row_values(cl, pat)):
            add(fid, "cloud_point_C", x, "C", np.nan, lab, "D2500", x, "C", "2.2", "pdf p40 (printed 18)")
    cn = p40.split("Cetane Number by Engine Method")[1].split("Cetane Index")[0]
    for lab, pat in (("CPChem", r"-CPChem C of A"), ("ORNL/SwRI", r"-ORNL/SwRI")):
        for fid, x in zip(FUELS, row_values(cn, pat)):
            add(fid, "cetane_number_D613", x, "-", np.nan, lab, "D613 engine", x, "-", "2.2", "pdf p40 (printed 18)")
    for ln in p53.splitlines():
        t = ln.split()
        if len(t) == 7 and t[0] in FUELS:
            add(t[0], "dcn_D6890", float(t[3]), "-", np.nan, "NREL", "D6890 IQT", float(t[3]), "-", "3.5", "pdf p53 (printed 31)")
    pd.DataFrame(prop).to_csv(V5 / "diesel_face_properties.csv", index=False)

    # ---- NIST TN 1949 surrogates: Table 1 (pdf p7), Table 11 (pdf p56)
    SM = {"n-hexadecane": "CCCCCCCCCCCCCCCC", "n-octadecane": "C" * 18, "n-eicosane": "C" * 20,
          "isohexadecane": "CC(CC(C)(C)C)CC(C)(C)CC(C)(C)C",            # = 2,2,4,4,6,8,8-heptamethylnonane (isocetane), see caveat in README results
          "2-methylheptadecane": "CC(C)" + "C" * 15, "n-butylcyclohexane": "CCCCC1CCCCC1",
          "1,3,5-triisopropylcyclohexane": "CC(C)C1CC(C(C)C)CC(C(C)C)C1", "trans-decalin": "C1CC[C@H]2CCCC[C@@H]2C1",
          "perhydrophenanthrene": "C1CCC2C(C1)CCC1CCCCC12", "1,2,4-trimethylbenzene": "Cc1ccc(C)c(C)c1",
          "1,3,5-triisopropylbenzene": "CC(C)c1cc(C(C)C)cc(C(C)C)c1", "tetralin": "c1ccc2c(c1)CCCC2", "1-methylnaphthalene": "Cc1cccc2ccccc12"}
    from rdkit import Chem
    from rdkit.Chem import rdinchi
    sur = []
    names = [n for n in SM]
    for ln in tn[6].splitlines():
        low = ln.lower()
        for nm in sorted(names, key=len, reverse=True):
            if low.startswith(nm):
                rest = ln[len(nm):].split()
                if len(rest) == 4:
                    for sid, tok in zip(("V0a", "V0b", "V1", "V2"), rest):
                        if tok != "-":
                            sur.append(dict(fuel_id=sid, component=nm, smiles=SM[nm], mole_fraction=float(tok)))
                break
    s = pd.DataFrame(sur)
    s["inchikey"] = [rdinchi.MolToInchiKey(Chem.MolFromSmiles(x)) for x in s.smiles]
    vals = None
    for ln in tn[55].splitlines():                  # the line that carries the 8 numbers of the viscosity row
        t = ln.split()
        if "Kinematic" in ln and len(t) >= 8 and all(NUM.match(x) for x in t[-8:]):
            vals = [float(x) for x in t[-8:]]
    assert vals is not None, "Table 11 viscosity row not found"
    exp = dict(zip(("V0a", "V0b", "V1", "V2"), vals[0::2]))
    mod = dict(zip(("V0a", "V0b", "V1", "V2"), vals[1::2]))
    s["nu40_exp_cSt_table11"] = s.fuel_id.map(exp)
    s["nu40_nist_model_cSt_table11"] = s.fuel_id.map(mod)
    s["source_table"], s["source_page"] = "NIST TN 1949 Table 1 (composition), Table 11 (viscosity)", "pdf p7 (printed 2) / pdf p56 (printed 51)"
    s["source_kind"], s["provenance"], s["data_tier"], s["double_count_group"] = "table", "nist_tn1949_tables", "T1-surrogate", "NIST_TN1949"
    s.to_csv(V5 / "diesel_avfl18a_surrogates.csv", index=False)
    print("built:", len(comp), "composition rows;", len(prop), "property rows;", len(s), "surrogate rows")


# ----------------------------------------------------------------------------------------------- loader (tables only)
def load_diesel_composition() -> pd.DataFrame:
    """FACE composition rows.  Refuses any row that is not table-sourced: there is no reconstruction path."""
    c = pd.read_csv(V5 / "diesel_face_composition.csv")
    if not (c.source_kind == "table").all():
        raise ValueError("diesel composition contains a non-table row; reconstructed diesel compositions are not allowed (R-S6-1)")
    return c


def load_diesel_properties() -> pd.DataFrame:
    return pd.read_csv(V5 / "diesel_face_properties.csv")


def load_surrogates() -> pd.DataFrame:
    s = pd.read_csv(V5 / "diesel_avfl18a_surrogates.csv")
    if not (s.source_kind == "table").all():
        raise ValueError("surrogate composition contains a non-table row")
    return s


# ----------------------------------------------------------------------------------------------- gate (R-S6-1)
TARGETS = ("density_kg_m3", "viscosity_kin_cSt", "flash_point_C", "cloud_point_C", "cetane_number_D613", "dcn_D6890")


def gate_table() -> pd.DataFrame:
    comp, prop = load_diesel_composition(), load_diesel_properties()
    rows = []
    for fid in FUELS:
        c = comp[comp.fuel_id == fid]
        total = float(c.mass_pct.sum())
        tg = sorted(set(prop[(prop.fuel_id == fid) & prop.property.isin(TARGETS) & prop.value.notna()].property))
        detail = (f"{c.type_raw.nunique()} hydrocarbon types x carbon number C{int(c.carbon_number.min())}-C{int(c.carbon_number.max())}; "
                  f"{(c.boiling_part == '>200C').sum()} GC-FIMS (>200 C) + {(c.boiling_part == '<200C').sum()} PIONA (<200 C) cells")
        ok_src = bool((c.source_kind == "table").all()) and len(c) > 0
        ok_level = c.type_raw.nunique() >= 4
        ok_sum = 95.0 <= total <= 105.0
        ok_tg = len(tg) > 0
        qual = ok_src and ok_level and ok_sum and ok_tg
        why = ("qualifies: table-sourced, class x carbon-number, sum in 95-105, measured targets present" if qual else
               "; ".join(r for r, ok in (("not table-sourced", ok_src), ("less than class-level", ok_level),
                                         (f"sum {total:.2f} outside 95-105", ok_sum), ("no measured target", ok_tg)) if not ok))
        rows.append(dict(fuel_id=fid, source="CRC FACE-1 (2010) Appendix B + C tables", composition_source_kind="table",
                         composition_detail=detail, composition_sum_pct=round(total, 3), targets_present="; ".join(tg),
                         in_gate_count=qual, decision="in" if qual else "out", reason=why))
    s = load_surrogates()
    for sid, g in s.groupby("fuel_id"):
        rows.append(dict(fuel_id=sid, source="NIST TN 1949 Table 1 / Table 11", composition_source_kind="table",
                         composition_detail=f"{len(g)} named pure compounds, exact mole fractions",
                         composition_sum_pct=round(100 * g.mole_fraction.sum(), 3), targets_present="viscosity_kin_cSt (40 C, Table 11)",
                         in_gate_count=False, decision="out of gate count",
                         reason="blend of pure compounds, not a real diesel fuel (README question A asks for real fuels); scored under R-S6-3"))
    return pd.DataFrame(rows)


def gate() -> pd.DataFrame:
    OUT.mkdir(parents=True, exist_ok=True)
    g = gate_table()
    g.to_csv(OUT / "gate_decision.csv", index=False)
    n = int(g.in_gate_count.sum())
    print(g[["fuel_id", "composition_sum_pct", "decision"]].to_string(index=False))
    print(f"gate count = {n} (>=5 -> {'MET' if n >= 5 else 'NOT met'})")
    return g


def gate_count(g: pd.DataFrame) -> int:
    """Only table-sourced real fuels count; surrogates never do."""
    return int(g[(g.composition_source_kind == "table") & g.in_gate_count & g.source.str.startswith("CRC FACE-1")].shape[0])


# ----------------------------------------------------------------------------------------------- scoring (R-S6-2, R-S6-3)
def _import(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def _models():
    from fuelblend.baselines.flash_unifac import _lookup
    from fuelblend.data.reconstruction_r3 import SubReps
    from fuelblend.mixture.common import CompoundRegistry, PurePropertyProvider
    vlb = _import(ROOT / "scripts" / "v4" / "build_viscosity_leaderboard_v4.py", "vlb")
    g1 = _import(ROOT / "scripts" / "v5" / "eval_ajftd_gcxgc_v5.py", "g1")
    reg = CompoundRegistry(ROOT / "data" / "unified" / "compounds.csv")
    prov = PurePropertyProvider(reg, ROOT / "data" / "cache" / "pure_grid.npz", verbose=False)
    reps = SubReps(reg, prov, vp_lookup=lambda s: _lookup(s, None))
    return reg, prov, reps, vlb.refutas_nu, g1


def stat(e) -> dict:
    e = np.asarray([x for x in e if np.isfinite(x)], float)
    return dict(n=len(e), bias=float(e.mean()) if len(e) else np.nan, mae=float(np.abs(e).mean()) if len(e) else np.nan,
                worst_over=float(e.max()) if len(e) else np.nan, worst_under=float(e.min()) if len(e) else np.nan)


def score() -> None:
    from rdkit.Chem import Descriptors
    from rdkit import Chem
    from fuelblend.baselines.flash_unifac import UnifacFallbackWarning, flash_point_unifac_detail
    import warnings
    g = gate_table()
    n_gate = gate_count(g)
    OUT.mkdir(parents=True, exist_ok=True)
    reg, prov, reps, refutas_nu, g1 = _models()
    margin = 0.1855823            # configs/uncertainty_v4.yaml viscosity_kin_cSt.lab_mixtures.raw_refutas.T_ge_283K applied_margin (read-only)
    prop = load_diesel_properties()
    comp = load_diesel_composition()

    # ---- R-S6-3 surrogates ---------------------------------------------------------------------------
    sur = load_surrogates()
    rows3 = []
    ik2i = reg.ik2i
    for sid, d in sur.groupby("fuel_id"):
        idx, x, name = [], [], []
        miss = 0.0
        for r in d.itertuples():
            i = ik2i.get(r.inchikey)
            ok = False
            if i is not None:
                prov.ensure([i])
                ok = bool(reg.usable[i]) and bool(prov.valid[i])
            if ok:
                idx.append(i)
                x.append(r.mole_fraction)
                name.append(r.component)
            else:
                miss += r.mole_fraction
        mw = {i: Descriptors.MolWt(Chem.MolFromSmiles(reg.smiles[i])) for i in idx}
        exp = float(d.nu40_exp_cSt_table11.iloc[0])
        nist_model = float(d.nu40_nist_model_cSt_table11.iloc[0])
        rec = dict(fuel_id=sid, n_components=len(d), mole_fraction_without_pure_viscosity_data=round(miss, 3), nu_exp_cSt=exp,
                   nist_model_cSt=nist_model, ln_err_nist_model=float(np.log(nist_model / exp)))
        if miss < 1e-9:
            w = np.array([xi * mw[i] for xi, i in zip(x, idx)])
            w = w / w.sum()
            nu, _, _, _ = refutas_nu(prov, np.array(idx), w, T_VISC_K)
            rec.update(status="evaluated", nu_refutas_cSt=nu, ln_err=float(np.log(nu / exp)), inside_margin=bool(abs(np.log(nu / exp)) <= margin))
        else:
            # sensitivity ONLY (not a score): drop the components without pure-viscosity data and renormalise
            w = np.array([xi * mw[i] for xi, i in zip(x, idx)])
            w = w / w.sum()
            nu, _, _, _ = refutas_nu(prov, np.array(idx), w, T_VISC_K)
            rec.update(status="not evaluable (pure eta(T) missing for: " + ", ".join(sorted(set(d.component) - set(name))) + ")",
                       nu_refutas_cSt=np.nan, ln_err=np.nan, inside_margin=np.nan,
                       sensitivity_covered_only_nu_cSt=nu, sensitivity_covered_only_ln_err=float(np.log(nu / exp)))
        rows3.append(rec)
    pd.DataFrame(rows3).to_csv(OUT / "surrogate_viscosity_R-S6-3.csv", index=False)

    # ---- R-S6-2 (only if the gate is met) -----------------------------------------------------------
    if n_gate < 5:
        (OUT / "R-S6-2_skipped.txt").write_text(f"gate count {n_gate} < 5: R-S6-2 not run (R-S6-4)\n")
        print("gate not met; R-S6-2 not run")
        return
    # mapping: reuse map_fuel / KIND of eval_ajftd_gcxgc_v5 (classes without a representative stay unmapped and are reported)
    map_rows, summ = [], []
    for fid in FUELS:
        c = comp[comp.fuel_id == fid]
        bins = c.groupby(["hc_class", "carbon_number"], as_index=False).mass_pct.sum().rename(columns={"carbon_number": "n_carbon"})
        agg, unmapped, lost = g1.map_fuel(bins, reps)
        tot_in = float(bins.mass_pct.sum())
        mapped_sum = sum(agg.values())
        for i, v in agg.items():
            map_rows.append(dict(fuel_id=fid, reg_index=i, inchikey=reps.ik[i], smiles=reps.smiles[i], mass_pct=v, fraction=v / mapped_sum))
        summ.append(dict(fuel_id=fid, input_sum_pct=tot_in, unmapped_fraction=unmapped, n_components=len(agg),
                         eligible=bool(agg) and unmapped <= MAX_UNMAPPED, unmapped_detail=json.dumps({k: round(v, 4) for k, v in lost.items()})))
    mp, sm = pd.DataFrame(map_rows), pd.DataFrame(summ)
    mp.to_csv(OUT / "mapped_representatives.csv", index=False)
    sm.to_csv(OUT / "mapping_summary.csv", index=False)

    # verified pure DCN (existing model: linear-by-volume of verified pure DCN, src/fuelblend/baselines/cetane.py + design/models_v4.py)
    pp = pd.read_csv(ROOT / "data" / "processed" / "v4" / "pure_props_v4.csv")
    pp = pp[(pp.prop == "DCN") & (pp.status == "verified")]
    from rdkit.Chem import rdinchi
    dcn_ik = {rdinchi.MolToInchiKey(Chem.MolFromSmiles(s)): float(v) for s, v in zip(pp.smiles, pp.value)}

    recs = []
    for s in sm.itertuples():
        gm = mp[mp.fuel_id == s.fuel_id]
        w = gm.fraction.to_numpy(float)
        idx = gm.reg_index.to_numpy(np.int64)
        iks, smi = list(gm.inchikey), list(gm.smiles)
        mws = np.array([Descriptors.MolWt(Chem.MolFromSmiles(x)) for x in smi])
        xm = w / mws
        xm = xm / xm.sum()
        t = prop[prop.fuel_id == s.fuel_id]
        base = dict(fuel_id=s.fuel_id, eligible=s.eligible, unmapped_fraction=s.unmapped_fraction)
        for r in t.itertuples():
            if r.property == "viscosity_kin_cSt":
                nu, _, _, _ = refutas_nu(prov, idx, w, T_VISC_K)
                recs.append(dict(base, property=r.property, lab=r.lab, model="refutas", measured=r.value, pred=nu, err=float(np.log(nu / r.value)),
                                 err_unit="ln(pred/meas)"))
            elif r.property == "density_kg_m3":
                _, rho, _, _ = refutas_nu(prov, idx, w, T_DENS_K)
                recs.append(dict(base, property=r.property, lab=r.lab, model="ideal_volume", measured=r.value, pred=rho, err=rho - r.value,
                                 err_unit="kg/m3"))
            elif r.property == "flash_point_C":
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore", UnifacFallbackWarning)
                    fp = flash_point_unifac_detail(smi, xm, use_unifac=False)
                recs.append(dict(base, property=r.property, lab=r.lab, model="ideal_le_chatelier", measured=r.value, pred=fp["flash_C"],
                                 err=fp["flash_C"] - r.value, err_unit="C", status=fp["status"]))
            elif r.property in ("cetane_number_D613", "dcn_D6890"):
                cov = float(sum(wi for wi, k in zip(w, iks) if k in dcn_ik))
                pred = np.nan
                if cov >= DCN_MIN_COVERAGE:
                    from fuelblend.baselines.cetane import blend_dcn_linear
                    rho_i = np.array([prov.lookup(np.array([[i]]), T_DENS_K)[1][0][0] for i in idx])
                    v = (w / rho_i)
                    v = v / v.sum()
                    pred = blend_dcn_linear(v, [dcn_ik.get(k, np.nan) for k in iks])
                recs.append(dict(base, property=r.property, lab=r.lab, model="linear_volume_verified_pure_DCN", measured=r.value, pred=pred,
                                 err=pred - r.value if np.isfinite(pred) else np.nan, err_unit="cetane units", verified_dcn_mass_coverage=cov))
    d = pd.DataFrame(recs)
    d.to_csv(OUT / "predictions_R-S6-2.csv", index=False)
    rows = []
    primary = {"viscosity_kin_cSt": "ORNL/SwRI", "density_kg_m3": "ORNL/SwRI", "flash_point_C": "ORNL/SwRI", "cetane_number_D613": "ORNL/SwRI", "dcn_D6890": "NREL"}
    for (pr, lab, mod), gg in d[d.eligible].groupby(["property", "lab", "model"]):
        st = stat(gg.err)
        row = dict(property=pr, lab=lab, model=mod, primary=(primary.get(pr) == lab), **st)
        if pr == "viscosity_kin_cSt":
            e = gg.err.to_numpy(float)
            row["fraction_inside_0.186_ln"] = float((np.abs(e) <= margin).mean()) if len(e) else np.nan
        if pr == "flash_point_C":
            e = gg.err.to_numpy(float)
            row["over_gt5C_rate"] = float((e > 5).mean()) if len(e) else np.nan
        rows.append(row)
    pd.DataFrame(rows).to_csv(OUT / "summary_R-S6-2.csv", index=False)
    print(pd.DataFrame(rows).to_string())
    # sensitivity (NOT the declared reading): all nine fuels regardless of unmapped mass; unmapped mass is dropped and the rest renormalised
    srows = []
    for (pr, lab, mod), gg in d.groupby(["property", "lab", "model"]):
        row = dict(property=pr, lab=lab, model=mod, **stat(gg.err))
        if pr == "viscosity_kin_cSt":
            e = gg.err.to_numpy(float)
            row["fraction_inside_0.186_ln"] = float((np.abs(e) <= margin).mean()) if len(e) else np.nan
        if pr == "flash_point_C":
            e = gg.err.to_numpy(float)
            row["over_gt5C_rate"] = float((e > 5).mean()) if len(e) else np.nan
        srows.append(row)
    pd.DataFrame(srows).to_csv(OUT / "summary_R-S6-2_sensitivity_all9.csv", index=False)
    print(pd.DataFrame(srows).to_string())
    print(sm.to_string())


# ----------------------------------------------------------------------------------------------- R-S6-5 transcription check
# Values below were read BY EYE from rendered page images (pypdfium2 at 2.2x), NOT from the text layer, and are compared here
# against the transcribed CSVs. Appendix B/C cells: (row label or column, carbon number, value on the page image).
IMG_B = {
    "B.1": [("isoparaffins", 12, 0.51), ("isoparaffins", 13, 3.00), ("isoparaffins", 14, 4.05), ("isoparaffins", 15, 3.85), ("n-paraffins", 16, 0.57),
            ("n-paraffins", 18, 0.77), ("monocycloparaffins", 14, 0.95), ("dicycloparaffins", 15, 1.18), ("polycycloparaffins", 15, 0.52),
            ("alkylbenzenes", 11, 0.80), ("naphthalenes", 11, 0.30), ("biphenyls", 15, 0.06), ("benzocycloalkanes", 10, 0.22)],
    "B.2": [("isoparaffins", 12, 1.14), ("isoparaffins", 13, 7.78), ("isoparaffins", 14, 5.48), ("n-paraffins", 12, 0.07), ("monocycloparaffins", 12, 2.21),
            ("alkylbenzenes", 10, 1.12), ("alkylbenzenes", 11, 2.86), ("naphthalenes", 11, 1.34), ("naphthalenes", 12, 0.73),
            ("dicycloparaffins", 21, 0.75), ("polycycloparaffins", 21, 1.39), ("benzothiophenes", 19, 0.16), ("biphenyls", 14, 0.12)],
    "B.3": [("isoparaffins", 12, 0.94), ("isoparaffins", 13, 3.13), ("isoparaffins", 14, 1.89), ("isoparaffins", 15, 2.81), ("n-paraffins", 12, 1.98),
            ("n-paraffins", 13, 1.34), ("monocycloparaffins", 13, 2.38), ("dicycloparaffins", 12, 2.13), ("alkylbenzenes", 10, 0.76),
            ("alkylbenzenes", 11, 2.58), ("alkylbenzenes", 12, 1.48), ("benzocycloalkanes", 12, 1.94)],
    "B.4": [("isoparaffins", 12, 0.17), ("isoparaffins", 13, 4.14), ("isoparaffins", 14, 7.75), ("isoparaffins", 15, 3.54), ("n-paraffins", 15, 0.01),
            ("monocycloparaffins", 13, 0.98), ("monocycloparaffins", 14, 2.64), ("monocycloparaffins", 15, 3.04), ("alkylbenzenes", 10, 1.78),
            ("alkylbenzenes", 11, 8.48), ("naphthalenes", 11, 2.25), ("polycycloparaffins", 21, 0.80)],
    "B.5": [("isoparaffins", 12, 2.30), ("isoparaffins", 13, 3.14), ("isoparaffins", 14, 2.25), ("n-paraffins", 12, 3.51), ("n-paraffins", 13, 2.35),
            ("n-paraffins", 14, 12.00), ("n-paraffins", 15, 7.01), ("monocycloparaffins", 12, 3.55), ("monocycloparaffins", 13, 3.48),
            ("alkylbenzenes", 11, 2.64), ("alkylbenzenes", 12, 1.61), ("polycycloparaffins", 16, 1.87)],
    "B.6": [("isoparaffins", 13, 1.57), ("isoparaffins", 14, 4.97), ("isoparaffins", 15, 8.16), ("isoparaffins", 16, 6.93), ("n-paraffins", 14, 1.42),
            ("n-paraffins", 15, 1.89), ("n-paraffins", 16, 1.34), ("polycycloparaffins", 16, 2.92), ("dicycloparaffins", 16, 2.27),
            ("biphenyls", 14, 0.36), ("alkylbenzenes", 14, 1.03), ("polycycloparaffins", 17, 2.75)],
    "B.7": [("n-paraffins", 12, 1.67), ("n-paraffins", 13, 1.22), ("n-paraffins", 14, 18.21), ("n-paraffins", 15, 10.09), ("n-paraffins", 16, 1.46),
            ("isoparaffins", 13, 1.45), ("biphenyls", 14, 2.24), ("biphenyls", 16, 2.79), ("biphenyls", 18, 0.97), ("naphthalenes", 16, 1.06),
            ("naphthalenes", 18, 1.33), ("alkylbenzenes", 14, 2.85)],
    "B.8": [("isoparaffins", 14, 1.83), ("isoparaffins", 15, 1.47), ("n-paraffins", 14, 16.44), ("n-paraffins", 15, 8.51), ("monocycloparaffins", 13, 0.88),
            ("monocycloparaffins", 19, 4.67), ("biphenyls", 14, 5.14), ("biphenyls", 16, 6.44), ("naphthalenes", 16, 2.53), ("naphthalenes", 18, 3.35),
            ("alkylbenzenes", 12, 1.36), ("alkylbenzenes", 14, 6.08)],
    "B.9": [("isoparaffins", 12, 0.14), ("isoparaffins", 13, 1.03), ("isoparaffins", 14, 2.99), ("isoparaffins", 15, 3.27), ("isoparaffins", 16, 3.33),
            ("isoparaffins", 17, 2.61), ("n-paraffins", 14, 1.22), ("monocycloparaffins", 13, 4.42), ("monocycloparaffins", 15, 5.12),
            ("polycycloparaffins", 16, 1.83), ("alkylbenzenes", 10, 0.23), ("alkylbenzenes", 11, 0.97)],
}
PIONA_NAME = {"monocyclo": "monocyclo paraffins", "iso": "iso paraffins", "normal": "normal paraffins", "cyclic": "cyclic olefins",
              "branched": "branched olefins", "arom": "monoaromatics"}
IMG_C = {   # (carbon number, column, value on the page image)
    "C.1": [(8, "monocyclo", .30), (9, "monocyclo", .57), (9, "iso", .52), (10, "iso", 1.96), (11, "iso", 4.67), (11, "normal", .41),
            (9, "arom", 26.99), (8, "arom", 5.02), (10, "arom", 4.21), (11, "branched", 1.63), (9, "normal", .29), (10, "monocyclo", .97)],
    "C.2": [(9, "arom", 34.54), (8, "arom", 7.14), (11, "iso", 4.52), (11, "branched", 2.43), (10, "arom", 4.45), (9, "normal", .22),
            (10, "iso", .97), (11, "monocyclo", .61), (10, "monocyclo", .10), (9, "iso", .09)],
    "C.3": [(9, "arom", 26.78), (8, "arom", 5.40), (10, "iso", 2.92), (11, "iso", 3.23), (10, "monocyclo", 1.65), (11, "normal", 1.40),
            (10, "normal", 1.66), (9, "monocyclo", .75), (9, "iso", .72), (8, "monocyclo", .41)],
    "C.4": [(9, "arom", 30.16), (8, "arom", 5.95), (10, "iso", 1.15), (11, "normal", 1.23), (10, "normal", .15), (9, "normal", .15),
            (8, "iso", .10), (9, "iso", .06), (11, "branched", .06), (10, "monocyclo", .02)],
    "C.5": [(9, "arom", 11.68), (8, "arom", 1.76), (10, "monocyclo", 2.93), (10, "iso", 3.26), (10, "normal", 2.91), (11, "normal", 2.92),
            (11, "iso", 2.04), (11, "monocyclo", 1.43), (9, "monocyclo", .98), (9, "iso", .76)],
    "C.6": [(9, "arom", 1.60), (8, "arom", .47), (10, "monocyclo", .79), (10, "iso", .94), (10, "normal", .84), (11, "normal", 1.32),
            (9, "monocyclo", .53), (9, "iso", .48), (8, "monocyclo", .37), (11, "iso", .25)],
    "C.7": [(9, "arom", 19.15), (10, "monocyclo", 1.40), (10, "iso", 1.63), (10, "normal", 1.50), (11, "normal", 1.99), (9, "monocyclo", .79),
            (9, "iso", .61), (9, "normal", .51), (11, "iso", .39), (8, "arom", .46)],
    "C.8": [(9, "arom", 2.92), (8, "arom", .13), (8, "iso", .10), (10, "normal", .05), (11, "normal", .05), (5, "iso", .06),
            (5, "cyclic", .01), (4, "normal", .02), (4, "iso", .01), (7, "iso", .03)],
    "C.9": [(9, "arom", 15.06), (8, "arom", 3.12), (9, "monocyclo", .90), (10, "monocyclo", .69), (10, "iso", 1.09), (10, "normal", .79),
            (11, "normal", 1.09), (9, "normal", .36), (8, "monocyclo", .32), (8, "iso", .16)],
}
IMG_T52 = [("FD1A", "n_paraffins", 6.78), ("FD1A", "paraffins", 59.06), ("FD1A", "isoparaffins", 52.28), ("FD1A", "cycloalkanes", 12.17),
           ("FD1A", "known_aromatics", 8.08), ("FD1A", "unknown_aromatics", 19.97), ("FD1A", "total_aromatics", 28.05), ("FD1A", "total", 99.28),
           ("FD3A", "n_paraffins", 10.49), ("FD3A", "total_aromatics", 54.56), ("FD5A", "n_paraffins", 31.39), ("FD5A", "total", 98.04),
           ("FD7A", "isoparaffins", 16.63), ("FD8A", "unknown_aromatics", 48.73), ("FD9A", "total", 98.30), ("FD6A", "cycloalkanes", 7.78),
           ("FD4A", "known_aromatics", 11.53), ("FD2A", "cycloalkanes", 14.66)]
_SG = {"CPChem": [.8084, .8037, .8401, .8355, .8086, .8411, .8375, .8682, .8465],
       "ORNL/SwRI": [.8084, .8037, .84, .8359, .8086, .8412, .8375, .8694, .8464],
       "NCUT": [.8076, .803, .8393, .8352, .8078, .8401, .8366, .8684, .8457]}
_VIS = {"ORNL/SwRI": [1.741, 1.96, 1.32, 1.785, 1.795, 3.097, 1.759, 3.278, 2.107],
        "NCUT": [1.56, 1.652, 1.317, 1.781, 1.587, 3.099, 1.581, 3.158, 2.104]}
_CLOUD = {"ORNL/SwRI": [-27, -51, -35, -50, -29, -25, -21, -20, -30], "NCUT": [-25, -59, -38, -59, -28, -28, -22, -19, -30]}
_FLASH = {"ORNL/SwRI": [135, 128, 141, 139, 131, 166, 149, 158, 140]}
_CN = {"CPChem": [29.9, 28.0, 32.0, 28.4, 54.2, 53.3, 44.3, 50.0, 45.0], "ORNL/SwRI": [30.7, 28.7, 30.7, 28.5, 55.0, 54.1, 45.9, 49.0, 43.5]}
_DCN = {"NREL": [35.3, 34.6, 33.8, 32.9, 54.9, 53.6, 45.4, 50.2, 44.6]}
IMG_PROP = (
    [(f, "density_kg_m3", lab, v) for lab, vs in _SG.items() for f, v in zip(FUELS, vs)]
    + [(f, "viscosity_kin_cSt", lab, v) for lab, vs in _VIS.items() for f, v in zip(FUELS, vs)]
    + [("FD6A", "viscosity_kin_cSt", "CPChem", 3.09), ("FD7A", "viscosity_kin_cSt", "CPChem", 1.55)]
    + [(f, "cloud_point_C", lab, v) for lab, vs in _CLOUD.items() for f, v in zip(FUELS, vs)]
    + [(f, "flash_point_C", lab, v) for lab, vs in _FLASH.items() for f, v in zip(FUELS, vs)]
    + [("FD1A", "flash_point_C", "CPChem", 140.4), ("FD2A", "flash_point_C", "CPChem", 134.4), ("FD4A", "flash_point_C", "CPChem", 130.4),
       ("FD5A", "flash_point_C", "CPChem", 138.1), ("FD7A", "flash_point_C", "CPChem", 150.3)]
    + [(f, "cetane_number_D613", lab, v) for lab, vs in _CN.items() for f, v in zip(FUELS, vs)]
    + [(f, "dcn_D6890", lab, v) for lab, vs in _DCN.items() for f, v in zip(FUELS, vs)]
)
_PROP_PAGE = {"density_kg_m3": ("Table 2.2 (specific gravity)", "pdf p39"), "viscosity_kin_cSt": ("Table 2.2", "pdf p40"),
              "cloud_point_C": ("Table 2.2", "pdf p40"), "flash_point_C": ("Table 2.2 (printed in F)", "pdf p40"),
              "cetane_number_D613": ("Table 3.5 / 2.2", "pdf p53 / p40"), "dcn_D6890": ("Table 3.5", "pdf p53")}
IMG_TN1 = [("V0a", "n-hexadecane", .278), ("V0a", "isohexadecane", .363), ("V0a", "trans-decalin", .148), ("V0a", "1-methylnaphthalene", .211),
           ("V0b", "n-octadecane", .235), ("V0b", "isohexadecane", .270), ("V0b", "tetralin", .209), ("V0b", "1,2,4-trimethylbenzene", .125),
           ("V0b", "1-methylnaphthalene", .161), ("V1", "n-hexadecane", .027), ("V1", "isohexadecane", .292), ("V1", "n-octadecane", .202),
           ("V1", "n-butylcyclohexane", .051), ("V1", "trans-decalin", .055), ("V1", "1,2,4-trimethylbenzene", .075), ("V1", "tetralin", .154),
           ("V2", "n-octadecane", .108), ("V2", "n-eicosane", .008), ("V2", "2-methylheptadecane", .073), ("V2", "n-butylcyclohexane", .191),
           ("V2", "1,3,5-triisopropylcyclohexane", .110), ("V2", "perhydrophenanthrene", .060), ("V2", "1,3,5-triisopropylbenzene", .147),
           ("V2", "tetralin", .164), ("V2", "1-methylnaphthalene", .139)]
IMG_TN11 = [("V0a", 2.452, 2.472), ("V0b", 2.303, 2.118), ("V1", 2.331, 2.172), ("V2", 2.378, 2.317)]


def transcription_check() -> pd.DataFrame:
    comp, prop, sur = load_diesel_composition(), load_diesel_properties(), load_surrogates()
    t52 = pd.read_csv(V5 / "diesel_face_class_table52.csv")
    rows = []

    def rec(table, page, item, img, csv):
        rows.append(dict(table=table, page=page, item=item, image_value=img, csv_value=round(float(csv), 6),
                         match=bool(abs(img - csv) < 5e-3 + 1e-9)))
    page_B = {t: 185 + i for i, t in enumerate(IMG_B)}
    for k, (tid, cells) in enumerate(IMG_B.items()):
        fid = FUELS[k]
        for lab, c, v in cells:
            x = comp[(comp.fuel_id == fid) & (comp.source_table == tid) & (comp.type_raw == lab) & (comp.carbon_number == c)].mass_pct.sum()
            rec(f"FACE-1 {tid}", f"pdf p{page_B[tid]}", f"{fid} {lab} C{c}", v, x)
    page_C = dict(zip(IMG_C, [197, 197, 198, 198, 199, 199, 200, 200, 201]))
    for k, (tid, cells) in enumerate(IMG_C.items()):
        fid = FUELS[k]
        for c, col, v in cells:
            q = comp[(comp.fuel_id == fid) & (comp.source_table == tid) & (comp.carbon_number == c) & (comp.type_raw == PIONA_NAME[col])]
            rec(f"FACE-1 {tid}", f"pdf p{page_C[tid]}", f"{fid} C{c} {col}", v, q.mass_pct.sum())
    for fid, col, v in IMG_T52:
        rec("FACE-1 Table 5.2", "pdf p80", f"{fid} {col}", v, t52[t52.fuel_id == fid][col].iloc[0])
    for fid, name, lab, v in IMG_PROP:
        q = prop[(prop.fuel_id == fid) & (prop.property == name) & (prop.lab == lab)]
        tb, pg = _PROP_PAGE[name]
        rec(f"FACE-1 {tb}", pg, f"{fid} {name} {lab} (as printed)", v, q.raw_value.iloc[0] if len(q) else np.nan)
    for sid, comp_name, v in IMG_TN1:
        q = sur[(sur.fuel_id == sid) & (sur.component == comp_name)]
        rec("NIST TN 1949 Table 1", "pdf p7", f"{sid} {comp_name}", v, q.mole_fraction.sum())
    for sid, e, m in IMG_TN11:
        s = sur[sur.fuel_id == sid]
        rec("NIST TN 1949 Table 11", "pdf p56", f"{sid} experimental nu40", e, s.nu40_exp_cSt_table11.iloc[0])
        rec("NIST TN 1949 Table 11", "pdf p56", f"{sid} model nu40", m, s.nu40_nist_model_cSt_table11.iloc[0])
    out = pd.DataFrame(rows)
    out.to_csv(OUT / "transcription_check.csv", index=False)
    print(out.groupby("table").agg(n=("match", "size"), n_match=("match", "sum")).to_string())
    print("total", len(out), "mismatches:", int((~out.match).sum()))
    print(out[~out.match].to_string())
    return out



if __name__ == "__main__":
    cmd = sys.argv[1] if len(sys.argv) > 1 else "gate"
    OUT.mkdir(parents=True, exist_ok=True)
    {"build": build, "gate": gate, "score": score, "check": transcription_check}[cmd]()
