"""Loader for literature smoke point / TSI / YSI data transcribed into data/raw/smoke_lit/ (v5).

Reads the five CSVs in ``data/raw/smoke_lit/`` (see SOURCES.md there) and returns one row per measured value with
InChIKeys resolved from SMILES via RDKit.  ``comp_inchikeys`` is the sorted InChIKey list (JSON) and ``fractions`` is
aligned to it, converted to fractions summing to 1 (percent / printed rounding handled here, never in the CSV).
``system_id`` = md5[:12] of the '|'-joined sorted InChIKeys (AGENTS.md section 6.2), or ``fuel:<source>:<sample>`` for
real fuels without a recipe.  ``is_hc`` is True when every component is a hydrocarbon (pure C/H); it is <NA> for real fuels.
Properties: smoke_point_mm, tsi, ysi.  ``qualifier`` keeps replicate / reference / reported_mean / literature_value /
anchor_assigned / text_conflict so a consumer can filter.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pandas as pd
from rdkit import Chem
from rdkit.Chem import inchi

ROOT = Path(__file__).resolve().parents[3]
RAW = ROOT / "data" / "raw" / "smoke_lit"
FILES = ["S3_mensch2009_smoke_tsi.csv", "S4_rainville2025_smoke.csv", "S5_botero_smoke_tsi.csv",
         "S6_boehm2022_smoke_tsi.csv", "S7_ysi_unified_hc.csv"]


def _mol(smiles: str):
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        raise ValueError(f"unparseable SMILES {smiles!r}")
    return m


def load_smoke_lit(raw_dir: Path | str = RAW, files=FILES) -> pd.DataFrame:
    raw_dir = Path(raw_dir)
    out = []
    for fn in files:
        df = pd.read_csv(raw_dir / fn, keep_default_na=False)
        for r in df.itertuples(index=False):
            smi = json.loads(r.smiles)
            names = json.loads(r.components)
            fr = [float(v) for v in json.loads(r.fractions)]
            mols = [_mol(s) for s in smi]
            iks = [inchi.MolToInchiKey(m) for m in mols]
            hc = all(a.GetSymbol() == "C" for m in mols for a in m.GetAtoms()) if mols else None
            order = sorted(range(len(iks)), key=lambda i: iks[i])
            if fr:
                tot = sum(fr)
                fr = [fr[i] / tot for i in order]
            iks_sorted = [iks[i] for i in order]
            if r.kind == "real_fuel":
                sid = f"fuel:{r.source_id}:{r.sample_id}"
            else:
                sid = hashlib.md5("|".join(iks_sorted).encode()).hexdigest()[:12]
            out.append(dict(source_id=r.source_id, doi=r.doi, table_or_page=r.table_or_page, sample_id=r.sample_id,
                            kind=r.kind, n_comp=len(iks), components=json.dumps([names[i] for i in order]),
                            comp_inchikeys=json.dumps(iks_sorted), fractions=json.dumps([round(v, 6) for v in fr]),
                            fraction_basis=r.fraction_basis, system_id=sid, is_hc=hc, property=r.property,
                            value=float(r.value), unit=r.unit,
                            T_C=(float(r.T_C) if str(r.T_C) != "" else float("nan")),
                            method=r.method, qualifier=r.qualifier, note=r.note))
    d = pd.DataFrame(out)
    d["is_hc"] = d["is_hc"].astype("boolean")
    return d
