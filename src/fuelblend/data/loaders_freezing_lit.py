"""Literature freezing-point and fuel-composition data (v5 paper ingest): Vozka 2018 SI, ASCENT 065A (2024 slides, 2025 report), Kittel 2023.

Reads data/raw/freezing_lit/ (see SOURCES.md there):
  vozka2018_gcxgc.csv        GCxGC class x carbon-number tables (wt %) of Jet A, CAME, TALL, MFAT
  vozka2018_properties.csv   measured freezing points (Jet A, three neat HEFAs, 18 blends) and three densities
  ascent065A_freezing.csv    pure and mixture freezing points of named hydrocarbons (mole basis)
  kittel2023_fuels.csv       seven real fuels, bulk properties
  oakley2017_liquidus.csv    DSC liquidus temperatures of hexane/heptane + hexadecane mixtures and pure hexane (mole basis; stored in C)

``load_properties`` returns one row per measured value with InChIKeys resolved through RDKit (``comp_inchikeys`` sorted, ``fractions`` aligned),
``system_id`` = md5[:12] of the sorted InChIKey set (named compounds) or ``fuel:<source>:<sample>`` (real fuels, blends of real fuels), ``is_hc``.
Values are returned as printed (the paper's unit in ``unit``); ``qualifier`` carries censoring ('<', '>') and caveats.

``vozka_as_ajftd_composition`` / ``blend_composition`` give the Vozka fuels in the layout of data/real_fuels/v5/ajftd_gcxgc_composition.csv so they
can be reconstructed exactly like the AJFTD GCxGC fuels (class names mapped in ``CLASS_MAP``).
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pandas as pd
from rdkit import Chem
from rdkit.Chem import inchi

ROOT = Path(__file__).resolve().parents[3]
RAW = ROOT / "data" / "raw" / "freezing_lit"
PROPERTY_FILES = ("vozka2018_properties.csv", "ascent065A_freezing.csv", "kittel2023_fuels.csv",
                  "oakley2017_liquidus.csv")

# paper class -> (AJFTD hc_class as in ajftd_gcxgc_composition.csv, mixture.common.CLASSES entry, reconstruction_r3 subkind)
CLASS_MAP = {
    "n-paraffins": ("n_paraffin", "n_paraffin", "n_paraffin"),
    "isoparaffins": ("iso_paraffin", "iso_paraffin", "iso_paraffin"),
    "monocycloparaffins": ("monocycloparaffin", "cycloalkane", "cyclo_mono"),
    "di- and tricycloparaffins": ("dicycloparaffin", "cycloalkane", "cyclo_di"),   # paper merges di and tri; filed under di
    "alkylbenzenes": ("alkylbenzene", "aromatic", "alkylbenzene"),
    "cycloaromatics": ("cycloaromatic", "aromatic", "indan_tetralin"),
    "alkylnaphthalenes": ("diaromatic", "aromatic", "aromatic_di"),               # AJFTD files naphthalenes under diaromatic
}
NEAT_FUELS = ("Jet A", "CAME", "TALL", "MFAT")


def _inchikey(smiles: str) -> str:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"unparseable SMILES {smiles!r}")
    return inchi.MolToInchiKey(mol)


def _is_hc(smiles: str) -> bool:
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    return all(a.GetSymbol() in ("C", "H") for a in mol.GetAtoms())


def _slug(s: str) -> str:
    return str(s).replace(" ", "_")


def _resolve(r: pd.Series) -> pd.Series:
    comps, smis, fr = json.loads(r["components"]), json.loads(r["smiles"]), json.loads(r["fractions"])
    if smis:
        iks = [_inchikey(s) for s in smis]
        order = sorted(range(len(iks)), key=lambda i: iks[i])
        f = [float(fr[i]) for i in order]
        if r["fraction_basis"] != "none" and sum(f) > 0:
            f = [round(v / sum(f), 6) for v in f]
        return pd.Series(dict(comp_inchikeys=json.dumps([iks[i] for i in order]), fractions_norm=json.dumps(f),
                              system_id=hashlib.md5("|".join(iks[i] for i in order).encode()).hexdigest()[:12],
                              is_hc=all(_is_hc(s) for s in smis)))
    return pd.Series(dict(comp_inchikeys="[]", fractions_norm=json.dumps([float(v) for v in fr]),
                          system_id=f"fuel:{r['source_id'].split('_')[0]}:{_slug(r['sample_id'])}", is_hc=True))


def load_properties(raw: Path = RAW) -> pd.DataFrame:
    d = pd.concat([pd.read_csv(raw / f, keep_default_na=False) for f in PROPERTY_FILES], ignore_index=True)
    d["value"] = pd.to_numeric(d["value"])
    d["T_C"] = pd.to_numeric(d["T_C"].replace("", float("nan")))
    return pd.concat([d, d.apply(_resolve, axis=1)], axis=1)


def load_vozka_gcxgc(raw: Path = RAW) -> pd.DataFrame:
    d = pd.read_csv(raw / "vozka2018_gcxgc.csv")
    m = d["class"].map(CLASS_MAP)
    d["hc_class"], d["common_class"], d["subkind"] = (m.str[0], m.str[1], m.str[2])
    return d


def vozka_as_ajftd_composition(raw: Path = RAW) -> pd.DataFrame:
    """The four neat fuels in the layout of data/real_fuels/v5/ajftd_gcxgc_composition.csv (zero bins kept, as printed)."""
    d = load_vozka_gcxgc(raw)
    return pd.DataFrame(dict(fuel_id=d.fuel_id, fuel_label=d.fuel_id, hc_class=d.hc_class, n_carbon=d.n_carbon, cn_open_ended=False,
                             metric_raw=d["class"] + " C" + d.n_carbon.astype(str), provenance="real_fuel_gcxgc_bins",
                             mass_pct=d.value))


def blend_composition(sample_id: str, raw: Path = RAW) -> pd.DataFrame:
    """Linear blend of two neat fuels (e.g. 'C-30' = 70 % Jet A + 30 % CAME) in the same layout.

    The SI does not state whether the blend percentages are volume or mass; this mixes the wt % tables linearly, which is exact for a
    mass blend and a ~2-3 % relative approximation for a volume blend (neat densities 0.757-0.761 g/cm3; Jet A density not tabulated)."""
    p = load_properties(raw)
    r = p[(p.source_id == "vozka2018_si") & (p.sample_id == sample_id) & (p.property == "freezing_point_C")]
    if r.empty:
        raise KeyError(sample_id)
    comps, fr = json.loads(r.iloc[0].components), json.loads(r.iloc[0].fractions)
    if not comps:
        return vozka_as_ajftd_composition(raw).query("fuel_id == @sample_id").reset_index(drop=True)
    base = vozka_as_ajftd_composition(raw)
    parts = [base[base.fuel_id == c].assign(mass_pct=lambda x, w=w: x.mass_pct * w) for c, w in zip(comps, fr)]
    out = pd.concat(parts).groupby(["hc_class", "n_carbon", "cn_open_ended", "provenance"], as_index=False).mass_pct.sum()
    out.insert(0, "fuel_id", sample_id)
    out.insert(1, "fuel_label", sample_id)
    out["metric_raw"] = out.hc_class + " C" + out.n_carbon.astype(str)
    return out[["fuel_id", "fuel_label", "hc_class", "n_carbon", "cn_open_ended", "metric_raw", "provenance", "mass_pct"]]
