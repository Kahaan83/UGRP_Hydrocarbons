#!/usr/bin/env python
"""v5: overlay on the v4 pure-constant audit (data/processed/v4/pure_props_v4.csv is read, never rewritten).

Two additions, both keep ruling R16 (`value` is always the LISTED value; a reference never overwrites it):

1. rho15 reference from MEASURED liquid densities.  The v4 audit compares the listed density with a value computed
   by the `thermo` library, which is itself wrong for some compounds (trans-decalin: 846.9 computed vs 869.8
   measured at 20 C).  Where data/unified/measurements_dedup.csv holds pure-liquid densities for the exact
   InChIKey (stereo included) the measured value becomes the reference: straight-line fit of density vs T over
   278-313 K at P <= 150 kPa, evaluated at 288.15 K.  Needs >= 3 temperatures from >= 2 distinct references and a
   point within 5 K of 288.15 K; rows flagged `conflict` are dropped.  Same 1 % threshold as v4.  The change can go
   either way (verified -> disputed too).  Rows that are disputed because two listings disagree are left alone.
2. Primary listings (configs/pure_constant_listings_v5.csv): a hand-entered literature value for a
   (compound, prop) that v4 marked `reference_only` (no listed value).  It becomes the listed value and is compared
   with the v4 reference under the v4 threshold.  Caveat: the library reference may trace back to the same
   measurement, so this is a transcription cross-check, not two independent measurements.

Outputs: data/processed/v5/pure_props_v5.csv (same columns as v4) and results/v5/P1/pure_props_changes_v5.csv.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

V4 = ROOT / "data/processed/v4/pure_props_v4.csv"
MEAS = ROOT / "data/unified/measurements_dedup.csv"
LISTINGS = ROOT / "configs/pure_constant_listings_v5.csv"
OUT = ROOT / "data/processed/v5/pure_props_v5.csv"
CHANGES = ROOT / "results/v5/P1/pure_props_changes_v5.csv"

T_RHO = 288.15
T_LO, T_HI = 278.0, 313.0
MAX_EXTRAP_K = 5.0
MIN_TEMPS, MIN_REFS = 3, 2
P_MAX_KPA = 150.0
THRESH = {"rho15": ("rel", 0.01), "dHfus": ("rel", 0.10), "Tm": ("abs", 2.0)}   # identical to v4 THRESH


def inchikey(smiles: str):
    from rdkit import Chem
    m = Chem.MolFromSmiles(smiles) if isinstance(smiles, str) else None
    return Chem.MolToInchiKey(m) if m is not None else None


def measured_rho15(meas: pd.DataFrame, ik: str):
    """(rho at 288.15 K, description) from measured pure-liquid densities, or None if the support rule fails."""
    q = meas[meas.comp_inchikeys == f'["{ik}"]']
    q = q[q.T_K.between(T_LO, T_HI) & (q.P_kPa.isna() | (q.P_kPa <= P_MAX_KPA)) & ~q.conflict.astype(bool)]
    if q.T_key.nunique() < MIN_TEMPS or q.ref.nunique() < MIN_REFS or (q.T_K - T_RHO).abs().min() > MAX_EXTRAP_K:
        return None
    slope, icpt = np.polyfit(q.T_K.to_numpy(float), q.value.to_numpy(float), 1)
    if not slope < 0:                       # a liquid density must fall with temperature
        return None
    how = "interpolated" if q.T_K.min() <= T_RHO <= q.T_K.max() else f"extrapolated {q.T_K.min() - T_RHO:.1f} K"
    return float(slope * T_RHO + icpt), (f"measured ThermoML liquid density, linear fit {how} to 288.15 K "
                                         f"({q.T_key.nunique()} temperatures {q.T_K.min():.2f}-{q.T_K.max():.2f} K, "
                                         f"{q.ref.nunique()} references)")


def _status(prop, listed, ref):
    kind, thr = THRESH[prop]
    d = (listed - ref) / ref if kind == "rel" else listed - ref
    return d, ("verified" if abs(d) <= thr else "disputed")


def build() -> tuple:
    pp = pd.read_csv(V4)
    meas = pd.read_csv(MEAS, low_memory=False)
    meas = meas[(meas.kind == "pure") & (meas.property == "density_kg_m3")]
    changes = []

    def log(i, why):
        changes.append(dict(compound=pp.at[i, "compound"], prop=pp.at[i, "prop"], change=why, **old,
                            new_status=pp.at[i, "status"], new_listed=pp.at[i, "listed_value"],
                            new_ref_value=pp.at[i, "ref_value"], new_ref_source=pp.at[i, "ref_source"],
                            new_delta=pp.at[i, "delta"]))

    for i in pp.index[(pp.prop == "rho15") & (pp.kind == "pure") & pp.listed_value.notna()]:
        if str(pp.at[i, "note"]).startswith("two listed values disagree"):
            continue
        ik = inchikey(pp.at[i, "smiles"])
        got = measured_rho15(meas, ik) if ik else None
        if got is None:
            continue
        old = dict(old_status=pp.at[i, "status"], old_ref_value=pp.at[i, "ref_value"], old_ref_source=pp.at[i, "ref_source"])
        d, st = _status("rho15", float(pp.at[i, "listed_value"]), got[0])
        pp.at[i, "ref_value"], pp.at[i, "ref_source"], pp.at[i, "delta"], pp.at[i, "status"] = got[0], got[1], d, st
        pp.at[i, "note"] = ("v5: reference is measured density (v4 used thermo-computed "
                            f"{old['old_ref_value']:.1f}); value kept as listed")
        log(i, "rho15 reference -> measured")

    for r in pd.read_csv(LISTINGS).itertuples():
        hit = pp.index[(pp.compound == r.compound) & (pp.prop == r.prop)]
        if len(hit) != 1 or pp.at[hit[0], "status"] != "reference_only":
            raise SystemExit(f"primary listing {r.compound}.{r.prop}: expected exactly one reference_only v4 row")
        i = hit[0]
        old = dict(old_status=pp.at[i, "status"], old_ref_value=pp.at[i, "ref_value"], old_ref_source=pp.at[i, "ref_source"])
        d, st = _status(r.prop, float(r.listed_value), float(pp.at[i, "ref_value"]))
        pp.at[i, "value"] = pp.at[i, "listed_value"] = float(r.listed_value)
        pp.at[i, "listed_source"] = pp.at[i, "source"] = r.listed_source
        pp.at[i, "delta"], pp.at[i, "status"] = d, st
        pp.at[i, "note"] = "v5 primary listing (configs/pure_constant_listings_v5.csv); " + str(r.note)
        log(i, "primary listing added")
    return pp, pd.DataFrame(changes)


def main() -> int:
    pp, ch = build()
    OUT.parent.mkdir(parents=True, exist_ok=True)
    CHANGES.parent.mkdir(parents=True, exist_ok=True)
    pp.to_csv(OUT, index=False)
    ch.to_csv(CHANGES, index=False)
    print(ch[["compound", "prop", "change", "old_status", "new_status", "new_listed", "old_ref_value", "new_ref_value",
              "new_delta"]].to_string())
    print(pd.crosstab(ch.old_status, ch.new_status).to_string())
    print(f"files: {OUT}\n       {CHANGES}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
