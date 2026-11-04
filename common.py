"""Shared helpers: hydrocarbon filter, class labels, composition conversions."""
from __future__ import annotations
import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors

RDLogger.DisableLog("rdApp.*")


def mol_from(smiles: str | None = None, inchi: str | None = None):
    if smiles:
        return Chem.MolFromSmiles(smiles)
    if inchi:
        return Chem.MolFromInchi(inchi)
    return None


def is_hydrocarbon(m) -> bool:
    """Only C/H, neutral, at least one carbon (excludes H2 etc.)."""
    if m is None:
        return False
    syms = [a.GetSymbol() for a in m.GetAtoms()]
    return (
        "C" in syms
        and all(s in ("C", "H") for s in syms)
        and all(a.GetFormalCharge() == 0 for a in m.GetAtoms())
    )


def hc_class(m) -> str:
    """aromatic | alkene_alkyne | cycloalkane | iso_paraffin | n_paraffin"""
    if any(a.GetIsAromatic() for a in m.GetAtoms()):
        return "aromatic"
    if any(b.GetBondType() != Chem.BondType.SINGLE for b in m.GetBonds()):
        return "alkene_alkyne"
    if m.GetRingInfo().NumRings() > 0:
        return "cycloalkane"
    branched = any(a.GetSymbol() == "C" and a.GetDegree() > 2 for a in m.GetAtoms())
    return "iso_paraffin" if branched else "n_paraffin"


def carbon_number(m) -> int:
    return sum(1 for a in m.GetAtoms() if a.GetSymbol() == "C")


def mol_weight(m) -> float:
    return float(Descriptors.MolWt(m))


# ---- composition conversions (match Section 0 of the methods document) ----
def mole_to_mass(x, M):
    x, M = np.asarray(x, float), np.asarray(M, float)
    return x * M / np.sum(x * M)


def mass_to_mole(w, M):
    w, M = np.asarray(w, float), np.asarray(M, float)
    return (w / M) / np.sum(w / M)


def mass_to_volume(w, rho):
    w, rho = np.asarray(w, float), np.asarray(rho, float)
    return (w / rho) / np.sum(w / rho)
