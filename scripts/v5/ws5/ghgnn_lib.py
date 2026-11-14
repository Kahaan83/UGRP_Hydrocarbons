"""WS5: batch inference and applicability-domain (AD) bookkeeping for GH-GNN (Sanchez Medina et al., Digital Discovery 2023).

Code and weights come from the cloned repository ``data/raw/ghgnn_code/repo`` (MIT licence, commit 0bbd2b3); two pure-python
dependencies are installed to ``data/raw/ghgnn_code/site`` (torch_geometric, mordredcommunity, tqdm, psutil, xxhash) plus a tiny pure-torch
``torch_scatter`` stand-in (no wheel exists for torch 2.15 CPU on Windows).  Nothing in the user's venv is modified.

What this module does differently from ``notebooks/GHGNN/ghgnn.py`` (the repo's own single-pair class), and why:
* the model and the AD pickles are loaded once, graphs are built once per molecule (the repo class reloads 30 000-row pickles per pair);
* the network itself, the atom/bond features and the data loader/collate are the repo's own code, unchanged;
* AD is the repo's definition (interpolation test on SMILES, otherwise the 10th-highest RDKFingerprint Tanimoto, and the ClassyFire
  class-pair count >= 25) except that the ClassyFire class is NOT fetched from the web: it is looked up in the repo's own
  ``data/interim/molecular_compounds.csv`` (InChIKey first block) and otherwise assigned by an explicit hydrocarbon rule, flagged ``class_source``.
All molecules here are hydrocarbons, so the hydrogen-bond counts (HBA/HBD) are 0 and ``inter_hb`` is 0 for every pair (asserted).
"""
from __future__ import annotations

import pickle
import sys
from functools import lru_cache
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
GH = ROOT / "data" / "raw" / "ghgnn_code"
REPO = GH / "repo"
for p in (str(Path(__file__).resolve().parent / "shims"), str(GH / "site"), str(REPO / "notebooks")):
    if p not in sys.path:
        sys.path.insert(0, p)

import torch  # noqa: E402
from rdkit import Chem, DataStructs  # noqa: E402
from rdkit.Chem import rdMolDescriptors  # noqa: E402
from mordred.Polarizability import APol, BPol  # noqa: E402
from mordred.TopoPSA import TopoPSA  # noqa: E402
from torch_geometric.data import Data  # noqa: E402

from GHGNN.GHGNN_architecture import GHGNN_model  # noqa: E402
from GHGNN.ghgnn import PairDataset_T, atom_features, bond_features, collate_T, get_bond_pair  # noqa: E402

NB = REPO / "notebooks" / "GHGNN"
HIDDEN, V_IN, E_IN, U_IN = 113, 37, 9, 3


def canon(smi: str) -> str:
    return Chem.MolToSmiles(Chem.MolFromSmiles(smi))


@lru_cache(maxsize=1)
def _model():
    m = GHGNN_model(V_IN, E_IN, U_IN, HIDDEN)
    m.load_state_dict(torch.load(NB / "GHGNN.pth", map_location="cpu"))
    return m.eval()


_GRAPHS: Dict[str, Data] = {}


def graph_of(smi: str) -> Data:
    """Per-molecule graph exactly as ``GH_GNN.sys2graphs`` builds it (inter_hb is set to 0 and asserted below)."""
    c = smi                      # graph built from the SMILES string as given (chiral-tag features depend on atom order, as in the repo)
    if c in _GRAPHS:
        return _GRAPHS[c]
    mol = Chem.MolFromSmiles(c)
    assert rdMolDescriptors.CalcNumHBA(mol) == 0 and rdMolDescriptors.CalcNumHBD(mol) == 0, "non-hydrocarbon passed to HC-only wrapper"
    ei = get_bond_pair(mol)
    ea = []
    for b in mol.GetBonds():
        ea.append(bond_features(b))
        ea.append(bond_features(b))
    g = Data(x=torch.tensor(np.array([atom_features(a) for a in mol.GetAtoms()]), dtype=torch.float),
             edge_index=torch.tensor(np.array(ei), dtype=torch.long),
             edge_attr=torch.tensor(np.array(ea), dtype=torch.float),
             ap=APol()(mol), bp=BPol()(mol), topopsa=TopoPSA()(mol), hb=0, inter_hb=0, y=torch.tensor([1]))
    _GRAPHS[c] = g
    return g


def predict_ln_gamma_inf(solute: Sequence[str], solvent: Sequence[str], T_C: Sequence[float], batch: int = 256) -> np.ndarray:
    """ln gamma_inf of ``solute`` in ``solvent`` at T (degC), the repo's own forward pass.  Training/validation T is degC."""
    model = _model()
    sv = [graph_of(s) for s in solvent]
    su = [graph_of(s) for s in solute]
    gt = [Data(x=torch.tensor(float(t), dtype=torch.float).reshape(1)) for t in T_C]
    loader = torch.utils.data.DataLoader(PairDataset_T(sv, su, gt), batch_size=batch, shuffle=False, collate_fn=collate_T)
    out = []
    with torch.no_grad():
        for b_solv, b_solu, b_T in loader:
            out.append(model(b_solv, b_solu, b_T).numpy().reshape(-1))
    return np.concatenate(out) if out else np.array([])


# --------------------------------------------------------------------------- applicability domain
@lru_cache(maxsize=1)
def _ad_data():
    with open(NB / "training_solutes.pickle", "rb") as f:
        solutes = set(pickle.load(f))
    with open(NB / "training_solvents.pickle", "rb") as f:
        solvents = set(pickle.load(f))
    with open(NB / "training_fp.pickle", "rb") as f:
        fps = pickle.load(f)
    with open(NB / "training_classes.pickle", "rb") as f:
        classes = pickle.load(f)
    comp = pd.read_csv(REPO / "data" / "interim" / "molecular_compounds.csv")
    comp["ik14"] = comp.inchikey.astype(str).str[:14]
    cmap = comp.dropna(subset=["class"]).drop_duplicates("ik14").set_index("ik14")["class"].to_dict()
    return solutes, solvents, fps, classes, cmap


def classyfire_class(smi: str) -> Tuple[str, str]:
    """(ClassyFire-style class name, source).  Source 'repo_table' = looked up in the GH-GNN repo's own molecular_compounds.csv;
    'rule' = hydrocarbon rule (no network call): no aromatic ring -> Saturated/Unsaturated hydrocarbons; one benzene ring -> 'Benzene and
    substituted derivatives'; naphthalene core -> 'Naphthalenes'; any other aromatic ring system -> 'Benzene and substituted derivatives'
    (ClassyFire would name a few of these, e.g. tetralins, differently: those pairs are then not matched to a training class and fail the AD)."""
    mol = Chem.MolFromSmiles(smi)
    ik = Chem.inchi.MolToInchiKey(mol)[:14]
    cmap = _ad_data()[4]
    if ik in cmap:
        return cmap[ik], "repo_table"
    ri = mol.GetRingInfo()
    arom_rings = [r for r in ri.AtomRings() if all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in r)]
    if not arom_rings:
        has_cc = any(b.GetBondTypeAsDouble() >= 2.0 for b in mol.GetBonds())
        return ("Unsaturated hydrocarbons" if has_cc else "Saturated hydrocarbons"), "rule"
    if Chem.MolFromSmiles("c1ccc2ccccc2c1") is not None and mol.HasSubstructMatch(Chem.MolFromSmiles("c1ccc2ccccc2c1")):
        return "Naphthalenes", "rule"
    return "Benzene and substituted derivatives", "rule"


@lru_cache(maxsize=None)
def _tanimoto10(smi: str) -> float:
    fps = _ad_data()[2]
    fp = Chem.RDKFingerprint(Chem.MolFromSmiles(smi))
    sims = sorted(DataStructs.BulkTanimotoSimilarity(fp, list(fps)))
    return float(sims[-10])


def ad_pair(solute: str, solvent: str) -> dict:
    """Repo AD definition for one ordered pair (solute in solvent): feasible iff interpolation (both SMILES seen in training as solute/solvent),
    or (10th-best Tanimoto >= 0.35 and ClassyFire class-pair count >= 25).  Unlike the repo (which returns 9999 for interpolation pairs) the
    real Tanimoto indicator and class-pair count are always computed and stored.  The Tanimoto indicator follows the repo's branch logic
    (min over both if neither molecule is in training, else the molecule that is not in training... as coded in the repo)."""
    solutes, solvents, _, classes, _ = _ad_data()
    cu, cv = canon(solute), canon(solvent)
    u_in, v_in = cu in solutes, cv in solvents
    interp = u_in and v_in
    if not u_in and not v_in:
        sim = min(_tanimoto10(cu), _tanimoto10(cv))
    elif v_in:                      # repo: solvent in training -> use the solute similarity
        sim = _tanimoto10(cu)
    else:
        sim = _tanimoto10(cv)
    if interp:                      # repo skips the indicator; we still compute it for reporting
        sim = min(_tanimoto10(cu), _tanimoto10(cv))
    cls_u, src_u = classyfire_class(cu)
    cls_v, src_v = classyfire_class(cv)
    n_cls = int(classes.get(cls_v + "_" + cls_u, 0))
    feasible = interp or ((sim >= 0.35) and (n_cls >= 25))
    return dict(solute_in_training=u_in, solvent_in_training=v_in, interpolation=interp, tanimoto10=sim, class_solute=cls_u,
                class_solvent=cls_v, class_source=f"{src_u}/{src_v}", n_class_pair=n_cls, ad_feasible=bool(feasible))
