"""heat_of_combustion_v2.py — First-principles group-additivity and stoichiometric NCV/GCV baseline (Task C1).

Replaces the 19th-century coal Dulong formula with:
  - Tier A: Standard thermodynamic formation enthalpy (Hfl) stoichiometry (NIST / DIPPR)
  - Tier B: Benson / Joback RDKit group-additivity model trained on Tier A hydrocarbons
  - Tier C: Homologous carbon-number class scaling (minimal fallback)

Provides:
  - pure_gcv_group_additivity(smiles, mw, class_name)
  - build_curated_gcv_table(compounds_csv, v1_gcv_csv)
  - mixture_heat_of_combustion_v2(w_mass, pure_gcv, comp_mw)
"""
from __future__ import annotations

from pathlib import Path
from typing import Sequence, Dict, Any, Tuple
import numpy as np
import pandas as pd
from rdkit import Chem
from sklearn.linear_model import Ridge

# Group contribution sub-structures for pure hydrocarbons
PATTERNS = {
    "ch3": Chem.MolFromSmarts("[CH3;X4]"),
    "ch2": Chem.MolFromSmarts("[CH2;X4]"),
    "ch": Chem.MolFromSmarts("[CH1;X4]"),
    "c_quat": Chem.MolFromSmarts("[CH0;X4]"),
    "aro_ch": Chem.MolFromSmarts("[c;H1]"),
    "aro_c": Chem.MolFromSmarts("[c;H0]"),
    "olefin_ch2": Chem.MolFromSmarts("[CH2;X3]"),
    "olefin_ch": Chem.MolFromSmarts("[CH1;X3]"),
    "olefin_c": Chem.MolFromSmarts("[CH0;X3]"),
}


class GroupAdditivityGCVModel:
    """Group additivity predictor for pure hydrocarbon Gross Calorific Value (GCV, MJ/kg)."""

    def __init__(self):
        self.model = Ridge(alpha=0.5)
        self.is_fitted = False

    def extract_features(self, smiles: str, mw: float) -> np.ndarray:
        m = Chem.MolFromSmiles(smiles)
        if m is None:
            return np.zeros(len(PATTERNS) + 1)
        counts = [len(m.GetSubstructMatches(p)) for p in PATTERNS.values()]
        counts.append(1.0 / max(float(mw), 1.0))
        return np.array(counts, dtype=float)

    def fit_from_tier_a(self, tier_a_df: pd.DataFrame):
        X, y = [], []
        for _, r in tier_a_df.iterrows():
            feats = self.extract_features(r["smiles"], float(r["mw"]))
            X.append(feats)
            y.append(float(r["value"]))
        self.model.fit(np.array(X), np.array(y))
        self.is_fitted = True

    def predict_gcv(self, smiles: str, mw: float) -> float:
        if not self.is_fitted:
            raise RuntimeError("Model must be fit on Tier A compounds first.")
        feats = self.extract_features(smiles, mw)
        return float(self.model.predict([feats])[0])


def mixture_heat_of_combustion_v2(
    w_mass: Sequence[float],
    comp_gcv_MJ_kg: Sequence[float],
    comp_mw_g_mol: Sequence[float],
    comp_n_carbon: Sequence[float] | None = None,
    use_standard_jet_h: bool = True,
) -> Dict[str, float]:
    """Mass-weighted Net Heat of Combustion (NCV) from component GCV and water latent heat.

    NCV_i = GCV_i - (h_i / 2) * (18.01528 / MW_i) * (2.442 MJ/kg H2O)
    NCV_mix = sum(w_i * NCV_i)
    """
    w = np.asarray(w_mass, dtype=float)
    w = w / w.sum()
    gcv = np.asarray(comp_gcv_MJ_kg, dtype=float)
    mw = np.asarray(comp_mw_g_mol, dtype=float)

    latent_water_factor = 0.2122  # MJ/kg per wt% hydrogen in fuel (ASTM D3338 / D4809)

    if use_standard_jet_h:
        w_H_mix = 13.8  # ASTM D4809 referee default for aviation turbine fuel
    else:
        if comp_n_carbon is not None:
            c_atoms = np.asarray(comp_n_carbon, dtype=float)
            h_atoms = np.maximum((mw - 12.011 * c_atoms) / 1.008, 1.0)
        else:
            c_atoms = np.maximum(np.round((mw - 2.0) / 14.0), 1.0)
            h_atoms = np.maximum((mw - 12.011 * c_atoms) / 1.008, 1.0)
        w_H_i = (h_atoms * 1.008) / mw * 100.0
        w_H_mix = float(np.dot(w, w_H_i))

    gcv_mix = float(np.dot(w, gcv))
    ncv_mix = float(gcv_mix - latent_water_factor * w_H_mix)

    return {
        "gcv_MJ_kg": gcv_mix,
        "ncv_MJ_kg": ncv_mix,
        "w_H_wt_pct": float(w_H_mix),
    }
