"""
fuelblend.mixture.cold_flow
===========================
Composition-level helpers for T9 that need no trained model:

  * wax_appearance_C   ideal-solubility (Schroeder-van Laar) wax appearance temperature, differentiable.
                       PROXY for freezing / cloud point. Each solute i with melting data precipitates when
                           ln x_i = -(dHf_i / R) (1/T - 1/Tm_i);   WAT = smooth-max_i T_i.
                       Compounds WITHOUT melting data are treated as never precipitating (reported by T9).
  * thermal_arrays     Tm / dHfus per compound: built-in n-alkane table (approximate literature values,
                       C8-C22 tabulated, Broadhurst fit beyond; total dH ~ 4.0 n - 10.6 kJ/mol), a few
                       aromatic solids, and an optional override file data/unified/thermal.csv with columns
                       inchikey, tm_K, hfus_kJ_mol (RECOMMENDED: verify/extend against NIST or CRC).
  * class / carbon-number one-hot matrices and a reference-composition distance
    (class marginals + Wasserstein-1 on the carbon-number distribution).

Install at: fuelblend_data/src/fuelblend/mixture/cold_flow.py
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import torch

R_GAS = 8.314462618
CLASSES = ["aromatic", "alkene_alkyne", "cycloalkane", "iso_paraffin", "n_paraffin"]
NC_MIN, NC_MAX = 4, 40

ALKANE_TM = {8: 216.4, 9: 219.7, 10: 243.5, 11: 247.6, 12: 263.6, 13: 267.8, 14: 279.0, 15: 283.1,
             16: 291.3, 17: 295.1, 18: 301.3, 19: 305.0, 20: 309.6, 21: 313.2, 22: 317.2}   # K, approximate
SPECIAL_SOLIDS = {                      # canonical SMILES -> (Tm K, dHfus J/mol), approximate
    "Cc1cc(C)c(C)cc1C": (352.4, 21.0e3),    # 1,2,4,5-tetramethylbenzene (durene)
    "c1ccc2ccccc2c1": (353.4, 19.0e3),      # naphthalene
    "Cc1ccc(C)cc1": (286.4, 17.1e3),        # p-xylene
    "c1ccccc1": (278.7, 9.87e3),            # benzene
}


def compound_meta(compounds_csv):
    """Same ordering as CompoundRegistry (drop_duplicates on inchikey, reset_index)."""
    from rdkit import Chem
    df = pd.read_csv(compounds_csv).drop_duplicates("inchikey").reset_index(drop=True)
    cls = df["cls"].astype(str).tolist()
    nc = df["n_carbon"].to_numpy(float) if "n_carbon" in df else np.full(len(df), np.nan)
    for i in np.where(~np.isfinite(nc))[0]:
        m = Chem.MolFromSmiles(df["smiles"].iat[i])
        nc[i] = sum(a.GetAtomicNum() == 6 for a in m.GetAtoms()) if m is not None else np.nan
    return df, cls, nc


def thermal_arrays(df, cls, nc, override_csv=None):
    from rdkit import Chem
    n = len(df)
    tm, hf = np.full(n, np.nan), np.full(n, np.nan)
    special = {Chem.CanonSmiles(k): v for k, v in SPECIAL_SOLIDS.items()}
    for i in range(n):
        if cls[i] == "n_paraffin" and np.isfinite(nc[i]) and nc[i] >= 8:
            k = int(nc[i])
            tm[i] = ALKANE_TM.get(k, 414.3 * (k - 1.5) / (k + 5.0))
            hf[i] = (4.0 * k - 10.6) * 1e3
            continue
        try:
            v = special.get(Chem.CanonSmiles(df["smiles"].iat[i]))
        except Exception:
            v = None
        if v:
            tm[i], hf[i] = v
    if override_csv and Path(override_csv).exists():
        o = pd.read_csv(override_csv).set_index("inchikey")
        for i, ik in enumerate(df["inchikey"]):
            if ik in o.index:
                tm[i], hf[i] = float(o.at[ik, "tm_K"]), float(o.at[ik, "hfus_kJ_mol"]) * 1e3
    return tm, hf


def wax_appearance_C(x, tm, hf, tau=0.5):
    """x mole fractions (B,K); tm (K,) K, hf (K,) J/mol (NaN = no data). Returns WAT in degC (B,).
    Smooth max with temperature tau [K] keeps gradients alive (overestimates by < tau*ln(#solutes))."""
    avail = torch.isfinite(tm)
    tm_s = torch.where(avail, tm, torch.full_like(tm, 300.0))
    hf_s = torch.where(avail, hf, torch.full_like(hf, 1.0e4))
    invT = 1.0 / tm_s - R_GAS * torch.log(x.clamp_min(1e-9)) / hf_s
    T = torch.where(avail[None, :], 1.0 / invT.clamp_min(1e-4), torch.full_like(invT, 100.0))
    return tau * torch.logsumexp(T / tau, -1) - 273.15


def class_matrix(cls_list):
    return (np.asarray(cls_list)[:, None] == np.asarray(CLASSES)[None, :]).astype(np.float32)


def carbon_matrix(nc):
    k = np.clip(np.nan_to_num(nc, nan=NC_MIN), NC_MIN, NC_MAX).astype(int) - NC_MIN
    m = np.zeros((len(nc), NC_MAX - NC_MIN + 1), np.float32)
    m[np.arange(len(nc)), k] = 1.0
    return m


def load_reference(path, reg, cls, nc):
    """CSV with inchikey or smiles, and mass_frac or wt_pct. Returns (r_cls (5,), r_nc (Nb,)) as mass fractions
    over the matched compounds (unmatched mass is reported and dropped)."""
    from rdkit import Chem
    d = pd.read_csv(path)
    col = "mass_frac" if "mass_frac" in d else "wt_pct"
    wts = d[col].to_numpy(float) / (1.0 if col == "mass_frac" else 100.0)
    r_cls, r_nc, matched = np.zeros(len(CLASSES)), np.zeros(NC_MAX - NC_MIN + 1), 0.0
    for row, wv in zip(d.itertuples(index=False), wts):
        ik = getattr(row, "inchikey", None)
        if ik not in reg.ik2i and hasattr(row, "smiles"):
            m = Chem.MolFromSmiles(row.smiles)
            ik = Chem.MolToInchiKey(m) if m is not None else None
        i = reg.ik2i.get(ik)
        if i is None or cls[i] not in CLASSES:
            continue
        r_cls[CLASSES.index(cls[i])] += wv
        r_nc[int(np.clip(nc[i], NC_MIN, NC_MAX)) - NC_MIN] += wv
        matched += wv
    print(f"[reference] matched {matched / wts.sum():.1%} of reference mass to registry compounds")
    if matched <= 0:
        raise ValueError("No reference compounds matched compounds.csv")
    return r_cls / matched, r_nc / matched


def reference_distance(w, cls_oh, nc_oh, r_cls, r_nc):
    """100 * sum_c (class mass - ref)^2  +  W1 distance of carbon-number distributions (carbon-number units)."""
    d_cls = ((w @ cls_oh - r_cls) ** 2).sum(-1)
    w1 = (torch.cumsum(w @ nc_oh, -1) - torch.cumsum(r_nc, -1)).abs().sum(-1)
    return 100.0 * d_cls + w1
