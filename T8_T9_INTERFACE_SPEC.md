# Interface & Specification Document for Tasks T8 & T9

**Project:** Hydrocarbon Fuel Property Prediction & Inverse Blend Design (UGRP)  
**Date:** November 3, 2026  
**Target:** Implementation guide for Task T8 (DeepSets Architecture) & Task T9 (Inverse Blend Design)

---

## 1. Mixture Table Schema

Experimental mixture data is stored in **long format** (one row per experimental measurement condition). Each row contains the blend's components, composition fractions (mass and mole), thermodynamic state variables, and target property values.

### 1.1 Data Sources & File Locations

1. **Primary Model Table:** `data/processed/model_table.csv`
   - **42,662 rows** filtered strictly to all-hydrocarbon systems (`is_hc == True`).
   - Viscosity: 2,270 pure, 16,037 binary, 625 multicomponent rows.
   - Density: 12,747 binary, 897 multicomponent rows.
   - Flash point: 789 pure, 123 binary, 3 multicomponent rows (126 mixture rows total).
   - Contains split assignments from `splits/s1s2s3.json` mapped to `system_id`.

2. **Supplemental Leenhouts Mixture Flash Table (Task T7):** `data/processed/t7_leenhouts_flash_mix.csv`
   - **1,006 rows** from Leenhouts `flashpoint_mix_exp.csv` (681 binary + 325 multicomponent blends, 113 unique systems).
   - Column `is_hc`: **126 rows** are strictly all-hydrocarbon (identical systems as in `model_table.csv`), while **880 rows** contain non-hydrocarbon fuel components (e.g. oxygenates/alcohols/esters from the broader Leenhouts benchmark).
   - Component InChIKeys, mass fractions, mole fractions, and MD5 `system_id` hashes are fully canonicalized.

### 1.2 Column Layout & Types

| Column Name | Type | Description |
|---|---|---|
| `property` | `string` | Target property name: `viscosity_dyn_mPas`, `density_kg_m3`, `flash_point_C` |
| `system_id` | `string` | 12-char hex MD5 hash of `'\|'.join(sorted(comp_inchikeys))`. Grouping key for `GroupKFold` |
| `T_K` | `float` | Temperature in Kelvin (NaN if not reported or temperature-independent) |
| `P_kPa` | `float` | Pressure in kPa (typically 101.325 or NaN) |
| `kind` | `string` | System type: `'pure'`, `'binary'`, `'multi'` |
| `n_comp` | `int` | Number of distinct components ($\ge 1$) |
| `comp_inchikeys` | `string` (JSON) | JSON array of sorted InChIKeys, e.g. `["BGHCVCJVXZWKCC-UHFFFAOYSA-N", ...]` |
| `w_mass` | `string` (JSON) | JSON array of mass fractions aligned with `comp_inchikeys` (sums to 1.0) |
| `x_mole` | `string` (JSON) | JSON array of mole fractions aligned with `comp_inchikeys` (sums to 1.0) |
| `value` | `float` | Target experimental value: mPa·s for viscosity, kg/m³ for density, °C for flash point |
| `quality_weight` | `float` | Sample weight based on measurement reliability |
| `is_hc` | `bool` | True if all components are hydrocarbons (C, H only) |
| `source` / `sources` | `string` | Provenance database (e.g. `leenhouts`, `nist_tde`, `thermoml`) |
| `conflict` | `bool` | True if replicate observations diverge beyond tolerance |

### 1.3 Exact Sample Rows

#### Binary Flash Point (`data/processed/model_table.csv`)
```python
{
    'property': 'flash_point_C',
    'system_id': 'f70a95ecf8ed',
    'T_K': np.nan,
    'P_kPa': np.nan,
    'kind': 'binary',
    'n_comp': 2,
    'comp_inchikeys': '["DCAYPVUWAIABOU-UHFFFAOYSA-N", "LTEQMZWBSYACLV-UHFFFAOYSA-N"]',
    'w_mass': '[0.900183, 0.099817]',
    'x_mole': '[0.866, 0.134]',
    'value': 117.0,
    'quality_weight': 1.0,
    'is_hc': True,
    'source': 'leenhouts',
    'conflict': False
}
```

#### Multicomponent Mixture Flash (`data/processed/t7_leenhouts_flash_mix.csv`)
```python
{
    'property': 'flash_point_C',
    'system_id': 'f630521a0719',
    'T_K': 298.15,
    'P_kPa': 101.325,
    'kind': 'multi',
    'n_comp': 3,
    'comp_inchikeys': '["BGHCVCJVXZWKCC-UHFFFAOYSA-N", "DIOQZVSQGTUSAI-UHFFFAOYSA-N", "SNRUBQQJIBEYMU-UHFFFAOYSA-N"]',
    'w_mass': '[0.599342, 0.200167, 0.20049]',
    'x_mole': '[0.539, 0.251, 0.21]',
    'value': 71.85,
    'quality_weight': 1.0,
    'is_hc': True,
    'source': 'leenhouts',
    'conflict': False
}
```

### 1.4 Compound Registry (`data/unified/compounds.csv`)
Contains 1,038 pure compounds keyed by `inchikey`.
- **Key columns:** `inchikey`, `smiles`, `inchi`, `formula`, `mw`, `cls` (`aromatic`, `alkene_alkyne`, `cycloalkane`, `iso_paraffin`, `n_paraffin`, `other`), `n_carbon`, `antoine_A`, `antoine_B`, `antoine_C`, `antoine_Tlo`, `antoine_Thi`, `lfl_vol_pct`, `lfl_T_C`.
- **Antoine equation convention:** $\log_{10}(P / \text{kPa}) = A - \frac{B}{T / \text{K} + C}$.

---

## 2. Pure-Compound Features

Feature extraction is implemented in:
`fuelblend_data/src/fuelblend/features/pure_features.py`

### 2.1 Function Signatures

```python
import numpy as np
from fuelblend.features.pure_features import (
    pure_feature_row,
    pure_feature_frame,
    FEATURE_NAMES,
    N_FEATURES
)

# 1. Single compound feature row (returns 1D np.ndarray of shape (19,))
feat_1d = pure_feature_row(
    inchikey: str,
    smiles: str,
    T_K: float,
    cls: str,  # One of: 'aromatic', 'alkene_alkyne', 'cycloalkane', 'iso_paraffin', 'n_paraffin'
    classes: list[str] = ['aromatic', 'alkene_alkyne', 'cycloalkane', 'iso_paraffin', 'n_paraffin'],
    T_fallback_K: float = 298.15
) -> np.ndarray

# 2. Batch feature frame (returns 2D np.ndarray of shape (N, 19))
feat_2d = pure_feature_frame(
    rows: list[dict],  # list of dicts with keys: 'inchikey', 'smiles', 'T_K', 'cls'
    classes: list[str] = ['aromatic', 'alkene_alkyne', 'cycloalkane', 'iso_paraffin', 'n_paraffin'],
    T_fallback_K: float = 298.15
) -> np.ndarray
```

### 2.2 Feature Vector Layout (19 Features)

```python
FEATURE_NAMES = [
    # 1. RDKit Molecular Descriptors (11)
    'mol_wt', 'heavy_atom_count', 'n_aromatic_rings', 'ring_count',
    'formal_charge', 'hbond_donor', 'hbond_acceptor', 'rotatable_bonds',
    'labute_asa', 'tpsa', 'n_carbon',
    # 2. One-Hot Chemical Class Encoding (5)
    'cls_aromatic', 'cls_alkene_alkyne', 'cls_cycloalkane', 'cls_iso_paraffin', 'cls_n_paraffin',
    # 3. Temperature Transformations (3)
    'logT', 'T_K', 'inv_T_K'
]
```

---

## 3. `predict_blend.predict` API

Implemented in:
`fuelblend_data/src/fuelblend/api/predict_blend.py`

### 3.1 Function Signature & Input Contract

```python
from fuelblend.api.predict_blend import predict

result = predict(
    composition: list[tuple[str, float]],  # List of (SMILES or InChIKey, fraction)
    T_K: float = 298.15,                   # Temperature in Kelvin
    mole_fractions: bool = False,          # True if fractions are mole fractions; default False (mass fractions)
    T_C: float | None = None               # Optional: Temperature in Celsius (takes precedence over T_K if set)
) -> dict
```

### 3.2 Output Data Contract

```python
{
    'T_K': 320.0,
    'density_kg_m3': {
        'value': 720.68,                      # Predicted density in kg/m³
        'std': None,                          # Epistemic uncertainty (float or None)
        'method': 'residual_gbm+dens_idvol',  # Calculation path used
        'warnings': []                        # Property-specific warnings
    },
    'viscosity_dyn_mPas': {
        'value': 0.7615,                      # Dynamic viscosity in mPa·s
        'std': 0.00397,                       # Ensemble standard deviation
        'method': 'residual_gbm+refutas',     # Calculation path used
        'warnings': []
    },
    'flash_point_C': {
        'value': 55.60,                       # Flash point in °C
        'std': None,
        'interval': [41.75, 69.45],           # 90% split-conformal coverage interval [lo, hi]
        'method': 'chat',                     # 'chat' | 'gbm_mix' | 'none'
        'warnings': [
            'flash_data_thin: mixture flash training is sparse (123 binary + 3 multi rows)'
        ]
    },
    'warnings': []                            # Global domain warnings (e.g. extrapolation, T-bounds)
}
```

---

## 4. Execution Environment

- **Python Version:** `3.12.5` (Windows x64).
- **Virtual Environment Path:** `C:\Users\kahaa\OneDrive\Desktop\UGRP_hydrocarbons\.venv`
- **Installed Packages:**
  - `scipy` (includes `scipy.optimize`)
  - `scikit-learn`
  - `lightgbm`
  - `rdkit`
  - `pandas` (`3.0.6`)
  - `numpy` (`2.5.3`)
- **PyTorch & Optimization Packages:**
  - `torch` is not currently installed in `.venv`. **Installing PyTorch (`pip install torch`) is not an issue and can be done immediately if needed for T8.** Both CPU and CUDA / CS-lab compute can be utilized.
  - `cma` is not currently installed; it can be installed via `pip install cma` if CMA-ES is selected, or `scipy.optimize` can be used.

---

## 5. Scope & Requirements for Task T8 (DeepSets Architecture)

### 5.1 Objectives & Questions Specified by User

1. **Permutation Invariance:**
   - Architecture must be invariant to the ordering of components in the mixture:
     $$\mathbf{z}_{\text{mix}} = \sum_{i=1}^{n} w_i \cdot \phi(\mathbf{h}_i) \quad \text{or} \quad \rho\left(\sum_{i=1}^{n} \text{Attention}(\mathbf{h}_i, w_i)\right)$$
   - Component feature $\mathbf{h}_i$ is derived from molecular features (e.g., the 19D `pure_feature_row` vector or learned representations).

2. **Targets & Residual Formulation:**
   - User specification: Provide heads for the raw property and/or the residual over the physics baseline, with the physics baseline fed in as an input feature to ensure a fair comparison:
     - $\ln \eta_{\text{Refutas}}$ for viscosity
     - $\rho_{\text{idvol}}$ for density
     - $T_{\text{LeChatelier}}$ for flash point
   - Data coverage across properties is uneven (density: ~13k rows, viscosity: ~19k rows, flash point: ~126 all-hydrocarbon mixture rows, or up to 1,006 rows including oxygenates). Masked loss handles missing target properties per row.

3. **Validation & Benchmark:**
   - Evaluate against the fixed **S2 held-out split** (`splits/s1s2s3.json`) grouped strictly by `system_id`.
   - Direct benchmark comparison against the existing physics baselines (WP3) and HistGradientBoosting residual models (WP5).

---

## 6. Scope & Requirements for Task T9 (Inverse Blend Design)

### 6.1 Problem Formulation

Given target property specifications (e.g., matching Jet A / RP-1 targets):
- Target kinematic viscosity $\nu_{\text{target}}$ at $-20$ °C $\le 8.0$ cSt
- Target density $\rho_{\text{target}}$ at $15$ °C in $[775, 840]$ kg/m³
- Target flash point $T_{\text{flash}} \ge 38.0$ °C

Find component mass fractions $\mathbf{w} = [w_1, w_2, \dots, w_K]^T$ on the unit simplex:
$$\sum_{i=1}^{K} w_i = 1.0, \quad w_i \ge 0$$
minimizing the objective:
$$J(\mathbf{w}) = \sum_{p \in \text{props}} \lambda_p \cdot \mathcal{L}_p\big(\hat{y}_p(\mathbf{w}, T), y_{p,\text{target}}\big) + \lambda_{\text{unc}} \cdot \Omega_{\text{uncertainty}}(\mathbf{w}) + \lambda_{\text{domain}} \cdot \Omega_{\text{penalty}}(\mathbf{w})$$

### 6.2 Candidate Component Pools

The optimizer can be run over:

1. **Curated Surrogate Palette (10–15 Commercially Available Pure Hydrocarbons):**
   - *n-Paraffins:* n-octane, n-decane, n-dodecane, n-tetradecane
   - *iso-Paraffins:* 2,2,4-trimethylpentane (isooctane), 2,2,4,6,6-pentamethylheptane (isododecane)
   - *Cycloalkanes:* cyclohexane, methylcyclohexane, decalin (decahydronaphthalene), bicyclohexyl
   - *Aromatics:* toluene, ethylbenzene, m-xylene, mesitylene, 1-methylnaphthalene
   - *Property data status:* All have pure GBM models, verified Antoine constants, and LFL values in `compounds.csv`.

2. **Full Pure Compound Palette:**
   - All 170 pure hydrocarbon compounds in `data/unified/compounds.csv` that carry both Antoine vapor pressure constants and LFL values.

### 6.3 Sparsity & Practical Constraints

- Real fuel blends and certification surrogates typically utilize $\le 4$ or $\le 5$ active components.
- The optimizer should support enforcing or penalizing active component counts (e.g. via thresholding, $L_0$ penalty, or subset selection).
