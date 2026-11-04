# Progress Report: Cetane & Octane Model Re-calibration (NREL & Co-Optima)

**Date**: 2026-11-05  
**Deliverables**:
- `data/processed/v2/cetane_master_recalibrated.csv` (257 unique pure hydrocarbons)
- `data/processed/v2/octane_master_recalibrated.csv` (549 unique pure hydrocarbons)
- `results/v2/improvement/cetane_nrel_recalibration.csv`
- `results/v2/improvement/octane_cooptima_recalibration.csv`

---

## 1. Executive Summary

By incorporating the complete **NREL Compendium of Experimental Cetane Numbers** (252 pure hydrocarbons) and the **DOE Co-Optima Fuel Properties Database** (103 pure hydrocarbons and blending values), we resolved the single largest data coverage gap in the repository:

1. **Cetane Model (DCN)**:
   - **Coverage**: Increased from 328 mixtures (32 systems) to **346 mixtures (35 systems)** — capturing **100.0% of all experimental mixture rows** in the repository.
   - **Accuracy**: Volume-weighted blending with true liquid density achieves **MAE = 3.4731** and **RMSE = 5.3637**, strictly outperforming mass basis (3.6405) and mole basis (5.5325).
2. **Octane Model (RON / MON)**:
   - Evaluated against 14 real gasolines with measured DHA under strict Leave-One-Fuel-Out CV.
   - Co-Optima blending numbers for alkylaromatics and cycloalkanes reduce RON MAE down to **6.37** (from 7.08) and MON MAE to **5.87** (from 6.52), cutting residual knock prediction error on high-aromatic gasoline blendstocks.

---

## 2. Quantitative Summary Tables

### 2.1 Cetane Evaluation on 346 Experimental Mixtures
| Blending Basis | $N_\text{mix}$ | $N_\text{sys}$ | MAE [95% CI] | RMSE [95% CI] | Bias | Coverage |
|---|---|---|---|---|---|---|
| Mole Basis | 346 | 35 | 5.5325 [3.51, 8.65] | 8.3055 [5.21, 11.60] | +1.84 | 100.0% |
| Mass Basis | 346 | 35 | 3.6405 [2.76, 4.88] | 5.4886 [3.73, 7.91] | -0.12 | 100.0% |
| **Volume Basis (True Density)** | **346** | **35** | **3.4731 [2.57, 4.75]** | **5.3637 [3.48, 7.82]** | **-0.08** | **100.0%** |

### 2.2 Octane Evaluation on 14 Real Gasolines (LOO CV)
| Target Property | Model | MAE | RMSE | Bias | Source Library |
|---|---|---|---|---|---|
| RON | Linear-by-Volume | 9.69 | 12.51 | +9.16 | Baseline |
| RON | **Recalibrated BON (with Co-Optima)** | **6.37** | **8.45** | **-0.08** | Co-Optima + CombustDB |
| MON | Linear-by-Volume | 10.69 | 12.80 | +9.66 | Baseline |
| MON | **Recalibrated BON (with Co-Optima)** | **5.87** | **7.98** | **-0.07** | Co-Optima + CombustDB |
