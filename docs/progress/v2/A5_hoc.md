# Task A5 Progress Report: Heat of Combustion Reinterpretation & Out-of-Sample Verification

**Task ID**: A5  
**Date**: 2026-11-04  
**Deliverable**: `results/v2/A5/hoc_recheck.csv`  
**Test Suite**: `tests/test_v2_a5_hoc.py`  

---

## 1. Executive Summary & Root Cause Analysis

In `FULL_REPORT.md` Section 3.4 and `docs/progress/T2.4.md`, the Ridge residual model was reported to achieve:
- Physics Baseline: MAE = 0.689 MJ/kg (MAPE = 1.60%)
- Ridge Residual: MAE = 0.078 MJ/kg (MAPE = 0.18%)

However, the report bizarrely dismissed this dramatic improvement with the statement:
> *"the physics baseline is already remarkably accurate... and the learned residual adds virtually no improvement (reducing RMSE by only 0.6020 MJ/kg)."*

Our verification audit resolves three key questions:
1. **The In-Sample Evaluation Flaw**: In v1 `scripts/eval_heat_of_combustion_baseline.py`, Ridge was fit on all 57 fuels and predicted on the exact same 57 fuels without cross-validation (producing an in-sample bias of $-1.99 \times 10^{-15}$).
2. **True Out-of-Sample Generalization (Leave-One-Fuel-Out CV)**: Under strict leave-one-fuel-out cross-validation across all 57 CRC aviation fuels:
   - Physics Baseline: MAE = **0.6889 MJ/kg** [0.6466, 0.7335], Bias = **+0.6889 MJ/kg**
   - Ridge (Strict LOO CV): MAE = **0.0860 MJ/kg** [0.0664, 0.1065], Bias = **-0.0026 MJ/kg**
   - **Error Reduction**: **87.5% reduction in MAE** (from 0.6889 to 0.0860 MJ/kg).
3. **Scientific Mechanism**:
   - The physics baseline exhibits a persistent +0.689 MJ/kg positive offset across all fuels due to pure-component GCV estimation (Dulong formula) or latent heat of water subtraction.
   - The Ridge model utilizes two bulk physical features: `density_15C_kg_m3` and `aromatics_vol_pct`.
   - In aviation fuels, Net Heat of Combustion (NCV) is directly determined by hydrogen content ($H$ wt%). Higher density and aromaticity correspond directly to lower hydrogen mass fraction. Standard engineering practice (ASTM D3338 and ASTM D4529) computes heating value from density and aromatics. Ridge on `[density, aromatics]` effectively re-estimates this empirical ASTM calibration, eliminating the +0.689 MJ/kg physics bias and achieving accuracy within the ASTM D4809 experimental reproducibility (~0.20 MJ/kg).

---

## 2. Evaluation Table

| Model | Evaluation | $N_\text{fuels}$ | Features | MAE (MJ/kg) [95% CI] | RMSE (MJ/kg) [95% CI] | Bias (MJ/kg) | MAPE (%) |
|---|---|---|---|---|---|---|---|
| Physics Baseline | Zero Parameters | 57 | Pure GCV + $w_i$ | 0.6889 [0.6466, 0.7335] | 0.7076 [0.6588, 0.7571] | +0.6889 | 1.60% |
| Physics + Ridge (v1 In-Sample) | In-Sample Fit | 57 | Density, Aromatics | 0.0783 [0.0608, 0.0972] | 0.1056 [0.0823, 0.1286] | 0.0000 | 0.18% |
| **Physics + Ridge (Strict LOO CV)** | **Leave-One-Out** | 57 | Density, Aromatics | **0.0860 [0.0664, 0.1065]** | **0.1219 [0.0945, 0.1504]** | **-0.0026** | **0.20%** |

---

## 3. Claims Register Impact

Claim CLM-07 in `docs/claims_register.md` is updated:
The v1 conclusion that *"learned residual adds virtually no improvement"* is **RETRACTED**. The Ridge residual achieves a massive **87.5% error reduction** out-of-sample under strict leave-one-fuel-out cross-validation, resolving the +0.689 MJ/kg physics baseline bias via ASTM D3338-style density-aromatics hydrogen content recalibration.
