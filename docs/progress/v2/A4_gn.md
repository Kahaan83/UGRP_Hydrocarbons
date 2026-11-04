# Task A4 Progress Report: Grunberg-Nissan Interaction Baseline Audit

**Task ID**: A4  
**Date**: 2026-11-04  
**Deliverable**: `results/v2/A4/gn_baselines.csv`  
**Test Suite**: `tests/test_v2_a4_gn.py`  

---

## 1. Executive Summary & Root Cause Analysis

In the v1 report (`FULL_REPORT.md` Section 3.4 and `results/t5/literature_baselines_eval.csv`), Grunberg-Nissan was quoted with a log-RMSE of **0.0699**, virtually tied with Refutas at **0.0698**.

Our audit reveals two critical findings:
1. **The v1 Table Evaluated $G = 0$ (Zero Interaction)**: `literature_baselines.py` and `baselines_eval.py` explicitly hardcoded `G_interaction = 0.0`. In fact, `build_harder_splits_and_eval.py` simply hardcoded the value `0.0699` across both models.
2. **Actual Ideal GN Error is 0.1014**: When rigorously evaluated on the 15,507 mixture observations across 321 systems, ideal log-linear mixing ($G=0$) achieves a log-RMSE of **0.1014** (MAPE = 7.42%), which is **45.2% worse than Refutas** (0.0699)!
3. **Fitted Interaction Schemes**:
   - **Global $G$ (5-Fold Grouped CV)**: Fitting a single global $G$ yields $G = -0.1494$, reducing log-RMSE slightly to **0.0976**.
   - **Class-Pair $G_{c_1, c_2}$ (5-Fold Grouped CV)**: Fitting 15 hydrocarbon class-pair parameters reduces log-RMSE to **0.0890**.
   - **Refutas Blending Index**: Without fitting *any* mixture interaction parameters, Refutas achieves **0.0699** log-RMSE, outperforming even the 15-parameter class-pair Grunberg-Nissan model.
   - **Per-System $G$ (In-Sample Bound)**: Fitting an individual interaction parameter for each of the 321 systems achieves **0.0412** log-RMSE, showing the theoretical upper bound of single-parameter pairwise interaction modeling.

---

## 2. Quantitative Results

| Model | $N_\text{obs}$ | $N_\text{sys}$ | Fold Discipline | Log-RMSE | Log-MAE | MAPE (%) | Parameters |
|---|---|---|---|---|---|---|---|
| Grunberg-Nissan ($G=0$, Ideal) | 15,507 | 321 | Zero Parameters | 0.1014 | 0.0669 | 7.42% | 0 |
| Grunberg-Nissan (Global $G$) | 15,507 | 321 | 5-Fold Grouped CV | 0.0976 | 0.0589 | 6.55% | 1 ($G = -0.1494$) |
| Grunberg-Nissan (Class-Pair $G$) | 15,507 | 321 | 5-Fold Grouped CV | 0.0890 | 0.0555 | 6.13% | 15 class pairs |
| **Refutas Blending Index** | 15,507 | 321 | Zero Parameters | **0.0699** | **0.0526** | **5.43%** | **0 (Pure Physics)** |
| Grunberg-Nissan (Per-System $G$) | 15,507 | 321 | In-Sample Bound | 0.0412 | 0.0232 | 2.50% | 321 |

---

## 3. Scientific Significance

1. **Why Refutas Beats Uncalibrated / Class-Calibrated Grunberg-Nissan**: Viscosity blending exhibits a profound double-logarithmic curvature ($\ln \ln(
u + 0.8)$). Pure log-linear mixing ($\ln \eta$) fails to capture this asymmetric curvature. Refutas encodes this thermodynamic curvature purely through mass-weighted VBN blending without requiring empirical binary interaction parameters.
2. **Claims Register Impact**: Claim CLM-06 in `docs/claims_register.md` is updated: the v1 quote of GN 0.0699 was an artifact of copying the Refutas baseline number; unparameterized GN is 0.1014, and Refutas is strictly superior to generalized class-pair Grunberg-Nissan.
