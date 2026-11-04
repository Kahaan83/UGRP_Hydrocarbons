# Task A6 Progress Report: Uncertainty Growth Factor Provenance and Refit

**Task ID**: A6  
**Date**: 2026-11-04  
**Deliverable**: `results/v2/A6/uncertainty_growth_refit.csv`  
**Test Suite**: `tests/test_v2_a6_uq.py`  

---

## 1. Executive Summary & Provenance Audit

In `FULL_REPORT.md` Section 3.4 and `docs/progress/Phase4_many_components.md`, property-specific uncertainty growth factors were introduced:
- $g_{\text{visc}} = 1.377$ (claimed real fuel 90% coverage = 92.5%)
- $g_{\text{dens}} = 2.530$ (claimed real fuel 90% coverage = 91.8%)
- $g_{\text{flash}} = 1.541$ (claimed real fuel 90% coverage = 94.2%)

Our provenance audit reveals:
1. **The Numbers Were Carried Over From Synthetic Multicomponent Errors**:
   In `scripts/run_many_component_studies.py` line 174, the factors were hardcoded as the exact ratio of multicomponent ($N \ge 3$) test error to binary ($N=2$) test error on synthetic mixtures:
   $$g_{\text{visc}} = \frac0.06030.0438 = 1.3767 \approx 1.377$$
   $$g_{\text{dens}} = \frac6.0342.385 = 2.52998 \approx 2.530$$
   $$g_{\text{flash}} = \frac2.451.59 = 1.54088 \approx 1.541$$
   **None of these factors were fit on real fuels.**
2. **Catastrophic Failure of Real-Fuel Coverage**:
   When the claimed factors ($g_{\text{visc}}=1.377, g_{\text{dens}}=2.530$) are evaluated against the DeepSets neural ensemble standard deviations (`extrap_deepsets_ens_sd`) on the 57 CRC real aviation fuels:
   - For Kinematic Viscosity: Claimed $g=1.377$ yields an empirical coverage of **only 6.43%** (NOT 92.5%)! Over 93% of real-fuel measurements fall outside the 90% confidence envelope!
   - For Density: Claimed $g=2.530$ yields an empirical coverage of **only 15.79%** (NOT 91.8%)!
3. **Root Cause**:
   The 5-seed DeepSets ensemble standard deviation (`extrap_deepsets_ens_sd`) reflects only parameter initialization variance on clean, low-noise laboratory binary mixtures ($\sigma \approx 0.04$ cSt for viscosity, $\approx 0.6$ kg/m³ for density). On real fuels, domain shift, distillation reconstruction error, and low-temperature non-ideality produce actual errors of $\sim 1.74$ cSt and $\sim 7.5$ kg/m³. The neural ensemble is **severely overconfident by 10x to 50x**.
4. **Refit Factors Under Strict Leave-One-Fuel-Out CV**:
   - For DeepSets Viscosity: achieving true 90% coverage on held-out fuels requires an inflation factor of **$g = 49.88$**.
   - For DeepSets Density: achieving true 90% coverage on held-out fuels requires an inflation factor of **$g = 14.49$**.
   - For Recalibrated Refutas Physics (Task A2): because the physical model captures the proper $1000/T$ slope, a modest inflation factor of **$g = 1.94$** achieves true 90% coverage on held-out fuels.

---

## 2. Quantitative Evaluation Table

| Property | Predictor | Dataset | Claimed $g$ | Claimed Coverage | Actual Real-Fuel Coverage | Refit $g$ (LOO Mean) | Held-Out LOO Coverage |
|---|---|---|---|---|---|---|---|
| Viscosity | DeepSets Residual | CRC 647 (171 pts) | 1.377 | 92.5% | **6.43%** | **49.88** | 90.06% |
| Density | DeepSets Residual | CRC 647 (57 fuels) | 2.530 | 91.8% | **15.79%** | **14.49** | 89.47% |
| Flash Point | HGB Residual (LOSO) | Leenhouts 126 Mix | 1.541 | 94.2% | **82.54%** | **2.72** | 90.48% |
| **Viscosity (Recal Phys)** | **Recalibrated Refutas** | **CRC 647 (171 pts)** | — | — | — | **1.94** | **89.47%** |

---

## 3. Claims Register Impact

Claim CLM-08 in `docs/claims_register.md` is recorded:
The v1 claims that $g_{\text{visc}}=1.377$ and $g_{\text{dens}}=2.530$ provided $>90\%$ empirical coverage on real fuels are **RETRACTED**. The factors were unvalidated synthetic carry-overs that produced severe under-coverage ($6.4\%$ and $15.8\%$). Real-fuel neural predictions require $g \approx 50$ due to severe ensemble overconfidence, whereas recalibrated physics is well-behaved ($g = 1.94$).
