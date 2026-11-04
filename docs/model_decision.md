# Master Model Strategy Decision Note

**Project:** UGRP Hydrocarbons v2  
**Status:** Canonical Architecture Decision  
**Scope:** Model Selection Across Properties and Composition Regimes  

---

## 1. Unequivocal Model Architecture Decision

Based on empirical evidence across Phase A (triage & LOO recalibration), Phase C (property baseline fixes), and Phase D (stacked evaluation, regime analysis, truncation sensitivity), the canonical modeling strategy is decided as follows:

| Property | Primary Production Model | Secondary / Lab Interpolator | Empirical Evidence Table |
|---|---|---|---|
| **Kinematic Viscosity** ($\nu$) | **2-Parameter Recalibrated Refutas** ($1000/T$ Slope) | Stacked Refutas + Shrunk Residual | `results/v2/A2/recalibration_loo.csv`<br>`results/v2/D4/stacked_eval.csv` |
| **Density** ($\rho$) | **Ideal Volume Blending** ($\rho_{ideal} + \Delta\rho$) | Stacked Ideal Volume + Shrunk Ridge | `results/v2/D4/stacked_eval.csv` |
| **Flash Point** ($T_{flash}$) | **Transparent Offset Correlation / Le Chatelier** | Le Chatelier + Shrunk HGB Residual (LOSO) | `results/v2/C2/flash_loso_v2.csv`<br>`results/v2/D4/stacked_eval.csv` |
| **Freezing Point** ($T_{freeze}$) | **Transparent Offset SLE WAT** ($WAT - 14.58$ °C) | Linear Recalibrated SLE | `results/v2/C4/cold_flow_v2.csv`<br>`results/v2/D4/stacked_eval.csv` |
| **Heat of Combustion** ($NCV$) | **Group Additivity GCV Baseline (v2)** | Stacked GCV + Shrunk Ridge | `results/v2/C1/hoc_eval.csv`<br>`results/v2/D4/stacked_eval.csv` |

---

## 2. Evidence Across Regimes: Ambient Lab Binaries vs Cryogenic Real Fuels

The empirical data across blend complexity buckets $N \in [1, 40+]$ (`results/v2/D5/error_vs_n.csv`) proves an **architectural inversion** between regimes:

1. **Ambient Laboratory Mixtures ($N = 2, 3$; $T \in [293, 323]$ K)**:
   - Learned residual models succeed: DeepSets achieves log-MAE $0.0372$ ($[0.0289, 0.0469]$) and Pairwise achieves log-MAE $0.0342$ ($[0.0266, 0.0432]$), improving over raw Refutas ($0.0616$) and Grunberg-Nissan ($0.0669$, `results/v2/A4/gn_baselines.csv`).
   - *Scope for learning*: Neural and tree-based residuals are strictly justified here to interpolate excess thermodynamic mixing friction and activity coefficients.

2. **Cryogenic Real Multi-Component Fuels ($N \ge 20$; $T \in [-40, +20]$ °C)**:
   - **DeepSets and Pairwise neural residuals fail catastrophically**: DeepSets extrapolates to log-MAE $0.3159$ ($[0.2905, 0.3411]$) and Pairwise to $0.3661$ ($[0.3406, 0.3913]$), worse than raw physics at moderate temperatures (`results/v2/A2/recalibration_loo.csv`).
   - **Empirical uncertainty collapse**: In `results/v2/A6/uncertainty_growth_refit.csv`, DeepSets claimed $92.5\%$ coverage on real fuels, but collapsed to an actual coverage of **$6.43\%$** (viscosity) and **$15.79\%$** (density).
   - **Recalibrated Physics Dominates**: The 2-parameter physical recalibration ($\ln \nu = \ln \nu_{Refutas} + a_0 + a_1 \frac{1000}{T}$) achieves log-MAE **$0.1042$** ($[0.0890, 0.1205]$), reducing error by **$69.6\%$** over raw Refutas and outperforming DeepSets by more than $3\times$.

---

## 3. Strict Boundary Rules for Downstream Tasks (Phase E & F)

1. **Physics Dominance Everywhere in Inverse Design (Phase E)**:
   - All property evaluations for jet, diesel, and gasoline inverse design must evaluate primary physics baselines first.
   - For viscosity, use the 2-parameter recalibrated Refutas slope.
   - For freezing point, use the transparent offset SLE WAT baseline ($T_{freeze} = WAT - 14.58$ °C).
   - For net heat of combustion, use the group additivity GCV baseline (`results/v2/C1/hoc_eval.csv`, MAE $0.179$ MJ/kg).

2. **Shrinkage Prior on Learned Residuals**:
   - Where a learned residual is employed (the stacked model architecture in `results/v2/D4/stacked_eval.csv`), it must carry an explicit $L_2$ shrinkage prior ($\alpha \ge 50.0$) ensuring that in data-sparse or out-of-distribution regimes, the residual collapses to zero rather than extrapolating arbitrarily.

3. **Surrogate Cardinality Rule (from Task D6)**:
   - Inverse design blends must not use over-simplistic 2- or 3-component surrogates (e.g., dodecane + tetralin), which fail freezing point and boiling range constraints.
   - Per `results/v2/D6/truncation.csv`, a surrogate complexity of **$K = 10$ to $20$ components with a lumped heavy tail** captures $>95\%$ of blend mass and keeps bulk property drift below $0.3\%$.
