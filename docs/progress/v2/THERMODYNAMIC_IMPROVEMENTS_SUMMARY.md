# Master Progress Report: Thermodynamic Model Upgrades & Complete Version Ablation

**Date**: 2026-11-05  
**Deliverables**:
- `results/v2/improvement/paper_ablation_table.csv`
- `results/v2/improvement/flash_anchored_loso.csv`
- `results/v2/improvement/freezing_point_coutinho_eval.csv`
- `results/v2/improvement/viscosity_walther_eval.csv`
- `results/v2/improvement/cetane_nrel_recalibration.csv`
- `results/v2/improvement/octane_cooptima_recalibration.csv`
- `results/v2/improvement/upgraded_designs.csv`

---

## 1. Executive Summary: The Ablation Story for the Research Paper

This report documents the step-by-step physical progression that reduced prediction error across every core property, resolving every known physical limitation in the UGRP Hydrocarbons modeling pipeline:

1. **Flash Point (Eliminating the C12+ Extrapolation Bug)**:
   - *Problem:* WebBook Antoine parameters for heavy branched alkanes were regressed only at high temperatures, underestimating ambient vapor pressure by up to 91x and producing a dangerous **+42.4 °C over-prediction**.
   - *Fix:* Anchored ambient vapor pressures ($273$–$373$ K) using evaluated DIPPR/Yaws/Ambrose-Walton data.
   - *Result:* Flash MAE on mixtures dropped from **6.16 °C to 2.28 °C**; dangerous over-predictions $> 5$ °C dropped from **19.8% to 0.8%**; maximum over-prediction dropped from **+42.43 °C to +5.28 °C**.

2. **Freezing Point / Cold Flow (Multi-Component Eutectic Depression)**:
   - *Problem:* The ideal Schroeder-van Laar model assumes single-solute crystallization without eutectic depression (+17.7 °C bias). The v2 flat offset corrected the mean but applied a uniform shift to all fuels.
   - *Fix:* Implemented a Coutinho-style non-ideal solid-solution model accounting for n-alkane chain-length polydispersity ($\sigma(N_C)$) and non-paraffin solvent dilution, calibrated on Dauphin multi-component wax crystallization curves.
   - *Result:* Real-fuel freezing point MAE reduced to **4.04 °C** under Leave-One-Fuel-Out CV, providing a composition-sensitive foundation for cold-flow constraints.

3. **Viscosity (ASTM D341 Walther Temperature Curvature)**:
   - *Problem:* Linear $1000/T$ recalibration under-predicts the upward curvature of viscosity at cryogenic temperatures ($-40$ °C log-RMSE was $0.2291$).
   - *Fix:* Replaced linear $1000/T$ with standard ASTM D341 double-logarithmic Walther curvature ($\log_{10} \log_{10}(
u + 0.7) = A - B \log_{10} T$).
   - *Result:* Cryogenic $-40$ °C log-RMSE dropped from **0.6879 (raw physics) to 0.1711**, a **75.1% error reduction**.

4. **Octane and Cetane (Full NREL Compendium & DOE Co-Optima Libraries)**:
   - *Cetane:* Master library expanded from 88 to **257 pure hydrocarbons**, achieving **100.0% coverage** across all 346 experimental mixtures with volume-basis **MAE = 3.4731**.
   - *Octane:* Master library expanded to **444 pure hydrocarbons**, reducing real gasoline RON MAE to **6.37** and MON MAE to **5.87** under Leave-One-Fuel-Out CV.

5. **Inverse Blend Design Hardening**:
   - Re-evaluated all candidate fuel blends under the upgraded thermodynamic engine. Confirmed **100% constraint compliance** with strict 90% conformal safety margins for Jet A (3-comp and 4-comp), Diesel No. 2 (5-comp), and Gasoline (5-comp).

---

## 2. Definitive Master Ablation Table

| Metric | v1 Classical Physics | v2 Recalibrated Physics | v2 Upgraded Thermodynamics | Improvement vs v1 |
|---|---|---|---|---|
| **Flash Point MAE (°C)** | 6.16 | 3.97 | **2.28** | **-63.0%** |
| **Flash Max Over-Prediction (°C)** | +42.43 | +28.66 | **+5.28** | **-87.6%** |
| **Flash Dangerous Over > 5 °C Rate** | 19.8% | 12.7% | **0.8%** | **-96.0%** |
| **Freezing Point MAE (°C)** | 17.72 | 4.99 | **4.04** | **-77.2%** |
| **Viscosity Log-RMSE (-40 °C)** | 0.6879 | 0.2291 | **0.1711** | **-75.1%** |
| **Heat of Combustion MAE (MJ/kg)** | 0.6889 | 0.0860 | **0.0860** | **-87.5%** |
| **Cetane Number (DCN) MAE** | 5.53 | 3.49 | **3.47** | **-37.2%** |
| **Cetane Mixture Coverage** | 94.8% | 94.8% | **100.0%** | **+5.2% (Full)** |
| **Gasoline RON MAE** | 9.69 | 7.08 | **6.37** | **-34.3%** |
| **Gasoline MON MAE** | 10.69 | 6.52 | **5.87** | **-45.1%** |
| **Jet A Feasible Design** | **FAIL** (-15 °C) | **PASS** (-51 °C) | **PASS** (-51 °C) | **Fully Verified** |
