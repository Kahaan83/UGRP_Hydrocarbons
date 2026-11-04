# Progress Report: Coutinho Multi-Component Wax Crystallization Model

**Date**: 2026-11-05  
**Deliverable**: `results/v2/improvement/freezing_point_coutinho_eval.csv`  

---

## 1. Physical Mechanism

The v2 flat offset model (T_freeze = WAT - 14.58 °C) accurately corrected the mean bias, but applied an identical shift to all fuels regardless of their chemical composition.

In reality, real fuels exhibit varying degrees of **eutectic depression** depending on:
1. **n-Alkane Chain Length Polydispersity** (sigma(N_C)): When multiple n-alkanes of neighboring carbon numbers co-precipitate, they form an orthorhombic solid solution that stabilizes the solid phase at a depressed temperature relative to the highest-melting pure compound.
2. **Non-Paraffin Solvent Dilution**: High concentrations of aromatics and cycloalkanes depress the wax appearance temperature by lowering the liquid-phase thermodynamic activity gamma_i^L.

The upgraded Coutinho model in `src/fuelblend/thermo/coutinho_wax.py` computes:
Delta T_eutectic = gamma_0 + alpha * sigma(N_C) + beta * (1 - w_alkanes)

---

## 2. Evaluation on 57 CRC Aviation Real Fuels

| Model | Evaluation | MAE (°C) | RMSE (°C) | Bias (°C) | Pearson r |
|---|---|---|---|---|---|
| Ideal Single-Solute WAT | Zero Parameters | 17.72 | 18.25 | +17.72 | 0.425 |
| v2 Flat LOO Offset | LOO CV | 4.99 | 6.54 | +3.14 | 0.425 |
| **Coutinho Multi-Solid Wax Model** | **LOO CV** | **4.04** | **6.07** | **+2.32** | **0.402** |
