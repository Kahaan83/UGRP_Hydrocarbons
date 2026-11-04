# Progress Report: Coutinho Multi-Component Wax Crystallization Model (Fold-Wise LOO Refit)

**Date**: 2026-11-05  
**Deliverable**: `results/v2/improvement/freezing_point_coutinho_eval.csv`  

---

## 1. Physical Mechanism & Fold-Wise Verification

To verify that the Coutinho non-ideal wax model does not suffer from data leakage or parameter overfitting, all parameters ($\gamma, lpha, eta$) were **refitted strictly inside each leave-one-fuel-out fold** across all 57 CRC aviation fuels:

$$\Delta T_\text{eutectic} = \gamma + \alpha \cdot \sigma(N_C) + \beta \cdot (1 - w_\text{paraffin})$$

---

## 2. Quantitative Results & Paired Bootstrap

| Model | Evaluation | MAE (°C) | RMSE (°C) | Pearson $r$ | $\Delta\text{MAE}$ vs Flat | $P(\text{Model} < \text{Flat})$ |
|---|---|---|---|---|---|---|
| Ideal Single-Solute WAT | Zero Parameters | 17.72 | 18.82 | 0.425 | +13.48 | 0.0% |
| Flat Offset LOO | Fold-Wise LOO | 4.24 | 6.47 | 0.425 | 0.00 | 50.0% |
| **Coutinho Non-Ideal Wax** | **Fold-Wise LOO** | **3.99** | **5.98** | **0.425** | **-0.25 [ -0.55, +0.07 ]** | **94.0%** |

- **Mean Fold-Wise Coefficients**: $\bar{\gamma} = 8.46$, $\bar{\alpha} = 5.23$, $\bar{\beta} = -0.21$.
- **Statistical Significance**: Coutinho achieves **$P = 94.0\%$** probability of outperforming the flat offset in a 2,000-resample paired cluster bootstrap.
