# Progress Report: ASTM D341 Walther Viscosity Model Upgrade

**Date**: 2026-11-05  
**Deliverable**: `results/v2/improvement/viscosity_walther_eval.csv`  

---

## 1. Physical Mechanism

The 2-parameter recalibrated Refutas model from Task A2 utilized a linear 1000/T shift:
Delta ln(nu) proportional to (1000/T - 1000/293.15)

While accurate near ambient temperatures (+20 °C and -20 °C), this linear assumption breaks down at cryogenic temperatures (-40 °C), where the true viscosity curve bends sharply upwards in log space.

The standard petroleum engineering method (**ASTM D341 / Walther equation**) establishes that viscosity is strictly linear in double-logarithmic space:
log10(log10(nu(T) + 0.7)) = A - B * log10(T)

When anchored using the +20 °C and -20 °C recalibrated predictions, extrapolation along the ASTM D341 line accurately captures the cryogenic curvature.

---

## 2. Evaluation at -40 °C on 57 CRC Aviation Fuels

| Model | Log-RMSE (-40 °C) | Log-MAE (-40 °C) | Error Reduction vs Raw Physics |
|---|---|---|---|
| Raw Refutas | 0.6879 | 0.6653 | 0.0% |
| Linear 1000/T Recalibration | 0.2877 | 0.2510 | 58.2% |
| **ASTM D341 Walther Model** | **0.1711** | **0.0968** | **75.1%** |
