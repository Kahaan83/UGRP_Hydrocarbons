# Progress Report: Flash Point Thermodynamic Anchoring Upgrade

**Date**: 2026-11-05  
**Deliverables**:
- `data/processed/v2/compounds_anchored_vp.csv`
- `results/v2/improvement/flash_anchored_loso.csv`

---

## 1. Physical Mechanism & Problem Solved

In previous versions, WebBook Antoine parameter sets for heavy alkanes (>= C12) were regressed solely at high temperatures (95–220 °C). When extrapolated down to flash point temperatures (20–60 °C):
- For isododecane (2,2,4,6,6-pentamethylheptane), WebBook extrapolated vapor pressure was 0.0029 kPa (actual is 0.2682 kPa, an underestimation factor of **91x**!).
- In Le Chatelier's equation, this caused an artificial **+42.4 °C dangerous over-prediction** on hydrocarbon mixtures containing heavy branched alkanes.

By anchoring vapor pressures across ambient temperatures (273.15–373.15 K) using evaluated DIPPR/Yaws/Ambrose-Walton data:
- Over-prediction rate > 5 °C on C12+ mixtures drops from **32.1% to 1.3%**!
- Maximum over-prediction drops from **+42.43 °C down to +5.28 °C**!
- Overall flash MAE on C12+ mixtures drops from **8.98 °C down to 2.84 °C**!

---

## 2. Quantitative Comparison Table

| Evaluation Slice | Old WebBook Base MAE | New Anchored Base MAE | Old Max Over-Prediction | New Max Over-Prediction | Dangerous Over >5 °C Rate |
|---|---|---|---|---|---|
| All Mixtures (N=126) | 6.16 °C | **2.28 °C** | +42.43 °C | **+5.28 °C** | 19.8% -> **0.8%** |
| C12+ Mixtures (N=78) | 8.98 °C | **2.84 °C** | +42.43 °C | **+5.28 °C** | 32.1% -> **1.3%** |
| <C12 Mixtures (N=48) | 1.56 °C | **1.37 °C** | +2.04 °C | **+1.88 °C** | 0.0% -> **0.0%** |
