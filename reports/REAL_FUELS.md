# REAL_FUELS.md - WP7 Real-Fuel Surrogate Bridge & Bulk Distillation Model

## 1. Overview
Real aviation and diesel fuels (e.g. Jet A, JP-8, RP-1, CRC-647/AV-27-18 fuel sets, and 2,600+ scraped AJFTD fuels) lack explicit InChIKey/SMILES molecular lists in experimental databases. Instead, they are characterized by **bulk inspection metrics**:
- ASTM D86 Distillation curves ($T_{0}$, $T_{5}$, $T_{10}$, $T_{50}$, $T_{90}$, $T_{100}$)
- Bulk density ($ho_{15.56\ ^\circ\text{C}}$)
- Hydrocarbon class fractions ($w_{\text{aromatics}}$, $w_{\text{saturates}}$, $w_{\text{olefins}}$)

WP7 establishes the predictive link between bulk inspection features and core combustion/flow safety properties.

---

## 2. Real Fuel Flash Point Prediction from Distillation (D86)

In petroleum engineering, Nelson's empirical correlation estimates flash point from the 10 % ASTM D86 distillation temperature:
$$T_{\text{flash, Nelson}} (^\circ\text{C}) = 0.69 \cdot T_{10} (^\circ\text{C}) - 68.0$$

We evaluate this classical engineering baseline against a 5-fold cross-validated ridge model using bulk distillation and composition features on $N = 1959$ real aviation fuels from CRC 647 and the full AJFTD database.

### Results on Real Fuels ($N=1959$)

| Model | Input Features | CV MAE (°C) | CV RMSE (°C) |
|---|---|---|---|
| **Nelson Classical Baseline** | $T_{10}$ only | 8.26 °C | 14.34 °C |
| **Bulk Inspection Ridge Model** | $T_{10}, T_{50}, T_{90}$, Density, Aromatics % | **3.36 °C** | **6.13 °C** |

The multi-metric bulk model achieves a **59.3 % MAE reduction** over the standard Nelson refinery correlation across 1959 real fuel records.

---

## 3. Status of Real Fuel Datasets
- Scraped AJFTD dataset now covers **2,679 unique fuels** across 91,175 metric rows.
- Validated flash point records: 2,413; Kinematic viscosity records: 2,481; Density records: 2,670; Freezing point records: 2,672.
- All uncharacterized blends are routed to `data/processed/excluded_by_R1.csv` for pure/blend training safety, while available as empirical benchmarks for bulk surrogate models.
