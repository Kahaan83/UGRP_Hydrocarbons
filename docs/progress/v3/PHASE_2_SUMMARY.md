# Phase 2 Summary: Property Model Repairs across Transportation Fuels (Tasks G1–G11)

**Date**: 2026-11-05  
**Deliverables**:
- `results/v3/G2/octane_lump_model_eval.csv`
- `results/v3/G2/blending_octane_values.csv`
- `results/v3/G6/dvpe_activity_eval.csv`
- `configs/uncertainty_v3.yaml`
- `results/v3/G8/diesel_cetane_cold_flow_eval.csv`
- `results/v3/G10/thermoml_integration_eval.csv`
- `tests/test_v3_phase2_repairs.py`

---

## 1. Executive Summary

Phase 2 addresses the critical property modeling failures diagnosed in the baseline evaluation across the 153 expanded real fuels (Jet, Diesel, Gasoline):
1. **Gasoline Octane Number (G1–G4)**: Linear-by-mole blending failed with **RON MAE 6.37** and **MON MAE 5.87** due to extreme non-ideality in oxygenated blends (ethanol synergetic boost) and aromatic-olefin interactions. A **Ghosh-style 7-lump model** with L2 shrinkage regularizer was implemented and evaluated under Leave-One-Fuel-Out (LOO) CV, reducing error down to **RON MAE 3.10** and **Sensitivity MAE 0.98**, cutting error by >50%.
2. **Gasoline Vapor Pressure / DVPE (G5–G6)**: Ideal Raoult vapor pressure baseline exhibited a massive **18.50 kPa MAE** (and >50 kPa uncalibrated) due to positive ethanol-hydrocarbon azeotropy and omission of $C_4$–$C_5$ light ends. An **Activity-Corrected DVPE model** incorporating Margules activity coefficients ($\gamma_\\text{EtOH} \\approx 4.5$–$6.0$ at low dilution) and light-end speciation dropped error to **4.39 kPa MAE** under LOO CV.
3. **Conformal Flash Margin (G7)**: Updated `configs/uncertainty_v3.yaml` reflecting Task R2 lab mixture metrics (2.28 °C MAE, 5.28 °C worst-case) and Task R3 real-fuel segregation (+4.20 °C Pensky-Martens closed-cup compliance margin with 0.0% over-prediction on commercial aviation fuels).
4. **Diesel Cetane & Cold-Flow Hindcasts (G8–G9)**: Evaluated DCN under LOO CV on 16 real diesels (MAE 4.34 DCN) and Coutinho multi-solid wax precipitation model for cloud point (MAE 10.48 °C) and pour point (MAE 10.69 °C).
5. **Viscosity & Density Integration with Cleaned ThermoML (G10–G11)**: Integrated 44,575 verified liquid-phase ThermoML records (excess viscosities removed, phase strictly liquid). Density ideal volume blending achieved 12.94 kg/m³ MAE (~1.5% error) across all 153 fuels.

---

## 2. Octane Diagnosis & Ghosh Lump Model (G1–G4)

### Why Linear-by-Mole Blending Failed
Linear blending assumes zero excess chemical interaction:
$$ON_\\text{mix} = \\sum_i x_i ON_{0,i}$$
This fails drastically in real gasolines due to four mechanisms:
1. **Oxygenate Synergetic Boosting**: Pure ethanol has RON 108.6, but at 10–20 vol% in gasoline, its effective **Blending Octane Value (BOV)** reaches **125–160 ON** because ethanol scavenges low-temperature alkylperoxy radicals. Linear blending completely misses this +15 to +40 ON synergy.
2. **Aromatic-Olefin Antagonism**: Toluene and 1-hexene have antagonistic kinetic competition during intermediate-temperature autoignition, suppressing the expected octane response.
3. **Aromatic-Paraffin Synergy**: High-aromatic reformate blends synergistically with paraffinic alkylates.
4. **Heavy Ends Lumping**: Unresolved DHA tails (~2–5 wt%) contribute non-negligible knock resistance (~82 ON).

### Ghosh 7-Lump Formulation
The components are grouped into 7 lumps:
- n-paraffins ($RON_0 = 0, MON_0 = 0$)
- iso-paraffins ($RON_0 = 100, MON_0 = 100$)
- aromatics ($RON_0 = 120.1, MON_0 = 103.5$)
- olefins ($RON_0 = 76.4, MON_0 = 63.4$)
- naphthenes ($RON_0 = 80.0, MON_0 = 75.0$)
- ethanol ($RON_0 = 108.6, MON_0 = 89.7$)
- heavy tail ($RON_0 = 82.0, MON_0 = 78.0$)

The model couples pure lump values with regularized binary interaction terms:
$$\\Delta ON = k_{E,P} w_E w_P + k_{E,\\text{sat}} w_E (1 - w_E) + k_{A,O} w_A w_O + k_{A,P} w_A w_P + \\dots$$

### Evaluation Results (55 Gasolines, LOO CV)
| Model | RON MAE | RON RMSE | MON MAE | MON RMSE | Octane Sensitivity MAE |
|---|---:|---:|---:|---:|---:|
| Linear-by-Volume Baseline | 9.69 | 12.51 | 10.69 | 12.80 | 4.12 |
| Linear-by-Mole Baseline | 6.37 | 8.14 | 5.87 | 7.42 | 2.85 |
| **Ghosh-Style Lump Model (LOO)** | **3.10** | **6.59** | **2.63** | **5.89** | **0.98** |

---

## 3. DVPE Diagnosis & Activity Model (G5–G6)

### Why Ideal Raoult Failed (18.50 kPa Error)
1. **Positive Azeotropy of Ethanol**: Ethanol hydrogen-bonding clusters dissociate when diluted into non-polar hydrocarbons. The activity coefficient jumps to $\\gamma_\\text{EtOH} \\approx 4.5$–$6.0$, causing the partial vapor pressure $P_\\text{EtOH} = x_\\text{EtOH} \\gamma_\\text{EtOH} P_\\text{EtOH}^\\text{sat}$ to jump by +7 to +10 kPa in E10–E25 blends.
2. **Missing Light Ends**: $C_4$ butane ($P^\\text{sat} \\approx 356$ kPa) and $C_5$ isopentane ($P^\\text{sat} \\approx 105$ kPa) constitute 5–15 wt% of commercial gasoline but generate >50% of total vapor pressure.

### Margules Activity Formulation
$$\\ln \\gamma_\\text{EtOH} = A_\\text{EtOH} (1 - x_\\text{EtOH})^2, \\quad A_\\text{EtOH} = 1.80$$
ASTM D5191 correlation:
$$DVPE = 0.965 P_\\text{tot} - 3.78 \\text{ kPa}$$

### Evaluation Results (55 Gasolines, LOO CV)
| Model | DVPE MAE (kPa) | DVPE RMSE (kPa) | Bias (kPa) |
|---|---:|---:|---:|
| Ideal Raoult Baseline | 18.50 | 22.45 | -16.80 |
| **Activity-Corrected DVPE Model (LOO)** | **4.39** | **5.46** | **-0.04** |

---

## 4. Conformal Flash Margin Recomputed (G7)

In `configs/uncertainty_v3.yaml`, the flash point safety margin is established at **+4.20 °C**:
- Grounded in Task R2 strict LOSO (2.28 °C MAE, 5.28 °C worst-case on lab mixtures).
- Grounded in Task R3 method segregation on 70 real transportation fuels: ASTM D93 closed-cup testing exhibits a systematic conservative offset (-8.12 °C bias) due to stirring vapor dilution.
- A +4.20 °C margin guarantees >= 95% compliance against ASTM D93 closed-cup specifications with **0.0% dangerous over-predictions** across commercial aviation fuels.

---

## 5. Diesel Cetane & Cold-Flow Hindcasts (G8–G9)

Evaluated across all 16 real diesel fuels:
- **DCN Blending**: Linear-by-volume baseline achieves 5.62 DCN MAE; LOO recalibration reduces error to **4.34 DCN MAE**.
- **Cloud Point (Coutinho Wax Model)**: Eutectic solid-solution model achieves **10.48 °C MAE**. Correctly flags the cold-flow vulnerability in the Sandia 4-component surrogate (+2.43 °C WAT vs -12 °C real ULSD).
- **Pour Point**: Correlated model ($WAT - 4.5$ °C) achieves **10.69 °C MAE**.

---

## 6. Viscosity & Density Integration with Cleaned ThermoML (G10–G11)

Cleaned ThermoML archive (44,575 rows: 19,404 liquid density, 11,861 liquid viscosity, 13,310 vapor pressure) was integrated:
- **Density**: Ideal volume blending achieves **12.94 kg/m³ MAE** (~1.5% relative error) across all 153 real fuels (82 jet, 55 gasoline, 16 diesel).
- **Viscosity**: Refutas kinematic viscosity achieves ~10.5–11.2% relative error across diesel and jet fuels.
