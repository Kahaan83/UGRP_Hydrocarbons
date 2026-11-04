# Experimental Validation Protocol: Laboratory Testing of Designed Hydrocarbon Fuel Blends

**Document Version:** 2.0 (Post-Audit Master Plan)  
**Date:** November 5, 2026  
**Target Audience:** Chemical Engineering, Fuels & Combustion, or Petroleum Testing Laboratory  
**Governing Standard Frameworks:** ASTM D1655 (Aviation Turbine Fuel), ASTM D975 (Diesel Fuel Oils), ASTM D4814 (Automotive Spark-Ignition Engine Fuel)

---

## 1. Executive Summary and Validation Objective

This document provides an actionable, rigorous experimental testing protocol for laboratory validation of the model-designed hydrocarbon fuel blends formulated in Phase E of the UGRP Hydrocarbons initiative.

The objective is to validate that sparse, manufactured hydrocarbon mixtures (3 to 5 pure components) reproduce the complex thermophysical behavior, volatility profiles, and certification specifications of real petroleum fuels within the 90% confidence envelopes established by the project's conformal uncertainty framework (`configs/uncertainty_v2.yaml`).

Three high-priority formulations have been selected for physical compounding and laboratory testing:
1. **Blend 1 (Jet A Minimal 3-Component Formulation):** Demonstrates that adding naphthenic decalin to isododecane and tetralin satisfies freezing point and density constraints without exceeding aromatics limits.
2. **Blend 2 (Jet A Realistic 4-Component Formulation):** Introduces normal alkane ($n$-dodecane) to better mimic petroleum molecular class distributions and burning characteristics.
3. **Blend 3 (Diesel No. 2 Realistic 5-Component Formulation):** Balances high cetane ($n$-alkanes), low cloud point (iso-alkane HMN), high density (bicyclohexyl), and aromatic solvency (tetralin).

---

## 2. Compounding Specifications & Batch Preparation

### 2.1 Component Purity & Supplier Requirements
All blending stocks must meet analytical reagent grade requirements:
- **Purity:** $\ge 99.0\%$ GC-FID purity for all pure components.
- **Water Content:** $< 50$ ppm by Karl Fischer titration (ASTM D6304). Water contamination severely biases cold-flow cloud/freeze measurements.
- **Peroxides:** $< 5$ ppm (active oxygen) in olefinic and ether components (inhibited if necessary).

| Compound Name | IUPAC / Structure | CAS Registry | Formula | Molecular Weight (g/mol) | Pure Density @ 20 °C (kg/m³) | Supplier Example (Sigma / TCI / Fisher) |
|---|---|---|---|---|---|---|
| **Isododecane** | 2,2,4,6,6-pentamethylheptane | 13475-82-6 | C12H26 | 170.34 | 743.0 | Sigma-Aldrich #324418 ($\ge 99\%$) |
| **$n$-Dodecane** | Dodecane | 112-40-3 | C12H26 | 170.34 | 748.7 | Sigma-Aldrich #D221104 ($\ge 99\%$) |
| **Decalin** | Decahydronaphthalene (cis/trans mix) | 91-17-8 | C10H18 | 138.25 | 896.5 | Sigma-Aldrich #D1300 ($\ge 99\%$) |
| **Tetralin** | 1,2,3,4-Tetrahydronaphthalene | 119-64-2 | C10H12 | 132.20 | 970.2 | Sigma-Aldrich #T2680 ($\ge 99\%$) |
| **$n$-Hexadecane** | Cetane | 544-76-3 | C16H34 | 226.44 | 773.3 | Sigma-Aldrich #H6703 ($\ge 99\%$) |
| **$n$-Octadecane** | Octadecane | 593-45-3 | C18H38 | 254.50 | 776.8 (melted @ 30 °C) | Sigma-Aldrich #O652 ($\ge 99\%$) |
| **Heptamethylnonane (HMN)**| 2,2,4,4,6,8,8-heptamethylnonane | 4390-04-9 | C16H34 | 226.44 | 785.0 | TCI America #H0363 ($\ge 98\%$) |
| **Bicyclohexyl** | Dicyclohexyl | 92-51-3 | C12H22 | 166.30 | 880.0 | Sigma-Aldrich #B35405 ($\ge 99\%$) |

---

### 2.2 Gravimetric Batch Compounding Protocol (500 mL Nominal Batches)
Because volumetric dispensing is subject to thermal expansion and excess volume of mixing, **all compounding must be conducted gravimetrically** on an analytical balance calibrated to $\pm 0.01$ g:

```
                            GRAVIMETRIC COMPOUNDING WORKFLOW
   ┌─────────────────────────────────────────────────────────────────────────┐
   │ 1. Zero a clean, oven-dried 500 mL borosilicate glass volumetric flask   │
   │    or Boston round bottle with PTFE-lined cap on an analytical balance.  │
   │ 2. Dispense high-boiling non-volatile stocks first (Decalin, Tetralin,   │
   │    Bicyclohexyl) to minimize vapor loss.                                │
   │ 3. Add medium alkanes (n-Dodecane, Hexadecane, HMN).                     │
   │ 4. For n-Octadecane (mp = 28 °C), warm bottle gently in a 35 °C water    │
   │    bath until liquid, dispense rapidly via pre-warmed glass pipette.     │
   │ 5. Add volatile fractions (Isododecane) last; immediately seal cap.      │
   │ 6. Invert gently 20 times to homogenize; equilibrate at 20 °C for 1 h.   │
   └─────────────────────────────────────────────────────────────────────────┘
```

#### Batch Recipes for 500 mL Test Samples:

##### Blend 1: Jet A Minimal 3-Component Formulation
- Target Blend Density @ 15 °C: 794.7 kg/m³
- Total Target Batch Mass: **397.35 g**

| Component | Target Mass Fraction ($w_i$) | Required Mass ($\pm 0.02$ g) | Approximate Volume @ 20 °C |
|---|---|---|---|
| **Isododecane** | 0.6500 | **258.28 g** | 347.6 mL |
| **Decalin** (cis/trans) | 0.2000 | **79.47 g** | 88.6 mL |
| **Tetralin** | 0.1500 | **59.60 g** | 61.4 mL |
| **Total** | **1.0000** | **397.35 g** | **~497.6 mL** |

##### Blend 2: Jet A Realistic 4-Component Formulation
- Target Blend Density @ 15 °C: 794.6 kg/m³
- Total Target Batch Mass: **397.30 g**

| Component | Target Mass Fraction ($w_i$) | Required Mass ($\pm 0.02$ g) | Approximate Volume @ 20 °C |
|---|---|---|---|
| **Isododecane** | 0.5800 | **230.43 g** | 310.1 mL |
| **$n$-Dodecane** | 0.1200 | **47.68 g** | 63.7 mL |
| **Decalin** (cis/trans) | 0.1500 | **59.60 g** | 66.5 mL |
| **Tetralin** | 0.1500 | **59.60 g** | 61.4 mL |
| **Total** | **1.0000** | **397.30 g** | **~497.7 mL** |

##### Blend 3: Diesel No. 2 Realistic 5-Component Formulation
- Target Blend Density @ 15 °C: 834.7 kg/m³
- Total Target Batch Mass: **417.35 g**

| Component | Target Mass Fraction ($w_i$) | Required Mass ($\pm 0.02$ g) | Approximate Volume @ 20 °C |
|---|---|---|---|
| **$n$-Hexadecane** | 0.2100 | **87.64 g** | 113.3 mL |
| **$n$-Octadecane** | 0.0700 | **29.21 g** | 37.6 mL |
| **Heptamethylnonane (HMN)**| 0.2700 | **112.68 g** | 143.5 mL |
| **Bicyclohexyl** | 0.2500 | **104.34 g** | 118.6 mL |
| **Tetralin** | 0.2000 | **83.47 g** | 86.0 mL |
| **Total** | **1.0000** | **417.35 g** | **~499.0 mL** |

---

## 3. Standard Test Methods & Predicted 90% Confidence Envelopes

Every physical measurement must be performed in triplicate using standard ASTM/IP apparatus. The table below lists the model-predicted values alongside the **data-driven 90% confidence envelopes** from `configs/uncertainty_v2.yaml`.

### 3.1 Jet A Candidate Formulations (ASTM D1655 Specifications)

| Test Property | ASTM Standard Method | Target Specification (ASTM D1655) | Blend 1 Predicted Value | Blend 1 90% Confidence Envelope | Blend 2 Predicted Value | Blend 2 90% Confidence Envelope | Acceptance Verdict |
|---|---|---|---|---|---|---|---|
| **Density @ 15 °C** | ASTM D4052 (Digital Density Meter) | 775.0 to 840.0 kg/m³ | **794.7 kg/m³** | **[792.2, 797.2]** kg/m³ | **794.6 kg/m³** | **[792.1, 797.1]** kg/m³ | **PASS** (centered in spec) |
| **Flash Point** | ASTM D56 (Tag Closed Cup) | $\ge 38.0$ °C | **44.8 °C** | **[39.2, 50.4]** °C | **45.1 °C** | **[39.5, 50.7]** °C | **PASS** (safe $> 38$ °C) |
| **Freezing Point** | ASTM D5972 (Phase Transition) | $\le -40.0$ °C | **$-51.3$ °C** | **[$-57.8, -44.8$]** °C | **$-47.1$ °C** | **[$-53.6, -40.6$]** °C | **PASS** (clears $-40$ °C) |
| **Kin. Viscosity @ $-20$ °C** | ASTM D445 (Ubbelohde / Cannon) | $\le 8.00$ cSt | **4.88 cSt** | **[3.71, 6.42]** cSt | **4.86 cSt** | **[3.70, 6.39]** cSt | **PASS** (well below 8.0) |
| **Kin. Viscosity @ $-40$ °C** | ASTM D445 (Low-Temp Viscometry) | Report | **10.42 cSt** | **[7.92, 13.70]** cSt | **10.35 cSt** | **[7.87, 13.61]** cSt | **PASS** (fluid, no gel) |
| **Kin. Viscosity @ $+20$ °C** | ASTM D445 | Report | **1.62 cSt** | **[1.23, 2.13]** cSt | **1.61 cSt** | **[1.22, 2.12]** cSt | Baseline sanity check |
| **Kin. Viscosity @ $+40$ °C** | ASTM D445 | Report | **1.15 cSt** | **[0.87, 1.51]** cSt | **1.14 cSt** | **[0.87, 1.50]** cSt | Baseline sanity check |
| **Net Heat of Comb. (NHOC)**| ASTM D4809 (Bomb Calorimetry) | $\ge 42.80$ MJ/kg | **43.43 MJ/kg** | **[43.15, 43.71]** MJ/kg | **43.49 MJ/kg** | **[43.21, 43.77]** MJ/kg | **PASS** (high gravimetric) |
| **Aromatics Content** | ASTM D1319 (FIA) or ASTM D6379 | $\le 25.0$ vol% | **13.9 vol%** | **[13.0, 14.8]** vol% | **13.9 vol%** | **[13.0, 14.8]** vol% | **PASS** (15.0 wt% tetralin) |
| **Distillation T10** | ASTM D86 (Automated Batch Dist.) | $\le 205.0$ °C | **175.2 °C** | **[168.4, 182.0]** °C | **176.4 °C** | **[169.6, 183.2]** °C | **PASS** |
| **Distillation T90** | ASTM D86 (Automated Batch Dist.) | Report ($\le 300$ °C)| **206.5 °C** | **[199.3, 213.7]** °C | **206.8 °C** | **[199.6, 214.0]** °C | **PASS** |

---

### 3.2 Diesel No. 2 Candidate Formulation (ASTM D975 Specifications)

| Test Property | ASTM Standard Method | Target Specification (ASTM D975) | Blend 3 Predicted Value | Blend 3 90% Confidence Envelope | Acceptance Verdict |
|---|---|---|---|---|---|
| **Density @ 15 °C** | ASTM D4052 | 820.0 to 845.0 kg/m³ | **834.7 kg/m³** | **[832.2, 837.2]** kg/m³ | **PASS** (compliant) |
| **Flash Point** | ASTM D93 (Pensky-Martens Closed Cup) | $\ge 52.0$ °C | **83.2 °C** | **[77.6, 88.8]** °C | **PASS** (well above 52 °C) |
| **Cloud Point** | ASTM D2500 (Optical Detection) | $\le 0.0$ °C | **$-3.4$ °C** | **[$-6.4, -0.4$]** °C | **PASS** (upper bound $< 0$ °C) |
| **Kin. Viscosity @ 40 °C** | ASTM D445 | 1.90 to 4.10 cSt | **2.89 cSt** | **[2.51, 3.33]** cSt | **PASS** (centered in spec) |
| **Kin. Viscosity @ 20 °C** | ASTM D445 | Report | **4.85 cSt** | **[4.21, 5.59]** cSt | Expected diesel curve |
| **Kin. Viscosity @ $-20$ °C** | ASTM D445 | Report | **24.50 cSt** | **[21.2, 28.3]** cSt | Mobile fluid |
| **Derived Cetane Number (DCN)**| ASTM D6890 (IQT Combustion) | $\ge 40.0$ | **48.5** | **[42.7, 54.3]** | **PASS** (exceeds 40.0) |
| **Aromatics Content** | ASTM D1319 / ASTM D5186 | $\le 35.0$ wt% | **20.0 wt%** | **[18.5, 21.5]** wt% | **PASS** (20.0% tetralin) |
| **Distillation T90** | ASTM D86 | 282.0 to 338.0 °C | **294.5 °C** | **[287.3, 301.7]** °C | **PASS** (within window) |

---

## 4. Laboratory Execution & Handling Protocol

### 4.1 Execution Sequence & Aliquot Splitting
To prevent cross-contamination and evaporative weathering, follow this strict testing sequence:

```
                            ALIQUOT DISPENSING SEQUENCE
   ┌────────────────────────────────────────────────────────────────────────┐
   │ 500 mL Compounded Stock                                                │
   │  ├── Aliquot A (100 mL): Volatile tests (ASTM D56 Flash Point first,    │
   │  │                       followed immediately by ASTM D4052 Density)   │
   │  ├── Aliquot B (150 mL): Low-Temperature Viscometry (ASTM D445 at     │
   │  │                       -40, -20, 20, 40 °C)                          │
   │  ├── Aliquot C (50 mL):  Phase Transition Freezing Point (ASTM D5972)  │
   │  │                       and Cloud Point (ASTM D2500)                  │
   │  ├── Aliquot D (100 mL): Distillation Curve (ASTM D86 Automated)        │
   │  └── Aliquot E (100 mL): Combustion Calorimetry (ASTM D4809) / DCN     │
   └────────────────────────────────────────────────────────────────────────┘
```

1. **Flash Point First:** Perform Tag Closed Cup (ASTM D56) or Pensky-Martens (ASTM D93) immediately upon unsealing Aliquot A. Even 1% evaporative loss of isododecane will artifactually elevate the measured flash point by 2–4 °C.
2. **Cold-Flow Cryogenic Conditioning:** For ASTM D445 at $-40$ °C and ASTM D5972 freezing point, dry nitrogen purge must blanket the sample tube to prevent ambient moisture condensation/frost formation, which can seed heterogeneous ice nucleation.
3. **Viscosity Cleaning:** Rinse viscometer tubes twice with analytical heptane, followed by acetone, and dry with filtered nitrogen. Polyaromatic residue from tetralin will cause capillary drag if not thoroughly degreased.

### 4.2 Health, Safety, and Hazardous Waste Management
- **Personal Protective Equipment (PPE):**
  * Nitrile gloves ($\ge 8$ mil thickness) or fluorocarbon (Viton) gloves. Tetralin and decalin penetrate thin 4 mil disposable nitrile within 15–20 minutes.
  * Chemical splash goggles and lab coat.
  * Conduct all compounding and volatile transfers in a certified chemical fume hood (face velocity $\ge 100$ fpm).
- **Flammability Hazards:**
  * Isododecane has a flash point of ~43 °C (Class II combustible liquid). Eliminate all open flames, sparks, and ungrounded electrical heating mantles.
  * Ground all metal dispensing containers to prevent electrostatic discharge.
- **Toxicity & Exposure:**
  * **Tetralin (119-64-2):** Skin and respiratory irritant; may form organic peroxides upon prolonged contact with air. Inspect bottles for peroxide crystals prior to warming.
  * **Decalin (91-17-8):** Moderate skin irritant; toxic to aquatic life.
- **Waste Disposal:** Collect all solvent flushes, redundant sample aliquots, and bomb calorimetry washings into designated Non-Halogenated Hydrocarbon Flammable Waste carboys. Do not drain dispose.

---

## 5. Result Verification & Discrepancy Reporting

### 5.1 Success Acceptance Criteria
A candidate blend is formally considered **validated** if:
1. Every physical measurement falls within the model's 90% confidence envelope defined in Section 3.
2. All measured properties comply with the governing commercial fuel specifications (ASTM D1655 / D975).
3. Repeatability across triplicate runs satisfies the standard ASTM repeatability criteria ($r$):
   - Density: $\Delta \le 0.1$ kg/m³
   - Flash Point: $\Delta \le 1.0$ °C
   - Freezing Point: $\Delta \le 0.5$ °C
   - Kinematic Viscosity: $\Delta \le 0.35\%$

### 5.2 Discrepancy Protocol
If a measured value falls outside the predicted 90% envelope:
1. Inspect the chromatogram (GC-FID) of the blended sample to verify gravimetric compounding accuracy ($\pm 0.5\%$ relative).
2. Verify water content ($< 50$ ppm).
3. If compounding is verified, report the delta ($\Delta = y_{\text{meas}} - \hat{y}_{\text{pred}}$) to the modeling team to update the residual calibration matrices in `results/v2/`.

---
*Protocol approved by UGRP Hydrocarbons Modeling Team, November 5, 2026.*
