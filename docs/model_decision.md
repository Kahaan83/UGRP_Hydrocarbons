> v4 NOTICE: Sections 1-4 are v2/v3-era and rely on the v3 corpus where noted; Section 6 governs; see docs/FULL_REPORT_v4.md and docs/claims_register.md

# Master Model Strategy Decision Note

**Project:** UGRP Hydrocarbons v3  
**Status:** Canonical Architecture Decision  
**Scope:** Model Selection Across Properties and Composition Regimes  

---

## 1. Unequivocal Model Architecture Decision

Based on empirical evidence across Phase A (triage & LOO recalibration), Phase C (property baseline fixes), Phase D (stacked evaluation, regime analysis, truncation sensitivity), and Task R1 (canonical viscosity leaderboard on 57 CRC 647 real fuels), the canonical modeling strategy is decided as follows:

| Property | Primary Production Model | Secondary / High-Precision Predictive Model | Empirical Evidence Table |
|---|---|---|---|
| **Kinematic Viscosity** ($\nu$) | **2-Parameter Recalibrated Refutas** (Parsimonious recalibrated physics baseline) | **Bilinear Pairwise Interaction** (Highest LOO predictive accuracy) | `results/v3/R1/viscosity_leaderboard.csv`<br>`results/v2/D1/bilinear_hindcast.csv`<br>`results/v2/A2/recalibration_loo.csv` |
| **Density** ($\rho$) | **Ideal Volume Blending** ($\rho_{ideal} + \Delta\rho$) | Stacked Ideal Volume + Shrunk Ridge | `results/v2/D4/stacked_eval.csv` |
| **Flash Point** ($T_{flash}$) | **Transparent Offset Correlation / Le Chatelier** | Le Chatelier + Shrunk HGB Residual (LOSO) | `results/v2/C2/flash_loso_v2.csv`<br>`results/v2/D4/stacked_eval.csv` |
| **Freezing Point** ($T_{freeze}$) | **Transparent Offset SLE WAT** ($WAT - 14.58$ °C) | Linear Recalibrated SLE | `results/v2/C4/cold_flow_v2.csv`<br>`results/v2/D4/stacked_eval.csv` |
| **Heat of Combustion** ($NCV$) | **Group Additivity GCV Baseline (v2)** | Stacked GCV + Shrunk Ridge | `results/v2/C1/hoc_eval.csv`<br>`results/v2/D4/stacked_eval.csv` |

---

## 2. Viscosity Model Reconciliation: Bilinear Pairwise vs Recalibrated Refutas

To resolve any apparent tension in model selection for kinematic viscosity:

1. **Predictive Accuracy Winner:** **Bilinear Pairwise (LOO)** achieves the lowest overall prediction error across all 57 real aviation fuels:
   - **Overall (all T):** log-RMSE **0.1242** ($[0.1041, 0.1507]$), log-MAE **0.0949** ($[0.0835, 0.1081]$), MAPE **9.56%**.
   - **Cryogenic ($T = -40$ °C):** log-RMSE **0.1146** ($[0.0929, 0.1342]$), log-MAE **0.0894** ($[0.0716, 0.1078]$), MAPE **8.45%**.
   - By capturing low-rank pairwise intermolecular friction terms among constituent hydrocarbon classes, Bilinear Pairwise provides the most accurate numerical predictions when trained on real fuels under Leave-One-Fuel-Out CV.

2. **Physical Parsimony Winner:** **2-Parameter Recalibrated Refutas (LOO)** is the most parsimonious transparent physical model:
   - **Overall (all T):** log-RMSE **0.1411** ($[0.1174, 0.1667]$), log-MAE **0.1042** ($[0.0890, 0.1205]$), MAPE **10.35%**.
   - **Cryogenic ($T = -40$ °C):** log-RMSE **0.1705** ($[0.1388, 0.1982]$), log-MAE **0.1352** ($[0.1085, 0.1621]$), MAPE **12.85%**.
   - With only **2 parameters** ($\ln \nu = \ln \nu_{Refutas} + a_0 + a_1 \frac{1000}{T}$), it captures $68.1\%$ of the error reduction without requiring any interaction matrix weights, guaranteeing robust monotonicity and zero risk of matrix overfitting.

3. **ASTM D341 Walther Curvature Baseline:**
   - Standard ASTM D341 double-logarithmic curvature ($\log_{10}\log_{10}(\nu + 0.7) = A - B\log_{10}T$) achieves log-RMSE **0.4012** overall and **0.6082** at -40 °C (MAPE 43.72%) with **zero empirical parameters fitted on real fuels**, improving over raw Refutas (0.6879 at -40 °C, MAPE 48.12%).
   - Combined with recalibrated anchors, Walther achieves log-RMSE **0.1642** overall and **0.2243** at -40 °C.

4. **Catastrophic Cryogenic Domain Collapse of Neural Residuals:**
   - Neural residuals trained strictly on ambient laboratory data (**DeepSets** and **Pairwise**) suffer disastrous out-of-domain failure when transferred zero-shot to real fuels at cryogenic temperatures:
     - **DeepSets Residual:** log-RMSE **0.4037** overall, collapsing to **0.6310** at -40 °C (MAPE 45.00%).
     - **Pairwise Residual:** log-RMSE **0.4504** overall, collapsing to **0.6870** at -40 °C (MAPE 48.15%, worse than raw physics).
   - In `results/v2/A6/uncertainty_growth_refit.csv`, DeepSets claimed $92.5\%$ nominal coverage on real fuels, but actual empirical coverage collapsed to **$6.43\%$** (viscosity) and **$15.79\%$** (density). Zero-shot ambient neural models must **never** be used for cryogenic fuel predictions without in-domain calibration.

---

## 3. Evidence Across Regimes: Ambient Lab Binaries vs Cryogenic Real Fuels

The empirical data across blend complexity buckets $N \in [1, 40+]$ (`results/v2/D5/error_vs_n.csv`) proves an **architectural inversion** between regimes:

1. **Ambient Laboratory Mixtures ($N = 2, 3$; $T \in [293, 323]$ K)**:
   - Learned residual models succeed: DeepSets achieves log-MAE $0.0372$ ($[0.0289, 0.0469]$) and Pairwise achieves log-MAE $0.0342$ ($[0.0266, 0.0432]$), improving over raw Refutas ($0.0616$) and Grunberg-Nissan ($0.0669$, `results/v2/A4/gn_baselines.csv`).
   - *Scope for learning*: Neural and tree-based residuals are strictly justified here to interpolate excess thermodynamic mixing friction and activity coefficients.

2. **Cryogenic Real Multi-Component Fuels ($N \ge 20$; $T \in [-40, +20]$ °C)**:
   - Real fuel complexity and cryogenic conditions invert model behavior: unconstrained neural models extrapolate erratically, while physics-anchored baselines (Recalibrated Refutas and Bilinear Pairwise) maintain strict monotonic bounds and low error.

---

## 4. Strict Boundary Rules for Downstream Tasks (Phase E & F)

1. **Production Deployment Recommendations**:
   - **Regulatory Certification & Explainable Screening:** Deploy **2-Parameter Recalibrated Refutas** for inverse design and constraint checking where regulatory auditability and zero black-box weights are mandatory.
   - **Mission-Critical Precision:** Deploy **Bilinear Pairwise Interaction (LOO)** when high-precision operational margins at -40 °C are required for aviation fuel safety.
   - **Prohibited:** Zero-shot ambient neural network residuals (DeepSets / Pairwise) are strictly barred from production evaluation of cryogenic real fuels.

2. **Shrinkage Prior on Learned Residuals**:
   - Where a learned residual is employed (e.g. the stacked model architecture in `results/v2/D4/stacked_eval.csv`), it must carry an explicit $L_2$ shrinkage prior ($\alpha \ge 50.0$) ensuring that in data-sparse or out-of-distribution regimes, the residual collapses to zero rather than extrapolating arbitrarily.

3. **Surrogate Cardinality Rule (from Task D6)**:
   - Inverse design blends must not use over-simplistic 2- or 3-component surrogates (e.g., dodecane + tetralin), which fail freezing point and boiling range constraints.
   - Per `results/v2/D6/truncation.csv`, a surrogate complexity of **$K = 10$ to $20$ components with a lumped heavy tail** captures $>95\%$ of blend mass and keeps bulk property drift below $0.3\%$.

---

## 5. v4 per-property routing of the prediction API (docs only)

`src/fuelblend/api/predict_blend.py` routes as below; the table mirrors `MARGIN_KEYS` there. Intervals use `applied_margin` from `configs/uncertainty_v4.yaml` only; `inf` means no interval (NaN) and `in_domain=False`.

| Property | Model | Margin key(s) (`MARGIN_KEYS`) | Certified scope |
|---|---|---|---|
| Kinematic/dynamic viscosity | Recalibrated Refutas (mixtures only; raw Refutas for single components and outside the window) | `viscosity_kin_cSt.crc647.recalibrated_refutas.all_T` (ln-ratio) | jet fuels only; -40..20 C |
| Density | Ideal volume mixing | `density_kg_m3.crc647_ideal_volume.all_T` | jet fuels only; -0.2..69.5 C |
| Flash point | Le Chatelier (ideal) + Antoine + LFL | max of `flash_point_C.lab_mixtures.ideal` and `flash_point_C.real_fuels_D93_reconstructed.ideal` (never the 0-floored D93 alone) | jet fuels only; `lo` is the safe bound |

Uncertified (no valid margin; value None, interval NaN, `in_domain=False`): freezing, smoke point, RON/MON, DCN, DVPE (data units unconfirmed), pour point.
The mixture-residual GBM is opt-in (`use_mixture_residual=True`), has no margin, and is never certified. `fuel_type='jet'` is trusted from the caller.

---

## 6. v4 decision summary (Task 19; extends section 5, does not replace it)

Evidence: `docs/FULL_REPORT_v4.md` and `results/v4/*/README.md`. Sections 1-4 above are v2/v3-era decisions; where they disagree with this section, this section governs. The v3 corpus and design optimizer were formula-derived and are quarantined (`data/quarantine/`, `results/v4/P1/AUDIT.md`); decisions that rested on them are void. The routing table in section 5 is unchanged.

| Property | v4 decision | Evidence | What it does NOT say |
|---|---|---|---|
| Viscosity | Production: recalibrated Refutas (margin 0.347 ln). Learned alternates (structured ridge 0.0826, bilinear 0.0901, GP 0.0888 dev log-RMSE) beat it with paired CI excluding 0 but are NOT ranked among themselves; the ridge is kept behind an envelope guard | `results/v4/V/README.md` | GP under-covers (0.77 dev / 0.69 locked); lab zero-shot does not fix -40 C; encoder lost (0.0862 vs 0.0833, CI includes 0) |
| Flash | Ideal Le Chatelier (UNIFAC is worse: dMAE +0.45 [+0.11, +0.78]); margin = max(lab 3.70, D93 0.00) | `results/v4/F/README.md`, `results/v4/F/flash_input_sensitivity.csv` | Lab margin is set by the mesitylene system (0.31 C without it); D93 evidence is reconstructed compositions only; D56 n = 1 |
| Density | Ideal volume mixing, margin 22.25 kg/m3 | `results/v4/B/README.md` | Bias -11 kg/m3, cause unknown |
| Net HoC | Group additivity, margin 0.194 MJ/kg | `results/v4/B/README.md` | Possible error cancellation (measured-H diagnostic) |
| Freezing | No model demonstrated; offset/Coutinho do not beat a constant mean on CRC fuels; ideal SvL is good only for n-alkane blends (8 systems, margin inf) | `results/v4/B/README.md`, `results/v4/SLE/README.md` | Not certifiable; model form changes design outcomes |
| Smoke point | AJFTD bulk regression is a baseline only (1.29 vs 1.47 mm grouped; no skill on the broad set); margin inf | `results/v4/SP/README.md` | Spec limits secondary-source / UNVERIFIED |
| RON/MON | bRON (base + one blendstock) only; learned models beat base-unchanged; no CI-supported winner among them | `results/v4/O/README.md` | Not for multi-component pools; locked thin (12 rows, 4 groups) |
| DVPE, pour, cloud, DCN | Uncertified, no v4 model or margin | `configs/uncertainty_v4.yaml` | DVPE data units unconfirmed |
| Diesel (any property) | Nothing validated: 0 qualifying fuels | `results/v4/P2/coverage_matrix_summary.md` | |
| Designs | 0 of 13 certified; certification is jet-only | `results/v4/D/README.md` | 8 of 13 are undetermined, not failed |
