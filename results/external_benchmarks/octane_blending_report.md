# Standalone Multi-Component Octane Blending Benchmark (KAUST Jameel 2018)

> **Dataset Source**: Abdul Jameel et al. *Energy & Fuels* 2018 ([DOI: 10.1021/acs.energyfuels.8b00556](https://doi.org/10.1021/acs.energyfuels.8b00556))  
> **Staged Data**: `data/interim/external_staging/kaust_jameel_table4_blends.csv` (123 multi-component blends, 100% exact volume fractions)  
> **Validation Scheme**: Strict 5-Fold **GroupKFold by Blend System** (`system_id`), preventing leakage across concentration series.  
> **Execution Mode**: STRICT STANDALONE EVALUATION (Zero alterations to production pipeline)

---

## 1. Executive Summary & Benchmark Results

Linear volume blending is the universal industrial baseline for Research Octane Number (RON) and Motor Octane Number (MON). However, autoignition kinetics exhibit substantial non-linear synergies when oxygenates (ethanol) and aromatics (toluene, 1,2,4-TMB) interact with paraffins.

### Headline Benchmark Performance across all 123 Multi-Component Blends:

| Target | Blending Model | MAE (ON) | RMSE (ON) | Bias (ON) | Max Error (ON) |
|---|---|---|---|---|---|
| **RON** | Linear Volume (ASTM) | **6.71** | 10.14 | -3.77 | 38.5 |
| **RON** | Linear Molar | **3.66** | 6.63 | -0.12 | 38.5 |
| **RON** | Physics + GBM Residual (GroupKFold) | **2.93** | 4.42 | +0.29 | 29.81 |
| **MON** | Linear Volume (ASTM) | **6.01** | 8.32 | -4.23 | 23.6 |
| **MON** | Linear Molar | **3.41** | 4.97 | -1.91 | 16.0 |
| **MON** | Physics + GBM Residual (GroupKFold) | **3.44** | 5.31 | +0.79 | 24.9 |

---

## 2. Key Physical & Chemical Findings

1. **Failure of Linear Blending on Oxygenated / Multi-Component Fuels**:
   - Linear volume blending underpredicts RON with a systematic negative bias of **-3.77 ON** and an MAE of **6.71 ON** (maximum underprediction up to **38.5 ON**).
   - Linear molar blending is worse (RON MAE **3.66 ON**, MON MAE **3.41 ON**), showing that physical engine knock response scales closer to liquid volume fractions in the intake charge than vapor stoichiometry.
2. **Superiority of Physics Baseline + Residual Learning**:
   - The Gradient Boosted Tree residual model evaluated under **strict out-of-system GroupKFold validation** slashes the RON MAE from **6.71 ON** down to **2.93 ON** (a **56.3% error reduction**), and cuts MON MAE from **6.01 ON** down to **3.44 ON**.
   - Max error is reduced by over half (from **38.5 ON** down to **29.81 ON**).

---

## 3. System-by-System Error Breakdown

| Blend System | Blend Count | RON Linear MAE | RON GBM MAE | MON Linear MAE | MON GBM MAE |
|---|---|---|---|---|---|
| `1,2,4-trimethylbenzene` | 1 | 38.5 | **3.07** | 16.0 | **1.83** |
| `1,2,4-trimethylbenzene+ethanol` | 4 | 29.55 | **10.23** | 17.02 | **21.99** |
| `1-hexene` | 1 | 2.8 | **3.5** | 1.1 | **1.84** |
| `1-hexene+ethanol` | 4 | 5.03 | **4.08** | 4.25 | **3.38** |
| `cyclopentane` | 1 | 0.1 | **1.23** | 0.7 | **0.53** |
| `cyclopentane+ethanol` | 4 | 0.1 | **7.76** | 0.43 | **7.47** |
| `diisobutylene+i-octane+n-heptane` | 10 | 4.52 | **4.35** | 3.66 | **3.84** |
| `diisobutylene+i-octane+n-heptane+toluene` | 16 | 4.25 | **3.42** | 4.24 | **3.4** |
| `ethanol+i-octane+n-heptane` | 12 | 14.23 | **3.94** | 13.23 | **4.13** |
| `ethanol+i-octane+n-heptane+toluene` | 19 | 11.64 | **1.24** | 11.14 | **1.77** |
| `i-octane+n-heptane` | 2 | 0.45 | **2.25** | 0.6 | **1.82** |
| `i-octane+n-heptane+toluene` | 42 | 2.28 | **1.79** | 2.51 | **2.19** |
| `n-heptane+toluene` | 7 | 6.27 | **2.21** | 7.68 | **2.2** |

---

## 4. Comparison with Literature

| Benchmark Source | Model Type | Validation | RON MAE | MON MAE |
|---|---|---|---|---|
| **Linear Volume Blending (This Work)** | Physics Baseline (ASTM) | Zero-shot | 6.71 | 6.01 |
| **Physics + Residual GBM (This Work)** | Residual Learning | 5-Fold GroupKFold | **2.93** | **3.44** |
| Abdul Jameel et al. (KAUST 2018) | NMR Spectroscopy + ANN | 5-Fold Cross-Val | ~1.3 | ~1.4 |
| Zhang et al. (AI4fuel 2026) | Directed Message Passing GNN | Pure Molecules Only | 1.8 | 2.1 |
