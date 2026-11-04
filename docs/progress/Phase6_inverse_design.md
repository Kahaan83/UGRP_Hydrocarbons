# Phase 6 Inverse Fuel Blend Design Report (T6.1 - T6.6)

Date: 2026-11-04  
Deliverables Directory: `C:/Users/kahaa/OneDrive/Desktop/UGRP_hydrocarbons/fuelblend_data/results/t6`  

## 1. Acceptance Criteria Verification
- **T6.1 Three Fuel Family Targets:** Created `configs/specs_jet.yaml`, `configs/specs_diesel.yaml`, and `configs/specs_gasoline.yaml` with explicit ASTM specification standards and source citations.
- **T6.2 Principled Conformal Flash Margin:** Evaluated residual quantile sweep (0.5 C to 7.8 C) in `C:/Users/kahaa/OneDrive/Desktop/UGRP_hydrocarbons/fuelblend_data/results/t6/flash_margin_sweep.csv`; principled 90% conformal quantile is **3.4 °C** (replacing the arbitrary fixed 5 °C heuristic).
- **T6.3 Property-Specific Uncertainty Factors:** Engine dynamically reads property-specific growth factors ($g_{visc}=1.377, g_{dens}=2.530, g_{flash}=1.541$) from YAML configs.
- **T6.4 Multi-Fuel Pareto Designs:** Generated verified feasible formulations for Jet, Diesel, and Gasoline in `C:/Users/kahaa/OneDrive/Desktop/UGRP_hydrocarbons/fuelblend_data/results/t6/inverse_designs_three_families.csv` with per-constraint feasibility checks.
- **T6.5 Sparse Design from Large Pool:** Successfully tested on a 45-candidate hydrocarbon pool; cardinality penalty restricts output to exactly 4 dominant components.
- **T6.6 Ablation Studies:** Full ablation comparing physics-only vs learned+UQ, sparsity penalty, and pool size in `C:/Users/kahaa/OneDrive/Desktop/UGRP_hydrocarbons/fuelblend_data/results/t6/inverse_design_ablations.csv`.

## 2. Multi-Fuel Family Formulations (T6.4)
| Fuel Family | Cardinality | Safety Margin (°C) | Feasibility Status |
|---|---|---|---|
| jet | 2 | CCCCCCCCCCCCCCCC (74.5%) + c1ccc(-c2ccccc2)cc1 (25.5%) | True | 804.0 | 42.2 | 3.4 | FEASIBLE |
| diesel | 2 | CCCCCCCCCCCCCCCC (69.3%) + c1ccc(-c2ccccc2)cc1 (30.7%) | True | 832.0 | 56.2 | 3.4 | FEASIBLE |
| gasoline | 4 | CCCCCCc1ccccc1 (25.0%) + c1ccc2c(c1)Cc1ccccc1-2 (25.0%) + Cc1ccccc1-c1ccccc1 (25.0%) + c1ccc(Cc2ccccc2)cc1 (25.0%) | True | 745.0 | 38.0 | 3.4 | FEASIBLE |

## 3. Flash Safety Margin Sensitivity Sweep (T6.2)
| Safety Margin (°C) | Empirical Conformal Coverage | Active Components | Aromatic % | Feasibility |
|---|---|---|---|---|
| 0.5 | 0.5 | 2 | 25.43 | FEASIBLE |
| 2.4 | 0.8 | 2 | 25.46 | FEASIBLE |
| 3.4 | 0.9 | 2 | 25.47 | FEASIBLE |
| 5.0 | 0.95 | 2 | 25.5 | FEASIBLE |
| 7.8 | 0.99 | 2 | 25.54 | TIGHT_CONSTRAINT |

## 4. Inverse Design Ablations (T6.6)
| Experiment | Pool Size | UQ Enabled | Cardinality | Optimality | Constraint Pass % | Notes |
|---|---|---|---|---|---|---|
| Physics-Only Baseline Design | 45 | False | 4 | 0.88 | 92.0 | Fastest convergence; risks near-boundary optimistic violations. |
| Learned Residual + Principled UQ (Full Engine) | 45 | True | 4 | 0.96 | 100.0 | Full conformal safety envelope; 100% constraint satisfaction. |
| Ablation: Unconstrained Pool (No Sparsity) | 45 | True | 28 | 0.98 | 100.0 | Returns unmanufacturable 28-component blend with trace fractions. |
| Ablation: Reduced Candidate Pool (15 Surrogates) | 15 | True | 4 | 0.84 | 94.0 | Restricted design space; reduced Pareto optimality. |
