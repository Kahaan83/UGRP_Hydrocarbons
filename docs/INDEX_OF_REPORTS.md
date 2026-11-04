# Master Documentation Index & Catalog of Reports

**Project**: UGRP Hydrocarbons (`fuelblend`)  
**Scope**: Complete directory of all technical reports, progress summaries, audit logs, and implementation plans across v1, v2, and v2.1.  
**Last Updated**: 2026-11-05 (Release Tag `v2.1-2026-11-05`)  

---

## 1. Master Syntheses & Capstone Reports

These documents synthesize the complete findings, architectural decisions, and final deliverables across the entire project.

| Document Title | File Path | Generation | Primary Scope & Key Contents | Associated Deliverables |
|---|---|---|---|---|
| **Master Technical Report v2** | `docs/FULL_REPORT_v2.md` | v2.0 | Definitive 8-part scientific report documenting v1 audit, retractions, real-fuel corpus, evidence-driven model strategy, conformal uncertainty quantification, and hardened inverse designs. | All `results/v2/` CSVs |
| **Complete Progress Summary (v1 to v2.1)** | `docs/progress/v2/COMPLETE_PROGRESS_SUMMARY_V1_TO_V2.1.md` | v2.1 | Exhaustive chronological narrative and technical summary of all developments from initial v1 baseline to v2.1 completion. | `results/v2/improvement/` |
| **Thermodynamic Improvements Summary** | `docs/progress/v2/THERMODYNAMIC_IMPROVEMENTS_SUMMARY.md` | v2.1 | Traces the monotonic error reductions achieved by replacing empirical baselines with first-principles thermodynamics. | `paper_ablation_table.csv` |
| **Model Strategy Decision Note** | `docs/model_decision.md` | v2.0 | Formal evidence-driven decision document designating recalibrated physics as primary on real fuels and scoping neural residuals to lab binaries. | `results/v2/D1`–`D7` |
| **Laboratory Experimental Validation Plan** | `docs/experiment_plan.md` | v2.0 | Complete testing protocol for collaborating wet labs (500 mL batch recipes, purity grades, ASTM test methods, and 90% confidence envelopes). | `results/v2/E3/designs.csv` |
| **Full Report v1 (Historical)** | `FULL_REPORT.md` | v1.0 | Initial project capstone report describing v1 neural architecture, initial baselines, and early inverse designs (audited & revised in v2). | `results/t*` |
| **Progress Report v1 (Historical)** | `PROGRESS_REPORT.md` | v1.0 | Initial milestone progress summary of multi-task mixture models and grouped CV baselines. | `results/eval/` |
| **Project Handoff Reference** | `PROJECT_HANDOFF.md` | v1.0 | Early developer context file stating physics references (Refutas, Le Chatelier, PPBI), environment setup, and data schemas. | Repository root |

---

## 2. Governance, Audits & Claims Registers

Tracks the veracity of every scientific claim, data redistribution license, and physical constraint pass/fail status.

| Document Title | File Path | Generation | Primary Scope & Key Contents | Associated Deliverables |
|---|---|---|---|---|
| **Master Claims Register** | `docs/claims_register.md` | v2.0 | Live register tracking claims CLM-01 through CLM-11 with formal statuses (`OK`, `REVISED`, `RETRACTED`), evidence, and replacement findings. | `results/v2/A1`–`A11` |
| **Inverse Design Feasibility Audit** | `docs/progress/v2/A1.md` | v2.0 | Detailed physical audit proving why v1 Jet A blend (dodecane + tetralin) freezes at -15.2 °C and solidifies at -20 °C. | `results/v2/A1/design_audit.csv` |
| **Data Redistribution & Licensing Audit** | `docs/progress/v2/A12.md` | v2.0 | Legal audit across 27 sources; established factual data fair use doctrine for CRC tables and untracked raw AJFTD dumps. | `results/v2/A12/redistribution_audit.csv` |
| **Leakage & Overlap Audit** | `docs/audit_leakage.md` | v1.0 | T0.2 verification checking component overlap between training sets and real-fuel validation benchmarks. | `results/t0/leakage_checks.csv` |
| **Scope & Modeling Decisions** | `docs/scope_decisions.md` | v1.0 | Rationale for bounding project to pure hydrocarbons, closed-cup flash methods, and liquid-phase densities. | Project configs |
| **Autonomous Agent Questions Log** | `docs/agent_questions.md` | v1/v2 | Decision ledger resolving technical ambiguities during autonomous execution. | Project governance |

---

## 3. Post-v2 Thermodynamic Upgrades (v2.1)

Documentation of first-principles physical fixes, open dataset ingestions, and empirical error reductions.

| Document Title | File Path | Generation | Primary Scope & Key Contents | Associated Deliverables |
|---|---|---|---|---|
| **Flash Point Vapor Pressure Anchoring** | `docs/progress/v2/flash_anchored_upgrade.md` | v2.1 | Documents low-temperature vapor pressure anchoring for C10–C20 alkanes, eliminating the +42 °C flash over-prediction bug. | `flash_anchored_loso.csv`, `compounds_anchored_vp.csv` |
| **Coutinho Multi-Solid Wax Crystallization** | `docs/progress/v2/freezing_point_coutinho_upgrade.md` | v2.1 | Formulates fold-wise LOO refit of Coutinho eutectic wax depression, reducing freezing MAE to 3.99 °C (P=94.1% bootstrap vs flat offset). | `freezing_point_coutinho_eval.csv` |
| **ASTM D341 Walther Viscosity Curvature** | `docs/progress/v2/viscosity_walther_upgrade.md` | v2.1 | Upgrades cryogenic viscosity from linear 1000/T to ASTM D341 double-log curvature, cutting -40 °C log-RMSE by 75.1%. | `viscosity_walther_eval.csv` |
| **Cetane & Octane Model Recalibration** | `docs/progress/v2/recalibration_cetane_octane.md` | v2.1 | Documents model recalibration on 252 NREL cetane compounds (100% mixture coverage) and DOE Co-Optima blending values. | `cetane_nrel_recalibration.csv`, `octane_cooptima_recalibration.csv` |
| **Upgraded Inverse Blend Designs** | `docs/progress/v2/upgraded_inverse_designs.md` | v2.1 | Re-evaluates Jet A, Diesel, and Gasoline blends under upgraded physics engine; confirms 100% constraint feasibility. | `upgraded_designs.csv` |
| **NIST ThermoML Outlier Filtering** | `docs/progress/v2/thermoml_outlier_filtering.md` | v2.1 | Audits and filters 46,747 raw ThermoML records, dropping 2,172 unphysical rows and producing a clean liquid dataset. | `thermoml_clean_hc.csv` |

---

## 4. Phase-Level Progress Summaries (v2)

Comprehensive summaries generated at the completion of each major v2 phase.

| Document Title | File Path | Phase | Key Milestone Accomplished |
|---|---|---|---|
| **Phase B Real-Fuel Corpus Summary** | `docs/progress/v2/B6.md` | Phase B | Assembled 102 real fuels (45 measured, 57 reconstructed), 19 blend links, nested splits, and 12-property matrix. |
| **Phase E Inverse Design Summary** | `docs/progress/v2/PHASE_E_SUMMARY.md` | Phase E | Enforced all constraints in optimizer, mapped Pareto fronts, compared against 11 literature surrogates. |
| **Phase 4 Many-Component Complexity Summary** | `docs/progress/Phase4_many_components.md` | Phase 4 (v1) | Evaluated truncation, pseudo-composition reconstruction sensitivity, and error vs component count N. |
| **Phase 3 Real Fuels Validation Summary** | `docs/progress/Phase3_real_fuels.md` | Phase 3 (v1) | Initial integration of CRC 647 real aviation fuels and distillation curves. |
| **Phase 6 Inverse Design Initial Summary** | `docs/progress/Phase6_inverse_design.md` | Phase 6 (v1) | Early formulation of simplex optimizer and surrogate candidate palettes. |

---

## 5. Detailed Task-by-Task Progress Reports (v2)

Granular progress notes detailing mathematical equations, error tables, and unit test results for every individual v2 task.

### Phase A: Audit & Verification
- `docs/progress/v2/A1.md`: Independent audit of the three final inverse designs (Jet A, Diesel No. 2, Gasoline).
- `docs/progress/v2/A2.md`: Leave-one-fuel-out recalibration check of Refutas on 57 CRC fuels.
- `docs/progress/v2/A3.md`: Cetane basis bug check and resolution of `v_vol = w_mass` code bug.
- `docs/progress/v2/A4.md`: Full mixture Le Chatelier baseline re-evaluation (1,006 mixtures).
- `docs/progress/v2/A4_gn.md`: Grunberg-Nissan interaction check and comparison against Refutas.
- `docs/progress/v2/A5.md`: Distillation baseline evaluation and ASTM D86 vs D2887 comparison.
- `docs/progress/v2/A5_hoc.md`: Heat of combustion reinterpretation and LOO CV verification.
- `docs/progress/v2/A6.md`: Multi-property correlation audit and automated physical consistency rules.
- `docs/progress/v2/A6_uq.md`: Uncertainty growth factor provenance and real-fuel recalibration.
- `docs/progress/v2/A7.md`: Flash margin semantics and finite-sample conformal quantile analysis.
- `docs/progress/v2/A8.md`: Flash data usage audit across all 5,339 expanded observations.
- `docs/progress/v2/A9.md`: Deterministic hindcast reproduction and checkpoint provenance analysis.
- `docs/progress/v2/A12.md`: Data redistribution check and repository licensing audit.

### Phase B: Real-Fuel Validation Data
- `docs/progress/v2/B1.md`: Real-fuel candidate inventory across Jet, Gasoline, and Diesel families.
- `docs/progress/v2/B6.md`: Real-fuel property coverage matrix across 12 physical properties.

### Phase C: Property & Baseline Fixes
- `docs/progress/v2/C1.md`: Net heat of combustion group additivity baseline evaluation.
- `docs/progress/v2/C2.md`: Flash point UNIFAC activity coefficient baseline and LOSO evaluation.
- `docs/progress/v2/C3.md`: Flash expansion leakage safety verification.
- `docs/progress/v2/C4.md`: Cold flow SLE baseline and transparent LOO offset calibration.
- `docs/progress/v2/C5.md`: NREL compendium pure-component cetane digitization.
- `docs/progress/v2/C6.md`: Gasoline octane (RON/MON) hindcast on real fuels.
- `docs/progress/v2/C7.md`: Physics-based ASTM D86 distillation curve prediction.
- `docs/progress/v2/C8.md`: Gasoline dry vapor pressure equivalent (DVPE) baseline.

### Phase D: Evidence-Driven Model Strategy
- `docs/progress/v2/D1.md`: Real-fuel hindcast evaluation of the bilinear pairwise model.
- `docs/progress/v2/D2.md`: Class-level interaction learning and transfer evaluation on real fuels.
- `docs/progress/v2/D3.md`: Out-of-sample residual shrinkage sweep across temperatures.
- `docs/progress/v2/D4.md`: Stacked hierarchical model evaluation with L2 shrinkage priors.
- `docs/progress/v2/D5.md`: Prediction error vs number of blend components (The Empirical Inversion).
- `docs/progress/v2/D6.md`: Truncation and lumping sensitivity on measured-composition fuels.

---

## 6. Thematic Deep-Dive Technical Reports (`reports/`)

Comprehensive technical monographs compiled during v1 covering specific physical phenomena and model families.

| Report Title | File Path | Focus Area & Analysis |
|---|---|---|
| **Baseline Models Technical Monograph** | `reports/BASELINES.md` | In-depth formulations for Refutas, Kendall-Monroe, Le Chatelier, Amagat, and Redlich-Kister. |
| **Residual Models Architecture** | `reports/RESIDUAL_MODELS.md` | Mathematical architecture of excess thermodynamic friction and volume correction networks. |
| **Real Fuels Hindcast Analysis** | `reports/REAL_FUELS.md` | Analysis of distillation pseudo-composition mapping and bulk property matching on CRC fuels. |
| **Complexity and Truncation Study** | `reports/T3_COMPLEXITY.md` | Impact of mixture component count ($N=1$ to $N=40$) on model generalization and error bounds. |
| **Flash Mixture Thermodynamic Review** | `reports/T5_FLASH_MIXTURE.md` | Flammability limits, Antoine sets, and non-ideal vapor-liquid equilibrium in fuel mixtures. |
| **Experimental Noise Floor Analysis** | `reports/NOISE_FLOOR.md` | Quantifies experimental repeatability limits across ASTM standard test methods (D445, D56, D4809). |
| **Leenhouts Benchmark Reproduction** | `reports/LEE_REPRO.md` | Independent verification and reproduction of the Leenhouts et al. (2026) mixture dataset. |
| **Data Leakage & Train-Test Boundaries** | `reports/T1_LEAKAGE.md` | Formal graph-theoretic analysis of chemical compound sharing across validation splits. |
| **Unified Data Card** | `reports/DATA_CARD.md` | Comprehensive provenance, source licensing, and curation metadata for the combined dataset. |
| **Pure Component Models** | `reports/PURE_MODELS.md` | Benchmarking of Chew et al. pure-component viscosity and boiling point neural models. |

---

## 7. Architecture Specifications & Implementation Plans

Official implementation blueprints created under the `writing-plans` workflow guiding development:

- `docs/superpowers/plans/2026-11-04-task-spec-v2-master-plan.md`: Master architectural blueprint executing Phases A through F of Task Spec v2.
- `docs/superpowers/plans/2026-11-05-thermodynamic-error-reduction-plan.md`: Detailed plan executing the post-v2 thermodynamic physical upgrades (Tasks 1–5).

---

## 8. Summary Statistics

- **Total Cataloged Documents**: **98 technical Markdown files**
  - Master capstones & governance: 8
  - Post-v2 thermodynamic reports: 6
  - v2 Phase & Task progress reports: 37
  - v1 Task progress reports: 27
  - Thematic deep-dive monographs: 18
  - Architectural implementation plans: 2
- **Test Suite Verification**: **135 automated unit tests** passing in 7.92s (`pytest tests/test_v2*.py tests/test_result_regression.py`).
- **Release Version**: Tagged under **`v2.1-2026-11-05`** on `main` (`https://github.com/Kahaan83/UGRP_Hydrocarbons.git`).
