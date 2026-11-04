# SDD ledger - plan: C:\Users\kahaa\Downloads\Coding Agent Task Spec Property Coverage and Many-Component Prediction.md

Started: 2026-11-04. Execution: lead session, subagents for phase work per user request.

Task T0.1: complete (commits 241a08b, ddface6; tests: eval scripts re-run → PARTIAL PASS). Git init'd, tag baseline-2026-11-03, archive/baseline_manifest.json + baseline_2026-11-03_code.zip (5.4 MB). Repro: flash LOSO byte-identical, grouped CV and hindcast/base and §2.4/§3.3/§3.4 within tolerance; BUT §3.2 learned hindcast numbers in PROGRESS_REPORT.md do NOT reproduce from the frozen bundles (deepsets -20C log-RMSE claimed 0.1633, actual 0.2720; pairwise 0.1795 → 0.3400; deepsets density bias sign flip +3.64 → -0.23) -> headline "viscosity error cut in half" is not supported by tag artifacts. §3.1 extrapolation unchecked (would need destructive retrain). Details: results/t0_repro/comparison.csv, docs/progress/T0.1_repro_note.md.
Ruling (T0.1): model fallback — antigravity claude was quota-blocked (>90% usage), T0.2 assigned to antigravity-gemini-3.8-flash-high instead; cost: slightly weaker reasoning on the audit, mitigated by an existing 27 KB audit script left by the interrupted Sonnet run that the agent now verifies/completes.
Task T0.3: complete (commit ddface6; deliverable docs/data_sources.csv, 24 sources, unclear entries marked).
Task T0.4: complete (commit ddface6; deliverable data/antoine_sources.csv, 19 compounds — 11 cited, 8 flagged needs-removal-candidate, unverified provenance ACTIVE in flash baseline).
Task T0.2: complete (deliverables docs/audit_leakage.md, results/t0/leakage_*.csv; C1 PASS, C2 WARN physical coincidences, C3a WARN low-severity shared pure VP, C3b PASS zero blend leakage, C4 WARN 22/31 systems in S1, C5 PASS). Acceptance criteria met.
Phase 0 complete (T0.1, T0.2, T0.3, T0.4 all finished).

Phase 1 (Unified Data Schema):
Task T1.1: complete (deliverable data/schema/measurement_schema.md + src/fuelblend/data/validate.py validate_measurement_df; tests: 3 passed).
Task T1.2: complete (deliverable data/schema/composition_schema.md + src/fuelblend/data/validate.py validate_composition_df; tests: 3 passed).
Task T1.3: complete (deliverables data/compounds_canonical.csv, scripts/build_canonical_compounds.py, docs/progress/T1.3.md; 1038 unique compounds, 952 with CAS (91.7%), canonical RDKit SMILES and InChIKeys; tests: 1 passed).
Task T1.4: complete (deliverables data/processed/mixtures_tagged.csv, src/fuelblend/features/fuel_family.py, scripts/build_mixtures_tagged.py, docs/progress/T1.4.md; all 510 existing mixture systems tagged: 457 model-mixture, 38 gasoline, 15 diesel, 0 real_fuel; tests: 5 passed).
Task T1.5: complete (deliverables data/processed/measurements_migrated_core.csv, scripts/migrate_to_v2_schema.py, docs/progress/T1.5_migration_report.md; 35,338 rows before and after, exact value SHA-256 match de8c199824425c7e0422275c124c74bba212fa57027f71431962dd8272b481c2; strict schema validator passed; tests: 1 passed).
Phase 1 complete (all 13 tests green).

Phase 5 (Statistics and Recalibration Ablation):
Task T5.1: complete (deliverables results/t5/system_cluster_bootstrap_summary.csv, scripts/run_cluster_bootstrap_eval.py, docs/progress/T5.1.md, src/fuelblend/eval/eval_grouped_cv.py, src/fuelblend/eval/eval_flash_loso.py, src/fuelblend/mixture/fuel_hindcast.py; all CIs explicitly cluster-bootstrapped over system_id/fuel_id; tests: 1 passed).
Task T5.2: complete (deliverables results/t5/recalibrated_physics_ablation.csv, scripts/recalibrated_physics_ablation.py, docs/progress/T5.2.md; proved 2-parameter 1000/T slope recalibration recovers 68.5% of error reduction on cold-flow viscosity at -40C/-20C, outpacing complex neural residuals; tests: 1 passed).

Phase 2 (Data Collection by Property — Section 3.1 Heat of Combustion):
Task T2.1: complete (deliverables data/processed/gcv_MJ_kg.csv, scripts/check_gcv_coverage.py, docs/progress/T2.1.md; 100% pure compound coverage across 1038 compounds: 114 Tier A stoichiometric, 924 Tier C modified Dulong).
Task T2.2: complete (deliverables src/fuelblend/baselines/heat_of_combustion.py, tests/test_heat_of_combustion_baseline.py; mass vs mole basis conversions, gross vs net conversions, water latent heat accounting; tests: 3 passed).
Task T2.3: complete (deliverables data/processed/measurements_real_fuels_heat_of_combustion.csv, scripts/build_heat_of_combustion_measurements.py, docs/progress/T2.3.md; 1,083 real-fuel heat of combustion measurements from AJFTD and CRC 647 with test methods ASTM D4809 / D3338, schema validated).
Task T2.4: complete (deliverables results/t2/heat_of_combustion_eval.csv, scripts/eval_heat_of_combustion_baseline.py, docs/progress/T2.4.md, tests/test_heat_of_combustion_eval.py; evaluated across 57 real aviation fuels with cluster bootstrap CIs over fuel_id; physics baseline achieves 1.60% MAPE and 0.689 MJ/kg MAE).

Phase 2 (Data Collection by Property — Section 3.4 Cetane / DCN):
Task T2.13: complete (deliverables data/processed/cetane_pure_curated.csv, scripts/build_dcn_pure_compounds.py, docs/progress/T2.13.md; 80 pure hydrocarbons with standardized method ASTM D6890/D613 and reliability flags; remaining un-digitized NREL compounds logged in data/MISSING.md).
Task T2.14: complete (deliverables data/processed/measurements_dcn.csv, scripts/build_dcn_measurements.py, docs/progress/T2.14.md; 391 DCN measurements collected across 346 speciated mixtures and 45 real fuels from AJFTD, strict schema validated).
Task T2.15: complete (deliverables src/fuelblend/baselines/cetane.py, results/t2/cetane_baseline_grouped_cv.csv, scripts/eval_cetane_baseline_grouped_cv.py, docs/progress/T2.15.md, tests/test_cetane_baseline.py; mole vs mass vs volume basis tested under grouped CV across 31 systems; mass/volume basis proved optimal with MAE 3.731 vs mole 5.707).

Phase 2 (Data Collection by Property — Section 3.3 Low-Temperature Behavior):
Task T2.9: complete (deliverables data/processed/measurements_low_temp_flow.csv, scripts/build_low_temp_measurements.py, docs/progress/T2.9.md; 189 low-temp flow property measurements from CRC 647, AJFTD, NASA, Kirgina with test methods ASTM D5972/D2386/D5949/D2500, schema validated).
Task T2.10: complete (deliverables scripts/build_sle_cold_flow_coverage.py, docs/progress/T2.10.md; documented 905 pure hydrocarbons with experimental Tm and SLE parameters for n-alkanes C8-C22 and solid aromatics).
Task T2.11: complete (deliverables src/fuelblend/baselines/cold_flow.py, results/t2/cold_flow_baseline_eval.csv, scripts/eval_wat_baseline.py, docs/progress/T2.11.md, tests/test_cold_flow_baseline.py; implemented Schroeder-van Laar ideal SLE baseline for wax appearance temperature, evaluated across 57 CRC aviation fuels with r = 0.443 correlation to ASTM D5972 freezing points and cluster-bootstrap CIs).
Task T2.12: complete (deliverables docs/scope_decisions.md Decision SD-01, docs/progress/T2.12.md; audit confirmed 0 speciated molecular mixtures with measured pour point < 30 threshold; formally recorded decision to scope pour point out in favor of thermodynamic freezing point and WAT).

Phase 3 (Real-Fuel Dataset with Composition):
Task T3.1: complete (deliverable data/real_fuels/inventory.csv; 65 real fuels catalogued with source, fuel family, component count, unresolved fraction, and available properties).
Task T3.2: complete (deliverable data/real_fuels/compositions.csv; 2,871 component fraction rows linked to canonical InChIKeys and SMILES; unresolved fraction tracked and validated by strict composition validator).
Task T3.3: complete (deliverable data/real_fuels/blend_links.csv; links RP-1 surrogates and commercial blends with blend ratios and measured properties).
Task T3.4: complete (CRC 647 pseudo-compositions segregated into separate subset labeled reconstructed_composition and isolated from measured_dha_gcxgc fuels).
Task T3.5: complete (deliverable splits/real_fuels_holdout.json; holds out all 65 real fuels from training by default with nested split definitions; Phase 3 report docs/progress/Phase3_real_fuels.md; tests: 4 passed).

Phase 4 (Models for Up to 40+ Components):
Task T4.1: complete (deliverable results/t4/pair_coverage.csv, docs/progress/T4.1.md; pair-coverage diagnostic showed molecular pair coverage drops to 0.57% on real fuels, while class pair coverage remains >95%).
Task T4.2: complete (deliverable src/fuelblend/mixture/bilinear_pairwise.py, tests/test_bilinear_pairwise.py; implemented low-rank bilinear interaction model scaling in O(N*d*r) time; exact match to explicit O(N^2) double sum verified; tests: 4 passed).
Task T4.3: complete (deliverable ClassSharedInteraction in bilinear_pairwise.py; symmetric 5x5 class-level interaction matrix).
Task T4.4: complete (deliverable tests/test_model_invariances.py, docs/progress/T4.4.md; 5 physical invariance principles verified; documented why DeepSets attention fails component splitting while Pairwise passes; tests: 5 passed).
Task T4.5: complete (deliverable results/t4/truncation_study.csv; proved retaining top 20 components captures 96.8% mass and 98.4% variance).
Task T4.6: complete (deliverable results/t4/reconstruction_sensitivity.csv; measured prediction spread across pseudo-composition perturbations).
Task T4.7: complete (deliverable results/t4/error_vs_n_components.csv; error stratified by component count bucket with cluster bootstrap CIs).
Task T4.8: complete (deliverable results/t4/uncertainty_growth.csv; fitted property-specific inflation factors g_visc=1.377, g_dens=2.530, g_flash=1.541, validating >90% empirical coverage on real fuels).
Task T4.9: complete (deliverable results/t4/physics_first_regularization.csv, docs/progress/Phase4_many_components.md; evaluated L2 residual shrinkage ensuring physics dominance at high N; tests: 3 passed).

Phase 2 (Data Collection by Property — Section 3.2 Flash Point Expansion):
Task T2.5: complete (deliverable data/processed/measurements_flash_expanded.csv; 5,339 flash measurements standardized with test method ASTM D56/D93/D3828/IP170 filled for 100% of rows; strict schema validated).
Task T2.6: complete (deliverable docs/progress/T2.6.md; 100% pure compound coverage catalogued for Antoine saturation vapor pressure and LFL parameters).
Task T2.7: complete (deliverables results/t2/flash_loso_activity_coeff.csv, scripts/build_flash_expanded.py; implemented Scatchard-Hildebrand regular solution activity-coefficient baseline; reduced RMSE from 10.24 C to 8.92 C and dangerous over-prediction rate from 19.8% to 16.7% across 31 flash LOSO systems).
Task T2.8: complete (deliverable results/t2/flash_loso_by_method.csv, docs/progress/T2.5_T2.8.md, tests/test_flash_expansion.py; added method feature and evaluated results separately by standardized test apparatus; tests: 3 passed).

Phase 2 (Data Collection by Property — Section 3.5 Octane & Section 3.6 Consistency Checks):
Task T2.16: complete (deliverable data/processed/measurements_octane.csv; 1,146 octane measurements collected across 685 pure, 292 measured_mixture, and 169 blending_value records; strict schema validated).
Task T2.17: complete (deliverables src/fuelblend/baselines/octane.py, results/t2/octane_baseline_grouped_cv.csv; linear-by-volume and interaction-based BON evaluated under grouped CV; interaction BON improved MAE to 2.479 and RMSE to 3.808).
Task T2.18: complete (deliverable results/t2/octane_sensitivity_report.csv; modeled derived octane sensitivity RON - MON across multicomponent FACE gasolines, mean sensitivity = 5.14).
Task T2.19: complete (deliverable docs/octane_data_notes.md; documented CFR engine test repeatability limits and inferred vs experimental data quality flags).
Task T2.20: complete (deliverable data/processed/measurements_distillation.csv; 342 ASTM D86 distillation curve points across 57 CRC aviation fuels, strict schema validated).
Task T2.21: complete (deliverables results/t2/distillation_consistency_check.csv, docs/progress/T2.20_T2.21.md, tests/test_octane_phase2.py, tests/test_distillation_phase2.py; distillation slope and boiling point consistency checks implemented; tests: 6 passed).
Phase 2 complete across all 6 property sections (T2.1-T2.21).

Phase 5 (Part 2 — Literature Baselines, Harder Splits & Multi-Seed Runs):
Task T5.3: complete (deliverables src/fuelblend/baselines/literature_baselines.py, results/t5/literature_baselines_eval.csv; Grunberg-Nissan, Kendall-Monroe, Refutas, and Redlich-Kister excess volume evaluated on identical splits).
Task T5.4: complete (tuned HistGradientBoosting and Bilinear Pairwise ML baselines benchmarked).
Task T5.5: complete (deliverables splits/harder_splits.json, results/t5/harder_splits_eval.csv; constructed and evaluated 4 strict holdout splits: unseen component pairs, unseen N components, and unseen temperature ranges).
Task T5.6: complete (deliverable results/t5/multiseed_ensemble_eval.csv, docs/progress/T5.3_T5.6.md, tests/test_baselines_harder_splits_phase5.py; executed 5 independent seeds with mean and standard deviation spread reported; tests: 3 passed).
Phase 5 complete (all tasks T5.1-T5.6 finished).

Phase 6 (Inverse Design for Three Fuel Families):
Task T6.1: complete (deliverables configs/specs_jet.yaml, configs/specs_diesel.yaml, configs/specs_gasoline.yaml; explicit property targets, ASTM standards, source citations, and property-specific uncertainty factors).
Task T6.2: complete (deliverable results/t6/flash_margin_sweep.csv; evaluated flash LOSO residual quantile sweep; principled 90% conformal quantile is 3.4 C, replacing arbitrary 5 C heuristic).
Task T6.3: complete (engine dynamically reads property-specific uncertainty inflation factors g_visc=1.377, g_dens=2.530, g_flash=1.541 from YAML configs).
Task T6.4: complete (deliverable results/t6/inverse_designs_three_families.csv; generated verified feasible multi-fuel formulations for Jet, Diesel, and Gasoline with per-constraint feasibility checks).
Task T6.5: complete (cardinality penalty restricts blend to 2-4 dominant components even when optimized over candidate pool of 45 hydrocarbons).
Task T6.6: complete (deliverable results/t6/inverse_design_ablations.csv, docs/progress/Phase6_inverse_design.md, tests/test_inverse_design_phase6.py; evaluated physics-only vs learned+UQ, sparsity penalty, and pool size; tests: 5 passed).
Entire task spec (Phases 0 through 6) complete. All 65 tests passing.










