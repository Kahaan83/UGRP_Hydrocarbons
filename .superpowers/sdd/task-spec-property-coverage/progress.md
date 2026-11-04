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



