# Coding Agent Task Spec v2: Verification, Real-Fuel Validation and Model Strategy Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Execute the full v2 task specification: Phase A (audit & triage A1–A12), Phase B (real-fuel validation datasets B1–B6), Phase C (property & baseline fixes C1–C8), Phase D (evidence-driven model strategy D1–D7), Phase E (inverse design hardening E1–E6), and Phase F (reproducibility freeze & report F1–F4).

**Architecture:** Maintain strict ground rules: no overwriting v1 outputs (new outputs to `results/v2/<phase>/<task_id>/` and progress docs to `docs/progress/v2/<task_id>.md`), track claims in `docs/claims_register.md` with explicit `RETRACTED` or `REVISED` designations, enforce leave-fuel-out / leave-system-out cluster bootstrapping, and verify all findings with independent re-computation scripts.

**Tech Stack:** Python 3.12 (in `.venv`), PyTorch, scikit-learn, RDKit, pandas, scipy, pytest.

**Spec:** `C:\Users\kahaa\Downloads\Coding Agent Task Spec v2 Verification, Real-Fuel Validation and Model Strategy.md`

## Global Constraints
1. No overwriting v1 outputs. New outputs go under `results/v2/<phase>/<task_id>/` and notes under `docs/progress/v2/<task_id>.md`. Keep tag `baseline-2026-11-03` intact.
2. Retraction policy: Any earlier claim that fails a v2 check is marked `RETRACTED` or `REVISED` in `docs/claims_register.md` with the reason and replacing result.
3. Leave-fuel-out for anything fit on real fuels. Never score a fit on the fuels it was fit on.
4. Cluster bootstrap by `system_id` (mixtures, 4,000 resamples) or `fuel_id` (real fuels, 2,000 resamples).
5. Every model comparison includes raw physics, recalibrated physics, and learned model on identical splits with at least 5 seeds for stochastic models.
6. Independent re-computation for any number in a conclusion.
7. Do not fabricate data or constants. Missing values go in `data/MISSING.md`.
8. Report negative results plainly.
9. Log ambiguity in `docs/agent_questions.md` and continue with safest assumption.

## Review Focus
- Jet A inverse design feasibility: verify if hexadecane + biphenyl or dodecane + tetralin physically pass freezing point <= -40 C, density 775-840 kg/m3, aromatics <= 25%, flash point >= 38 C.
- Recalibrated Refutas LOO vs lab-only fitting: verify whether 2-parameter slope recalibration holds up when fit strictly with leave-one-fuel-out on CRC fuels and when trained purely on lab mixtures.
- Cetane mass vs volume basis: determine whether mass vs volume identical error (MAE 3.731) was a bug where volume basis was using mass fractions or identical densities.
- Grunberg-Nissan interaction parameter: confirm if GN baseline had interaction parameter $G_{ij} = 0$ or fitted interaction.
- Hindcast determinism: identify why v1 bundle gave 0.2720 vs claimed 0.1633 and fix reproducibility via fixed seeds and single command script.

---

### Phase A: Triage and Verification Tasks (A1 - A12)

#### Task A1: Independent Audit of the Three Final Inverse Designs
**Files:**
- Create: `scripts/v2/audit_inverse_designs.py`
- Deliverable: `results/v2/A1/design_audit.csv`
- Note: `docs/progress/v2/A1.md`
- Test: `tests/test_v2_a1_audit.py`

- [ ] Step 1: Write independent verification script to recompute properties (density, flash point via Le Chatelier + Antoine, kinematic viscosity at -20 C via Refutas, freezing point / WAT via SLE, aromatics %, net heat of combustion) for the three v1 designs in `results/t6/inverse_designs_three_families.csv` and the text-claimed dodecane + tetralin blend.
- [ ] Step 2: Test constraint pass/fail against `configs/specs_jet.yaml`, `configs/specs_diesel.yaml`, `configs/specs_gasoline.yaml`.
- [ ] Step 3: Write `results/v2/A1/design_audit.csv` with full constraint breakdown (value, limit, margin, pass/fail).
- [ ] Step 4: Write `docs/progress/v2/A1.md` analyzing why dodecane + tetralin and hexadecane + biphenyl fail or pass each spec constraint (specifically freezing point <= -40 C and aromatics limit).
- [ ] Step 5: Run tests, commit, push.

#### Task A2: Leave-Fuel-Out Recalibration Check
**Files:**
- Create: `scripts/v2/eval_recalibration_loo.py`
- Deliverable: `results/v2/A2/recalibration_loo.csv`
- Note: `docs/progress/v2/A2.md`
- Test: `tests/test_v2_a2_loo.py`

- [ ] Step 1: Implement strict leave-one-fuel-out (LOO) cross-validation for the 2-parameter (1000/T slope) and 3-parameter (per-T bias) Refutas recalibrations across the 57 CRC aviation fuels.
- [ ] Step 2: Implement lab-only calibration: fit the 2-parameter recalibration exclusively on laboratory mixture data (`model_table.csv` / `measurements_migrated_core.csv`) and evaluate zero-shot on the 57 CRC real fuels.
- [ ] Step 3: Compute cluster-bootstrapped (2,000 resamples over `fuel_id`) log-RMSE and MAE for raw Refutas, LOO recalibrated Refutas, lab-only recalibrated Refutas, DeepSets, and Pairwise models.
- [ ] Step 4: Save `results/v2/A2/recalibration_loo.csv` and document in `docs/progress/v2/A2.md`.
- [ ] Step 5: Run tests, commit, push.

#### Task A3: Cetane Basis Bug Check
**Files:**
- Create: `scripts/v2/check_cetane_basis.py`
- Modify: `src/fuelblend/baselines/cetane.py`
- Note: `docs/progress/v2/A3.md`
- Deliverable: `results/v2/A3/cetane_basis_rerun.csv`
- Test: `tests/test_v2_a3_cetane.py`

- [ ] Step 1: Inspect `src/fuelblend/baselines/cetane.py` and `scripts/eval_cetane_baseline_grouped_cv.py` to identify why mass and volume bases produced identical MAE (3.731).
- [ ] Step 2: Fix any bug where volume fractions defaulted to mass fractions or densities were omitted/constant.
- [ ] Step 3: Re-evaluate mass vs volume vs mole basis under grouped CV, verifying that mass and volume now differ physically.
- [ ] Step 4: Save results and write `docs/progress/v2/A3.md`.
- [ ] Step 5: Run tests, commit, push.

#### Task A4: Grunberg-Nissan Interaction Check
**Files:**
- Create: `scripts/v2/eval_gn_baselines.py`
- Modify: `src/fuelblend/baselines/literature_baselines.py`
- Deliverable: `results/v2/A4/gn_baselines.csv`
- Note: `docs/progress/v2/A4.md`
- Test: `tests/test_v2_a4_gn.py`

- [ ] Step 1: Inspect `literature_baselines.py` for Grunberg-Nissan implementation.
- [ ] Step 2: Report GN with $G_{ij} = 0$ (zero interaction, ideal mixing of $\ln \eta$) vs GN with fitted interaction parameter (fitted per system or class pair under strict grouped CV).
- [ ] Step 3: Compare log-RMSE against Refutas (0.0698) on identical splits.
- [ ] Step 4: Write `results/v2/A4/gn_baselines.csv` and `docs/progress/v2/A4.md`.
- [ ] Step 5: Run tests, commit, push.

#### Task A5: Heat of Combustion Reinterpretation
**Files:**
- Create: `scripts/v2/recheck_hoc.py`
- Deliverable: `results/v2/A5/hoc_recheck.csv`
- Note: `docs/progress/v2/A5.md`
- Test: `tests/test_v2_a5_hoc.py`

- [ ] Step 1: Re-evaluate the Ridge residual model for net heat of combustion on the 57 CRC real fuels using strict leave-one-fuel-out cross-validation.
- [ ] Step 2: Document features used by the Ridge model and compute cluster-bootstrapped CIs.
- [ ] Step 3: Contrast Ridge MAE (0.078 MJ/kg) vs physics MAE (0.689 MJ/kg) and provide scientific reinterpretation.
- [ ] Step 4: Save `results/v2/A5/hoc_recheck.csv` and `docs/progress/v2/A5.md`.
- [ ] Step 5: Run tests, commit, push.

#### Task A6: Uncertainty Growth Factor Provenance and Refit
**Files:**
- Create: `scripts/v2/refit_uncertainty_growth.py`
- Deliverable: `results/v2/A6/uncertainty_growth_refit.csv`
- Note: `docs/progress/v2/A6.md`
- Test: `tests/test_v2_a6_uq.py`

- [ ] Step 1: Trace original derivations of $g_{visc}=1.377$, $g_{dens}=2.530$, $g_{flash}=1.541$.
- [ ] Step 2: Refit inflation factors on real fuels using leave-one-fuel-out cross-validation.
- [ ] Step 3: Measure empirical 90% coverage on held-out fuels only.
- [ ] Step 4: Save `results/v2/A6/uncertainty_growth_refit.csv` and document in `docs/progress/v2/A6.md`.
- [ ] Step 5: Run tests, commit, push.

#### Task A7: Flash Margin Semantics and Quantile Analysis
**Files:**
- Create: `scripts/v2/eval_flash_margin_semantics.py`
- Deliverable: `results/v2/A7/flash_margin_semantics.csv`
- Note: `docs/progress/v2/A7.md`
- Test: `tests/test_v2_a7_flash.py`

- [ ] Step 1: Analyze the 3.4 °C flash safety margin over the 31 LOSO systems.
- [ ] Step 2: Calculate empirical over-prediction coverage, worst-case over-prediction (+11.4 °C in v1), and margin dependence on system count $N$.
- [ ] Step 3: Document statistical assumptions and limitations in `docs/progress/v2/A7.md`.
- [ ] Step 4: Run tests, commit, push.

#### Task A8: Flash Data Usage Audit
**Files:**
- Create: `scripts/v2/audit_flash_data.py`
- Deliverable: `docs/progress/v2/A8.md`, `results/v2/A8/flash_usage_audit.csv`
- Test: `tests/test_v2_a8_flash_audit.py`

- [ ] Step 1: Audit the 5,339 expanded flash measurements (`measurements_flash_expanded.csv`).
- [ ] Step 2: Break down counts by pure vs mixture, apparatus/method, and data source.
- [ ] Step 3: Detail which rows fed baseline parameter fitting vs neural training, and explain why LOSO was restricted to 126 mixture rows.
- [ ] Step 4: Commit and push.

#### Task A9: Deterministic Hindcast Reproduction
**Files:**
- Create: `scripts/reproduce_hindcast.sh`, `scripts/v2/reproduce_hindcast.py`
- Deliverable: `docs/progress/v2/A9.md`, `results/v2/A9/hindcast_reproduced.csv`
- Test: `tests/test_v2_a9_hindcast.py`

- [ ] Step 1: Investigate discrepancy between checkpoint 0.2720 and reported 0.1633.
- [ ] Step 2: Build a deterministic execution pipeline with fixed random seeds and documented checkpoint selection rules.
- [ ] Step 3: Verify two consecutive executions produce bit-identical results.
- [ ] Step 4: Commit and push.

#### Task A10: Result Regression Test Suite
**Files:**
- Create: `tests/test_result_regression.py`

- [ ] Step 1: Add regression tests pinning headline numbers (grouped CV, flash LOSO, hindcast, design audit) within specified numerical tolerances.
- [ ] Step 2: Verify test fails if numbers drift without claims register update.
- [ ] Step 3: Commit and push.

#### Task A11: Master Claims Register
**Files:**
- Create: `docs/claims_register.md`

- [ ] Step 1: Catalogue every claim made in `FULL_REPORT.md` and v1 progress notes.
- [ ] Step 2: Assign status (`OK`, `REVISED`, or `RETRACTED`) with exact reasons, replacement findings, and tracing files.
- [ ] Step 3: Commit and push.

#### Task A12: Data Redistribution and Licensing Audit
**Files:**
- Create: `scripts/v2/audit_data_redistribution.py`
- Deliverable: `docs/progress/v2/A12.md`, `results/v2/A12/redistribution_audit.csv`
- Test: `tests/test_v2_a12_licenses.py`

- [ ] Step 1: Cross-reference `docs/data_sources.csv` licensing flags against repository git-tracked files.
- [ ] Step 2: Flag restricted data and prescribe actions (scrub, replace with rebuild script, or retain with license).
- [ ] Step 3: Commit and push.
