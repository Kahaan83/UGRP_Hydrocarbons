# v4 — Data Integrity Repair, Re-evaluation and Upgrades: Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Remove fabricated and circular data from the fuel-blend pipeline, rebuild the real-fuel corpus from traceable sources, re-run every result built on it, then add three new properties (smoke point/YSI, SLE freezing, real-blend octane) and the model upgrades the clean data supports.

**Architecture:** New work lives in `scripts/v4/`, `src/fuelblend/` (new modules only where named), `results/v4/<task>/`, `data/real_fuels/v4/`, `data/quarantine/`, and `tests/test_v4_*.py`. Nothing in `results/v2`, `results/v3` or old docs is edited or deleted. Quarantined items get moved or marked, and retractions are recorded in `docs/claims_register.md` (Ground Rule 2). Every number a model or margin uses has to trace to a raw file or a cited document. Otherwise the value is NaN and gets logged in `data/MISSING.md`.

**Tech Stack:** Python 3.12 (`.venv\Scripts\python.exe`), pandas, numpy, scipy, scikit-learn, torch, RDKit, `thermo`/`chemicals` (UNIFAC, DIPPR-style correlations), pytest, PyYAML.

**Spec:** this document plus `docs/PROJECT_HANDOFF.md` (state), `AGENTS.md` §2/§6/§9/§13 (physics, schema, fixed bugs, conventions), `docs/model_decision.md`, `docs/claims_register.md`.

## Why this plan exists (findings from the 2026-11-09 review — verify each in Task 1)

1. **Fabricated "measured" fuel compositions.** In `scripts/v3/ingest_phase1_expansion.py` (and the same code in `scripts/v2/ingest_real_fuels_v2.py` ~lines 300–330), the FACE gasoline compositions are computed from each fuel's own RON/MON/density (`aro_frac=(dens-680)/200`, `ole_frac=(RON-MON)/20`, `iso_frac∝RON`, ~lines 216–223). They are still labelled `composition_type="measured_dha_gcxgc"` with "78 resolved components". The octane model is then scored on predicting RON from these compositions, which is circular.
2. **Formula-derived targets.** Ethanol-blend DVPE is `base + 7.0 (≤E20) else +3.0`. Every FACE heat of combustion is the constant `43.5`. FACE diesel flash is `52 + (ρ−800)·0.25`. Diesel compositions are back-solved from DCN and density. Freezing is `cloud − 3`, pour is `freeze − 4`, NASA flash is a constant `58.0`, and the RP-1/POSF densities are `804.0 if "RP1" else 795.0`. All five POSF additional fuels get the same 48/32/20 split.
3. **Wrong identities.** The NASA paraffinic-diesel rows pair InChIKeys with the wrong SMILES. `DCAYPVUWAIABOU` (hexadecane) is used for decane, `NHTMVDHEPJAVLT` (isooctane) for tridecane, and `IMNFDUFMRHMDMM` (heptane) for hexadecane.
4. **The v3 design optimizer bypasses the library physics** (`scripts/v3/run_v3_inverse_designs.py`):
   - Flash uses `1/T = Σ x_i/T_flash,i` (~89–103), which is not Le Chatelier.
   - "Coutinho" freeze is an invented `0.60 + 3.20σ + 14.50·(non-paraffin fraction)` offset (~105–147). That offset is why the diesel designs show pour ≈ −22 °C; ideal SLE gives a WAT of about −5 °C for 20 wt % n-C16.
   - DVPE uses a hand-typed Psat dict and its own γ (~165–203).
   - The margins (~38–49) disagree with `configs/uncertainty_v3.yaml` (density 8.0 vs 2.5, DVPE 5.0 vs 7.2, RON 3.5 vs 4.5).
   - The "5 model-form robustness, 30/30 pass" check (~522–545) only adds constant biases (+0.2 °C, ×1.02) to one prediction.
   - There is no DVPE minimum, so gasoline designs at 17.8–19.4 kPa pass. Real gasoline volatility specs would reject these. Check the actual limits in IS 2796 and the ASTM D4814 classes.
5. **Internal inconsistency.** The density MAE on real fuels is 12.94 kg/m³ (`results/v3/G10`), but the v3 design margin is ±2.50 kg/m³.
6. **Real data already on disk but unused:**
   - `extras/fuelblend_data/data/raw/octane/face_gasoline_set_{I,II}.csv` (measured PRF/TPRF/FACE+ethanol RON/MON)
   - `raw/cooptima/{bron_data_with_cas.xls,bob_data.xlsx,prop_db_amr.xls,fuel_engine_co_optimization.csv}`
   - `raw/gc_composition/fuel_gc_compositions.csv` (cited class-level jet compositions)
   - `raw/sooting/*YSI*`, `processed/{tsi_table.csv,smoke_point_ajftd.csv,sooting_real_fuels.csv}`
   - `raw/cetane/nrel_cn_merged_210831.csv`
   - ThermoML SLE data, which `thermoml_extract.py` drops
7. **Hygiene.** `requirements.txt` is missing scikit-learn, torch, joblib, pytest and PyYAML (all imported). The bilinear model's init is unseeded (CLM-19). Open items F1–F6 remain (PROJECT_HANDOFF §4 Stage 5).

**Trustworthy and kept as-is (re-verify only):** the lab pipeline (Leenhouts, NIST/TDE, cleaned ThermoML), the 126-mixture flash LOSO (CLM-12), the CRC 647 *measured properties* (ingested from the PDF), and the CRC 647 reconstructed compositions. The reconstruction is derived from the measured aromatics and D86; it is legitimate but must keep its `reconstructed` label.

## Global Constraints

- **No fabricated values.** A property value or composition enters `data/real_fuels/v4/` only with a `provenance` that points at a raw file path + row, or at a document + table/page. Otherwise it is NaN and gets a row in `data/MISSING.md`.
- **No circularity.** A fuel's composition may not be derived from any property the same fuel is evaluated on.
- Never delete or silently edit v2/v3 results or reports; retractions go in `docs/claims_register.md` as new rows CLM-20 onward.
- Splits are grouped by `system_id` (lab) or `fuel_id` (real). Any hyperparameter or model choice is selected inside the outer fold (nested CV).
- Units are converted only in loaders. Antoine is `log10(P/kPa)=A−B/(T/K+C)`. Viscosity is modelled in ln space. Flash is reported with a one-sided over-prediction metric.
- Default to hydrocarbon-only (`is_hc`); oxygenates are allowed only in the gasoline octane/DVPE tasks.
- Every new loader gets a test whose fixture mirrors the **real** file's header (copy the first rows of the real file into `tests/data/`).
- Never commit `ajftd_session.json`, raw AJFTD pages, or licensed raw data.
- Work on branch `v4-integrity`; commit after every task; run with `.venv\Scripts\python.exe`.
- Training time is not a constraint: overnight laptop runs are fine, and GPU jobs go in a job pack the user runs on the CS lab machine.
- The user may not be able to get the paywalled documents in Task 8. Every task must run and report without them.

## Review Focus

1. **Circular derivation that hides behind a label** (composition computed from a target, then that target predicted). The ledger must detect it from code and data, not trust `composition_type`. Test: Task 1 `test_ledger_flags_face_composition_as_target_derived`.
2. **SMILES/InChIKey mismatch** in any composition table. Expect it flagged, with the row excluded from v4. Test: Task 2 `test_nasa_diesel_identity_mismatch_detected`.
3. **Too few groups for a conformal margin** (n_groups < 10, e.g. 6 D56 fuels). The margin must come back `inf` with a warning, not a small optimistic number. Test: Task 10 `test_margin_infinite_when_too_few_groups`.
4. **Design optimum outside the model's validated domain** (component mass fraction, temperature or class mix outside the training/validation envelope). Expect the design marked `extrapolating=True` and excluded from PASS. Test: Task 17 `test_design_flags_extrapolation`.
5. **Scale/unit mix-ups in new sources** (YSI scales, mole vs volume octane blending, kPa vs bar vs psi). Expect loaders to normalise and tests to pin one known value per source. Tests: Task 14 `test_ysi_unified_scale_anchor`, Task 15 `test_bron_units_and_basis`.

---

## Phase 0 — Setup

### Task 0: Branch, dependencies, determinism

**Files:** Modify `requirements.txt`, `src/fuelblend/mixture/bilinear_pairwise.py:28-30`, `Makefile`. Test `tests/test_v4_setup.py`.

- [ ] Create branch `v4-integrity` from `main`.
- [ ] Run `.\.venv\Scripts\python.exe -m pytest tests -q` and save the summary to `results/v4/baseline_pytest.txt` (expected 290 passed per FULL_REPORT_v3 §11; record what you actually get).
- [ ] Add pinned lower bounds to `requirements.txt`: `scikit-learn>=1.3`, `torch>=2.1`, `joblib>=1.3`, `pytest>=7.4`, `PyYAML>=6.0`.
- [ ] Write `test_bilinear_init_is_deterministic`: two `BilinearPairwiseInteraction(...)` constructions with `seed=0` have identical parameters.
- [ ] Add a `seed: int = 0` constructor argument that builds a local `torch.Generator` for the `randn` init. Re-run `scripts/v3/build_viscosity_leaderboard.py` twice and confirm the CSVs are byte-identical. Record the CLM-19 fix in the claims register.
- [ ] Add a `make v4` target that runs `scripts/v4/reproduce_all.py` (Task 19).
- [ ] Commit.

## Phase 1 — Integrity audit and quarantine

### Task 1: Provenance ledger

**Files:** Create `scripts/v4/audit_provenance.py`, `src/fuelblend/data/provenance.py`. Output `results/v4/P1/provenance_ledger.csv`, `results/v4/P1/AUDIT.md`. Test `tests/test_v4_provenance.py`.

**Interfaces — Produces:**
- `classify_fuel_rows(inventory: pd.DataFrame, compositions: pd.DataFrame, source_code: dict[str, str]) -> pd.DataFrame` with columns `fuel_id, field, value, provenance_class, evidence`.
- `provenance_class ∈ {raw_file, cited_typed, formula_derived, target_derived_composition, placeholder_constant, unknown}`.

The ledger covers every row of `data/real_fuels/inventory*.csv`, `compositions*.csv` and `data/real_fuels/v3/*`.

- [ ] Write the tests. Each one runs the classifier on the real current files and asserts it detects the finding:
  - `test_ledger_flags_face_composition_as_target_derived`: every `CRC_FACE_[A-J]` composition is `target_derived_composition`.
  - `test_ledger_flags_ethanol_dvpe_formula`: `CRC_FACE_A_E10.dvpe_kpa` is `formula_derived`.
  - `test_ledger_flags_constant_hoc`: FACE `heat_of_combustion_MJ_kg == 43.5` is `placeholder_constant`.
  - `test_ledger_flags_diesel_flash_formula`: `CRC_FACE_Diesel_1.flash_point_C` is `formula_derived`.
  - `test_ledger_flags_posf_48_32_20`.
  - `test_ledger_keeps_crc647_measured`: CRC647 `freezing_point_C` is `raw_file`, with evidence pointing at `crc647_fuel_properties.csv`.
- [ ] Implement detection in two ways:
  - (a) Static: AST-scan the two ingest scripts for literal tuples and arithmetic assigned to property or fraction keys.
  - (b) Data: test for exact functional relations, e.g. `aro_frac == clip((dens-680)/200)`, `freeze == cloud-3`, `pour == freeze-4`, and identical fraction vectors across fuels. `raw_file` requires a value match against a file under `extras/fuelblend_data/data/raw` or `data/raw`.
- [ ] Write `AUDIT.md`: counts per class per family, and the list of v2/v3 deliverables that consume non-`raw_file` rows (grep the scripts that read `inventory_v3.csv`/`compositions_v3.csv`).
- [ ] Run the script and the tests (they must PASS); commit.

### Task 2: Identity consistency check

**Files:** Create `scripts/v4/check_identities.py`. Output `results/v4/P1/identity_mismatches.csv`. Test `tests/test_v4_identities.py`.

**Produces:** `inchikey_mismatches(df: pd.DataFrame, smiles_col="smiles", key_col="inchikey") -> pd.DataFrame`, which recomputes the key with RDKit and skips lumps and blank SMILES.

- [ ] Write `test_nasa_diesel_identity_mismatch_detected` (≥3 mismatches in `compositions_v3.csv` among the `NASA_LFP*` systems) and `test_canonical_registry_clean` (0 mismatches in `compounds.csv`).
- [ ] Implement, run over every composition table and the `compounds.csv` files; commit.

### Task 3: Quarantine and retractions

**Files:** Create `data/quarantine/README.md`, `scripts/v4/quarantine.py`. Modify `docs/claims_register.md`, `data/MISSING.md`. Move tests to `tests/quarantine/` (add it to `norecursedirs` in `pytest.ini`).

- [ ] The script copies every ledger row with a class other than `raw_file`/`cited_typed`, plus every identity-mismatch row, into `data/quarantine/{inventory,compositions}_quarantined.csv` with a `reason`. It writes `data/real_fuels/v4/{inventory_v4.csv,compositions_v4.csv}` holding only clean rows, plus a `provenance` column. Fuels left with no clean composition are dropped from v4 and listed in `MISSING.md`.
- [ ] `cited_typed` rows (values typed from a named paper but not in any local file) stay in v4 with `verified=False`. They are excluded from headline metrics until a raw file confirms them.
- [ ] Add claims rows. Status RETRACTED for: G2 octane, G6 DVPE, G8 diesel cetane/cold flow, the U6 designs and the U2 "30/30 robustness" check. Status REVISED-PENDING for: G10 density/viscosity, R3 real-fuel flash (subset), R4 surrogates (only if a surrogate's evaluation used quarantined data), and the `uncertainty_v3.yaml` margins. Each row cites the ledger and the code lines.
- [ ] Move v3 tests that pin quarantined numbers (inspect `test_v3_phase1_expansion.py`, `test_v3_phase2_repairs.py`, `test_v3_phase3_designs.py`, `test_v3_r4_surrogates.py`, the relevant cases in `test_v3_phase0_regression.py`) into `tests/quarantine/`, each with a module docstring naming the claim ID. Keep the CRC-only cases active.
- [ ] Run `pytest tests -q`: all active tests pass. Commit.

### Task 4: Pure-component constant audit (F1)

**Files:** Create `scripts/v4/audit_pure_constants.py`. Output `results/v4/P1/pure_constant_diffs.csv`. Test `tests/test_v4_pure_constants.py`.

**Produces:** `reference_value(smiles: str, prop: str) -> tuple[float, str] | None`, using the `chemicals`/`thermo` libraries and `compounds.csv`. Properties: Tm, Tb, ρ(15 °C), MW, flash, ΔHfus, Psat(37.8 °C); RON/MON/DCN come from the Co-Optima and NREL cetane files on disk.

- [ ] Compare `data/processed/v3/purchasable_compounds.csv` and the hard-coded table in `src/fuelblend/design/inverse_design_v2.py`. Flag |Δ| > 2 °C (temperatures), > 1 % (density), > 3 ON/DCN, > 5 % (Psat).
- [ ] Test `test_known_reference_values`: n-dodecane Tm within 0.5 °C of −9.6 °C; n-octane RON flagged if its listed value differs from the Co-Optima/NREL value by more than 3 ON (currently listed −10).
- [ ] Write `data/processed/v4/pure_props_v4.csv`: one value per compound and property, with `source`. This is what Phase 4 consumes. Commit.

## Phase 2 — Rebuild the real-fuel corpus from traceable data

### Task 5: Gasoline blends with measured octane (FACE set I/II + Co-Optima)

**Files:** Create `src/fuelblend/data/loaders_octane.py`, `scripts/v4/build_gasoline_octane.py`. Output `data/real_fuels/v4/gasoline_octane_blends.csv`. Test `tests/test_v4_octane_loaders.py` + fixtures copied from the real headers.

**Produces:** `load_face_sets(paths) -> pd.DataFrame` and `load_cooptima_bron(path) -> pd.DataFrame`, both returning `blend_id, base_fuel, components (JSON [{name, inchikey|LUMP:<base>, frac, basis}]), ron, mon, provenance`.

- [ ] Inspect the real files first. The FACE CSV header (`Name,n-N,C_1,n-P,C_6,…,RON,MON`) alternates component-name and fraction columns. Document the decoded meaning in the loader docstring.
- [ ] A FACE base gasoline is a pseudo-component `LUMP:FACE_A` etc. Its RON/MON is used only as a model input, never as a target, for blends containing it.
- [ ] Write `test_face_prf10_parses`: the row `PRF+10` contains n-heptane, iso-octane and ethanol with fractions summing to 1 ± 0.01, RON 55.6, MON 53.
- [ ] Write `test_bron_units_and_basis` (Review Focus 5): blend fractions are volume-basis, RON is in 0–130, and the row count matches the file.
- [ ] Implement, run, report counts in `results/v4/P2/gasoline_corpus.md`; commit.

### Task 6: Jet fuels with cited class-level compositions

**Files:** Create `scripts/v4/build_jet_corpus.py`. Output `data/real_fuels/v4/jet_class_compositions.csv`. Test `tests/test_v4_jet_corpus.py`.

**Consumes:** the `fuel_hindcast.py` bulk route (class × carbon-number pseudo-components).

- [ ] Map `raw/gc_composition/fuel_gc_compositions.csv` (POSF-10325/13511/13516, JetA_API14avg, JP-5/JP-8 Potter-Simmons, UDRI HTA) to class × carbon-number pseudo-components. Where the carbon-number split isn't given, use D86 from a raw file for that same fuel only. If there is none, mark it `assumed_cn_split` and exclude it from headline metrics.
- [ ] Take property targets only from raw files (`crc647_fuel_properties.csv`, `av2718_fuel_properties.csv`, `rp1_properties.csv`, the AJFTD processed tables). For POSF fuels, if no local file holds the property, record NaN and add a `MISSING.md` row naming the document to fetch (Task 8).
- [ ] Test `test_posf_10325_class_sums`: class fractions sum to 1 ± 0.02, and no two POSF fuels have identical fraction vectors.
- [ ] Commit.

### Task 7: Diesel corpus decision and splits v4

**Files:** Create `scripts/v4/build_splits_v4.py`. Output `splits/v4/real_fuels_v4_splits.json`, `results/v4/P2/coverage_matrix.csv`. Test `tests/test_v4_splits.py`.

- [ ] Diesel: include only fuels with raw-file composition (the Kirgina/NASA data, if a raw file exists in `extras/.../data/raw/pour_point` or `wax`; check first). If fewer than 5 qualify, record in `docs/scope_decisions.md` (SD-02) that the diesel real-fuel evaluation waits on Task 8 documents.
- [ ] Splits: 5-fold GroupKFold by `fuel_id`, plus a locked test set (20 %, seed 20261030) that no tuning touches. Test `test_no_fuel_in_two_folds` and `test_locked_set_disjoint`.
- [ ] Commit.

### Task 8: Optional manual-download checklist (user-assisted)

**Files:** Create `docs/v4_manual_downloads.md`. Ingesters for these are written only once a file appears in `data/raw/<name>/`.

- [ ] List each document with: what it contains, which task it unblocks, and the target folder. Candidates:
  - AFRL-RQ-WP-TR-2020-0017 (Edwards, *Jet Fuel Properties*). DTIC blocks scripts, so a browser download is needed.
  - CRC AV-07-02 / FACE diesel characterization.
  - FACE gasoline DHA supplementary data (CRC AVFL-24 / Sarathy et al.).
  - Purdue GCxGC jet-fuel studies (Vozka, Kilaz et al.) with measured properties.
  - Liaw et al. hydrocarbon mixture flash-point papers (ternary/quaternary data).
  - Any newer ThermoML archive release than 2020-09-30 on data.nist.gov. Scripts can check this.
- [ ] For each, state the fallback if the user cannot get it (task runs on current data; claim scope narrowed).
- [ ] Commit.

## Phase 3 — Re-evaluate on clean data

### Task 9: Flash point (activity coefficients + real-fuel re-hindcast)

**Files:** Create `scripts/v4/eval_flash_v4.py`, `src/fuelblend/baselines/flash_unifac.py`. Output `results/v4/F/flash_lab_loso.csv`, `results/v4/F/flash_real_fuels.csv`. Test `tests/test_v4_flash.py`.

**Produces:** `flash_point_unifac(smiles: list[str], x: np.ndarray, P_kPa=101.325) -> float` (°C). It uses the existing anchored Psat + LFL and `thermo` UNIFAC γ, with the same brentq solver as `baselines/flash.py`.

- [ ] Test `test_unifac_reduces_to_ideal_for_alkane_pair` (|Δ| < 0.5 °C for n-decane/n-dodecane) and `test_unifac_lowers_flash_for_aromatic_alkane` (γ > 1 → flash ≤ ideal).
- [ ] Evaluate ideal vs UNIFAC on the 126-mixture LOSO set (31 systems). Report MAE, RMSE, the over-prediction > 5 °C rate, and the worst over-prediction.
- [ ] Re-run the real-fuel hindcast only on v4 fuels with `raw_file` flash. Split by D56/D93 method. Report per-family N honestly.
- [ ] Commit.

### Task 10: Conformal margins computed from residuals (replaces hand-typed margins)

**Files:** Create `src/fuelblend/uq/conformal.py`, `scripts/v4/build_uncertainty_v4.py`. Output `configs/uncertainty_v4.yaml` (generated; header says "do not edit"). Test `tests/test_v4_conformal.py`.

**Produces:**
- `one_sided_margin(residuals: np.ndarray, groups: np.ndarray, alpha: float = 0.10, direction: str = "over") -> float`. It takes the per-group max residual and returns the finite-sample conformal quantile ⌈(n+1)(1−α)⌉/n over n groups.
- `two_sided_margin(...)` for density and viscosity (log space).

Steps:
- [ ] Tests:
  - `test_margin_matches_hand_example`: 9 groups with max residuals 1..9, α = 0.1 → 9.0.
  - `test_margin_infinite_when_too_few_groups`: n < 10 at α = 0.1 → `inf` + warning (Review Focus 3).
  - `test_margin_uses_groups_not_rows`: duplicating a group's rows doesn't change the result.
- [ ] Build the yaml from the out-of-fold residual CSVs of Tasks 9, 11, 12, 13, 14 and 15. Each entry carries `{value, n_groups, residual_file, model}`. Density uses the real 153→v4 residuals (this resolves the 12.94 vs 2.50 inconsistency).
- [ ] Commit.

### Task 11: Viscosity — nested CV and physics-structured residuals

**Files:** Create `scripts/v4/build_viscosity_leaderboard_v4.py`, `src/fuelblend/mixture/structured_residual.py`. Output `results/v4/V/viscosity_leaderboard_v4.csv`. Test `tests/test_v4_viscosity.py`.

**Produces:** `StructuredResidual(features: list[str], kind: "ridge"|"gp").fit(X, T_K, y_resid)`. The residual model is `Δlnν = a0(z) + a1(z)·1000/T`, where `z` holds the fuel descriptors (class fractions, mean/std carbon number, aromatics). Its T-dependence is linear in 1000/T by construction, so it cannot blow up at −40 °C.

- [ ] Tests:
  - `test_structured_residual_is_linear_in_invT`: the second difference over 1000/T is ≈ 0.
  - `test_nested_cv_selects_inside_fold`: the selected hyperparameters vary by outer fold and are never fit on the held-out fuel.
- [ ] Leaderboard rows on the 57 CRC fuels, LOO by fuel:
  - raw Refutas
  - recalibrated Refutas
  - Walther
  - bilinear pairwise (seeded; rank and λ chosen by inner CV)
  - structured ridge
  - structured GP (report the GP's 90 % interval coverage)
  - a lab-trained structured residual applied zero-shot (to test the cryogenic collapse fix)
- [ ] Include cluster-bootstrap CIs (existing `cluster_bootstrap_ci`) and a paired bootstrap of each model vs recalibrated Refutas.
- [ ] Commit.

### Task 12: Density, heat of combustion, freezing — re-run on v4 corpus

**Files:** Create `scripts/v4/eval_bulk_props_v4.py`. Output `results/v4/B/{density,hoc,freezing}_v4.csv`. Test `tests/test_v4_bulk_props.py`.

- [ ] Re-run the ideal-volume density, group-additivity NHOC, ideal SLE and Coutinho (`src/fuelblend/thermo/coutinho_wax.py`, the library version only) on v4 fuels with `raw_file` targets. Use the same metrics as v3 so the results can be compared side by side.
- [ ] Test `test_v4_inputs_are_clean`: every evaluated fuel has `provenance` other than quarantined.
- [ ] Commit.

## Phase 4 — New properties

### Task 13: Freezing/WAT from ThermoML solid–liquid equilibrium

**Files:** Modify `extras/fuelblend_data/thermoml_extract.py`, adding an `--sle` mode that leaves the default whitelist untouched; copy the result to `data/processed/v4/thermoml_sle_hc.csv`. Create `scripts/v4/eval_sle.py`, output `results/v4/SLE/sle_eval.csv`. Test `tests/test_v4_sle.py` with a real ThermoML SLE record copied to `tests/data/` (find one in the tgz with a hydrocarbon binary, e.g. n-alkane + aromatic).

- [ ] Extract liquidus temperatures (phase-equilibrium T with phases Liquid + Crystal), with liquid composition, for all-hydrocarbon systems. Keep the existing phase rules (AGENTS.md §9 item 1).
- [ ] Test `test_sle_record_parsed`: the fixture gives the expected T and x. Test `test_default_mode_unchanged`: the existing `test_thermoml.py` still passes.
- [ ] Evaluate ideal Schroeder–van Laar vs Coutinho on the SLE binaries/ternaries (GroupKFold by `system_id`). Use the result to set the freezing margin (Task 10) for simple blends.
- [ ] Commit.

### Task 14: Smoke point / sooting (D1655 constraint)

**Files:** Create `src/fuelblend/baselines/sooting.py`, `scripts/v4/eval_smoke_point.py`. Output `results/v4/SP/smoke_point_eval.csv`. Test `tests/test_v4_sooting.py`.

**Produces:**
- `ysi_mix(smiles, x_mole) -> float`: linear-by-mole blending on the unified YSI scale from `raw/sooting/ysi_compounds.csv`.
- `smoke_point_from_ysi(ysi, mw) -> float` (mm), fit on `tsi_table.csv` pure compounds that have both values.
- `smoke_point_from_bulk(aromatics_vol, diaromatics_vol) -> float`, fit on `smoke_point_ajftd.csv`.

Steps:
- [ ] Tests:
  - `test_ysi_unified_scale_anchor`: benzene ≈ 100.3 and n-heptane ≈ 36.0 on the stored scale (Review Focus 5).
  - `test_ysi_mix_linear`: a 50/50 mix gives the mean.
  - `test_smoke_point_censoring_respected`: rows with `censored=True` (> 50 mm) are not used as exact targets.
- [ ] Evaluate:
  - (a) the bulk model with GroupKFold by `fuel_id` on AJFTD;
  - (b) the YSI route on any fuel with both speciation and a measured smoke point (`sooting_real_fuels.csv`).
- [ ] Add a smoke-point constraint to `configs/specs_astm.yaml` (D1655: ≥ 25 mm, or ≥ 18 mm with naphthalenes ≤ 3.0 vol %) and to `specs_bis.yaml`. Verify both limits against the standard text before adding them; leave a note if they can't be confirmed.
- [ ] Commit.

### Task 15: Octane refit on real blends

**Files:** Modify `src/fuelblend/models/octane_ghosh_lumps.py` (add a `fit(blends)` path that accepts base-fuel lumps). Create `scripts/v4/eval_octane_v4.py`. Output `results/v4/O/octane_eval_v4.csv`. Test `tests/test_v4_octane.py`.

- [ ] Models: linear-by-volume, linear-by-mole, Ghosh lumps, and a GBM/GP on (blendstock bRON, volume fraction, base-fuel RON/MON/sensitivity). CV is GroupKFold by `base_fuel` (and by blendstock for a leave-component-out view).
- [ ] Test `test_prf_linear_by_volume_exact`: an n-heptane/iso-octane PRF gives RON equal to vol % iso-octane within 1 ON.
- [ ] Re-derive the DVPE evaluation **only** if a raw file with measured DVPE exists (check `gasoline_fuel_properties.csv`, Co-Optima). Otherwise record SD-03 (DVPE evaluation pending data).
- [ ] Commit.

## Phase 5 — Learned models

### Task 16: Mixture encoder retrain (laptop overnight / GPU job pack)

**Files:** Create `src/fuelblend/mixture/set_encoder_v4.py`, `scripts/v4/train_set_encoder.py` (CLI: `--device`, `--seeds`, `--epochs`, `--out`), `jobs/v4_gpu/README.md` + `run_all.sh`. Output `results/v4/NN/*.csv`, `runs/v4/*.pt` (gitignored). Test `tests/test_v4_set_encoder.py`.

**Produces:** `SetEncoder(component_featurizer, pool="attention"|"sum", heads: dict[str, Head])`. Each head outputs a residual on top of its physics baseline. The temperature enters only through the structured form from Task 11 (`a0 + a1·1000/T`).

- [ ] Tests:
  - `test_permutation_invariance`: shuffling components leaves the output unchanged.
  - `test_zero_fraction_component_no_effect`.
  - `test_residual_zero_init_returns_physics`.
  - a CPU smoke run (`--epochs 2 --seeds 1`) finishes in < 2 min.
- [ ] Training data: lab mixtures (`model_table.csv`, GroupKFold by `system_id`) + v4 real fuels (GroupKFold by `fuel_id`), using a masked multi-task loss. Deep ensemble of 5 seeds, then split-conformal calibration with Task 10 functions.
- [ ] Report against Task 11's best model. The encoder is promoted to `model_decision` only if it wins under the outer CV, with paired-bootstrap CI excluding 0, and its 90 % coverage is ≥ 85 % on real fuels. Otherwise document that it lost.
- [ ] Commit.

## Phase 6 — Inverse design rebuilt on the library

### Task 17: Design engine using library physics and generated margins

**Files:** Create `src/fuelblend/design/engine_v4.py`, `scripts/v4/run_designs_v4.py`, `scripts/v4/verify_designs_v4.py`. Output `results/v4/D/designs_v4.csv`, `results/v4/D/robustness_v4.csv`. Test `tests/test_v4_design.py`.

**Consumes:** `pure_props_v4.csv` (Task 4), the flash functions (Task 9), the viscosity models (Task 11), the Coutinho library function (Task 12/13), `ysi_mix` and smoke point (Task 14), the octane model (Task 15), `uncertainty_v4.yaml` (Task 10), and `configs/specs_*.yaml` (the only spec source; this fixes F5).

**Produces:** `design(spec_name: str, pool: list[str], models: dict, n_starts: int = 30, seed: int = 42) -> DesignResult` with fields `weights, predictions, margins, passed, extrapolating, binding_constraints`.

- [ ] Tests:
  - `test_design_flags_extrapolation` (Review Focus 4): a pool forcing > 60 wt % of a class absent from validation data gives `extrapolating=True` and `passed=False`.
  - `test_specs_read_from_yaml_only`: no numeric spec literals in `engine_v4.py` (AST check).
  - `test_volatility_minimum_enforced`: a gasoline below the yaml DVPE minimum fails.
  - `test_robustness_runs_each_model`: the robustness table calls each model function (mock counts) instead of adding offsets.
- [ ] The verifier recomputes Psat, flash and SLE with the `thermo`/`chemicals` libraries and shares no property code with the engine. Report disagreements > margin.
- [ ] Re-run jet, diesel and gasoline (ASTM + BIS). Designs that no longer pass are reported as such.
- [ ] Commit.

### Task 18: Prediction API cleanup

**Files:** Modify `src/fuelblend/api/predict_blend.py`. Test `tests/test_predict_blend.py` (extend).

- [ ] Remove the lazy retrain at load. A missing bundle raises `FileNotFoundError` with the command to build it. Route each property to the `docs/model_decision.md` v4 choice, and return `(value, lo, hi, model_name, in_domain)` using `uncertainty_v4.yaml`.
- [ ] Test `test_missing_bundle_raises_not_trains` and `test_interval_from_v4_yaml`. Commit.

## Phase 7 — Reporting and reproducibility

### Task 19: v4 orchestrator, report, handoff

**Files:** Create `scripts/v4/reproduce_all.py` (mirrors v3: steps in order, pytest, deliverable row-count verification, `SHA256SUMS_v4.txt`, `--verify-only`). Create `docs/FULL_REPORT_v4.md`. Modify `docs/model_decision.md` (v4 section appended), `docs/PROJECT_HANDOFF.md`, `AGENTS.md` §12.

- [ ] The report states, per property: data N and provenance mix, best model under outer CV with CI, margin with `n_groups`, what changed vs v3, and what is still unverified. Every number is transcribed from a `results/v4` CSV and cited.
- [ ] Run `python scripts/v4/reproduce_all.py` end to end twice; all deliverables must be byte-identical across the runs (determinism check).
- [ ] Final `pytest tests -q` all pass; commit; do not push (user decides).

---

## Self-review notes

- Coverage of user requests:
  - script improvements: Tasks 0, 4, 11, 17, 18
  - missed issues: Tasks 1–3
  - upgrades: Tasks 9, 11, 16
  - new data: Tasks 5, 6, 8, 13
  - existing data for new properties (user picked smoke point, SLE freezing, bRON octane): Tasks 13–15
  - retraining: Tasks 11, 15, 16
- Order matters: Phase 1 must finish before any Phase 3+ number is reported. Tasks 5, 6, 13 and 14 are independent of each other and can run in parallel after Task 3.
