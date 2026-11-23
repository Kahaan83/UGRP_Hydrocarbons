# v8: ML-filled component constants and design certification Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:executing-plans (native, in one session). No subagents. Steps use checkbox (`- [ ]`) syntax for tracking.

## User decisions (recorded 2026-11-23, before any v8 result)

* G1 tolerances: flash 2 C, freezing 2 C, density 2 kg/m3, viscosity 0.02 ln.
* `splits/locked.json`: restore only this file to its committed version in Task 0 (BASE files untouched); locked readings (G4) run.
* Git: commit after each task; push to `origin feat/external-datasets-ingest` once, after Task 7.
* Execution: native, no subagents.

## Unattended execution rules

The plan runs start to finish without asking the user. All decisions are above.
* A failed gate is a reading, not a stop: write it in that task's `READING.md`, commit, and continue. Later tasks use only what passed (a constant kind that fails G1 is not used in Tasks 3 to 5; a property that fails G2 stays `uncertified`).
* If the Task 3 first pass shows freezing `no_skill_shown`, continue anyway; G5 is then expected to move only among P0 to P0y-b, and the reading says so.
* If P4 does not form (Task 5 Step 1), skip to Task 6 with 16 rows.
* Stop only on a broken invariant: `reproduce_base.py` failing, a test that cannot be made to pass without changing a gate, or any write into `results/v5/BASE/`. Then commit what is done, write `results/v8/STOPPED.md` with the cause, and do not push.
* Long runs (D-MPNN refits, design runs) go in the background with resumable outputs; logs go to `scratch/`.

**Goal:** Use the v7 pure-compound models to supply missing component constants to the v5 physics blend routes, measure exactly how much accuracy that costs, earn real-fuel margins for GCxGC-represented blendstocks, and re-run inverse design under unchanged gates.

**Architecture:** Physics stays the blending backbone (v5 routes). ML enters at the component level: a leak-free table of pure-compound predictions fills constants a component lacks. Every ML-filled or route constant keeps a non-verified status, and a constraint that uses one is certifiable only with a margin derived on the same input pipeline (new file `configs/uncertainty_v8.yaml`). `configs/uncertainty_v5.yaml` and `results/v5/BASE/` are never changed; with the v8 margin set off the engine reproduces v5 exactly.

**Tech Stack:** Python 3.12 (`.venv`), pandas, numpy, scikit-learn, RDKit, existing `fuelblend` packages; the GPU venv only for D-MPNN refits in Task 1.

**Spec:** this plan and the decisions recorded in Task 0 (`docs/v8_predeclaration.md`). Evidence it argues from: `docs/final_numbers_v5.md` (sections C, E, F), `docs/final_numbers_v7.md`, `docs/v7_errors_by_property.md` (section E), `results/v5/BASE/delta_adopt_designs/log_P1_astm_jet_a.txt`, `results/v5/BASE/delta_ws6_gcxgc/proposed_margins.yaml`.

## Why designs are not certified today (from committed files)

* P0 to P0y-b (known compounds): fail on slack (flash, freezing, D86), not on missing margins. 3 of 16 certified.
* P1 to P3 (real blendstocks): up to 11 constraints `uncertified` for two reasons in the log: (a) stand-in compounds use non-verified constants (`rho15=route_grid`, `flash=route_anchored`, `Tm/dHfus=coutinho_table`); (b) `represented_blendstock_margin_not_validated` (margins exist only for known compositions; the WS6 real-fuel margins were proposed, not adopted).
* So the honest lever is (b) plus (a): margins earned on real GCxGC fuels with the same constants pipeline. ML-filled constants add components to that pipeline.

Expectation context (committed numbers, not targets): ML pure errors are large next to blend errors (Andrade ln viscosity ridge 0.284, log10 vapour pressure ridge 0.220, melting point random forest 27.38 K). Task 2 measures whether they survive substitution; a fail is a reportable result.

**Likely ceiling.** A design row certifies only if every constraint does, so freezing decides every P1 to P3 row. Real-fuel UNIFAC SvL is worse than a constant (8.425 vs 6.413 C pooled) with bias +4.3 to +5.9 C on the safe side, so its coverage is uninformative. If freezing does not earn a skill-backed margin in Task 3 (G2), the certified count can only move among P0 to P0y-b and will likely stay near 3; v8 still delivers the substitution result (G1).

**Execution order:** Task 0, then Task 3 first pass (existing stand-ins only, no ML constants; it needs nothing from Tasks 1 and 2) as the go/no-go, then Tasks 1, 2, a Task 3 second pass with ML constants for kinds that passed G1, then Tasks 4 to 8. The Task 3 first-pass reading is committed and execution continues (see Unattended execution rules).

**P4 is conditional.** None of the 7 candidates in `results/v5/BASE/delta_representatives/candidate_gaps.csv` has measured Antoine constants (`has_antoine` False for all), and trap 3 forbids ML or un-refit Antoine. P4 exists only for candidates whose Antoine can be refit from measured vapour pressure already on disk (source `refit_measured`).

## Global Constraints

* Never invent or alter numbers; every figure traces to a committed file. Gate tolerances are fixed in Task 0 before any v8 result exists.
* Do not change `data/unified/`, `scripts/v5/`, `configs/uncertainty_v5.yaml`, `results/v5/BASE/**` (manifest hashes exact bytes; never `git checkout` or normalise them). New code goes in `src/fuelblend/v8/`, `scripts/v8/`, results in `results/v8/`.
* `GroupKFold` on `system_id`; train and evaluate from `measurements_dedup.csv`; hydrocarbon-only (`is_hc`) default; `data/quarantine/` and `DatasetsFuelProperties.xlsx` never used.
* Margins: one-sided 90 % conformal, `k=ceil((n+1)(1-alpha))` of per-group max unsafe-side residuals, `inf` if `k>n` (same rule as `configs/uncertainty_v5.yaml`). Unsafe side: flash over-prediction, freezing under-prediction, viscosity as v5 (two-sided ln), density as v5.
* Trap 3: no un-refit Antoine for C10+ alkanes in the 40 to 100 C window. Antoine constants predicted by parameter heads are not used (v7: +16.4 C over-prediction).
* Locked sets are scored once, after dev readings are committed. Dev-only gains cannot reopen a v5 or v6 verdict.
* Real-fuel freezing and cold viscosity keep "status: no skill shown" unless a Task 3 margin passes its transfer reading.
* Design protocol unchanged: `n_starts=10`, `seed=42`, `strict_constraints=True`, `synthetic_cap` 0.50 for P2/P3, flash overlay on.
* No AI attribution in commits. No em dashes in manuscript text. Tests: `.venv\Scripts\python.exe -m pytest -q <file>`; full suite only at the end (baseline 741 passed, 11 skipped).

## Review Focus

1. A compound that has a measured constant and an ML prediction: the measured value must always win; ML fills gaps only.
2. Leakage: a compound inside the v7 pure training pool must get its out-of-fold prediction in Task 2, never the deployed model's in-sample value.
3. Out-of-domain components (oxygenates, heavier than the pure pool's carbon range, unparseable SMILES): status `ml_out_of_domain`, value not used, constraint stays `uncertified`.
4. v8 margin set off (the default): every v5 design and `reproduce_base.py` output is byte-identical.
5. A margin that comes out `inf` (`k>n`) or fails its transfer reading must leave the constraint `uncertified`, never certified with a finite fallback.

---

### Task 0: Predeclaration and user decisions

**Files:**
- Create: `docs/v8_predeclaration.md`

- [ ] **Step 1: Write `docs/v8_predeclaration.md`** with: the five gates below, the sets each uses, the reading text for pass and fail, and the date. Commit before any v8 script exists.

| gate | question | pass rule (fixed; values from the User decisions section) |
|---|---|---|
| G1 substitution | Does replacing a measured constant by its out-of-fold ML value keep blend accuracy? Per property: density (`rho15`), viscosity (Andrade from ML log viscosity at 293 K plus the v7 descriptor-ridge slope), flash (LFL from ML pure flash point via `Psat(Tf)=LFL*P`, only where measured Antoine exists), freezing (`Tm`, `dHfus`) | upper end of the paired group-bootstrap 90 % CI of dMAE (ML minus measured) at or below: density 2 kg/m3, viscosity 0.02 ln (both traced to the `delta_representatives` reading); flash 2 C, freezing 2 C (user decision 2026-11-23) |
| G2 real-fuel margin | Can a margin be earned for represented blendstocks on real GCxGC fuels (AJFTD + AV-33-22)? | all three: (a) the route beats the leave-one-group-out constant on the same fuels (paired group bootstrap, CI of the gain excludes 0); (b) margin finite; (c) the WS6 R-WS6-3 leave-source-out coverage reading holds in both directions. Coverage alone never passes (freezing coverage is uninformative under safe-side bias) |
| G3 no regression | With the v8 margin set off, are v5 outputs unchanged? | `reproduce_base.py` passes and re-evaluating the committed v5 design weights with `evaluate_recipe` gives identical constraint statuses and values |
| G4 locked | Do G1 passes hold on locked systems? | same rule as G1, through `src/fuelblend/v7/guard.py`. Viscosity and density lab locked systems were already read by v5 R-MLV-2 / R-MLD-2, so their G4 rows are labelled `already_read`; only unspent locks give a fresh reading |
| G5 design | How many of the 16 v5 rows (plus P4 x 2 specs if P4 forms) certify under v8? | descriptive count; no target. Each certified row lists which v8 margins it used |

- [ ] **Step 2: Copy the "User decisions" section of this plan into the file verbatim.**
- [ ] **Step 2b: Restore `splits/locked.json` only:** `git restore splits/locked.json`, then `git status --short splits/ results/v5/BASE/` to confirm the BASE files are still listed as modified (untouched). Run `.venv\Scripts\python.exe -m pytest -q tests -k s4_splits`; expect PASS.
- [ ] **Step 3: Commit** `git add docs/v8_predeclaration.md; git commit -m "docs(v8): predeclare ML-constant substitution and real-fuel margin gates"`

### Task 1: Leak-free ML constants table

**Files:**
- Create: `src/fuelblend/v8/__init__.py`, `src/fuelblend/v8/ml_constants.py`, `scripts/v8/build_oof_constants.py`
- Test: `tests/v8/test_ml_constants.py`

**Interfaces:**
- Consumes: `fuelblend.api.predict_pure.predict_pure(smiles, props)` (deployed values); v7 fold assignments and OOF predictions in `results/v7/track_b/descriptor_gbm/` and `results/v7/track_b/curves_v1_long/`.
- Produces:
  - `ml_constant(smiles: str, prop: str) -> tuple[float, str]` with prop in `{"rho15", "Tm", "dHfus", "logvisc293", "flash_pure"}`; status `ml_rf_rdkit2d`, `ml_dmpnn_long`, `ml_out_of_domain`, or `invalid_smiles`; Tm returned in C, dHfus in kJ/mol, rho15 in kg/m3.
  - `oof_constants(prop: str) -> pd.DataFrame` columns `inchikey, smiles, value, fold, method` (out-of-fold for every compound in the pure pool).
  - `results/v8/oof_constants.csv` (all props, long format).

- [ ] **Step 1: Confirm the temperature basis of the Track B density label** from the pool loader used by `scripts/v7/track_b/` and write it in the module docstring. If it is not 288.15 K, `ml_constant(..., "rho15")` converts with the pure-property grid slope where the compound has one, else returns `ml_out_of_domain`.
- [ ] **Step 2: Write failing tests**

```python
def test_measured_value_not_returned_by_ml_layer():
    v, st = ml_constant("CCCCCCCCCCCC", "Tm")      # dodecane
    assert st.startswith("ml_") and np.isfinite(v)

def test_out_of_domain_oxygenate():
    assert ml_constant("CCO", "rho15")[1] == "ml_out_of_domain"

def test_invalid_smiles():
    assert ml_constant("not_a_smiles", "Tm")[1] == "invalid_smiles"

def test_oof_has_no_in_sample_rows():
    df = oof_constants("Tm")
    assert df["fold"].notna().all() and df["inchikey"].is_unique
```

- [ ] **Step 3: Run** `.venv\Scripts\python.exe -m pytest -q tests/v8/test_ml_constants.py`; expect FAIL (module missing).
- [ ] **Step 4: Implement.** Domain = hydrocarbon (C and H only) and carbon count within the min..max of the property's training pool. `oof_constants` reads committed OOF predictions where they exist; otherwise `build_oof_constants.py` refits the deployed method per v7 scaffold fold (D-MPNN refits in the GPU venv) and writes `results/v8/oof_constants.csv`.
- [ ] **Step 5: Run tests; expect PASS.**
- [ ] **Step 6: Commit** `feat(v8): leak-free ML pure-constant table`

### Task 2: Substitution test on lab mixtures (gate G1)

**Files:**
- Create: `src/fuelblend/v8/substitution.py`, `scripts/v8/substitution_test.py`
- Test: `tests/v8/test_substitution.py`
- Output: `results/v8/substitution/{rows.csv, summary.csv, READING.md}`

**Interfaces:**
- Consumes: `oof_constants(prop)` (Task 1); the v5 designated routes in `fuelblend.api.predict_route` (ideal volume, Refutas, Le Chatelier with LFL overlay, UNIFAC SvL).
- Produces: `score_substitution(prop: str, set_name: str = "dev") -> pd.DataFrame` with one row per system: `system_id, status, n_rows, mae_measured, mae_ml, bias_ml, n_substituted` (`status` in `ok`, `skipped_no_measured_antoine`, `skipped_no_ml_value`); `summary.csv` with dMAE and paired group-bootstrap 90 % CI per property.

- [ ] **Step 1: Write failing tests**

```python
def test_identity_when_ml_equals_measured(monkeypatch):
    # if the ML table returns the measured constants, dMAE is exactly 0
    df = score_substitution("density", set_name="dev", _ml_override="measured")
    assert np.allclose(df["mae_ml"], df["mae_measured"])

def test_only_systems_with_measured_truth_scored():
    df = score_substitution("freezing", set_name="dev")
    assert (df["n_substituted"] >= 1).all()

def test_flash_skips_components_without_measured_antoine():
    df = score_substitution("flash", set_name="dev")
    skipped = df[df["status"] == "skipped_no_measured_antoine"]
    assert skipped["mae_ml"].isna().all()                          # never scored with a substitute Antoine
    assert (df.loc[df["status"] == "ok", "n_substituted"] >= 1).all()
```

- [ ] **Step 2: Run; expect FAIL.**
- [ ] **Step 3: Implement.** One constant kind substituted at a time, all components of the system substituted together. Flash LFL from ML flash point: `LFL_i = 100 * Psat_i(Tf_i) / P_atm` with the component's measured Antoine. Viscosity: `ln nu(T) = a + b/T`, b from the v7 descriptor ridge (`src/fuelblend/ml/param_heads.py`), a anchored to the ML log viscosity at 293.15 K.
- [ ] **Step 4: Run tests; expect PASS.**
- [ ] **Step 5: Run** `.venv\Scripts\python.exe scripts/v8/substitution_test.py --set dev > scratch/v8_sub.log`; write `READING.md` with the G1 verdict per property using the Task 0 rule verbatim.
- [ ] **Step 6: Commit** `exp(v8): ML-constant substitution test on dev mixtures`
- [ ] **Step 7: Only for properties that passed on dev:** run `--set locked` once through `src/fuelblend/v7/guard.py` (G4); label viscosity and density rows `already_read`; append to `READING.md`, commit `exp(v8): locked substitution reading`.

### Task 3: Real-fuel margins for represented blendstocks (gate G2)

**Files:**
- Create: `src/fuelblend/v8/represented_margins.py`, `scripts/v8/build_represented_margins.py`
- Test: `tests/v8/test_represented_margins.py`
- Output: `results/v8/represented_margins/{oof.csv, margins.csv, READING.md}`, `configs/uncertainty_v8.yaml`

**Interfaces:**
- Consumes: the real GCxGC fuel set and representation used in `results/v5/BASE/delta_ws6_gcxgc/` (AJFTD + AV-33-22) and `fuelblend.design.blendstock_v5` stand-in construction; `ml_constant` for stand-ins lacking constants (only property kinds that passed G1).
- Produces: `represented_margin(prop: str) -> dict` with keys `value, n_groups, direction, coverage_a_by_b, coverage_b_by_a, reading`; yaml entries keyed `<property>.real_fuels_represented.<route>` in the v5 yaml schema (`value`, `applied_margin`, `unit`, `direction`, `alpha`, `n_groups`, `group_unit`, `conformal_rank_k`).

- [ ] **Step 1: Write failing tests**

```python
def test_inf_when_too_few_groups():
    m = represented_margin_from_residuals(np.array([1.0, 2.0]), alpha=0.1)
    assert m == float("inf")

def test_no_skill_route_gets_no_margin_even_if_covered():
    # safe-side bias: every residual on the safe side gives coverage 1.0, but the route loses to the constant
    r = g2_reading(pred_err=np.full(30, 8.0), const_err=np.full(30, 6.0), coverage=(1.0, 1.0), margin=0.0)
    assert r == "no_skill_shown"

def test_failed_transfer_not_written():
    cfg = yaml.safe_load(open(V8_YAML))
    for k, e in cfg["margins"].items():
        assert e["reading"] == "held" and np.isfinite(e["applied_margin"])
```

- [ ] **Step 2: Run; expect FAIL.**
- [ ] **Step 3: Implement.** Properties: density 15 C, flash, viscosity -20 C, freezing, NHOC, D86 (T10, FBP, T50-T10, T90-T10). Predictions use exactly the design pipeline (same stand-ins, same constant statuses). Residuals per fuel (group = fuel). `g2_reading(pred_err, const_err, coverage, margin) -> str` returns `held`, `no_skill_shown`, `inf_margin` or `transfer_not_held`, checking skill against the leave-one-group-out constant first. Leave-source-out coverage as WS6 R-WS6-3. Only entries whose reading is `held` go into `configs/uncertainty_v8.yaml`; all entries, held or not, go to `margins.csv`.
- [ ] **Step 4: Run tests; expect PASS.**
- [ ] **Step 5: Run the script, write `READING.md` (verdict per property, verbatim rule).**
- [ ] **Step 6: Commit** `exp(v8): real-fuel margins for represented blendstocks`

### Task 4: Engine wiring behind a margin-set switch (gate G3)

**Files:**
- Modify: `src/fuelblend/design/library_route.py` (`build_route_models`, `_flag`), `src/fuelblend/design/blendstock_v5.py` (where `REPRESENTED_REASON` is applied)
- Test: `tests/v8/test_engine_v8_switch.py`

**Interfaces:**
- Consumes: `configs/uncertainty_v8.yaml` (Task 3), `ml_constant` (Task 1).
- Produces: `build_route_models(pool_ids, pure_props=PURE_PROPS_V5, require=..., margin_set: str = "v5") -> dict`; with `margin_set="v8"`, a constraint whose non-verified statuses are all covered by a held v8 entry gets that entry's margin key; ML-filled constants carry status `ml_*` and fill `own_data_gaps` only for kinds that passed G1.

- [ ] **Step 1: Write failing tests**

```python
def test_v5_default_unchanged():
    # re-evaluate the committed v5 weights (no optimiser run) and compare statuses and values
    w = load_committed_weights("P1", "astm:jet_a")
    a = evaluate_recipe("astm:jet_a", build_route_models(pool("P1"), margin_set="v5"), w)
    assert constraint_table(a).equals(load_committed_constraints("P1", "astm:jet_a"))

def test_v8_uses_v8_key_for_represented_density():
    r = design_row("P1", "astm:jet_a", margin_set="v8")
    assert r.constraint("density_15C_kg_m3").margin_key.endswith("real_fuels_represented.ideal_volume")

def test_property_without_held_v8_margin_stays_uncertified():
    r = design_row("P1", "astm:jet_a", margin_set="v8")
    assert r.constraint("freezing_point_C").status == "uncertified"   # unless Task 3 READING says held

def test_measured_constant_wins_over_ml():
    L = RouteLibrary(["n-dodecane"])
    assert L.status("n-dodecane", "Tm") == "verified"
```

`design_row`, `load_committed_weights`, `load_committed_constraints`, `constraint_table` and `pool` are test helpers in the test file; they read `results/v5/BASE/delta_adopt_designs/designs/` and write only under `scratch/`. `design_row` runs `n_starts=2` and is used only for v8 assertions, never compared with a committed 10-start row.

- [ ] **Step 2: Run; expect FAIL.**
- [ ] **Step 3: Implement the `margin_set` argument (default `"v5"`).**
- [ ] **Step 4: Run tests; expect PASS. Run** `.venv\Scripts\python.exe scripts/v5/base/reproduce_base.py`; expect pass (G3).
- [ ] **Step 5: Commit** `feat(v8): margin-set switch for ML-filled and represented constants`

### Task 5: Pool P4 (components admitted through ML constants)

**Files:**
- Create: `scripts/v8/pools.py`
- Test: `tests/v8/test_pools.py`

**Interfaces:**
- Produces: `make_pool(name: str) -> BlendPool | None` for `P4`: P0y plus each `results/v5/BASE/delta_representatives/candidate_gaps.csv` compound that (a) gets a `refit_measured` Antoine from measured vapour pressure already on disk (ThermoML), via the existing refit with clipped starts and per-compound `failed:` status, and (b) has every other missing constant kind filled by a kind that passed G1. Returns `None` if fewer than one candidate qualifies; then G5 covers the 16 v5 rows only.

- [ ] **Step 1: Count measured vapour-pressure rows on disk per candidate** (short Python summary, capped output) and record it in `results/v8/pools/README.md`. If no candidate has enough points for a refit inside its flash window, write that reading, skip Steps 2 to 5, and commit the README.
- [ ] **Step 2: Failing test:** every P4 member has no `absent` status under `margin_set="v8"`; a candidate whose refit status is `failed:` is excluded with its reason logged.
- [ ] **Step 3: Run; expect FAIL. Step 4: implement. Step 5: run; expect PASS.**
- [ ] **Step 6: Commit** `feat(v8): pool P4 from refit Antoine and ML-admitted components`

### Task 6: Design re-run (gate G5)

**Files:**
- Create: `scripts/v8/run_designs.py` (wraps `delta_adopt_designs.run_one` with `margin_set="v8"`, output root `results/v8/designs/`)
- Output: `results/v8/designs/{designs.csv, compare_v5.csv, READING.md}`

- [ ] **Step 1: Write `results/v8/designs/README.md` declaring the protocol (pools P0..P3, plus P4 if it formed, x `astm:jet_a`, `bis:jet_a1`; protocol constants from Global Constraints) and commit it before running.**
- [ ] **Step 2: Smoke run** with `n_starts=2` into `scratch/v8_smoke/` (the runner refuses other values in `results/v8/designs/`).
- [ ] **Step 3: Full run** in the background, resumable per pool; logs to `results/v8/designs/log_*.txt`.
- [ ] **Step 4: Write `READING.md`:** certified count out of 16 (18 if P4 formed), per row the v8 margin keys used, labels (`d86_model_dependent`, `freezing_route_dependent`, `overlay_dependent`, `ml_constant_dependent`), and any component fraction outside AJFTD GCxGC class ranges.
- [ ] **Step 5: Commit** `exp(v8): design re-run under v8 margins`

### Task 7: Numbers, methodology, full suite

**Files:**
- Create: `scripts/v8/build_final_numbers_v8.py`, `docs/final_numbers_v8.md` (generated)
- Modify: `docs/methodology_v7.md` (add section "v8: ML-filled constants"), `docs/GUIDE.md` (status), `docs/v5_negative_results.md` table (failed G1/G2 items)

- [ ] **Step 1: Generator reads only committed v8 CSVs; tables: G1 per property, G2 margins with readings, G5 design table vs v5.**
- [ ] **Step 2: Run it; commit generated doc and docs edits** `docs(v8): final numbers and methodology`.
- [ ] **Step 3: Full suite** `.venv\Scripts\python.exe -m pytest -q > scratch/v8_full.log`; expect 741 + new tests passed, 11 skipped, 0 failed (`splits/locked.json` was restored in Task 0). A failure here is fixed before Task 8; BASE files are never normalised.

### Task 8: Push

- [ ] **Step 1:** `git status --short` and confirm no v8 change is uncommitted and no `ajftd_session.json`, raw AJFTD page or licensed data is staged.
- [ ] **Step 2:** `git push origin feat/external-datasets-ingest` (skipped if `results/v8/STOPPED.md` exists).
- [ ] **Step 3:** Report to the user: G1 to G5 verdicts, certified count vs v5 (3 of 16), and every gate that failed.
