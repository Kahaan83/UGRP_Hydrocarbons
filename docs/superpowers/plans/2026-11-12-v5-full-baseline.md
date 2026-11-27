# v5 Full Baseline Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Produce one scoreboard that says how every property model performs today on the v5 data, with the audit fixes applied, and then measure each further fix or data addition as a change against that frozen scoreboard.

**Architecture:** Three parts. Part A moves the composition-to-prediction path into `src/fuelblend` so every property is scored through one mapping and one statistics module. Part B scores or aggregates every property x evidence-tier cell into `results/v5/BASE/scoreboard.csv`. Part C is a ranked list of improvements, each scored as a delta against the frozen scoreboard. Nothing in `configs/`, `results/v4/` or any existing `results/v5/<task>/` folder is rewritten.

**Tech Stack:** Python 3.12 (`.venv/Scripts/python.exe`), pandas 3, numpy, rdkit, scikit-learn, pytest, PyYAML.

**Spec:** No separate spec file. This plan implements the open items of:
- `C:\Users\kahaa\.claude\plans\i-want-to-move-joyful-platypus.md` (the v5 plan; Phases 0-2 done, Phase 3 held by the user)
- `docs/handoff_data_audit_to_v5_plan.md` and `docs/data_audit_2026-11-11.md` (data audit session)
- `docs/audit_fable_v5_2026-11-11.md` (findings B1-B6, code audit, proposals 1-12)

## Background (state on 2026-11-12, commit `271fdcc`, branch `v5`)

Already done, do not redo:

| Item | Where | Status |
|---|---|---|
| B1 NaN Refutas on zero-mass representative | `results/v5/S3/corrected_2026-11-12_*` (commit `7476658`) | fixed in S3 only; `refutas_nu` itself still unguarded |
| B2 censored freezing counted as measured | same commit; `targets(exclude_censored=True)` | fixed, default is still `False` |
| B3 heat-of-combustion circularity | `results/v5/S4b/` | skill survives on measured rows (n = 172 fuels, 16 groups) |
| Cold viscosity correction | `results/v5/C1/` | K1 per-temperature offset works both directions at -20 and -40 C; not adopted |
| ThermoML 2020 files | `data/real_fuels/v5/thermoml2020_lab.csv`, `results/v5/S1b/` | viscosity held-out test spent (n = 3, not evaluable); 1,440 mixture density rows unused |

Still open and owned by this plan:

1. No single table says how each property is doing. Numbers are spread over 18 result folders with different fuel sets, group rules and censoring rules.
2. Two mapping implementations (`map_fuel` vs S2's `rep_for`, which parses a markdown report); S2's flash-specific "unmapped" rule feeds viscosity and density.
3. Four separate paired-bootstrap implementations (S3, S4, S5, S7).
4. AV-33-22 was scored at 2 of 9 temperatures in S3 (C1 used all 9 for viscosity only); density at 9 temperatures, heat of combustion, smoke point and vapour pressure are unread.
5. Properties with no v5 number: heat of combustion on measured composition, smoke point on measured composition, distillation, pour point, cetane.
6. Pending adoptions (user decisions): LFL overlay, corrected mesitylene Antoine set, K1 cold offset, SLE2 constants tier.
7. B4 (S4 flash correlations unverified), B5 (no v5 claims, no v5 manifest), B6 (declarations not provable from git).

## Global Constraints

- Run everything with `.venv/Scripts/python.exe`; system Python lacks rdkit.
- Result files go only to `results/v5/BASE/` and `results/v5/BASE/delta_<name>/`. Never modify `results/v4/**`, existing files in `configs/**`, `splits/**`, or any file already committed under `results/v5/{S1..S7,S1b,S4b,C1,G1,L1,SLE2,D3,D4,P1}`. Old READMEs get a dated note only if a number in them is superseded. The only permitted writes elsewhere are the ones a task lists under **Files**: new modules under `src/fuelblend/`, new scripts under `scripts/v5/base/`, new tests, the function moves in Tasks 3 and 6 (bodies unchanged, re-exported from the old script), `SHA256SUMS_v5_base.txt`, `data/real_fuels/v5/nistir6659_density.csv` (Task 15), a new `configs/uncertainty_v5.yaml` (Task 10, only after the user adopts), and the docs in Task 17.
- **Designated model per tier** (what "our model" means in each column). T1 and T2: the shipped API's `composition_kind="known"` route in `src/fuelblend/api/predict_blend.py` (raw Refutas, ideal-volume density, ideal Le Chatelier flash, lab-mixture margins from `configs/uncertainty_v4.yaml`). CRC 647 r3: the `composition_kind="reconstructed"` route (recalibrated Refutas inside -40..20 C, CRC margins). T3: no shipped model; `hgb` and `ridge` from S4 are reported against the constant. Freezing, smoke, pour, cetane: no shipped model (API returns uncertified); the physics or regression baseline scored is named in the `model` column.
- **Locked stages are never re-invoked.** Do not run `eval_bulk_descriptor.py locked`, `eval_lab_thermoml2020_holdout_v5.py`, or any script whose README says "scored once". Locked cells are aggregated from the committed CSVs and labelled `already_read`.
- The CRC 647 locked 12 fuels have been read five times. Report them as `fold_b`, not as a locked test.
- Every T2 fuel (AV-33-22, AJFTD GCxGC, FuelLib) has been read. The baseline is labelled `descriptive`; no cell may be called held-out unless Part C Task 15 creates it.
- One fuel per double-count group in every statistic. Groups come from `results/v5/S3/corrected_2026-11-12_groups.csv` (source order AV-33-22, FuelLib, AJFTD).
- Censored values (`qualifier` not null, e.g. "<") are bounds, not measurements: excluded and counted.
- A prediction with status `ok` must be finite. Non-finite predictions are excluded and counted in `n_excluded_nonfinite`.
- Error is `predicted - measured`. Viscosity error is `ln(pred) - ln(meas)`. Flash reports the unsafe side: fraction over-predicted by more than 5 C. Viscosity at T < 283 K reports the unsafe side: fraction under-predicted by more than the S1 margin (0.4154 below 263 K, 0.6084 at 263-283 K).
- Pending adoptions are scored as side-by-side variants in memory. The `designated` variant stays the v4 routing (`docs/model_decision.md` section 5) until the user decides.
- Scripts read committed tables under `data/` only. A script that needs `extras/` is skip-guarded and reports `SKIPPED`.
- Readings and metric definitions are committed in `results/v5/BASE/README.md` in a commit that contains no result file (fixes B6 for this work).
- No download without the user's go-ahead. AJFTD raw rows are never published.
- Out of scope (user said "not yet"): `FULL_REPORT_v5.md`, paper draft, PR to `main`, `splits/locked.json` replacement.
- pandas 3: assign mixed-type values per column.
- Commit after each task. No attribution lines in commit messages.

## Review Focus

1. A fuel whose composition maps under 90 % of its mass: scored apart as `low_coverage`, never dropped silently and never in the headline. (Test in Task 3.)
2. A (fuel, property, T) with several measured values from different methods: one median target, and the method family recorded so flash and heat can be stratified. (Test in Task 3.)
3. A cell with fewer than 8 groups: the scoreboard prints the numbers with verdict `n_too_small`, not a skill verdict. (Test in Task 9.)
4. A property with no model or no data in a tier: the scoreboard carries an explicit `no_evaluation` row with the reason. A missing row is a test failure. (Test in Task 9.)
5. A variant (LFL overlay, K1 offset) that changes a designated number: it must appear only in the variant column; the designated column must reproduce the committed number to the stated tolerance. (Tests in Tasks 5, 6 and 9.)

## The scoreboard (what "full baseline" means)

Rows = property x tier x set x model x variant. Columns:
`property, tier, set, model, variant, status, n_groups, n_rows, bias, mae, rmse, unsafe_rate, trivial_model, mae_trivial, diff_trivial_minus_model, ci_lo, ci_hi, verdict, margin90, coverage90, n_excluded_censored, n_excluded_nonfinite, n_low_coverage, evidence_path, label`

- `tier` in `T1_lab | T2_measured | CRC647_r3 | T3_ajftd | diesel | designs`.
- `status` in `recomputed_v5 | aggregated_v5 | carried_v4 | no_evaluation`.
- `label` in `descriptive | dev | fold_b | already_read | exploratory`.
- `verdict` in `beats_trivial | no_skill_shown | worse_than_trivial | n_too_small | not_applicable`. Rule: paired group bootstrap (5000, seed 0) 95 % interval of `MAE(trivial) - MAE(model)`; above 0 = beats, below 0 = worse, else no skill shown; fewer than 8 groups = n too small.
- `trivial_model`: leave-one-group-out constant. For flash and heat of combustion, the constant is fitted within method family.

Cells that must exist (Task 9 asserts every one):

| Property | T1 lab mixtures | T2 measured composition | CRC 647 r3 (reconstructed) | T3 AJFTD bulk | Diesel |
|---|---|---|---|---|---|
| Viscosity, by T bin (<263, 263-283, >=283 K) | Task 4 (aggregate S1) | Task 6 (recompute, 9 T) | Task 8 (carry v4 R3) | Task 7 (aggregate S4, -20 C) | Task 8 (aggregate S6) |
| Density, by T | Task 4 (v4 L + ThermoML 2020 new systems) | Task 6 (recompute, 9 T) | Task 8 | no_evaluation (input, not target) | Task 8 |
| Flash point | Task 4 (v4 F, 126 rows) | Task 5 (recompute) | Task 8 | Task 7 | Task 8 |
| Freezing point | Task 4 (SLE2 tiers) | Task 6 (aggregate S5: C, I, F1, P) | Task 8 | Task 7 | no_evaluation (cloud point only) |
| Net heat of combustion | no_evaluation | Task 6 (new) | Task 8 | Task 7 (S4 + S4b) | no_evaluation |
| Smoke point | no_evaluation | Task 6 (new, YSI route) | no_evaluation | Task 8 (carry v4 SP) | no_evaluation |
| Distillation (T10, T50, T90) | no_evaluation | Task 6 (new) | no_evaluation (input) | no_evaluation (input) | no_evaluation |
| Pour point | no_evaluation | no_evaluation (no data) | Task 7 part 4 (bulk baseline, n = 57) | no_evaluation | no_evaluation |
| Cetane / DCN | no_evaluation (v4: uncertified, no model scored) | no_evaluation | no_evaluation | no_evaluation | Task 8 (S6: not evaluable) |
| RON / MON, DVPE | Task 8 (carry v4 O) | - | - | - | - |
| Designs certified | Task 8 (aggregate D4: 0 of 4) | | | | |

---

## Part A: one path for composition, prediction and statistics

### Task 0: Record the starting state (observed 2026-11-12 at `271fdcc`)

Observed while writing this plan, to be re-checked before Task 1:

- `.venv/Scripts/python.exe -m pytest tests -q --tb=no -p no:cacheprovider` -> `433 passed` (94 s).
- `.venv/Scripts/python.exe scripts/v4/reproduce_all.py --verify-only` -> **exit 1**: `[MANIFEST] DRIFT: 1 changed, 0 missing, 20 new`. The changed file is `docs/claims_register.md` (edited after `SHA256SUMS_v4.txt` was last written in `65212cc`); the new files are R3 outputs added after it. This drift predates this plan. Do not run `--update-manifest` for v4 here; it is listed in Task 17 for the user to approve. Because of this, no step in this plan uses `--verify-only` as a pass/fail check; v4 integrity is checked with `git diff` (see Verification).

- [ ] **Step 1:** Re-run both commands. Expected: `433 passed` and the same single CHANGED line. If either differs, stop and report before Task 1.
- [ ] **Step 2:** Record the commit hash: `git rev-parse HEAD > results/v5/BASE/START_COMMIT.txt` is written in Task 1 Step 1 together with the README (so Task 1's commit still holds no result file).

### Task 1: Declare the baseline (README only)

**Files:**
- Create: `results/v5/BASE/README.md`

- [ ] **Step 1: Write `results/v5/BASE/START_COMMIT.txt`** (output of `git rev-parse HEAD`) **and `results/v5/BASE/README.md`** with: status line `DECLARED 2026-11-12, not scored`; the scoreboard column definitions, verdict rule, unsafe-side definitions, group rule, censoring rule and `low_coverage` rule copied from "Global Constraints" and "The scoreboard" above; the cell table above; the variant list (`designated`, `lfl_overlay`, `mesitylene_fixed`, `lfl_overlay+mesitylene_fixed`, `k1_cold_offset_loso`, `sle2_two_compilations`, `sle2_single_compilation`); the sentence "Every T2 fuel was read before this declaration; the baseline is descriptive." and these five readings:
  - R-BASE-1: every cell in the table exists with a status.
  - R-BASE-2: each `recomputed_v5` cell that has a committed counterpart reproduces it within tolerance (flash 0.05 C, ln viscosity 0.005, density 0.5 kg/m3, freezing 0.05 C) or the difference is listed with its cause.
  - R-BASE-3: per property, name the tier where the designated model beats the trivial model and the tier where it does not.
  - R-BASE-4: per property, name the largest unsafe-side rate and the set it occurs on.
  - R-BASE-5: list every cell whose verdict is `n_too_small` or `no_evaluation`; this is the data-gap list.
- [ ] **Step 2: Commit** (this commit must contain only the README and the start-commit file)

```bash
git add results/v5/BASE/README.md results/v5/BASE/START_COMMIT.txt
git commit -m "v5 BASE: declare full-baseline scoreboard definitions and readings (no results)"
```

### Task 2: Shared composition route in `src/fuelblend`

**Files:**
- Create: `src/fuelblend/eval/composition_route.py`
- Test: `tests/test_base_composition_route.py`
- Not modified: `scripts/v4/build_viscosity_leaderboard_v4.py::refutas_nu` stays as it is (v4 evidence); the guard lives in the new `refutas_nu_safe`.

**How the pieces are built today (pin these, do not guess):**
- Context: `scripts/v5/s3_composition_vs_physics.py::main` lines 271-281: `reg = CompoundRegistry(ROOT/"data/unified/compounds.csv")`, `prov = PurePropertyProvider(reg, ROOT/"data/cache/pure_grid.npz", verbose=False)`, `reps = SubReps(reg, prov, vp_lookup=lambda s: _lookup(s, None))`. The test fixtures `route_ctx` and `reps` are built the same way (module-scoped).
- Flash: `fuelblend.baselines.flash_unifac.flash_point_unifac_detail(smiles, x_mole, registry=reg_df, use_unifac=False)`; the default table is `flash_unifac._default_registry()`.
- In-memory input overrides: copy that registry DataFrame and overwrite columns, as `scripts/v5/eval_flash_lfl_overlay_v5.py:35-46` does for `lfl_vol_pct` and `scripts/v5/s2_mesitylene.py:182-186` does for `antoine_A/B/C`; pass the copy through `registry=`. `RouteCtx` builds this copy once in `__post_init__`.
- `score_like_s3()` in the test: for each entity in `corrected_2026-11-12_fuels.csv` with `M_ok`, map its bins with `map_bins` and score the (property, T_K) keys present in `corrected_2026-11-12_predictions.csv`.

**Interfaces:**
- Consumes: `fuelblend.data.reconstruction_r3.SubReps`; `fuelblend.mixture.common.CompoundRegistry`, `PurePropertyProvider`; `map_fuel`, `hc_class`, `ideal_sle_C` bodies from `scripts/v5/eval_ajftd_gcxgc_v5.py`; `predict_comp`, `drop_zero_mass` bodies from `scripts/v5/s3_composition_vs_physics.py`.
- Produces:
  - `map_bins(bins: pd.DataFrame, reps: SubReps) -> MappedFuel` where `bins` has columns `hc_class, n_carbon, mass_pct` and `MappedFuel` is a dataclass `idx: np.ndarray[int64], w: np.ndarray[float] (sums to 1 over mapped mass), unmapped_fraction: float, unmapped_by_class: dict[str, float]`. Property-independent: a bin is unmapped only when `reps.get` returns `None`.
  - `CLASS_MAPS: dict[str, dict[str, str]]` with keys `"ajftd"`, `"av3322"`, `"fuellib"`, `"face_diesel"`, each mapping the source's class label to the `SubReps` kind.
  - `refutas_nu_safe(provider, idx, w, T_K) -> tuple[float, float]` returning `(nu_cSt, rho_kg_m3)`; drops `w <= 0`; returns `(nan, rho)` when any remaining component has `nu_i + 0.8 <= 1`.
  - `predict_properties(ctx: RouteCtx, fuel: MappedFuel, targets: pd.DataFrame) -> pd.DataFrame` adding `pred, status`; `targets` has `property, T_K`; properties handled: `flash_point_C`, `viscosity_kin_cSt`, `density_kg_m3`, `freezing_point_C`. Status `ok` implies finite.
  - `RouteCtx(reg, prov, tmap, lfl_override: dict[str, float] | None = None, antoine_override: dict[str, tuple[float, float, float]] | None = None)`; overrides are keyed by InChIKey and applied in memory only.
  - `COVERAGE_MIN = 0.90`.

- [ ] **Step 1: Write the failing tests**

```python
def test_zero_mass_rep_does_not_poison_refutas(route_ctx, propane_idx, decane_idx):
    nu, rho = refutas_nu_safe(route_ctx.prov, np.array([propane_idx, decane_idx]), np.array([0.0, 1.0]), 313.15)
    assert np.isfinite(nu) and np.isfinite(rho)

def test_status_ok_implies_finite(route_ctx, av3322_fuel_01):
    out = predict_properties(route_ctx, av3322_fuel_01, ALL_TARGETS_9T)
    assert np.isfinite(out.loc[out.status == "ok", "pred"]).all()

def test_map_bins_is_property_independent(reps, bins_with_diaromatic):
    m = map_bins(bins_with_diaromatic, reps)          # fixture: 3 % of mass in a diaromatic bin
    assert m.unmapped_by_class["aromatic_di"] == pytest.approx(0.03)
    assert m.w.sum() == pytest.approx(1.0)

def test_reproduces_s3_corrected_M_predictions():
    ref = pd.read_csv("results/v5/S3/corrected_2026-11-12_predictions.csv")
    new = score_like_s3()            # helper in the test: same fuels, same targets, new route
    j = ref.merge(new, on=["entity", "property", "T_K"])
    ok = j.status_M == "ok"
    assert np.allclose(j.pred_M[ok], j.pred[ok], rtol=1e-9, atol=1e-9)

def test_lfl_override_changes_only_flash(route_ctx, av3322_fuel_01):
    base = predict_properties(route_ctx, av3322_fuel_01, ALL_TARGETS_9T)
    over = predict_properties(replace(route_ctx, lfl_override=ONE_LFL_CHANGE), av3322_fuel_01, ALL_TARGETS_9T)
    nf = base.property != "flash_point_C"
    assert np.array_equal(base.pred[nf].values, over.pred[nf].values)
    assert base.pred[~nf].iloc[0] != over.pred[~nf].iloc[0]
```

- [ ] **Step 2: Run** `.venv/Scripts/python.exe -m pytest tests/test_base_composition_route.py -q` - expected: FAIL (module missing).
- [ ] **Step 3: Implement `composition_route.py`** by moving the function bodies named under "Consumes" (do not rewrite the science). If `test_reproduces_s3_corrected_M_predictions` cannot pass for AV-33-22 or FuelLib because S3 took their mapping from S2's `rep_for` (which marks bins unmapped for missing Antoine/LFL), do not bend the new route: write the differing rows to `results/v5/BASE/route_vs_s3_differences.csv` (columns `entity, property, T_K, pred_s3, pred_new, cause`) and narrow that test to AJFTD entities plus an assertion that every difference row has a non-empty `cause`.
- [ ] **Step 4: Run the tests** - expected: PASS.
- [ ] **Step 5: Run the full suite** `.venv/Scripts/python.exe -m pytest tests -q -p no:cacheprovider` - expected: 0 failed, at least 433 + the new tests passed. Record the count in the commit message.
- [ ] **Step 6: Commit** `git commit -m "v5 BASE: shared composition route (map_bins, guarded Refutas, in-memory input overrides)"`

### Task 3: Shared statistics and the T2 target table

**Files:**
- Create: `src/fuelblend/eval/group_stats.py`
- Create: `scripts/v5/base/build_t2_targets.py`
- Test: `tests/test_base_stats_targets.py`
- Output: `results/v5/BASE/t2_targets.csv`, `results/v5/BASE/t2_fuels.csv`

**Interfaces:**
- Produces (`group_stats.py`):
  - `paired_group_bootstrap(abs_err_ref, abs_err_model, groups, n_boot: int = 5000, seed: int = 0) -> tuple[float, float, float]` returning `(diff, lo, hi)` of `MAE(ref) - MAE(model)`.
  - `logo_constant(y, groups, strata=None) -> np.ndarray` leave-one-group-out mean, within `strata` when given.
  - `verdict(n_groups: int, lo: float, hi: float, min_groups: int = 8) -> str` returning the scoreboard verdict strings.
  - `cell_metrics(err, groups, unsafe: Callable[[np.ndarray], np.ndarray] | None) -> dict` with keys `n_groups, n_rows, bias, mae, rmse, unsafe_rate`.
- Produces (`build_t2_targets.py`):
  - `t2_targets.csv` columns: `entity, source, fuel_id, label, group, property, T_K, measured, n_values, method, method_family, qualifier_dropped, already_read_locked` - one row per (entity, property, T_K), median over values.
  - `t2_fuels.csv` columns: `entity, source, group, chosen_for_group, unmapped_fraction, low_coverage, n_components`.
  - Sources: `crc_av3322_properties.csv` (viscosity and density at all 9 temperatures, flash, freezing, `heat_of_combustion_net_MJ_kg`, `smoke_point_mm`, `vapor_pressure_kPa`, `distillation_T10_C/T50_C/T90_C`), `ajftd_fuels_long.csv` joined to `ajftd_gcxgc_fuels.csv` (flash, freezing, viscosity with `T_source == "reported"`, density, heat, plus smoke and D86 where present), `fuellib_properties.csv` (rows with finite `T_C` only; the 48 HEFA-blend rows with `T_C = NaN` are counted and skipped).
  - `method_family(property, method) -> str`: flash `D56 | D93_IP170 | coded_220x | other`; heat `measured (D4809, D240) | calculated (D3338) | unclassified`. Move `_flash_family` from `scripts/v5/eval_bulk_descriptor.py:66` and `classify_method` from `scripts/v5/s4b_heat_measured.py:64` into `group_stats.py` and re-export them from the scripts.

- [ ] **Step 1: Write the failing tests**

```python
def test_bootstrap_matches_s5_committed_row():
    # recompute S5 reading "F1 vs C, headline" from results/v5/S5/s5_predictions.csv
    d, lo, hi = paired_group_bootstrap(abs_c, abs_f1, groups)
    assert d == pytest.approx(2.47, abs=0.01)
    assert lo > 0

def test_logo_constant_never_sees_own_group():
    y = np.array([1.0, 1.0, 9.0]); g = np.array(["a", "a", "b"])
    assert np.allclose(logo_constant(y, g), [9.0, 9.0, 1.0])

def test_verdict_small_n():
    assert verdict(7, 0.5, 1.0) == "n_too_small"

def test_targets_exclude_censored_and_count_them():
    t = pd.read_csv("results/v5/BASE/t2_targets.csv")
    assert t[t.property == "freezing_point_C"].qualifier_dropped.sum() >= 5
    assert not ((t.entity == "AJFTD:293") & (t.property == "freezing_point_C")).any()

def test_av3322_has_nine_temperatures():
    t = pd.read_csv("results/v5/BASE/t2_targets.csv")
    a = t[(t.source == "av3322") & (t.property == "viscosity_kin_cSt")]
    assert a.groupby("entity").T_K.nunique().min() == 9

def test_one_median_target_per_key_with_method_family():          # Review Focus 2
    t = pd.read_csv("results/v5/BASE/t2_targets.csv")
    assert not t.duplicated(["entity", "property", "T_K"]).any()
    assert t[t.property.isin(["flash_point_C", "heat_of_combustion_net_MJ_kg"])].method_family.notna().all()

def test_one_chosen_fuel_per_group():
    f = pd.read_csv("results/v5/BASE/t2_fuels.csv")
    assert (f.groupby("group").chosen_for_group.sum() == 1).all()

def test_low_coverage_flag_uses_090():                            # Review Focus 1
    f = pd.read_csv("results/v5/BASE/t2_fuels.csv")
    assert (f.low_coverage == (f.unmapped_fraction > 0.10)).all()
```

If the S5 bootstrap used a different resampling unit and the interval cannot be matched, keep `paired_group_bootstrap` as specified; the first test asserts the point estimate within 0.01 and `lo > 0` only.

- [ ] **Step 2: Run** the tests - expected: FAIL.
- [ ] **Step 3: Implement** both files; run `build_t2_targets.py`.
- [ ] **Step 4: Run** the tests - expected: PASS.
- [ ] **Step 5: Commit** `git commit -m "v5 BASE: shared group statistics and unified T2 target table (9-temperature AV-33-22, censored excluded)"`

---

## Part B: the baseline

Every scoring script in Part B lives in `scripts/v5/base/`, writes `results/v5/BASE/cells_<name>.csv` in the scoreboard column schema plus a per-row file `results/v5/BASE/rows_<name>.csv`, and is tested in `tests/test_base_cells.py` (one test function per task, named `test_cells_<name>`).

### Task 4: T1 lab-mixture cells

**Files:**
- Create: `scripts/v5/base/cells_t1_lab.py`
- Test: `tests/test_base_cells.py::test_cells_t1_lab`

**Interfaces:**
- Consumes: `results/v5/S1/viscosity_lab_residuals_v5.csv` (`method == "refutas"`, `split_role`, `t_bin`), `results/v5/S1/uncertainty_v5.yaml`, `results/v4/L/density_lab_residuals.csv`, `results/v4/F/flash_lab_loso.csv`, `results/v5/L1/comparison.csv`, `results/v5/S2/mesitylene_flash_rerun.csv`, `results/v5/SLE2/summary_by_tier.csv`, `data/real_fuels/v5/thermoml2020_lab.csv`, `data/unified/measurements_dedup.csv`; `fuelblend.uq.conformal.two_sided_margin, one_sided_margin, coverage`.
- Produces: `cells_t1_lab.csv`, `rows_t1_density_thermoml2020.csv`.

Cells:
1. Viscosity, 3 temperature bins x (`dev`, `already_read` for locked): aggregate; margin and coverage from the S1 yaml. `status = aggregated_v5`. Trivial model: `not_applicable`.
2. Flash lab (31 systems, 126 rows): designated ideal; variants `lfl_overlay`, `mesitylene_fixed` read from the L1 and S2 files. `status = aggregated_v5`.
3. Freezing (SLE): three rows from `summary_by_tier.csv` as variants `designated` (v4 verified, margin inf), `sle2_two_compilations`, `sle2_single_compilation`.
4. Density lab, v4: aggregate `density_lab_residuals.csv` by `t_bin` and `split_role`.
5. **Density lab, new systems (recomputed).** Select `thermoml2020_lab.csv` rows with `property == "density_kg_m3"`, `kind != "pure"`, whose `system_id` is absent from `measurements_dedup.csv` and from `data/processed/v2/thermoml_clean_hc.csv`. Predict with ideal volume mixing from the provider's pure densities at the same T. These systems have never been scored for density, so this one cell is labelled `dev` with the note "first read"; before running, append to `results/v5/BASE/README.md` under a dated note the reading: "ideal-volume density on unseen systems: abs bias <= 2 kg/m3 and 90 % of systems inside the v4 lab density margin". Commit that note before Step 3 of this task.

- [ ] **Step 1: Write the failing test**

```python
def test_cells_t1_lab():
    c = pd.read_csv("results/v5/BASE/cells_t1_lab.csv")
    v = c[(c.property == "viscosity_kin_cSt") & (c.set == "T_lt_263K") & (c.label == "dev")].iloc[0]
    assert v.n_groups == 20 and v.bias == pytest.approx(-0.1135, abs=1e-3) and v.margin90 == pytest.approx(0.4154, abs=1e-3)
    f = c[(c.property == "flash_point_C") & (c.variant == "designated")].iloc[0]
    assert f.n_groups == 31 and f.n_rows == 126 and f.bias == pytest.approx(-3.30, abs=0.01)
    s = c[(c.property == "freezing_point_C") & (c.variant == "designated")].iloc[0]
    assert np.isinf(s.margin90) and s.n_groups == 8
    d = c[(c.property == "density_kg_m3") & (c.set == "thermoml2020_new_systems")]
    assert len(d) >= 1 and d.n_rows.sum() > 1000 and d.label.eq("dev").all()
```

- [ ] **Step 2: Run** - expected: FAIL.
- [ ] **Step 3:** Commit the dated README note, then implement and run the script.
- [ ] **Step 4: Run** the test - expected: PASS. If the new-systems density row count is below 1,000, print the count by drop reason and set the assertion to the observed count (the audit figure of 1,440 was not run through this filter).
- [ ] **Step 5: Commit** `git commit -m "v5 BASE: T1 lab cells (viscosity bins, flash, SLE tiers, density incl. unseen ThermoML 2020 systems)"`

### Task 5: T2 flash cells with input variants

**Files:**
- Create: `scripts/v5/base/cells_t2_flash.py`
- Test: `tests/test_base_cells.py::test_cells_t2_flash`

**Interfaces:**
- Consumes: Task 2 `map_bins, predict_properties, RouteCtx`; Task 3 tables and `group_stats`; `data/processed/v5/lfl_v5.csv` (overlay); `results/v5/S2/mesitylene_best_input.json` (corrected Antoine set); `results/v5/S2/t2_flash_predictions.csv` and `results/v5/G1/predictions.csv` (reproduction check only).
- Produces: `cells_t2_flash.csv`, `rows_t2_flash.csv` (columns `entity, source, fuel_id, group, method_family, measured, pred_designated, pred_lfl_overlay, pred_mesitylene_fixed, pred_both, status, low_coverage, already_read_locked`), `rows_t2_flash_outliers.csv`.

Sets: `pooled`, `av3322`, `ajftd_gcxgc`, `fuellib`, and each `method_family` on the pooled set. Variants: the four flash variants, plus one row with `model = api_as_shipped`: call `fuelblend.api.predict_blend.predict(composition=[(inchikey, mass_fraction), ...], T_K=298.15, fuel_type="jet", composition_kind="known")` on the mapped composition, take its flash value and interval, and report `margin90` (the applied margin the API used) and `coverage90` (fraction of groups whose measured flash is at or above the API's lower bound). If the API's flash differs from `pred_designated` by more than 0.05 C for any fuel, list those fuels in `rows_t2_flash_api_differences.csv`; do not reconcile them here. Unsafe rate: fraction with error > +5 C. Trivial: LOGO constant within method family. Write the four AJFTD fuels the audit names (labels 7695, 9324, 8092, 8438) to `rows_t2_flash_outliers.csv` with mass % below C8 by class; they stay in every statistic.

- [ ] **Step 1: Write the failing test**

```python
def test_cells_t2_flash():
    r = pd.read_csv("results/v5/BASE/rows_t2_flash.csv")
    g1 = pd.read_csv("results/v5/G1/predictions.csv").query("property == 'flash_point_C' and model == 'ideal'")
    j = r[r.source == "ajftd_gcxgc"].merge(g1, on="fuel_id")
    assert len(j) == 37 and np.allclose(j.pred_designated, j.pred, atol=0.05)          # R-BASE-2
    c = pd.read_csv("results/v5/BASE/cells_t2_flash.csv")
    assert set(c.variant) == {"designated", "lfl_overlay", "mesitylene_fixed", "lfl_overlay+mesitylene_fixed"}
    assert (c[c.n_groups < 8].verdict == "n_too_small").all()
    assert c.unsafe_rate.between(0, 1).all()
```

(Check the `model` label G1 uses for the ideal route in `predictions.csv` and use that exact string in the query.)

- [ ] **Step 2: Run** - FAIL.
- [ ] **Step 3: Implement and run.**
- [ ] **Step 4: Run** - PASS. If the 37-row reproduction fails because Task 2 recorded route differences, restrict the `allclose` to fuels absent from `route_vs_s3_differences.csv` and assert the remainder is listed there.
- [ ] **Step 5: Commit** `git commit -m "v5 BASE: T2 flash cells with LFL-overlay and mesitylene variants, method-stratified constant"`

### Task 6: T2 viscosity, density, freezing, heat, smoke, distillation cells

**Files:**
- Create: `scripts/v5/base/cells_t2_bulk.py`, `src/fuelblend/eval/cold_offset.py`
- Modify: `scripts/v5/c1_cold_viscosity_correction.py` (import `fit_k1, apply_k1` from the new module instead of defining them)
- Test: `tests/test_base_cells.py::test_cells_t2_bulk`

**Interfaces:**
- Consumes: Tasks 2-3; `results/v5/S5/s5_predictions.csv` (models `C, I, F1, P`); `fuelblend.baselines.heat_of_combustion_v2.mixture_heat_of_combustion_v2`; `fuelblend.baselines.sooting.ysi_mix, mixture_mw, smoke_point_from_ysi`; `fuelblend.baselines.distillation.predict_d86_from_composition`.
- Produces: `cold_offset.fit_k1(T, resid) -> dict[float, float]`, `cold_offset.apply_k1(k1, T) -> np.ndarray` (bodies moved unchanged from the C1 script); `cells_t2_bulk.csv` (scoreboard columns plus `T_K`, NaN for properties without temperature); `rows_t2_bulk.csv`; `viscosity_error_vs_T.csv` (columns `source, T_K, n_groups, bias_ln, mae_ln, sd_ln`); `density_error_vs_T.csv` (same, in kg/m3).

Cells:
1. Viscosity at each measured temperature and per S1 bin, sets `pooled, av3322, ajftd_gcxgc, fuellib`. Temperatures are snapped to the nine AV-33-22 temperatures within 0.3 K (the C1 `snap` rule) before grouping, so AJFTD rows at 253.2 K land on 253.15. Rows:
   - `designated` (raw Refutas), with `margin90`/`coverage90` from the shipped `configs/uncertainty_v4.yaml` lab keys (`viscosity_kin_cSt.lab_mixtures.raw_refutas.{bin}`; `inf` below 283 K means coverage is `not_applicable`), and a second pair of columns `margin90_s1`, `coverage90_s1` from `results/v5/S1/uncertainty_v5.yaml` (0.4154 / 0.6084), which is not adopted in `configs/`.
   - `model = api_as_shipped`, `variant = known`: `predict(..., composition_kind="known")` at each temperature, with the API's own interval coverage.
   - `model = api_as_shipped`, `variant = reconstructed_route`: `predict(..., composition_kind="reconstructed")`, i.e. recalibrated Refutas with the CRC 647 margin, applied to measured compositions; this shows what the CRC-fitted recalibration does on fuels it was not fitted on.
   - `variant = k1_cold_offset_loso` (offset fitted with the scored fuel's source held out: AV-33-22 offsets applied to AJFTD and the reverse; FuelLib gets the AV-33-22 offsets).
   - `model = structured_ridge`, `status = no_evaluation`, reason "trained inside `build_viscosity_leaderboard_v4.py::main` on CRC 647 fuels; no saved model to apply to T2 (Part C Task 16f)". Same for `bilinear_pairwise` and `structured_gp`.
   Unsafe: under-prediction beyond the S1 margin.
2. Density at each temperature, same sets, `designated` (ideal volume) plus `api_as_shipped` with the API's side-specific lab density margins. Trivial: LOGO constant per temperature.
3. Freezing: aggregate `s5_predictions.csv` headline rows (do not refit); models `C` (trivial), `I` (designated: ideal SvL), `F1`, `P` as variants `f1_recalibrated`, `p_ridge`. Locked AJFTD rows reported as `already_read`.
4. Net heat of combustion (new): `mixture_heat_of_combustion_v2(w_mass, comp_gcv_MJ_kg, comp_mw_g_mol, comp_n_carbon=None, use_standard_jet_h=True) -> dict` on the mapped composition. Component gross heats come from `GroupAdditivityGCVModel` in the same module, built the way `scripts/v4/eval_bulk_props_v4.py` builds it for `hoc_v4.csv` (copy that call; the two columns there, `ncv_v2_const_H` and `ncv_comp_H`, become models `ga_const_H` and `ga_comp_H`; `ga_comp_H` is the designated one). Stratified by `method_family`; headline on `measured` rows; AV-33-22 CoA rows as their own set.
5. Smoke point (new, exploratory): `ysi_mix(smiles, x_mole)` then `smoke_point_from_ysi(ysi, mixture_mw(smiles, x_mole))`; `label = exploratory`; exclude and count fuels whose YSI-table mass coverage is under 0.90.
6. Distillation T10, T50, T90 (new, exploratory): `predict_d86_from_composition(boiling_points_C, fractions)`. Boiling point per representative = the temperature where the registry Antoine set gives 101.325 kPa (solve with `scipy.optimize.brentq` on 250-800 K; a representative with no Antoine set is dropped and its mass counted; fuels under 0.90 mass coverage are excluded and counted). `fractions` are volume fractions from mass fractions and the provider's 15 C densities. `label = exploratory`.

Reproduction fallback (same rule as Tasks 2 and 5): the C1 numbers asserted below came from S3's mapping. If Task 2 wrote `route_vs_s3_differences.csv` rows for AV-33-22, the two `av3322` bias assertions use tolerance 0.05 instead of 0.02, and the observed values are written to `results/v5/BASE/route_vs_c1_bias.csv` (`T_K, bias_c1, bias_base`) so the shift is on record.

For cells 4-6 append a dated note with these readings to the BASE README and commit it before scoring: "heat: group additivity beats a within-method constant on measured rows (n >= 8 groups)"; "smoke and distillation: reported, no skill claim; the verdict comes from the rule".

- [ ] **Step 1: Write the failing test**

```python
def test_cells_t2_bulk():
    vt = pd.read_csv("results/v5/BASE/viscosity_error_vs_T.csv")
    a = vt[vt.source == "av3322"].set_index("T_K")
    assert len(a) == 9
    assert a.loc[233.15, "bias_ln"] == pytest.approx(-0.637, abs=0.02)      # C1 R-C1-5
    assert a.loc[313.15, "bias_ln"] == pytest.approx(0.066, abs=0.02)
    c = pd.read_csv("results/v5/BASE/cells_t2_bulk.csv")
    k1 = c[(c.property == "viscosity_kin_cSt") & (c.variant == "k1_cold_offset_loso") & (c.set == "ajftd_gcxgc") & (c.T_K == 253.15)]
    assert abs(k1.bias.iloc[0]) <= 0.10                                     # C1 R-C1-1
    fz = c[(c.property == "freezing_point_C") & (c.set == "pooled")]
    assert fz[fz.variant == "f1_recalibrated"].mae.iloc[0] == pytest.approx(3.46, abs=0.02)
    for p in ["heat_of_combustion_net_MJ_kg", "smoke_point_mm", "distillation_T10_C", "density_kg_m3"]:
        assert (c.property == p).any()
    r = pd.read_csv("results/v5/BASE/rows_t2_bulk.csv")
    assert np.isfinite(r.loc[r.status == "ok", "pred"]).all()
```

- [ ] **Step 2: Run** - FAIL.
- [ ] **Step 3:** Commit the README notes, implement, run; then run `.venv/Scripts/python.exe -m pytest tests/test_c1_cold_viscosity.py -q` - expected: PASS (the move did not change C1).
- [ ] **Step 4: Run** the new test - PASS.
- [ ] **Step 5: Commit** `git commit -m "v5 BASE: T2 viscosity and density at 9 temperatures, freezing, heat, smoke, distillation cells"`

### Task 7: T3 AJFTD bulk cells and the pour-point baseline

**Files:**
- Create: `scripts/v5/base/cells_t3_bulk.py`
- Test: `tests/test_base_cells.py::test_cells_t3_bulk`

**Interfaces:**
- Consumes: `results/v5/S4/dev_scores.csv`, `locked_scores.csv`, `dev_predictions.csv`; `results/v5/S4b/s4b_metrics.csv`, `s4b_bootstrap.csv`; `data/real_fuels/v4/jet_properties_v4.csv`; Task 3 `group_stats`.
- Produces: `cells_t3_bulk.csv`, `rows_pour_point.csv`.

Cells:
1. Flash, freezing, viscosity at -20 C, heat: `dev / logo` rows aggregated from `dev_scores.csv` for models `const, ridge, hgb`, and the two flash correlations with `variant = unverified_coefficients` (B4). Locked rows aggregated from `locked_scores.csv` with `label = already_read`. Never call `cmd_locked`.
2. Flash with a method-stratified constant (dev only, new): from `dev_predictions.csv` (`scheme == "logo"`), recompute the constant with `logo_constant(strata=flash_family)` and re-run the verdict for `ridge` and `hgb` against it; `trivial_model = logo_constant_within_method`.
3. Heat on measured rows: aggregate S4b.
4. Pour point (new, n = 57, `label = dev`, CRC 647 only). Target `pour_point_C` from `jet_properties_v4.csv`. Models, leave-one-out: `const`; `ridge` on D86 T10/T50/T90, aromatics vol %, density (alpha by inner leave-one-out over `[0.1, 1, 10]`); `freezing_offset` (measured freezing minus the leave-one-out mean of freezing-minus-pour; it uses another measured property, so `variant = uses_measured_freezing` and `verdict = not_applicable`). The PPBI rule cannot run: the registry holds no pure-component pour points; say so in the T1/T2 `no_evaluation` rows.

- [ ] **Step 1: Write the failing test**

```python
def test_cells_t3_bulk():
    c = pd.read_csv("results/v5/BASE/cells_t3_bulk.csv")
    f = c[(c.property == "flash_point_C") & (c.model == "const") & (c.label == "dev") & (c.variant == "designated")].iloc[0]
    assert f.n_rows == 1728 and f.mae == pytest.approx(3.943, abs=1e-3)
    assert (c[c.label == "already_read"].status == "aggregated_v5").all()
    assert ((c.property == "flash_point_C") & (c.trivial_model == "logo_constant_within_method")).any()
    p = c[c.property == "pour_point_C"]
    assert set(p.model) >= {"const", "ridge", "freezing_offset"} and (p.n_rows == 57).all()
```

- [ ] **Step 2: Run** - FAIL.
- [ ] **Step 3: Implement and run.**
- [ ] **Step 4: Run** - PASS.
- [ ] **Step 5: Commit** `git commit -m "v5 BASE: T3 AJFTD cells (aggregated, method-stratified flash constant) and CRC 647 pour-point baseline"`

### Task 8: Carried cells (CRC 647 r3, diesel, designs, octane, smoke, cetane)

**Files:**
- Create: `scripts/v5/base/cells_carried.py`
- Test: `tests/test_base_cells.py::test_cells_carried`

**Interfaces:**
- Consumes: `results/v4/R3/F/flash_real_fuels.csv`, `results/v4/R3/B/{density_v4,hoc_v4,freezing_v4}.csv`, `results/v4/R3/V/viscosity_predictions_v4.csv`, `results/v5/L1/comparison.csv`, `results/v5/S6/summary_R-S6-2.csv`, `results/v5/S6/surrogate_viscosity_R-S6-3.csv`, `results/v5/D4/designs_v4.csv`, `results/v4/SP/smoke_point_eval.csv`, `results/v4/O/octane_eval_v4.csv`, `configs/uncertainty_v4.yaml`.
- Produces: `cells_carried.csv`.

Rules: CRC 647 rows are recomputed from the per-fuel files with `group_stats` (so they carry every scoreboard column), split by `split_role` into `dev` (45 fuels) and `fold_b` (12, the former locked set), `status = carried_v4`. Viscosity models carried: every `pred_*` column of `viscosity_predictions_v4.csv` (`Raw Refutas`, `Recalibrated Refutas` (designated), `Walther (D341)`, `Walther + recalibration`, `Bilinear pairwise`, `Structured ridge`, `Structured GP`, `Lab-trained structured residual (zero-shot)`), each per temperature (-40, -20, 20 C) and overall, with the designated row's margin and coverage from `viscosity_kin_cSt.crc647.recalibrated_refutas.all_T`. Flash gets variant `lfl_overlay` from L1. The mesitylene-corrected set was never run on CRC 647: add a `no_evaluation` row with reason "mesitylene fix not run on CRC 647 (Part C Task 10)". Diesel rows are `n_too_small` (n = 4). Designs: one row, `property = designs_certified`, `n_rows = 4`, `bias = 0` (number certified), evidence `results/v5/D4`.

- [ ] **Step 1: Write the failing test**

```python
def test_cells_carried():
    c = pd.read_csv("results/v5/BASE/cells_carried.csv")
    fl = c[(c.property == "flash_point_C") & (c.set == "crc647_r3") & (c.label == "dev") & (c.variant == "designated")].iloc[0]
    assert fl.n_groups == 45 and fl.bias == pytest.approx(0.44, abs=0.02) and fl.unsafe_rate == pytest.approx(0.133, abs=0.005)
    assert c[c.set == "crc647_r3"].label.isin(["dev", "fold_b"]).all()
    assert (c[c.set == "face_diesel"].verdict == "n_too_small").all()
```

- [ ] **Step 2: Run** - FAIL.
- [ ] **Step 3: Implement and run.**
- [ ] **Step 4: Run** - PASS.
- [ ] **Step 5: Commit** `git commit -m "v5 BASE: carried cells (CRC 647 r3 dev and fold_b, diesel, designs, octane, smoke)"`

### Task 9: Assemble, verify, freeze

**Files:**
- Create: `scripts/v5/base/build_scoreboard.py`, `scripts/v5/base/reproduce_base.py`
- Output: `results/v5/BASE/scoreboard.csv`, `results/v5/BASE/SCOREBOARD.md`, `results/v5/BASE/data_gaps.csv`, `SHA256SUMS_v5_base.txt`
- Test: `tests/test_base_scoreboard.py`

**Interfaces:**
- Consumes: every `cells_*.csv`; `compute_hashes, diff_manifest, read_manifest, write_manifest` imported from `scripts/v4/reproduce_all.py`.
- Produces:
  - `REQUIRED_CELLS: list[tuple[str, str]]` in `build_scoreboard.py`, the `(property, tier)` pairs transcribed from the cell table in this plan (every table cell, including the `no_evaluation` ones).
  - `scoreboard.csv` (all cells, sorted by property, tier, set, model, variant).
  - `SCOREBOARD.md`: one table per property with the designated row first, then variants; then the five R-BASE readings answered from the data; then "Decisions waiting on the user" (LFL overlay, mesitylene set, K1 offset, SLE2 tier), each with its designated and variant numbers side by side.
  - `data_gaps.csv`: every `no_evaluation` and `n_too_small` cell with `reason` and `what_would_fill_it` (a Part C task number or a user action U1-U8).
  - `reproduce_base.py [--update-manifest]`: runs the Task 3-9 scripts in order, then pytest on `tests/test_base_*.py`, then diffs `results/v5/BASE/**` against `SHA256SUMS_v5_base.txt`; exits 1 on drift; only `--update-manifest` writes the manifest.

- [ ] **Step 1: Write the failing tests**

```python
def test_every_required_cell_exists():                      # Review Focus 4, R-BASE-1
    s = pd.read_csv("results/v5/BASE/scoreboard.csv")
    have = set(zip(s.property, s.tier))
    assert not [c for c in REQUIRED_CELLS if c not in have]

def test_no_evaluation_rows_have_reasons():
    g = pd.read_csv("results/v5/BASE/data_gaps.csv")
    assert g.reason.notna().all() and g.what_would_fill_it.notna().all()

def test_small_n_never_gets_skill_verdict():                # Review Focus 3
    s = pd.read_csv("results/v5/BASE/scoreboard.csv")
    assert s[(s.n_groups < 8) & s.verdict.isin(["beats_trivial", "worse_than_trivial", "no_skill_shown"])].empty

def test_designated_columns_untouched_by_variants():        # Review Focus 5
    s = pd.read_csv("results/v5/BASE/scoreboard.csv")
    d = s[(s.property == "flash_point_C") & (s.set == "ajftd_gcxgc") & (s.variant == "designated")].iloc[0]
    assert d.bias == pytest.approx(-0.43, abs=0.05) and d.mae == pytest.approx(3.54, abs=0.05)

def test_tampered_file_makes_reproduce_exit_1(tmp_path):    # same pattern as tests/test_reproduce_drift.py
    ...  # copy BASE to tmp_path, write manifest, edit one byte, assert check_manifest(...) is True
```

Note on the last-but-one test: G1's n = 37 includes locked and duplicate fuels; if the BASE `ajftd_gcxgc` headline set (one per group) differs, add a set `ajftd_gcxgc_g1_definition` that uses G1's 37 fuels and point the test at it.

- [ ] **Step 2: Run** - FAIL.
- [ ] **Step 3: Implement; run `reproduce_base.py --update-manifest` once.**
- [ ] **Step 4: Run** `.venv/Scripts/python.exe -m pytest tests -q -p no:cacheprovider` - expected: 0 failed; then `.venv/Scripts/python.exe scripts/v5/base/reproduce_base.py` - expected: exit 0; then `git status --porcelain` - expected: changes only under `results/v5/BASE/`, `scripts/v5/base/`, `tests/` and `SHA256SUMS_v5_base.txt`.
- [ ] **Step 5: Append a dated `SCORED` note to the BASE README** (the declaration text stays as committed in Task 1) and answer R-BASE-1 to R-BASE-5 in `SCOREBOARD.md`.
- [ ] **Step 6: Commit** `git commit -m "v5 BASE: full-baseline scoreboard, data-gap list, drift-checked manifest"`
- [ ] **Step 7: Stop and report to the user.** Show `SCOREBOARD.md` and ask the four adoption decisions (U1). Record the Task 9 commit hash in a dated README note as the frozen reference. Part C starts only after that.

---

## Part C: improvements, each measured against the frozen scoreboard

Shared rule for every Part C task: declare the reading in `results/v5/BASE/delta_<name>/README.md` and commit it alone; run; write `delta.csv` with the scoreboard columns plus `baseline_mae, baseline_bias, baseline_unsafe_rate`; never edit `scoreboard.csv`. A change becomes the designated model only by the user's decision. Tasks are ordered by expected value over effort and are independent unless a dependency is named.

### Task 10: Flash input package (uses the user's decision U1)

**Files:** Create `scripts/v5/base/delta_flash_inputs.py`, `results/v5/BASE/delta_flash_inputs/`; Test `tests/test_base_delta_flash_inputs.py`.

Run the mesitylene-corrected set and the LFL overlay, separately and together, on CRC 647 r3 (dev and fold_b), lab mixtures and T2, through `RouteCtx` overrides. Recompute the lab flash margin in memory with `one_sided_margin` and report it beside the v4 value 3.70 C, with and without the mesitylene system. Reading: "with both inputs corrected, no set's unsafe rate rises, and the leave-mesitylene-system-out margin is within 0.5 C of the full margin". If the user adopts: a follow-up commit writes the flash margin once into a new `configs/uncertainty_v5.yaml` (v4 yaml untouched), so flash numbers change a single time.

Tests: `test_v4_designated_reproduced` (CRC 647 dev bias +0.44 with no overrides); `test_overrides_do_not_touch_registry_files` (sha256 of `data/processed/v5/pure_props_v5.csv` and `data/processed/v5/lfl_v5.csv` unchanged after the run).

### Task 11: Cold-viscosity correction with its own margin (depends on Task 6)

**Files:** Create `scripts/v5/base/delta_cold_viscosity.py`, `results/v5/BASE/delta_cold_viscosity/`; Test `tests/test_base_delta_cold_viscosity.py`.

C1 showed K1 works but gave it no margin. Cross-fit: for each group, fit K1 with that group's source held out, take the corrected residual, and compute a two-sided conformal margin per temperature on those out-of-source residuals with `two_sided_margin` (AV-33-22 35 groups at -20 and -40 C; AJFTD 19 and 8). Write `proposed_margins.yaml` in the delta folder (not `configs/`). Reading: "the corrected-model margin at -20 C is below the raw S1 margin 0.4154 and 90 % coverage holds on the other source".

Then the cause check (audit proposal 2): for every representative with more than 2 % mean mass, tabulate whether the provider is inside its measured range at 253 and 233 K and its error against Boehm pure data where present (the `c1_pure_check.csv` logic, extended from 5 compounds to all); output `pure_cold_extrapolation.csv` with `smiles, mean_mass_pct, data_T_min_K, source_at_253K (fit|ml_extrapolation), source_at_233K, err_vs_boehm_253K, err_vs_boehm_233K`. If more than half of the mass is on `ml_extrapolation`, the README states that the fix belongs in the pure provider and names Walther extrapolation from each compound's own measured range as the next delta.

Tests: the offset applied to a group never used that group's source; margin is `inf` when groups < 9.

### Task 12: Better representatives (unblocks diesel)

**Files:** Create `src/fuelblend/data/representatives_v5.py`, `scripts/v5/base/delta_representatives.py`, `results/v5/BASE/delta_representatives/`; Test `tests/test_base_delta_representatives.py`.

**Interfaces:** `class SubRepsV5(SubReps)` with the same `get(kind, n_carbon, max_gap=3)`; `REPLACED_BINS: set[tuple[str, int]]`; `candidate_table(reg, prov) -> pd.DataFrame` with columns `smiles, kind, n_carbon, has_antoine, has_lfl, has_viscosity, has_density, has_fusion, admitted, reason`.

Candidates (audit proposal 8): 1-methylnaphthalene and 2,6-dimethylnaphthalene (diaromatics); perhydrophenanthrene or perhydrofluorene (tricycloparaffins); alkyl-tetralins and alkyl-decalins at C12-C14; a C9 aromatic other than mesitylene (1,2,4-trimethylbenzene, n-propylbenzene). A candidate is admitted only if the registry already holds verified Antoine, LFL and eta(T)/rho(T) for it. Re-score Tasks 5 and 6 and the S6 diesel cells with `SubRepsV5`. Reading: "unmapped mass falls below 10 % for at least 7 of 9 FACE diesels, and T2 viscosity at 40 C and density bias do not worsen by more than 0.02 ln / 2 kg/m3". Candidates with missing constants become `data_gaps.csv`-style rows in the delta folder naming the source (ThermoML on disk first; fusion data needs U6).

Tests: `SubRepsV5.get` equals `SubReps.get` for every bin not in `REPLACED_BINS`; a candidate lacking a verified constant has `admitted == False`.

### Task 13: Pooled T2 residual model, leave-source-out

**Files:** Create `src/fuelblend/mixture/t2_residual.py`, `scripts/v5/base/delta_t2_residual.py`, `results/v5/BASE/delta_t2_residual/`; Test `tests/test_base_delta_t2_residual.py`.

**Interfaces:** `FEATURE_NAMES: list[str]` (20 names, fixed order); `bin_features(bins: pd.DataFrame) -> np.ndarray` (10 class mass fractions; mean and sd of carbon number overall and for n-, iso-, cyclo-paraffins and aromatics as four groups; mass fraction below C8; mass fraction above C15); `fit_predict_residual(X_train, r_train, g_train, X_test, alpha_grid=(0.1, 1, 10, 100)) -> np.ndarray` (ridge on features standardised on the training rows, alpha by inner leave-one-group-out).

This is the "physics baseline + learned residual" on measured compositions, which so far exists only on reconstructed CRC 647 fuels. Targets: residual of the designated physics model for flash, density at 15 C, ln viscosity at 40, -20 and -40 C, freezing, heat. Validation: train on AV-33-22 and test on AJFTD GCxGC, and the reverse; FuelLib is test-only. `label = exploratory`. Reading per property: "the residual model beats raw physics in both directions with the bootstrap interval above 0 and does not raise the unsafe rate". Low capacity is deliberate: about 79 labelled fuels for flash and fewer elsewhere, so no deep model here and no GPU.

Tests: the test source's rows never enter the scaler or the alpha search; a toy set with a constant residual is recovered within 1e-6.

### Task 14: Freezing follow-ups

**Files:** Create `scripts/v5/base/delta_freezing.py`, `results/v5/BASE/delta_freezing/`; Test `tests/test_base_delta_freezing.py`.

(a) F1 with the complete Coutinho n-alkane table (`F1x` columns already in `s5_predictions.csv`; sensitivity MAE 3.07 C) as a declared variant, scored with `group_stats`. (b) `data/processed/v2/wax_dauphin_coutinho.csv` (57 wax-appearance temperatures) as a second T1 SLE source under the SLE2 constants, method label `wax_appearance`, reported apart from liquidus data. (c) The unreported lead from the audit: the D86-reconstructed route gave freezing MAE 3.5 C against 5.9 C for the constant on 34 AV-33-22 fuels; test whether D86 T90/FBP plus measured n-paraffin share explains it (ridge, LOGO). Reading: "F1x beats the constant on AV-33-22 and on the pooled set; AJFTD stays not evaluable below 20 groups".

Tests: wax-appearance rows never enter a cell labelled liquidus; (a) reproduces 3.07 within 0.02.

### Task 15: A fresh holdout (the only non-descriptive T2 number available without new measurements)

**Files:** Create `src/fuelblend/data/loaders_nistir6659.py`, `scripts/v5/base/delta_fresh_holdout.py`, `results/v5/BASE/delta_fresh_holdout/`; Test `tests/test_nistir6659_loader.py`, `tests/test_base_delta_fresh_holdout.py`.

**Interfaces:** `load_nistir6659(path: Path) -> pd.DataFrame` with columns `fuel_id, posf, property, value, T_K, P_kPa, source_file, provenance`; writes `data/real_fuels/v5/nistir6659_density.csv`.

Step order is binding. (1) Declare readings and name the Task 9 commit hash as the frozen model version. (2) Load the POSF 4658 rows (121) from `extras/fuelblend_data/data/processed/nist_density_compressed_liquid.csv` (skip-guarded if `extras/` is absent); keep P <= 1000 kPa for the ideal-volume comparison. (3) List the AJFTD GCxGC fuels with a complete composition and no property row, by POSF label, and search only tables already under `data/` for properties with the same POSF label. (4) Check every found fuel against the double-count groups. (5) Score once. Reading: "density bias on the fresh fuels lies inside the T2 descriptive 95 % interval". The task runs on what is on disk and reports its n; U3 and U5 would enlarge it.

Tests: loader fixture built from 10 real rows of the source file with a unit check (kg/m3, K, kPa); the scoring script refuses to run a second time once `scored_once.json` exists.

### Task 16: Smaller data wires (each about half a day, independent)

| Sub-task | Data | Output | Blocker |
|---|---|---|---|
| 16a | AV-33-22 vapour pressure (266 rows) against Raoult + Antoine on the mapped composition | `delta_vapor_pressure/`: separates Psat error from LFL error in the flash model | none |
| 16b | Leenhouts published splits (`splits.zip`) for lab flash and viscosity | `delta_leenhouts_splits/`: our numbers on the benchmark's own splits | file is in `extras/`; skip-guarded |
| 16c | FuelLib 18 Jet-A/HEFA blends (48 rows with `T_C = NaN`) | `delta_fuellib_blends/`: real-fuel blending-rule test | U7 |
| 16d | Supervised reconstruction: map (D86, simdis, aromatics, density) to light-end fraction and n-paraffin share on 71 fuels, leave-source-out, then re-score CRC 647 dev | `delta_supervised_recon/` | none; about 2 days |
| 16e | Check the two S4 flash correlation forms against their sources (B4) | dated note in `results/v5/S4/README.md`; variant renamed from `unverified_coefficients` | source documents |
| 16f | Refit the v4 structured ridge and bilinear viscosity models on all CRC 647 dev fuels, save them, and apply them to T2 mapped compositions (transfer test of the learned models) | `delta_learned_viscosity_transfer/`: fills the three `no_evaluation` rows from Task 6 | needs the fit code lifted out of `build_viscosity_leaderboard_v4.py::main` into a function; about 1 day |

Each sub-task follows the shared Part C rule and adds one test that its designated-model numbers equal the frozen scoreboard's.

### Task 17: Bookkeeping owed (B5 and handoff section 4)

**Files:** Modify `docs/data_sources.csv`, `docs/claims_register.md`, `AGENTS.md` sections 8 and 12.

- [ ] Add `data_sources.csv` rows for the AJFTD GCxGC table, the LFL overlay and the extended SLE constants; remove `fill_lfl.py` as a user of `lfl_chemicals.csv`.
- [ ] Add one claim per R-BASE reading and per adopted delta to `claims_register.md` (new CLM numbers after CLM-41; retire superseded claims, do not edit them), each citing its scoreboard row.
- [ ] Update `AGENTS.md` 8.0 to point at `results/v5/BASE/SCOREBOARD.md` and this plan; correct "81 fuels" in 8.1.3.
- [ ] Ask the user before this step: the v4 manifest has drifted since `65212cc` (Task 0). With their go-ahead, run `.venv/Scripts/python.exe scripts/v4/reproduce_all.py --verify-only --update-manifest` after the claims register edit, so `SHA256SUMS_v4.txt` matches the tree again; say in the commit message which files were added to it.
- [ ] Commit `git commit -m "docs: v5 baseline claims, data-source rows, AGENTS.md refresh"`

## User actions (none of these is done by the agent)

| # | Action | Unblocks |
|---|---|---|
| U1 | Decide after Task 9: adopt the LFL overlay, the corrected mesitylene set, the K1 cold offset, the SLE2 two-compilations tier | Task 10, a design re-run with real margins |
| U2 | Browser download: Li et al. 2014, doi 10.1155/2014/239023 (n-alkane binary and ternary flash points) | more T1 flash systems (31 today) |
| U3 | Browser download: Puetz et al. 2025 Supporting Information E (8 fuels incl. 4 IsoPar solvents, GCxGC and properties) | Task 15 holdout; iso-paraffin representatives |
| U4 | Library access: Mueller et al. 2016, doi 10.1021/acs.energyfuels.5b02879 | diesel surrogate flash, density, cloud point (n = 4) |
| U5 | DTIC browser download: AFRL-RQ-WP-TR-2017-0091 | measured composition for the 57 CRC 647 fuels (moves them from reconstructed to T2) |
| U6 | Library access: Acree-Chickos fusion compendia | 9 missing enthalpies of fusion (21 blocked SLE systems), dimethylnaphthalene representatives |
| U7 | Confirm FuelLib HEFA-blend measurement temperatures from the Vozka paper | Task 16c |
| U8 | Delete `extras/ajftd_session.json` and `extras/fuelblend_data/ajftd_session.json`; move `extras/` out of OneDrive | Phase 0 leftover |

## Verification

- `.venv/Scripts/python.exe -m pytest tests -q -p no:cacheprovider` reports 0 failed after every task.
- v4 evidence untouched: `git diff --stat $(cat results/v5/BASE/START_COMMIT.txt)..HEAD -- results/v4 configs splits data/real_fuels/v4 scripts/v4` prints nothing after Task 9. (`reproduce_all.py --verify-only` is not used as a check: it already exits 1 at the starting commit, see Task 0.)
- `.venv/Scripts/python.exe scripts/v5/base/reproduce_base.py` exits 0 on a clean checkout and 1 after a hand edit of any file in `results/v5/BASE/`.
- `git log --diff-filter=A --format=%H -- results/v5/BASE/README.md` returns a commit older than the first commit that adds any `cells_*.csv`.
