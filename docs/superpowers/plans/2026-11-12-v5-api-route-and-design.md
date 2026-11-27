# v5 Part C delta: API on the scoring route, any-size blends, inverse design on a user pool — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** After this work, `predict()` gives physics-based viscosity, density, flash, freezing, net heat, smoke point and aromatics for a blend of any size (50+ components).
* It returns exactly the numbers the frozen BASE scoreboard scored.
* It reports per-property coverage.
* No fallback depends on temperature.

`design()` accepts a pool given as SMILES or InChIKeys, and its recipe re-predicts identically through `predict()`.

**Architecture:** A new `src/fuelblend/api/predict_route.py` reuses the BASE scoring machinery:
* `CompoundRegistry`;
* a read-only `PurePropertyProvider` over the frozen `data/cache/pure_grid.npz`;
* `SubReps`;
* the flash registry `compounds_anchored_vp.csv` with `flash_point_unifac_detail`;
* `composition_route.refutas_nu_safe` and `ideal_sle_C`.

Each property uses a component's own data where it exists. Otherwise it substitutes the SubReps representative of the same subkind and carbon number, explicitly and counted. Anything else is unmapped and reported. `predict_blend.predict` becomes a thin wrapper over the new module. Inverse design gets `RouteLibrary(Library)`, which builds constants from the same sources; the engine is unchanged.

**Tech Stack:** Python 3.12, numpy, pandas, RDKit, scipy, pytest (venv `.venv\Scripts\python.exe`).

**Spec:** the Fable decision memo of 2026-11-12 (summarised below) and the user's request of 2026-11-12:

> "add the fix and make sure the inverse works as well … fix the problems 1 2 3 or at least make it predict properly"

The three problems:
1. All-or-nothing physics with a silent fallback and no uncertainty.
2. A flash fallback that depends on the request temperature.
3. Only 3 properties.

The frozen reference is commit `5a62e20` (`results/v5/BASE`). The parent plan is `docs/superpowers/plans/2026-11-12-v5-full-baseline.md`; Part C deltas live in `results/v5/BASE/delta_<name>/` (its line 498).

## Decisions taken (Fable memo; spot-checks verified)

**D1. Re-route the API onto the scoring route.**
* A component missing pure data for a property is substituted by `SubReps.get(subkind, nC, max_gap=99 if subkind in NEAREST_ANY else 3)`, for that property only.
* Non-hydrocarbon, unknown subkind, or no representative → unmapped. If unmapped mass ≤ 0.10, drop and renormalise. Above that, value `None` and status `low_coverage`.
* Unparseable SMILES → `ValueError`.
* Lab yaml margins apply only when substituted mass = unmapped mass = 0.

**D2. Delete the `gbm_mix` flash fallback and `_flash_lechatelier`.** If flash still fails after substitution, return value `None`, method `none`, and a reason.

**D3. Flash = `flash_point_unifac_detail(..., registry=ctx.flash_reg, use_unifac=False)`** with a 150–650 K bracket.
* A root is always bracketed, because `psat_antoine` returns 0 when T + C ≤ 1.
* Flag `antoine_nominal_window_exceeded` when the flash is outside 0–100 °C (the anchored Antoine Tlo/Thi are the constant 273.15/373.15 K).
* **Revised after BASE Task 10 (U1 decided: LFL overlay adopted, mesitylene set rejected).**
  * The API's designated flash uses the overlay: `RouteCtx(reg, prov, lfl_override=cells_t2_flash.lfl_overlay_dict())`, the same call that produced `pred_lfl_overlay`.
  * The old route value is returned as `variants['v4_designated']`, from a plain `RouteCtx(reg, prov)`.
  * `scoreboard.csv` is not edited. Its `designated` row stays the v4 routing, and the API README says the API follows U1.
  * The lab margin stays the v4 yaml value of 3.70 °C, because the overlay leaves it unchanged.
  * Every flash result with a finite `lo`/`hi` carries the warning `lab_margin_set_by_one_system`. Task 10 found that this margin is the group score of the single mesitylene system, and that it falls to 0.31 °C when that system is left out.
  * A blend containing mesitylene also gets the warning `antoine_suspect:mesitylene`.

**D3b. What BASE Tasks 11–16 change in this plan.**
* **No change:**
  * SubReps stays as it is (Task 12: the C9 swap made flash worse).
  * No learned residual correction ships (Task 13: no transfer except a cold offset that K1 already gives).
  * The indicative interval stays labelled not held-out (Task 15: there are no fresh fuels).
* **Information only:** Task 16a showed that vapour pressure does not explain the flash error.
* **K1:** stays a variant with no margin. The −20 °C margin did not transfer between sources: the AV-33-22 margin covers 75 % of AJFTD fuels. The candidate margins in `delta_cold_viscosity/proposed_margins.yaml` are not read by the API.

**D4. Viscosity and density come from the frozen provider grid.**
* Never call `ensure()`; pin the npz SHA-256.
* Designated output: kinematic cSt from `refutas_nu_safe`, and dynamic = ν·ρ/1000.
* The K1 cold offset appears only under `variants`, and only at its fitted temperatures.

**D5. New properties.**

| Property | Designated model | Variant / label |
|---|---|---|
| Freezing | ideal SvL in mole fraction; label "no skill shown, bias +8 °C", uncertified | `f1x_recalibrated` = a + b·T_first, computed with the complete Coutinho n-alkane table (`F1x` in `S5/s5_predictions.csv`), a and b fitted on all 49 S5 headline groups. BASE Task 14 found that it beats the constant pooled and on AV-33-22, at MAE 3.07 °C. Quote only that out-of-sample number; never quote an in-sample one. Jet only; needs n-alkanes. Uncertified. |
| Net heat | `ga_comp_H`; yaml margin `heat_of_combustion_net_MJ_kg.crc647.group_additivity_composition_H` | jet only |
| Smoke point | YSI, with a 0.90 YSI-coverage gate | status `no_evaluation`; never in domain |
| Aromatics wt % and vol % | | `by_construction` |
| Distillation | cut | |

**D6. Uncertainty.**
* `lo`/`hi` stay NaN and `in_domain` False whenever there is no yaml margin.
* A separate `indicative_interval` holds the 90 % abs-error quantile of the frozen T2 residuals per property and T bin, labelled `empirical_T2_descriptive_not_held_out`.
* It is returned only for `fuel_type='jet'` with unmapped mass ≤ 0.10. Design never reads it.

**D7. Inverse design on an arbitrary pool: `RouteLibrary(Library)` accepts SMILES and InChIKeys.**
* Constants come from the verified table first, then from route sources, with a status per constant.
* A model that touches a non-verified constant becomes `uncertified` (via `margin_missing_reason`).
* A pool member without its own data for a constrained property is refused by name. No representatives inside a pool.
* Freezing in the adapter is mole-fraction SvL, the same as the forward model.
* `n_starts` defaults to 8 when the pool has more than 20 components.

**D8. Definition of "inverse works".** It finds blends that meet every limit at the point estimate, re-predicts them identically through `predict()`, and names what blocks certification.
* Freezing stays uncertified: 8 systems support at most an 88.9 % level, and lowering alpha after the failure would be post hoc.
* An optional declared 4 K freezing guard band is labelled "indicative guard band, not certified".

**D9.** See the execution order at the end.

**Spot-checks verified:**
* `pure_grid.npz` holds 1042 compounds; 114 have a full finite grid; the grid starts at T = 230 K; the file is gitignored.
* `fuelblend.eval.cold_offset` has `fit_k1(T, resid)` and `apply_k1(k1, T)`.

## Global Constraints

- Do not modify any of:
  - `results/v4/**`, `configs/**`, `splits/**`;
  - committed `results/v5/<task>/` files;
  - `results/v5/BASE/scoreboard.csv`;
  - `src/fuelblend/design/engine_v4.py`.
- Never call `PurePropertyProvider.ensure()` from the API or the design adapter. `data/cache/pure_grid.npz` is read-only and its SHA-256 is pinned.
- Keep `predict_blend._predict_viscosity` and `_predict_density` callable with their current signatures, because `mixture/common._fetch_via_api` imports them.
- No downloads and no AJFTD raw rows in any output.
- New data enters only in Tasks 10–11, and only by transcription from the two PDFs the user already put in `~/Downloads`.
  - Copy each PDF to `extras/staging/` (gitignored) first; never commit a PDF or quoted text.
  - Commit only the transcribed numbers, with the DOI, table number and page in every row.
- Commit the delta README with its declared readings alone, before any delta result is generated.
- A property without a valid yaml margin never gets a finite `lo`/`hi` or `in_domain=True`.
- Commit on branch `v5` with no attribution lines. Push only when the user asks.

## Review Focus

1. **A 60-component blend where no component has own data.** Expected: values with `substituted_mass ≈ 1`, `in_domain` False, no exception. Test: `test_all_substituted_blend` (Task 2).
2. **Fractions not summing to 1, duplicate components, or zero weights.** Expected: normalise, merge duplicates by InChIKey, drop zero weights, never return NaN. Test: `test_duplicates_and_zero_weights` (Task 2).
3. **Temperature outside 230–470 K.** Viscosity and density return `None` with status `T_out_of_grid` (today the grid clips silently). Flash, heat and aromatics are unaffected. Test: `test_T_outside_grid` (Task 2).
4. **A pool given as SMILES for compounds already in the verified table.** The verified constants win, and the design is identical to the name-based pool. Test: `test_smiles_pool_equals_name_pool` (Task 6).
5. **A non-hydrocarbon (ethanol) in a blend.** It is unmapped for physics, with a warning, and lab margins are off. Above 10 % mass the status is `low_coverage`. Test: `test_non_hc_component` (Task 2).

---

### Task 0: Declare the delta (README first, grid pinned)

**Files:**
- Create: `results/v5/BASE/delta_api_route/README.md`
- Create: `tests/test_delta_api_route.py`

**Interfaces:** Produces `GRID_SHA256` (in the README and the test), used by Task 1.

- [ ] **Step 1:** Compute the SHA-256 of `data/cache/pure_grid.npz`. Write the README with:
  - **Purpose.**
  - **Baseline:** `api_as_shipped` flash bias −21.66 and MAE 21.92 on 70 groups (`cells_t2_flash.csv`).
  - **Declared readings**, written before any result exists:
    - R-API-1: on all 79 T2 rows, the API's designated flash equals the route's `pred_lfl_overlay`, and `variants['v4_designated']` equals `pred_designated`, both to 1e-6 °C (columns of `rows_t2_flash.csv`).
    - R-API-2: API kinematic viscosity and density equal the route's values on all T2 rows at the 9 BASE temperatures, to 1e-9 relative.
    - R-API-3: flash, heat and aromatics are bit-identical across T_K ∈ {233.15, 288.15, 313.15}.
    - R-API-4: both D4 v5pool recipes re-predict through `predict()` within 0.01 °C for flash, 1e-6 relative for viscosity, and exactly for heat. Density and freezing differences are reported with their cause.
  - **Grid SHA-256 and disclosure:** the grid is an unpinned, gitignored binary built 2026-11-11. Its values differ from today's `_fetch_via_api` by up to 3.25 in ln η at 230 K. This is disclosed, not fixed here.
  - **Out of scope:** grid rebuild; U1 promotions (K1, F1, LFL overlay); distillation.
- [ ] **Step 2:** Write `test_grid_pinned`: the npz SHA-256 equals `GRID_SHA256`; skip with a reason if the file is absent.
- [ ] **Step 3:** Run `.venv\Scripts\python.exe -m pytest tests/test_delta_api_route.py -v`. Expected: PASS.
- [ ] **Step 4:** Commit the README and test alone: `v5 delta_api_route: declared readings and pinned grid (before results)`.

### Task 1: Route context and per-property component mapping

**Files:**
- Create: `src/fuelblend/api/predict_route.py`
- Test: `tests/test_predict_route.py`

**Interfaces (all produced here):**
- `class ReadOnlyProvider(PurePropertyProvider)`: `ensure()` is a no-op; the constructor raises `FileNotFoundError` if the npz is missing.
- `get_ctx() -> RouteCtx` (cached). It builds:
  - `CompoundRegistry(data/unified/compounds.csv)`;
  - `ReadOnlyProvider`;
  - `SubReps(reg, prov, vp_lookup=lambda s: _lookup(s, None))`, exposed as the attribute `reps`;
  - `RouteCtx(reg, prov)`.
- `@dataclass Component(ik: str, smiles: str, w: float, reg_idx: int | None, subkind: str | None, n_carbon: int, is_hc: bool)`.
- `resolve(composition, mole_fractions=None) -> list[Component]`:
  - accepts InChIKey or SMILES;
  - an InChIKey is looked up in the registry, and an InChIKey not in the registry raises `ValueError`;
  - a SMILES is canonicalised, then converted to an InChIKey with RDKit;
  - duplicates are merged, zero weights dropped, and weights normalised.
- `@dataclass PropMap(idx: np.ndarray, w: np.ndarray, own_mass: float, substituted_mass: float, unmapped_mass: float, substitutions: list[tuple[str, str, int]])`. Each substitution is (input InChIKey, representative InChIKey, carbon gap).
- `map_for(comps, has_own: Callable[[int], bool], reps) -> PropMap`. For each component:
  - if `reg_idx` is not None and `has_own(reg_idx)`, use the component itself;
  - otherwise use `reps.get(subkind, n_carbon, max_gap=99 if subkind in NEAREST_ANY else 3)` as the representative;
  - otherwise mark it unmapped. Non-hydrocarbons are always unmapped.
  - Representatives are aggregated as in `map_bins`.
- Own-data predicates:
  - `HAS_GRID = lambda i: prov.valid[i]`
  - `HAS_FLASH = lambda i: _lookup(reg.smiles[i], None) is not None`
  - `HAS_SLE`: the tmap entry is finite
  - `HAS_YSI`: the canonical SMILES is in the YSI table

- [ ] **Step 1: Write the failing tests:**
  - `test_resolve_smiles_and_inchikey_agree`: decane given as SMILES and as an InChIKey gives the same `Component`.
  - `test_resolve_merges_duplicates`.
  - `test_resolve_bad_smiles_raises`.
  - `test_map_for_own_data_no_substitution`: n-decane + n-dodecane + toluene with `HAS_GRID` gives `substituted_mass == 0` and `unmapped_mass == 0`.
  - `test_map_for_substitutes_missing`: a hydrocarbon without a grid gets one substitution with `gap <= 3`.
  - `test_map_for_route_equivalence`: for 5 T2 fuels from `results/v5/BASE/t2_fuels.csv`, the `map_bins` result fed through `resolve` and `map_for(HAS_FLASH)` gives the same `idx` and `w`.
  - `test_provider_never_writes`: the npz mtime is unchanged after `get_ctx()` and a lookup.
- [ ] **Step 2:** Run `pytest tests/test_predict_route.py -v`. Expected: FAIL (module missing).
- [ ] **Step 3:** Implement. Subkind comes from `reconstruction_r3.subkind(smiles)`; the carbon count from RDKit.
- [ ] **Step 4:** Run. Expected: PASS.
- [ ] **Step 5:** Commit: `v5 delta_api_route: route context and per-property component mapping`.

### Task 2: `predict_fuel` for viscosity, density and flash; rewire `predict()`

**Files:**
- Modify: `src/fuelblend/api/predict_route.py`.
- Modify: `src/fuelblend/api/predict_blend.py`:
  - `predict()` (lines 907–995) delegates to `predict_fuel` and keeps the old keys `viscosity_dyn_mPas`, `density_kg_m3` and `flash_point_C`;
  - delete `_flash_lechatelier` (lines 483–541) and the `gbm_mix` branch of `_predict_flash` (lines 708–752);
  - `_predict_flash` calls `predict_route`;
  - update the module docstring.
- Modify only these existing tests, which pin the old behaviour:
  - `tests/test_predict_blend.py::test_flash_antoine_out_of_range_warns`, and any flash value pinned to the `compounds.csv` Antoine constants;
  - the cold `gn_mass`/`gn_fallback` asserts in `tests/test_v4_api_known_composition.py::test_binary_is_continuous_with_pure`.
  - Every assert that a missing margin gives NaN `lo`/`hi` and `in_domain` False stays.
- Test: `tests/test_predict_route.py`.

**Interfaces (produced here):**
- `predict_fuel(composition, T_K: float = 298.15, fuel_type: str | None = None, mole_fractions=None, composition_kind: str = "known") -> dict`.
  - Top-level keys: `T_K`, `warnings`, `components` (the resolved list).
  - Each property maps to `{value, unit, method, status, coverage: {own_mass, substituted_mass, unmapped_mass, substitutions}, lo, hi, margin, in_domain, indicative_interval: None, variants: {}, warnings}`.
  - Properties: `viscosity_kin_cSt`, `viscosity_dyn_mPas`, `density_kg_m3`, `flash_point_C`.
  - `status` is one of `ok | low_coverage | T_out_of_grid | non_finite | no_evaluation`.
- Margins: reuse `_attach_uncertainty` unchanged.
  - The lab path applies only when `substituted_mass == unmapped_mass == 0`.
  - Otherwise `lo`/`hi` are NaN, `in_domain` is False, and the reason is `substituted_or_unmapped_mass`.

- [ ] **Step 1: Write the failing tests:**
  - `test_three_component_physics`: decane 0.5 / dodecane 0.3 / toluene 0.2 at 313.15 K gives flash method `ideal_le_chatelier_anchored`, status `ok` and `substituted_mass == 0`; density and viscosity equal a direct `refutas_nu_safe` call to 1e-9.
  - `test_flash_overlay_designated`: on that same blend, the flash value equals `flash_point_unifac_detail` with the overlay registry, and `variants['v4_designated']` equals it with the plain registry.
  - `test_mesitylene_warning`: a blend with 10 % mesitylene carries `antoine_suspect:mesitylene`, and every finite-margin flash carries `lab_margin_set_by_one_system`.
  - `test_flash_independent_of_T`: the same blend at 233.15, 288.15 and 313.15 K gives the identical flash (`==`).
  - `test_sixty_component_blend_runs_physics`: 60 jet-range hydrocarbons from the registry (seed 0; n-C8–C16, iso, cyclo and aromatics) give a physics method on every property (never `none`), status `ok`, finite values, and run in under 20 s.
  - `test_all_substituted_blend`, `test_duplicates_and_zero_weights`, `test_T_outside_grid`, `test_non_hc_component` (see Review Focus).
  - `test_no_gbm_mix_anywhere`: the string `"gbm_mix"` does not occur in the source of `predict_blend.py`.
  - `test_predict_wrapper_keys`: `predict()` still returns the three old keys plus `uncertified`.
- [ ] **Step 2:** Run. Expected: FAIL.
- [ ] **Step 3:** Implement:
  - viscosity and density: `map_for(HAS_GRID)`, then `refutas_nu_safe(prov, idx, w, T_K)`; T outside `[prov.T[0], prov.T[-1]]` gives `T_out_of_grid`;
  - flash: `map_for(HAS_FLASH)`, mole fractions from RDKit molecular weights, then `flash_point_unifac_detail(smiles, x, registry=ctx.flash_reg, use_unifac=False)`;
  - the `antoine_nominal_window_exceeded` warning;
  - then the rewiring in `predict_blend.py`.
- [ ] **Step 4:** Run `pytest tests/test_predict_route.py tests/test_predict_blend.py tests/test_v4_api_known_composition.py -v`. Expected: PASS.
- [ ] **Step 5:** Commit: `v5 delta_api_route: predict() on the scoring route; gbm_mix and compounds.csv flash removed`.

### Task 3: Freezing, heat, smoke, aromatics; K1 and F1 as labelled variants

**Files:**
- Create: `src/fuelblend/eval/freezing_f1.py`.
- Modify: `src/fuelblend/api/predict_route.py`.
- Modify: `predict_blend.py`, removing `freezing_point_C` and `smoke_point_mm` from `UNCERTIFIED_PROPERTIES`.
- Test: `tests/test_predict_route.py`.

**Interfaces:**
- Consumes:
  - `composition_route.ideal_sle_C` and `build_tmap`;
  - `s5_freezing_nalkane.t_first` and `audited_nalkane_constants`, imported via `importlib` as `cells_t2_bulk` does;
  - `mixture_heat_of_combustion_v2`;
  - the `PureData.comp_info` logic of `scripts/v5/base/cells_t2_bulk.py`, copied into `predict_route` as `_comp_info(ctx, i)` (the script stays untouched);
  - `sooting.ysi_mix` and `smoke_point_from_ysi`;
  - `cold_offset.fit_k1` and `apply_k1`.
- Produces:
  - `freezing_f1.fit_f1x() -> tuple[float, float]`: least squares of `freezing_C` on the complete-table first-species temperature, using the same column the S5 `F1x` model reads (check this in `s5_freezing_nalkane.py` around line 340). It fits over all S5 headline groups, one row per group, from `results/v5/S5/s5_predictions.csv`.
  - The constant `F1X_LOGO_MAE`: read it from `results/v5/BASE/delta_freezing/delta.csv` (pooled F1x, n = 49, about 3.07), never type it in. It is the only freezing score that may be quoted.
  - New `predict_fuel` keys:
    - `freezing_point_C`: ideal SvL is designated; `variants['f1x_recalibrated']` appears only when `fuel_type == 'jet'` and the n-alkane mass is above 0.
    - `heat_net_MJ_kg`: `ga_comp_H` with the yaml margin key above; jet only.
    - `smoke_point_mm`: YSI. Status `no_evaluation` when YSI own-coverage is below 0.90; never in domain.
    - `aromatics_wt_pct` and `aromatics_vol_pct`: `method='by_construction'`; the volume basis uses the grid ρ at 288.15 K.
  - `variants['k1_cold_offset_loso']` on `viscosity_kin_cSt`, only when T_K is within 0.01 K of a fitted BASE temperature. K1 is fitted once with `fit_k1` on all designated viscosity rows of `rows_t2_bulk.csv` and labelled `in_sample_T2`.

- [ ] **Step 1: Write the failing tests:**
  - `test_freezing_svl_matches_route`: on 5 T2 fuels, API freezing equals the route's `ideal_sle_C` to 1e-9.
  - `test_f1x_coefficients_reproduce`: `fit_f1x()` is finite, and a LOGO refit reproduces the `F1x` column of `s5_predictions.csv` to 1e-6 for 3 groups.
  - `test_heat_matches_bulk_rows`: on 5 T2 fuels, API heat equals the `ga_comp_H` prediction in `rows_t2_bulk.csv` to 1e-9.
  - `test_smoke_coverage_gate`: a blend that is 50 % compounds without YSI gives `no_evaluation`.
  - `test_aromatics_by_construction`: decane/toluene 0.8/0.2 gives `aromatics_wt_pct == 20.0`.
  - `test_k1_only_at_fitted_T`: the variant is present at 233.15 K and absent at 240.0 K.
  - `test_uncertified_props_have_nan_interval`: freezing and smoke have NaN `lo`/`hi` and `in_domain` False.
- [ ] **Step 2:** Run. Expected: FAIL.
- [ ] **Step 3:** Implement.
- [ ] **Step 4:** Run. Expected: PASS.
- [ ] **Step 5:** Commit: `v5 delta_api_route: freezing, heat, smoke, aromatics in predict(); K1/F1 as labelled variants`.

### Task 4: Indicative interval

**Files:** Modify `src/fuelblend/api/predict_route.py`; test in `tests/test_predict_route.py`.

**Interfaces:** Produces `_indicative_interval(prop: str, T_K: float, value: float) -> dict | None`, returning `{lo, hi, q90_abs, n_groups, kind: 'empirical_T2_descriptive_not_held_out'}`.
- Built once from the frozen residuals:
  - flash: `rows_t2_flash.csv`, error = `pred_designated - measured`;
  - bulk properties: `rows_t2_bulk.csv`, designated variant, `status == 'ok'`, `~low_coverage`, one row per group, nearest BASE T bin.
- Viscosity is in ln space; the other properties are absolute.
- Returned only when `fuel_type == 'jet'` and `unmapped_mass <= 0.10`.
- Never written into `lo`/`hi`.

- [ ] **Step 1: Write the failing tests:**
  - `test_indicative_flash_q90`: `q90_abs` equals `np.quantile(abs(err), 0.9)` recomputed from the CSV inside the test.
  - `test_indicative_absent_without_jet`.
  - `test_lo_hi_untouched_by_indicative`.
- [ ] **Step 2:** Run. Expected: FAIL.
- [ ] **Step 3:** Implement.
- [ ] **Step 4:** Run. Expected: PASS.
- [ ] **Step 5:** Commit: `v5 delta_api_route: indicative T2 interval, labelled, separate from margins`.

### Task 5: Delta result — API vs route on the 79 T2 fuels, plus a 60-component demo

**Files:**
- Create: `scripts/v5/base/delta_api_route.py`, writing to `results/v5/BASE/delta_api_route/`:
  - `delta.csv`
  - `api_vs_route_rows.csv`
  - `forward_demo_60.csv`
- Modify: `results/v5/BASE/delta_api_route/README.md`. Append a "Result" section only; the declared text stays unchanged.
- Test: `tests/test_delta_api_route.py`.

**Interfaces:**
- Uses the same fuel loop as `cells_t2_flash.main` (`kind_bins`, then `map_bins`), then calls `predict_fuel` on `[(reg.ik[i], w)]`.
- `delta.csv` columns: `reading, n_rows, max_abs_diff, pass, before_bias, after_bias, after_mae`.
- `forward_demo_60.csv`: the Task 2 60-component blend at 233.15, 253.15, 288.15 and 313.15 K, with every property, method, status and coverage.

- [ ] **Step 1: Write the failing test** `test_delta_readings_pass`. It reads `delta.csv` and asserts:
  - R-API-1, R-API-2 and R-API-3 all have `pass == True`;
  - the pooled flash `after_bias` equals the bias of the frozen `lfl_overlay` T2 cell to 1e-6, with the value read from `cells_t2_flash.csv`, not typed in.
- [ ] **Step 2:** Run. Expected: FAIL (no CSV yet).
- [ ] **Step 3:** Run `.venv\Scripts\python.exe scripts/v5/base/delta_api_route.py`. It prints pass or fail for each reading.
- [ ] **Step 4:** Run the test. Expected: PASS. If a reading fails, record it and its cause in the README's Result section; never edit the declared text.
- [ ] **Step 5:** Commit: `v5 delta_api_route: API equals scored route on T2 (R-API-1..3)`.

### Task 6: `RouteLibrary` — design on a SMILES/InChIKey pool

**Files:**
- Create: `src/fuelblend/design/library_route.py`
- Test: `tests/test_v5_design_forward.py`

**Interfaces:**
- Consumes: `models_v4.Library`, `make_fns` and `build_models`; `predict_route.get_ctx`, `map_for` and `_comp_info`.
- Produces `RouteLibrary(pool_ids: list[str], pure_props=data/processed/v5/pure_props_v5.csv)`, a subclass of `Library`:
  - A pool id is a table name, a SMILES or an InChIKey. A SMILES or InChIKey whose InChIKey matches a table compound uses that table name.
  - Compounds not in the table get constants with a status each:

    | Constant | Source | Status |
    |---|---|---|
    | MW | RDKit | `verified` |
    | rho15 | grid ρ at 288.15 K | `route_grid` |
    | Tm, dHfus | `build_tmap` | `coutinho_table` |
    | flash | anchored registry | `route_anchored` |
    | RON, MON, DCN, Psat37_8 | none | absent |

  - `v(name, prop)` returns table-verified values as before. It also returns route-sourced values, recording each use in `self.nonverified_used`.
- Produces `build_route_models(pool_ids) -> dict`:
  - starts from `build_models(RouteLibrary(pool_ids))`;
  - swaps freezing to `freezing_ideal_mole`;
  - gives any `PropertyModel` whose call touched `nonverified_used` the `margin_missing_reason` "non-verified constant(s): …", which the unchanged engine reports as `uncertified`;
  - raises `ModelUnavailable` naming the member when a pool member has no own data for a constrained property.
- Produces `design_forward_check(result: DesignResult, T_C_visc: float = -20.0) -> pd.DataFrame`. It re-predicts `result.weights` through `predict_fuel` and returns the columns `constraint, design_pred, forward_pred, abs_diff, tol, pass, cause`.

- [ ] **Step 1: Write the failing tests:**
  - `test_name_pool_unchanged`: the D4 v5pool given by name gives the weights in `results/v5/D4/designs_v4.csv` to 1e-9.
  - `test_smiles_pool_equals_name_pool`.
  - `test_nontable_compound_uncertifies`: add the first n-alkane absent from `pure_props_v5.csv` to the pool. The design runs, and every constraint that used that compound is `uncertified` with a reason.
  - `test_member_without_own_data_refused`.
  - `test_forward_check_d4`: the R-API-4 tolerances hold for both v5pool recipes, and the density and freezing rows have a non-empty `cause`.
- [ ] **Step 2:** Run. Expected: FAIL.
- [ ] **Step 3:** Implement.
- [ ] **Step 4:** Run the new tests, then `pytest tests/test_v4_design.py -v`. Expected: PASS for both, with `test_v4_design.py` unchanged.
- [ ] **Step 5:** Commit: `v5 delta_api_route: design on SMILES/InChIKey pools with honest constant statuses; forward check`.

### Task 7: User-pool design runs (presentation evidence)

**Files:**
- Create: `scripts/v5/run_design_userpool_v5.py`, writing to `results/v5/BASE/delta_api_route/designs/`.
- Modify: the README's "Result" section.

**Interfaces:**
- Command-line options:
  - `--spec astm:jet_a|bis:jet_a1`;
  - `--pool` (comma-separated ids, or a .txt file);
  - `--n-starts N` (default 30, or 8 above 20 components; the value used is written to the outputs);
  - `--seed 42`;
  - `--freezing-guard-K 0|4`;
  - `--out DIR`.
- Outputs: `designs.csv`, `constraints.csv`, `forward_check.csv`, `runtime.txt`.

- [ ] **Step 1:** Run the D4 v5pool for both specs. The recipes match D4; the script prints the difference.
- [ ] **Step 2:** Write the 12-component pool rule into the README before any run. The pool is the v5pool plus, for each of these (subkind, nC) bins, the first registry compound (lowest InChIKey) with own grid, flash and SLE data:

  | Subkind | nC |
  |---|---|
  | n_paraffin | 13 |
  | n_paraffin | 14 |
  | cyclo_mono | 9 |
  | cyclo_mono | 10 |
  | indan_tetralin | 10 |
  | aromatic_di | 11 |

  Then run both specs. A bin with no qualifying compound is recorded and skipped.
- [ ] **Step 3:** Run once with `--freezing-guard-K 4`, labelled "indicative guard band, not certified".
- [ ] **Step 4:** Add a README table. Per design:
  - category: A (fails at the point estimate), B (meets the point estimate but cannot be certified), or certified;
  - the named blocking gaps;
  - whether the forward check passes;
  - runtime;
  - n_starts sensitivity (8 vs 30 starts on the 12-component pool).
- [ ] **Step 5:** Commit: `v5 delta_api_route: user-pool designs, guard-band run, forward checks`.

### Task 8: Full suite, usage notes, handoff

**Files:**
- Modify: `results/v5/BASE/delta_api_route/README.md`. Add a "How to call it" section with short examples: `predict()` on a 60-component list, and `run_design_userpool_v5.py`.
- Modify: `AGENTS.md` §12, adding one paragraph that points to the delta.
- `docs/FULL_REPORT_v4.md` is not touched.

- [ ] **Step 1:** Run `.venv\Scripts\python.exe -m pytest -q`. Expected: the count after BASE Tasks 10–17 (483 passed at `0a85f20`) plus the new tests, with no new failures. List every changed test and its reason in the README.
- [ ] **Step 2:** Run `.venv\Scripts\python.exe scripts/v5/base/reproduce_base.py` (default mode). Expected:
  - the drift check passes;
  - `scoreboard.csv` is byte-identical to commit `5a62e20`;
  - `delta_*` folders are outside the manifest (commit `23d3914`), so `delta_api_route/` cannot fail the check.
- [ ] **Step 3:** Commit: `v5 delta_api_route: suite green, usage docs`.
- [ ] **Step 4: Whole-branch review.** Run one fresh reviewer (superpowers:requesting-code-review, most capable model) over `0a85f20..HEAD`. The BASE session skipped this review; it is part of this plan so the branch gets one before the presentation.
  - Fix what it confirms.
  - List what is left unfixed in the README with a reason.
  - Commit the fixes as `v5 delta_api_route: review fixes`.

## Part B: what the BASE session (Tasks 10–17) left open

The BASE session committed through `0a85f20`. It left these items open. Tasks 9–11 run in this plan; the rest are listed with what blocks them.

### Task 9: Bookkeeping that is waiting for the user (gated)

- [ ] **Step 1 (user go-ahead needed: Q1).** Run `.venv\Scripts\python.exe scripts/v4/reproduce_all.py --verify-only --update-manifest`, so that `SHA256SUMS_v4.txt` matches the tree again.
  - Before running it, list the files whose hashes change and show them to the user.
  - Commit only the manifest: `v4 manifest: re-hash after v5 BASE (no results/v4 content changed)`.
- [ ] **Step 2 (user go-ahead needed: Q3).** Commit both plan files on their own, before Task 0: this plan and `2026-11-12-v5-full-baseline.md`.
  - Commit message: `docs: v5 plans (full baseline; API route and design)`.
  - Leave `docs/data_audit_2026-11-11.md` alone; it belongs to the BASE session.

### Task 10: U2, more lab flash systems from Li et al. 2014 (n-alkane binary and ternary flash points)

Why: Task 10 showed that the lab flash margin of 3.70 °C rests on one system out of 31. More measured hydrocarbon mixture systems are the honest way to settle Q2. Choosing between 3.70 and 0.31 °C after seeing them is not.

**Files:**
- Copy the PDFs from `~/Downloads` into `extras/staging/`:
  - `Advances in Materials Science and Engineering - 2014 - Li - Prediction of the Flash Point of Binary and Ternary.pdf` (doi 10.1155/2014/239023);
  - `Measurement of Lower Flash Point for Ternary Mixture, n-Nonane+n-Decane+n-Dodecane System.pdf` (the same kind of data; check whether it is already among the lab rows).
- Create:
  - `data/raw/flash_lab_li2014/li2014_flash_mixtures.csv`, with columns `system_label, component_name, component_smiles, x_mole, measured_C, method, doi, table, page`;
  - `src/fuelblend/data/loaders_li2014.py::load_li2014() -> pd.DataFrame`, which outputs rows in the layout of `results/v4/F/flash_lab_loso.csv` (`system_id` = md5[:12] of the sorted InChIKey set, `comp_smiles` and `x_mole` as JSON);
  - `results/v5/BASE/delta_lab_flash_li2014/README.md`, then `delta.csv` and `lab_margin.csv`;
  - `scripts/v5/base/delta_lab_flash_li2014.py`;
  - `tests/test_li2014_loader.py` and `tests/test_delta_lab_flash_li2014.py`.

- [ ] **Step 1:** Open the PDF.
  - Record the licence (Hindawi journals are usually CC BY; confirm it on the first page) and the tables that hold measured mixture flash points.
  - Transcribe only measured values; skip predicted columns.
  - Write a fixture test that copies 3 rows exactly as printed.
- [ ] **Step 2:** Run the de-duplication check.
  - Compare against `results/v4/F/flash_lab_loso.csv`: same `system_id`, x within 0.005 and flash within 0.5 °C counts as the same point. These are already in the 31 systems (Leenhouts compiled literature), so they are counted and dropped, not added twice.
  - Do the same check for the nonane/decane/dodecane ternary.
  - Write the counts (new systems, new points, duplicates) in the README.
- [ ] **Step 3: Declare readings in the README and commit it alone, before any score.**
  - R-LI-1: on the new points, the overlay flash has |bias| ≤ 2 °C and an over-prediction rate (> 5 °C) ≤ 10 %.
  - R-LI-2: the lab one-sided margin (`one_sided_margin(direction="over", alpha=0.10)`) over the old plus new systems, computed with and without the mesitylene system, differs by less than 0.5 °C.
  - If R-LI-2 holds, that margin is proposed for `configs/uncertainty_v5.yaml`. Writing that file stays the user's decision (Q2).
  - The new points are labelled `dev` (first read), never `held_out`.
- [ ] **Step 4:** Implement the loader and the delta script. The script uses the same per-row prediction as `delta_flash_inputs.lab_predictions`, with both the overlay and the designated variant. Run the tests and the script.
- [ ] **Step 5:** Append a Result section to the README and commit: `v5 delta_lab_flash_li2014: n-alkane mixture flash points, margin sensitivity`.

### Task 11: U4, diesel surrogates from Mueller et al. 2016 (known-composition test of `predict()`)

Why: a surrogate is a recipe of pure compounds with measured properties. That makes it the cleanest real test of the any-size API on a known composition, and the diesel cell is `no_evaluation` today.

**Files:**
- Copy `~/Downloads/Mueller et al. 2016, doi.pdf` to `extras/staging/`. Check on page 1 that it is doi 10.1021/acs.energyfuels.5b02879; if not, stop and tell the user. This is an ACS paper, so never commit the PDF or any text from it.
- Append the recipes to `data/raw/surrogates/fuel_surrogates.csv` in its existing columns (`fuel_code = Mueller2016_<name>`, full reference with the DOI and table).
- Create `data/raw/surrogates/mueller2016_properties.csv`, with columns `fuel_code, property, value, unit, T_C, method, table, page`. It covers flash, density and cloud point, and the target diesel's measured values where the paper gives them.
- Create `results/v5/BASE/delta_diesel_surrogates/` (README first), `scripts/v5/base/delta_diesel_surrogates.py` and `tests/test_delta_diesel_surrogates.py`.

- [ ] **Step 1:** Transcribe.
  - Each recipe's mole or mass fractions must sum to 100 ± 0.5; the test asserts this.
  - Every component SMILES is checked with RDKit, and its InChIKey must exist in `data/unified/compounds.csv`. A component that is not there is listed, never guessed.
- [ ] **Step 2: Declare in the README before scoring, and commit it alone.**
  - The cell is descriptive and `n_too_small` (n = 4), so no skill claim is made.
  - The readings are the per-surrogate flash, density and cloud-point errors, and the coverage dict.
  - Cloud point is scored with the designated ideal SvL freezing estimate, labelled "SvL vs cloud point, different property, indicative only".
- [ ] **Step 3:** Run `predict_fuel(..., fuel_type='diesel', composition_kind='known')` on each recipe at the measurement temperature.
  - Write `rows.csv` with value, status, coverage and measured.
  - Write a second table with the surrogate's own predicted values against the target diesel's measured values.
- [ ] **Step 4:** Add the Result section and commit: `v5 delta_diesel_surrogates: Mueller 2016 surrogates through predict()`.

### Left open: blocked or after the 17th

| Item | What | Blocked by | When |
|---|---|---|---|
| BASE 16c | FuelLib 18 Jet-A/HEFA blends as a real-fuel blending-rule test | U7: FuelLib HEFA-blend measurement temperatures from the Vozka paper | after the user answers U7 |
| BASE 16d | Supervised reconstruction (D86, simdis, aromatics, density → light end and n-paraffin share), leave-source-out, re-score CRC 647 dev | about 2 days | after the 17th |
| BASE 16e | Check the two S4 flash correlation forms against their sources | the source documents are not on disk | when the user supplies them |
| BASE 16f | Refit the v4 structured ridge and bilinear viscosity models on CRC 647 dev and apply them to T2 (fills three `no_evaluation` rows) | about 1 day; the fit code must be lifted out of `build_viscosity_leaderboard_v4.py::main` | after the 17th |
| K1 pure-viscosity fix | Extrapolate each compound from its own measured range (51–60 % of fuel mass is extrapolated at −20/−40 °C) | goes together with the grid rebuild | after the 17th |
| Freezing D86 lead | Ridge on D86 T90, FBP and n-paraffin share (MAE 2.98 vs 6.49, exploratory) | needs a second source to be more than exploratory | after the 17th |
| LFL and light-end check | Task 16a pointed at LFL values or the light end as flash suspects | nothing; it is analysis | after the 17th |

## Execution order and cut list

- **Before Task 0:** Task 9, once the user answers Q1 and Q3. Neither blocks Task 0 if the answer is no.
- **13 Nov morning:** Tasks 0–2.
- **13 Nov afternoon:** Tasks 3–5.
- **14 Nov morning:** Tasks 6–7.
- **14 Nov afternoon:** Task 8, including the whole-branch review.
- **15 Nov morning:** Task 10 (U2), then Task 11 (U4), each with its README committed before it is scored.
- **15–17 Nov (already agreed):**
  - repo simplification on a branch, with a tag on the full state;
  - the methodology document;
  - the manuscript;
  - the slides.

If time runs short, cut in this order:
1. distillation (already cut);
2. Task 11 (U4 diesel surrogates);
3. Task 10 (U2 lab flash), except Step 2: the duplicate count is quick and tells the user whether the PDF adds anything;
4. smoke point (Task 3 smoke only);
5. the 12-component design runs (Task 7 Step 2; the v5pool runs stay);
6. the K1 and F1x variants;
7. the indicative interval (Task 4).

Tasks 0–2, 5 and 6 are the fix itself and are never cut.

## User decisions this plan does not take

- **U1 (decided 2026-11-12):**
  - The LFL overlay is adopted, so it is the API's designated flash (D3).
  - The mesitylene Antoine set is rejected.
  - K1 is a variant below 263 K, with no margin.
  - Whether F1x ever becomes the designated freezing model is still open; here it is a variant.
- **Q1:** may Task 9 Step 1 re-hash `SHA256SUMS_v4.txt`?
- **Q2:** the flash margin. Keep 3.70 °C, which is conservative but set by the suspect mesitylene system, or wait for Task 10's evidence. Until then the API uses 3.70 °C with the `lab_margin_set_by_one_system` warning, and no `configs/uncertainty_v5.yaml` is written.
- **Q3:** may Task 9 Step 2 commit the two plan files?
- **U7:** the FuelLib HEFA-blend measurement temperatures, which unblock BASE 16c.
- **Grid disclosure:** whether to mention the `pure_grid.npz` finding (an unpinned binary behind the frozen BASE) to the faculty. Recommended: yes. The delta README discloses it either way.
- **Grid rebuild:** whether to rebuild the grid as a separate delta after the 17th. A rebuild changes every BASE viscosity and density number.
