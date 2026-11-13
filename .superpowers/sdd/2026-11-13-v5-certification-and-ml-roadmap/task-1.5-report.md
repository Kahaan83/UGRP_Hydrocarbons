# Task 1.5 report: recipe-dependent margins in the optimiser

## Baseline (Rule 0)
Full `pytest -q` at HEAD e8b7819 (before edits): 1 failed, 643 passed, 8 skipped (19 min). The one failure, `tests/test_v4_sle.py::test_default_mode_unchanged` (NotADirectoryError), is environmental and pre-existing. Output: `task-1.5-pytest-base.txt`. Final run: `task-1.5-pytest-final.txt`.

## What changed
- `engine_v4.py`: `PropertyModel` gains `margin_fn` (plus `margin_ln`, `margin_direction`, `margin_source`). A `None` margin from `margin_fn` makes that constraint side infeasible (`BAD_SLACK`) and the status uncertified, never margin 0. `_Sides` reads the margin per weight vector inside the existing `_k/_v` memo (d1b2b60 pattern). `_Sides` is alt-aware when `spec["use_alt_limits"]`: smoke slack = max(primary, min(alt smoke slack, naphthalenes slack)); if naphthalenes is unavailable the alt branch is infeasible. `evaluate_blend` computes the same quantity (`_alt_combo_slack`), so optimiser and evaluator agree. `robustness_table` accepts `alt_limits` / `synthetic_blend`.
- `library_route.py`: `with_recipe_margins(models, margins)` attaches `margin_fn` to viscosity (ML R-MLV-6 margin when `visc_ml_status` allows, else the Refutas lab margin of the T bin) and smoke (7.8 mm lab key, or the v5 real-fuel key with the class-mix gate). Models with `margin_missing_reason` are skipped and stay uncertified. Flag `recipe_margins=False` by default on `design_routed` and `evaluate_recipe`. `_rejudge_smoke` is alt-aware; `_consistency` compares the engine and rejudge on status, margin and slack to 1e-9.
- Rule 3: if `smoke_point_mm.ajftd_v5_in_domain.bulk_regression` is present in the yaml it replaces the lab key for in-domain, in-class-mix recipes; outside the class-mix range the margin is None, reason `class_mix_outside_real_fuel_range`. Key absent: 7.8 mm behaviour unchanged. The same logic is in `predict_route._smoke`. The yaml was not edited; tested with an in-memory config.
- `smoke_margin_v5.py`: `recipe_class_fractions` / `class_mix_ok_for_recipe` use the GCxGC basis (cyclo = mono+di+tri, same normalisation, decalins in cyclo); `_EPS` tolerance at range edges.
- `scripts/v5/run_design_userpool_v5.py`: `--recipe-margins` (default off).

## Tests (`tests/test_v5_margin_fn.py`, 25 pass)
`test_margin_fn_matches_rejudge` (8 stored recipes), `test_margin_fn_absent_unchanged` (baseline JSON captured before edits by `tests/data/gen_margin_fn_baseline.py`), `test_optimiser_sees_smoke_margin` (stub 3.0 mm, seed 42, n_starts=5, alt on), `test_none_margin_is_infeasible`, `test_class_mix_basis` (a real GCxGC fuel's own fractions pass the gate), in-memory v5-key tests, and predict_route tests.

## Concerns
1. Under the alt flag the reported smoke slack is now max(primary, alt combo) per Rule 4. This is intended and differs from the old alt=True number.
2. The v5 pool has no iso-alkane, so the real-fuel margin gate would always reject its recipes (iso-paraffin range lower bound 0.0068 > 0). The class-mix tests use isooctane because isododecane is unavailable (disputed rho15) in the routed models.
3. The engine reports FAIL for a recipe that satisfies only the alt limit when the spec is `secondary_sources_only`, as in 1.2; the legacy rejudge leaves it uncertified, so the status part of the consistency check is one-way for smoke.
4. `test_v4_sle` failure is environmental (also in baseline).
5. Full `design_routed(recipe_margins=True)` runs were not executed (10-20 min each); tests use `evaluate_recipe` and stub models.
