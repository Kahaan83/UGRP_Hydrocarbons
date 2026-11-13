# delta_smoke_margin: smoke-point margin 7.8 mm

## Decision

User decision 2026-11-13: "Yes, add it", answer to "Add the 7.8 mm smoke margin to uncertainty_v5.yaml?". The user was told beforehand that the
margin comes from binaries and surrogates and probably certifies nothing (designs predict about 27.6 mm against a limit of 25, so 27.6 - 7.8 < 25).

Source row: `results/v5/BASE/delta_newdata_freeze_smoke/summary.csv`, scope `smoke_mixtures`, model `bulk_designated`,
`margin_alpha10` = 7.8009965628100595 (direction over, 12 systems, 56 scored points). This exact float is used.

## Scope caveats

* The 12 systems are binaries and surrogates, not real jet fuels.
* 37 of the 56 scored points lie outside the regression's training domain (aromatics 8-23 vol %).
* The margin (7.8 mm) is wider than the 19-25 mm range the regression is trusted on.
* The user chose it after seeing these numbers. It is a bound from out-of-domain mixtures, not a coverage statement for jet fuels.

## What changes

* `configs/uncertainty_v5.yaml` (regenerated) gains `smoke_point_mm.lab_mixtures.bulk_regression`. The flash entry stays 2.206932808388615.
* API `predict_route._smoke`: margin applied only when the aromatics and result are inside the regression domain.
* Design side: the smoke model reads the new key; the existing OutOfDomainWarning gate keeps out-of-domain designs uncertified.

## Declared readings

* **R-SM-1**: `predict_fuel` smoke on an in-domain blend (aromatics 8-23 vol %, result 19-25 mm) returns margin 7.8009965628100595,
  lo = value - margin, hi = value + margin, `safe_bound` = lo, in_domain True.
* **R-SM-2**: out-of-domain blends (aromatics outside 8-23 vol % or result outside 19-25 mm) keep NaN lo/hi and in_domain False.
* **R-SM-3**: all non-smoke outputs of `predict_fuel` are unchanged (same comparison method as delta_flash_margin_q2 R-Q2-2).
* **R-SM-4**: the four design runs (seed 42, as in `../delta_flash_margin_q2/`) are re-run; per design the smoke constraint status and the overall
  category are reported. Expectation: smoke is still not certified.

## Result

Full pytest: 581 passed, 0 failed (was 575 + the delta test; +2 new tests, existing smoke tests updated).

* **R-SM-1: pass.** Blend decane 0.5 / dodecane 0.3 / toluene 0.2 at 288.15 K: smoke 21.89 mm, aromatics 17.5 vol %, margin 7.8009965628100595 mm,
  lo = value - margin = 14.08, hi = 29.69, `safe_bound` = lo, `in_domain` True, warning `smoke_margin_from_binaries_and_surrogates (user decision)`
  (test `test_smoke_margin_in_domain_user_decision`). The API applies it only when no `out_of_domain` warning was raised and nothing was substituted or unmapped.
* **R-SM-2: pass.** Decane/toluene 50/50 (45.7 vol % aromatics) and decane/dodecane (0 vol %, 27.6 mm) keep NaN lo/hi, `in_domain` False, `uncertified` True
  (test `test_smoke_margin_not_applied_out_of_domain`).
* **R-SM-3: pass.** `check_r_sm3.py` / `check_r_sm3.txt`: the 60-component demo (seed 0) and a 5-compound blend at 233.15, 288.15, 313.15 K, 72 records;
  only the 6 `smoke_point_mm` records differ; 0 non-smoke records differ. "Before" is the same code with the smoke key removed from the margins dict
  (equivalent to the previous behaviour), not a separate checkout.
* **R-SM-4: designs re-run, and the expectation "nothing else changes" did not hold.** Smoke is not certified in any design (status `fail`, not `uncertified`),
  but all 8 designs moved from category B to category A (`A_fails_a_limit_on_design_or_API_point_estimate`), 0 certified overall.

| run | design | smoke (pred 27.56, margin 7.801) | other non-pass constraints | category (was) |
|---|---|---|---|---|
| v5pool | Jet A | fail | density, freezing | A (B) |
| v5pool | Jet A-1 | fail | density, freezing | A (B) |
| ext_n30 | Jet A | fail | freezing | A (B) |
| ext_n30 | Jet A-1 | fail | density, freezing | A (B) |
| ext_n8 | Jet A | fail | freezing | A (B) |
| ext_n8 | Jet A-1 | fail | density, freezing | A (B) |
| guard4_v5pool | Jet A | fail | density, freezing | A (B) |
| guard4_v5pool | Jet A-1 | fail | density, freezing | A (B) |

Why: 27.56 - 7.80 = 19.76 mm is below the 25 mm minimum, so smoke can never pass with the margin, and the engine ranks "limit violated after applying the margin"
(fail) above the out-of-domain block (uncertified). With smoke infeasible the optimiser returns a compromise: the recipes changed substantially (butylcyclohexane
0.61-0.67 before, 0.10-0.21 now; decalin 0 before, 0.15-0.31 now) and several point-estimate limits (density at 15 C, freezing) are now violated at the point
estimate, on both the design and the API side (forward check passes: design and API agree). Before the margin the optimiser could satisfy all limits at the point
estimate. The smoke result is a weak extrapolation (0 vol % aromatics), so the optimiser is now being driven by a margin that was never validated there.
`engine_v4.py` was not changed; making smoke `uncertified` out of domain instead of `fail` would need an engine rule change (not done, out of scope).

### Domain fix

What went wrong: the first rerun (`designs/`, kept as the record) gave the optimiser the smoke margin through the engine's yaml key. The engine's margin does
not depend on the recipe, so the margin applied at 0 vol % aromatics (outside the regression domain), smoke (27.56 mm) could never reach 25 + 7.8, and the
optimiser gave up density and freezing. All 8 designs fell from category B to A. This contradicted brief step 4 (smoke stays uncertified out of domain).

Fix (design side only; `engine_v4.py` unchanged): `library_route` no longer hands the smoke key to the optimiser, so the optimiser and the engine see smoke as
before the margin (OutOfDomainWarning gives `uncertified`, margin None). After the engine, `annotate` / `evaluate_recipe` re-judge smoke constraints that are
IN the domain (no OutOfDomainWarning in the reason): margin 7.8009965628100595 applied, slack = (prediction - margin - limit) / scale, `fail` if negative.
The optimiser is therefore not margin-aware in domain; an in-domain recipe is only judged correctly after the fact. Test:
`tests/test_v5_design_forward.py::test_smoke_margin_design_side_domain_gated` (0 vol % recipe: uncertified, margin None; in-domain recipe: margin applied).

Rerun into `designs_domain_fix/` (seed 42, same arguments). Finished: v5pool and ext_n30. NOT finished at the time of writing: ext_n8 and guard4_v5pool
(still running; their folders may be incomplete or absent and are not part of this result).

| run | design | smoke status | smoke margin | category | vs 9ee34cc |
|---|---|---|---|---|---|
| v5pool | Jet A | uncertified | none | B | recipe identical (max weight change 0.0000) |
| v5pool | Jet A-1 | uncertified | none | B | recipe identical |
| ext_n30 | Jet A | uncertified | none | B | recipe identical |
| ext_n30 | Jet A-1 | uncertified | none | B | recipe identical |

Tests run for this fix: `tests/test_v5_design_forward.py` and `tests/test_predict_route.py` only (the full suite was not rerun, by instruction; the earlier
full run of the first version gave 581 passed).

Update (reruns completed): ext_n8 and guard4_v5pool also finished: all four designs in each are category B (8 of 8 designs B, smoke uncertified, 0 certified).
Full pytest after the fix: 582 passed, 0 failed.
