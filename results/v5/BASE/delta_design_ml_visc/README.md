# delta_design_ml_visc

Declared 2026-11-13 before any design is re-run; this commit contains no result file. Frozen reference: branch v5 at commit `5d57504`
(delta_freeze_margin). The declared text below is never edited after results; results are appended under a "Result" heading.

## Why
Since commit `cc5f5ef` (user decision "Make it main") `predict_fuel` designates the ML pair-residual viscosity (Refutas x exp(sum x_i x_j g_ij))
inside its domain (known composition, own grid data, every class pair with >= 5 training systems, 233-393 K), with the R-MLV-6 ln margin of its
temperature bin (0.287 below 263 K, 0.608 for 263-283 K, 0.134 at or above 283 K; `results/v5/BASE/delta_ml_viscosity/readings.json`). The design
side (`library_route`) still predicted and certified viscosity with raw Refutas and the yaml Refutas margins, so the design and the API disagreed: the
forward-check row `kin_visc_m20C_cSt[tol]` failed in all 8 designs of `delta_freeze_margin` (tolerance 1e-6 relative).

## What changes (design side only; `engine_v4.py` unchanged)
* `library_route.build_route_models`: the viscosity model returns the API designated value: Refutas x exp(ML correction) when
  `visc_ml_status` (same rule as the API for a known composition) says in domain, else Refutas.
* `library_route.annotate` -> `_rejudge_viscosity`: after the engine, an in-domain viscosity constraint is re-judged with the ML margin of its T bin
  (two-sided ln-ratio; slack = (limit - pred x exp(m)) / scale for a maximum limit; `fail` if negative). The optimiser itself runs with the Refutas
  yaml margin of the bin (the engine's margin cannot depend on the recipe's domain), so the optimiser is not ML-margin aware; the final recipe is
  judged correctly after the fact. Out of domain the engine's Refutas judgement is left unchanged, with an `ml_not_designated` note.
* Forward check: the `[tol]` row now compares like with like (both designated values); the old report-only `[ml_pair_residual]` row becomes an
  informational `[refutas_raw]` row (the Refutas value the ML value replaced).

## Declared readings (descriptive; no gate)
* **R-DV-1**: the four design runs of `delta_freeze_margin/designs/` (`v5pool`, `ext_n30`, `ext_n8`, `guard4_v5pool`; seed 42,
  `scripts/v5/run_design_userpool_v5.py`, same arguments) re-run into `designs/`. Per design: viscosity prediction, model, margin, slack and status
  before (`delta_freeze_margin`) and after; recipe change (max weight difference); category; forward-check pass.
* **R-DV-2**: the forward-check `kin_visc_m20C_cSt[tol]` row passes in every design (expected, by construction).
* Expectation stated now: viscosity at -20 C is in the ML domain for these pools (n-alkanes, butylcyclohexane, decalin, trimethylbenzene), so the
  margin changes from 0.415 (Refutas, T < 263 K bin) to 0.287; viscosity stays `pass` (predictions near 3 cSt vs 8 cSt). Categories stay B and
  0 designs are certified, because smoke blocks every design (aromatics 0 vol %, out of the regression domain). Nothing here is a new certification
  claim beyond the per-constraint statuses.

## Result (appended 2026-11-13; declared text above unchanged)
Runs re-done into `designs/<run>/` (seed 42, same arguments). v5pool, ext_n8, guard4_v5pool ran before commit `d1b2b60` (memoised ML pair matrix,
same values); ext_n30 was stopped and re-run after it (speed only).

R-DV-1, viscosity constraint `kin_visc_m20C_cSt` (limit <= 8 cSt), before (`delta_freeze_margin/designs`) and after:

| run | design | before: model, pred (cSt), margin (ln), slack, status | after: model, pred, margin, slack, status | max weight change | category | `[tol]` row | forward check |
|---|---|---|---|---|---|---|---|
| v5pool | Jet A | refutas_raw 2.918, 0.415, 0.447, pass | refutas+ml_pair_gn 3.172, 0.287, 0.471, pass | 0.0000 | B | pass | pass |
| v5pool | Jet A-1 | refutas_raw 2.994, 0.415, 0.433, pass | refutas+ml_pair_gn 3.156, 0.287, 0.474, pass | 0.0000 | B | pass | pass |
| ext_n30 | Jet A | refutas_raw 2.968, 0.415, 0.438, pass | refutas+ml_pair_gn 3.255, 0.287, 0.458, pass | 0.0000 | B | pass | pass |
| ext_n30 | Jet A-1 | refutas_raw 3.015, 0.415, 0.429, pass | refutas+ml_pair_gn 3.197, 0.287, 0.467, pass | 0.0001 | B | pass | pass |
| ext_n8 | Jet A | refutas_raw 2.968, 0.415, 0.438, pass | refutas+ml_pair_gn 3.255, 0.287, 0.458, pass | 0.0000 | B | pass | pass |
| ext_n8 | Jet A-1 | refutas_raw 3.015, 0.415, 0.429, pass | refutas+ml_pair_gn 3.197, 0.287, 0.467, pass | 0.0001 | B | pass | pass |
| guard4_v5pool | Jet A | refutas_raw 2.965, 0.415, 0.438, pass | refutas+ml_pair_gn 3.164, 0.287, 0.473, pass | 0.0000 | B | pass | pass |
| guard4_v5pool | Jet A-1 | refutas_raw 3.025, 0.415, 0.427, pass | refutas+ml_pair_gn 3.146, 0.287, 0.476, pass | 0.0000 | B | pass | pass |

R-DV-2: the `[tol]` row passes in all 8 designs (it failed in all 8 before), so every forward check now passes. The ML value is 5-10 % above Refutas
at -20 C; viscosity stays `pass` with the smaller ML margin; recipes are unchanged to 1e-4. As expected, all 8 designs stay category B and 0 are
certified: smoke blocks every design (aromatics 0 vol %, out of the regression domain, no margin). Freezing `pass` in these designs uses the 16.98 C
margin measured on binary liquidus data (delta_freeze_margin), transferred to these multi-component recipes, not validated for them.
Tests: `tests/test_v5_design_forward.py` (incl. `test_design_viscosity_is_api_designated_with_ml_margin`), `tests/test_predict_route.py`,
`tests/test_delta_freeze_margin.py`: 77 passed (targeted only, as instructed).
