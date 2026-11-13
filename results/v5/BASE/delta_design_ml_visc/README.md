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
