# Overnight run 2026-11-15: what changed, what it showed, what you need to decide

Branch `v5`. Everything is committed locally; **nothing was pushed** (the push approval in memory covered the 12 Nov run only).
Rules kept: every reading declared in its README before it was computed; no gate, margin, spec limit, spec `verification:` flag, judge
or `configs/uncertainty_v5.yaml` changed; all new behaviour opt-in with the default bit-identical; negative results kept.

## Headline

* **Certification: 2 designs certified by the unchanged judge, both conditional; 0 unconditional.** Pool P0y-b (pure pool P0x-b +
  bicyclohexyl) gives optima that pass every constraint of ASTM Jet A and BIS Jet A-1 (CLM-72). They need the opt-in flash-registry
  overlay (bicyclohexyl Antoine, not adopted) and pass D86 only under the designated Engler route; under the Rayleigh model both miss the
  D7566 slopes. **Jet A is the stronger design** (passes both data sensitivities); Jet A-1 also fails if bicyclohexyl's solid-solid
  transition enters SvL, and loses its flash slack at the lower measured flash cluster. Details and the condition table:
  `results/v5/BASE/delta_pool_bicyclohexyl/README.md`.
* **Two improvements that are ready for adoption (your call):**
  1. **Freezing: UNIFAC-SvL** lowers the freezing margin from **16.98 C to 13.84 C** and the S_F MAE from 3.33 to 2.62 C, with no UNIFAC
     fallbacks. It meets the declared preference rule (`results/v5/BASE/delta_freeze_unifac/README.md`).
  2. **`predict_fuel` now reports ASTM D86** (T10, T50, T90, FBP, T50-T10, T90-T10) with WS7 margins and one-sided safe bounds, a
     Rayleigh variant, and a warning when the two models disagree by more than 20 C on T90-T10. This is reporting only: it does not change
     any other output.
* **Two problems found and documented** (no silent fixes):
  1. The optimiser silently drops a constraint when its model is unavailable for a pool member (`engine_v4._Sides.excluded`). The P0y run
     optimised without flash and drove isooctane to 37-39 wt %. An **opt-in guard** now exists: `design_pool(..., strict_constraints=True)`
     raises instead (default False, unchanged).
  2. The frozen Coutinho lookup table gives 112 aromatics another molecule's constants (97 as benzene); fluorene and phenanthrene land on
     the unsafe side. Audit only (`results/v5/BASE/delta_coutinho_lookup_audit/`).
* A committed reading was slightly wrong: the delta_freeze_margin S_F table had 10 duplicate rows. The 16.98 C margin is unchanged; the
  true ideal-model MAE is 3.33 C, not 3.01 C.

## Deltas (one folder each under `results/v5/BASE/`)

### delta_d86_rayleigh (commit e92d5c6)
A zero-parameter ideal single-stage Rayleigh D86 model (`src/fuelblend/baselines/distillation_rayleigh.py`), opt-in through
`add_d86_models(d86_model="rayleigh" | {quantity: model})`; the default stays Engler.

R-DR-2, 65 real fuels (Engler vs Rayleigh):

| quantity | margin Engler | margin Rayleigh | MAE Engler | MAE Rayleigh |
|---|---|---|---|---|
| T10 | 7.93 | 0 | 4.52 | 11.63 (bias +11.6) |
| T50 | | | 8.75 | 6.18 |
| T90 | | | 13.06 | 7.82 |
| FBP | 3.26 | 6.10 | | |
| T50-T10 | 18.68 | 1.24 | 7.67 | 10.88 |
| T90-T10 | 18.25 | 5.31 | 12.56 | 13.79 |

Not adopted by the declared rule (margin smaller AND MAE not larger). **R-DR-3:** on the pure pools P0, P0x-a and P0x-b, Engler reads
T90-T10 of about 50-55 C and Rayleigh 17-19 C. The Engler number comes from linear interpolation through 4-5 wt % isooctane. The pure pools
are therefore physically far from the D7566 slope requirement (T90-T10 >= 40 C); P1, P2 and P3 meet the D86 limits under both models.

### API: predict_fuel D86 block (commit 4fab6dd)
`out["distillation_D86_C"]`: designated method `engler_tbp_offsets`, WS7 margins, per-quantity `safe_bound`, `in_domain`, `tbp_spread_C`,
`monotone`, `variants["rayleigh"]` with the R-DR-2 margins, `d86_models_disagree` warning, `adopted_in_uncertainty_yaml: False`.
79 API tests pass.

### delta_freeze_unifac (commits 8feb2e7, ff1b2f2)
* R-FU-1 (S_F, 9 systems, 67 rows): ideal MAE 3.33 / margin 16.98; **UNIFAC MAE 2.62 / margin 13.84**. Naphthalene + trans-decalin
  -16.98 -> +0.06. 1-methylnaphthalene + 2,6-dimethylnaphthalene (-13.84, aromatic-aromatic, gamma about 1) now sets the margin.
  Labelled exploratory / dev (the ideal residuals were seen before the declaration; no locked SLE set exists).
* R-FU-2 (40 systems, 643 rows, descriptive): ideal 4.24 / 17.81; UNIFAC 3.92 / 17.23.
* R-FU-3 (stored pure-pool optima, counterfactual): no recipe flips. P0x-b Jet A would be 0.4 C short under UNIFAC; pool C passes under
  both.

### delta_coutinho_lookup_audit (commit 4fab6dd)
112 misassigned aromatics in the frozen table, listed in `misassigned.csv`. No code change.

### delta_pool_bicyclohexyl (commit 96d8ad6 and the final commit)
* R-BC-1: Antoine for bicyclohexyl (thermo LANDOLT, extrapolated below its 424-577 K range) and cis-decalin; Tb 232.5 vs 239 and 187.7 vs
  195.8 C (gate 10 C passed). Stored in the opt-in overlay `data/processed/v5/flash_registry_additions.csv`; default registry unchanged
  (test now checks frame equality with the base build).
* R-BC-2, pool P0y (+ both): **invalid**. cis-decalin's MW tier is unverified, so flash and NHOC raised `ModelUnavailable`; the optimiser
  silently dropped them and put 37-39 wt % isooctane in. Led to the `strict_constraints` guard.
* R-BC-4, pool P0y-b (+ bicyclohexyl only), 10 starts, seed 42, about 25 min per spec:

| | Jet A | Jet A-1 |
|---|---|---|
| judge | certified, min slack 0.0115 | certified, min slack 0.0058 |
| recipe (wt %) | butylcyclohexane 27.7, TMB 17.8, bicyclohexyl 17.3, n-C10 15.1, n-C11 13.6, isooctane 3.0, n-C12 2.9, n-C13 2.6 | bicyclohexyl 35.9, butylcyclohexane 24.7, TMB 14.4, n-C10 9.3, n-C11 7.9, isooctane 3.5, trans-decalin 1.6, n-C12 1.4, n-C13 1.3 |
| flash / freezing (C) | 40.64 / -57.4 | 40.43 / -64.2 |
| Rayleigh T50-T10 / T90-T10 (needs 16.2 / 45.3) | 9.5 / 35.4: `d86_model_dependent` | 14.6 / 44.0: `d86_model_dependent` |
| R-BC-3, dH with transition (13.94 kJ/mol) | -57.4, unchanged | -40.7 (+16.98 > -47): `freezing_dH_choice_dependent` |
| R-BC-5, pure flash at 74 C cluster | bound 38.25, passes | bound 37.82: `flash_data_dependent` |
| UNIFAC-SvL + 13.84 (counterfactual) | -43.1, passes | -49.7, passes |

* Usability check: `predict_fuel` on the Jet A recipe (overlay on) returns the `d86_models_disagree` warning
  ("T90-T10 Engler 70.4 vs Rayleigh 35.4 C"), so a user of the API sees the D86 caveat without reading this report.
* Data-consistency finding: with the registry (Jones) LFL, bicyclohexyl's 74 C measured flash needs Tb 210 C (literature 239 C); the 92 C
  cluster is the one consistent with Tb.
* Opt-in guard `design_pool(..., strict_constraints=True)`: raises when `_Sides` dropped a constraint as unavailable; default off.

## Data search
* Lead, not downloaded: Sequeira et al. 2024, Int. J. Thermophys., doi 10.1007/s10765-023-03317-9 (n-octane/n-decane and n-decane/n-dodecane
  SLE). Open copy: repositorioaberto.uab.pt. Two alkane-alkane systems would not move the freezing margin, which is set by an
  aromatic-aromatic system.
* bicyclohexyl Antoine: the NIST WebBook has no Antoine set (and returned an internal error); thermo's LANDOLT data were used instead
  (see delta_pool_bicyclohexyl).

## Tests run tonight
* Design suites (`test_delta_d86_rayleigh`, `test_delta_smoke_branch_opt`, `test_v5_sak_blendstock`, `test_v5_ws7_design_realism`,
  `test_v5_design_forward`): 42 passed, with the guard in place.
* API and flash suites after the `flash_unifac.py` overlay edit (`test_delta_api_route`, `test_ml_api_variants`, `test_predict_blend`,
  `test_predict_route`, `test_v4_api_known_composition`, `test_base_delta_flash_inputs`, `test_delta_lab_flash_li2014`, `test_s2_flash`,
  `test_v2_a4_flash`, `test_v4_flash`, `test_delta_d86_rayleigh`): 109 passed.
* `test_delta_pool_bicyclohexyl`: 3 passed (overlay off = base frame exactly; overlay on adds the entry; guard raises).
* The full suite was not re-run.

## Decisions needed from you
These are coupled: **adopting the overlay (3) is necessary but not sufficient** for the P0y-b certification to stand. It also rests on
keeping Engler as the designated D86 route (2) and, for Jet A-1, on the fusion-only bicyclohexyl dHfus (4).
1. **Adopt UNIFAC-SvL as the designated freezing model with the 13.84 C margin?** Meets the declared preference rule; labelled exploratory.
   It does not decide either P0y-b design (both pass freezing under it).
2. **D86 route.** Rayleigh is not preferred by the declared rule (T10 bias), so Engler stays designated. But on pure-compound pools Engler
   inflates the slopes by 35-45 C, and under Rayleigh both P0y-b recipes miss D7566 T50-T10 and T90-T10. Do you accept an Engler-only
   D86 pass for a pure-compound recipe, or require both models to pass (which would leave 0 certified)?
3. **Adopt the flash-registry overlay** into the default registry? It could move GCxGC representatives and so frozen T2 / WS7 / P2 / P3
   numbers; it would need a re-score. Without it, bicyclohexyl cannot be in any pool.
4. **Bicyclohexyl dHfus for SvL:** fusion only (6.86 kJ/mol, current) or fusion + solid-solid transition (13.94 kJ/mol)? The transition
   sits about 3 K below the melt and both liquidi are far below it, which favours the total; under it the Jet A-1 design fails freezing.
5. **Fix the Coutinho lookup table** (112 aromatics)? A frozen v5 artefact, so a fix changes frozen numbers.
6. **Make `strict_constraints=True` the default in design runs?** Changes no certified result; turns invalid optima like P0y into errors.
7. **cis-decalin MW tier:** stored as unverified, which blocks flash and NHOC for any pool containing it. Verify it (C10H18, 138.25)?
8. **Push approval** for the overnight commits on `v5` (e92d5c6 .. final).
