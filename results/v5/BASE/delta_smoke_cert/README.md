# delta_smoke_cert: smoke point path to certification (WS1)

Plan: `docs/superpowers/plans/2026-11-13-v5-certification-and-ml-roadmap.md` (WS1, Tasks 1.1-1.6).
Written 2026-11-13 BEFORE any result. Gates below are never moved after results are seen.

## Declared readings

- **R-SC-1 (descriptive):** number of AJFTD fuels with smoke point, aromatics and naphthalenes all present.
- **R-SC-2:** smoke margin, alpha = 0.10, direction `over` (pred - measured, over-prediction is unsafe), refit on in-domain real jet fuels (aromatics 8-23 vol %, uncensored smoke).
  - Grouping is chosen now, knowing the levels. `source_group` has 10 levels on fuels with aromatics (JAA 672, JA 452, JFK 348, LAB_REPORT 271, NUMERIC_ID 190, OTHER 170, PAREN_BATCH 65, HR 49, CRC 2, JETSCREEN 2), so it cannot carry a >= 30-group gate.
  - (i) Split conformal over individual fuels (unit `fuel_id`), with leave-one-`source_group`-out fits, so no fuel is scored by a model that saw its group. (ii) The margin is also computed per held-out `source_group`.
  - **Certifying margin = max(fuel-level margin, largest per-group margin over groups with >= 40 fuels).**
  - Gate: at least 200 in-domain fuels in at least 5 groups with >= 40 fuels, else the margin stays `inf`.
- **R-SC-2b (user decision, not automatic):** if R-SC-2 passes its gate, the user decides whether the real-fuel margin REPLACES the 7.8 mm lab-mixture margin for in-domain jet recipes. `_combined_margin` takes the max over keys, so an extra key alongside 7.8 mm changes nothing. The 7.8 mm margin was the user's own choice, so replacing it needs a new quoted answer. Without that answer the smaller margin is not used.
- **R-SC-3:** the 4 existing design runs (seed 42, same pools) plus 2 new aromatic pools, re-run with `--synthetic-blend --alt-limits` (aromatics min 8 vol %, D1655 alternative smoke limit 18 mm with naphthalenes <= 3.0 vol %). Per design: aromatics vol %, smoke pred, margin, max margin that certifies (pred - 18), smoke status, freezing slack, category. A reading, not a unit test.
  - Stated now: smoke passes only where the replacing margin is below pred - 18 at the recipe's aromatics (6.97 mm at 8 vol %, 4.71 mm at 15 vol %) and the class-mix gate holds. Butylcyclohexane-heavy recipes likely fail the gate; freezing may fail or block (feeds WS2).
- **R-SC-4:** with every new option off, new code gives identical constraint results to pre-change code on the same in-memory weights (stored CSV weights are rounded to 4 decimals, so CSV comparison is not bit-exact).
- **Stated now:** certification also needs the user's `verification` flag upgrade. Without it, best outcome = "would certify; blocked only by `secondary_sources_only`". The agent never edits the flag.

## Declared data rules (written before any result, after reading the raw value strings, not any fit)

- Source: `Desktop/extras/fuelblend_data/data/raw/ajftd/ajftd_all_metrics.csv` (91,175 rows, 2,679 fuels; outside the repo, never committed). It already holds `smoke_point` (2,512 rows, 2,379 fuels, mm) and `naphthalenes` (58 rows, 56 fuels; vol% by D1840 or mass% by D2425). No re-scrape and no HTML parse needed.
- Smoke value strings: a leading `>` or `<` is censored (`smoke_qualifier`, value NaN). `No Test`, `Not Found`, blank give no row. Plain numbers are kept as measured. Plain values above 40 mm and the odd method label `420` are kept and recorded, not filtered by hand; the in-domain fit uses fuels with aromatics 8-23 vol % as declared.
- Repeats: one value per `fuel_id`, the median of its uncensored rows. Fuels with only censored rows are counted and excluded.
- `naphthalenes` (class total, vol% only) is not merged with `naphthalene` (single compound); mass% rows are not converted. Unit is kept in the output.
- Sensitivity only (cannot change the certifying margin): the margin is also reported with fuels at smoke == 25.0 exactly removed, because 25 may be a reporting cap.

## Declared before Task 1.4 results (added 2026-11-13, no fit run yet)

- `smoke_point_from_bulk` docstring says it was trained on 100 AJFTD jet fuels, so scoring those fuels with its stored coefficients would leak. Declared: the certifying margin uses the same model form refit on the training folds only (leave-one-`source_group`-out), predicting each held-out fuel from aromatics only (diaromatics unknown, as in the design call), residual = pred - measured. The stored-coefficient residuals are reported as a leaky reference only.
- Per-group margin: groups with >= 40 in-domain fuels; fuel-level margin = conformal (1-alpha) quantile, alpha = 0.10, of out-of-fold residuals over all in-domain fuels, `ceil((n+1)(1-alpha))/n` rule.
- Class-mix gate (Review Focus 6): real-fuel margin applies only if the recipe's cycloalkane and n-/iso-alkane mass fractions lie inside the [min, max] range of the AJFTD GCxGC fuels (113 mapped); else smoke is `uncertified`, reason `class_mix_outside_real_fuel_range`. Range taken as observed min/max over those fuels, fixed before any design is scored.

## User decisions (quoted verbatim)

R-SC-2b and the yaml key, answered 2026-11-13 in the parent session (relayed to this session by that session; the user chose from the options shown):

- Q: "Should 4.11 mm REPLACE your earlier 7.8 mm lab-mixture margin for in-domain jet recipes inside the real-fuel class-mix range? (Butylcyclohexane-heavy recipes stay uncertified. With 4.11 the 18 mm limit can certify only at aromatics <= 15 vol %. With 7.8 every design fails smoke.)" A: "Yes, replace (Recommended)" (option text: "Use 4.11 mm for in-domain jet recipes; 7.8 mm stays for everything else.")
- Q: "May WS1 add the key smoke_point_mm.ajftd_v5_in_domain.bulk_regression = 4.11 to configs/uncertainty_v5.yaml?" A: "Yes, add it (Recommended)"

Done: key added (value 4.107069860688682, generated by `scripts/v5/base/build_uncertainty_v5_q2.py`, `add_smoke_cert`). Task 1.0 is still the user's: spec `verification:` flags are unchanged and were never edited by this work.

## Results

### Task 1.1 ingest (R-SC-1 counts; `r_sc_1_counts.json`, table `data/real_fuels/v5/ajftd_smoke_naph.csv`)

- Table has 2314 fuels (smoke, naphthalenes or aromatics). 2379 fuels have a smoke row in the raw file; 2183 have an uncensored value, 14 are censored-only, the rest are `No Test`/`Not Found`.
- smoke + aromatics: 2092; smoke + aromatics + naphthalenes vol% (R-SC-1): 36; in domain (aromatics 8-23, smoke uncensored): 1990.
- `source_group` counts among fuels with aromatics match the declared table exactly (JAA 672, JA 452, JFK 348, LAB_REPORT 271, NUMERIC_ID 190, OTHER 170, PAREN_BATCH 65, HR 49, CRC 2, JETSCREEN 2).
- Plain smoke values are kept as recorded, including odd ones (e.g. fuel 156 has a row of -24 mm; >40 mm values).


### Task 1.4 R-SC-2 results (script `scripts/v5/base/delta_smoke_cert.py`; `margin_fit.csv`, `readings.json`)

- Model: production form (`fit_bulk`, OLS) refit on the training folds with **aromatics as the only feature**, because diaromatics is absent from the AJFTD smoke table (the stored model imputes it at its training mean; a single-feature OLS is the matching fold-level model). Leave-one-`source_group`-out, residual = pred - measured, 14 censored-only fuels excluded.
- **Gate (>= 200 in-domain fuels in >= 5 groups with >= 40 fuels): PASSED.** 1990 in-domain fuels (1989 in eligible groups), 7 eligible groups: JAA 665, JA 432, JFK 348, LAB_REPORT 269, OTHER 140, NUMERIC_ID 70, PAREN_BATCH 65. JETSCREEN has 1 in-domain fuel (margin inf, not eligible). HR (49) and CRC (2) have no in-domain fuels (aromatics outside 8-23 or absent).
- Fuel-level margin (alpha 0.10): **2.95 mm**. Per-group: JFK 4.11, LAB_REPORT 2.77, JAA 1.79, OTHER 1.55, PAREN_BATCH 1.18, JA 0.16, NUMERIC_ID 0.07.
- **Certifying margin = max = 4.11 mm (JFK)**, 0.53 x the 7.8 mm lab margin. Out-of-fold coverage of that margin: 0.981. OOF residual mean -0.06, median 0.03 mm.
- Against the max margin that certifies (pred - 18): 4.11 < 6.97 (8 vol %), 6.33 (10), 5.68 (12), 4.71 (15): certifies; 4.11 > 3.09 (20 vol %): does not.
- Sensitivity (smoke == 25.0 removed, n = 1791): certifying 3.86 mm (fuel-level 2.73, JFK 3.86). Cannot change the certifying margin by declared rule.
- Leaky reference (stored coefficients on the same 1990 fuels, trained on 100 of them): fuel-level 2.32, certifying 3.16 mm (JFK). The leak-free margin is larger, as expected.
- Class-mix gate (`fuelblend.design.smoke_margin_v5`), observed range over the 113 mapped GCxGC fuels, class mass fractions normalised to 1 per fuel (cycloalkane = mono+di+tri-cycloparaffin): cyclo 0.000-0.709, n-paraffin 0.000-0.508, iso-paraffin 0.007-0.995. A 75 % butylcyclohexane recipe is outside the cyclo range -> `uncertified`, reason `class_mix_outside_real_fuel_range`.
- Gate interpretation used: the 200-fuel count is taken over fuels in groups with >= 40 fuels (the stricter reading); the total-count reading gives the same pass.
- R-SC-2b not answered: `configs/uncertainty_v5.yaml` NOT edited; the 7.8 mm margin stays in force until a quoted user decision.

### R-SC-2c (consistency check, declared before result; added after Task 1.4 review, before this check was run)

Residuals of `smoke_point_from_bulk(aromatics, None)` (the predictor the design and `predict_route` use) on in-domain fuels NOT among the 100 training fuels of `data/processed/v4/sooting/smoke_point_ajftd.csv` (matched by `fuel_id`; the file has ids; the number matched is reported), as fuel-level conformal margin and per-group margin (groups >= 40), same rule as R-SC-2 (alpha 0.10, over, certifying = max, same gate). If the production-predictor certifying margin is LARGER than the refit margin, the certifying margin used downstream (`certifying_margin_final`) is the larger of the two (the margin must belong to the predictor the design uses); if smaller, the refit margin stays. Also reported: max abs difference between refit-OLS and `smoke_point_from_bulk` predictions over aromatics 8..23 step 0.5. The max-margin-that-certifies levels are derived (pred - 18 at the aromatics, with pred from `smoke_point_from_bulk`) instead of hard-coded.

#### R-SC-2c results

- Matched by `fuel_id`: all 100 training fuels are among the 1990 in-domain fuels; 1890 scored. Eligible groups (>= 40): JAA 665, JA 432, JFK 348, LAB_REPORT 268, OTHER 106, NUMERIC_ID 70 (PAREN_BATCH drops below 40 as most of it is training data); gate passed.
- Production predictor (`smoke_point_from_bulk(aromatics, None)`): fuel-level margin 2.41 mm, per-group JFK 3.16, LAB_REPORT 2.19, JAA 1.74, OTHER 1.37, JA 0.38, NUMERIC_ID -0.05; certifying margin **3.16 mm**.
- Refit margin 4.11 mm is larger, so `certifying_margin_final` = **4.11 mm** (refit stays).
- Max abs difference refit-OLS vs `smoke_point_from_bulk` over aromatics 8..23 step 0.5: **1.87 mm** (the two predictors differ in slope/offset; the refit margin is not calibrated to the stored predictor's bias, but the stored predictor's own margin is smaller, so the larger refit value covers both).
- Max-margin-that-certifies levels in `readings.json` are now derived as `smoke_point_from_bulk(a) - 18`.

### Task 1.5 notes (recipe-dependent margins in the optimiser; appended 2026-11-13)

- **Change to `alt_limits` output.** With the alt flag on (`use_alt_limits`) the reported smoke slack is now max(primary slack, alternative-branch slack), where the alternative branch = min(slack of the 18 mm limit with the margin, slack of every requirement, here naphthalenes <= 3.0). This is the number the optimiser maximises, and the engine, `_Sides`, the rejudge and `robustness_table` all report it (agreement to 1e-9). Alt-on smoke slacks in earlier outputs were primary-only, so they are not directly comparable; alt-off outputs are unchanged (R-SC-4, frozen baseline test).
- **Status rule under the alt flag.** Overall status = better of (primary, alternative branch); the alternative branch = worse of (alternative limit, requirements). A recipe that meets the alternative limit with its margin on a spec whose `verification` is `secondary_sources_only` is `uncertified` on every path, never `fail`. `fail` remains for a violated limit (with margin) on both branches.
- **v5 pool class-mix rejection.** The real-fuel smoke margin (4.107 mm, user decision R-SC-2b 'Yes, replace (Recommended)') applies only to recipes inside the GCxGC class-mix range (cyclo 0-0.709, n-paraffin 0-0.508, iso-paraffin 0.007-0.995). The v5 pool has no iso-alkane, so every v5-pool recipe has iso-paraffin fraction 0 and is outside the range: its smoke margin is None, reason `class_mix_outside_real_fuel_range`, status `uncertified`. That is a property of the pool, not of the optimiser; a pool with an iso-alkane (for example isooctane) can reach the in-range case.
- A requirement that is not `by_construction` and has no margin is no longer treated as margin 0 in `_Sides` (the alternative branch is infeasible); naphthalenes is `by_construction`, so margin 0 is legitimate there and in the rejudge.

### R-SC-3 run protocol (declared before results; written 2026-11-13 before any Task 1.6 run)

Common arguments for every run: `--synthetic-blend --alt-limits --recipe-margins`, seed 42 (default), `--spec astm:jet_a,bis:jet_a1`, default n_starts (30 unless stated), margins `configs/uncertainty_v5.yaml` (contains the user-decided real-fuel smoke key 4.107 mm). Python `.venv/Scripts/python.exe`, `PYTHONPATH=src`. Output root `results/v5/BASE/delta_smoke_cert/designs/<name>/`.

The 4 existing pools (same arguments as `delta_design_ml_visc`, plus the three flags):

```
python scripts/v5/run_design_userpool_v5.py --pool v5pool --synthetic-blend --alt-limits --recipe-margins --out results/v5/BASE/delta_smoke_cert/designs/v5pool
python scripts/v5/run_design_userpool_v5.py --pool v5pool --freezing-guard-K 4 --synthetic-blend --alt-limits --recipe-margins --out results/v5/BASE/delta_smoke_cert/designs/guard4_v5pool
python scripts/v5/run_design_userpool_v5.py --pool ext --n-starts 30 --synthetic-blend --alt-limits --recipe-margins --out results/v5/BASE/delta_smoke_cert/designs/ext_n30
python scripts/v5/run_design_userpool_v5.py --pool ext --n-starts 8 --synthetic-blend --alt-limits --recipe-margins --out results/v5/BASE/delta_smoke_cert/designs/ext_n8
```

The 2 new pools (v5pool members: n-decane, n-undecane, n-dodecane, butylcyclohexane, decalin_trans, 1,2,4-trimethylbenzene):

```
python scripts/v5/run_design_userpool_v5.py --pool "n-decane,n-undecane,n-dodecane,butylcyclohexane,decalin_trans,1,2,4-trimethylbenzene,toluene,isooctane" ...  --out .../designs/poolA_toluene_isooctane
python scripts/v5/run_design_userpool_v5.py --pool "<v5pool members>,hexylbenzene,isooctane" ...  --out .../designs/poolB_hexylbenzene_isooctane
```
(`1,2,4-trimethylbenzene` contains commas, so the pool is passed as a `.txt` file, one id per line, written to `designs/poolA_ids.txt` / `designs/poolB_ids.txt`; the run uses `--pool results/v5/BASE/delta_smoke_cert/designs/poolA_ids.txt`.)

Pool choices, fixed before running (from `data/processed/v5/pure_props_v5.csv` constant statuses and `library_route.build_route_models`):
* **Correction to the expectation given in the task brief.** The existing v5pool is NOT aromatic-free: it already contains 1,2,4-trimethylbenzene (class `aromatic`), so `blocked_aromatics_min_infeasible` will NOT fire for the 4 old pools. They will be optimised. Their expected outcome is instead smoke `uncertified`, reason `class_mix_outside_real_fuel_range` (no iso-alkane, so iso-paraffin fraction 0 < 0.007). This is stated before any run.
* **Pool A = v5pool + toluene + isooctane.** 1,2,4-trimethylbenzene is already in v5pool, so the new aromatic blendstock is toluene (rho15, flash `verified`; CORRECTED 2026-11-14: toluene is NOT in `data/raw/sle_constants/tm_dhfus_verified.csv`, so its Tm/dHfus come from the Coutinho table, which is non-verified: freezing is blocked for pool A as well, not only pool B). Isooctane (2,2,4-trimethylpentane) is added because the class-mix gate needs an iso-alkane fraction >= 0.007 (the minimum over the 113 GCxGC fuels) and v5pool has none; isooctane has `verified` rho15, Tm, dHfus, MW (it is in `tm_dhfus_verified.csv`; flash is supplied by the registry as `verified` for the 2,2,4-trimethylpentane entry). Isododecane is NOT used (rho15 `disputed`).
* **Pool B = v5pool + hexylbenzene + isooctane.** Hexylbenzene is the only C10-C12 alkylbenzene in the verified-constants table that also passes the own-data gate (grid, flash, SLE). Its constants are NOT all verified: rho15, Tm, MW, Tb `verified`; dHfus `disputed`, flash `reference_only`. Consequence declared now: constraints that touch hexylbenzene flash or freezing may be `uncertified` through non-verified constants, and freezing for a recipe with hexylbenzene uses the Coutinho entry not the verified csv. Reported as found, not corrected. The nearest alternative (a registry alkylbenzene from the route) would carry only `route_*` / `coutinho_table` statuses, which is weaker, so hexylbenzene is the pick.
* Run plan: 4 old pools first (each could take 15-40 min per spec; the old pools are no longer instant), the two new pools in the background, concurrently if CPU allows.
* Table columns (r_sc_3.csv): design, category, aromatics vol %, smoke pred, margin used, max margin that certifies (pred - 18), smoke status/reason, naphthalenes, class-mix gate, freezing slack/status, flash, density, viscosity -20 C, NHOC status, and whether the only blocker is `secondary_sources_only` verification. No gate is moved after results.

### Task 1.6 R-SC-3 results (appended 2026-11-13; protocol above unchanged)

Runs: 6 pools x 2 specs, seed 42, `--synthetic-blend --alt-limits --recipe-margins`; each design takes 4-19 s (not 20+ min). Outputs `designs/<run>/`, table `r_sc_3.csv`, builder `scripts/v5/base/delta_smoke_cert_r_sc_3.py`. Smoke margin is the user-decided 4.107 mm where the class-mix gate holds, otherwise none.

| run | design | cat | arom vol % | smoke pred (mm) | margin used (mm) | max margin that certifies (pred - 18) | smoke status | freezing | flash | density | visc -20 C / NHOC |
|---|---|---|---|---|---|---|---|---|---|---|---|
| v5pool | Jet A | A | 15.0 | 22.70 | none (class-mix) | 4.70 | uncertified | fail, slack -0.495 | pass 0.300 | pass | pass/pass |
| v5pool | Jet A-1 | A | 15.0 | 22.70 | none | 4.70 | uncertified | fail -0.570 | pass 0.300 | pass | pass/pass |
| guard4_v5pool | Jet A | A | 15.0 | 22.70 | none | 4.70 | uncertified | fail -0.495 | pass 0.300 | pass | pass/pass |
| guard4_v5pool | Jet A-1 | A | 15.0 | 22.70 | none | 4.70 | uncertified | fail -0.570 | pass 0.300 | pass | pass/pass |
| ext_n30 | Jet A | A | 11.2 | 23.95 | none | 5.95 | uncertified | fail -0.821 | pass 0.399 | pass | pass/pass |
| ext_n30 | Jet A-1 | A | 11.2 | 23.95 | none | 5.95 | uncertified | fail -0.848 | pass 0.399 | pass | pass/pass |
| ext_n8 | Jet A | A | 11.2 | 23.95 | none | 5.95 | uncertified | fail -0.821 | pass 0.399 | pass | pass/pass |
| ext_n8 | Jet A-1 | A | 11.2 | 23.95 | none | 5.95 | uncertified | fail -0.848 | pass 0.399 | pass | pass/pass |
| poolA (toluene, isooctane) | Jet A | A | 20.6 | 20.89 | 4.107 | 2.89 | fail | uncertified, slack +0.155 | fail -0.188 | pass | pass/pass |
| poolA | Jet A-1 | A | 20.6 | 20.89 | 4.107 | 2.89 | fail | fail -0.017 | fail -0.188 | pass | pass/pass |
| poolB (hexylbenzene, isooctane) | Jet A | A | 20.7 | 20.88 | 4.107 | 2.88 | fail | uncertified +0.146 | fail -0.039 | pass | pass/pass |
| poolB | Jet A-1 | A | 20.7 | 20.88 | 4.107 | 2.88 | fail | fail -0.025 | fail -0.039 | pass | pass/pass |

(max margin that certifies = smoke pred - 18; the 4.107 margin is only applied where the class-mix gate holds.)

- **Correction to the stated expectation**: the 4 old pools were not blocked (`blocked_aromatics_min_infeasible` did not fire): v5pool contains 1,2,4-trimethylbenzene. They were optimised; smoke is `uncertified`, reason `class_mix_outside_real_fuel_range` (no iso-alkane), as the Task 1.5 notes predicted. Their smoke predictions (22.7 / 23.95 mm) would pass the 18 mm alternative limit even with 4.107 mm (max margins 4.70 / 5.95 > 4.107) if the class-mix gate held; the gate is what blocks them. On old pools freezing FAILS (all 8 designs), because the 8 vol % aromatics minimum forces trimethylbenzene (high freezing) into the recipe.
- **New pools (with isooctane the class-mix gate holds)**: the optimiser's max-min recipe sits at 20.6-20.7 vol % aromatics, where pred - 18 = 2.88-2.89 mm < 4.107 mm, so smoke is `fail` (limit violated after the margin; alt_min 18 not passed). Pool A flash fails with margin (slack -0.188, API point estimate also fails 33.1 C vs 38 C); pool B flash fails with the margin (-0.039) but not at the point estimate. Freezing: Jet A `uncertified` (passes on the margin transferred from binaries, constants include non-verified entries: pool A has isooctane/toluene on route/table statuses, pool B hexylbenzene dHfus disputed); Jet A-1 `fail` (-0.017/-0.025).
- Aromatics vol % is `uncertified` in every design only because its spec limit has `verification: secondary_sources_only`. No design is "blocked only by secondary_sources_only": every one also has at least one `fail` (freezing, flash and/or smoke) or a class-mix `uncertified`.
- Nothing certifies. 0 of 12 designs are certified; all 12 are category A.
- Recipes in pool A and pool B are identical to 4 decimals except the aromatic blendstock name (toluene vs hexylbenzene weight 0.0394); not investigated (optimiser result at a shared binding set).
- The viscosity constraint model label in pools A/B is `refutas_raw` with margin 0.287 ln, whereas the `delta_design_ml_visc` v5pool designs reported `refutas+ml_pair_gn`; not investigated here (viscosity passes in all 12 designs).

### Task 1.6 review checks (appended 2026-11-13 after the coordinator's three checks; facts first)

- **Result still: 0 of 12 designs certified.** Category `A_fails_a_limit_on_design_or_API_point_estimate` (R4.category, `scripts/v4/run_designs_v4.py:82`, and `cat_final` in `run_design_userpool_v5.py`) means "a limit is violated at the point estimate or by design/API"; it is a FAILURE category, not a pass. Only `C_certified` is a pass, and none has it.
- **CHECK 1 (viscosity label).** Real finding, label only. In `recipe_margins=True` runs the engine reports the property-model name `refutas_raw`, but the value is the ML designated value (`_visc_designated_fn`) and the margin is the ML one: evaluate_recipe on the pool A recipe gives 2.6306 cSt, margin 0.28748, identical with recipe_margins off (label `refutas+ml_pair_gn`) and on (label `refutas_raw`); `visc_ml_status` says in_domain, g = +0.0909 ln (Refutas alone would be 2.40 cSt). The ML model is still called in recipe_margins runs. Not an unsafe-side bug; fixed the label (`_label_viscosity_ml`, `library_route.py`, test `test_recipe_margins_viscosity_label_is_ml_when_in_domain`). The viscosity row in `delta_design_ml_visc/designs/v5pool` (3.172 cSt, 0.2875) vs here (3.163 cSt, 0.2875) differs only because the recipe differs. Run time 4-19 s vs 20+ min: the ML pair matrix is memoised since `d1b2b60`, and the old runs predate it; arguments (n_starts 30/8, seed 42, same pools) match the old README.
- **CHECK 2 (A and B weights identical).** Models differ (`models['pool']` ends in toluene vs hexylbenzene; the .txt ids parsed correctly). The recipe equals the Dirichlet start no. 27 of seed 42 (max difference 4.6e-5): SLSQP did not move off a random start, in either pool, so the two pools return the same weights for different compounds. The returned recipes are NOT optimised optima.
- **CHECK 3 (optimiser).** `best_min_slack` is -100 (`BAD_SLACK`) for all 8 old-pool designs, i.e. no start was feasible: without an iso-alkane the class-mix gate makes the smoke margin None for every recipe, so the smoke side is a constant -100 and the objective is flat. Those old-pool designs are therefore the equal-weight start (all weights 1/6 or 1/8), not optimised, and their freezing failure says nothing about the best recipe. For pools A and B, 26 of 30 starts hit `BAD_SLACK` (class-mix or the smoke domain), the rest are feasible with min slack -0.64 .. -0.19; best -0.1878 (A, flash) and -0.0682 (B, smoke). `n_converged` 30/30 is the SLSQP success flag and does not mean an optimum. Internal vs post-hoc smoke slack agree (A: -0.0677 vs -0.067639).
- **Conclusion limits.** "0 of 12 certified" is true for the returned recipes. It is NOT shown that no recipe in these pools could certify, because the optimiser is stuck at plateaus (`BAD_SLACK` cliffs); that would need a different search (not run; any probe would be declared here first).
- **Pool B is not usable for certification:** hexylbenzene has dHfus `disputed` and flash `reference_only`; its results are reported, but are not a basis for any certification claim.
- Outside the class-mix range the recipe is `uncertified` (the 7.8 mm lab margin is not used as a fallback), per the user's question text and Review Focus 6; 7.8 mm could not certify anything in domain anyway (the largest margin that certifies is 6.97 mm at 8 vol % aromatics).

### Optimiser plateau: declared change and diagnostic (written 2026-11-14 BEFORE any re-run; review IMPORTANT-2)

Finding being addressed: `engine_v4._Sides.vector` gives `BAD_SLACK` (-100) when a recipe-dependent margin is `None`. 26 of 30 starts in pools A/B and every start in the 4 old pools sat on that constant, SLSQP had no gradient and returned a Dirichlet start (start 27) or the centroid.

Change (OPTIMISER ONLY, active only with `--recipe-margins`; judge unchanged):
- New optional fields on `PropertyModel`: `margin_stand_in` (float) and `margin_penalty_fn(pool, w, T_C) -> float >= 0`. New argument `graded_objective` (default False) of `_Sides` and `engine_v4.design`; `library_route.design_routed` and the runner set it True iff `recipe_margins` is on. With it off, or for a model without those fields, the code path is the old one (flag-off results bit-identical; the frozen baseline tests must still pass).
- Smoke model: stand-in margin = the real-fuel key `smoke_point_mm.ajftd_v5_in_domain.bulk_regression` (4.107 mm; the margin that applies inside the certifiable region, so the objective is continuous at its edge; if the key is absent, the lab key). Penalty = (a) distance to the class-mix box (sum over cyclo / n-paraffin / iso-paraffin of how far the mass fraction lies outside the observed GCxGC range), plus (b) aromatics vol % outside 8-23 (distance / 23) plus smoke prediction outside 19-25 mm (distance / 25), with a floor of 0.01 so a `None`-margin recipe is never scored as clean. When the margin is `None`, the slack of the smoke sides is computed with the stand-in and reduced by the penalty.
- Viscosity keeps `BAD_SLACK` for a `None` margin (no stand-in): its margin is `None` only when no margin entry exists.
- The judge is untouched: `engine_v4.evaluate_blend` / `_evaluate_primary` still treat `None` as `uncertified`, never margin 0. The reported `best_min_slack` becomes the graded objective value (flagged `graded_objective: true` in `designs.csv`); the post-hoc constraint statuses and slacks are the judge's.
- Re-run: the same six pools and arguments as the protocol above. Expected (stated now): old pools stop returning the centroid; pools A/B move off start 27; nothing certifies (specs `secondary_sources_only`/unverified, class-mix, freezing/flash constants). If the optimiser still does not move, that will be said plainly.
- Not a probe: no new pool or cap is added.

- Flag-off outputs changed through the yaml key, not through code (review MINOR-6): adding `smoke_point_mm.ajftd_v5_in_domain.bulk_regression` to `configs/uncertainty_v5.yaml` (user decision R-SC-2b) changes in-domain smoke results of runs without the new flags; the code with every new option off is unchanged (R-SC-4).
- Corrections of 2026-11-14 (review): (1) pool A freezing is blocked by non-verified Coutinho Tm/dHfus for toluene as well as pool B (hexylbenzene); the `uncertified` freezing status of the Jet A rows of both pools is therefore expected, not a margin result. (2) `_rejudge_smoke` no longer turns an alternative-met smoke row on a verified spec into `pass` when the engine reported another blocker (extrapolation, non-verified constants, margin direction, validated range). (3) The stale 7.8 mm wording in the smoke reason, the API note and `docs/model_decision.md` now names the margin actually applied. (4) `margin_fn_consistency` and `graded_objective` are written to `designs.csv`.

### Task 1.6 results after the graded objective (appended 2026-11-14; replaces the 2026-11-13 table as the current reading)

Same six pools and arguments as the protocol; `graded_objective: true` and `margin_fn_consistency: []` (consistent) are in every `designs.csv`; `r_sc_3.csv` is regenerated. Design times are now 2-11 min per spec (the graded objective evaluates the ML viscosity everywhere).

- **The optimiser moves and finds feasible-in-the-objective recipes in all 12 designs** (`best_min_slack` +0.0136 .. +0.0209, was -100 / a random start). Recipes are now butylcyclohexane-rich with 8.2-8.4 vol % aromatics (the 8 vol % minimum is tight) and are different between pools; pools A and B are again identical to each other, because toluene, hexylbenzene and isooctane all end at weight 0.
- **Which constraint binds:** the objective is a max-min, so the binding set ties at the equalised slack: freezing (old pools; Jet A-1 slightly lower), NHOC, aromatics vol % minimum (and density in pools A/B). Pools A/B: density, NHOC and freezing/aromatics tie at 0.0205-0.0209.
- Per design (Jet A / Jet A-1 give the same smoke row): aromatics 8.2-8.4 vol %, smoke pred 24.86-24.90 mm, margin none (class-mix gate: iso-alkane fraction 0 in every recipe), max margin that certifies 6.86-6.90 mm (pred - 18), smoke `uncertified`; freezing `pass` with slack 0.014-0.018 in the old pools (0.100-0.118 in guard4), flash `pass`, density, viscosity -20 C and NHOC `pass`; pools A/B freezing `uncertified` (toluene / hexylbenzene constants from the Coutinho table), slack +0.021.
- **No design is certified: 0 of 12.** Blockers: (i) smoke `uncertified` on every row, because the recipe is outside the class-mix range (iso-paraffin fraction 0 < 0.007) and the spec limit is `secondary_sources_only` / `UNVERIFIED`; (ii) aromatics vol % `uncertified` (`secondary_sources_only`); (iii) pools A/B freezing `uncertified` (non-verified constants). A design is "blocked only by `secondary_sources_only`" in none of the 12 rows: the class-mix reason is always present too.
- **Why the optimiser stays outside the class-mix box (plainly stated):** it moves, but not to a certifiable recipe. The smoke stand-in margin (4.107 mm) is the same inside and outside the box, the penalty is the 0.01 floor plus the distance, and the floor makes a jump of 0.01 in slack at the box edge, which SLSQP cannot cross; pools A/B drop isooctane to 0 although 0.7 % would cost almost nothing. So "0 of 12 certified" is still not "no certifiable recipe exists". The follow-up that is NOT run here: a continuous penalty without the jump (or an inward bound on iso-paraffin >= 0.007 in pools that contain an iso-alkane), declared before running.
- **Category A:** now caused only by the point-estimate check against the PRIMARY smoke limit (24.9 < 25 mm by 0.1 mm); `R4.point_estimate_view` ignores the D1655 alternative 18 mm limit, which is nominally met. It is a failure label for the point-estimate view, not a certified or passing state.
- ext_n8 reports 0/8 SLSQP successes (the flag only), recipes are the same as ext_n30 to 1e-3.

### Declared before running (2026-11-14, review 3): label fix, continuous penalty, nudge diagnostic

Declared BEFORE any of the following is run.

1. **Why slack > 0 and category A coexisted.** The optimiser (and the reported `slack` of the smoke row, 0.38) uses the D1655 alternative branch: slack = max(primary, min(alt 18 mm - margin, naphthalenes)). `R4.point_estimate_view` and the forward check's `design_meets` / `api_meets` test only the PRIMARY 25 mm limit on the point estimate (24.9 < 25), so the same recipe is feasible for the optimiser and "fails" in the category. That is a labelling defect: with the alt branch met at the point estimate (smoke >= 18 and naphthalenes <= 3.0 vol %, both design and API values) the smoke result is `alt met, uncertified` (class-mix gate plus an unverified spec), not failed. Fix in `scripts/v5/run_design_userpool_v5.py` only (v4 `point_estimate_view` untouched): when `--alt-limits` is on, an alt-met smoke row counts as meeting its limit in the point-estimate view and in `design_meets` / `api_meets`; `designs.csv` gets `smoke_alt_branch_met` and the category of a design whose only point-estimate "violation" was alt-met becomes `B_meets_both_not_certifiable` (smoke label `alt met, uncertified`). Nothing becomes `C_certified` by this: `res.passed` is unchanged.
2. **Continuous penalty (replaces the jump).** In the graded objective the penalty is now applied always (inside and outside the box), `pen = d_classmix + d_arom + d_smoke`, where `d_classmix` is the distance to the class-mix box shrunk by an inner buffer of 0.005 mass fraction at every NON-ZERO edge (iso-paraffin lower edge 0.007 becomes 0.012; the 0 edges of cyclo / n-paraffin are not shrunk since a fraction cannot be negative), `d_arom` the aromatics distance to 8-23 vol % (/ 23) and `d_smoke` the distance of the smoke prediction to 19-25 mm (/ 25). The 0.01 floor now applies only when the margin is None and no distance explains it. The buffer is there because an exact-penalty edge lets SLSQP stop on the wrong side by round-off. Judge unchanged (None => uncertified), flag-off unchanged. Then all six designs are re-run with the same arguments as the protocol.
3. **Diagnostic (nudge).** `scripts/v5/base/delta_smoke_cert_nudge.py`: for each returned optimum (the committed 2026-11-14 run first, then the continuous re-run) set every iso-paraffin pool member to 1 wt % (others renormalised) and call `evaluate_recipe` with the same flags as the designs (`--synthetic-blend --alt-limits --recipe-margins`). Reported per design: smoke prediction, status, margin, reason, and the change of every other constraint's slack. A pool with no iso-paraffin member (v5pool, guard4) cannot meet the class-mix gate at all: reported as a real blocker, not an optimiser artifact. Reading: if smoke becomes pass-or-uncertified-only-by-`secondary_sources_only` with the other slacks barely moved, the gate was an optimiser artifact in that pool; if other constraints fail or smoke prediction leaves 19-25, it is a real trade-off.

### Amendment declared before the second re-run (2026-11-14, review 3)

The first continuous-penalty re-run returned the SAME recipes as the committed run (pools A/B again isooctane = 0, best_min_slack 0.0209, identical weights) and the nudge diagnostic ran on those. Cause found: the penalty was subtracted from the smoke rows only, and the optimiser maximises the MINIMUM slack; the smoke rows sit at slack 0.38 (alternative branch) while the minimum is about 0.02, so no penalty below about 0.36 can change the optimum. (The earlier 0.01-floor version had the same defect: the movement of the first graded run came from removing the BAD_SLACK plateau, not from the class-mix penalty.) Fix, optimiser only: the total graded penalty of a vector is subtracted from EVERY slack entry (it lowers the whole max-min objective), so `best_min_slack` of a graded run is "min slack minus total penalty" (equal to the plain min slack where the penalty is 0). Test added (a penalty must move the optimum although the penalised row is not binding). All six designs are re-run again with the same arguments; the nudge diagnostic is repeated on the new optima (`nudge_final.csv`; the diagnostic on the committed run is kept as `nudge_committed.csv`).

### Final results (2026-11-14, review 3; replaces the two earlier result sections as the current reading)

Run: continuous graded penalty applied to the whole objective, label fix, same six pools and arguments. `r_sc_3.csv` regenerated (new columns `smoke_label`, `min_slack_constraints`, `binding`). Under the graded objective `best_min_slack` is "min slack minus total penalty", so for the pools without an iso-paraffin member it is lower than the plain minimum slack (a constant class-mix penalty cannot be removed there); `min_slack_constraints` is the plain minimum over the reported constraint slacks.

| pool | recipe | min slack (plain) | binding set | smoke | category |
|---|---|---|---|---|---|
| v5pool, guard4, ext_n30, ext_n8 (no iso-paraffin member) | butylcyclohexane-rich, aromatics 8.3-14.4 vol %, smoke pred 22.9-24.9 mm | 0.010-0.018 | freezing, NHOC, aromatics vol % (the Jet A-1 guard4 recipe: NHOC alone) | `uncertified`: class-mix gate outside the real-fuel range (iso-paraffin 0) AND spec `secondary_sources_only` / `UNVERIFIED` | `B_meets_both_not_certifiable`, smoke label `alt met, uncertified` |
| pool A (toluene + isooctane) | isooctane 1.18 wt %, toluene 0, aromatics 8.35 vol %, smoke pred 24.86 mm | 0.0204-0.0207 | density, freezing, NHOC, aromatics vol % (all tied) | `uncertified`: gate INSIDE (margin 4.107 mm applied, 24.86 - 4.11 = 20.75 >= 18), blocked only by the unverified spec | `B_meets_both_not_certifiable`, `alt met, uncertified` |
| pool B (hexylbenzene + isooctane) | identical to pool A (hexylbenzene 0); NOT usable for certification (hexylbenzene dHfus disputed, flash reference_only) | same | same | same | same |

- **Category meaning.** `B_meets_both_not_certifiable` = every limit is met at the point estimate on the design side and the API side (smoke through the D1655 alternative branch: >= 18 mm with naphthalenes <= 3.0 vol %, naphthalenes 0 by construction) but the design cannot be certified. It is NOT a pass. The earlier category A in all 12 rows was a labelling defect (primary-limit-only point-estimate view against an alt-branch optimiser objective); fixed in the runner via `LR.apply_alt_point_view` (no status and no `passed` changed).
- **Nudge diagnostic** (`nudge_committed.csv`: the optimum of the committed run, iso-paraffin set to 1 wt %; `nudge_final.csv`: the final optimum, already at 1.18 wt %, so a no-op): in pools A/B the nudge moves the smoke row from "class mix outside range, margin None" to "class mix inside, margin 4.107 mm applied" with the prediction 24.86 -> 24.89 mm; slack changes: aromatics vol % -0.0055, density -0.0145, flash -0.12, freezing +0.004, NHOC +0.0002; every other constraint stays `pass`. Reading: in pools A/B the class-mix gate was an OPTIMISER ARTEFACT (it cost about 0.015 slack to satisfy it); the smoke blocker that remains is the unverified spec. In the pools with no iso-paraffin member (v5pool, guard4, ext) the gate cannot be met at all: a REAL blocker for those pools (add an iso-paraffin to the pool).
- **Still uncertified: 0 of 12.** Remaining blockers: smoke and aromatics vol % limits are `secondary_sources_only` / `UNVERIFIED` (never certified, spec flags untouched); pools A/B freezing `uncertified` because the pool contains toluene / hexylbenzene with non-verified Tm, dHfus (the model flags a non-verified member even at weight 0; dropping zero-weight members from the judged pool would clear this and was not done here); the old pools additionally have the class-mix gate.
- Why the first continuous-penalty run did not move: see the amendment above (penalty on a non-binding row).
