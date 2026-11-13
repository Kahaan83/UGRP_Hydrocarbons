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
