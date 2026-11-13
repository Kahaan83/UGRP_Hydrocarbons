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

(none yet)

## Results

### Task 1.1 ingest (R-SC-1 counts; `r_sc_1_counts.json`, table `data/real_fuels/v5/ajftd_smoke_naph.csv`)

- Table has 2314 fuels (smoke, naphthalenes or aromatics). 2379 fuels have a smoke row in the raw file; 2183 have an uncensored value, 14 are censored-only, the rest are `No Test`/`Not Found`.
- smoke + aromatics: 2092; smoke + aromatics + naphthalenes vol% (R-SC-1): 36; in domain (aromatics 8-23, smoke uncensored): 1990.
- `source_group` counts among fuels with aromatics match the declared table exactly (JAA 672, JA 452, JFK 348, LAB_REPORT 271, NUMERIC_ID 190, OTHER 170, PAREN_BATCH 65, HR 49, CRC 2, JETSCREEN 2).
- Plain smoke values are kept as recorded, including odd ones (e.g. fuel 156 has a row of -24 mm; >40 mm values).

