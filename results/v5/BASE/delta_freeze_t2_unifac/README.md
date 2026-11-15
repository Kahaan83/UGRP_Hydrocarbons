# v5 delta_freeze_t2_unifac: designated freezing route (UNIFAC-SvL) on the T2 real-fuel freezing set

Status: DECLARATION, written and committed before the scoring script was run (Task 14 of
`docs/superpowers/plans/2026-11-15-apply-decisions-and-audit-fixes.md`). The Results section is appended after the run.

## Why

Decision 1 (2026-11-15, CLM-73) made UNIFAC Schroeder-van Laar (`svl_unifac`, `DESIGNATED_FREEZING_MODEL` in
`src/fuelblend/eval/composition_route.py`) the designated freezing model of the API and the design judge. Its margin
(13.8427 C, `freezing_point_C.lab_liquidus.svl_unifac` in `configs/uncertainty_v5.yaml`) was derived on 12 lab liquidus
systems (S_F, `results/v5/BASE/delta_freeze_adopt/`). Every real-fuel freezing cell in the scoreboard is still ideal SvL:
`scripts/v5/base/cells_t2_bulk.py::freezing_cells()` aggregates column `I` of the committed `results/v5/S5/s5_predictions.csv`
(rerun: `results/v5/BASE/delta_adopt_rescore/rerun/cells_t2_bulk.csv`, `freezing_point_C,T2_measured,pooled,ideal_svl`:
n = 49 groups, bias = MAE = 7.986 C, trivial 5.936 C, `no_skill_shown`). UNIFAC-SvL has never been scored on real fuels.

## What is scored

Same fuels and rows as the scoreboard ideal-SvL freezing cells: the entities in `results/v5/S5/s5_predictions.csv`
(one per double-count group, S5 selection; headline 50 = AJFTD 13 + AV-33-22 34 + FuelLib 3; locked 8 AJFTD), with S5's
measured `freezing_C`, `group` and `role`. Compositions: S5 `build_entities()` (AJFTD GCxGC G1 mapping, AV-33-22 and FuelLib
S2 representative maps; the same `mapped` table S5 used for `I`), imported by path; `s5_freezing_nalkane.main()` is not called.
Mole fractions built exactly as `s5_freezing_nalkane.ideal_svl_I` (mass fraction / RDKit MolWt, normalised).

Columns (each computed with `composition_route.freezing_C`, `ideal_sle_C` underneath for the ideal model):

| column | model | constants | role |
|---|---|---|---|
| `I_s5` | ideal SvL | S5's own map (`get_pure_thermal_properties`, default substring aromatic match, as `ideal_svl_I`) | reproduction check only |
| `unifac_s5` | UNIFAC SvL | same S5 map | isolates the activity-coefficient effect (sensitivity) |
| `ideal_desig` | ideal SvL | designated map `fuelblend.api.predict_route._tmap_freezing()` (verified where present, else exact-match Coutinho) | reference on the API constants |
| `unifac_desig` | UNIFAC SvL | designated map | **primary: the designated API route** |

Reported per row: `unifac_fallback` (count and mole fraction of members with x > 1e-4 and no UNIFAC groups, gamma = 1),
mapped InChIKeys missing from the designated map.

## Metrics, comparator, verdict

Per cell: n_groups, n_rows, bias, MAE, RMSE (pred - measured, C). Sets copied from `freezing_cells()`: role headline ->
`pooled`, `AV-33-22`, `AJFTD`, `FuelLib`; role locked -> the same names with suffix `_locked`; rows dropped with
`dropna(subset=[col, "C", "freezing_C"])` and counted as `n_excluded_nonfinite`.

Comparator: the scoreboard's trivial model `logo_constant` = S5 column `C` (leave-one-group-out mean of headline training
freezing points; locked fuels scored from the all-headline mean), `err_trivial = C - freezing_C`. Not recomputed.

Verdict: `fuelblend.eval.cellkit.make_row` is called with exactly the `freezing_cells()` arguments, so the paired group
bootstrap (B = 5000, seed 0, 95 % percentile CI of MAE(trivial) - MAE(model)) and the rule are the scoreboard's. Copied
verbatim from `src/fuelblend/eval/group_stats.py`:

```python
def verdict(n_groups: int, lo: float, hi: float, min_groups: int = MIN_GROUPS) -> str:
    if n_groups < min_groups:
        return "n_too_small"
    if lo > 0:
        return "beats_trivial"
    if hi < 0:
        return "worse_than_trivial"
    return "no_skill_shown"
```
with `MIN_GROUPS = 8`, and from `src/fuelblend/eval/cellkit.py::make_row`:
```python
    if row["status"] != "no_evaluation" and np.isfinite(row["n_groups"]) and row["n_groups"] < 8 and row["verdict"] != "not_applicable":
        row["verdict"] = "n_too_small"
```
If a model has a different non-finite count than ideal, an intersection-rows MAE/bias (rows finite for all four columns)
is also reported, descriptive, no verdict.

Margin coverage (one-sided, direction `under`): covered iff measured <= pred + 13.8427 C; coverage = mean over rows of the
cell. For `ideal_*` the same is reported with the ideal margin 16.9755 C (`svl_verified`) for context. Caveats declared now:
(a) this is a *would-be* coverage: the API withholds the margin for every one of these fuels (representative-mapped, not a
`known` composition, not all-verified constants); (b) the margin comes from 12 lab systems; with a positive bias one-sided
coverage is near-trivially high and says little, so bias is printed next to coverage.

Sanity check (must pass before any reading): `I_s5` reproduces S5 `I` to |diff| <= 1e-6 C on every row, and the `I_s5`
pooled cell reproduces the rerun cell (n 49, bias 7.98592, MAE 7.98592, trivial 5.93576, `no_skill_shown`) to 1e-4 C. If it
fails, nothing is read.

## Expectations (stated before the run)

* No claim of improvement is made in advance. Ideal SvL scores `no_skill_shown` pooled (bias +7.99 C, every pooled error
  positive) and `worse_than_trivial` locked; UNIFAC for alkane/aromatic mixtures typically has gamma > 1 for the
  crystallising n-alkanes, which raises the predicted liquidus, so UNIFAC is expected to be no better and probably worse
  (larger positive bias) than ideal on these fuels. This is a reading, not a gate.
* Nothing is adopted or changed: no scoreboard row, frozen BASE file, margin, gate or verdict rule is edited; the
  designated model stays as decided. Results feed a proposed Table A addition reported to the controller.

## Files (written by `scripts/v5/base/delta_freeze_t2_unifac.py`)

`rows.csv` (per fuel: entity, source, group, role, measured, C, I_s5_check, the four columns, fallbacks, missing constants),
`summary.csv` (per set x column: make_row fields, coverage of the margin, intersection MAE/bias), `readings.json`.

## Results (run 2026-11-15, `python scripts/v5/base/delta_freeze_t2_unifac.py`, a few seconds)

Sanity check PASSED: `I_s5` reproduces S5 `I` on all 58 rows (max |diff| 7e-15 C, NaN pattern identical); the `I_s5` cells
reproduce the rerun scoreboard ideal-SvL cells for pooled, pooled_locked, AV-33-22, AJFTD, FuelLib, AJFTD_locked (max |diff|
5e-8 on n, bias, MAE, trivial MAE, CI; verdicts identical). No UNIFAC fallback on any row; no mapped InChIKey lacks a
designated constant; no fuel has all-verified constants (so the API would withhold the margin for all 58).

Main cells (C; pred - measured; trivial = S5 `C`, logo_constant; CI = 95 % paired group bootstrap of MAE(trivial) - MAE(model)):

| set | column | n | bias | MAE | trivial MAE | CI | verdict | coverage (margin) |
|---|---|---|---|---|---|---|---|---|
| pooled | **unifac_desig (designated route)** | 50 | +5.878 | 8.425 | 6.413 | [-4.33, +0.44] | no_skill_shown | 0.98 (13.84) |
| pooled | ideal_desig | 50 | +6.263 | 8.825 | 6.413 | [-4.75, +0.04] | no_skill_shown | 0.98 (16.98) |
| pooled | unifac_s5 | 49 | +7.586 | 7.586 | 5.936 | [-3.51, +0.50] | no_skill_shown | 1.00 (13.84) |
| pooled | I_s5 (= scoreboard) | 49 | +7.986 | 7.986 | 5.936 | [-3.92, +0.12] | no_skill_shown | 1.00 (16.98) |
| AV-33-22 | unifac_desig | 34 | +6.907 | 6.907 | 5.858 | [-2.90, +0.99] | no_skill_shown | 1.00 |
| AJFTD | unifac_desig | 13 | +2.824 | 12.620 | 9.126 | [-10.69, +4.05] | no_skill_shown | 0.92 |
| FuelLib | unifac_desig | 3 | +7.447 | 7.447 | 0.951 | [-7.74, -4.70] | n_too_small | 1.00 |
| pooled_locked | **unifac_desig** | 8 | +4.286 | 9.310 | 5.176 | [-7.71, -0.04] | worse_than_trivial | 0.875 |
| pooled_locked | I_s5 (= scoreboard) | 8 | +4.639 | 9.555 | 5.176 | [-7.88, -0.31] | worse_than_trivial | 0.875 |

Like-for-like (the 49 pooled rows finite for all four columns; descriptive, no verdict): MAE unifac_desig 7.297, ideal_desig
7.698, unifac_s5 7.586, I_s5 7.986, constant 5.936. Every one of these 49 errors is positive for all four columns (the
models over-predict the freezing point of every headline fuel but AJFTD:325).

Reading (honest, no spin):
* **The designated UNIFAC-SvL route does not beat the constant on real fuels.** Pooled: `no_skill_shown`, MAE 8.42 C vs 6.41 C
  for the constant, the CI almost entirely on the worse side. Locked AJFTD fuels: `worse_than_trivial` (9.31 vs 5.18 C), as for
  ideal SvL. The real-fuel verdicts of the scoreboard are unchanged by the model switch.
* UNIFAC vs ideal: UNIFAC lowers the predictions by about 0.4 C on average (MAE 7.30 vs 7.70 on the same 49 fuels, designated
  constants; 7.59 vs 7.99 on S5's constants). The declared expectation (UNIFAC no better, likely a larger positive bias) is
  refuted on its second half: the bias is slightly smaller. The difference is small next to the +7 C bias common to both and
  carries no verdict (UNIFAC vs ideal was not a declared comparison with a rule).
* The pooled n rises from 49 to 50 with the designated constants because AJFTD:325 (measured -81 C, S5 constants gave no
  value) now gets a value, -144.7 C (error -63.7 C). That single fuel pulls the pooled bias from +7.3 to +5.9 C and the AJFTD
  MAE from 8.4 to 12.6 C; the constant's pooled MAE moves from 5.94 to 6.41 C on the 50 rows.
* Margin coverage 0.98 pooled (49 of 50; the uncovered fuel is AJFTD:325, under-predicted by 63.7 C) and 0.875 locked
  (AJFTD:97, -61 C measured, under-predicted by 19.7 C). As declared, with a +6 C bias this coverage is near-trivial and does not
  validate the 13.84 C margin on real fuels; the API withholds the margin for all these fuels anyway.
* Nothing adopted or changed. Designated model, margins, gates, verdict rules, frozen BASE files and the scoreboard are untouched.
