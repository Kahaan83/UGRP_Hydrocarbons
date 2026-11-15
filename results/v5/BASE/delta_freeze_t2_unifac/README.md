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
