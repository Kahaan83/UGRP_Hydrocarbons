# delta_flash_margin_q2: flash lab margin 2.21 C

## Why

User decision 2026-11-13, made after seeing the Task 10 numbers (`results/v5/BASE/delta_lab_flash_li2014/`): "use 2.2 for Q2".
The flash lab margin used by the API and by the designs was 3.698768343688585 C (all 31 old Leenhouts systems, one-sided over-prediction,
alpha 0.10). Task 10 showed that this number is set by one system, mesitylene, and that the margin falls to 0.31 C without it.
The user excluded mesitylene because its Antoine data are suspect. The replacement is the row `designated,old_plus_new_minus_mesitylene` of
`results/v5/BASE/delta_lab_flash_li2014/lab_margin.csv`: **2.206932808388615 C** (35 systems, 150 points, one-sided over, alpha 0.10).

Stated plainly: this is a choice made after the numbers were seen. Removing the system that sets the margin lowers the margin
(3.70 to 2.21 C), so the new value is less conservative and carries a cherry-picking risk. It is still set by one system: the new
ternary n-alkane system, whose maximum over-prediction is +2.2 C. The margin is not a validated coverage statement for novel blends,
and mesitylene-containing blends still carry the warning `antoine_suspect:mesitylene`.

## What changes

* New generated file `configs/uncertainty_v5.yaml` (script `scripts/v5/base/build_uncertainty_v5_q2.py`): a copy of
  `results/v5/S1/uncertainty_v5.yaml` with only `flash_point_C.lab_mixtures.ideal.applied_margin` set to 2.206932808388615
  (plus `applied_margin_v4`, `applied_margin_source`, `n_groups_applied: 35`, `user_decision`).
  `configs/uncertainty_v4.yaml` and `results/v5/S1/uncertainty_v5.yaml` are not modified.
* `src/fuelblend/api/predict_route.py`: the flash margin (and only it) is read from `configs/uncertainty_v5.yaml`
  (max of the lab ideal and the D93 reconstructed ideal applied margins, unchanged keys). Viscosity, density, heat keep their yaml.
  Legacy `predict_blend` v4 paths stay on v4. Warning `lab_margin_set_by_one_system` is kept; `flash_margin_excludes_mesitylene (user decision Q2)`
  is added whenever a finite margin applies.
* `scripts/v5/run_design_userpool_v5.py` and `scripts/v5/base/delta_api_route.py` read the new yaml for the design margins.

## Declared readings

* **R-Q2-1**: `predict_fuel` flash `margin` on an in-domain known-composition hydrocarbon blend equals 2.206932808388615, and lo = value - margin.
* **R-Q2-2**: viscosity, density and all non-flash outputs of `predict_fuel` are bit-identical to before. Input set: the forward-demo 60-component
  blend and a 5-compound known blend, each at 233.15, 288.15 and 313.15 K (snapshot before and after the change; every field of every non-flash property).
* **R-Q2-3**: the designs are re-run (same arguments, seed 42; runs v5pool, ext_n30, ext_n8, guard4_v5pool). Per design: is flash now certified, and the
  overall category. Expectation: still 0 certified overall, because freezing and smoke have no margin.

## Result

(filled in after the runs)
