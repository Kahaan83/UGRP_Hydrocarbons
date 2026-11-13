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
