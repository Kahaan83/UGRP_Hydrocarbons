# v5 S1b: one-time held-out test of raw Refutas and its lab margins on 2020 ThermoML systems

Declared 2026-11-12, before any residual of the new rows was computed. An earlier same-day draft of this file described using
these rows as dev (training) data; the user cancelled that plan before anything was scored. The readings below replace it. Do not
edit the readings after scoring; add dated notes.

## Data and rule

* New table: `data/real_fuels/v5/thermoml2020_lab.csv` (built by `src/fuelblend/data/loaders_thermoml2020.py`, `sources` =
  `thermoml_jced_2020_unlisted`, measurements_dedup schema).
* Test systems = systems of that table that pass `lab_selection` (S1) and are NOT in the S1 split
  (`results/v5/S1/lab_systems_splits_v5.json`, neither dev nor locked). Systems of the new table that are already in the S1 split are
  not test systems; their new rows are reported separately and labelled "not independent".
* Nothing is fitted. Margins and bias offsets come only from S1 dev systems (new rows excluded). The S1 locked systems are not used
  for any choice here. The test set is scored once. Nothing is written to `configs/`, `splits/` or `results/v5/S1/`.
* Scoring code: the functions of `scripts/v5/eval_lab_known_composition_v5.py` (`score_rows`, `add_bias_corrected`, `method_key`,
  same temperature bins) and `fuelblend.uq.conformal` (`two_sided_margin`, `coverage`, `group_scores`, alpha 0.1, group = system,
  two-sided |ln(pred/meas)|).

## Readings (fixed before scoring)

* **R-S1b-1 (regression).** Margins recomputed here from S1 dev systems only (new rows excluded) equal
  `results/v5/S1/uncertainty_v5.yaml` for raw Refutas in the two cold bins and for bias-corrected Refutas in all three bins
  (absolute difference < 1e-6). If not, stop; nothing else is read.
* **R-S1b-2 (primary, coverage).** Per temperature bin, a test system is "covered" if its largest |ln(pred/meas)| under raw Refutas
  is at most the margin in use for that bin. Warm bin `T_ge_283K`: the v4 value in `configs/uncertainty_v4.yaml`
  (`viscosity_kin_cSt.lab_mixtures.raw_refutas.T_ge_283K`, read from the file); cold bins: the S1 values. With n test systems in a bin,
  the 90 % margin "holds" if covered >= the smallest k such that P(Binomial(n, 0.9) <= k) >= 0.05 (for n = 17 this is 13); it "fails"
  otherwise; "not evaluable" if n < 9 (counts still reported). k is computed with scipy.stats.binom and printed.
* **R-S1b-3 (bias-corrected).** The same coverage test for `refutas_biascorr` with the S1 all-dev offset and the S1 bias-corrected
  margins. Reported next to R-S1b-2; no choice between the two is made from the test.
* **R-S1b-4 (accuracy).** For the test systems, per bin: mean ln residual, RMSE of ln residual, number of systems and rows, next to
  the same numbers for S1 dev and S1 locked systems (from `results/v5/S1/summary_by_bin_method.csv`). No pass/fail.
* **R-S1b-5 (sensitivity, declared now).** Repeat R-S1b-2 without the systems containing an adamantane (tricyclic cage compounds,
  not typical of jet fuel). The headline is the all-systems result; the sensitivity is reported beside it, never instead of it.
* **R-S1b-6 (per-system table).** Components, bin, rows, mean ln residual, max |ln residual|, covered yes/no. Systems dropped by
  `score_rows` (component not usable in the provider/registry or class not in CLASSES) are listed with the reason category
  "not scorable" and counted; they are not counted as covered or uncovered.
* **R-S1b-7.** After this run these systems count as read; they must not be presented again as an untouched test.

A reading that cannot be evaluated is stated "not evaluable" with the reason.

## Dated notes (added after scoring; readings above unchanged)

* 2026-11-12: scored once (`scripts/v5/eval_lab_thermoml2020_holdout_v5.py`; a second run only fixed a label in the generated
  bullet lines, and all numeric outputs were byte-identical). R-S1b-1 held (5/5 margins equal to the S1 yaml). Of the 19 lab systems in
  the new table, 15 were already in the S1 split because the same binaries are in `measurements_dedup.csv` through NIST TDE /
  Leenhouts (two also through the 2019 ThermoML archive); only 4 systems were test systems, 1 of them not scorable (a
  diethylbenzene component has no pure-property provider data), leaving 3 scored test systems, all in the warm bin. Every coverage reading is
  therefore "not evaluable" (n < 9). The earlier statement that 17 of 19 component sets are absent from the 2019 archive is true for
  the ThermoML archive only, not for the project's unified data.
* 2026-11-12 (second note): 738 of the new table's 1,278 mixture viscosity rows match a NIST TDE row of `measurements_dedup.csv` on system,
  composition and temperature, all within 0.5 % in value, so the overlap is the same measurements, not repeat measurements. Six of the
  15 not-independent systems are locked in the S1 split (`6b00f1bbc8c1`, `6fb947bc968c`, `b30d7680e486`, `bedb14fc8474`, `fd5cabac54de`,
  `ff0678aca0d0`); their new rows were scored and tabulated here under R-S1b-6 but no choice was made from them. All 27 rows at 253 K
  belong to three of these locked systems.
