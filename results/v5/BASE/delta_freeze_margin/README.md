# delta_freeze_margin

Declared 2026-11-13 before `scripts/v5/base/delta_freeze_margin.py` exists; this commit contains no result file.
Frozen reference: branch v5 at commit `d3e13d1` (plus the verified constants of commit `3ecf5d3`). User request: "make sure that all of them get
certified, however possible". Controller ruling: by honest means only (no gate loosening, no post-hoc scope cut). The declared text below is never edited
after results; results are appended under a "Result" heading.

## Why
`results/v5/BASE/delta_newdata_freeze_smoke/` could not give a freezing margin: the jet-cold ThermoML scope shrank from 10 to 4 systems because the
Coutinho Tm / dHfus map has no constants for six of the ten systems. `data/raw/sle_constants/tm_dhfus_verified.csv` (39 compounds, primary-source Tm and
dHfus, SOURCES.md) now supplies constants for those compounds. This delta re-scores the same kind of liquidus data with those constants and tests whether
a finite freezing margin exists.

## Disclosure: the data were scored once already
The ThermoML liquidus set (`data/processed/v5/thermoml_sle_hc.csv`, `point_kind == liquidus`, eutectic rows excluded) was already scored in
`delta_newdata_freeze_smoke` with the old Coutinho constants (37 rows / 4 systems in its jet-cold scope; 629 finite rows / 45 systems overall). The
verified constants were selected AFTER seeing that scoring (the 31 list-B compounds are exactly the components of the unscored systems). The scope S_F below
is therefore not a fresh holdout; it is the same data re-scored with the missing constants filled in. The constants themselves come from primary
literature (NIST WebBook calorimetry rows), not from fitting to these liquidus points; but their selection was informed by which systems failed.

## Declared choices
* **Constants, `tmap_verified`.** A map InChIKey -> {tm_K, dh_fus_J_mol} built from `tm_dhfus_verified.csv`: `tm_K` = column `tm_K`, `dh_fus_J_mol` =
  1000 x `svl_dhfus_kJ_mol` (the "total incl. solid-solid" convention for odd n-alkanes, as recommended in SOURCES.md). A row with a blank Tm or dHfus
  (stereo-unspecified decahydronaphthalene, 1-hexyne) has no entry. For SCORING a component without an entry has no constant and does not crystallise in the SvL
  sum (the same treatment as a compound absent from the Coutinho table today). Coutinho values are NOT mixed into scoring. Builder:
  `fuelblend.eval.composition_route.build_tmap_verified` (new function; `get_pure_thermal_properties` is frozen and not edited).
* **Known defect of the old table (not fixed; follow-up).** `get_pure_thermal_properties` identifies aromatics by SMILES substring: `"c1ccccc1"` matches
  fluorene and phenanthrene (both receive benzene's 278.7 K / 9.87 kJ/mol; true 387.94 K / 19.58 and 372.38 K / 16.46), and `"c1ccc2ccccc2c1"` matches
  2-methylnaphthalene (receives naphthalene's 353.4 K / 19.0; true 307.73 K / 12.13). The verified map bypasses it. The old-constant numbers reported below
  still contain this defect.
* **Scope S_F (fixed now).** ThermoML liquidus rows (eutectic excluded) such that: every component is a hydrocarbon with 8 <= n_carbon <= 16 (carbon atoms
  counted from the SMILES), every component has an entry in `tmap_verified`, and row T <= 273.15 K. PAH and naphthalene systems are NOT removed.
  group = `system_id`. Prediction: ideal Schroeder-van Laar in mole fraction at the stored liquid composition (`ideal_sle_C`: maximum over components with
  x > 1e-4), the function the API and designs use. A row with a non-finite prediction is listed with its status and counted, never dropped silently.
* **Constant comparator.** Leave-one-system-out mean of the measured liquidus temperatures (degC) over S_F, per row; row-level MAE.
* **Margin.** `fuelblend.uq.conformal.one_sided_margin(resid, groups, alpha=0.10, direction="under")`, resid = predicted - measured (degC), group = `system_id`
  (under-prediction of a freezing point is the unsafe side against a maximum limit; the margin is added to the prediction, `safe_bound = hi = value + margin`).

## Declared readings
* **R-FM-1**: designated SvL with `tmap_verified` on S_F has n_systems >= 9 (systems with at least one scored row).
* **R-FM-2**: MAE on S_F (SvL) < leave-one-system-out constant MAE (same rows). A system-bootstrap interval (2000 resamples) of (constant MAE - SvL MAE) is
  reported; a pass whose interval includes 0 is called "point estimate only". The decision uses the point estimate, as declared in `delta_newdata_freeze_smoke`.
* **R-FM-3**: the margin above is finite on S_F.
* **Rule (declared now).** If R-FM-1, R-FM-2 and R-FM-3 all pass, that margin is written into `configs/uncertainty_v5.yaml` as
  `freezing_point_C.lab_liquidus.svl_verified` (direction `under`, unit C, alpha 0.1, n_groups, n_rows, residual file, scope text, caveats) and used by the API
  and the designs. If any fails, NOTHING is written to the yaml, the API and designs are left as they are, and the result says "not written".
* **R-FM-4** (only if a margin is written): the four design runs of `results/v5/BASE/delta_smoke_margin/designs_domain_fix/` (`v5pool`, `ext_n30`, `ext_n8`,
  `guard4_v5pool`; seed 42, `scripts/v5/run_design_userpool_v5.py`, the same arguments) are re-run with freezing judged as `pred + margin <= limit` and the
  members' freezing constants marked verified (source csv) instead of `coutinho_table`. Per design, all constraint statuses and the category (A / B / certified)
  are reported. Stated expectation: freezing may certify; smoke cannot (the in-domain regression result is limited to 19-25 mm, so it can never reach the
  25 mm limit plus the 7.8 mm margin; and some smoke limits are unverified), so 0 designs are expected to be certified overall unless smoke changes.
  A design is `C_certified` only if the engine passes all constraints and the forward check meets every limit (existing definition in the run script).

## Implementation choices (declared now; applied only if the rule writes the margin)
* **API** (`predict_route._freezing`): the designated SvL value switches to constants "verified where present, else the Coutinho entry (flagged non-verified),
  else none". The old value (all-Coutinho constants) stays as `variants["coutinho_constants"]`. The margin (`hi = value + margin`, `safe_bound = hi`) is applied
  only when every mapped hydrocarbon component has a verified constant and nothing was substituted or unmapped (`_finish` rule); otherwise the property is
  `uncertified` as today. (This is stricter than "every crystallising component": a component without a verified constant cannot be shown not to crystallise.)
  If the rule fails, `_freezing` is not changed at all.
* **Designs** (`library_route`, freezing part only): same constants; members with verified constants are not recorded as non-verified, the yaml key is attached to
  the freezing model (the engine then applies the margin on the maximum limit, as for flash). Members without a verified constant fall back to the Coutinho entry
  and are recorded `coutinho_table` (freezing then stays uncertified for that design), as today. `engine_v4.py` is not edited. Viscosity code is not touched.
* `scripts/v5/base/build_uncertainty_v5_q2.py` gains the freezing entry (asserts: flash margin 2.206932808388615 and smoke margin 7.8009965628100595 unchanged).

## Descriptive extra reports (no decision attached)
* Vozka 2018 real-fuel freezing points (22 rows, reconstructed GCxGC composition, same route as in `delta_newdata_freeze_smoke`) and ASCENT 065A mixtures, with
  `tmap_verified` for the SvL constants (verified-only as in scoring), next to the old-constant numbers of `delta_newdata_freeze_smoke/summary.csv`.
* The old-constant numbers on the same S_F rows (scored with the Coutinho map) for comparison.

## Outputs
`scripts/v5/base/delta_freeze_margin.py` -> `rows.csv` (one row per scored point: system_id, sample, T_K, measured, pred_verified, pred_old, residual, status,
scope flags), `summary.csv` (scope x constants: n_rows, n_systems, MAE, bias, constant MAE, margin), `readings.json`. Tests: `tests/test_delta_freeze_margin.py`.
