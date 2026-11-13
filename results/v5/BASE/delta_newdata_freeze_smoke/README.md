# delta_newdata_freeze_smoke

Declared 2026-11-13 before `scripts/v5/base/delta_newdata_freeze_smoke.py` exists; this commit contains no result file.
Frozen reference: commit `9ee34cc` (branch v5). User request: "Yes, score them" (score the shipped freezing and smoke models on the newly
ingested literature data and test whether a safety margin is possible). Nothing in the models, configs or yaml is changed; margins are PROPOSED only.

## Why
Freezing and smoke have no yaml margin today (freezing: no skill over a constant; smoke: the AJFTD bulk regression was fitted on jet fuels with
aromatics 8-23 vol %). The ingests of 2026-11-13 (ThermoML liquidus rows, ASCENT 065A, Vozka 2018, Kittel 2023, Mensch 2009, Rainville 2025,
Botero, Boehm 2022) are the first new measured data that none of these models saw. This delta scores the shipped models once, reports what
happens, and states whether the group-conformal margin (`fuelblend.uq.conformal.one_sided_margin`, alpha = 0.10) is even defined.

## Models scored (shipped code is called, not changed)
* Freezing designated: `predict_fuel(...)["freezing_point_C"]` (`predict_route._freezing`: ideal SvL in mole fraction, Coutinho Tm / dHfus map, max over components).
  Variant F1x (`freezing_point_C.variants.f1x_recalibrated`, `fuel_type="jet"`): reported only where it returns a value (needs an n-alkane in the input; it was
  fitted on jet-fuel S5 groups, so on binaries and ASCENT it is outside its intended scope and is labelled so).
* Smoke designated: `predict_fuel(...)["smoke_point_mm"]` = `smoke_point_from_bulk(aromatics vol %)` with `aromatics_vol_pct` computed from the recipe by the API
  (volume basis: pure-grid density at 288.15 K).
* Smoke variant YSI route: `smoke_point_mm.variants.ysi_route` (`ysi_mix` + `smoke_point_from_ysi`), only where its own-coverage gate (>= 0.90) passes; coverage reported.
* Trivial baseline in every scope: the leave-one-group-out mean of the measured values of the scope (group = `system_id`, or the parent-pair group for Vozka, or the fuel for Kittel).

## Data scopes and handling (all fixed now)
Compositions are passed as SMILES. Mole-basis recipes use `mole_fractions=`; mass-basis recipes are mass fractions; volume-basis recipes are converted to mass with the
pure-component grid density at 288.15 K (`ctx.prov.lookup`). A scored row needs a finite prediction; rows without one (coverage, non-finite) are listed with their status
and counted, never dropped silently. Censored rows (`<`, `>`) are kept out of MAE and counted separately. Rows tagged `reference` (ASTM standard values) and
`reported_mean` are not measurements of the paper and are excluded. Rows tagged `replicate` and rows sharing (source, `system_id`, composition, basis) are averaged into
one point before scoring (averaging is within a source, never across labs). `anchor_assigned` TSI values play no part because neither shipped model predicts TSI (TSI rows are
not scored; S7 YSI is an input of the YSI route, not a held-out target, so it is not scored either).

Freezing
* `thermoml`: `data/processed/v5/thermoml_sle_hc.csv`, `point_kind == liquidus` only (838 rows, 59 systems; the 11 eutectic rows are excluded). Prediction at the row's liquid
  composition (mole fractions as stored). The SvL liquidus is what is measured, but a liquidus is not an ASTM D5972 freezing point; no D5972 claim is made from these rows.
  * `thermoml_jet_cold` (R-NF-1 and R-NF-2 scope): every component has 8 <= n_carbon <= 16 (computed from SMILES; the loader flag "any jet-range compound" is not used) and row T_K <= 273.15.
    Pre-run count (composition inspection only): 74 rows, 10 systems.
  * `thermoml_all`: all liquidus rows (descriptive; also the all-systems margin of R-NF-2).
* `ascent`: `ascent065A_freezing.csv` known mixtures (binary and multi, plus the 100 % n-alkane fuel; mole basis) = `ascent_mixtures`; the 10 neat compounds = `ascent_pure` (descriptive only).
* `vozka`: 22 freezing points of Jet A, CAME, TALL, MFAT and 18 blends. Composition reconstructed like the AJFTD GCxGC fuels in BASE:
  `loaders_freezing_lit.blend_composition` / `vozka_as_ajftd_composition` -> `composition_route.kind_bins("ajftd", ...)` -> `map_bins(..., ctx.reps)` (fuel scored only if unmapped mass <= 0.10) ->
  `predict_fuel` on the mapped representatives with `composition_kind="reconstructed"`, `fuel_type="jet"`. Blend basis is UNCONFIRMED (the SI never says volume or mass; the linear mass mix is used): every blend row carries `blend_basis_unconfirmed`.
  Group = parent pair (Jet A|CAME, Jet A|TALL, Jet A|MFAT); the single neat Jet A value is scored once and assigned to Jet A|CAME.
* `kittel2023`: seven neat fuels with no composition (aromatics vol % only): NOT scorable for the freezing models (composition models); listed as such, 7 rows, 2 censored `<`.

Smoke (`data/raw/smoke_lit/`, property `smoke_point_mm`)
* R-NF-3 scope `smoke_mixtures`: known-composition hydrocarbon mixtures = Mensch (binary Appendix A replicates and Appendix B, surrogates Table 7/9), Rainville surrogates (kind surrogate; the recipe-less Jet-A real fuel is excluded), Botero binaries (H90T10 ... H30T70), Boehm binary and 6-component surrogate.
  Excluded from the scope: all pure compounds, all `real_fuel` rows (no recipe: Mensch JP-8, Rainville Jet-A, Botero gasoline), Rainville iC8/toluene calibration binaries (extra scope below).
* Extra scopes (descriptive): `smoke_pure` (pure compounds; flags `text_conflict` iso-octane and `literature_value` n-heptane kept visible); `smoke_rainville_calib` (iC8/toluene binaries; mean of manual, FURTI, VBA; `reference` excluded).
* `smoke_kittel`: the seven Kittel fuels, designated bulk regression on the measured total aromatics vol % (its actual input); the two `>` smoke values (HEFA) are censored and counted apart; FCC HN has no smoke value. The YSI route does not apply (no composition).
* In-domain flag for the bulk regression: aromatics vol % in [8, 23] and result in [19, 25] mm (its own `OutOfDomainWarning` conditions). MAE is reported in-domain and out-of-domain.
* The YSI route's pure-compound smoke-point fit (`smoke_point_pure_citedtyped.csv`) was hand-typed partly from Mensch 2010 and Boehm 2022, so pure-compound rows are partly in-sample for that variant; the mixtures are not.

## Declared readings
* R-NF-1 (freezing skill): on `thermoml_jet_cold`, designated SvL row-level MAE < constant-baseline row-level MAE (same rows, leave-one-system-out constant). Per-system macro MAE and a system-bootstrap interval (2000 resamples) of (constant MAE - SvL MAE) are reported; a PASS whose interval includes 0 is called "point estimate only".
* R-NF-2 (freezing margin feasibility): `one_sided_margin(resid, groups, alpha=0.10, direction="under")`, resid = predicted - measured, group = `system_id`, is finite on the R-NF-1 scope (needs >= 9 systems). Also reported on `thermoml_all`, `ascent_mixtures`, `vozka`. Semantics (read from `uq/conformal.py`): for freezing, a maximum limit, the unsafe side is under-prediction; `direction="under"` uses score = -residual, group score = max over its rows, margin = k-th smallest group score, k = ceil((n+1)(1-alpha)); the margin would be added to the predicted freezing point. `inf` means no valid margin at that alpha.
* R-NF-3 (smoke skill): on `smoke_mixtures` (after replicate averaging), designated bulk-regression row-level MAE < constant-baseline MAE; the same for the YSI variant on the rows where its gate passes (constant on the same rows). Gate coverage is reported.
* R-NF-4 (smoke margin feasibility): `one_sided_margin(resid, groups, alpha=0.10, direction="over")` (unsafe = over-predicting the smoke point against a minimum limit), per model, on the R-NF-3 scope, group = `system_id`, finite needs >= 9 systems. Margins are PROPOSED only: no yaml is written (`configs/uncertainty_v5.yaml` is the user's decision). A finite margin is also read against the data's own spread: if it exceeds the span between measured values of the scope it is called uninformative.

## Scope mismatch (stated now)
Binaries and surrogates are not 100+-component jet fuels, and the ThermoML systems are 2-3 component binaries (many PAH or waxy), not fuels; a margin from them is evidence about those mixtures, not a certified jet-fuel margin. The bulk regression was fitted on AJFTD fuels with aromatics 8-23 vol % and 19-25 mm and cannot respond outside; rows outside its domain are counted and reported separately.

## Outputs
`scripts/v5/base/delta_newdata_freeze_smoke.py` -> `rows_freezing.csv`, `rows_smoke.csv` (one row per scored point: source, system_id, sample, measured, predicted per model, residual, qualifiers, in_domain flags, status), `summary.csv` (scope x model: n_rows, n_systems, MAE, bias, constant MAE, margin), `readings.json`. Test: `tests/test_delta_newdata_freeze_smoke.py`.
