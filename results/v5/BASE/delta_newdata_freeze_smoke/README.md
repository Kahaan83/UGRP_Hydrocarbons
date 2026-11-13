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

## Result 2026-11-13

Run with `scripts/v5/base/delta_newdata_freeze_smoke.py` (`PYTHONPATH=src`); files `rows_freezing.csv`, `rows_smoke.csv`, `rows_kittel_freezing_not_scorable.csv`, `summary.csv`, `readings.json`.
Plumbing note: the first run of the script scored 22 pure compounds as errors ("composition has no positive weight") because the loader stores no fractions for a pure compound; the script was fixed (pure = fraction 1.0) and the extra columns `max_group_unsafe`, in/out-of-domain constants and a post-hoc smoke bootstrap were added before any result was committed. No declared scope, rule or reading changed. "Source" in the replicate rule was applied as the paper family (ASCENT F2 and F3 are one lab and were averaged).

### Row accounting
* ThermoML: 838 liquidus rows (eutectic excluded); 629 have a finite SvL value, 209 do not (status `non_finite`: the Coutinho Tm / dHfus map has no constants for any component, mass 1.0 "liquid solvent only"; typically PAH, bicyclohexyl, decalins). Jet-cold scope: 74 rows / 10 systems declared, 37 rows / 4 systems scored, 37 rows / 6 systems unscored (bicyclohexyl + butyl/ethylcyclohexane, 1-methylnaphthalene + 1,5/2,6/2,7-dimethylnaphthalene, 1,6- + 2,6-dimethylnaphthalene; plus 2 unscored naphthalene + decalin points).
* ASCENT: 17 averaged points (10 pure, 7 mixtures); iso-octane and bicyclohexyl pure have no SvL value (2 unscored). Vozka: 22 rows, all reconstructed (unmapped mass <= 0.017) and scored. Kittel: 7 rows, not scorable (no composition), 2 censored `<`.
* Smoke: 106 averaged points; `smoke_mixtures` 68 points (56 scored by the bulk regression, 12 `low_coverage`: Mensch iso-cetane/TMB 6, MURI 3-comp A and B, Rainville S6 UM2, S7 and S8, Boehm surrogate 1 contain compounds absent from the pure-property grid), 55 by the YSI variant (gate pass 98.2 %; Eddings S2 fails the 0.90 gate).
  Excluded without recipe: Mensch JP-8 drum, Rainville Jet-A, Botero gasoline. Kittel: 4 scored, 2 censored `>` (HEFA, predictions 27.5 and 27.6 mm, below 50), FCC HN has no smoke value.

### Summary (from `summary.csv`; residual = predicted - measured)
| scope | model | rows | systems | MAE | bias | constant MAE | margin (alpha 0.10) |
|---|---|---|---|---|---|---|---|
| thermoml_jet_cold | SvL | 37 | 4 | 6.42 C | -4.83 C | 10.16 C | inf (4 systems) |
| thermoml_jet_cold | F1x (variant) | 27 | 2 | 12.41 C | -12.41 C | 4.60 C | inf |
| thermoml_all | SvL | 629 | 45 | 15.19 C | -9.19 C | 34.63 C | 76.18 C |
| ascent_mixtures | SvL | 7 | 6 | 4.02 C | +1.46 C | 5.63 C | inf (6 systems) |
| ascent_mixtures | F1x | 7 | 6 | 6.91 C | -6.42 C | 5.63 C | inf |
| ascent_pure | SvL | 8 | 8 | 0.35 C | +0.18 C | 22.0 C | inf |
| vozka | SvL | 22 | 3 | 8.01 C | +8.01 C | 3.09 C | inf (3 groups) |
| vozka | F1x | 22 | 3 | 0.95 C | -0.07 C | 3.09 C | inf |
| smoke_mixtures | bulk (designated) | 56 | 12 | 4.49 mm | -2.38 mm | 7.80 mm | 7.80 mm |
| smoke_mixtures | YSI variant | 55 | 12 | 7.56 mm | -7.56 mm | 7.83 mm | -0.32 mm |
| smoke_pure | bulk | 20 | 14 | 11.07 mm | -10.59 mm | 18.07 mm | 4.86 mm (descriptive) |
| smoke_pure | YSI | 21 | 15 | 9.87 mm | -9.87 mm | 17.60 mm | -0.02 mm |
| smoke_rainville_calib | bulk | 6 | 1 | 3.12 mm | -3.12 mm | n/a (1 system) | inf |
| smoke_kittel | bulk (measured aromatics) | 4 | 4 | 1.89 mm | -0.71 mm | 5.90 mm | inf (4 fuels) |

### R-NF-1 (freezing skill): PASS on the point estimate only
On `thermoml_jet_cold` the designated SvL has row MAE 6.42 C against 10.16 C for the leave-one-system-out constant (macro: 9.89 C against 13.43 C). The system-bootstrap interval of (constant MAE - SvL MAE) is -3.45 to +9.44 C, which includes 0.
Honest reading: the scope that was declared (10 systems, 74 rows) shrank to 4 systems / 37 rows because the model has no melting constants for six of the ten systems, and the four that remain are two n-alkane binaries (SvL group-mean offsets +1.8 and -1.5 C, individual points -4.1 to +5.6 C) and the two naphthalene + decalin systems (SvL under-predicts by 13 and 23 C on average, worst point 35.7 C: the unsafe direction). The pass is therefore a statement about four systems and does not show skill. The designated model does show skill, trivially, on pure n-alkanes (ASCENT pure MAE 0.35 C), where SvL is the melting point itself.
Other freezing scopes (descriptive): on Vozka, real jet-fuel blends with reconstructed GCxGC composition, the designated SvL over-predicts all 22 points by 6.3 to 10.8 C (MAE 8.01 C, constant 3.09 C: no skill; the over-prediction is the safe direction for a maximum limit), while the F1x variant has MAE 0.95 C and bias -0.07 C (residuals -1.7 to +1.8 C; 3 groups, blend basis unconfirmed, a and b were not fitted on these fuels). On ASCENT mixtures SvL gives one value (-47.2 C) for six solvent mixtures that differ by 12 C, MAE 4.02 vs 5.63 (6 systems). F1x on binaries is outside its intended scope and under-predicts (jet-cold scope bias -12.4 C, ASCENT mixtures -6.4 C).

### R-NF-2 (freezing margin feasibility): FAIL on the declared scope
`one_sided_margin(..., "under")`: R-NF-1 scope 4 systems -> inf (needs >= 9); ASCENT mixtures 6 systems -> inf; Vozka 3 parent-pair groups -> inf. Over all liquidus systems (45 scored) the margin is finite, 76.18 C, which is added to a predicted freezing point; it is dominated by PAH and wax binaries (worst group 321 C), has no jet-fuel meaning and is uninformative. The pooled F1x margin on all systems is 82.1 C for the same reason. Direction confirmed in `uq/conformal.py`: `"under"` scores -residual (under-prediction of freezing is unsafe against a maximum). No freezing margin is proposed; the data are insufficient (6 of the 10 jet-range ThermoML systems are unscorable by the model, and no multi-component jet-fuel freezing data with composition beyond 3 Vozka series exist).

### R-NF-3 (smoke skill): PASS for the designated model, point estimate only for the YSI variant
Designated bulk regression on `smoke_mixtures`: MAE 4.49 mm against 7.80 mm for the constant (macro 4.73 vs 7.61; 56 points, 12 systems; post-hoc system-bootstrap interval of the difference +2.26 to +4.16 mm, excludes 0). Scope mismatch: 37 of 56 points are outside the regression's domain (aromatics 8-23 vol % or result 19-25 mm). In-domain (19 points) MAE is 5.61 mm against a constant of 8.38 mm; out-of-domain (37 points) 3.92 mm against 7.51 mm. The regression is biased -2.4 mm and misses the high-smoke-point end (Botero H90T10 39 mm vs 24.3 mm, Honnet 32.9 vs 22.0, m-xylene/n-dodecane 36.9 vs 24.3). Pure compounds (all 20 scored are out of domain): MAE 11.07 vs constant 18.07 mm.
YSI variant (gate pass 98.2 %, 55 points, 12 systems): MAE 7.56 mm against 7.83 mm for the constant; the post-hoc interval of the difference is -1.77 to +3.15 mm, so it shows no skill. It under-predicts every one of the 55 points (bias -7.56 mm = MAE): a monotone offset, not noise. The pure-compound fit behind it is partly the same compounds scored in `smoke_pure`, which flatters that scope. Rainville calibration binaries (extra): bulk MAE 3.12 mm, YSI 11.65 mm, 1 system.
Kittel real fuels (bulk regression on measured aromatics, its true input): 4 scored, MAE 1.89 mm against a constant of 5.90 mm; the 3 in-domain fuels give 2.17 mm. The two censored HEFA values (>50 mm) are predicted at 27.5 mm, i.e. the regression under-predicts smoke point at zero aromatics by more than 22 mm.

### R-NF-4 (smoke margin feasibility): MET in form for both models, but no margin is recommended
`one_sided_margin(..., "over")`, 12 systems (k = ceil(13 x 0.9) = 12, so the margin is the largest group over-prediction):
* bulk regression: 7.80 mm (finite). Subtracted from a predicted smoke point it would be 7.8 mm, wider than the whole trusted range of the regression (19-25 mm, 6 mm) and a quarter of the measured span (6.1-39 mm); it is not "uninformative" by the declared span test, but it exceeds what the model can resolve. It comes from binaries and surrogates, not 100+-component jet fuels, 37 of 56 points are out of domain, and the Kittel real fuels (n = 4) give an over-prediction of at most 1.32 mm in the 4 real-fuel groups. It is proposed as an upper bound from the new data only, to be placed next to the user's decision in `configs/uncertainty_v5.yaml`; nothing was written there.
* YSI variant: the "margin" is -0.32 mm: no group is ever over-predicted because the variant is biased low by 7.6 mm. A negative margin is not a licence to add 0.3 mm to a prediction; the model has no skill (R-NF-3), so no margin is proposed.

### Proposed (not applied)
| property | model | proposed margin | status |
|---|---|---|---|
| freezing | SvL designated | none | not feasible (inf at 4/6/3 groups; 76 C on all systems is uninformative) |
| freezing | F1x variant | none | inf on every jet or mixture scope |
| smoke | bulk regression | 7.8 mm (subtract from predicted; upper bound from binaries / surrogates, 12 systems) | user decision; exceeds the model's trusted 6 mm range |
| smoke | YSI variant | none | no skill, systematic -7.6 mm bias |
