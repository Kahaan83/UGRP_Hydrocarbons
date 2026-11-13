# delta_ml_ysi: learned group-contribution model of pure-compound YSI (gap-filler)

Declared 2026-11-13 before `scripts/v5/base/delta_ml_ysi.py` exists; this commit contains no result file and no model has been fitted on the evaluation data.
Plan source: `.superpowers/sdd/2026-11-12-v5-api-route-and-design/ml-plan-fable.md` section 2B (advisor addendum applied: the shipped model is linear).
User requirement (verbatim): "NO i need to have some ML part in the system today, that was the main point of the whole thing".
Ruling: real ML, honestly validated, shipped as a labelled gap-filler variant; nothing becomes designated or certified unless its declared readings pass.
Frozen reference: branch v5, commit 3ecf5d3. Seed 42 everywhere.

## Why
The YSI smoke route has 0 evaluable real fuels: own-YSI mass coverage on the 42 BASE smoke fuels (38 CRC AV-33-22 + 4 FuelLib; there is no AJFTD fuel in that set) is 0.01-0.56, gate 0.90.
The missing molecules are 23 mapped representatives, all saturated (n-alkanes C13-C18, single-methyl iso-alkanes C10-C17, alkylcyclohexanes C11-C18, decalin). The YSI database stops at
n-dodecane (C12) for linear alkanes and about C10 for saturated cyclics, so every fill is an extrapolation in carbon number. A tree model predicts flat beyond its last training
carbon number; the shipped gap-filler is therefore LINEAR in group counts (a learned group-contribution model). LightGBM is fitted and reported for comparison only and never fills a gap.

## Data
`data/processed/v4/sooting/ysi_hc.csv` (224 hydrocarbons; exactly what `fuelblend.baselines.sooting.ysi_table` reads; all 224 are in the Harvard Dataverse Volume 2 workbook, whose 229 hydrocarbon rows differ only by duplicates).
Canonical SMILES WITHOUT stereo (`isomericSmiles=False`); duplicates collapsed to the median. A 447-compound all-molecule model (oxygenates included) is descriptive only and does not ship.
Disclosure: the YSI -> smoke-point map `smoke_point_from_ysi` was fitted on 28 pure compounds partly from the same papers as the literature mixtures scored in R-MLY-4; it is not changed.

## Features (16, RDKit counts)
C count, H count, H/C, MW, aromatic C, aromatic ring count, aliphatic (non-aromatic) ring count, fused-ring atoms (atoms in >= 2 rings), non-aromatic sp2 C, sp C, CH3,
tertiary CH (sp3 C with 3 heavy neighbours), quaternary C (sp3 C with 4 heavy neighbours), benzylic C (sp3 C bonded to an aromatic atom), rotatable bonds, double-bond equivalents.
No folded fingerprints.

## Models
* Shipped: ridge regression on the 16 features after a StandardScaler fitted on the training fold (so it stays linear in the counts), target RAW YSI, sample weight 1/(YSI+50)^2,
  alpha in {0.1, 1, 10} chosen by inner GroupKFold(4) (inner criterion: weighted mean squared error with the same weights); predictions clipped at >= 0.
* Comparison only (never fills a gap): LightGBM n_estimators=300, learning_rate=0.05, num_leaves=7, min_child_samples=5, subsample=0.8, subsample_freq=1, colsample_bytree=0.8,
  reg_lambda=1.0, seed=42, deterministic=True, force_row_wise=True, n_jobs=1; same features, target and weight.
* Baselines: (a) global constant (training mean); (b) class mean, classes = {paraffin: acyclic saturated unbranched; iso: acyclic saturated branched; other: acyclic unsaturated; cyclo: ring(s), no aromatic ring;
  mono-aromatic: 1 aromatic ring; poly-aromatic: >= 2 aromatic rings}; a class absent from the training fold falls back to the global mean; (c) class-wise linear rule: for the paraffin, iso and cyclo classes an OLS line
  YSI ~ a + b * n_carbon fitted per class on the training fold (unweighted; fewer than 3 training points -> class mean), all other classes use the class mean.

## Split
Group key: generic Murcko scaffold (`MakeScaffoldGeneric`, canonical SMILES) for ring compounds, `acyclic_C{n_carbon}` for acyclics. Outer GroupKFold(5) (deterministic order, no shuffle); alpha selection is nested
(inner GroupKFold(4) on the outer-training groups); no hyper-parameter is chosen on a test fold. Bootstrap = resample groups with replacement, 2000 resamples, numpy default_rng(42), 95 % percentile interval.
Relative error = |pred - YSI| / YSI.

## Declared readings
* R-MLY-1 (interpolation, gate): outer-CV out-of-fold median relative error of the ridge model <= 20 %, AND ridge MAE < class-mean MAE with the group-bootstrap 95 % interval of (MAE_classmean - MAE_ridge) above 0.
  Also reported (no gate): ridge vs the global constant, the class-wise linear rule and LightGBM, each with a paired group-bootstrap interval; median relative error per class.
* R-MLY-1x (carbon-number extrapolation, gate): train on all compounds with n_carbon <= 10 (alpha by inner GroupKFold(4)); test on the union of (i) saturated compounds (no double/triple bond, no aromatic ring) with
  n_carbon 11-12 and (ii) non-aromatic ring compounds with n_carbon above the maximum n_carbon of non-aromatic ring compounds in that training set. Pass if the median relative error <= 15 %.
  Reported against the class-wise linear rule (paired bootstrap over test compounds); if the linear rule is better with the interval excluding 0 this is stated. Disclosure: this tests 1-2 carbons of
  extrapolation (n_carbon 11-12) while the fills need up to 6 (C18); n-alkane test points will be few. The test-set size is reported.
* R-MLY-2 (coverage on the 42 fuels): the 42 entities with a smoke_point_mm target in `results/v5/BASE/t2_targets.csv`; components = mapped representatives (`map_bins` over the BASE composition files, as in `cells_t2_bulk.py`).
  The final model (all 224) fills every component without measured YSI. Tier "full": a fuel is evaluable when own+ML mass coverage >= 0.90 (every fill counted). Tier "strict": counts only measured components and fills
  that are (a) within 2 carbons of the maximum n_carbon of their class (classes as above) among the 224 training compounds and (b) for ring compounds, whose generic Murcko scaffold is present among the training scaffolds.
  Pass if the full tier reaches >= 32 of 42; the strict-tier count is printed beside it (expected low). The full tier is nearly guaranteed by construction (every hydrocarbon gets a number); the strict tier and R-MLY-5 carry the information.
* R-MLY-3 (first scoring, descriptive, these fuels were already read): the BASE scoring subset (`one_per_group` of non-low-coverage fuels from `cells_t2_bulk.py`, measured smoke from the same rows): smoke MAE (mm) of the
  ML-filled YSI route vs the leave-one-group-out constant (`fuelblend.eval.group_stats.logo_constant`) vs the designated bulk regression (`predict_fuel(...)["smoke_point_mm"]["value"]` on the mapped representatives,
  `composition_kind="reconstructed"`, `fuel_type="jet"`), with signed bias. Expectation stated now: no skill (the YSI route under-predicts all 55 literature mixtures by 7.6 mm).
* R-MLY-4 (mixtures): scope `smoke_mixtures` of delta_newdata_freeze_smoke (`rows_smoke.csv`; recipes rebuilt from `loaders_smoke_lit` with that script's basis handling). Report: points evaluable without ML fill (existing gate)
  vs with ML fill; MAE of the YSI route with/without (without on its own subset; with on all and on the same subset); the designated bulk MAE on the same rows; and the split-conformal 90 % interval on YSI from OOF errors:
  score |y - yhat|/(yhat + 50), q = the ceil((n+1)*0.9)-th smallest of the n OOF scores; interval yhat -/+ q*(yhat+50). For `ml_predicted_extrapolated` fills the interval uses max(q, q_x) where q_x is the same quantile of the
  R-MLY-1x test scores (descriptive; calibration for extrapolation is NOT claimed). Coverage of the 224 measured values (OOF) is reported as counts.
* R-MLY-5 (sensitivity): for each of the 42 fuels recompute the ML-filled YSI-route smoke point with every ML-filled component YSI multiplied by 0.7 and by 1.3; report the largest absolute change (mm) per fuel. If <= 1 mm for >= 90 % of fuels,
  the fills are immaterial to the smoke conclusion (that is the claim made); otherwise the route stays "not evaluable on real fuels".

## Adoption rule
Both gates (R-MLY-1 and R-MLY-1x) must pass -> a new variant `ysi_route_ml_filled` is offered; the existing `ysi_route` variant and the `HAS_YSI` predicate are unchanged; no margin and no certification; the designated smoke model
stays the bulk regression (user decision). If R-MLY-1x fails, the module and results stay as a scored negative and the API exposes nothing ("sooting index for heavy saturates cannot be extrapolated from this database; the sooting route stays unscored on real fuels").
Source labels per component: `measured` | `ml_predicted` (strict tier) | `ml_predicted_extrapolated` (beyond it). Never overwrites a measured value.

## Wiring plan (for the controller; I do not edit predict_route.py)
Call site: `src/fuelblend/api/predict_route.py::_smoke_ysi_variant`, after the existing `ysi_route` entry. Call `fuelblend.ml.ysi_model.predict_ysi([smiles...])` (DataFrame columns `smiles, ysi, lo, hi, source`) for the components, then
`prop["variants"]["ysi_route_ml_filled"] = {value, status, ysi_measured_mass, ysi_ml_mass, ysi_ml_extrapolated_mass, ysi_source_by_component, ysi_interval, method: "ysi_mix+ridge_group_contribution"}`;
`status = "ok_extrapolated_inputs"` whenever `ysi_ml_extrapolated_mass > 0.10`, else `"ok"`. Emit only when the adoption rule is met (`readings.json: adoption`).

## Files
`src/fuelblend/ml/ysi_model.py`, `scripts/v5/base/delta_ml_ysi.py`, `models/v5/ml_ysi/{model.json,meta.json}`, `results/v5/BASE/delta_ml_ysi/*`, `tests/test_ml_ysi.py`. Report: `.superpowers/sdd/2026-11-12-v5-api-route-and-design/ml-ysi-report.md`.
