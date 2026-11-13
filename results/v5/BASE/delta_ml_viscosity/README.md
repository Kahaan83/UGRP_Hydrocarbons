# delta_ml_viscosity: learned pair-interaction residual on the Refutas mixture viscosity (declared BEFORE any fit or score)

Requirement source: user ("NO i need to have some ML part in the system today"), controller ruling, Fable-reviewed plan section 2A
(`.superpowers/sdd/2026-11-12-v5-api-route-and-design/ml-plan-fable.md`). Real ML, honestly validated. It ships only as a labelled
variant; nothing becomes designated or certified unless the readings below pass. Promotion to designated is a user decision.

## Why
AGENTS.md 12.1: property = physics baseline + learned correction. Earlier real-fuel residual attempts (16-35 fuels) did not transfer
(`delta_t2_residual`). The lab mixture viscosity table is far larger (15,827 refutas rows, 334 systems), enough to learn a pair-interaction
(Grunberg-Nissan style) correction between the pair of components.

## Model
ln eta_ML = ln eta_Refutas + sum_{i<j} x_i x_j g(d_i, d_j, T), x = mole fractions, g = LightGBM on pair-symmetric features.
No composition feature inside g: the correction is a quadratic form in x, zero at the pure components, permutation-invariant, and defined
for any number of components. Scoring is on the lab target (dynamic viscosity, mPa s, ln); because the Refutas density is the same ideal-volume
density in numerator and denominator, the additive ln correction is identical on kinematic viscosity.

## Data (all read-only)
- Rows: `results/v5/S1/viscosity_lab_residuals_v5.csv`, `method == "refutas"` only (15,827 rows, 334 systems). `resid` = ln(pred/meas)
  so the target is r = -resid = ln(meas) - ln(pred).
- That table has no compositions. The rows are regenerated with x_mole and InChIKeys by loading `scripts/v5/eval_lab_known_composition_v5.py`
  by path (`lab_selection`, same measurement concatenation, same provider). Before any fit the script asserts positional agreement with the
  committed predicted/measured to 1e-6 and stops if it fails; it asserts the pure-grid SHA-256 =
  b931e969926698dff52d69ff51dd99a990c24ff12eca545aa98e970afd3e7daf.
- Roles: `results/v5/S1/lab_systems_splits_v5.json` (264 dev / 70 locked systems; group key `system_id`).

## Split and fit
- Train: binary rows of dev systems with x1*x2 >= 0.01 (rows with smaller x1*x2 are scored, never trained on). Ternary rows are never trained on.
- Dev estimate: GroupKFold(5) on dev system_id (out-of-fold). Dev ternary systems are scored with the model of the fold that held them out.
- Test: the 70 locked systems, scored ONCE by the model fitted on all dev training rows, after this README + the hyperparameters are committed.
  Disclosure: raw Refutas was already scored on the locked systems (S1); no fitted model was. A locked ternary may contain a dev binary pair (intended).
- Target for g: g = r/(x1 x2), with r winsorised at +/-0.5 ln before the division (declared reading of "training targets winsorised at +/-0.5 ln");
  sample weight (x1 x2)^2 normalised to sum 1 per system. Scoring always uses raw values.
- Hyperparameters (fixed, no tuning): LightGBM n_estimators=400, learning_rate=0.05, num_leaves=15, min_child_samples=40, subsample=0.8,
  subsample_freq=1, colsample_bytree=0.8, reg_lambda=1.0, seed=42, deterministic=True, force_row_wise=True, n_jobs=1.
- Features (pair-symmetric mean and |difference| of each): n_carbon, MW, ring count, aromatic ring count, CH3 count, quaternary-C count,
  pure ln eta at T (frozen grid), pure molar volume at T (MW/rho from the same grid). Plus: unordered class pair (categorical, classes from
  `compound_meta`), 1000/T.
- Baselines: raw Refutas (g = 0, the gate); one global constant g (weighted mean on training rows); class-pair mean table of g (reported).

## Metric
System-macro MAE of |ln error| (mean over a system's rows, then mean over systems). Intervals: paired system bootstrap, 2000 resamples,
seed 42, 95 % percentile. "T >= 283 K" means T >= 283.15 K (bin T_ge_283K); cold = T < 283.15 K (bins T_lt_263K and T_263_283K).

## Declared readings
- R-MLV-1 (dev out-of-fold, T >= 283 K): MAE_ML < MAE_Refutas and the interval of (MAE_Refutas - MAE_ML) is above 0.
- R-MLV-2 (locked 70 systems, one run, T >= 283 K): same rule, and at most 25 % of systems worse by > 0.02 ln (system MAE_ML - MAE_Refutas > 0.02).
- R-MLV-3 (ternaries): on all ternary rows MAE_ML <= MAE_Refutas (dev out-of-fold + locked with the final model); descriptive split by whether all
  three unordered InChIKey pairs occur in the training binaries used for that model.
- R-MLV-4 (cold gate, T < 283.15 K): dev out-of-fold + locked systems with cold rows: MAE_ML < MAE_Refutas with the interval above 0 and the
  under-prediction rate (share of cold rows with pred < meas) not raised. If it fails the correction is defined only for T >= 283 K.
  Disclosure: all cold rows are Boehm 2022, so a cold pass is about one source.
- R-MLV-5 (real-fuel transfer, descriptive, already read by BASE): T2 rows of `results/v5/BASE/rows_t2_bulk.csv` (ajftd_gcxgc, av3322, fuellib),
  correction as the quadratic form over mapped representatives; MAE per temperature vs Refutas. Expected weak.
- R-MLV-6 (margin): conformal margin per T bin from dev out-of-fold residuals, same rule as the yaml (per-system max |ln ratio| score,
  k = ceil((n+1)(1-alpha)), alpha 0.10, two-sided); a bin with < 9 dev systems stays inf. Locked coverage reported as counts only.

## Adoption rule (declared)
R-MLV-1 and R-MLV-2 pass -> ships as variant `ml_pair_residual` with the R-MLV-6 margin, in domain only (every component has its own grid data,
every class pair has >= 5 training systems, T inside the R-MLV-4 gate, known composition; reconstructed real fuels get the value with
in_domain False and no margin). Designated stays raw Refutas; promotion is offered to the user only if R-MLV-3 passes and the ML margin for
T >= 283 K is <= 0.186 ln. If a gate fails the module and results stay as a scored negative and the API exposes nothing.

## Leakage traps (known before fitting)
Multiple method rows per point (filtered to refutas); residual sign (target = -resid); endpoint mismatch (paper pure values vs frozen grid give a
residual linear in x that x_i x_j g cannot represent: report the share of systems with endpoint rows and their endpoint residual; the target is
not corrected); cold rows are all Boehm; no tuning after seeing locked results.

## Wiring plan (for the controller; I do not edit predict_route.py)
In `predict_route.predict_fuel`, after the Refutas value, call
`fuelblend.ml.viscosity_pair.ln_correction(smiles_list, x_mole, T_K)` and set
`res["viscosity_kin_cSt"]["variants"]["ml_pair_residual"] = {value, ln_correction, status, in_domain, margin, margin_unit: "ln",
method: "refutas+ml_pair_gn", training_scope: "lab binaries, 264 dev systems"}`; dynamic viscosity mirrors it; not stacked with K1.
Design: `design_forward_check` gains a column `ml_variant_pred`.

## Files
src/fuelblend/ml/{__init__,viscosity_pair}.py, scripts/v5/base/delta_ml_viscosity.py, models/v5/ml_viscosity_pair/{model.txt,meta.json},
results/v5/BASE/delta_ml_viscosity/**, tests/test_ml_viscosity_pair.py.

---

## Outcome (appended AFTER scoring; the declaration above is unchanged)

Reproduce: `.venv/Scripts/python.exe scripts/v5/base/delta_ml_viscosity.py` (full run; `--dev-only` skips the locked systems),
then `scripts/v5/base/delta_ml_viscosity_realfuels.py` (R-MLV-5). Numbers: `readings.json`, `readings_R-MLV-5.json`, `rows_scored.csv`, `run_log.txt`.
Positional regeneration check against the committed S1 table passed (15,827 rows, 334 systems). The locked systems were scored once.
MAE below = system-macro mean |ln error|, Refutas -> ML; interval = paired system bootstrap of the gain.

| reading | result | numbers |
|---|---|---|
| R-MLV-1 dev OOF, T >= 283 K (264 systems) | PASS | 0.0824 -> 0.0599, gain interval [0.0180, 0.0269] (global-constant g 0.0692, class-pair table 0.0649) |
| R-MLV-2 locked, T >= 283 K (70 systems) | PASS | 0.0544 -> 0.0311, gain [0.0163, 0.0302]; 1.4 % of systems worse by > 0.02 ln (const 0.0431, class-pair 0.0390) |
| R-MLV-3 ternaries (664 rows, 18 systems) | PASS | row MAE 0.0841 -> 0.0398; all-pairs-in-training 5 systems 0.0834 -> 0.0410; some pair unseen 13 systems 0.0844 -> 0.0393 |
| R-MLV-4 cold T < 283.15 K (26 systems, Boehm 2022 + thermoml) | PASS | 0.0994 -> 0.0612, gain [0.0225, 0.0547]; under-prediction rate 0.631 -> 0.571 |
| R-MLV-5 real fuels (descriptive) | MIXED | 394 rows; Refutas reproduced to 4.5e-8; MAE 0.160 -> 0.086 ln pooled, but ML is worse at warm cells (e.g. av3322 273-313 K, ajftd 313 K, fuellib 313 K): not worse in 13 of 21 source x T cells; large gain only at 233/253 K |
| R-MLV-6 margin (ln, alpha 0.10) | computed | ML: T < 263 K 0.287 (20 dev systems), 263-283 K 0.608 (16), >= 283 K 0.134 (264); Refutas 0.415 / 0.608 / 0.185. Locked coverage at >= 283 K: ML 64/70, Refutas 66/70 |

Adoption: R-MLV-1 and R-MLV-2 pass, so the variant `ml_pair_residual` ships with the margin above, in domain only. R-MLV-3 passes and the ML
>= 283 K margin 0.134 <= 0.186, so promotion to designated is OFFERED to the user; it is not made here. Because R-MLV-3 rests on 18 ternary
systems and R-MLV-5 is mixed on 100-component fuels, no claim beyond lab mixtures is supported.

Caveats (measured): 35 % of binary systems (112) carry endpoint rows with x1 x2 < 0.01 (mean endpoint residual +0.012 ln, mean |.| 0.034 ln),
which the pair form cannot represent. Training systems per class pair are < 5 for alkene_alkyne|{cycloalkane, n_paraffin} (1 each) and
iso_paraffin|iso_paraffin (2); such mixtures are out of domain. Cold evidence is essentially one source (Boehm 2022). The cold gate passed, so the
correction is defined at all T in the training range (233-393 K). Locked ternaries may contain dev binary pairs (intended).
