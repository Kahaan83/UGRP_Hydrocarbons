# delta_ws4_mixture_ml: neural mixture encoders versus the LightGBM pair residuals (declared BEFORE any fit or score)

Requirement source: roadmap `docs/superpowers/plans/2026-11-13-v5-certification-and-ml-roadmap.md` WS4; detailed plan
`docs/superpowers/plans/2026-11-14-v5-ws4-mixture-ml.md`. This is a declared comparison. Nothing here is designated, exposed in
`predict_fuel`, or certified. A candidate that passes the gate becomes *eligible* for designation, and designation is a **user decision**.

**Expected outcome (stated now): a loss or a tie is likely.** The evidence is 18 viscosity ternary systems and 21 density ternary
systems (2 locked), with 264 / 167 dev systems. Density ternaries already beat the pair model with ideal mixing (R-MLD-3 FAIL), and
Leenhouts et al. report that about 10^3 unique mixtures are needed to beat ideal mixing. A negative result is a valid result.

## Baselines (identical, read, not refitted)
- Viscosity: `ml_pair_residual` (delta_ml_viscosity). ln eta = ln eta_Refutas + sum_{i<j} x_i x_j g_ij (LightGBM).
- Density: `ml_pair_excess_volume` (delta_ml_density). V = V_ideal + sum_{i<j} x_i x_j g_ij (LightGBM).
- The baseline errors are taken from the committed `results/v5/BASE/delta_ml_viscosity/rows_scored.csv` (`resid_ml`, pred - meas, ln)
  and `results/v5/BASE/delta_ml_density/rows_scored.csv` (`err_ml`, kg/m3). Those are the dev out-of-fold and locked errors behind R-MLV-1/2/3
  and R-MLD-1/2/3. The baseline models are not refitted, and `models/v5/ml_*` is not touched.
- Rows are regenerated with the baseline scripts' own functions, loaded by path: `regenerate_rows` gives 15,827 rows and 334 systems;
  `regenerate` gives 10,900 rows and 208 systems. The script asserts positional agreement with those `rows_scored.csv` files (system_id,
  T_K to 1e-4, n_comp, role) and stops if it fails.
- Folds come from the identical call `GroupKFold(n_splits=5).split(dev_rows, groups=system_id[dev_rows])` on the same row order. The
  test-system lists are written to `folds_viscosity.json` / `folds_density.json`.
- Locked set: the same 70 viscosity / 41 density locked systems (`results/v5/S1/lab_systems_splits_v5.json`; v4 L `split_role`). The
  viscosity and density tables carry the same role for every density system (checked before fitting: 167 dev / 41 locked, 0 disagreements).

## Candidates (all residuals on the same physics value; never raw targets)
Plain PyTorch (`src/fuelblend/ml/set_models.py`), CPU, no PyG.
- Per-component features (14): n_carbon, MW, ring count, aromatic ring count, CH3 count, quaternary-C count, pure ln eta(T) from the frozen grid,
  molar volume M/rho_i(T), the 5-class one-hot, and an ln-eta-missing flag. Viscosity takes rho_i from the frozen grid, as the baseline does.
  Density takes rho_i from the tier-1 fit that defines its baseline, and ln eta from the grid where the compound has grid data, else 0 with the
  flag set. The row-level input is 1000/T.
- Standardisation uses training-fold statistics only.
1. **pds, pairwise DeepSets.** correction = sum_{i<j} x_i x_j phi(h_i + h_j, |h_i - h_j|, 1000/T), with h = MLP(features).
   This is the neural version of the pair term. It is symmetric, zero for a pure component, and defined for any n.
2. **gnn, MolSets-style.** Molecular graphs come from RDKit (atom features: C/other, degree 0-4, total H 0-3, aromatic, in ring).
   - 3 mean-aggregation message-passing layers with sum readout, run once per unique molecule per step.
   - The readout is concatenated with the pure properties at T and passed through an MLP to give z_i.
   - Fraction-weighted attention pooling: a_i = softmax_i(w . z_i + ln x_i). The MLP head output is multiplied by (1 - sum_i x_i^2), so the
     correction is zero for a pure component. That factor is a declared physics guard added to MolSets.
3. **st, set transformer, multi-task.**
   - Each token is MLP(features, x_i, 1000/T), followed by 2 SAB blocks (4 heads) and PMA with 1 seed.
   - There is one head for the viscosity residual and one for the density excess volume, each output times (1 - sum x^2).
   - Training is joint on viscosity and density rows with a masked loss. Each property's loss is system-macro and scaled by its target's
     training std, and the two losses are summed with equal weight.
   - Leakage rule: in a fold run for property P, every system in P's test fold is excluded from training in both tables.
   - The flash head named in the roadmap is NOT built. The flash lab mixtures (126 hydrocarbon rows, about 31 systems) have their own split
     protocol and are too few for a third head. This is stated now, not after results.
- **Training rows (like-for-like with the baseline):** binary dev rows with x_i x_j >= 0.01. Ternaries are never trained on, so ternary
  evidence is extrapolation for every model.
- **Targets:**
  - viscosity: r = ln meas - ln Refutas, winsorised at +/-0.5 ln;
  - density: VE = V_meas - V_ideal (cm3/mol), winsorised at +/-2.
  - Loss: MSE on the target scaled by its training std, with row weight 1/(training rows of its system), so the loss is system-macro.
  - Scoring always uses raw values. The viscosity error is resid_Refutas + correction. The density error is
    sum(x M)/(V_ideal + correction)*1000 - meas.
- **Fixed hyper-parameters (no tuning):**
  - width 64; AdamW with lr 2e-3 and weight decay 1e-4; batch 512; cosine schedule; gradient clip 1.0; torch threads 4;
  - **5-seed deep ensemble** (seeds 0-4). The prediction is the ensemble mean, and the ensemble std is kept.
  - **Declared CPU reading config: 60 epochs.** The GPU pack can run 300 epochs as a separate, labelled reading (`jobs/v5_gpu/`).
  - The CPU epoch count may be lowered only if a run proves infeasible before ANY dev number is printed, and the change is recorded here.
  - If a candidate crashes, it is reported as "not run".

## Metric
- The metric is system-macro MAE: the mean |error| over a system's rows, then the mean over systems. Viscosity is in ln and density in kg/m3.
- Intervals come from a paired system bootstrap of gain = MAE_pair - MAE_candidate: 2000 resamples, seed 42, 95 % percentile.
- The bootstrap code is the baseline scripts' own `paired_boot`.
- Viscosity masks mirror R-MLV-1/2: T >= 283.15 K. Density masks mirror R-MLD-1/2: all T.

## Declared readings (per candidate x property: pds, gnn, st for viscosity and density = 6 comparisons)
- **R-WS4-1 (dev out-of-fold):** gain interval lower bound > 0. A CI entirely below 0 is a significant loss.
- **R-WS4-2 (locked, scored ONCE by the final ensemble fitted on all dev training rows):** gain interval lower bound > 0.
- **R-WS4-3 (ternaries):** row MAE_candidate <= row MAE_pair on all ternary rows, combining dev OOF and locked from the final ensemble.
  This is the same row-MAE form as R-MLV-3 / R-MLD-3. Density is also reported against ideal mixing, because the pair baseline already
  loses there.
- **Gate (WS4):** R-WS4-1 AND R-WS4-2 AND R-WS4-3. Passing makes the candidate eligible for designation, which is a user decision.
  - There is no multiplicity correction across the 6 comparisons. A pass whose lower bound sits barely above 0 is to be read with that in mind.
- **R-WS4-4 (descriptive, viscosity cold T < 283.15 K):**
  - gain against the pair model and against Refutas;
  - the under-prediction rate, the share of cold rows with pred < meas, which is the unsafe side for cold viscosity.
- **R-WS4-5 (descriptive, uncertainty):**
  - Split-conformal margins come from dev OOF ensemble errors with the R-MLV-6 rule (per-system max |error|, alpha 0.10, per T bin, inf if
    fewer than 9 systems) for viscosity, and two-sided by group for density. Locked coverage is reported as counts.
  - The Spearman correlation of ensemble std with |error| is reported. These margins are not exposed anywhere.
- **R-WS4-6 (Leenhouts benchmark splits):** status recorded now. The benchmark's `splits.zip` / `datasets.zip` are not on disk: the
  `extras/` path used by `delta_leenhouts_splits.py` is absent, and `data/interim/leenhouts_benchmark` holds only one `compounds.csv`. The
  flash numbers on the benchmark folds already exist in `delta_leenhouts_splits`, and viscosity row matching was judged unworkable there.
  -> **not run, data absent.** Zenodo 19914760 is open; fetching it and matching viscosity rows is a follow-up, not part of this reading.

## Procedure
1. This README and the plan are committed first.
2. `--stage dev`: OOF only, and the locked systems are not scored. The output goes to `readings_dev.json`, which is committed.
3. `--stage locked`: run once, giving `readings.json` and `rows_scored.csv`. The results are appended below a separator, and this
   declaration is not edited.

## Files
src/fuelblend/ml/set_models.py, scripts/v5/base/delta_ws4_mixture_ml.py, tests/test_v5_ws4_set_models.py, jobs/v5_gpu/{README.md,run_ws4.sh},
results/v5/BASE/delta_ws4_mixture_ml/**. Checkpoints are not saved; they would go under the gitignored `runs/` if needed.
