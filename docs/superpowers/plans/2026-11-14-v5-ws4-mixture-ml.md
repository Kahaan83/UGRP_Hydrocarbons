# v5 WS4: Mixture ML (pairwise DeepSets, MolSets-style GNN, set transformer) versus the LightGBM pair residual

> For agentic workers: executed inline by the WS4 agent (no subagents). Subagents and reviewers: Opus at most, never Fable.
> Roadmap source: `docs/superpowers/plans/2026-11-13-v5-certification-and-ml-roadmap.md`, section WS4. Common rules: WS2-7 rules file (side-session notes).

**Goal:** a declared, like-for-like comparison of three neural mixture encoders with the two shipped LightGBM pair models,
`ml_pair_residual` (viscosity, delta_ml_viscosity) and `ml_pair_excess_volume` (density, delta_ml_density). Every model is a residual on the
same physics value. The comparison uses the same rows, the same GroupKFold(5) folds on dev `system_id`, the same locked systems, and the
same metric. **Expected outcome, stated before any fit: a loss or a tie is likely.** The reasons are 18 viscosity ternary systems, 21 density
ternary systems (2 locked), 264/167 dev systems, and Leenhouts' finding that about 10^3 unique mixtures are needed to beat ideal mixing.
The value of the work is the evidence either way.

**Not in scope (user decisions):** designating or exposing any variant in `predict_fuel` / `predict_route.py`; touching `models/v5/ml_*`;
pushing. A candidate that passes its gate only becomes *eligible* for designation, and designation is the user's decision.

**Architecture:** one script `scripts/v5/base/delta_ws4_mixture_ml.py` plus one module `src/fuelblend/ml/set_models.py` (plain PyTorch, no PyG).
- The script loads `delta_ml_viscosity.py` / `delta_ml_density.py` by path and reuses their row regeneration (with positional asserts
  against the committed S1 / v4-L tables), `train_set`, `macro`, `paired_boot` and `conformal`.
- The baseline errors are READ from the committed `rows_scored.csv` (`resid_ml`, `err_ml`). They are never refitted, so the committed
  baseline models and results are not touched.
- The regenerated rows are asserted to align positionally with those files.
- The fold test-system lists are recomputed with the identical `GroupKFold(5).split(dev_rows, groups=...)` call and written to the delta folder.

**Tech stack:** Python 3.12, torch 2.15 CPU (already in `.venv`), RDKit, numpy/pandas, sklearn GroupKFold, pytest.

## Tasks

### Task 0: Declare (no code, no fit)
- [ ] Write this plan and `results/v5/BASE/delta_ws4_mixture_ml/README.md`: data, models, fixed hyper-parameters, CPU reading config, gates, expected outcome, Leenhouts status.
- [ ] Commit both, before any training code runs.

### Task 1: Module `src/fuelblend/ml/set_models.py` (TDD)
- [ ] Tests first, `tests/test_v5_ws4_set_models.py`:
  1. permutation invariance of all three models (shuffle component order, same output to 1e-6);
  2. zero correction for a pure component (one component with x = 1, or a padded set with one live component);
  3. padding invariance (an extra masked slot does not change the output);
  4. the deterministic seed gives the same output twice;
  5. the graph featuriser yields the right atom count and edges for propane and toluene.
- [ ] Implement:
  - a `PairDeepSets` encoder h = MLP(features); correction = sum_{i<j} x_i x_j phi(h_i + h_j, |h_i - h_j|, 1000/T);
  - a `MolSetsGNN` (plain-torch mean-aggregation message passing over unique molecules, sum readout, concatenated with pure properties at T; fraction-weighted attention pooling, logits + log x_i; output times (1 - sum x_i^2));
  - a `SetTransformer` (2 SAB blocks, PMA with 1 seed, per-property heads, output times (1 - sum x_i^2));
  - a `train_ensemble(...)` helper: AdamW, cosine schedule, fixed epochs, system-macro weighted MSE, 5 seeds, mean prediction and std.

### Task 2: Script `scripts/v5/base/delta_ws4_mixture_ml.py` (TDD for the data path)
- [ ] Tests (fast, synthetic fixtures; no real data): fold-list reproduction from a toy frame equals `GroupKFold` on the same call; the alignment assert fails on a shifted frame; the gate function returns pass only when lo > 0.
- [ ] Data path: viscosity rows from `delta_ml_viscosity.regenerate_rows` (15,827 rows, 334 systems); density rows from `delta_ml_density.regenerate` (10,900 rows, 208 systems).
  - Assert alignment with the committed `rows_scored.csv` (system_id, T_K, n_comp, role).
  - Write `folds_viscosity.json` / `folds_density.json`.
- [ ] Per-component features: static descriptors, class one-hot, pure ln eta(T) (viscosity), molar volume M/rho_i(T), 1000/T.
  - Standardise with training-fold statistics only.
  - Graphs from the registry SMILES (viscosity) or from the InChI (density).
- [ ] `--stage dev` produces OOF readings only (locked untouched) -> `readings_dev.json`. `--stage locked` fits the final 5-seed ensemble on all dev training rows, scores the locked systems once -> `readings.json`, `rows_scored.csv`.
- [ ] `--epochs/--seeds/--width/--device` flags exist for the GPU pack. The declared CPU reading uses the defaults in the README.

### Task 3: Run (background), dev first, then locked once
- [ ] Run `--stage dev` for all candidates in the background. Commit `readings_dev.json` and the run log.
- [ ] Run `--stage locked` once. Append the results to the README below a separator; the declaration stays unchanged.

### Task 4: GPU pack `jobs/v5_gpu/`
- [ ] `README.md` (inputs incl. gitignored data, rules, explicit "not run on a GPU", what a GPU run adds and does not add) and `run_ws4.sh [cuda|cpu]`.

### Task 5: Records
- [ ] `docs/claims_register.md`: provisional claim ids `CLM-WS4-*`, which the controller renumbers after CLM-49 when the WS branches merge.
- [ ] `docs/model_decision.md` "v5 ML components": one bullet with the outcome. No designation.
- [ ] One full pytest at the end (baseline: 1 known failure `tests/test_v4_sle.py::test_default_mode_unchanged`).
- [ ] Commit on the worktree branch (no push, no attribution lines).

## Review focus
1. Baseline identity: the candidate and the baseline must be scored on the same rows, folds and locked set. The asserts are the evidence.
2. Leakage in the multi-task set transformer: a system that is in the test fold of the scored property must not enter training through the other head. The roles of the viscosity and density tables agree (checked: 167 dev / 41 locked density systems carry the same role in `lab_systems_splits_v5.json`). The per-property fold runs exclude that fold's test systems from both tables.
3. Gate direction: gain = MAE_pair - MAE_candidate. A CI that lies entirely below 0 is a significant loss, not a pass.
4. No tuning after locked results: the hyper-parameters are fixed in code and in the README.
