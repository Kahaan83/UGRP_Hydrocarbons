# v6 flash arm: oxygenate mixtures (DRAFT pre-declaration, written before any table build, lock, fit or score)

Status: draft for the author's review. Nothing here has been run beyond the loader test and a row count. Date: 2026-11-21.
Model-level only; nothing is exposed in `predict_fuel` or used for certification. A negative or tied result is valid and goes to
`docs/v6_negative_results.md`.

## Question
Does adding oxygenate mixture data (1) give a mixture model that beats the strongest comparator on oxygenate blends, and (2) leave
hydrocarbon-blend accuracy unharmed? Motivation: Leenhouts et al. (Fuel 381, 2025) report ML gains over blending laws mainly for
alcohols, ketones and ethers; our unified tables hold hydrocarbon mixtures only, so no model of ours has seen oxygenate flash data.

## Data (loader `src/fuelblend/data/loaders_leenhouts_flash.py`, test `tests/test_leenhouts_flash_loader.py`, 6 passed)
Leenhouts benchmark `flashpoint_mix_exp.csv` (Zenodo, from `extras/.../datasets.zip`), 1,006 rows, all molar, unit string `Celcius`,
no T or P column (P set to 101.325 kPa, closed cup). Real-file counts from the loader (no repeated compositions, no conflicts):

| class | rows | systems | by components (systems) |
|---|---|---|---|
| hydrocarbon (`hc`) | 126 | 31 | 2: 30, 3: 1 |
| oxygenate (C, H, O only, `oxy`) | 880 | 82 | 2: 66, 3: 11, 4: 5 |
| N, S, halogen | 0 | 0 | n/a |

Staged output only (`data/interim/external_staging/`); `data/unified/` is not touched. `system_id` = md5[:12] of sorted InChIKeys
(production helper). The 31 hydrocarbon systems are the same Leenhouts rows already in the v5 data.
Pure inputs: measured pure flash points from the benchmark's `flashpoint_pure_exp.csv`, used by every method alike (as pure MON was
in Track C). Rule: rows whose measured-pure coverage (mole fraction of components with a pure flash value) is < 0.9 are excluded for
all methods and logged; they are never filled with model values and never used for any claim.

## Split (outcome-blind, fixed now)
Per chemistry class separately (`hc`, `oxy`), groups = `system_id`, ordered by md5(`v6f:` + system_id); groups are locked in that order
while locked rows stay <= 25 % of the class's rows. Remaining groups are dev. Locked sets are scored once. Dev tuning uses inner
GroupKFold(3). Expected size: about 8 hc systems and about 20 oxy systems locked, so the hc locked set is small (see gate 4).

## Arms (same frozen Track B encoder, same v1 MolPool configuration, no v2 settings, 5-seed ensemble)
- **Arm A**: train on dev hc systems only. **Arm B**: train on dev hc plus dev oxy systems.
- Candidate: frozen Track B D-MPNN embeddings, per-component projection, mole-fraction pooling, pair term, small MLP, output a residual
  over the best blending baseline (the baseline with the lowest dev OOF MAE). The encoder is not retrained and not fine-tuned.
  Caveat: the Track B pool kept oxygenates out of the flash head, so the encoder never saw an oxygenate flash label.
- Comparators, all fitted inside the same folds and arm: (a) linear molar mixing of pure flash points; (b) the Leenhouts blending law
  (flash-point blending index, mole fractions; their Eq. 2); (c) v5 ideal Le Chatelier + Raoult + Antoine, evaluated only where
  Antoine constants exist (not evaluable elsewhere, never imputed); (d) best of (a) to (c) plus an in-fold GBM on mole-weighted
  descriptors. Gate comparator = lowest dev OOF MAE, chosen on dev before the locked set is scored.

## Primary head and gate (one head; no multiplicity)
Primary claim: arm B MolPool vs the gate comparator on **oxy** systems. Pass requires all of:
1. dev out-of-fold paired group bootstrap (resample systems, 2000 draws), 90 % CI lower bound of (MAE_comparator minus MAE_candidate) > 0;
2. the same on the locked oxy systems, scored once;
3. MAE on rows with n_comp >= 3 no worse than the comparator on dev and locked (not evaluable if < 5 systems in that stratum).

Secondary claim (needed before any use on hydrocarbon fuels): **non-inferiority on hc**. 90 % bootstrap upper bound of
(MAE_armB minus MAE_armA) on locked hc systems <= 1.0 degC. Tolerance rationale: under half of the ASTM D3828 reproducibility
0.0152 x (T + 110) degC, about 2.3 degC at 40 degC (Flora et al. 2019, Table 1). If the interval is too wide to meet it, the result
is "not established", and the oxygenate model is not used for hydrocarbon fuels.

Safety reporting (trap 6, no pass criterion): signed error summary per class, mean over-prediction, 90th percentile of
(predicted minus measured). Veto: if the candidate's 90th percentile over-prediction exceeds the gate comparator's by more than
1.0 degC, the model is labelled not recommended for design use even if the gate passes.

## Reporting-only columns (never used for a pass or a skill claim)
- Benchmark's own 5-fold split definitions (as in `scripts/v5/base/delta_leenhouts_splits.py`), to compare with the paper's protocol.
- Random 10-fold over rows, to measure the gap to system-grouped scoring.
- Leave-component-out (components of a held-out system never in training), to test the stricter case.
- Per-class and per-n_comp error breakdown, as in Leenhouts Figs 5 and 6.
Each is run once, labelled "reporting-only" in every table, strict column first. Note: AGENTS.md says always GroupKFold, so the first
two need an explicit exception recorded there (not yet made).

## Expected outcome (stated now)
MolPool beats linear molar and the blending law on oxygenate blends (Leenhouts report 5.72 to 4.86 degC MAE on their own test), but a
tie with the descriptor-GBM comparator is the most likely result, as in Track C. Hydrocarbon non-inferiority is plausible but the
locked hc set is too small to confirm it with a 1.0 degC tolerance. Failure to meet either is a result, not a reason to retune.

## Files to be created after review (not yet written)
`scripts/v6/flash_oxy/build_tables.py` (stage table, coverage log, lock; writes `lock.json`, `rows.csv`, `build_log.json`),
`run_cv.py`, readings and README, each committed before the locked set is scored.

## Amendment 1 (2026-11-21): build outcome and implementation choices, written after `build_tables.py` and BEFORE any scoring
Build (`scripts/v6/flash_oxy/build_tables.py`; `rows.csv`, `lock.json`, `build_log.json`). Rows with measured-pure coverage < 0.9: 26
dropped for every method (logged in `build_log.json`); 980 kept. Five InChIKeys (stereoisomer keys) have no pure flash value; they are
never filled. Roles (locked while locked rows <= 25 % of the class, md5 order `v6f:`):

| class | dev rows / systems | locked rows / systems | n_comp >= 3 systems (dev / locked) |
|---|---|---|---|
| hc | 95 / 17 | 31 / 14 | 1 / 0 (not evaluable) |
| oxy | 641 / 57 | 213 / 21 | 8 / 7 |

The hc locked set is 14 systems with only 31 rows, so hc non-inferiority will be wide; the rule above still applies unchanged.
Implementation choices fixed now (the paper's Eq. 2 text was not available to me, so these are stated, not claimed to match it):
- **Blending-law comparator (b)** is implemented as a mole-fraction blending index, T = -theta ln(sum x_i exp(-T_i/theta)), with
  theta chosen from {10, 20, 40, 80, 160, 1e4} (1e4 = linear molar) by lowest training MAE inside each fold and arm. It is our
  implementation of a flash-point blending index, not a transcription of Leenhouts Eq. 2.
- **Comparator (c)** (ideal Le Chatelier, v5 registry constants, `lfl_temp=True`) is reported where every component has constants. It is a
  gate candidate only if it covers >= 90 % of dev oxy rows (smoke run: 0 %, so it will be reported as not evaluable on oxy).
- **Residual base**: the better of linear molar and the blending index by dev out-of-fold MAE, chosen per arm (arm B on all dev rows,
  arm A on hc rows), before the candidate is fitted. Descriptor GBM uses the same base.
- **Folds**: dev GroupKFold(5) over hc and oxy systems jointly; arm B trains on hc+oxy of the training folds, arm A on hc only. Locked:
  fit on all dev once. Descriptors, GBM, MolPool configs (v1), 5 seeds and inner GroupKFold(3) as in Track C. Locked scoring is guarded
  by an assert so it cannot be run twice.
- Reporting-only columns (benchmark splits, random 10-fold, leave-component-out) are not run in this pass; they will be added after the
  primary result and labelled reporting-only, with the AGENTS.md exception recorded then.
