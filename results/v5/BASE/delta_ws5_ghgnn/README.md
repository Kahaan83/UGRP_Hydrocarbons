# delta_ws5_ghgnn: does GH-GNN activity-coefficient non-ideality improve Le Chatelier flash? (a test, not a fix)

Plan: `docs/superpowers/plans/2026-11-14-v5-ws5-ghgnn.md`. Roadmap source: WS5 of `docs/superpowers/plans/2026-11-13-v5-certification-and-ml-roadmap.md` (not in this worktree; read from the main checkout).

Why: original UNIFAC made flash worse (v4 F: dMAE +0.45 C). GH-GNN (Sanchez Medina et al., Digital Discovery 2023, 2, 781-798, doi 10.1039/D2DD00142J) predicts ln gamma_inf(T) for a solute in a solvent from the two molecular graphs, and was trained on measured infinite-dilution data (Brouwer et al. 2021 set). This delta tests whether using it in place of gamma = 1 helps. **No live model is changed.** The live route stays ideal Le Chatelier with the LFL overlay.

Code was cloned (MIT licence, commit 0bbd2b3, `data/raw/ghgnn_code/`, see `SOURCES.md` there), read before running, and reproduces the repo's own stored predictions on hydrocarbon rows (`validation_repo_predictions.json`: 400 rows, max |diff| 2.4e-4 in ln gamma, mean 4.3e-5; the repo's stored CSV values are rounded). GH-GNN runs, so the UNIFAC-Dortmund fallback is NOT used.

## Declared before any flash score was computed

Written and committed before `ws5_flash.py` produced a single prediction. Nothing below is changed after results.

### Model (primary, one only)

* **GH-GNN-Wilson.** For every ordered pair (solute i in solvent j) the network gives ln gamma_ij_inf(T) = A_ij + B_ij / T_K exactly (the architecture's last line). A and B come from two network evaluations (25 and 75 C), and a third (50 C) must reproduce them to 1e-4 (unit test). For a mixture, each unordered pair {i, j} is inverted to Wilson parameters (Lambda_ij, Lambda_ji) at the T the solver is asking about, from ln gamma_i_inf = 1 - ln Lambda_ij - Lambda_ji and ln gamma_j_inf = 1 - ln Lambda_ji - Lambda_ij. The multicomponent Wilson expression then gives gamma_i(T, x) with Lambda_ii = 1. gamma_fn is handed to `predict_flash_le_chatelier` (same solver, Antoine, LFL, `lfl_temp=True`, `lfl_eff` scaling as `flash_point_unifac_detail`). `flash_unifac.py` is not modified.
* **Wilson inversion failure rule (fixed now):** if a pair has no solution (residual > 1e-8 after the solver and a fsolve restart), that pair gets Lambda_ij = Lambda_ji = 1 (no interaction). The number of such pairs is counted and reported per set.
* **Two-parameter Margules is a secondary, binary-only sensitivity** (lab binaries; A12 = ln gamma_1_inf, A21 = ln gamma_2_inf at T). It has no multicomponent form here, so it cannot carry CRC or T2 results and is not part of the gate.
* **Applicability domain (AD).** Repo definition: a pair is in domain if both molecules were seen in training (solute as solute, solvent as solvent), or if the 10th-highest RDKFingerprint Tanimoto >= 0.35 AND the ClassyFire class-pair count >= 25. The class names come from the repo's own `molecular_compounds.csv` (InChIKey first block), otherwise from an explicit hydrocarbon rule (flagged `class_source`); no web call. For EVERY ordered pair used, the table stores the real Tanimoto indicator, the real class-pair count, the interpolation flag and the verdict (not 9999). Per mixture: pair-weight in domain = sum over i != j of x_i x_j * [pair in domain in both directions] / sum over i != j of x_i x_j.
* **Primary variant = all pairs as predicted by the network (gated).** Secondary variant (reported, not gated) = out-of-domain pairs set to Lambda = 1.

### Baseline and test sets (all re-predicted, nothing taken from a table unchecked)

Baseline = the live ideal Le Chatelier route with the LFL overlay registry (`RouteCtx(lfl_override=...)`), gamma = 1. The first test of the code is that gamma = 1 through the new wrapper equals the stored values: lab = `err_lfl_overlay` in `delta_lab_flash_li2014/points.csv`; CRC 647 r3 = `ideal_C` in `results/v5/L1/F_r3/flash_real_fuels.csv`; T2 = `pred_lfl_overlay` in `results/v5/BASE/rows_t2_flash.csv`. A mismatch above 0.01 C stops the run.

| id | set | unit of analysis (cluster) | role in the gate |
|---|---|---|---|
| S-LAB | lab mixtures, `results/v4/F/flash_lab_loso.csv`, 126 rows, 31 systems | `system_id` | gated |
| S-CRC | CRC 647 r3, 57 fuels with `statistics_allowed`, split_role != locked | `fuel_id` | gated |
| S-CRCL | the 12 locked CRC 647 fuels | `fuel_id` | read once, descriptive; if it shows a rise in the unsafe rate it is reported as contradicting the gate |
| S-T2 | T2 measured flash (chosen, not low-coverage, status ok), pooled over av3322, fuellib, ajftd_gcxgc | T2 `group` | gated |
| S-LI | 33 Li 2014 / Ha-Lee 2016 n-alkane points, 6 systems; with S-LAB gives 36 systems | `system_id` | margin reading only (R-4) |

### Readings and gate

* **R-1 gate, per gated set (S-LAB, S-CRC, S-T2).** Paired cluster bootstrap (resample clusters with replacement, B = 5000, seed 20261109, the v4 values) of dMAE = MAE(GH-GNN-Wilson) - MAE(ideal). PASS on a set needs BOTH (a) the 95 % percentile CI of dMAE entirely below 0 AND (b) the over-prediction rate (error = predicted - measured > 5 C) of GH-GNN <= that of ideal, row level, no tolerance. **The gate passes overall only if all three gated sets pass.** Anything else is reported as "gate not met" with the sets that failed. Over-prediction is the unsafe direction.
* **R-2 reported with R-1, not gated:** bias, MAE, RMSE, worst over-prediction, over-prediction rate at > 2 C and > 5 C, per set; one-sided margin (`one_sided_margin`, direction over, alpha 0.10) per set where it is finite; per-source split of S-T2 and per-method (D93, D56) split of S-CRC.
* **R-3 applicability domain:** counts of ordered and unordered pairs used, share in domain, the class-pair table, the T range of the Brouwer hydrocarbon rows against the flash temperatures solved, and the per-set distribution (median, minimum, 10th percentile) of the mixture-level in-domain pair weight.
* **R-4 mesitylene and the 2.21 C margin (corrected premise).** The yaml margin 2.21 C is the 35-system margin (36 systems without the mesitylene system), set by the ternary nonane + decane + dodecane; the 3.70 C margin with mesitylene included is the one mesitylene sets. Declared reading: the one-sided margin (overlay) for ideal and for GH-GNN-Wilson in the four scopes old-31, old minus mesitylene, old+new (36), old+new minus mesitylene (`delta_lab_flash_li2014/lab_margin.csv` scope names), the difference with versus without mesitylene (flag < 0.5 C as in R-LI-2), the system that sets the quantile in each, and the mean error of the mesitylene system. Mesitylene's error was traced earlier to its Antoine set (Psat ratio 0.52), not to non-ideality; so if GH-GNN appears to "fix" the mesitylene system, the test is repeated with the corrected Antoine (`mesitylene_best_input.json`, same override as `mesitylene_fixed`) for both ideal and GH-GNN, to tell a real effect from gamma compensating a vapour-pressure error. No margin is written to any yaml.
* **R-5 Margules sensitivity (binaries of S-LAB only, secondary).** MAE and over-prediction rate of Margules vs ideal vs Wilson on the same binary rows. Descriptive.
* **Not run:** SVL freezing with GH-GNN (the WS5 freezing item). It is not attempted here unless time remains after R-1 to R-5; if it is not run, none is claimed.

### Stated expectation (before results)

A loss or a tie is the likely outcome. The ideal route already has bias about -3.3 C on the lab set (it under-predicts, the safe side), and most hydrocarbon pairs have gamma_inf > 1, which pushes the predicted flash point lower still, so GH-GNN is expected to help only where the ideal model over-predicts. The UNIFAC precedent was dMAE +0.45 C. A pass would be a surprise and would still not designate anything.

### User decisions (not made here)

Designating GH-GNN (or any gamma model) in the live flash route is a user decision. Writing any margin for it to `configs/uncertainty_v5.yaml` is a user decision. This delta does neither.

## Results

(appended below after the run)
