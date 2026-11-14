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

### R-1 gate (declared, applied unchanged): NOT MET

| set | rows / clusters | MAE ideal | MAE GH-GNN-Wilson | dMAE [95 % CI] | over > 5 C ideal -> GH | pass |
|---|---|---|---|---|---|---|
| S-LAB | 126 / 31 | 3.692 | 4.393 | +0.701 [+0.145, +1.333] | 0.8 % -> 0.8 % | no |
| S-CRC | 45 / 45 | 4.094 | 4.298 | +0.204 [-0.302, +0.721] | 4.4 % -> 2.2 % | no (CI includes 0) |
| S-T2 | 70 / 70 | 2.576 | 3.665 | +1.089 [+0.849, +1.314] | 5.7 % -> 5.7 % | no |

Two of three gated sets are significantly WORSE. The unsafe-side rate did not rise anywhere, but it did not fall on the lab or T2 sets either. Locked CRC (12 fuels, descriptive, read once): dMAE -0.106 [-1.016, +0.809], over > 5 C 8.3 % -> 8.3 %; no unsafe rise. Bias moves further to the safe side everywhere: S-LAB -3.25 -> -4.20 C, S-CRC -0.73 -> -2.52, S-T2 -0.69 -> -2.18. The stated expectation (a loss is likely) was right: gamma > 1 for most pairs lowers the predicted flash point on top of an ideal model that already under-predicts.

### R-2 splits (descriptive; full numbers in `summary.csv`)

* S-T2 by source: D56 dMAE +1.40 [1.05, 1.68]; D93/IP170 +0.68 [0.33, 1.02]; other +1.54 [1.03, 1.91]. S-CRC is all D93 (same as above).
* S-LAB binaries Wilson +0.70 [0.12, 1.33]; Margules +1.08 [0.32, 1.92] (R-5: Margules is worse than Wilson). Ternaries: 1 system, 3 rows, +0.86.
* S-LI (6 systems, margin reading only) +0.30 [0.02, 0.71]. S-LAB + S-LI (36 systems) +0.62 [0.16, 1.12].
* Wilson inversion failures (declared rule L = 1): 0 pairs in every set (about 66 000 pairs). Solved flash temperatures: lab -5 to 151 C, CRC 40-60 C, T2 16-67 C.

### R-3 applicability domain

* 4970 ordered hydrocarbon pairs examined (about 58 unique molecules); 95.2 % in domain; 1867 are interpolation pairs (solute and solvent both seen in training). All 240 out-of-domain pairs fail the class-pair count (< 25 training systems), none fail Tanimoto (> 0.35 always).
* Mixture-level in-domain pair weight: lab and Li 1.00 (every mixture fully in domain); CRC 647 median 0.69 (minimum 0.54, 0 of 45 fully in domain); T2 median 0.75 (3 of 70 fully in domain). The real-fuel reconstructions lean on the heavy-hydrocarbon class pairs, which GH-GNN's training set barely covers.
* Training temperature range (Brouwer hydrocarbon rows, 2372 rows, 805 pairs): 0-180 C, 5th-95th percentile 20-139 C; the flash temperatures solved (-5 to 151 C) sit in range except a few lab points below 0 C.
* AD-gated variant (out-of-domain pairs set to L = 1): S-CRC +0.059 [-0.146, +0.270], S-T2 +0.670 [0.48, 0.85], S-LAB unchanged (all in domain). AD gating reduces the loss but does not turn it into a gain; it also is not the gated variant.

### R-4 mesitylene and the lab margin

* Ideal route reproduces the yaml numbers: 3.70 C with the mesitylene system (36 systems), 2.21 C without it (set by another system); the mesitylene-Antoine-fixed ideal also gives 2.21 C. So the caller's 2.21 C is the no-mesitylene figure.
* GH-GNN-Wilson one-sided lab margin is 0.70 C with or without mesitylene (difference 0.0, within the 0.5 C flag); on the 31 old systems 0.26 C. **This is not evidence that GH-GNN fixes the mesitylene system.** The margin falls only because every prediction moves 1 C lower (bias -3.31 -> -4.20 C), so the over-prediction tail the margin is the quantile of shrinks. MAE got worse. With the corrected mesitylene Antoine, the mesitylene system itself goes from -4.74 (ideal) to -7.34 C (Wilson): GH-GNN pushes it further into under-prediction. A smaller margin bought by uniformly lowering predictions is a conservatism trade, not a model gain, and the gate (R-1) already failed.
* No margin is written to any yaml. Any use of the 0.70 C figure is a user decision.

### R-5 Margules

Binaries of S-LAB: Margules MAE 4.77 vs Wilson 4.39 vs ideal 3.69; over > 5 C identical (0.8 %). Margules is worse than Wilson and is not recommended.

### Not run

SVL freezing with GH-GNN was not run (time was spent on the flash gate); nothing is claimed about it. UNIFAC-Dortmund fallback not needed (GH-GNN ran and reproduces the repo's stored predictions to 2.4e-4 ln units).

### Reproduce

`python scripts/v5/ws5/ws5_flash.py --stage all` then `python scripts/v5/ws5/ws5_report.py` (GH-GNN clone and weights under `data/raw/ghgnn_code/`, see its `SOURCES.md`; MIT licence).

### Needs user decision

Designating GH-GNN: not recommended, gate failed. Writing any margin: none written.
