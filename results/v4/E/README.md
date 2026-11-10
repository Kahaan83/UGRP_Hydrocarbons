# v4 Task 16 - residual set encoder vs Task 11's best viscosity model

**Outcome: the encoder LOST (not promoted to `model_decision`).** It ties Task 11's structured ridge: the paired fuel-cluster bootstrap CI of the
RMSE difference includes 0 under both protocols. (`promotion_decision.json`: `PROMOTE_TO_MODEL_DECISION = false`; the script now gates on fuel-level dev coverage and records locked fuel coverage.) **Coverage: row coverage (0.970) is inflated** (3 rows per fuel, but the conformal score is per fuel); the fuel-level numbers are
dev cross-conformal 41/45 = 0.911 and **locked 8/12 = 0.67, below the 0.85 bar**. Any future winner must be gated on fuel-level coverage (dev and locked), not row coverage. Numbers below are copied from the CSVs in this folder (`viscosity_encoder_vs_task11.csv`, `conformal_calibration.csv`,
`flash_encoder_real_d93.csv`, `lab_groupkfold_system.csv`); the sum-pool ablation is in `../E_sum_ablation/`.

## What was run (CPU only; no GPU result exists)
`python scripts/v4/train_set_encoder.py --device cpu --seeds 5 --epochs 40 --threads 3 --pools attention --lab-folds 5 --out results/v4/E`
(sum-pool ablation: same with `--pools sum --lab-folds 0`, run concurrently, 3 threads each; torch 2.15.0+cpu). Wall time 1605 s (attention), 800 s (sum).
* 5-seed deep ensemble per fold (seeds 20261110-14), 40 epochs, AdamW lr 3e-3 wd 1e-3 cosine, hidden 64, batch 512 (lab) / 32 (real), Huber loss.
  Each epoch = ceil(14237/512) = 28 steps, so ~1100 optimiser steps per model.
* **Hyper-parameters were fixed a priori; nothing was selected on CV or locked rows.** Epoch count came from the compute budget. This is a small, under-tuned,
  under-trained CPU run (no width/lr/epoch search, no early stopping; real-fuel rows are oversampled ~7x per epoch relative to a single pass: 28 steps x 32 rows / ~120 real viscosity training rows per fold, a further over-fit risk); a different configuration could do better or worse. No tuning on locked.
* Heads: viscosity (residual on ln nu vs raw Refutas, structured `a0 + a1*(1000/T - 3.7)`, exactly linear in 1000/T) and flash (degC residual on the Task 9 ideal
  Le Chatelier value). Domain flag (lab=0/real=1) is a head input. Masked multi-task Huber loss: one minibatch per (domain, head) table per step.
* Training data: 14,237 lab dynamic-viscosity rows (322 systems, same filter as Task 11's lab row), 126 lab flash records (31 systems), 171 CRC647
  viscosity rows and 57 CRC647 D93 flash rows (reconstructed compositions, so errors mix reconstruction and physics).
* Outer CV: 9-fold grouped by fuel over the 45 dev CRC fuels. Locked (12 fuels, 36 viscosity rows) never enter a training/calibration index
  (`split_audit.json`, 9 fold checks, plus `tests/test_v4_set_encoder.py`); final model = dev-trained, scored once on locked.
* Reference models: Task 11 structured ridge as published (LOO-by-fuel, 44 training fuels per fit: more data per fit than the 9-fold encoder) and the same
  ridge recomputed on the identical 9 folds with inner-CV alpha (the fair row).

## Viscosity (ln nu RMSE, all T; 95% cluster-bootstrap CI)
| model | dev (45 fuels, 135 obs) | locked (12 fuels, 36 obs) |
|---|---|---|
| Raw Refutas | 0.4267 [0.3959, 0.4549] | 0.3806 [0.3135, 0.4564] |
| Recalibrated Refutas (T11 LOO) | 0.1351 [0.1157, 0.1550] | 0.1939 [0.1341, 0.2594] |
| Structured ridge, Task 11 (LOO) | 0.0826 [0.0731, 0.0949] | 0.1497 [0.0810, 0.2289] |
| Structured ridge, same 9 folds | 0.0833 [0.0739, 0.0953] | 0.1497 [0.0810, 0.2289] |
| **Set encoder, attention pool (5 seeds)** | 0.0862 [0.0758, 0.0983] | 0.1500 [0.0852, 0.2248] |
| Set encoder, sum pool (ablation) | 0.0852 [0.0753, 0.0966] | 0.1504 [0.0879, 0.2240] |

Paired dRMSE (encoder - ridge), dev overall: vs T11 LOO +0.0036 [-0.0004, +0.0088]; vs same-fold ridge +0.0029 [-0.0022, +0.0100] -> CI includes 0, no win
(point estimates favour ridge). Sum pool: +0.0026 [-0.0016, +0.0078] / +0.0020 [-0.0034, +0.0086]. Locked: +0.0003 [-0.0060, +0.0145] (12 fuels, uninformative).
Per-T slices are in the CSV; only locked at -40 C has a CI excluding 0, in the direction ridge better (+0.0119 [+0.0028, +0.0218]) - a post-hoc slice of 12 fuels, not used for the decision.

## Uncertainty (Task 10 split-conformal functions, dev groups only; alpha 0.10)
Two-sided |ln residual| of the ensemble mean, group = fuel, n = 45, k = 42: margin 0.170 ln (x1.185). Dev cross-conformal coverage (margin from the other folds'
out-of-fold residuals): 131/135 rows = 0.970 (inflated, see headline), 41/45 fuels (all three temperatures covered) = 0.911. Locked, reported as counts, margin from all dev groups:
31/36 rows, 8/12 fuels. The 12-fuel locked count is low-powered, but is below the dev figure. Ensemble seed spread (mean sd 0.0076 ln) is far smaller than the error,
so deep-ensemble variance alone is NOT a usable uncertainty here; coverage comes from the conformal margin. Flash (one-sided over-prediction, D93 fuels): margin
+2.82 C (applied = max(q,0)), locked 8/12 fuels covered. Margins computed here are NOT written to `configs/uncertainty_v4.yaml`.

## Flash and lab (secondary)
* CRC D93 (57 fuels): ideal MAE 9.44 C dev / 8.44 locked (all under-predicted, -9.2 C bias); fold-trained constant offset 1.67 / 1.92 C; encoder 1.70 / 1.99 C
  -> the encoder only relearns the constant offset; max over-prediction dev 4.9 C (offset: 7.0 C), locked 3.6 C.
* Lab GroupKFold(system_id, 5 folds), attention encoder: viscosity ln-RMSE 0.3081 vs lab ridge 0.3093 vs Refutas 0.3096 (RMSE dominated by a few bad systems; CI wide
  [0.04, 0.55]); MAE 0.0451 vs 0.0585 vs 0.0673; the ln-RMSE differences (0.3081 vs 0.3093 vs 0.3096) are ties within the CI, and the MAE gap has no CI computed, so treat lab viscosity as a tie. Flash lab (126 records, 31 systems, system-grouped folds): MAE ideal 3.74 C, ideal + constant offset fitted in training folds 2.68 C (RMSE 3.28), encoder 2.50 C (RMSE 3.43). So most of the 3.74 -> 2.50 gain is bias removal; the encoder adds ~0.2 C MAE over the constant offset and is worse on RMSE - not distinguishable from bias removal with 31 systems.

## Caveats
* Single CPU configuration; undertrained/untuned. Conclusion is "no evidence the encoder beats the structured ridge on 45 dev fuels", not "set encoders cannot".
* Dev = 45 fuels -> CIs are wide (+-0.005 ln); a difference below ~0.01 ln is not resolvable.
* Real-fuel compositions are reconstructed (formula-derived); component features are the registry's static descriptors of class x carbon representatives.
* GPU pack in `jobs/v4_gpu/` is untested on GPU; only the CPU result above is claimed.
* Claims register (not edited): CLM-19/19a (neural residual seed jitter - this run uses explicit seeds), CLM-28 (conformal margin - viscosity encoder margin would not replace the ridge's).
