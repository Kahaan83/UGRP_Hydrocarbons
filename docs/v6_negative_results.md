# v6 negative results (KEEP in methodology and code)

Same convention as `docs/v5_negative_results.md`: each row is a hypothesis with a gate declared in the results README before the locked
holdout was scored. Numbers come from the CSVs named in the last column.

| Track | What was tried | Outcome (declared gate, not moved) | Where |
|---|---|---|---|
| A freezing | Monotone LightGBM residual on UNIFAC-SvL with isomer-tile constants (best dev base), 4 coarse features, fixed settings | Dev LOGO (52 groups): MAE 3.31 vs constant 6.64, CI of the gain [1.56, 5.25] (skill). Locked holdout (12 groups, scored once): MAE 9.30 vs 6.86, CI [-6.31, +1.31], no skill. The dev gain did not replicate. Status stays "no skill shown". | `results/v6/freezing/dev_logo_summary.csv`, `holdout_summary.csv` |
| A freezing | Ridge on coarse composition (c2); Bell-inspired single-species control curve (d) | No skill on dev (c2 6.43, d 5.92 vs constant 6.64; CIs include 0) and on holdout (6.06, 7.08 vs 6.86). Method details of Shi 2017 and Bell 2025 were not verified, so these are not replications. | same |
| A freezing | Fixed mass-weighted average of representative Tm (c1) | Fails: dev MAE 16.72, holdout 23.99, bias +16.7 and +24.0 C. Linear averaging of Tm is not a liquidus. | same |
| A freezing | v5 UNIFAC-SvL (b) on the extended set | Dev MAE 8.83 vs constant 6.64 (CI [-4.75, +0.28]); holdout 13.48 vs 6.86, driven by AJFTD 325 (-144.7 C predicted vs -81 measured) and 97 (-80.7 vs -61). Also -176 C for AJFTD 64 (all iso-paraffin, no constants). Reproduces v5 on the 50 headline rows (MAE 8.4246, bias +5.8778). | `baseline_b_preds.csv`, `repro_gate.json` |
| A freezing | Isomer-averaged tile constants (Hall 2022 idea) | Not testable with current data: only 86 of 1042 designated-tmap entries have both Tm and dHfus, so only 4 tile representatives change, and results equal (b) to 2 decimals. Needs Tm and dHfus for the missing isomers (Track B). | `src/fuelblend/data/tile_isomers.py`, `dev_logo_summary.csv` |

## Reading
- Dev LOGO is dominated by 34 AV-33-22 fuels (65 % of dev). On the 15 dev AJFTD fuels the LightGBM residual had MAE 4.93 vs constant 9.37; on
  the 12 AJFTD holdout fuels it was 9.30 vs 6.86. A 52-group LOGO with a CI that excludes 0 was not enough evidence; the holdout is the gate.
- One-sided margin from dev out-of-fold residuals (90th percentile of measured minus predicted, applied upward): 3.98 C for the LightGBM
  residual. On the holdout it covered 8 of 12 fuels (0.67), below the 0.90 target, so it is not adopted. Live freezing margin stays 13.84 C.
- Nothing in v6 Track A changes `predict_fuel`. Do not re-score the holdout; `holdout_scored.flag` blocks it.

## Track B: pure-compound D-MPNN (2026-11-21)
- Tm head (MAE 30.4 K) does not beat the v5 fusion_estimate (28.4 K) on the same scaffold-grouped folds; not adopted for Track A.
- Tb and LHV heads not trained (no Tb source; LHV units unverified). Details: `results/v6/gnn/README.md`.

## Track C mixture ML retry (2026-11-21, `results/v6/mixture/`)
Frozen D-MPNN encoder + MolPool + pair term, residual over physics, gates declared before fitting. Failed: KAUST RON and MON (no gain
over GBM or linear molar), Boehm viscosity (dev and locked loss, CI < 0) and Boehm density (dev loss; head reverted to no correction).
Passed: Leenhouts MON only (dev gain CI [0.04, 1.33], locked [0.07, 0.98] MON units vs descriptor GBM), single head of five, no
multiplicity correction. The "corrected-pure" baseline was worse than linear molar everywhere.
