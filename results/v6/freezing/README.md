# v6 Track A: freezing point (declared before any model beyond baseline (b) was scored)

Status of the v5 reference: "status: no skill shown" (`docs/final_numbers_v5.md`): pooled n 50, MAE 8.425 C vs constant 6.413 C.

## Data and lock
- `dataset.csv`: 64 groups = v5 T2 entities plus AJFTD 21, 22, 64, 154, 155 (dev) and 20 (holdout). AJFTD 32 is the same fuel as FuelLib
  posf10264 and is not added. One entity per `dup_group`. Rule and hashes in `holdout_lock.json` (committed as `dd65d08`).
- Dev = 52 groups (leave-one-group-out, LOGO). Holdout = 12 groups, scored once, after the model list is frozen. Holdout = v5 S4 locked 8, plus
  AJFTD 86, 195, 325 and 20, which are in the older `splits/locked.json` that the v5 headline did not honour.
- Fuel 64 (iso-paraffin, 99.6 %): uncensored reading -90 C, censored `<-61` readings dropped. Baseline (b) gives -176 C for it (no constants for
  the heavy iso-paraffin tile). Every result is reported with and without fuel 64.
- Reproduction gate (`repro_gate.json`): baseline (b) reproduces v5 on the 50 headline rows (MAE 8.4246, bias +5.8778, max abs diff 3e-6 C, the
  rounding of the stored v5 rows).
- Baseline (b) predictions were also computed for the 12 holdout groups at step 2 (physics, no fitting, nothing selected on them).

## Baselines (all scored LOGO on dev, then once on holdout)
- (a) constant: mean of the other dev groups.
- (b) v5 UNIFAC-SvL, designated constants, unchanged code.
- (c1) fixed mass-weighted average of representative-tile Tm (no fitting). (c2) ridge on coarse features (n-paraffin %, n-paraffin C14+ %,
  aromatic %, mean carbon number), alpha by inner LOO. Neither is claimed to be Shi 2017 or Watanasiri 2026 (method details not verified).
- (d) single-species control curve, Bell-inspired (method details unverified): isotonic fit of freezing point on one n-paraffin carbon-number
  mass % (or aromatic %); species chosen by inner LOO MAE inside each fold.
- (e) UNIFAC-SvL with isomer-averaged tiles: each representative's Tm and dHfus replaced by the mean over the registry isomers of its tile
  (`SubReps.by[(kind, n)]`) that have constants in the designated tmap (Hall 2022 idea; `src/fuelblend/data/tile_isomers.py`).

## Residual model (declared a priori)
- LightGBM on the residual (measured minus base), base = whichever of (b), (e) has the lower dev OOF MAE, base clipped to [-100, -20] C. Fixed
  settings, no tuning: num_leaves 4, min_child_samples 8, n_estimators 60, learning_rate 0.05, no subsampling.
- Features and monotone constraints (on the residual added to the base): n-paraffin C14+ % (+1), n-paraffin % (+1), mean carbon number (+1),
  aromatic % (0).

## Gates
LOGO on dev; one-sided margin = 90th percentile of out-of-fold (measured minus predicted), added to the prediction; beat-trivial = paired
bootstrap CI (over fuels) of (MAE constant minus MAE model) above 0; holdout scored once. A model that fails these is logged in
`docs/v6_negative_results.md`.

## Disclosures made after the dev run (before holdout scoring)
- (e) first used only the 77 compounds that SubReps accepts (4 tiles changed, all C6 to C10) and equalled (b) exactly. It was then widened to
  every registry compound with constants (564 structures in sub-kind tiles), which still changes only 4 representatives: only 86 of 1042
  designated-tmap entries have both Tm and dHfus. (e) therefore equals (b) on dev to 2 decimals. Isomer averaging is data-limited, not
  refuted; it needs Tm and dHfus for the missing isomers (Track B).
- Dev LOGO results are in `dev_logo_summary.csv`; the model list and settings above were frozen at this commit, before `step4_holdout.py`.

## Results (dev LOGO n 52, holdout n 12; MAE in C; full tables in dev_logo_summary.csv and holdout_summary.csv)
| model | dev MAE | dev bias | holdout MAE | holdout bias | holdout vs constant |
|---|---|---|---|---|---|
| a constant | 6.64 | 0.00 | 6.86 | +3.54 | trivial |
| b UNIFAC-SvL (v5) | 8.83 | +5.51 | 13.48 | -0.69 | worse |
| c1 weighted Tm | 16.72 | +16.72 | 23.99 | +23.99 | worse |
| c2 ridge | 6.43 | +2.66 | 6.06 | -0.95 | no skill (CI [-1.90, +4.07]) |
| d control curve | 5.92 | -1.50 | 7.08 | -1.50 | no skill |
| e isomer tiles | 8.83 | +5.51 | 13.48 | -0.69 | equals b (data-limited) |
| f LightGBM residual | 3.31 | -0.16 | 9.30 | -1.77 | no skill (CI [-6.31, +1.31]) |

Verdict: status "no skill shown" is unchanged. The dev LOGO gain of (f) did not replicate on the locked holdout. See `docs/v6_negative_results.md`.
Holdout was scored once (`holdout_scored.flag`, entry in `reports/LOCKED_ACCESS.log`).

## Post-hoc step 5: tiles filled with Track B constants (exploratory, DEV ONLY; `scripts/v6/freezing/step5_tiles_gnn.py`, `tiles_gnn/`)
Written after the holdout was scored; the holdout is not used and these variants are not part of the declared result. The registry has 1,042
compounds, of which 172 constant slots (Tm and dHfus) were designated, 452 filled from measured pool values and 1,460 estimated. e2 uses the
Track B D-MPNN for Tm and dHfus; e3 uses the v5 fusion_estimate for Tm and the D-MPNN for dHfus (D-MPNN dHfus MAE 12.8 kJ/mol, n 81; Tm MAE
about 30 K). 31 representatives change (4 before). Dev (n 52, no fitting, so no LOGO needed):

| model | dev MAE | bias | CI of gain vs constant | no fuel 64 MAE |
|---|---|---|---|---|
| a constant | 6.64 | 0.00 | trivial | |
| b / e (v5 tiles) | 8.83 | +5.51 | [-4.75, 0.28] | 7.31 |
| e2 tiles, D-MPNN fill | 7.88 | +7.88 | [-2.98, 0.85] | 7.31 |
| e3 tiles, fusion Tm fill | 8.09 | +8.09 | [-3.21, 0.66] | 7.32 |

Reading: the whole improvement (8.83 -> 7.88) is fuel 64, whose tile had no constants; without it e2 and e3 equal the v5 value (7.31), and every
fuel is over-predicted (frac_under 0). Neither beats the constant. Status unchanged: no skill shown. Filling tiles with estimated constants
does not help the freezing point on the other 51 dev fuels.
