# L1: flash evaluation re-run on the LFL overlay (v5)

Status: READINGS DECLARED 2026-11-11, before the re-run. Results are written by `scripts/v5/eval_flash_lfl_overlay_v5.py`
to `RESULTS.md`; this file is not edited afterwards (add a dated note).

## Change under test

Only the lower flammability limit of the pure compounds changes: `data/processed/v5/lfl_v5.csv` replaces a Jones-rule estimate
by a tabulated value where one exists (IEC 60079-20-1 table first, Co-Optima compiled value second). Measured Zabetakis values
and every Antoine set stay as they are. The overlay is applied in memory; `compounds_anchored_vp.csv`, `results/v4/` and
`configs/` are not modified.

## Runs

A. v4 flash evaluation, unchanged code (`scripts/v4/eval_flash_v4.py`) on the r3 compositions: 126 lab mixtures (31 systems)
   and the CRC 647 fuels. Baseline for comparison: `results/v4/R3/F/flash_summary.csv`.
B. The AJFTD GCxGC fuels of `results/v5/G1` (ideal and UNIFAC), baseline `results/v5/G1/predictions.csv`.

## Readings (fixed before the re-run)

R1. Report, per slice and model: MAE, bias, fraction over-predicted by more than 5 C, before and after, and how many pure
    inputs changed in each run.
R2. "The overlay is safe to adopt" holds only if, for the ideal model, on each of the three sets (lab mixtures all 126, CRC 647
    D93 dev, AJFTD G1 all eligible): the fraction over-predicted by more than 5 C does not rise, and MAE rises by at most 0.2 C.
R3. "The overlay improves accuracy" holds only if ideal MAE falls by more than 0.2 C on at least two of the three sets and
    rises on none.
R4. The CRC 647 locked slice is shown but takes no part in R2 or R3.

If R2 does not hold the overlay stays a side table. Either way adoption as the default input is the user's decision.
