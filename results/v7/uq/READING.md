# Task 19 reading: family-conditional and shift-weighted one-sided margins (dev rows only)

Method: `scripts/v7/uq/coverage_by_family.py`. For each held-out group the 90 % one-sided margin is calibrated on all other groups
(global, Mondrian per family with fallback to global below 10 rows, or density-ratio weighted on the model prediction, plus nC for YSI),
and coverage is the share of held-out rows whose unsafe-side score is within it. Unsafe side: flash over-prediction, freezing and YSI
under-prediction (YSI direction assumed, not predeclared). Locked rows excluded (11 flash rows, 53 YSI rows). Numbers in `coverage_by_family.csv`.

Readings (exploratory, dev rows already seen; no v5 or v6 verdict is reopened):
- Flash (68 rows, 59 groups): global coverage 0.912 overall but 0.839 for D93/IP170 (24 groups); Mondrian lifts it to 0.871 (designated) and 0.903 (overlay) with overall 0.897 and 0.912. Weighted gives 0.956 to 0.971 overall, i.e. conservative.
- Freezing (52 groups): the global margin covers 0.826 of crc_jeta fuels for the constant baseline; Mondrian gives 0.957. Ridge: ajftd_posf 0.800 global, 0.933 Mondrian. Families have 3 to 23 groups; fuellib_reference (3) always falls back to the global margin.
- YSI (169 rows, 21 groups): the global margin covers only 0.28 (ridge) and 0.24 (LightGBM) of poly_arom; Mondrian restores 0.96 and 0.92. The margins themselves (about 300 and 700 YSI units globally) are uninformative: the YSI errors are heavy tailed, so coverage is bought with width.
- The weighted variant over-covers almost everywhere (0.92 to 1.00), so the shift correction as specified is conservative here; its tight per-family cells (3 to 6 groups) are not reliable.
- Caveats: per-family cells with few groups (stated per row as n_groups) carry large sampling error; Mondrian margins are calibrated on in-sample dev groups and are not validated on unseen systems.
