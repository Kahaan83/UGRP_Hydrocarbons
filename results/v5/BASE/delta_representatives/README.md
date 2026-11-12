# delta_representatives (Task 12)

Declared 2026-11-12 before `delta_representatives.py` exists; this commit contains no result file. Frozen reference: `5a62e20`.

Candidates (audit proposal 8): 1-methylnaphthalene, 2,6-dimethylnaphthalene, perhydrophenanthrene, perhydrofluorene, 1-methyltetralin,
n-butyltetralin, 2-methyldecalin, 1,2,4-trimethylbenzene, n-propylbenzene. A candidate is admitted only if the registry already holds:
Antoine constants with source `refit_measured` or `webbook`; an LFL that is not a Jones-rule estimate (measured, or tabulated in the adopted v5 overlay);
a measured pure viscosity fit; and a valid provider density. Fusion data is recorded but not required. Admitted candidates replace their own exact
(subkind, carbon number) bin only (`REPLACED_BINS`); every other bin returns what `SubReps` returns.

Re-scoring: only fuels whose mapped composition changes are re-predicted (flash designated and LFL overlay; viscosity at 313.15 K; density at 288.15 K);
all other fuels keep their frozen predictions. Diesel: unmapped mass of the 9 FACE diesels with `SubReps` and `SubRepsV5` (S6 mapping).

Reading (plan text, fixed now): "unmapped mass falls below 10 % for at least 7 of 9 FACE diesels, and T2 viscosity at 40 C and density bias do not worsen
by more than 0.02 ln / 2 kg/m3". Candidates with missing constants are listed in `candidate_gaps.csv` with the source that would fill them.

## Result 2026-11-12

Run with `scripts/v5/base/delta_representatives.py`. Re-predicted "old" values reproduce the frozen flash rows to 5e-5 C (69 scored fuels checked).
* Admitted: 1,2,4-trimethylbenzene and n-propylbenzene (refit Antoine, Co-Optima LFL, measured viscosity fit). Both sit in the C9 alkylbenzene bin; 1,2,4-trimethylbenzene
  replaces mesitylene there (`REPLACED_BINS = {(alkylbenzene, 9)}`). 140 T2 fuel entries change composition. Every diaromatic, tricycloparaffin, tetralin and decalin candidate is
  refused: missing Antoine/LFL in the registry, or not in it at all (`candidate_gaps.csv`).
* Reading FAILS on diesel: unmapped mass below 10 % for 4 of 9 FACE diesels with SubRepsV5, unchanged from SubReps; the unmapped mass is diaromatics and tricycloparaffins,
  for which no candidate was admitted.
* T2 half passes: pooled viscosity at 40 C |bias| +0.009 ln (0.042 -> 0.051), density at 15 C |bias| -0.25 kg/m3 (-3.80 -> -3.55).
* Not in the reading but material: flash gets worse. Pooled T2 bias -0.73 -> -0.99 C, MAE 2.58 -> 2.74 C (with the LFL overlay: -0.69 -> -0.98 C, MAE 2.58 -> 2.77 C).
  Unsafe rate unchanged (5.7 %). Same direction as the mesitylene Antoine fix, which was also rejected.
Recommendation: do not adopt SubRepsV5. Diesel needs diaromatic (1-methylnaphthalene Antoine + LFL) and tricycloparaffin constants added to the registry first.
