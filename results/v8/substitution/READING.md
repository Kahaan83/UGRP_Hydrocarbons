# v8 G1 reading: ML-constant substitution on dev lab mixtures

Rule (docs/v8_predeclaration.md): pass if the upper end of the one-sided 90 % paired group-bootstrap bound of dMAE (ML minus measured) is at or below the tolerance. Out-of-fold ML constants; one kind at a time; all components substituted.

| property | systems | skipped | dMAE | 90 % upper | tolerance | verdict |
|---|---|---|---|---|---|---|
| density (kg/m3) | 167 | 0 | 10.91 | 11.79 | 2 | fail |
| viscosity (ln) | 236 | 19 (ml_route_failed) | 0.130 | 0.144 | 0.02 | fail |
| flash (C) | 31 | 0 | 1.33 | 2.07 | 2 | fail |
| freezing (C) | 8 | 0 | 3.82 | 4.81 | 2 | fail |

G1 failed on dev for all four constant kinds. None is used in Tasks 3 to 5; constraints that need an ML-filled constant stay `uncertified`.
Flash fails by 0.07 C on the upper bound only; the rule is not tuned.
Limits: viscosity keeps the measured T-curve shape (optimistic); 19 viscosity systems gave an undefined Refutas blend after substitution (VBN needs nu + 0.8 > 1); freezing has 8 scored systems with ML values; flash and freezing have no locked lab systems.
G4 (locked): not run, because no property passed G1 on dev.
Run: `scripts/v8/substitution_test.py --set dev` (rows.csv, summary.csv).
