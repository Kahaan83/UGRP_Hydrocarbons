# ThermoML phase-bug re-extraction (v5)

Fix: `phase_of()` now delegates to `any_phase()` so real-record tags `eVarPhase` / `eConstraintPhase` are read
(old code read only `ePhase`, so vapour-phase VLE compositions passed `is_liquid(None)` and were used as liquid composition).
Extractor moved to `src/fuelblend/data/thermoml_extract.py` (self-contained; `--out PATH` added). `extras/` untouched.
Regression test: `tests/test_thermoml_phase.py` with REAL trimmed record `tests/data/thermoml_vle_gas_comp_fluid2005_03_032.xml`
(methane+butane vapour pressure, composition variable is eVarPhase=Gas). Fails on the old phase_of (2 failures), passes on the fix.

Run: `PYTHONPATH=src python -m fuelblend.data.thermoml_extract <tgz> --out results/v5/thermoml_reextract/thermoml_hc_v5.csv`
Output `thermoml_hc_v5.csv` (new); old-parser rerun reproduced the existing `extras/.../thermoml_hc.csv` byte-for-byte (rows equal).
Logs: `log_old.txt`, `log_new.txt`. Same 11,923 files, 2 unparseable.

## Row diff (old -> new), hydrocarbon rows kept: 46,747 -> 46,493 (-254)
| prop | n_comp | old | new |
|---|---|---|---|
| Vapor or sublimation pressure | 2 | 5570 | 5376 (-194) |
| Vapor or sublimation pressure | 3 | 174 | 114 (-60) |
| all density, viscosity, kinematic viscosity rows | | unchanged | unchanged |

Dropped rows come from 15 files; 0 rows added. By reason (not-kept counts):
| reason | old | new |
|---|---|---|
| not all-hydrocarbon / no InChI | 612551 | 606244 |
| mixture: unsupported/ambiguous composition scale | 288653 | 295331 |
| mixture: composition incomplete | 5336 | 5387 |
| mixture: inconsistent composition | 2588 | 2420 |
| all other reasons | identical | identical |

(The 254 spurious rows were vapour pressures whose only composition variable was gas-phase; they now fall into the
"ambiguous composition scale"/"incomplete" buckets or are dropped earlier.)

## Downstream
`scripts/v4/eval_lab_known_composition.py` reads `data/processed/v2/thermoml_clean_hc.csv` (same schema, a cleaned copy), not
`thermoml_hc.csv`, and uses liquid `Mass density` n_comp>=2 only; density/viscosity rows are identical old vs new, so it and the
tier-1 density fits are unaffected. NOT rerun. Note `thermoml_clean_hc.csv` contains all 254 now-removed vapour-pressure rows
(any vapour-pressure / Antoine refit consumers of it are the ones affected); not regenerated here.
