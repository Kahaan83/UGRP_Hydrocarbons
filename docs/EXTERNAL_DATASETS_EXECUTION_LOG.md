# External Fuel Datasets: Execution Log (corrected)

> **Branch**: `feat/external-datasets-ingest`. Staging only; production paths (`data/unified/`, `src/fuelblend/`, `scripts/v5/`) are not edited by the ingest scripts.

## Correction note (supersedes the earlier version of this log)

An audit of the earlier version found these defects. All are corrected below or removed.

1. **In-sample "conformal" ±24.6% was invalid.** The margin (1.246x) was computed from the same Boehm rows used to fit the model, so it was not a held-out conformal margin. The claim that it can replace the production cold-viscosity margin is withdrawn. Boehm is also already in production (`src/fuelblend/data/loaders_boehm2022.py`), so the Boehm staging scripts, CSVs and reports were deleted.
2. **Octane headline claims were wrong.** The earlier "linear molar MAE 3.66 / GBM slashes error" numbers used typed pure RON/MON literals (1,2,4-TMB 148/124 versus the dataset's own measured 109.5/108.0; cyclopentane and 1-hexene also differed), dropped diisobutylene from all 26 MC blends and renormalised the rest, and grouped by a system string that put a hydrocarbon and its ethanol blends in different folds (leakage). The literature comparison row (~1.3 / ~1.4 ON) was unsourced and is removed.
3. **The stale "inf" margin story**: the log presented the production cold margin as a live `inf` defect fixed by this work. That was not verified against the current production code and is withdrawn.
4. **"Zero production modifications" was an unverified claim.** `git status` on this branch shows modified tracked files outside the ingest area (`results/v5/BASE/*`, `splits/locked.json`, `data/real_fuels/v5/*`). This log makes no such claim. The ingest scripts themselves write only to `data/interim/external_staging/` and `results/external_benchmarks/`.
5. Other defects fixed: only 7 of 128 Table 3 compounds had a structure (now 128 via a curated SMILES table and RDKit); Table 5 (FACE + ethanol) was not ingested; single-component rows were labelled `multi`; invented InChIKeys for TMB, 1-hexene and diisobutylene in the FACE resolver; unsourced `FACE_BASE_PIONA` literals; hardcoded 114+30 counts (the run gives 120+24); AI4fuel density at 293.15 K for 15-30 C values, melting labelled as freezing, CN labelled D613, ExactMolWt used for MW; `system_id` not following the production convention.
6. **Removed as duplicates of production**: Boehm cold viscosity, CRC FACE-1 diesel, ternary flash staging (scripts, CSVs, reports; recoverable from git history).

## 1. Decisions in force

Kept: KAUST (octane mixture validation; oxygenates allowed for octane only, tagged by `is_hc`), FACE surrogate rows, AI4fuel (pure-compound pool for a future multi-task model). Design certification remains model-level only.

## 2. Scripts (`scripts/external_ingest/`)

| Script | Role |
|---|---|
| `_common.py` | Imports the production `system_id` (`fuelblend.data.thermoml_sle`, md5 of `"|".join(sorted InChIKeys)`, 12 chars); RDKit helpers |
| `kaust_registry.py` | Curated name-to-SMILES for all 128 Table 3 compounds and the blend components; pure RON/MON read from measured rows; densities from `data/unified` or CRC (AI4fuel) |
| `parse_kaust_jameel.py` | Parses Tables 3, 4, 5 from the PDF |
| `ingest_kaust_jameel.py` | Canonical KAUST rows (pure, Table 4, Table 5), compounds table |
| `eval_octane_blending_benchmarks.py` | Benchmark and computed report |
| `resolve_face_gasoline_dha.py` | Resolves 144 legacy FACE rows; computed report |
| `parse_sjtu_ai4fuel.py`, `ingest_sjtu_ai4fuel.py` | AI4fuel ingest with provenance tier, T ranges, dedupe flags |
| `ingest_leenhouts_mon_multicomponent.py` | Leenhouts MON mixtures (production `system_id`) |

## 3. Staged holdings (computed row counts)

| File | Rows | Systems | Notes |
|---|---|---|---|
| `kaust_jameel_canonical_measurements.csv` | 562 (281 RON, 281 MON) | 145 | 274 pure, 90 binary, 198 multi; 68 ethanol-containing blends (`is_hc` false); tier `literature_compiled` |
| `canonical_ai4fuel_measurements.csv` | 7,709 | 1,763 | tier `compiled`; 978 rows duplicate a (property, system_id) already in `data/unified`; compounds overlapping `compounds.csv`: 603 |
| `canonical_leenhouts_mon_mix_measurements.csv` | 292 | 67 | MON mixtures |
| `canonical_resolved_face_gasoline.csv` | 144 | 16 | 120 surrogate + 24 FACE + ethanol |

## 4. Octane blending benchmark (KAUST Table 4, computed)

Full report: `results/external_benchmarks/octane_blending_report.md`. 120 blends scored (TMB+0, CP+0, Hex+0 define the pure values and are excluded). GroupKFold, 5 splits, groups = component set without ethanol (8 groups). Pure values: TMB 109.5/108.0, cyclopentane 100.0/85.6, 1-hexene 73.6/64.5 (+0 rows); diisobutylene 100.6/86.5 (Table 3 row 72 isomer proxy); n-heptane 0/0 (PRF definition).

| Target | Model | MAE | RMSE | Bias | Max |
|---|---|---|---|---|---|
| RON | Constant (train mean) | 12.476 | 14.825 | +1.541 | 49.38 |
| RON | Linear volume | 5.748 | 7.986 | -5.127 | 23.9 |
| RON | Linear molar | 2.791 | 3.783 | -1.141 | 13.87 |
| RON | Linear molar + GBM residual | 2.453 | 3.161 | +0.505 | 11.376 |
| MON | Constant (train mean) | 10.021 | 12.113 | +1.140 | 43.544 |
| MON | Linear volume | 5.644 | 7.751 | -4.767 | 23.6 |
| MON | Linear molar | 3.101 | 4.364 | -2.303 | 12.799 |
| MON | Linear molar + GBM residual | 2.092 | 2.815 | +0.133 | 10.156 |

FACE + ethanol (Table 5, 24 blends), linear volume from measured FACE+0 and ethanol: RON MAE 7.044, MON MAE 4.960 (both biased low).
The GBM residual is cross-validated across component sets, but with 8 groups and 120 rows it is a small-sample result; per-family it is worse than linear molar for TMB and cyclopentane (see report). Diisobutylene density sensitivity (±5%) changes linear-molar MAE by at most 0.013 ON.

## 5. Isolation audit

The ingest scripts write only under `data/interim/external_staging/`, `results/external_benchmarks/` and (tests) `tests/external_ingest/`. Tests in `tests/external_ingest/` check InChIKeys against RDKit, KAUST pure values against the +0 rows, DIB presence in all 26 MC blends, the `system_id` convention, and the AI4fuel corrections.
