# ThermoML solid-liquid equilibrium (liquidus) rows for hydrocarbon mixtures

Extraction and counts only. No scoring, no model fitting.

- Code: `src/fuelblend/data/thermoml_sle.py` (run `PYTHONPATH=src python -m fuelblend.data.thermoml_sle`)
- Summary counts: `scripts/v5/sle_summary.py`
- Output: `data/processed/v5/thermoml_sle_hc.csv` (849 rows, 22 source files / DOIs)
- Test: `tests/test_thermoml_sle.py` with a real record block, `tests/data/thermoml_sle_benzene_cyclohexane_fluid2007_01_029.xml`
- Inputs (read-only): `extras/fuelblend_data/data/raw/thermoml/ThermoML.v2020-09-30.tgz` (11,923 XML files; the 11,932 in the run log include the 9 of the second archive) and `data/raw/thermoml_2020/thermoml_jced_2020_unlisted.tgz` (9 files, 0 contain SLE data).

## How ThermoML encodes SLE (from real records)

- Property `Solid-liquid equilibrium temperature, K` (PhaseTransition group), `ePropPhase` = Liquid.
- `PhaseID` blocks: one Liquid, one Crystal; the Crystal block's `RegNum/nOrgNum` is the component that crystallises (several Crystal blocks = several solids; a Crystal block without RegNum = solid unidentified).
- The composition `Variable` carries `VarPhaseID/eVarPhase` = Liquid and is the liquid composition at the liquidus point. The last component is implied. The parser mirrors the existing ThermoML extractor rules (Constraint > ConstraintID > ConstraintType, ePropPhase, VarPhaseID).
- Eutectic points are separate datasets: `Eutectic temperature, K` (no composition variable) and `Eutectic composition: mole|mass fraction` (one component's fraction, nOrgNum under Property-MethodID). We pair them only when exactly one composition dataset with the same component set exists in the file.

## Filters

- Keep only liquid-phase composition variables (non-Liquid composition variables are ignored; a point with no complete liquid composition is dropped).
- Mole fractions as given; mass fractions converted with RDKit molecular weights; volume fractions, molality and ratios dropped.
- Pressure > 1000 kPa dropped. Rows with no stated pressure kept (P_kPa empty).
- Hydrocarbon-only: all components C/H, neutral (same rule as `thermoml_extract.is_hydrocarbon`). `is_hc` is therefore always True in the main CSV.
- Pure-component rows (n_comp 1, including 2,747 `Normal melting temperature` rows) are counted but not written.
- Exact duplicates (same file, point kind, components, x, T, P, solid) removed (4).

Rows dropped, hydrocarbon systems only: 431 pressure > 1000 kPa, 43 molality (one file, naphthalene solubility, 10.1021/acs.jced.9b00064), 9 with no composition variable (10.1016/j.jct.2015.06.021, phenanthrene + anthracene, solid unidentified). Non-hydrocarbon dropped: 25,604 point-level rows (counted before filtering by property type, mostly salts, polar organics, solubility data), plus the non-HC share of the parse-stage reasons in the run log.

## Counts (`scripts/v5/sle_summary.py`)

| item | value |
|---|---|
| rows | 849 (838 liquidus, 11 eutectic) |
| rows by n_comp | 2: 842, 3: 7 |
| distinct systems (`system_id`) | 60 (59 with liquidus rows, 11 eutectic) |
| systems by n_comp | 2: 58, 3: 2 (>= 3 components: 2) |
| systems with a jet-range compound (C8-C16 saturated/cyclo or aromatic) | 47 systems, 759 rows |
| systems with an n-alkane C10+ | 33 systems, 459 rows |
| systems with all components in `data/unified/compounds.csv` (InChIKey match) | 60 of 60, 849 rows |
| liquidus rows at T <= 253.15 K / <= 233.15 K | 96 (15 systems) / 15 |
| liquidus T range | 204.75 to 490 K |

The two ternary systems: 1,4-dimethylbenzene + hexane + hexadecane and 1-methylnaphthalene + decahydronaphthalene + hexatriacontane. There is no jet-fuel-like multicomponent (>= 4) SLE data in the archives.

## Caveats

- Liquidus temperature is the temperature at which the first crystal appears on cooling (or the last disappears on heating) in a binary/ternary mixture of defined composition. ASTM D5972 freezing point and D2386 are measured on real fuel in a defined protocol with stirring; they are close to the liquidus but not identical (supercooling, wax haze versus first crystals).
- Eutectic rows are not liquidus points of a free composition; exclude them from any freezing-margin use. Two eutectic rows have no unique paired composition (`x_liquid` empty).
- The data are dominated by n-alkane/aromatic or n-alkane/cycloalkane binaries; many involve C20+ waxes (e.g. hexatriacontane, tetratetracontane) far above jet range, and several naphthalene-type PAH binaries. Only 15 systems reach T <= 253.15 K.
- `solid_phase_desc` says which component the authors report as crystallising; for mixtures with solid solutions or several Crystal phases all are listed. The SLE column `T_K` is the measured temperature, the composition columns refer to the liquid.
- The CSV has system-level repeats (one row per T/x point); use `system_id` as the group key for splits.
- Anything in the archives but outside the fixed property whitelist (e.g. `Equilibrium temperature of phase transition`) is not extracted.
