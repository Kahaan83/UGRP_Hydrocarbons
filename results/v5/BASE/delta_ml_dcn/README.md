# delta_ml_dcn: linear blend + learned pair residual for cetane number (DCN) - NEGATIVE FEASIBILITY, NOT TRAINED

Requirement source: user ("Yes build density and cetane number"). Pattern: `results/v5/BASE/delta_ml_viscosity/README.md`.
The feasibility gate was declared before looking at the data: if fewer than about 40 mixture systems have every component with a measured pure
DCN, stop and report, with no model. That gate failed. No model, no split, no readings were produced.

## Data (Leenhouts benchmark, Zenodo 19914760, CC BY; read from `extras/fuelblend_data/data/raw/zenodo/leenhouts_benchmark/{datasets,splits}.zip`,
extracted copy at `data/interim/leenhouts_dcn/`, gitignored)
- `dcn_pure_exp.csv`: 229 pure rows (229 distinct InChIKeys, no repeats), 80 of them hydrocarbon-only (all atoms C/H).
- `dcn_mix_exp.csv`: 484 rows, all molar fractions (each row sums to 1), 2 to 4 components (265 binary, 134 ternary, 85 quaternary).
- Hydrocarbon-only mixture rows: 346 rows, 38 systems (system_id = md5 of the sorted InChIKey set).
- Hydrocarbon-only mixture systems with every component having a measured hydrocarbon pure DCN: **31 systems, 324 rows**
  (n_comp 2: 16 systems / 202 rows, n_comp 3: 6 / 71, n_comp 4: 9 / 51; 7 hydrocarbon systems / 22 rows have a component with no pure DCN).
- Mixture components: 25 distinct hydrocarbons, 21 with a measured pure DCN.
- Even without the hydrocarbon filter (all 57 mixture systems incl. oxygenates): 36 systems with all pure DCN known, still below 40.
- Leenhouts also ships a 5-fold mixture-grouped DCN split (`splits/fuel_ignition_numbers/dcn_mix_exp/split_definitions/mixture/cross_validation_5fold`);
  it was not used because the gate stopped the work.

## Declared plan (unexecuted)
Baseline: linear blend of pure DCN by mole fraction (no densities in this table, so a volume-fraction variant is not available from this source).
Learned form: DCN = DCN_linear + sum_{i<j} x_i x_j g(pair features), LightGBM with the viscosity model's fixed hyperparameters, seed 42.
Split: deterministic locked 20 % of systems by md5(system_id) hash. Readings R-MLC-1 to R-MLC-4 as in the task brief.

## Why it is stopped (reading)
31 systems is below the 40-system floor, and a 20 % locked split would leave about 6 locked systems, below the 9 systems that a finite conformal
margin at alpha 0.10 needs. A pair-residual model with 25 components and 31 systems (the viscosity model had 264 dev systems) cannot be
validated honestly, so nothing is trained. Scored: NEGATIVE FEASIBILITY (R-MLC-1 to R-MLC-4 not run, data).

## What would unblock it
More cetane mixture data with pure-component DCN for each component (the NREL cetane compendium is the Leenhouts source; ASTM D6890 mixture
measurements for GC-characterised fuels), or measured DCN for the 4 components that lack a pure value (this would add 7 systems, giving 38).
No `src/`, `scripts/`, `models/`, or tests were added for DCN.
