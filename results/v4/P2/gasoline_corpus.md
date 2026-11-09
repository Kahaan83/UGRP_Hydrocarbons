# Gasoline octane corpus (v4 Task 5)

Output: `data/real_fuels/v4/gasoline_octane_blends.csv` (210 rows).  Built by `scripts/v4/build_gasoline_octane.py`; raw inputs under `extras/fuelblend_data/data/raw/` (gitignored).

## Rows per source

| source | rows | with RON | with MON | composition_resolved | target_ok |
|---|---|---|---|---|---|
| cooptima_bron | 66 | 65 | 66 | 0 | 0 |
| face | 144 | 144 | 144 | 0 | 0 |
| **all** | 210 | 209 | 210 | 0 | 0 |

## Composition coverage

* FACE rows: 144; with an ethanol volume fraction parsed from the Name suffix: 57; all other FACE rows are a single unresolved lump (fraction 1.0).
* FACE rows whose hydrocarbon fractions are known: 0 (the CSVs contain component NAMES only, no fractions; see the loader docstring).
* FACE rows on a FACE A-J base lump: 24.
* bRON rows: 66 from 21 blendstocks; the blendstock fraction is a measured volume fraction in every row, the base gasoline (distinct base RON values: [87.7, 90.3]) is a lump.
* `target_ok` is False for every row: each blend contains a lump pseudo-component (R: lump RON/MON are model inputs only).  `base_lump_input_only` marks the bRON rows where the base RON/MON are recorded as inputs for a base-conditioned model (Task 15).

## Provenance

* `provenance` (ron/mon values): raw_file for all rows.
* `composition_provenance`: raw_file for bRON (numeric vol fractions in file); unknown for FACE (ethanol % inferred from the Name convention; hydrocarbon ratios absent).

## Skipped / not ingested

* Nothing skipped inside the two files (every data row yields a row).
* `TPRF 4+N` (N=1..15) not interpreted as ethanol percentages (series index): single lump.
* Co-Optima `fuel_engine_co_optimization.csv` Blend_RON1..3/percent columns and `bob_data.xlsx` are NOT ingested here (outside Task 5 scope).
* v3's typed FACE A-J base RON/MON and compositions (`ingest_phase1_expansion.py`) are NOT used: 47 of them conflict with the raw file (P1 audit) and no raw GC composition exists.
* Diisobutylene (CAS 107-39-1) carries the InChIKey of the `chemicals` database entry; the commercial material is an isomer mixture.  'Ketone Mixture' and 'Furan Mixture' (no CAS) are blendstock lumps.
