# delta_fresh_holdout (Task 15)

Declared 2026-11-12 before `delta_fresh_holdout.py` exists; this commit contains no result file.
Frozen model version: scoreboard commit `5a62e20` (designated ideal-volume density on the mapped measured composition, raw provider densities).

Steps (binding order): (1) this declaration. (2) Load the NISTIR 6659 / Outcalt et al. 2009 POSF 4658 rows (121) from
`extras/fuelblend_data/data/processed/nist_density_compressed_liquid.csv` (skip-guarded), write `data/real_fuels/v5/nistir6659_density.csv`; keep P <= 1000 kPa
for the ideal-volume comparison (the model has no pressure term). (3) List the AJFTD GCxGC fuels with complete composition and no property row by POSF label, and
search only tables already under `data/` for rows carrying the same POSF label. (4) Check every found fuel against the double-count groups. (5) Score once;
`scored_once.json` then blocks any second run.

Honest scope: POSF 4658 is AJFTD:23, already a T2 fuel (its 288.15 K density 806.5 kg/m3 was read). The NIST rows are new temperatures for a read fuel, not a fresh fuel.
Reading (fixed now): "density bias on the fresh fuels lies inside the T2 descriptive 95 % interval", where the interval is the 2.5-97.5 % range of the
frozen T2 pooled density errors at 288.15 K (one fuel per group, from `rows_t2_bulk.csv`).

## Result 2026-11-12 (scored once; `scored_once.json` written)

* Step 2: 121 POSF 4658 rows loaded to `data/real_fuels/v5/nistir6659_density.csv`; 11 have P <= 1000 kPa (about 0.5 MPa, 270-470 K) and were scored.
* Step 3: 73 AJFTD GCxGC fuels have complete composition and no property row. A search of every CSV under `data/real_fuels`, `data/processed` and `data/unified`
  (AJFTD tables excluded) for their POSF labels found 0 rows. No fresh fuel exists on disk; U3 and U5 are what would add one.
* Step 4: AJFTD:23 is its own double-count group (no other member).
* Step 5: density bias -12.3 kg/m3 (MAE 12.3) over 11 temperatures, roughly constant from 270 to 470 K (-9.2 to -13.7). The T2 descriptive 95 % interval at 288.15 K is
  -14.2 to +2.4 kg/m3, so the reading PASSES, but only just: the bias sits near the 3rd percentile of the T2 errors, and this is one already-read fuel (its frozen 288.15 K
  error was about the same). This is not a held-out claim for the T2 population; n of fresh fuels = 0.
