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
