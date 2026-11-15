# delta_coutinho_lookup_audit: which aromatics the frozen Coutinho table gives another molecule's Tm / dHfus

Audit only, 2026-11-15. No code, constant, margin or result is changed. `src/fuelblend/thermo/coutinho_wax.get_pure_thermal_properties` is
frozen (see `composition_route.build_tmap_verified` docstring); this note records the full extent of its known SMILES-substring defect,
which was found case by case (toluene CLM-65, hexylbenzene delta_smoke_branch_opt).

The aromatic branch tests `"c1ccccc1" in smiles`, `"Cc1ccc(C)cc1" in smiles`, `"c1ccc2ccccc2c1" in smiles` and the durene pattern, so any
SMILES that **contains** one of these strings gets benzene's (278.7 K, 9.87 kJ/mol), p-xylene's, naphthalene's or durene's constants.

Over the 1,000+ registry compounds in `data/unified/compounds.csv`: **112 aromatics** carry another molecule's constants
(`misassigned.csv`): 97 as benzene, 8 as p-xylene, 6 as naphthalene, 1 as durene. Six of them have verified values:

| compound | assigned as (Tm K) | verified Tm K | direction for a freezing maximum |
|---|---|---|---|
| 2-methylnaphthalene | naphthalene 353.4 | 307.73 | conservative (predicts too warm) |
| hexylbenzene | benzene 278.7 | 209.70 | conservative |
| toluene | benzene 278.7 | 178.15 | conservative |
| o-xylene | benzene 278.7 | 247.82 | conservative |
| **fluorene** | benzene 278.7 | **387.94** | **unsafe** (predicts too cold) |
| **phenanthrene** | benzene 278.7 | **372.38** | **unsafe** |

Where it matters:
* The designated freezing route (`tm_dhfus_verified.csv` only, `build_tmap_verified`) does not use this table. Every verified compound above
  is scored with its verified value there.
* The table is the `coutinho_table` fallback of `library_route` / `predict_route` for members without verified constants; such freezing
  constraints are `uncertified` by rule, so no certificate depends on it. Point predictions printed on that path (e.g. the P2/P3 freezing
  values in delta_smoke_branch_opt) carry the defect.
* Fix, if wanted (user decision, changes frozen BASE numbers): exact canonical-SMILES identity for the four special cases, and NaN (no
  constants) for every other aromatic.
