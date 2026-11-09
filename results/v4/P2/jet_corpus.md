# Jet corpus (v4 Task 6)

Built by `scripts/v4/build_jet_corpus.py` (loader `src/fuelblend/data/loaders_jet.py`).  Outputs: `data/real_fuels/v4/jet_class_compositions.csv` (7 fuels), `data/real_fuels/v4/jet_properties_v4.csv` (2768 rows).

## Fuels (class x carbon-number compositions)

* fuels: 7; composition_complete: 5; excluded_partial_composition: 2
* **headline_ok: 1**; not headline: 6 (assumed_cn_split: 4, partial composition: 2)
* link to inventory_v4: none of the 7 GC fuels is in inventory_v4 (those POSF ids were quarantined in Task 3); ids are disjoint (asserted), `inventory_fuel_id` is empty.  The 57 `CRC647_*` inventory fuels are linked on the property side only.

| fuel_id | complete | cn_split_method | headline_ok | pseudo-components | target props | reason |
|---|---|---|---|---|---|---|
| JetA_API14avg | True | assumed_cn_split | False | 6 | 0 | assumed_cn_split: ['aromatic', 'cycloalkane', 'paraffin_unsplit'] |
| JP5_PotterSimmons | False | source | False | 27 | 0 | partial_composition (listed 0.409 of the fuel; unresolved 0.591) |
| JP8_PotterSimmons | False | source | False | 29 | 0 | partial_composition (listed 0.423 of the fuel; unresolved 0.577) |
| UDRI_HTA_repfuel | True | assumed_cn_split | False | 20 | 0 | assumed_cn_split: ['aromatic', 'cycloalkane', 'n_paraffin'] |
| POSF-10325 | True | D2887_same_fuel | True | 48 | 3 |  |
| POSF-13511 | True | assumed_cn_split | False | 4 | 0 | assumed_cn_split: ['aromatic', 'cycloalkane', 'iso_paraffin', 'n_paraffin'] |
| POSF-13516 | True | assumed_cn_split | False | 3 | 0 | assumed_cn_split: ['cycloalkane', 'iso_paraffin', 'n_paraffin'] |

## Properties present per source (raw files only)

| source | fuels | rows | target rows | targets present (fuels with it) |
|---|---|---|---|---|
| av2718_fuel_properties.csv | 1 | 33 | 11 | density_kg_m3 (1), flash_point_C (1), viscosity_kin_cSt (1) |
| crc647_fuel_properties.csv | 57 | 2735 | 797 | density_kg_m3 (57), flash_point_C (57), freezing_point_C (57), net_heat_of_combustion_MJ_kg (57), pour_point_C (57), viscosity_kin_cSt (57) |

* inventory_v4 CRC647 fuels with properties: 57 of 57
* `data/MISSING.md` rows added (block `v4-task6`): 7

## Notes / limits

* The raw GC file lists class totals AND sub-totals (POSF-10325 raw rows sum to 152.28 %); the loader sums only a named non-overlapping partition per fuel (POSF-10325: 0.9999).
* POSF-10325's carbon-number split comes from the same fuel's av2718 ASTM **D2887** simulated distillation (no D86 exists for it in a raw file) through `fuel_hindcast.carbon_distribution`; that function applies ONE carbon-number distribution to every class (class-independent), which is an approximation.  Column `cn_split_method` records it.
* Property values were never interpolated; nothing is NaN-filled (absence = no row).  Btu/lb -> MJ/kg is the exact factor 0.002326; the density `(60F)` unit label gives T = 15.56 C.
