# External Datasets: Parsers & Staged Ingestion Suite

> **Isolation Notice:** These scripts and datasets are strictly isolated from the active production pipeline (`data/unified/`, `data/real_fuels/v4/`, `scripts/v4/`, `scripts/v5/`). They download, parse, and stage experimental measurements in `data/interim/external_staging/` without modifying or triggering any active production pipeline workflows.

---

## 1. Overview of Parsers & Ingestion Scripts

| Script | Source & Reference | Raw Data Input | Staged Output in `data/interim/external_staging/` | Records |
|---|---|---|---|---|
| `parse_sjtu_ai4fuel.py`<br/>`ingest_sjtu_ai4fuel.py` | **SJTU AI4fuel / LG-Transformer**<br/>Zhang, Huang et al. (*Nat. Commun.* 2026)<br/>DOI: `10.1038/s41467-026-73853-z`<br/>Figshare: `10.6084/m9.figshare.29571449` | `data/raw/git/AI4fuel/all_data.xlsx` | `canonical_ai4fuel_measurements.csv`<br/>`canonical_ai4fuel_compounds.csv`<br/>`ai4fuel/*.csv` (15 sheets) | **7,709 measurements**<br/>**1,763 unique compounds**<br/>(RON 379, MON 344, CN 637, DCN 178, FP 1054, LFL 813, MP 1068, Dens 1187, Visc 308, VP 692, LHV 1049) |
| `parse_kaust_jameel.py`<br/>`ingest_kaust_jameel.py` | **KAUST CCRC / KFUPM**<br/>Abdul Jameel et al. (*Energy & Fuels* 2018)<br/>DOI: `10.1021/acs.energyfuels.8b00556`<br/>Handle: `10754/627695` | `data/raw/kaust_jameel2018/acs.energyfuels.8b00556.pdf` | `kaust_jameel_canonical_measurements.csv`<br/>`kaust_jameel_table3_pure.csv` (128 pure)<br/>`kaust_jameel_table4_blends.csv` (123 blends)<br/>`kaust_jameel_table5_face.csv` (30 FACE) | **124 canonical measurements** (62 RON, 62 MON; 82 multi, 28 binary, 14 pure) |
| `parse_crc_face1_diesel.py`<br/>`ingest_crc_face1_diesel.py` | **CRC Report FACE-1**<br/>Cannella et al. (CRC AVFL Committee, 2014)<br/>9 FACE Research Diesel Fuels (FD1A–FD9A) | `data/raw/diesel/CRC_FACE-1.pdf` | `canonical_crc_face1_diesel_measurements.csv`<br/>`crc_face1_diesel_properties.csv`<br/>`crc_face1_diesel_composition_table52.csv` | **63 real diesel measurements**<br/>(9 Cetane D613, 9 DCN D6890, 9 Visc at 40°C, 9 Cloud point, 9 Flash point, 9 Density, 9 Net HOC) |
| `parse_ternary_flash.py`<br/>`ingest_ternary_flash.py` | **Liaw et al. & Li et al.**<br/>*Korean J. Chem. Eng.* & *Adv. Mater. Sci. Eng.* (2014)<br/>DOI: `10.1155/2014/239023` | `data/raw/ternary_flash/liaw_ternary_nonane_decane_dodecane.pdf`<br/>`data/raw/ternary_flash/li2014_ternary_alkanes.pdf` | `canonical_ternary_flash_measurements.csv`<br/>`ternary_flash_liaw_korean.csv`<br/>`ternary_flash_li2014.csv` | **24 canonical flash records**<br/>(17 multi-component ternary, 7 pure calibration) |
| `ingest_leenhouts_mon_multicomponent.py` | **Leenhouts / Kuzhagaliyeva et al.**<br/>*Commun. Chem.* (2022)<br/>Zenodo `19914760` | `data/interim/leenhouts_dcn/.../mon_mix_exp.csv` | `canonical_leenhouts_mon_mix_measurements.csv`<br/>`canonical_leenhouts_mon_compounds.csv` | **292 multi-component blend rows**<br/>(up to 80 components per blend across 67 unique systems) |
| `parse_boehm2022_cold_viscosity.py`<br/>`ingest_boehm2022_cold_viscosity.py` | **Boehm et al. HEAT Lab**<br/>*Front. Energy Res.* (2022)<br/>DOI: `10.3389/fenrg.2022.1074699` | `data/raw/boehm2022/HEAT_LAB_ViscosityData_2022.xlsx` | `canonical_boehm2022_cold_viscosity.csv`<br/>`boehm2022_raw_extracted.csv` | **2,162 cold-T measurements**<br/>(1,081 kinematic viscosity + 1,081 density from -40 °C to 15 °C) |

---

## 2. Summary of Staged Canonical Data

Across all external datasets parsed, a total of **10,374 experimental measurement rows** have been validated, converted to standard SI units, mapped to InChIKeys and chemical structures, and staged in `data/interim/external_staging/`.

```
data/interim/external_staging/
├── ai4fuel/                                      # Extracted per-property CSVs (15 sheets)
├── canonical_ai4fuel_measurements.csv            # 7,709 pure property measurements
├── canonical_ai4fuel_compounds.csv               # 1,763 unique compound structures
├── canonical_boehm2022_cold_viscosity.csv        # 2,162 cold viscosity & density rows
├── canonical_crc_face1_diesel_measurements.csv   # 63 real diesel property measurements
├── canonical_leenhouts_mon_mix_measurements.csv  # 292 multi-component MON blend rows
├── canonical_leenhouts_mon_compounds.csv         # 131 unique blend components
├── canonical_ternary_flash_measurements.csv      # 24 ternary & binary flash points
└── kaust_jameel_canonical_measurements.csv       # 124 pure & multi-component octane rows
```

All staged measurement files strictly follow the canonical `data/unified/measurements.csv` schema:
`source, property, value, T_K, P_kPa, kind, n_comp, is_hc, system_id, comp_inchikeys, x_mole, w_mass, fuel_id, ref, raw_value, raw_unit`.
