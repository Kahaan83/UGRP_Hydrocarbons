# Data audit, 2026-11-11: what is on disk, what is unused, and what was recovered

Scope: `data/` (repo root) and `extras/fuelblend_data/data/`. "Used" means read by the current pipeline (`src/fuelblend`,
`scripts/v4`, `scripts/v5`, `tests`). The method was a file-by-file scan of 331 data files against that code, followed by
opening each unused file. Counts below were computed on 2026-11-11 from the files named.

## 1. How the two folders relate

* `data/` is what v4/v5 reads. `extras/fuelblend_data/data/` is the raw staging area (gitignored) plus the v1 outputs.
* `data/unified/measurements_dedup.csv` is identical in both trees (76,307 rows). The root `data/processed/*.csv` files are
  copies of the extras ones.
* 33,636 of those rows are real-fuel rows with no composition (`data/processed/excluded_by_R1.csv`). Rule R1 keeps them out
  of mixture training. That exclusion is correct; they are used as bulk descriptors (v5 S4).

## 2. Recovered in this session

### 2.1 AJFTD GCxGC compositions (largest find)

* File: `extras/fuelblend_data/data/raw/ajftd/ajftd_fuel_details.csv` (162,943 cells, 2,657 fuels). No v4/v5 code read it.
* Why it was unused: `docs/ajftd_v5_note.md` concluded "no hydrocarbon-type composition suitable for blend modelling". That
  conclusion came from `ajftd_all_metrics.csv` only. The detail file holds the page's composition tree: weight % by
  hydrocarbon class and carbon number (source label on the page: "AFRL GCxGC").
* Fix: `src/fuelblend/data/loaders_ajftd_gcxgc.py`, outputs in `data/real_fuels/v5/ajftd_gcxgc_{composition,fuels}.csv` and
  `ajftd_gcxgc_build_log.json`; tests in `tests/test_v5_data_audit.py`.
* Result: 116 fuels have bins; 113 have a complete composition (leaf sum 98 to 102 wt %; median 99.85); 6,785 bins.
  Of the 113, the AJFTD bulk table has a measured flash point for 37, freezing point for 26, kinematic viscosity for 32,
  density for 30 and net heat of combustion for 36. 38 have at least one of these; 21 have flash, freezing, viscosity and
  density together. 18 of the 113 have under 1 wt % aromatics (synthetic blend stocks).
* Fuel labels are POSF numbers. Four fuels that v4 dropped for "no clean composition" now have one: POSF 10289, 10325,
  12344, 12345. Three overlap FuelLib (4658, 10289, 10325) and three overlap CRC AV-27-18 (10289, 10325, 11769): these must
  share a double-count group before any scoring.
* AJFTD fuel 2 (label 10325) carries this AFRL composition with 18.5 wt % aromatics, which supports the open identity
  question (v4 item M7) but does not settle the 48 C vs 55 C flash conflict.
* Three fuels fail closure (labels 13350: 97.9 %, 9968: 78.9 %, 3941: 110.9 %) and are flagged, not repaired.
  74 complete compositions have no scoreable property row in the AJFTD bulk table.
* Scored the same day (`results/v5/G1/`, readings declared in its README before the run; nothing fitted). Bins were mapped
  to r3 representative compounds; 112 of 113 fuels keep at least 90 % of their mass (median unmapped 1.4 %), and 39 have a
  measured target. Error is predicted minus measured.

| property | model | n fuels | bias | MAE | note |
|---|---|---|---|---|---|
| flash point | ideal Le Chatelier | 37 | -0.43 C | 3.54 C | 11 % over-predicted by more than 5 C; worst +17.2 C |
| flash point | UNIFAC | 37 | -0.54 C | 3.77 C | no better than ideal |
| density, 15 C | ideal volume | 30 | -6.1 kg/m3 | 8.1 kg/m3 | biased low, as on CRC 647 |
| viscosity, 40 C | Refutas | 18 | -0.2 % | 5.3 % | |
| viscosity, -20 C | Refutas | 32 | -26 % | 34 % (ln 0.342) | under-predicts |
| viscosity, -40 C | Refutas | 12 | -49 % | ln 0.683 | under-predicts |
| freezing point | ideal SvL | 31 | +0.4 C | 15.2 C | constant mean gives 7.3 C: no skill |

  Readings: "measured composition removes the unsafe flash bias" holds (bias at most +2 C and at most 20 % over by 5 C);
  "UNIFAC helps" does not hold; "the freezing model has skill" does not hold. Four fuels over-predict flash by 8 to 17 C
  (labels 7695, 9324, 8092, 8438, each with one or two D93 tests); they are not explained. The 11 fuels in the S4 locked
  slice were read once here (flash bias -2.1 C, MAE 2.8 C, none over by 5 C).
* AJFTD is login-gated. The new CSVs sit beside the existing AJFTD v5 tables; do not publish raw rows.

### 2.2 Freezing (SLE) systems blocked by missing constants

* File: `data/processed/v4/thermoml_sle_hc.csv`, 59 scored systems, 778 points. v4 scored 8 systems (103 points), so the
  freezing margin was `inf` (9 systems needed).
* Why: the model needs a `verified` melting point and enthalpy of fusion for every component, and only design compounds had
  been audited. The data were never the limit.
* Fix: `scripts/v5/eval_sle_extended_v5.py` (same model and margin rule, wider constants in declared tiers), outputs in
  `results/v5/SLE2/`. The first tier reproduces v4 exactly.

| constants tier (cumulative) | systems | points | bias K | MAE K | 90 % margin K |
|---|---|---|---|---|---|
| v4 verified | 8 | 103 | -1.10 | 2.02 | inf |
| + n-alkanes confirmed by the Broadhurst table | 11 | 136 | -1.41 | 2.13 | 12.34 |
| + CRC-table enthalpy with no second source | 38 | 501 | -4.02 | 7.03 | 24.25 |

* Reading: a finite margin now exists, but it is 12 K at best. The v5 jet designs sit 0.5 to 0.7 C inside the freezing limit,
  so this does not certify them. The wide tier is a diagnostic: its errors are dominated by aromatic pairs and long n-alkanes,
  where the ideal model lacks solid-solution and rotator-transition terms.
* Still blocked: 21 systems, 277 points. Nine components have no tabulated enthalpy of fusion (seven dimethylnaphthalenes,
  1-hexyne, bicyclohexyl), three have disputed melting points (1-methylnaphthalene, tetracosane, octacosane), one is
  stereo-unspecified decalin.
* `configs/` was not changed.

### 2.3 Tabulated LFL values that were never merged

* Files: `extras/.../raw/lfl/lfl_chemicals.csv` (99 IEC 60079-20-1 values) and Co-Optima `Pure_LFL_LEL` (106 values, 39
  hydrocarbons). `lfl.csv` used 25 measured values and the Jones estimate for the other 670.
* Fix: `scripts/v5/build_lfl_overlay_v5.py` writes `data/processed/v5/lfl_v5.csv` with priority Zabetakis measured, IEC
  table, Co-Optima compiled, Jones estimate. 91 estimates are replaced (71 IEC, 20 Co-Optima); mean absolute change 6.5 %,
  largest 23 %. IEC and Co-Optima agree to 5 % on the 13 compounds both list.
* In the v4 flash evaluation 82 of 98 pure inputs use the estimate; 22 of those have an IEC value.
* The Co-Optima values are compiled; their primary source is not stated per row.
* Re-run on the overlay the same day (`results/v5/L1/`, readings declared first; overlay applied in memory, so
  `results/v4/` and `configs/` are unchanged). 25 of 79 pure inputs change in the v4 run and 10 of 45 in the AJFTD run.

| set (ideal model) | n | MAE before | MAE after | bias before | bias after | over by >5 C before | after |
|---|---|---|---|---|---|---|---|
| lab mixtures | 126 | 3.74 | 3.69 | -3.30 | -3.25 | 0.008 | 0.008 |
| CRC 647 D93, dev | 45 | 4.04 | 4.09 | +0.44 | -0.73 | 0.133 | 0.044 |
| CRC 647 D93, locked (shown only) | 12 | 3.72 | 3.26 | +2.46 | +1.41 | 0.250 | 0.083 |
| AJFTD GCxGC | 37 | 3.54 | 3.48 | -0.43 | -0.28 | 0.108 | 0.108 |

  Readings: "the overlay is safe to adopt" holds (no set has more over-prediction, no MAE rises by more than 0.2 C);
  "the overlay improves accuracy" does not hold (MAE moves by less than 0.1 C everywhere). Its effect is on the safety side:
  CRC 647 over-predictions fall from 6 to 2 of 45 dev fuels. Adopting it as the default input is still your decision.

## 3. Unused data that can be wired next (checked, not yet done)

| data | what it holds | why unused | what it takes |
|---|---|---|---|
| `extras/.../processed/nist_density_compressed_liquid.csv` (NISTIR 6659) | 634 liquid densities, 270 to 470 K, 0.5 to 40 MPa, for Jet A 3602/3638/4658, JP-8 3773, S-8 | no composition was linked | POSF 4658 has a GCxGC composition in FuelLib and in AJFTD (section 2.1) but no property rows: link its 121 density rows. The other four need AJFTD labels checked. |
| `extras/.../raw/rp1/rp1_blends_properties.csv` (NISTIR 6646) | 1,465 rows for the RP-1 "original sample": 486 densities, 34 kinematic viscosities, thermal conductivity | v4 names the inventory fuel "SRM 4572" and treats the property sample as a different, unconfirmed fuel | composition (p. 36) and properties come from the same report and sample, so they can be linked under a NISTIR 6646 id. Limit: the composition lists 80.6 wt %, 65.6 % with a structure. |
| `data/processed/v2/wax_dauphin_coutinho.csv` | 57 measured wax-appearance temperatures: n-heptane with C23/C25/C28/C32/C36 (53 points) and 4 multicomponent mixtures | v4 SLE evaluation reads ThermoML only | add as a second SLE source using the section 2.2 constants. Appearance temperature is not the liquidus; keep the method label. |
| `extras/.../raw/vozka2019/` | freezing point of 25 Jet A blends, GCxGC class composition for 11 | class totals only, no carbon numbers | usable for a class-level freezing check; cannot enter the bin-based route. |
| Co-Optima `fuel_engine_co_optimization.csv` | measured pure flash point (259 rows, 102 hydrocarbons), melting point (336) | only octane and DVPE columns were ingested | cross-check the pure flash table; second source for melting points. |
| Leenhouts `splits.zip` (200 MB) | the benchmark's published train/test splits | v4 built its own grouped splits | lets our flash and viscosity numbers be compared on the benchmark's own splits. |

## 4. Unused and should stay unused

| data | reason |
|---|---|
| `data/quarantine/*`, `data/real_fuels/v3/` | formula-derived v3 values (v4 Task 3). No recovery possible without a raw source. |
| `extras/.../raw/pires2018/` class compositions | GC mass closure 48 to 180 % (stated in the paper). The 130 property rows are fine as bulk descriptors. |
| `NASA_CR_179472_19860018789.txt` | n-alkane fractions only, not a full composition, so no fuel passes the 0.95 gate. Useful only as an n-paraffin descriptor against pour and freezing point. |
| `chew2024/si_viscosity_dataset-predictions.csv`, `zenodo_model.zip`, Bilodeau `pretrained_models.zip` (models) | model outputs and weights, not measurements. |
| `extras/.../processed/tsi_table.csv` derived rows | circular (v4 Task 14). |
| DCN, MON, cetane tables | measured and clean, but not target properties. |
| Fortin 2012 distillation curves for the 11 RP-1 blends | recipes are confidential, so there is no composition to pair them with. |

## 5. New sources found online (2026-11-11)

Status words: CONFIRMED = document opened and the table seen; LISTED = found, contents not opened.

| # | source | fills | status |
|---|---|---|---|
| 1 | NIST Technical Note 1949 (Huber), https://doi.org/10.6028/NIST.TN.1949 | Diesel gap (SD-02). Table 1 gives exact mole fractions of the four CRC diesel surrogates V0a, V0b, V1, V2 (4 to 9 pure compounds); Table 11 gives measured kinematic viscosity at 40 C and surface tension. | CONFIRMED. Public NIST report. Only 4 fuels, 1 viscosity each. |
| 2 | Mueller et al., Energy & Fuels 30 (2016), the paper NIST TN 1949 cites for the measured values (CRC AVFL-18a) | Same four surrogates: density, cloud point, flash point, distillation, cetane | NOT OPENED. Citation recalled from memory; verify the DOI before relying on it. |
| 3 | CRC Report FACE-1 (2010), https://crcao.org/wp-content/uploads/2019/06/FACE-1-Chem-and-Phys-Props-of-FACE-Research-Diesel-Fuels.pdf | Nine FACE diesel fuels: 2-D GC group types (Table 5.2), PIONA and SOAP tables (7.1, 8.1 to 8.3), DHA (10.1), plus ASTM properties. Covers the 9 dropped `CRC_FACE_Diesel_*` fuels. | CONFIRMED that the tables exist (242 pages). Not yet checked whether composition is per carbon number. |
| 4 | CRC Report AVFL-24 (2014), https://crcao.org/wp-content/uploads/2019/05/AVFL-24-FACE-Gasolines-Report-071414.pdf | FACE gasolines and ethanol blends (v4 item M9; 34 dropped fuels) | LISTED. 60-page report opened to its title page only. |
| 5 | Li, Duan, Sun, Yan, Adv. Mater. Sci. Eng. 2014, 239023, doi 10.1155/2014/239023 | Measured flash points of binary and ternary n-alkane mixtures | LISTED (abstract). Systems and counts not confirmed. Hindawi journals are normally open access. |
| 6 | Hwang and In, Fire Sci. Eng. 31(3) 19 (2017), doi 10.7731/KIFSE.2017.31.3.019 | Flash points of toluene, methylcyclohexane, n-heptane each with ethylbenzene (Seta closed cup) | CONFIRMED abstract; values are in Table 4 as an image, so they need typing. Binary only. |
| 7 | Acree and Chickos fusion-enthalpy compendia (J. Phys. Chem. Ref. Data 2016/2017; Thermochim. Acta 2009 supplement) | The 9 missing enthalpies of fusion in section 2.2 | LISTED. Access terms restrict reuse of the JPCRD parts. |

Not found: an open Zenodo-style dataset of jet fuels with GCxGC composition and measured properties. The AJFTD detail file
(section 2.1) is that dataset, and it was already on disk.

## 6. Housekeeping noticed

* `AGENTS.md` section 8.1.3 still says AJFTD has 81 fuels; the v5 scrape has 2,601 in the metrics file and 2,657 in the detail file.
* `docs/ajftd_v5_note.md` said no composition exists; a correction line was added.
* `extras/fuelblend_data/data/unified/BASELINE_VISCOSEITY_REPORT.md` is a 2-byte misspelled duplicate.
* `extras/ajftd_session.json` and `extras/fuelblend_data/ajftd_session.json` hold login cookies. `extras/` is gitignored; keep it so.
* `docs/data_sources.csv` lists `fill_lfl.py` as a user of `lfl_chemicals.csv`; it never read that file.
