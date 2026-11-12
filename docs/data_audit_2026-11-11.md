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

## 5b. Second online pass (2026-11-12)

Status words as in section 5. Nothing in this section has been downloaded into the repo.

### NIST ThermoML: files on the server that are not in the bulk archive (largest find of this pass)

The bulk archive `ThermoML.v2020-09-30.tgz` and NIST's own browse index (`/ThermoML/jsons/journal_issue.json`, 11,922 DOIs) stop at
J. Chem. Eng. Data 2019. The server still answers `https://trc.nist.gov/ThermoML/<doi>.json` and `.xml` for later articles.
Every J. Chem. Eng. Data DOI of the forms `acs.jced.9b00001-9b01250` and `acs.jced.0c00001-0c01150` that is absent from the index was
checked for existence (2,185 checks): **367 unlisted files exist** (179 with `9b` DOIs, 188 with `0c` DOIs, the highest being 0c00724; all those opened are volume 65, 2020).
Each file's metadata was read; 52 hold at least one hydrocarbon-only data set. The ones that matter here:

| DOI (10.1021/acs.jced.) | first author, J. Chem. Eng. Data 65 (2020) | hydrocarbon-only content | T range, K |
|---|---|---|---|
| 0c00020 | Luning Prak, p 2625 | iso-dodecane (2,2,4,6,6-pentamethylheptane) + n-propylbenzene, 1,3,5- and 1,2,4-trimethylbenzene: 198 binary viscosity, 264 density, 231 speed of sound | 253.15-333.15 |
| 0c00280 | Luning Prak, p 3941 | n-dodecane + five C10 aromatics: 220 binary viscosity, 308 density | 288.15-333.15 |
| 0c00581 | Luning Prak, p 4699 | iso-dodecane + butylbenzene isomers and 1,3-diethylbenzene: 185 binary kinematic viscosity, 259 density | 288.15-333.15 |
| 0c00309 | Lei, p 3977 | isopropylcyclohexane + n-tridecane: 99 binary viscosity, 99 density | 293.15-333.15 |
| 9b01171 | Zhao, p 2527 | exo-tetrahydrodicyclopentadiene + n-decane: 99 binary viscosity, 99 density | 293.15-333.15 |
| 9b01167 | Qin, p 2512 | three alkyladamantanes + n-nonane or n-undecane: 726 binary viscosity, 726 density | 293.15-343.15 |
| 9b01026 | Gahlyan, p 2964 | benzene + n-heptane + n-octane: 108 ternary density | 298.15-318.15 |
| 9b01208 | Touronen, p 3270 | vapour pressure of toluene, hexylbenzene, 2-ethylnaphthalene: 40 points | 400-620 |
| 0c00111 | Troni, p 3334 | n-undecane, n-dodecane, n-tridecane boiling temperatures (pure and binary) | 377-449 |

Counts are data-set sizes read from the JSON (`NumValues`), not rows passed through our extractor.
Why it matters: about 1,500 new hydrocarbon binary viscosity points and 1,850 density points, mostly iso-paraffin + aromatic and
cycloparaffin + n-paraffin pairs, which looked thin when compared with the ThermoML archive alone (see the correction below; the earlier reference here to the S7 pairwise model was wrong, S7 re-analyses real-fuel predictions and has no lab-pair training). 0c00020 reaches -20 C, the temperature
of the jet viscosity limit, where the lab margin is widest (S1).
Limits: (a) ThermoML has no flash-point property, so the flash points these papers report are still not available this way;
(b) the other four journals (Elsevier and Springer DOIs) cannot be enumerated the same way, so their 2020 files, if any, are unchecked;
(c) the search API (`/ThermoML-API`) returned HTTP 503 on both days, so a query by compound was not possible;
(d) no J. Chem. Eng. Data 2021+ file was found (75 spot checks of `1c`, `2c`, `3c` DOIs, all absent);
(e) these files carry the same NIST liability statement as the archive and are not listed by NIST, so they may be less checked.
Total size of the nine JSON files is about 4 MB.

**Update, same day:** the nine XML files were downloaded with the user's go-ahead and run through the unchanged extractor:
3,768 hydrocarbon rows, 0 file errors; 1,527 mixture viscosity rows in 19 systems, 17 of them absent from the 2019 archive.
See `results/v5/thermoml_2020/README.md`.

**Correction, 2026-11-12 (after the S1b run):** "17 systems absent from the 2019 archive" compared against the ThermoML archive only.
Against the project's unified data the files add much less viscosity: 15 of the 19 mixture viscosity systems are already in
`measurements_dedup.csv` through NIST TDE and Leenhouts. They are the same measurements: 738 rows match a NIST TDE row on system, composition and temperature, all within 0.5 % in value. After end points are removed and repeats merged
there are 1,278 mixture viscosity rows: 756 sit at a (system, temperature) the unified data already has, 342 are new temperatures of
known systems (including all 27 rows at 253 K), and 180 belong to 4 systems that are new (iso-dodecane with butylbenzene isomers or
1,3-diethylbenzene). Density is mostly new: 1,440 of 1,599 mixture density rows are in systems found neither in the unified data nor in `data/processed/v2/thermoml_clean_hc.csv`.
All 27 rows at 253 K belong to three systems that are locked in the S1 split, so they cannot enter a dev margin.
The "about 1,500 new viscosity points" in the paragraph above is therefore an overstatement; read it as 522 new rows.
Used once as a held-out test (`results/v5/S1b/`): only 3 independent scorable systems, so every coverage reading is "not evaluable".
Lesson for the next source: check overlap against `measurements_dedup.csv`, not against one upstream archive.

### Other sources checked in this pass

| # | source | fills | status |
|---|---|---|---|
| 8 | NIST public data mds2-2289 (Bell), https://doi.org/10.18434/mds2-2289, file `MIDAS_alkanes.csv` | Literature viscosity of the n-alkanes from NIST ThermoData Engine 10.4.2, with a reference code per point | CONFIRMED record and file list (catalog.data.gov). Carbon-number range, point count and T range not stated on the record; the CSV was not opened. NIST open licence. Pure fluids only. |
| 9 | Puetz et al., Energy & Fuels 39 (2025) 18001, doi 10.1021/acs.energyfuels.5c02762 (open copy at elib.dlr.de/216631) | 8 samples (Jet A-1, ATJ-SPK, two FT-SPK, IsoPar G/H/J/L) with exact GCxGC composition in Supporting Information E and measured density, viscosity, distillation, heat of combustion, cetane | CONFIRMED from the article text. The Supporting Information is on the ACS site, which blocks scripted access; needs a browser download. Small, but the four IsoPar solvents are nearly pure iso-paraffin, where our representatives are weakest. |
| 10 | Mueller et al., Energy & Fuels 30 (2016) 1445, doi 10.1021/acs.energyfuels.5b02879 | Measured properties of the four CRC diesel surrogates (section 5, item 2) | Citation now CONFIRMED (PNNL, Sandia, OSTI records). No free full text found: OSTI record 1237000 has no file, PMC has no copy. Needs library access. |
| 11 | Li et al. 2014, doi 10.1155/2014/239023 (section 5, item 5) | Binary and ternary n-alkane flash points | Still not opened: Hindawi and Wiley both refuse scripted access (HTTP 403). Needs a browser download. |

Checked and closed:
* CRC AV-25-16 (2021, 93 airport jet fuel samples): 39 pages opened. Composition is class totals only (Appendix II, one page); the rest is trace
  heteroatoms, microbes and ions. No flash, freezing or viscosity tables. Not useful here.
* Luedtke et al., Energy & Fuels 39 (2025) 4434: a composition-method paper; its supplement has no property data.
* arXiv 2408.01530 (FTIR and machine learning for jet fuel properties): its data are AJFTD, which we already hold.
* No open diesel dataset with GCxGC composition and flash, viscosity and cloud point was found beyond CRC FACE-1 (already used in S6).
* No open table of pure iso-paraffin or alkylcyclohexane viscosity at -40 C was found.

## 6. Housekeeping noticed

* `AGENTS.md` section 8.1.3 still says AJFTD has 81 fuels; the v5 scrape has 2,601 in the metrics file and 2,657 in the detail file.
* `docs/ajftd_v5_note.md` said no composition exists; a correction line was added.
* `extras/fuelblend_data/data/unified/BASELINE_VISCOSEITY_REPORT.md` is a 2-byte misspelled duplicate.
* `extras/ajftd_session.json` and `extras/fuelblend_data/ajftd_session.json` hold login cookies. `extras/` is gitignored; keep it so.
* `docs/data_sources.csv` lists `fill_lfl.py` as a user of `lfl_chemicals.csv`; it never read that file.
