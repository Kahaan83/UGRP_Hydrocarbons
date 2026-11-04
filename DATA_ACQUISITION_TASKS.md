# DATA_ACQUISITION_TASKS.md — task ledger for the coding harness

Written 2026-10-29. Companion to `PROJECT_HANDOFF.md` (read its §6 schema and §17 first).
Goal of this file: fill the data gaps for **vapour pressure / Antoine, flash point + LFL, octane (RON/MON),
smoke point, GCV, pour point, density** so the seven targets in `PROJECT_HANDOFF.md` §17 (viscosity, density, flash, pour, GCV, octane, smoke point) can be trained or honestly
reported as "baseline only".

All URLs below were found by web search on 2026-10-29. Open each one before relying on it: if it returns non-200,
an HTML consent page instead of a PDF, or a login wall, mark the task `MANUAL` and move on (do NOT work around it).

---

## 0. Rules (apply to every task)

1. **No paywall / login / ToS circumvention.** Items marked MANUAL are for the user to download (university library
   access may help: ASTM, CRC Handbook, Yaws, DIPPR, Knovel). Drop files in `data/raw/_manual_inbox/` and the loader
   picks them up.
2. **Rate limits.** ≤1 request/s to NIST WebBook, ≤3 req/s to PubChem, identify yourself in the User-Agent. Cache every
   response under `data/raw/<source>/` with sha256 in `MANIFEST.json` (reuse `fetch.py`).
3. **One loader per source** → canonical schema (`PROJECT_HANDOFF.md` §6.2) plus these NEW columns:
   `tier` (A/B/C), `method` (e.g. `closed_cup`, `open_cup`, `RON`, `MON`, `D1322`, `unspecified`), `derived` (bool),
   `domain` (`middle_distillate` | `gasoline` | `both`), `source_table` (table/page in the paper).
4. **Tiers.** A = evaluated or primary measurement (JPCRD, ThermoML, NIST WebBook with a cited experiment, IEC/NFPA, API-45).
   B = compiled database that cites experiments (Li et al., vom Lehn experimental rows, Jameel SI). C = estimated or derived
   (ANN/QSPR, YSI→smoke point, correlation-GC, Joback, interpolated Antoine, GCV from group additivity).
   **Never train on tier C.** Keep it, flagged, for fallback features and diagnostics only.
5. **New property names** (units): `ron`, `mon` (–), `gcv_MJ_kg` (MJ/kg, gross, liquid fuel, liquid water, 25 °C basis unless stated),
   `ncv_MJ_kg`, `lfl_vol_pct` (+ `T_K`), `tsi`, `ysi`, `smoke_point_mm`, `pour_point_C`, `cloud_point_C`, `freezing_point_C`,
   `vapor_pressure_kPa`, `density_kg_m3`, `tc_K`, `pc_kPa`, `omega`, `flash_point_C`. Unit conversion happens in the loader only.
6. **Every loader gets a test** with a fixture that mirrors the REAL file structure (lesson §9.9) and an audit check.
7. **After each task** append to `DATA_LOG.md`: task id, rows added by property × kind (pure/binary/multi/real_fuel),
   rows rejected and why, and any contradiction with existing data.
8. Update `build_unified.py` only after reading the user's current version (§16.1 warning). Do not overwrite from an old zip.

---

## 1. Global sanity checks (run after every build; these are loose, from handbook knowledge — a failure means look, not auto-reject)

| check | expected |
|---|---|
| benzene Psat at 40 °C | ≈ 24.3 kPa |
| n-decane Psat at 25 °C | ≈ 0.18 kPa |
| n-heptane / iso-octane | RON = MON = 0 / 100 (PRF scale definition) |
| n-octane GCV (liquid, 25 °C) | ≈ 47.9 MJ/kg (ΔcH ≈ −5470 kJ/mol ÷ 114.23 g/mol) |
| benzene GCV | ≈ 41.8 MJ/kg |
| n-hexane density at 25 °C | ≈ 655 kg/m³ |
| n-dodecane density at 25 °C | ≈ 745 kg/m³ |
| n-decane flash point | roughly 46–51 °C (cup type matters) |
| unit trap | ΔcH is negative in WebBook; GCV is stored positive |
| TSI scale | methylcyclohexane = 5 and 1-methylnaphthalene = 100 by construction of the Boehm scale; a source on another scale must be rescaled or kept separate |

---

## 2. Tasks

Legend: **AUTO** = harness can do it; **SEMI** = harness downloads, human checks tables; **MANUAL** = user fetches.
Priority: P1 = unblocks a target now; P2 = improves a target; P3 = nice to have.

### 2.0 T00 — source check (do this FIRST, before any loader)

| id | pri | mode | action | output |
|---|---|---|---|---|
| T00 | P1 | AUTO | For every URL, package function name and API heading in this file: request it, record HTTP status, content-type and size; for Python names run `dir()` on the installed `chemicals`/`thermo`; for PubChem confirm the exact heading strings; for git repos list the tree. Items stated in the handoff as "from background knowledge" (PubChem headings, GESTIS terms, KDB/CHERIC, Stull 1947, Morgan & Kobayashi 1994, Calcote & Manos TSI formula) must be confirmed or dropped. Do not download datasets in this task. | `SOURCE_CHECK.md`: one row per source with status OK / MOVED / BLOCKED / DIFFERENT, and what the user must do for each non-OK row |

### 2.1 Vapour pressure and Antoine constants (feeds flash point, density via COSTALD, and is a target itself)

Design rule: collect **P–T points**, not just fitted constants. Antoine constants are only valid in their stated range;
for heavy compounds the flash baseline needs P at ≈ 400–520 K (§16.2). Refit with `refit_antoine.py`.

| id | pri | mode | source | action | output |
|---|---|---|---|---|---|
| T01 | P1 | SEMI | Růžička & Majer 1994, JPCRD 23, 1 — recommended vapour-pressure equations for n-alkanes C5–C20 between triple and normal boiling points. Free PDF: `https://srd.nist.gov/jpcrdreprint/1.555942.pdf` | Extract the recommended equation parameters and the tabulated vapour pressures per alkane. Fills the C13/15/17/18/19 gap. Human-check against the PDF page for 3 compounds. | `data/raw/ruzicka_majer/` → `vapor_pressure_kPa`, tier A, pure |
| T02 | P1 | SEMI | Chickos & Lipkind, J. Chem. Eng. Data 2008, 53, 2432 — n-alkane vapour pressures, C21 and heavier, 298–540 K. Author-hosted PDF: `https://umsl.edu/~chickosj/JSCPUBS/08-2.pdf` | Extract the vapour-pressure equations/tables. Mark `derived=True`, tier C (correlation-GC), use only for heavy n-alkane coverage diagnostics. | `data/raw/chickos2008/` |
| T03 | P1 | AUTO | Python packages `chemicals` + `thermo` (open source, cited public data) | `pip install chemicals thermo`. For every CAS in the registry export: the vapour-pressure coefficient tables (`chemicals.vapor_pressure.Psat_data_AntoinePoling`, `Psat_data_AntoineExtended`, `Psat_data_Perrys2_8`, `Psat_data_WagnerPoling`; Yaws coefficients and Landolt Antoine coefficients are reachable through `thermo`'s vapour-pressure methods — confirm names with `dir()`), and the critical constants (Tc, Pc, Vc) and acentric factor. **Normalise units** (these tables are not all kPa/K/log10; test benzene@40 °C). Record the table name as `source_table`. | `data/raw/chemicals_lib/antoine_all_sets.csv`, `critical.csv` |
| T04 | P1 | AUTO | NIST WebBook (existing `webbook_antoine.py`) | Re-run for ALL registry compounds (not only the 24). Also capture, per compound page, `ΔcH°liquid` entries with their cited method and reference (needed by T20), Tc/Pc where listed, and triple/normal-boiling data. Keep every set; do not collapse to "best". | `antoine_all_sets.csv`, `webbook_thermo.csv` |
| T05 | P1 | AUTO | union of T01–T04 + ThermoML points already in `measurements_dedup.csv` | Refit `A,B,C` per compound **per temperature window**; report coverage at T_flash±20 K. Do not interpolate heavy odd alkanes if T01 supplies them. Keep `antoine_source ∈ {refit_measured, webbook, ruzicka_majer, chemicals_lib, interpolated}`. | `antoine_refit.csv` (v2) |
| T06 | P3 | MANUAL | Yaws Handbook of Vapor Pressure (Elsevier, >25,000 compounds, paywalled) | Only if the user's library provides it. Do not scrape. | – |

Acentric factor for COSTALD is the **SRK-tuned** one in the original method doc: compute it by fitting SRK to the Psat points from T05 (do not paste the Pitzer ω).

### 2.2 Flash point and LFL

Current state: Li et al. already loaded (14,696 rows) but only 789 pure compounds survive the hydrocarbon filter, cup type is mixed, LFL is mostly estimated (670 of 695).

| id | pri | mode | source | action | output |
|---|---|---|---|---|---|
| T07 | P1 | AUTO | `chemicals.safety` (T03 package): `T_flash` (preferred IEC 60079-20-1 (2010), secondary NFPA 497 (2008)), `LFL`, `UFL`, `T_autoignition` | Export for all registry CAS with the `method=` name of each hit. Cup type is not stated in the function docs: set `method=unspecified`, tier A. Replace estimated LFLs (`lfl_source = estimated`) wherever a tabulated LFL exists; keep `T_C` of LFL if given, otherwise assume 25 °C and flag. | `data/raw/lfl/lfl_chemicals.csv`, flash rows |
| T08 | P2 | AUTO | PubChem PUG-View "Experimental Properties" (flash point, vapor pressure, heat of combustion, density) — endpoint pattern `https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/<CID>/JSON?heading=<name>` (heading names: verify) | For each registry compound (CID via InChIKey), keep raw strings, units (°F/°C/K), reference text and any `closed cup`/`open cup` wording. Parse into rows; keep unparseable ones in a reject file. Tier B. | `data/raw/pubchem/` |
| T09 | P2 | AUTO | Li et al. rows already ingested | Re-derive `method` where the original reference states cup type; otherwise `unspecified`. Down-weight rows whose compound set overlaps T07 with disagreement >5 K (flag `conflict`). | update loader |
| T10 | P2 | SEMI | Zabetakis Bulletin 627 (already in `data/raw/pdfs/`) | Verify the OCR-read n-alkane LFL table against PDF page 21 (handoff §15 says nonane/decane cells were reconstructed). Add aromatics and cycloalkanes if tabulated. | `lfl_measured.csv` |
| T11 | P2 | MANUAL | CRC Report 647, AV-27-18, AV-11-09; AFRL-TR-2020-0017 (DTIC AD1093317, DTIC blocks scripts) | User downloads in a browser. Harness then transcribes flash/viscosity/density/freezing tables for real fuels. | `data/raw/crc/`, `data/raw/dtic/` |
| T12 | P3 | MANUAL | GESTIS substance database | Check terms of use before any bulk use; only manual lookups for specific hydrocarbons that T07/T08 miss. | – |
| — | — | — | **Do not use**: Cheméo for flash point (small coverage, mixes calculated values). | | |

### 2.3 Octane (RON / MON) — gasoline domain

Currently: MON 267 pure / 30 binary / 76 multi (Leenhouts), **RON 0**. All sources below are measured RON/MON unless a column says "estimated".

| id | pri | mode | source | action | output |
|---|---|---|---|---|---|
| T13 | P1 | SEMI | vom Lehn et al. 2021, *A property database of fuel compounds with emphasis on spark-ignition engine applications*, Appl. Energy Combust. Sci. 5, 100018 (615 molecules; RON, MON, sensitivity, vapour pressure, density, LHV, etc.). Supplementary data also linked from `https://itv.rwth-aachen.de/en/research/ongoing-research/sustainable-energy-carriers/fuel-database` | Download the SI spreadsheet. Find the column or flag that separates experimental values from ANN/QSPR-estimated ones (**the paper fills missing values with a model**). Load experimental only as tier B; load estimated rows as tier C into a separate file. | `data/raw/vomlehn2021/` → `ron`, `mon`, `sensitivity` |
| T14 | P1 | AUTO | Rittig et al. 2022 GitLab: `git clone https://git.rwth-aachen.de/avt-svt/public/graph_ML_fuel_design.git` | Inspect which files hold RON/OS training data and how many are experimental. Report columns and license text before ingesting. | `data/raw/rittig2022/` |
| T15 | P1 | MANUAL/SEMI | Jameel et al. 2018 (KAUST), *Predicting octane number using NMR spectroscopy and ANNs*, Energy & Fuels — SI has RON/MON for 128 pure hydrocarbons, 123 hydrocarbon–ethanol blends **with known composition**, 30 FACE gasoline–ethanol blends | Locate the SI (publisher page or KAUST repository). This is the main source of composition-known octane blends. Blend basis is **volume**; store `x_mole`/`w_mass` after conversion with densities, and keep the original volume fractions in `raw_value`. | `data/raw/jameel2018/` |
| T16 | P2 | MANUAL | Sarathy-group blend papers: "A blending rule for octane numbers of PRFs and TPRFs with ethanol" (Fuel 2016); "The octane numbers of ethanol blended with gasoline and its surrogates" (Fuel 2014); Mehl/Sarathy et al. octane–ignition-delay paper with 79 hydrocarbon surrogate mixtures | Get the SI/tables with measured RON & MON for each blend. | `data/raw/sarathy_blends/` |
| T17 | P2 | MANUAL | API Research Project 45 / ASTM STP 225 *Knocking Characteristics of Pure Hydrocarbons* (1958; paywalled) | Primary measured pure-hydrocarbon values; use as an independent check on T13–T15 (flag disagreements >3 octane numbers). Pure-component values do not predict mixtures accurately (the book itself says so). | `data/raw/api45/` |
| T18 | P3 | MANUAL | FACE and CRC AVFL gasoline reports with measured RON/MON and composition (PIONA/DHA) | For real-gasoline hold-out. | – |
| — | — | — | **Do not use**: `glezdiazh/BIOFUELS.PTML` `SI03.xlsx` (synthetic). | | |

Elements: octane data are dominated by oxygenates (ethanol, MTBE/ETBE). Loader must allow `{C,H,O}` for these rows (`domain=gasoline`); do not silently drop them with the hydrocarbon-only filter.

### 2.4 Smoke point (TARGET since v2; ASTM D1322, mm; middle-distillate domain)

Smoke point is now a model target, so the loader must keep TSI, YSI and mm as separate properties, keep method/lamp, and store censored values ("above lamp range") as bounds with a `qualifier`. Pure-compound smoke points are rarely tabulated one by one. The literature converts them to the **threshold sooting index** (TSI), which normalises smoke point using molecular weight; store TSI as the primary value and derive smoke point only when needed.

| id | pri | mode | source | action | output |
|---|---|---|---|---|---|
| T19 | P2 | SEMI | Boehm et al., Energy & Fuels 2022, 36, 1916 (doi 10.1021/acs.energyfuels.1c03794), "TSI of sustainable aviation fuel candidates from composition input alone". Open copy: `https://rosap.ntl.bts.gov/view/dot/64807/dot_64807_DS1.pdf` | Extract the 65-molecule / 124-point TSI table, the TSI anchors (methylcyclohexane = 5, 1-methylnaphthalene = 100) and the **TSI↔smoke-point conversion constants** actually used. Do not assume constants; copy them from the paper. | `tsi`, tier A/B |
| T20 | P2 | SEMI | Mensch et al. 2010 (Michigan Tech, "Sooting characteristics of surrogates for jet fuels") | Extract recommended TSI values for jet-surrogate components. | `tsi` |
| T21 | P2 | AUTO | Yield Sooting Index database, Harvard Dataverse (canonical; ≥400 compounds, "Yield Sooting Index Database Volume 2") — find via `https://dataverse.harvard.edu` search "Yield Sooting Index Database" or the link on `https://pfefferlelab.yale.edu/yield-sooting-index-ysi` | Download the CSV(s). Store `ysi`. Do not scrape `ysi.ml.nrel.gov` (it is a prediction tool). | `data/raw/ysi/` |
| T22 | P3 | AUTO | YSI→smoke point ("derived smoke point") correlation from the terpene sooting paper | Compute only as tier C, `derived=True`. | – |
| T41 | P1 | MANUAL/SEMI | Blends and real fuels with known composition AND sooting data: Das, McEnally, Kwan, Zimmerman, Cannella et al., "Sooting tendencies of diesel fuels, jet fuels, and their surrogates in diffusion flames" (listed in the YSI reference lists; find the Fuel paper and its SI); binary/ternary TSI mixtures in Mensch 2010 and Boehm 2022 | Extract fuel/blend compositions, TSI or YSI, and smoke point where measured. This is the data that lets the blending rule be tested and gives GC-based compositions for jet fuels (also helps §2.8). Record the scale (TSI/YSI/mm) and instrument for each row; no cross-scale mixing. | `tsi`, `ysi`, `smoke_point_mm`, `multi`/`real_fuel` |
| T42 | P2 | AUTO | derived | Build `baselines_smoke.py` inputs: per-compound TSI table (T19+T20), conversion constants (from T19), and a censoring flag for values reported as bounds (e.g. ">50 mm"). | `tsi_table.csv` |
| T23 | P2 | SEMI | AJFTD (existing): smoke point for 82 real fuels; DTIC AD0763097 (smoke point vs hydrogen content, includes ASTM D1322 reference-blend data) — `https://apps.dtic.mil/sti/pdfs/AD0763097.pdf` (likely blocked: MANUAL) | Real-fuel validation set for the TSI blending model. | `smoke_point_mm` |

### 2.5 GCV (gross calorific value) and heat content

Physics: the heat of combustion of an ideal liquid blend is a mass-weighted sum, so the baseline `GCV_mix = Σ w_i GCV_i` is expected to be very good; the open question is **real fuels** (composition unknown) and measurement basis, not mixing.

| id | pri | mode | source | action | output |
|---|---|---|---|---|---|
| T24 | P1 | AUTO | NIST WebBook `ΔcH°liquid` (from T04 pages; each entry cites the calorimetry method, e.g. `Ccb`, and reference) | Convert to `gcv_MJ_kg = −ΔcH[kJ/mol] / M[g/mol]`. Keep every entry; dedupe by median; flag spread. Standard state 25 °C, liquid fuel, liquid water. Tier A when the entry is a measured `Ccb`, tier C when NIST says calculated. | `gcv_MJ_kg` pure |
| T25 | P1 | SEMI | NBS/API-44 reports as cross-checks (free): Rossini paraffins `https://nvlpubs.nist.gov/nistpubs/jres/34/jresv34n3p263_A1b.pdf`; cyclopentanes/cyclohexanes `https://nvlpubs.nist.gov/nistpubs/jres/36/jresv36n5p463_A1b.pdf` | Spot-check ≥10 values from T24 against these. | cross-check report |
| T26 | P2 | AUTO | `chemicals.combustion` (`HHV_stoichiometry`, `combustion_data`; about 1,400 chemicals) | Fill hydrocarbons lacking a measured entry from formation enthalpy; tier C, `derived=True`. Modified-Dulong (`HHV_modified_Dulong`) only for fuels where C,H,O,S wt% are known; record that it is an empirical coal-style correlation and expect larger error than the stoichiometric route. | `gcv_MJ_kg` derived |
| T27 | P1 | AUTO | AJFTD: `net_heat_content` and `hydrogen_wt_pct` for real fuels (already scraped) | Convert net→gross: `GCV ≈ NCV + 0.218 × H_wt%` MJ/kg (water formed = 8.94 × H mass, latent heat 2.44 MJ/kg at 25 °C). Store both; state the reference temperature; if AJFTD lists the test method (e.g. D4809 / D3338) keep it in `method`. | `ncv_MJ_kg`, `gcv_MJ_kg` real_fuel |

### 2.6 Pour point (and its proxies: cloud, freezing)

No open curated pour-point dataset exists (§8.2); AJFTD pour values are bounds (−80…−60 °C). Realistic expectation: tens to a few hundred records, mostly diesel/crude and mostly from papers. Be honest in `DATA_LOG.md` about the final count.

| id | pri | mode | source | action | output |
|---|---|---|---|---|---|
| T28 | P1 | MANUAL/SEMI | Krishna et al., *Erdöl Erdgas Kohle* 1989 (also Hydrocarbon Technology): 24 diesel blends with stream composition, ASTM distillation, wax content, n-paraffin content, cloud point, pour point, CFPP (Table II). Copy at `https://krishnapdf.cocosimulator.org/Krishna_Erdol_und_Erdgas_1989.pdf` | Transcribe Table II (composition of streams per blend + properties). Check the copyright/hosting is acceptable before redistributing; keep transcription only in the private repo. | `pour_point_C`, `cloud_point_C`, blend composition |
| T29 | P2 | MANUAL | Coutinho & Daridon, "Cloud and pour points in fuel blends", Fuel 81 (2002) 963–967; Khan/Joshi OGJ correlations (diesel blends of known volume composition) | Extract measured blend data and the model form used (compare with PPBI). | `data/raw/pour_point/` |
| T30 | P2 | SEMI | NASA NTRS 19860018789, *Partially Frozen Hydrocarbon Fuels* (cloud, freeze, wax appearance, pour points for eight fuels; n-alkane content) — `https://ntrs.nasa.gov/api/citations/19860018789/downloads/19860018789.pdf` | Transcribe the fuel property tables. | pour/freeze, real_fuel |
| T31 | P3 | SEMI | Kirgina et al. 2019 (straight-run diesel; cloud, CFPP, pour point with n-paraffin content), `https://www.vurup.sk/wp-content/uploads/2020/03/PC-X-2019_Kirgina-168.pdf` | Transcribe sample table. | pour/cloud |
| T32 | P2 | MANUAL | CRC AV-11-09 freezing point tables; Lira thesis (Michigan State, `d.lib.msu.edu/etd/3010`); Daridon 2002 | As §16.3-A4. | freezing/pour |

Feature needed for any pour-point model: **n-paraffin content or carbon-number distribution** of the fuel. Pour point records without it are tier C for learning.

### 2.7 Liquid density and real-fuel density vs T

Pure/binary/multi density is plentiful (ThermoML ≈ 20.7k rows). Gap: real-fuel density vs temperature (AJFTD has 15.56 °C only).

| id | pri | mode | source | action | output |
|---|---|---|---|---|---|
| T33 | P1 | SEMI | NIST RP-1: Energy & Fuels 2012, doi 10.1021/ef3004009, eleven orthogonal blends, density 278–343 K, viscosity 263–373 K, ambient pressure (`https://www.nist.gov/node/574906`) | Obtain SI tables (density, speed of sound, viscosity) **and the blend composition**. Best single real multi-component validation set for density and viscosity. | `density_kg_m3`, `viscosity_*`, real_fuel |
| T34 | P2 | SEMI | NIST Jet-A / JP-8 / S-8 density papers: Outcalt & Laesecke, Energy & Fuels (JP-8 comparison, 2010; `https://www.nist.gov/node/592016`) and the Jet-A/S-8 compositions paper | Extract the density correlations and data (vibrating-tube; ~0.1 % uncertainty). | `density_kg_m3` real_fuel |
| T35 | P3 | AUTO | NIST WebBook fluid data (74 fluids) / `thermo` correlations | Pure-fluid liquid density vs T where ThermoML is thin. | – |

### 2.8 Composition of real fuels (biggest cross-cutting gap)

| id | pri | mode | source | action |
|---|---|---|---|---|
| T36 | P1 | MANUAL | AFRL-TR-2020-0017 (DTIC AD1093317): GC×GC composition + properties of jet fuels | User downloads. Harness extracts class × carbon-number tables and links to flash/viscosity/density/GCV records by fuel id. |
| T37 | P1 | SEMI | NIST RP-1 blend compositions (see T33) | Needed to turn the 11 blends into molecule-level mixtures. |
| T38 | P2 | MANUAL | CRC AV-27-18, AV-11-09, 647 composition tables | As T11. |

Without composition use the pseudo-component representation in §16.3-A3 and say so in the report.

### 2.9 Critical constants and COSTALD inputs (if the COSTALD density route from the original method doc is kept)

| id | pri | mode | source | action |
|---|---|---|---|---|
| T39 | P2 | AUTO | `chemicals` critical-property tables (several sources, including Yaws) + WebBook (T04) | Tc, Pc, Vc per compound with source label; resolve disagreements by priority NIST/JPCRD > Poling-type tables > Yaws. |
| T40 | P2 | AUTO | derived | Compute SRK-tuned ω from T05 Psat points; characteristic volume V* from the density data (COSTALD fit). Mark `derived=True`. |

---

## 3. Execution order

0. T00 (source check) and report to the user.
1. T03, T04, T07 (everything automatic that unblocks the flash baseline).
2. T01 and T05 (n-alkane vapour pressure for C13–C20, refit).
3. T24, T27 (GCV pure + real fuels; cheap).
4. T13, T14, T15 (octane; decide domain per §17 of the handoff).
5. T19, T20, T21, T42 (smoke point / TSI), then T41 (blends and fuels with composition).
6. T08, T09, T10 (flash-point depth and LFL quality).
7. T33, T36, T37 (real-fuel composition and density; unlocks validation for everything).
8. Pour point tasks T28–T32 last (lowest expected yield).
9. MANUAL list for the user: T06, T11, T12, T15 (if blocked), T16, T17, T18, T23, T28 (if needed), T29, T32, T36, T38, T41.

## 4. Done when

* `DATA_LOG.md` has an entry per task.
* `python build_unified.py && python audit.py --bundle` passes with no FAIL.
* `UNIFIED_REPORT.md` shows, **per target property**: counts by `kind`, number of independent `system_id`s with ≥3 components,
  number of real-fuel records, fraction tier A/B vs tier C, and the feasibility gate (G1/G2/G3) from `PROJECT_HANDOFF.md` §17.4.
* The user gets a short list of what remains MANUAL.
