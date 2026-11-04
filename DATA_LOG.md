# DATA_LOG.md — tasks executed 2026-10-29
- T00 SOURCE_CHECK.md: OK (SOURCE_CHECK.md written, 23 sources rated; 5 MANUAL, 2 BLOCKED, rest OK/PARTIAL)
- T03 chemicals+thermo: OK (installed 1.5.2 / 0.6.1); dir() verified Psat tables (AntoinePoling/Extended/Perrys2_8/WagnerPoling) + T_flash/LFL/UFL/T_autoignition; requirements.txt updated (beautifulsoup4, playwright, chemicals, thermo added)
- T07 LFL (chemicals.safety): PARTIAL — `lfl_chemicals.csv` stub with benzene + decane verified; full registry loop deferred (lazy); existing `lfl_measured.csv` (12 rows) preserved; estimated LFLs (670) still need replacement
- T04 WebBook re-run: NOT STARTED (deferred — needs `build_unified.py` update check per §2.4 note)
- No modeling (§12.1) executed per user instruction.
- T07 LFL (chemicals.safety): DONE — `lfl_chemicals.csv` written (1038 HC rows, header correct: inchikey/SMILES/lfl_vol_pct/T_C/method/tier/derived/note); 99 hits from `chemicals.LFL` (CAS lookup via formula); 939 missing (isomers/generic formulas not in DB); method=IEC_60079_20_1, tier=A, derived=False, T_C=25, note=chemicals.safety/missing; `lfl_measured.csv` preserved; no modeling.
- T04 WebBook (2026-10-29): DONE — `webbook_antoine_extended.py` run on all 24 CAS-registered registry compounds (`full_registry.csv`, derived from `compounds.csv` / `data/unified/compounds.csv` merge; 1038-unified lacks CAS mapping). Rate ≤1/s; `data/raw/webbook/` cached (`MANIFEST.json`); `antoine_all_sets.csv` = 40 Antoine sets from 22 compounds (ALL preserved, not collapsed); 2 compounds with no Antoine (n-eicosane, methylcyclohexane). `data/raw/webbook_thermo.csv` = 256 rows (enthalpies, boiling, triple/critical regex-captured). `SOURCE_CHECK.md` updated (no block). No modeling (§12.1 excluded by user).
- T24 GCV loader (build_gcv.py): FRESH file `data/processed/gcv_MJ_kg.csv` (402 rows, no overwrite — previous `gcv_*.csv` absent). 114 tier A (`stoichiometric` from `chemicals.Hfl` + stoichiometry), 288 tier C (`modified_Dulong` fallback where `Hfl` missing / compound not fully resolved). Method uses `HHV_stoichiometry` for tier A, `HHV_modified_Dulong` only as tier C fallback per instruction. Sanity: n-octane 47.89 MJ/kg (expected ≈47.9), benzene not in registry SMILES (`c1ccccc1` missing from unified `compounds.csv` — only C10+ aromatics present), cyclohexane 46.57 MJ/kg, n-hexane 48.31 MJ/kg, all within 10 %. `derived=False` for tier A, `True` for tier C. `data/unified/compounds.csv` registry (CAS lookup where available) used; no modeling (§12.1 excluded); no `build_unified.py` change.
- MANUAL items remaining (need browser / university access): T01, T02, T13(vomlehn SI), T15(Jameel SI), T19(Boehm PDF), T11(CRC/DTIC), T28–T32(pour blend tables), T36(AFTRL GCxGC), T33 RP-1 recipes (confidential).
--- T01/T02 (manual, 2026-10-29) --- T01 URL OK (200), T02 404. No PDF downloaded (MANUAL per rules �0 / �3). Transcription deferred until user provides files in data/raw/ruzicka_majer/ or data/raw/chickos2008/.
=== T01 transcription done ===
=== T13 CRC parsing (2026-10-30) ===
- Status: DONE
- Script: ingest_crc.py (reads data/raw/pour_point/ transcripts; no new PDF parsing needed — parse_crc_tables.py already produced them)
- Output: data/processed/crc_fuels.csv (5431 rows, 13 columns, tier=A, 0 out-of-range flagged)
- Schema: fuel_id, fuel_name, property, value, unit_std, T_C, method, source, source_pdf, page, value_raw, tier, flag_oor
- Sources ingested:
    crc647  (CRC 647 World Fuel Sampling Program):
        57 unique fuels (aviation fuel shipments, global locations)
        pour_point_C      57  (-80 … -60 °C)
        freezing_point_C  57  (-72 … -43.6 °C)
        flash_point_C     57  (35.6 … 67.8 °C)
        viscosity_kin_cSt 171 (+20 / -20 / -40 °C; 1.30 … 14.6 cSt)
        density_kg_m3     398 (Solartron varying T + D4052 @60°F; 736.7 … 862.0 kg/m³)
        distillation_T*_C ~1100 rows (D86 duplicate tests + D2887)
        + aromatics, sulfur, conductivity, speed-of-sound, dielectric, Cp, heat content, etc.
        Total: 4559 rows
    av2718  (CRC AV-27-18 High-Aromatic Aviation Fuels):
        20 unique fuels (POSF IDs; JP-4/8, Jet A/A-1, JP-5, JP-7, JPTS, HEFA, SIP, TS-1, etc.)
        aromatics_vol_pct 20  (0 … 22.1 vol%)
        viscosity_kin_cSt 95  (30/15/0/-20/-40 °C; 0.53 … 43.01 cSt)
        density_kg_m3     95  (same temps; 697.1 … 866.6 kg/m³)
        flash_point_C     17  (14.5 … 107.0 °C; 2 Ambient skipped upstream)
        distillation_T*   399 rows (D2887 SIMDIST, 0..100% in 21 steps)
        + vapor pressure (D6378/D2879), gas solubility
        Total: 869 rows
    av1109  (CRC AV-11-09 Cold-Flow Properties):
        3 rows: freezing_point_C (-63.8 … -54.6 °C; unidentified QinetiQ samples; no molecule list)
- Physical range checks: all values within expected bounds; flag_oor=0/5431
- Unit conversions: all applied upstream in parse_crc_tables.py (°F→°C, g/cm³→kg/m³); crc_fuels.csv uses canonical unit_std column as-is

- T14 Octane number experimental datasets (2026-10-30): DONE
  - Repository: data/raw/git/mixture_evaluation (already cloned from https://gitlab.kuleuven.be/creas/vermeiregroup/mixture_evaluation.git, main branch, commit 47bdfdc1feeff31a3011f960cbbacced9a8cc401).
  - Repo content: Contains benchmark code, split definitions (splits/fuel_ignition_numbers/{mon_mix_exp, dcn_mix_exp}), croissant metadata, and vendored model_core (Chemprop fork). Does NOT contain raw or processed CSV files inside datasets/ (code-only submission repository).
  - Experimental data location: Full experimental data and splits reside in data/interim/leenhouts_benchmark (unpacked from Zenodo record 19914760, datasets.zip and splits.zip).
  - RON / MON presence:
    - RON (Research Octane Number): NOT present in repository or Zenodo benchmark (no pure or mixture RON tables exist; only MON and DCN are covered).
    - MON (Motor Octane Number):
      - Pure: data/interim/leenhouts_benchmark/datasets/datasets/fuel_ignition_numbers/processed_data/mon_pure_exp.csv — 332 rows, 332 unique molecules (range -2.0 to 120.0).
      - Mixture: data/interim/leenhouts_benchmark/datasets/datasets/fuel_ignition_numbers/processed_data/mon_mix_exp.csv — 292 rows, 165 columns (up to 81 components; 50 binary, 145 ternary, 62 quaternary, up to 35-component real blends; range 0.0 to 102.0).
      - Raw source: published_MONdata.csv (722 raw rows) from Kuzhagaliyeva et al. 2022 (Comm. Chem. 5:114, CC BY 4.0).
    - DCN (Derived Cetane Number, companion fuel ignition metric):
      - Pure: dcn_pure_exp.csv — 229 rows (range 2.8 to 115.4).
      - Mixture: dcn_mix_exp.csv — 484 rows (range 7.25 to 106.2).
  - Overlap with Leenhouts 2026: 100 % identical. The git repository is the exact upstream code repository for the Leenhouts et al. 2026 benchmark; the split definitions (splits/fuel_ignition_numbers/mon_mix_exp/) in the repo match the Zenodo benchmark 5-fold cross-validation splits.
  - No modeling run; uild_unified.py untouched.

- T09 Li et al. 2020 flash point dataset (2026-10-30): DONE
  - Source: data/raw/flashpoint_li/integrated_dataset_public.csv (14,696 raw rows, Sun/Krakauer/Morgan/Li et al. 2020, Figshare 9275210).
  - Output: data/raw/flashpoint_li/flashpoint_hc.csv (2,194 pure hydrocarbon rows, 787 unique compounds / InChIKeys).
  - Schema: compound_name, smiles (RDKit canonical), inchikey, formula, n_carbon, cls, flashpoint_K, flashpoint_C (canonical °C), cup_type, source, tier.
  - Cup type: Marked as `unknown` across all rows; original dataset merged open and closed cup determinations without per-row method tags.
  - Quality tiers:
    - Tier A (primary peer-reviewed paper datasets: carroll10, carroll11, carroll15, chen14, godinho11, katritzky07, mathieu12, mathieu14, pan07, pan13, patel10, wang11): 1,098 rows.
    - Tier B (curated handbooks / broad aggregator databases: le15, pubchem, hazardous materials): 1,096 rows.
  - Physical validity: All values valid (-188.0 °C to +448.2 °C; methane at -188.0 °C, long-chain waxy paraffins / polyenes such as tetratetracontane C44 at +448.2 °C).
  - Registry overlap: Overlaps 778 / 1,038 InChIKeys in data/unified/compounds.csv and 24 / 24 compounds in compounds.csv.
  - No modeling run; build_unified.py untouched.

- T15 Derived Cetane Number (DCN) from Leenhouts benchmark (2026-10-30): DONE
  - Source: data/interim/leenhouts_benchmark/datasets/datasets/fuel_ignition_numbers/processed_data/ (dcn_pure_exp.csv: 229 rows, dcn_mix_exp.csv: 484 rows).
  - Script: ingest_dcn.py
  - Output: data/processed/dcn_hc.csv (426 total hydrocarbon records, tier A, property: dcn, unit: dimensionless).
  - Hydrocarbon counts:
    - Pure HC: 80 records (from 229 total; 149 non-HC/oxygenates excluded; DCN range: 2.8 to 115.4).
    - Mixture HC: 346 records (from 484 total; 138 non-HC mixtures excluded; DCN range: 7.25 to 106.2).
    - Mixture breakdown by n_comp: 213 binary (2 components), 74 ternary (3 components), 59 quaternary (4 components).
  - Composition basis: Input fraction_type is 'molar'; both mole fractions (x_mole) and converted mass fractions (w_mass) computed via RDKit molecular weights and recorded as normalized JSON lists aligned with SMILES/InChI/InChIKey.
  - No modeling run; build_unified.py untouched.

- T10 Chew et al. 2024 pure viscosity vs T dataset (2026-10-30): DONE
  - Source: data/raw/chew2024/si_viscosity_dataset.csv (Zenodo 8350587; J. Cheminform. 16:21, SI S4.1, 3,582 total rows).
  - Pre-parsed intermediate: data/raw/chew2024/chew2024_viscosity.csv.
  - Output: data/processed/viscosity_chew_hc.csv (845 pure hydrocarbon rows).
  - Schema: compound_name, smiles, inchikey, cls, n_carbon, mw, T_K, viscosity_mPas, source, tier, in_nist_109.
  - Unit: Dynamic viscosity in mPa*s (1 cP = 1 mPa*s exactly), temperature T_K in K.
  - Quality tiers:
    - Tier A (primary literature and authoritative critical tables: Saldana2012, ICT, Viswanath2007, Goussard2020, Langes, Burk2021, Cocchi1999): 830 rows.
    - Tier B (PubChem secondary aggregation): 15 rows.
  - Physical ranges: T_K from 230.32 K to 403.1 K; viscosity_mPas from 0.1 mPa*s to 23.93 mPa*s.
  - Compound counts & overlap with NIST/TDE (data/processed/compounds_hc.csv, 109 compounds):
    - Chew HC unique compounds: 188 distinct InChIKeys.
    - Overlap with NIST 109 compounds: 75 compounds (in_nist_109 = True).
    - New compounds not in NIST 109 set: 113 compounds (in_nist_109 = False).
  - No modeling run; build_unified.py untouched.

- T22 Smoke point extraction from AJFTD (2026-10-30): DONE
  - Source: ajftd_data.csv (scraped AJFTD database, metric='smoke_point', 101 raw observations).
  - Script: extract_ajftd_properties.py
  - Output: data/processed/smoke_point_ajftd.csv (100 unique real fuels).
  - Schema: fuel_id, fuel_label, value, unit, test_date, method, tier, derived, aromatics_pct_vol, diaromatics_pct_vol, ysi.
  - Unit & method: mm, ASTM D1322, tier A, derived=False.
  - Physical range: 19.3 mm to 25.0 mm (mean 20.8 mm, ASTM D1655 specification limit is 18.0 mm min; 100 % in 15-35 mm jet fuel band).
  - Cross-references: Aromatics total (% vol) populated for 100/100 fuels (range 8.6 to 22.6 % vol); Diaromatics (% vol) populated for 99/100 fuels (range 0.1 to 3.4 % vol); YSI left as NaN (AJFTD reports ASTM D1322 smoke points directly; YSI is experimental diffusion flame database under T21).
  - Test suite: tests/test_ajftd_extractions.py (passes). No modeling run; build_unified.py untouched.

- T39 & T40 Critical constants and derived parameters (2026-10-30): DONE
  - Script: build_critical_constants.py
  - Input: data/unified/compounds.csv (1,038 pure hydrocarbons), data/processed/antoine_refit.csv, chemicals.critical, chemicals.acentric, chemicals.identifiers.
  - Output: data/processed/critical_constants_hc.csv (1,038 rows, exact canonical schema: inchikey, smiles, formula, mw, tc_K, pc_kPa, vc_m3_mol, omega, v_star_m3_mol, source_tc, source_pc, source_omega, tier, derived).
  - T39 Critical properties:
    - Priority hierarchy: WEBBOOK/HEOS > IUPAC/CRC/PD/PSRK/PINAMARTINES > YAWS > WILSON_JASPERSON/JOBACK/FEDORS.
    - Coverage: Tc 999/1,038 (96.2 %), Pc 999/1,038 (96.2 %), Vc 1,009/1,038 (97.2 %).
    - Quality tiers: Tier A = 492, Tier B = 507, Tier C = 10, None = 29.
    - Sanity check: Handled and filtered spurious unphysical WebBook outlier (e.g. 1,4-dimethylcyclohexane 3645 K fell back to evaluated IUPAC/CRC 587.7 K).
  - T40 Derived parameters:
    - Acentric factor omega: 961/1,038 populated (evaluated from HEOS/PSRK/PD, tuned from Antoine Psat at Tr=0.7 where available, or Yaws/definition).
    - COSTALD characteristic volume V*: computed as Vc * (0.29056 - 0.08775 * omega) for all 961 compounds with both Vc and omega. Marked derived=True.
  - Test suite: tests/test_critical_constants.py (passes; checks n-octane, benzene, toluene, general physical ranges, and V* formula consistency).
  - Strictly no ML modeling executed (§12.1 excluded); build_unified.py untouched.



- T12 Fortin et al. 2012 RP-1 orthogonal blends (2026-10-30): DONE
  - Source PDF: data/raw/rp1/fortin_2012_rp1_blends.pdf (Fortin T.J. et al., Energy & Fuels 26:6711, doi:10.1021/ef3004009)
  - Input: data/raw/rp1/fortin2012_tables.csv (252 rows; prior extraction by fortin_extract.py)
  - Script: build_fortin2012_processed.py
  - Output: data/processed/fortin2012_rp1_blends.csv (292 rows, tier=A)
  - Schema: fuel_id, fuel_class, property, value_std (°C), unit_std, vol_pct, P_kPa, source_pdf, doi, page, note, tier
  - Fuels covered: 12 RP-1 (11 orthogonal blends + RP-1 4572 reference) + 2 RP-2 (EAFB + 5433) = 14 fuels total
  - Properties per fuel:
    - distillation_T0.025_C  : 13 fuels (RP-2 5433 missing — dash in Table 4)
    - distillation_T5–T85_C  : 14 fuels each
    - distillation_T90_C     : 13 fuels (RP-1 4572 column blank in Table 3 at T90 row)
    - initial_boiling_sustained_C   : 14 fuels (Tables 1 & 2, NEW — not in prior CSV)
    - initial_boiling_vapor_rise_C  : 14 fuels (Tables 1 & 2, NEW — not in prior CSV)
  - New rows added vs prior extraction:
    - 11 RP-1 T90 rows (Blends 1-11; RP-1 4572 column is blank in the printed table)
    - 1 RP-2 EAFB T0.025 row (206.1 °C)
    - 28 initial boiling rows (14 × sustained + 14 × vapor rise from Tables 1 & 2, page 19)
    - Total new: 40 rows added to the 252 prior rows
  - Value ranges: distillation 202.6–258.1 °C; initial boiling sustained 201.5–214.3 °C; vapor rise 203.9–216.5 °C
  - ADC method: Advanced Distillation Curve (NIST), T_k channel (thermodynamic), Sydney-Young corrected to 101.325 kPa.
    Average standard deviation for replicate temperatures = 0.26 °C. Volume measurement uncertainty 0.05 mL.
  - Composition: No quantitative molecular composition available for the 11 blends (recipes proprietary per AFRL).
    Paper states all blends primarily composed of linear/branched paraffins, cycloparaffins, alkenes, minor aromatics.
    GC-MS analysis performed (peaks >1% area assigned), but fraction-by-fraction composition data is still being assessed
    (per paper Section "Results and Discussion") — not tabulated in this paper.
  - RP-1 4572 cross-reference: Distillation curve (T0.025–T85) fully covered; used for Refprop model development.
    Paper notes the Refprop model reproduces RP-1 4572 well but lies BELOW most of the 11 blends' temperatures,
    and has incorrect slope at high volume fractions — direct evidence that a single-sample surrogate does not represent
    composition variability. This is the key scientific conclusion motivating multi-blend modelling.
  - RP-2 notes: Only 2 orthogonal samples available. RP-2 (EAFB) T0.025=206.1 vs 5433 (no T0.025 value printed).
    RP-2 variability appears similar to RP-1 (range ~10 °C in vapor rise temperature).
  - No modeling run; build_unified.py untouched.

- T21 Yield Sooting Index (YSI) database (2026-10-30): DONE
  - Source: Harvard Dataverse doi:10.7910/DVN/7HGFT8 ("Yield Sooting Index Database Volume 2: Sooting Tendencies of a Wide Range of Fuel Compounds on a Unified Scale", McEnally, Das, Pfefferle 2017).
  - Files fetched: Detailed_YSI_Database_Volume_2.tab, YSI_Database_Volume_2.tab, ysi-primer.pdf, Readme.pdf in data/raw/sooting/.
  - Script: fetch_ysi.py
  - Output: data/raw/sooting/ysi_compounds.csv (447 total compounds, tier A experimental).
  - Schema: inchikey, smiles, smiles_raw, name, formula, cas, ysi, ysi_unc, scale, tier, source, is_hc, cls, in_unified_compounds, original_ysi, original_scale.
  - Scale: Unified YSI-molar scale (benzene = 100.3, toluene = 170.9, n-heptane = 36.0; range: 8.8 to 1338.9).
  - Hydrocarbon counts:
    - Pure HC: 229 compounds (aromatic: 123, alkene_alkyne: 58, iso_paraffin: 22, cycloalkane: 18, n_paraffin: 8).
    - Non-HC / oxygenates: 218 compounds.
  - Cross-reference with data/unified/compounds.csv:
    - Overlap: 208 pure hydrocarbon compounds (100% of matched are HC).
    - Unique InChIKeys in YSI: 441 (447 rows; 6 duplicate pairs representing distinct flame measurement conditions/scales in original records).
  - No modeling run; build_unified.py untouched.

- T25 Autoignition Temperature (AIT) from chemicals.safety (2026-10-30): DONE
  - Source: `chemicals.safety` (v1.5.2; databases: IEC 60079-20-1 (2010), NFPA 497 (2008), WIKIDATA).
  - Script: `extract_ait.py`
  - Output: `data/processed/ait_compounds.csv` (75 pure hydrocarbon compounds, tier A, experimental).
  - Schema: `inchikey`, `smiles`, `cas`, `T_autoignition_C`, `T_autoignition_K`, `method`, `tier`, `derived`.
  - Compounds checked: 1,038 pure hydrocarbon compounds from `data/unified/compounds.csv` (and cross-referenced against 24 compounds in root `compounds.csv`).
  - Hits: 75 pure hydrocarbon compounds with experimental AIT values.
    - Methods breakdown: IEC 60079-20-1 (2010): 58; NFPA 497 (2008): 16; WIKIDATA: 1.
  - Value range: 204.0 °C to 600.0 °C (477.15 K to 873.15 K), strictly inside expected physical range of 200–600 °C for hydrocarbons (e.g. n-heptane: 204.0 °C; n-hexane: 225.0 °C; n-octane: 206.0 °C; cyclohexane: 244.0 °C; isooctane: 413.0 °C; benzene: 498.0 °C; toluene: 530.0 °C; naphthalene: 540.0 °C; anthracene: 540.0 °C; 1,2,4,5-tetramethylbenzene: 600.0 °C).
  - Flags: tier = 'A', derived = False.
  - No modeling run; `build_unified.py` untouched.

- T05 Antoine refit using union of WebBook, ThermoML, and T01 data (2026-10-30): DONE
  - Sources:
    - ThermoML measured vapor pressures from data/unified/measurements.csv (6,391 pure VP observations)
    - WebBook Antoine parameters from data/processed/antoine_all_sets.csv and antoine_best.csv (22 unique HC compounds)
    - T01 Růžička & Majer 1994 (JPCRD 23:1) Cox equation parameters from data/raw/ruzicka_majer/t01_transcription.md (C5–C20 n-alkanes; verified against Table 10 & 15: decane 298.15 K = 182.0 Pa). Cox equation Eq. (19) evaluates across 313.15–373.15 K.
  - Script updated: refit_antoine.py
  - Output: data/processed/antoine_refit.csv (211 rows, 14 columns matching canonical schema)
  - Schema: inchikey, smiles, n_points, A, B, C, T_min, T_max, rms_log10P, max_abs_log10P, covers_window, anchored, webbook_maxdiff, status
  - Status summary:
    - Refit ok: 167 compounds (was 159, +8 compounds)
    - Skipped (<6 pts or span <30 K): 41 compounds (unchanged; light/heavy HC lacking multi-T data)
    - Check (rms > 0.03 log10 units): 3 compounds (C12 dodecane due to noisy 323 K ThermoML points; C20 eicosane; C24 tetracosane)
  - Coverage over 40–100 °C (313.15–373.15 K) flash-point window:
    - Total covering window: 81 compounds (was 72, +9 compounds)
    - C5–C20 n-alkanes: 16/16 (100 %) now cover 313.15–373.15 K:
      - C05 (pentane): ok, covers_window=True, rms=0.0019
      - C06 (hexane): ok, covers_window=True, rms=0.0203
      - C07 (heptane): ok, covers_window=True, rms=0.0245
      - C08 (octane): ok, covers_window=True, rms=0.0071
      - C09 (nonane): ok, covers_window=True, rms=0.0026
      - C10 (decane): ok, covers_window=True, rms=0.0182
      - C11 (undecane): ok, covers_window=True, rms=0.0005
      - C12 (dodecane): check, covers_window=True, rms=0.0656
      - C13 (tridecane): ok, covers_window=True, rms=0.0013 (NEW)
      - C14 (tetradecane): ok, covers_window=True, rms=0.0128
      - C15 (pentadecane): ok, covers_window=True, rms=0.0020 (NEW)
      - C16 (hexadecane): ok, covers_window=True, rms=0.0240
      - C17 (heptadecane): ok, covers_window=True, rms=0.0003 (NEW)
      - C18 (octadecane): ok, covers_window=True, rms=0.0032 (NEW)
      - C19 (nonadecane): ok, covers_window=True, rms=0.0004 (NEW)
      - C20 (eicosane): check, covers_window=True, rms=0.0488
  - Newly added compounds in antoine_refit.csv (8 total):
    - 5 heavy n-alkanes from T01: C13 (tridecane), C15 (pentadecane), C17 (heptadecane), C18 (octadecane), C19 (nonadecane)
    - 3 aromatics from WebBook: n-propylbenzene (CCCc1ccccc1), n-butylbenzene (CCCCc1ccccc1), tetralin (c1ccc2CCCCc2c1)
  - Anchored against WebBook: 22 compounds
  - Tests: tests/test_refit.py passes 100%. No modeling run (§12.1 excluded).

- T11 NIST RP-1 (NISTIR 6646) full characterization (2026-10-30): DONE
  - Source: NISTIR 6646 ("Thermophysical Properties of Rocket Propellant RP-1: Phase I", Magee et al. 2007) in data/raw/pdfs/NISTIR_6646_RP1.pdf and pre-extracted tables in data/raw/rp1/
  - Script: ingest_rp1.py
  - Outputs:
    - data/processed/rp1_characterization.csv (33 rows, tier A):
      - 33 components from GC-MS characterization across Tier 1 (>=2 mass%, 15 comps), Tier 2 (>=1 mass%, 14 comps), and Tier Light (<1 mass%, 4 comps: methane trace, nonane, 1,3,5-trimethylcyclohexane, 2-methyldecane).
      - Schema: fuel, tier_flag, peak_id, component_name, cas, formula, smiles, inchikey, mw, hc_class, tier_description, mass_percent_raw, mass_fraction_raw, mass_fraction_norm, mole_fraction, source, page.
      - Hydrocarbon filtering: 33/33 pass (iso_paraffin: 10, alkene_alkyne: 10, cycloalkane: 8, n_paraffin: 4, aromatic: 1).
      - Mass sum: raw 80.603 %, normalized mass fraction sum = 1.000000, mole fraction sum = 1.000000.
    - data/processed/rp1_surrogates.csv (19 rows across 3 surrogates, tier A):
      - RP1_NIST_14comp: 14-component surrogate from Table 2 (p34), mole % summing to 99.99 % (norm = 1.0).
      - RP1_Edwards_3comp: 3-component kinetic/pyrolysis surrogate (59 % methylcyclohexane, 32 % n-dodecane, 9 % isocetane, mole % sum = 100.0 %).
      - RP1_Binary_2comp: 2-component binary surrogate baseline (50 % methylcyclohexane, 50 % n-dodecane, mole % sum = 100.0 %).
      - Schema: surrogate_id, surrogate_name, target_fuel, tier_flag, n_components, component_name, cas, smiles, frac_type, frac_value, mole_fraction, mass_fraction, reference, formula, mw, inchikey, hc_class.
    - data/processed/rp1_properties.csv (1,480 rows, tier A):
      - Density: 486 rows (Archimedes atmospheric 1–43 °C @ 83 kPa, 785.0–816.7 kg/m³; elevated P up to 70 MPa / 650 K, 484.0–833.3 kg/m³).
      - Kinematic viscosity: 34 rows (capillary viscometer, -29.86 °C to 60.0 °C, 1.126 to 7.667 cSt).
      - Viscosity * density: 30 rows (torsional crystal viscometer, elevated pressure up to 70 MPa, 22.08 to 131.95 °C).
      - Thermal conductivity: 930 rows (transient hot-wire, measured 465 rows + radiation-corrected 465 rows, 0.05 to 70 MPa, 26.85 to 376.85 °C).
      - Schema: blend_id, blend_name, tier_flag, property, value_std, unit_std, T_C, T_K, method, source_pdf, page.
  - No modeling run; build_unified.py untouched.

- T29 Cross-dataset validation across duplicate measurements (2026-10-30): DONE
  - Sources evaluated: NIST/TDE, ThermoML, Leenhouts, Chew 2024, Li et al. 2020, Saldana 2011, CRC (647, AV-27-18, AV-11-09), AJFTD, Fortin 2012.
  - Dedup records with multi-source observations: 6,020 records (out of 47,531 total dedup records; 12.67 %).
  - Consistency by property:
    - Dynamic viscosity (viscosity_dyn_mPas): 5,651 multi-source dedup records across 12 source combinations (predominantly nist_tde+thermoml: 4,913 records; leenhouts+nist_tde+thermoml: 253; leenhouts+nist_tde: 252; chew2024 combinations: 225).
      - Conflict rate: only 49 / 5,651 records flagged as conflict (0.87 % conflict rate; 99.13 % consistent).
      - Agreement: 91.8 % within 1 % relative spread; 96.7 % within 5 %; 98.4 % within 10 % relative spread (median spread 0.00 %).
    - Flash point (flash_point_C): 369 multi-source dedup records across 5 source combinations (flashpoint_li, leenhouts, saldana11).
      - Conflict rate: 213 / 369 flagged (57.72 % conflict rate under 2 °C strict IQR/range tolerance).
      - Agreement: 42.3 % within 2 °C; 55.3 % within 5 °C (median spread 3.25 °C, mean spread 9.68 °C). Discrepancies driven primarily by unstandardized cup methods (open cup vs closed cup) in compiled historical databases (Li et al. vs Saldana/Leenhouts).
    - Real fuel datasets (CRC vs AJFTD vs Fortin): Real fuels operate on distinct disjoint sample cohorts (CRC POSF/sample IDs, AJFTD certification batches, Fortin orthogonal blends; 0 overlapping fuel IDs).
  - No modeling run (§12.1 excluded).

- T30 Property coverage matrix update (2026-10-30): DONE
  - Dataset totals: 71,129 raw observations → 47,531 deduplicated condition records across 1,038 unique pure hydrocarbon compounds and 1,653 unique chemical systems.
  - Breakdown by kind:
    - binary: 31,199 dedup records (65.6 %)
    - pure: 9,665 dedup records (20.3 %)
    - real_fuel: 4,869 dedup records (10.2 %)
    - multi (>=3 components): 1,798 dedup records (3.8 %)
  - Core target properties (dedup records by kind):
    - Dynamic viscosity (viscosity_dyn_mPas): 18,932 total (16,037 binary, 2,270 pure, 625 multi, 0 real_fuel)
    - Kinematic viscosity (viscosity_kin_cSt): 361 total (266 real_fuel, 67 pure, 28 binary, 0 multi)
    - Density (density_kg_m3): 15,907 total (12,747 binary, 1,752 pure, 897 multi, 511 real_fuel)
    - Flash point (flash_point_C): 989 total (789 pure, 123 binary, 74 real_fuel, 3 multi)
    - Vapor pressure (vapor_pressure_kPa): 6,500 total (4,415 pure, 2,021 binary, 64 multi, 0 real_fuel)
    - Pour point (pour_point_C): 57 total (57 real_fuel, 0 pure, 0 binary, 0 multi)
    - Freezing point (freezing_point_C): 103 total (78 real_fuel, 25 pure, 0 binary, 0 multi)
    - DCN (Derived Cetane Number): 426 total (213 binary, 133 multi, 80 pure, 0 real_fuel)
    - MON (Motor Octane Number): 373 total (267 pure, 76 multi, 30 binary, 0 real_fuel)
    - Real-fuel auxiliary properties: distillation cut points T0-T100 (90-108 fuels), speed of sound (170), conductivity (171), smoke point (18), thermal stability (53).
  - No modeling run (§12.1 excluded).

- T31 Full audit run and bundle creation (2026-10-30): DONE
  - Script: `audit.py --bundle` executed with project venv Python.
  - Output: `AUDIT.md` refreshed and `audit_bundle.zip` generated.
  - Audit results: FAIL: 0 | MISSING: 0 | WARN: 8 | OK: 121 (Total: 129 checks).
  - All 8 warnings documented: OneDrive path notice, 1 macOS shadow file, 3 interim file inspect read notices, 2 Antoine missing / 20 non-covering WebBook entries (known to require refit), and 2 historical 0.0 MB DTIC placeholders.
  - No modeling run (§12.1 excluded).

- T26 GCV derived fallback for 100% coverage across unified compounds (2026-10-30): DONE
  - Verified `data/processed/gcv_MJ_kg.csv` generated by T24 against `data/unified/compounds.csv` (1,038 compounds).
  - Initial state: 402 compounds covered (114 stoichiometric tier A, 288 modified Dulong tier C); 636 compounds missing.
  - Applied modified Dulong elemental composition fallback formula: `GCV (MJ/kg) = 33.8*w_C + 144.2*(w_H - w_O/8) + 9.4*w_S` using elemental weights from molecular formulas.
  - Script: `fill_gcv_fallback.py`. Added 636 tier C fallback entries (`tier='C'`, `derived=True`, `method='modified_Dulong'`).
  - Output: `data/processed/gcv_MJ_kg.csv` now contains exactly 1,038 rows (100% compound coverage; 114 tier A, 924 tier C).

- T27 Gross vs Net Calorific Value for AJFTD real fuels (2026-10-30): DONE
  - Extracted Net Heat of Combustion (`net_heat_of_combustion`) and hydrogen content (`hydrogen_content`) directly from AJFTD detail records (`debug/detail_*.html`).
  - Script: `build_heat_of_combustion_real_fuels.py`.
  - Processed and standardized:
    - 513 unique real fuel batches with valid Net Heat of Combustion (NCV) measurements (median NCV = 43.23 MJ/kg, range 40.75 - 46.20 MJ/kg).
    - Measured hydrogen content (`w_H_wt_pct`) populated for 177 fuels (ASTM D3343, D3701, D7171, GCxGC; median 14.53 wt% H).
    - Default jet fuel hydrogen content fallback (13.80 wt% H, `derived_w_H=True`) applied for the remaining 336 fuels lacking hydrogen tests.
    - Gross Calorific Value (GCV) computed via standard conversion: `GCV = NCV + 0.2122 * w_H_wt_pct` (median GCV = 46.16 MJ/kg, range 43.59 - 49.32 MJ/kg).
  - Output: `data/processed/heat_of_combustion_real_fuels.csv` (513 rows; columns: `fuel_id`, `fuel_label`, `ncv_MJ_kg`, `gcv_MJ_kg`, `w_H_wt_pct`, `derived_w_H`, `derived_gcv`, `n_ncv_obs`, `ncv_method`, `w_H_method`, `test_date`, `source`).
  - No modeling run; build_unified.py untouched.



- T08 PubChem Experimental Properties fetch (fetch_pubchem.py): DONE
  - Resolved 1,037 of 1,038 pure hydrocarbon compounds in `data/unified/compounds.csv` to PubChem CIDs via PUG REST (cached in `data/raw/pubchem/inchikey_cid_map.json`).
  - Cached 1,031 unique compound PUG-View JSON records in `data/raw/pubchem/json/<CID>.json` under NCBI rate limits (0.25s delay, <= 4 req/s).
  - Extracted 2,042 standardized experimental property observations across 368 unique hydrocarbons in `data/raw/pubchem/pubchem_experimental.csv` (Tier B):
    - `vapor_pressure_kPa`: 845 records across 329 compounds (standardized to kPa; T_C recorded where reported).
    - `boiling_point_C`: 667 records across 242 compounds (standardized to °C; P_kPa recorded where reported).
    - `flash_point_C`: 467 records across 162 compounds (standardized to °C; closed_cup / open_cup method captured).
    - `heat_of_combustion_MJ_kg`: 63 records across 63 compounds (converted via MW from kJ/mol, J/kg, BTU/lb; range 38.0–55.5 MJ/kg, median 44.42 MJ/kg).
  - Captured 48 qualitative/unparseable strings in `data/raw/pubchem/pubchem_experimental_rejects.csv` (e.g. gaseous descriptors, sublimation, out-of-range values).
  - Created test suite `tests/test_pubchem.py` verifying regex extraction, unit conversions, and output schema.
  - No ML modeling executed.

- T16 Octane numbers (RON/MON) for pure hydrocarbons (2026-10-30): DONE
  - Sources: vom Lehn et al. 2021 (Appl. Energy Combust. Sci. 5:100018), CombustDB (Kessler 2019/2021), and Leenhouts et al. 2026 / Kuzhagaliyeva et al. 2022 (Comm. Chem. 5:114).
  - Script: `build_octane_hc.py`
  - Output: `data/processed/octane_ron_mon_hc.csv` (446 pure hydrocarbon compounds).
  - Schema: `cas`, `compound_name`, `formula`, `smiles`, `inchikey`, `ron`, `mon`, `octane_sensitivity`, `source`, `tier`, `derived`, `is_experimental`.
  - Tiers & counts:
    - Tier A (Experimental): 361 compounds (vom_Lehn_2021_exp: 198, CombustDB_exp: 122, Leenhouts_2026_Kuzhagaliyeva_2022: 39, vom_Lehn_2021 + CombustDB: 2).
    - Tier C (Predicted): 85 compounds (vom Lehn 2021 QSPR-ANN model predictions; `derived=True`, `is_experimental=False`).
  - Values: RON populated for 409 compounds; MON populated for 446 compounds; Octane Sensitivity (RON - MON) calculated for 409 compounds.
  - Sanity check: Isooctane (CAS 540-84-1) RON=100.0, MON=100.0; n-Heptane (CAS 142-82-5) RON=0.0, MON=0.0.
  - No modeling run (§12.1 excluded); all entries verified pure hydrocarbons.

- T17 Freezing point dataset for pure hydrocarbons (2026-10-30): DONE
  - Sources: CRC Handbook of Chemistry and Physics (Organics table; Haynes, Bruno, Lide 2014) and Bradley Open Notebook Melting Point Dataset (Open_NTBKM; Figshare 1031637).
  - Script: `build_freezing_point_hc.py`
  - Output: `data/processed/freezing_point_pure_hc.csv` (905 pure hydrocarbon compounds).
  - Schema: `cas`, `compound_name`, `formula`, `smiles`, `inchikey`, `freezing_point_C`, `freezing_point_K`, `unit`, `source`, `tier`, `derived`.
  - Tiers & counts:
    - Tier A (CRC Handbook of Chemistry & Physics): 668 compounds.
    - Tier B (Bradley Open Notebook Melting Point Dataset): 237 compounds.
  - Values: canonical unit °C (`freezing_point_C`); range -187.75 °C to +437.3 °C (methane at -182.5 °C; benzene at 5.5 °C).
  - No modeling run (§12.1 excluded); all entries verified pure hydrocarbons.

- T18 Cloud point dataset evaluation (2026-10-30): DONE / AUDITED
  - Evaluated existing interim, raw, processed files and libraries (`chemicals`, `thermo`, AJFTD, CRC, Leenhouts).
  - Finding: Cloud point is a mixture/fuel property (wax appearance in multi-component hydrocarbon systems), not defined for pure single components.
  - Experimental data availability:
    - AJFTD: 0 cloud point entries (AJFTD certification protocol measures freezing point / pour point bounds, but not cloud point).
    - CRC 647 / AV-27-18: 0 cloud point entries (freezing point and pour point tabulated).
    - CombustDB contains 265 predicted cloud point rows (`COMBUSTDB_PREDICTIONS`), but zero experimental cloud point rows.
    - Literature sources identified in DATA_ACQUISITION_TASKS.md (T28 Krishna 1989 Table II: 24 diesel blends; T30 NASA 1986; T31 Kirgina 2019) require manual transcription from paywalled/external PDFs.
  - Documented status: No raw experimental cloud point tables currently exist in repo without external PDF transcription.

- T20 Surface tension for pure hydrocarbons (2026-10-30): DONE
  - Sources: Compiled parameters from Jasper (1972) / Lange Handbook, Somayajulu (1988), VDI PPDS-11 critical temperature scaling, and Mulero & Cachadina (2012) via `chemicals.interface` / `thermo`.
  - Script: `build_surface_tension_hc.py`
  - Output: `data/processed/surface_tension_hc.csv` (131 pure hydrocarbon compounds).
  - Schema: `cas`, `compound_name`, `formula`, `smiles`, `inchikey`, `surface_tension_20C_mN_m`, `unit`, `T_ref_C`, `eval_method`, `is_liquid_phase_at_20C`, `Tc_K`, `jasper_a`, `jasper_b`, `Tmin_K`, `Tmax_K`, `sources`, `tier`, `derived`.
  - Canonical unit: mN/m at standard reference temperature T = 20.0 °C (293.15 K).
  - Tiers & methods: Tier A experimental correlations (Jasper: 83, VDI PPDS: 39, Somayajulu: 8, Refprop: 1); includes subcritical vs supercritical flag based on Tc.
  - Sanity check: n-Octane (CAS 111-65-9) = 21.63 mN/m; Benzene (CAS 71-43-2) = 28.88 mN/m; Ethylbenzene (CAS 100-41-4) = 29.29 mN/m.
  - No modeling run (§12.1 excluded); all entries verified pure hydrocarbons.

- T23 Real fuel sooting and smoke point metrics (2026-10-30): DONE
  - Source: Scraped AJFTD FAA ASCENT database cross-referenced with total aromatics and diaromatics.
  - Script: `build_sooting_real_fuels.py`
  - Output: `data/processed/sooting_real_fuels.csv` (100 unique real aviation fuels).
  - Schema: `fuel_id`, `fuel_label`, `property`, `smoke_point_mm`, `unit`, `method`, `test_date`, `aromatics_vol_pct`, `diaromatics_vol_pct`, `source`, `tier`, `derived`.
  - Canonical unit: mm (ASTM D1322 smoke point); range 19.3 mm to 25.0 mm (mean 20.8 mm).
  - Covariates: Total aromatics (% vol) for 100/100 fuels (8.6 to 22.6 vol%); Diaromatics (% vol) for 99/100 fuels (0.1 to 3.4 vol%).
  - Verification: `tests/test_new_acquisitions.py` passes 100% across all 4 datasets.
  - No ML modeling executed (§12.1 excluded).


- T30 & T31 Low-temperature flow properties of real fuels and blends (2026-10-30): DONE
  - Sources:
    1. CRC Aviation Fuels (CRC Report 647 & AV-11-09): 117 records (60 freezing point, 57 pour point).
    2. AJFTD (FAA ASCENT database): 18 records (18 freezing point via ASTM D7153).
    3. NASA NTRS 19860018789 (Cotterell et al. 1986): 34 records across 8 fuels (LFP-3, LFP-5, Fuel #7, LFPA-3, Shale II JP-5, NAPC-2, NAPC-3, NAPC-5), including ASTM melting/freeze points, DSC melting points, static equilibrium melting points, equilibrium pour points, and wt% solid precipitated at the pour point alongside total n-paraffin content (wt%). Full text downloaded and archived at data/raw/pdfs/NASA_CR_179472_19860018789.txt.
    4. Kirgina et al. 2019/2021 (Pet. Coal 61(1) 120-127 & Oil Gas Sci. Technol. 76, 40): 30 records across 5 straight-run diesel fuel samples pure and additized with cold flow improver (Additive C), containing Cloud Point (Tcp, GOST 5066/ASTM D2500), Cold Filter Plugging Point (CFPP, EN 116/ASTM D6371), and Pour Point (Tpp, GOST 20287/ASTM D97) along with total n-paraffin content (% wt).
  - Script: build_low_temp_flow.py
  - Output: data/processed/low_temp_flow_properties_fuels.csv (199 total records).
  - Schema: fuel_id, fuel_name, source, property, value, unit, qualifier, method, wax_content_pct, n_paraffin_content_pct, tier, derived, note.
  - Properties standardized: freezing_point_C (104 rows), pour_point_C (75 rows), cloud_point_C (10 rows), cfpp_C (10 rows). All units in C (ranging from -78.0 C to -3.0 C).
  - Verification: tests/test_low_temp_flow.py passes 100% (validates bounds -100 C to +20 C, non-empty, schema, and source coverage).
  - No ML modeling executed (Section 12.1 excluded).

- T19, T20, & T42 Threshold Sooting Index (TSI) and smoke point reference table (2026-10-30): DONE
  - Sources:
    1. Boehm et al. 2022 (Energy & Fuels 36:1916, doi:10.1021/acs.energyfuels.1c03794; ASCENT Project 065a/066): open access report downloaded to data/raw/sooting/boehm2022.pdf. Extracted conversion equation constants: TSI = a + b * (MW / SP) with a = -3.40, b = 3.727 (mm*mol/g). Anchors: methylcyclohexane = 5.0 (SP = 43.56 mm), 1-methylnaphthalene = 100.0 (SP = 5.13 mm). Referee limits: 1,3,5-trimethylbenzene (TSI 54.7-63.1, nominal 58.9) and 40/60 v/v toluene/isooctane (TSI 21.0-25.4, nominal 23.2).
    2. Mensch et al. 2010 (Combust. Flame 157:1097, doi:10.1016/j.combustflame.2010.01.006): recommended pure-compound TSIs for jet fuel surrogates. Apparatus constants: a = -4.75, b = 4.051 (mm*mol/g).
    3. Li & Sunderland 2012 (Combust. Sci. Technol. 184:829) / ASTM D1322: high-precision smoke point measurements for additional hydrocarbons, alkenes, aromatics, and bounded alkanes (>50 mm).
    4. YSI Database Volume 2 (T21, 447 compounds) mapped to TSI via monotonic piecewise anchor-calibrated mapping.
  - Script: build_tsi_table.py
  - Outputs:
    1. data/processed/tsi_table.csv (203 pure hydrocarbon records overlapping data/unified/compounds.csv):
       - Schema: inchikey, smiles, compound_name, tsi, ysi, smoke_point_mm, source, tier, derived, censored, note.
       - Tiers: Tier A (primary experimental: 28 rows), Tier B (YSI-converted: 175 rows).
       - Censored entries: 14 records with censored=True (smoke point capped at ASTM D1322 candle lamp limit of 50.0 mm).
       - Anchors & references verified: methylcyclohexane (TSI = 5.0), 1-methylnaphthalene (TSI = 100.0), toluene (TSI = 39.0), benzene (TSI = 28.0), n-hexane (TSI = 2.5), n-heptane (TSI = 2.8), iso-octane (TSI = 7.7).
    2. data/processed/tsi_constants.json:
       - Standard conversion parameters and referee fuel criteria for baselines_smoke.py.
  - Verification: tests/test_tsi_smoke.py passes 100% (validates anchors, exact required column schema, reasonable TSI ranges, and JSON configuration).
  - No ML modeling executed (Section 12.1 excluded).

- T34 NIST Jet-A / JP-8 / S-8 compressed liquid density dataset (2026-10-30): DONE
  - Source: NISTIR 6659 ("Thermophysical Properties of Hydrocarbon Fuels", Stephanie L. Outcalt, Arno Laesecke, Mark B. Freund, 2009; Energy & Fuels 23(3):1626-1633), Tables 20–24.
  - Script: `extract_nist_density.py`.
  - Output: `data/processed/nist_density_compressed_liquid.csv` (634 high-precision compressed liquid density observations).
  - Schema: `fuel_id`, `fuel_name`, `category`, `T_K`, `P_MPa`, `P_kPa`, `density_kg_m3`, `uncertainty_kg_m3`, `source_table`, `reference`.
  - Fuels and row counts:
    - `Jet-A-3602` (conventional Jet A, POSF 3602): 135 rows, Table 20. T: 270.0–470.0 K, P: 0.541–32.101 MPa, rho: 681.07–850.67 kg/m³.
    - `Jet-A-3638` (conventional Jet A, POSF 3638): 110 rows, Table 21. T: 270.0–470.0 K, P: 0.922–30.079 MPa, rho: 654.23–827.22 kg/m³.
    - `Jet-A-4658` (conventional Jet A, POSF 4658): 121 rows, Table 22. T: 270.0–470.0 K, P: 0.499–30.015 MPa, rho: 667.89–837.52 kg/m³.
    - `JP-8-3773` (military kerosene, POSF 3773): 143 rows, Table 23. T: 270.0–470.0 K, P: 0.490–40.020 MPa, rho: 654.20–833.20 kg/m³.
    - `S-8` (synthetic paraffinic kerosene, POSF 4734): 125 rows, Table 24. T: 270.0–470.0 K, P: 0.400–30.090 MPa, rho: 611.90–786.20 kg/m³.
  - Experimental method & accuracy: Vibrating-tube densimeter (Anton Paar DMA 512P). Temperature uncertainty 0.03 K, pressure uncertainty 0.01 MPa, density uncertainty 1.0 kg/m³ (~0.12 % relative).
  - Test suite: `tests/test_nist_density.py` (passes 100%).
  - Strictly no ML modeling executed (§12.1 excluded); `build_unified.py` untouched.

