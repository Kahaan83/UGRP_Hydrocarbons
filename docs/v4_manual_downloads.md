# v4 manual-download checklist (user-assisted)

Written 2026-11-10 for v4 Task 8. Nothing here is required: every task runs on the data already in the repo. A document only widens what the report can claim.

## Honesty conventions in this file

- **In repo** = the identifier appears in the plan, the Task 8 brief, `data/MISSING.md`, `docs/`, or a raw file. It is still not proof that the document holds the rows we hope for.
- **UNVERIFIED - confirm before relying** = I could not check it. Treat the title, author and year as a search hint only. I have not invented URLs, DOIs or report numbers. Where none is given, none was available.
- The one scripted check (ThermoML release) is in section 4.
- "Expected filename" is the name the ingester will look for. **No ingester exists yet** for any item below. Each is written only after the file appears (Task 8 rule), so the filename is a contract I will code to, not an existing code path.
- Task numbers 9-19 are those in `docs/superpowers/plans/2026-11-09-v4-integrity-and-upgrade.md`.

## Data roots (canonical paths)

Two raw-data roots exist, checked on disk 2026-11-10:
- **Repo root `data/raw/`** holds only `gc_composition/`, `pour_point/`, `surrogates/`, `git/`. The v4 scripts and loaders read these.
- **`extras/fuelblend_data/data/raw/`** holds everything older (ajftd, flash_mixtures, vozka2019, pires2018, rp1, octane, pdfs, thermoml, ...) and the only `MANIFEST.json`.

Rule used in this file: a *new* document goes in a new folder under repo-root `data/raw/<name>/`. Where an existing extras folder is the natural home (M10, M11, M12), the path is written in full as `extras/fuelblend_data/data/raw/...`. Every path below is complete from the repo root. Plain-language terms: "headline" = the fuel set allowed into headline statistics; "same-fuel D86/D2887 check" = comparing a measured D86 distillation curve with the simulated-distillation (D2887) curve the carbon-number split was built from.

## 1. Summary table

Priority: **A** = closes a named MISSING row or SD-02 and unblocks a headline claim; **B** = widens a claim; **C** = nice to have.

| ID | Document | Unblocks | Folder (`data/raw/...`) | Fallback if you cannot get it | Pri |
|---|---|---|---|---|---|
| M1 | Edwards, *Jet Fuel Properties*, AFRL-RQ-WP-TR-2020-0017 (DTIC AD1093317) | MISSING rows POSF-10325 / 13511 / 13516 / UDRI; Tasks 9, 11, 12 (more jet fuels with raw-file properties) | `afrl_jet_fuel_properties/` | Tasks run on 57 CRC-647 fuels + POSF-10325. No claim about alternative-jet or POSF-13511/13516 | A |
| M2 | ASCENT Project 31a report, April 2019 (property and composition slides) | UDRI_HTA_repfuel, POSF-13511 properties and D86; Tasks 9, 11, 12 | `ascent_31a/` | UDRI_HTA stays property-less (it is the locked jet GC fuel, so locked-set jet metrics stay empty) | A |
| M3 | ATSDR TP-121 Ch. 4 and the API 2010 Jet A property source | JetA_API14avg properties; Tasks 9, 11, 12 | `atsdr_tp121/` | JetA_API14avg remains composition-only and is not scored | B |
| M4 | POSF property sheets (Edwards AIAA 2017-0146, UNVERIFIED) | POSF-10325 freezing/pour/heat; POSF-13511/13516 all properties; D86 for their carbon-number split | `posf_property_sheets/` | POSF-13511/13516 keep `assumed_cn_split` and stay out of headline metrics | A |
| M5 | Potter & Simmons (1998) full GC tables for JP-5 and JP-8 (UNVERIFIED citation) | JP5/JP8 composition (now `excluded_partial_composition`, 41-42 % mass identified); Tasks 9, 11, 12, 16 | `potter_simmons_1998/` | JP-5 and JP-8 stay excluded; no JP-5/JP-8 claim | B |
| M6 | SRM 4572 certificate and Fortin et al. 2012 property table (RP-1) | RP-1 inventory fuels' property targets; Tasks 9, 11, 12 | `rp1_srm4572/` | RP-1 stays distillation-only; no RP-1 property claim | B |
| M7 | AJFTD fuel_id 2 identity check (browser, login) | settles POSF-10325 vs AJFTD "10325"; an independent same-fuel D86 vs D2887 check (the POSF-10325 carbon-number split currently comes from D2887 of the same fuel); Tasks 9, 11, 12 | `ajftd_identity/` | Fuel 2 stays unlinked; the n=1 headline stays single-source | A |
| M8 | Diesel: Kirgina et al. 2019/2021 tables; NASA CR-179472 composition tables; CRC AV-07-02 / FACE diesel (UNVERIFIED) | SD-02 (diesel corpus is empty); Tasks 9, 11, 12, 16, 17 | `diesel_composition/` | **Diesel real-fuel evaluation stays nil; no diesel claim; Task 17 diesel designs reported as unvalidated** | A |
| M9 | FACE gasoline detailed hydrocarbon analysis (CRC AVFL-24 / Sarathy et al., UNVERIFIED) | FACE rows `target_ok=False` (no fractions); Tasks 9, 15, 16, 17 | `face_gasoline_dha/` | Octane refit uses only the 66 bRON blends (21 blendstock groups). No claim about FACE octane from composition | B |
| M10 | Vozka et al. 2019 main article, plus any Purdue GCxGC jet studies (UNVERIFIED) | blend density/viscosity/flash for 25 Jet A blends; Tasks 9, 11, 12 | `extras/fuelblend_data/data/raw/vozka2019/` (exists) | Only freezing and heat of combustion (SI) are usable; no density, viscosity or flash claim on these blends | B |
| M11 | Liaw et al. and other measured ternary/quaternary flash tables (JCED 2017, Korean Gas Society 2016) | more multi-component HC flash rows; Task 9 | `extras/fuelblend_data/data/raw/flash_mixtures/` (exists) | Task 9 evaluates the 126 mixtures (123 binary, 3 ternary) plus the 9 Ha & Lee ternaries. Claim: binary only | B |
| M12 | AJFTD full database with the filters cleared | "of N entries" check; the >160 GCxGC fuel records; Tasks 9, 11, 14 | `extras/fuelblend_data/data/raw/ajftd/` | The 81 scraped fuels only; no composition-based AJFTD claim | B |
| M13 | SLE sources: ACS Omega 2022 quaternary wax data (open); Coutinho 1995/1996 (paywalled) | Task 13 (and the Task 12 freezing margin) | `sle_papers/` | Task 13 runs on the ThermoML SLE records only | C |
| M14 | JPCRD 29, 1 (2000) n-alkane vapour pressure to C36 | Antoine/Cox for C21-C24, whose refit status is disputed (see M14 body); Task 9 Psat | `psat_alkanes/` | C21-C24 keep the `check` flag; Task 9 reports the heavy-alkane Psat as low confidence | C |
| M15 | Newer ThermoML release | Task 13 | n/a | None exists (checked, section 4) | closed |
| M16 | No source identified in repo (MISSING.md fuels listed in M16) | nothing specific; Tasks 9, 11, 12, 16 | `unsourced_fuels/` if you find any | They stay dropped; claims cover only the fuels in `data/real_fuels/v4/` | C |

## 2. Entries

Format: **contains / unblocks / folder + expected filename / identifiers / fallback**.

### M1. AFRL-RQ-WP-TR-2020-0017, Edwards, *Jet Fuel Properties*
- **Contains (per the brief and leads doc, SNIPPET only, not opened):** property tables for reference fuels A-1/A-2/A-3 (A-2 as POSF 10325 is **UNVERIFIED**) and alternative fuels: viscosity, flash, freezing. Whether it lists POSF-13511 or POSF-13516 is **UNVERIFIED**.
- **Unblocks:** MISSING rows for POSF-10325 (freezing, pour, heat of combustion), POSF-13511 and POSF-13516 (all properties), UDRI_HTA_repfuel. The other dropped POSF ids (POSF_10264, 10289, 11498, 12344, 12345) are named in no MISSING row's suggested document; whether the report lists them is **UNVERIFIED** (M16). Feeds Tasks 9, 11, 12 (`raw_file` targets) and the Task 6 corpus rebuild.
- **Folder / filename:** `data/raw/afrl_jet_fuel_properties/AFRL-RQ-WP-TR-2020-0017.pdf`
- **Identifiers in repo:** the title and report number above (brief); `https://apps.dtic.mil/sti/pdfs/AD1093317.pdf` (leads doc, SNIPPET). Stanford HyChem's A2 page `https://web.stanford.edu/group/haiwanglab/HyChem/fuels/A2_spec.html` has a key-properties table for POSF 10325 (leads doc, SNIPPET). It is a web page, not a source document: save it as `hychem_A2_spec.html` in the same folder.
- **How:** DTIC blocks scripts (HTTP 403, AGENTS.md §7). Download in a normal browser. Delete any file under 100 KB: it is likely an HTML consent page.
- **Fallback:** Tasks 9, 11, 12 run on the 57 CRC-647 fuels plus POSF-10325 (flash, density, viscosity from AV-27-18). Narrowed claims: no result for alternative or synthetic jet fuels, and none for POSF-10325 freezing/pour/heat.

### M2. ASCENT Project 31a public technical report (UDRI, April 2019)
- **Contains:** the source already cited in `data/raw/gc_composition/fuel_gc_compositions.csv` for POSF-10325 (slides 16, 18), POSF-13511 (slide 16) and UDRI_HTA_repfuel (slide 12). The Task 6 MISSING row asks for the "fuel's specification-property table from the same ASCENT 31a report and a D86/D2887 curve". Whether such a table exists in the deck is **UNVERIFIED**.
- **Unblocks:** UDRI_HTA_repfuel properties (it is the locked jet GC fuel from Task 7); POSF-13511 properties and D86 (carbon-number split); Tasks 9, 11, 12.
- **Folder / filename:** `data/raw/ascent_31a/april-2019-ascent-031.pdf`
- **Identifier in repo:** `https://ascent.aero/documents/2019/04/april-2019-ascent-031.pdf` (composition CSV). Not opened by me.
- **Fallback:** UDRI_HTA stays composition-only, so locked-set jet metrics have no scored GC fuel. POSF-13511 stays `assumed_cn_split` and outside headline metrics.

### M3. ATSDR TP-121 Ch. 4 and the API (2010b) Jet A source
- **Contains:** ATSDR Table 4-2 compiles 14 API Jet A samples (this is the JetA_API14avg composition source). The MISSING row asks for the original API fuel-property table and a mean D86 curve. The exact API document title is **not in the repo: UNVERIFIED**. Search "API 2010 Jet A kerosene category assessment" and confirm.
- **Unblocks:** JetA_API14avg density, viscosity, flash, freezing, pour, heat of combustion; Tasks 9, 11, 12.
- **Folder / filename:** `data/raw/atsdr_tp121/tp121-c4.pdf` (ATSDR chapter, URL in repo: `https://www.atsdr.cdc.gov/toxprofiles/tp121-c4.pdf`); the API document as `data/raw/atsdr_tp121/api_2010b_jet_a.pdf`.
- **Fallback:** JetA_API14avg stays composition-only and unscored. The v4 jet claim loses one more fuel.

### M4. POSF property sheets
- **Contains:** MISSING rows name "Edwards, 'Reference Jet Fuels for Combustion Testing', AIAA 2017-0146" and "AFRL 'Alternative Jet Fuels' tables" as an **UNVERIFIED SUGGESTION**. Confirm the title, number and that it lists POSF 10325/13511/13516 before relying on it.
- **Unblocks:** POSF-10325 freezing/pour/heat; POSF-13511/13516 density, viscosity, flash, freezing, pour, heat; their D86 or D2887 curve for the carbon-number split. Tasks 9, 11, 12.
- **Folder / filename:** `data/raw/posf_property_sheets/edwards_2017_aiaa_2017-0146.pdf` (default name; if the document is a different one, keep its own name and tell me which POSF ids it covers).
- **Fallback:** POSF-13511/13516 keep `assumed_cn_split` and no properties, so they are excluded from headline metrics.

### M5. Potter & Simmons (1998), JP-5 / JP-8 GC tables
- **Contains:** the MISSING row says ATSDR TP-121 Tables 3/4 list identified compounds only (JP-5 40.9 %, JP-8 42.3 % of mass). It asks for Potter & Simmons (1998) full GC tables, or a DHA/GCxGC-type analysis of the same two fuels. The journal and title are **not in the repo: UNVERIFIED**.
- **Unblocks:** JP-5/JP-8 composition (currently `excluded_partial_composition`, never normalised) and anything that needs it (Tasks 9, 11, 12, 16).
- **Folder / filename:** `data/raw/potter_simmons_1998/potter_simmons_1998.pdf`
- **Also dropped for incomplete composition (kept fraction below 0.95): BrunoHuber2010_3comp [0.69], Chen2016_JI [0.341], Honnet2009_Aachen [0.80].** These are surrogate definitions from papers; the repo names no source for the missing fractions (the original papers by those authors and years are search hints only, **UNVERIFIED**). Listed in M16.
- **MISSING.md ids:** JP5_PotterSimmons, JP8_PotterSimmons.
- **Fallback:** JP-5 and JP-8 stay excluded. No JP-5/JP-8 claim.

### M6. RP-1: SRM 4572 certificate and Fortin et al. 2012 property table
- **Contains:** MISSING row: `rp1_properties.csv` (NISTIR 6646) describes an "RP-1 original sample" whose identity with SRM 4572 (`NIST_RP1`) is unconfirmed. Fortin et al. 2012 is local (`extras/fuelblend_data/data/raw/rp1/fortin_2012_rp1_blends.pdf`, distillation Tables 1-4 only). The certificate or the paper's property table would settle identity and give targets. The certificate's NIST document id is **not in the repo: UNVERIFIED**.
- **Unblocks:** property targets for NIST_RP1, RP1_Binary_2comp, RP1_Edwards_3comp, RP1_NIST_14comp (the surrogates have no measurements, so only NIST_RP1 can ever be scored). Tasks 9, 11, 12.
- **Folder / filename:** `data/raw/rp1_srm4572/srm4572_certificate.pdf` (the Fortin 2012 PDF stays in `extras/fuelblend_data/data/raw/rp1/`)
- **Fallback:** RP-1 stays distillation-only; no RP-1 property claim. The 11 orthogonal blend recipes are confidential (MISSING row dated 2026-11-04) and not obtainable by download.

### M7. AJFTD fuel_id 2 identity check (needs your AJFTD login)
- **What to do:** log in to AJFTD, open the fuel list, find the fuel labelled `10325` (internal fuel_id 2; the scraper stored its label and page title but **not its URL**: `record_url` is empty in `ajftd_all_metrics.csv`), open its View Fuel page and read the fuel name, aliases and any POSF/CRC identifier on the page. Save the page as HTML (File, Save as, "Webpage, HTML only"). **Do not save or share `ajftd_session.json`** (login cookies).
- **Why:** fuel 2 has flash_point, freezing_point, density, viscosity_kinematic, net_heat_of_combustion, `distillation` and `micro_distillation` in the local metrics (checked). If it is POSF-10325 it would supply the missing freezing and heat values and an independent D86 against the D2887 split. Known conflict: its flash (48 C) disagrees with AV-27-18 (55 C), so a match is not guaranteed.
- **Unblocks:** MISSING row POSF-10325; the n=1 headline; Tasks 9, 11, 12.
- **Folder / filename:** `data/raw/ajftd_identity/ajftd_fuel_2_detail.html`, plus one line in chat: "id 2 is / is not POSF 10325, because ...".
- **Fallback:** fuel 2 stays unlinked and the headline stays single-source.

### M8. Diesel (SD-02: the v4 diesel corpus is empty)
The **NASA CR-179472 text is already local** (`extras/fuelblend_data/data/raw/pdfs/NASA_CR_179472_19860018789.txt`), so that source needs a transcription step, not a download. The other two candidate sources need a download. Nothing in the repo shows which has the most usable raw composition tables. MISSING.md fuels this would cover: CRC_FACE_Diesel_1, CRC_FACE_Diesel_2, CRC_FACE_Diesel_3, CRC_FACE_Diesel_4, CRC_FACE_Diesel_5, CRC_FACE_Diesel_6, CRC_FACE_Diesel_7, CRC_FACE_Diesel_8, CRC_FACE_Diesel_9, NASA_LFP3, NASA_LFP5 (no clean composition), and ASTM_D975_ULSD_Reference, SURR_EN590_European_Diesel, SURR_IDEAL_4comp, SURR_Mueller_4comp, SURR_Peugeot_4comp (kept fraction 0.50-0.58, below the 0.95 gate). A source for those last five full compositions is **not identified in the repo (UNVERIFIED)**; they need their own GC/DHA tables.
- **Kirgina et al. 2019/2021**: 10 diesel fuels with pour points and bulk boiling cuts (SD-01). `docs/data_sources.csv` gives doi:10.2516/ogst/2021021 (Oil & Gas Science and Technology). Whether the composition is speciated or only bulk is **UNVERIFIED**; SD-01 says bulk boiling cuts.
- **NASA CR-179472 / NTRS 19860018789**: n-alkane mass fractions and pour points for 8 fuels incl. LFP-3/LFP-5. **The text is already local** (`extras/fuelblend_data/data/raw/pdfs/NASA_CR_179472_19860018789.txt`). This needs transcription into a raw CSV, not a download. Say so and I will do it as a separate task.
- **CRC AV-07-02 / FACE diesel characterization**: from the brief only. The title and content are **UNVERIFIED**; AV-07-02 may be a different report. Search "CRC FACE diesel fuels characterization".
- **Unblocks:** SD-02 (needs at least 5 diesel fuels with raw-file composition), Tasks 9, 11, 12 (diesel rows), 16 (training data), 17 (diesel designs).
- **Folder / filename:** `data/raw/diesel_composition/kirgina_2021.pdf`, `data/raw/diesel_composition/crc_face_diesel.pdf`. After transcription, `data/raw/diesel_composition/diesel_fuels.csv` (columns as `fuel_gc_compositions.csv`).
- **Fallback:** the diesel real-fuel evaluation stays nil. **Claim narrowed:** no diesel result of any kind; the BIS/ASTM diesel designs in Task 17 are reported as unvalidated against real fuels.

### M9. FACE gasoline detailed hydrocarbon analysis
- **Contains:** FACE set I/II CSVs hold component names but no fractions (Task 5), so all 144 FACE rows are `target_ok=False` (`data/real_fuels/v4/gasoline_octane_blends.csv`). A DHA table (CRC "AVFL-24" report and/or a Sarathy et al. paper, both **from the brief only: UNVERIFIED**) would give per-component fractions. Ethanol volume percentages in FACE rows are inferred from the "+N" suffix and have unknown provenance (Task 5 ledger). The gasoline MISSING.md fuels beyond FACE A-J (CoOptima_4comp and its E10/E20/E30/M10/M20 variants, NREL_BOB_01, PRF_E10 to PRF_E90, SURR_TRF_60/70/80/90, TIGAS_Green_Gasoline) have no source identified in the repo and are listed in M16.
- **MISSING.md fuels covered (gasoline, no clean composition):** CRC_FACE_A, CRC_FACE_B, CRC_FACE_C, CRC_FACE_D, CRC_FACE_E, CRC_FACE_F, CRC_FACE_G, CRC_FACE_H, CRC_FACE_I, CRC_FACE_J, CRC_FACE_A_E10, CRC_FACE_A_E25, CRC_FACE_A_E40, CRC_FACE_A_E60, CRC_FACE_C_E10, CRC_FACE_C_E25, CRC_FACE_C_E40, CRC_FACE_C_E60, CRC_FACE_F_E10, CRC_FACE_F_E25, CRC_FACE_F_E40, CRC_FACE_F_E60, CRC_FACE_G_E10, CRC_FACE_G_E25, CRC_FACE_G_E40, CRC_FACE_G_E60, CRC_FACE_I_E10, CRC_FACE_I_E25, CRC_FACE_I_E40, CRC_FACE_I_E60, CRC_FACE_J_E10, CRC_FACE_J_E25, CRC_FACE_J_E40, CRC_FACE_J_E60.
- **Unblocks:** FACE rows as octane targets (Task 15), their use in Tasks 9 and 17 (gasoline), 16.
- **Folder / filename:** `data/raw/face_gasoline_dha/face_dha.pdf` (and any supplementary spreadsheet as `face_dha_supplement.xlsx`).
- **Fallback:** Task 15 uses only the bRON blends (66 blends, 21 blendstock groups, 65 RON / 66 MON targets; `data/real_fuels/v4/gasoline_octane_blends.csv`), and the locked octane set stays thin (12 fuels in 4 groups; `splits/v4/real_fuels_v4_splits.json`). **Claim narrowed:** no octane-from-composition claim for FACE gasoline.

### M10. Vozka et al. 2019 main article and Purdue GCxGC jet studies
- **Contains:** *Energy & Fuels* 33(4) 3275 (2019). Only the SI is ingested (freezing for 25 Jet A blends; heat of combustion for the neat blendstocks; class compositions for 11 blends). The main article has density, viscosity and flash of the blends (leads doc: still needed). The ACS figshare record `https://acs.figshare.com/articles/journal_contribution/Impact_of_Alternative_Fuel_Blending_Components_on_Fuel_Composition_and_Properties_in_Blends_with_Jet_A/7844741` may hold free SI (leads doc, SNIPPET). The brief also names "Purdue GCxGC jet-fuel studies (Vozka, Kilaz et al.)". No further paper is identified in the repo, so any other Purdue paper is **UNVERIFIED**; the HEFA-feedstock companion paper (E&F 2018) is mentioned in the leads doc without a citation.
- **Unblocks:** blend density/viscosity/flash for 11 composition-known blends; Tasks 9, 11, 12.
- **Folder / filename:** `extras/fuelblend_data/data/raw/vozka2019/vozka2019_main_article.pdf` (folder exists; `vozka2019_ef_si.pdf` is already there).
- **Fallback:** freezing and heat of combustion only. No density, viscosity or flash claim on alternative-blend Jet A.

### M11. Measured ternary/quaternary flash tables
- **Contains:** the plan asks Task 9 to evaluate the 126 mixtures (123 binary, 3 ternary); Ha & Lee 2016 added 12 rows (9 ternary). Leads doc candidates (all SNIPPET or paywalled, none opened):
  - Liaw et al., *Combust. Flame* 138 (2004) 308 and similar: mostly model papers, use only measured tables.
  - JCED 2017, doi 10.1021/acs.jced.7b00466 (and the 2016 binary doi 10.1021/acs.jced.6b00516): ternary n-butylcyclohexane + n-C12/n-C16 + aromatic. The flash values are **not** in ThermoML (0 of 11,923 XML files have a flash property), so they must be typed from the paper.
  - *J. Korean Gas Society* (2016), n-nonane + n-decane + n-dodecane: `https://koreascience.or.kr/article/JAKO201610364969655.page`.
- **Unblocks:** a ternary-or-higher claim in Task 9 (LOSO needs several systems).
- **Folder / filename:** `extras/fuelblend_data/data/raw/flash_mixtures/<first_author>_<year>.pdf` (e.g. `liaw_2004_combust_flame.pdf`)
- **Fallback:** Task 9 claims binary mixtures only; ternary flash is reported as 12 rows, "not validated at scale".

### M12. AJFTD full database (filters cleared)
- **What to do:** open the AJFTD list page with every filter cleared, read the "Showing x to y of N entries" text, and read `/about/database-access` for the terms before any re-scrape (leads doc). We hold 81 fuels (ids 1 and 365-462). The About page claims over 160 GCxGC fuel records (CONFIRMED in leads doc; "records" may mean test rows).
- **Unblocks:** more real fuels with composition (Tasks 9, 11, 12, 14 smoke point). Whether the GCxGC data is exposed per fuel is **UNVERIFIED**.
- **Folder / filename:** `extras/fuelblend_data/data/raw/ajftd/` (existing); tell me N in chat.
- **Fallback:** the 81 scraped fuels. No class-composition claim from AJFTD.

### M13. SLE papers (Task 13)
- **Contains:** Coutinho et al. 1995 (FPE, about 60 binary alkane SLE systems) and Coutinho & Stenby 1996 (IECR) are paywalled (leads doc, SNIPPET). The open ACS Omega 2022 paper on quaternary n-C11/C14/C16/C18 wax **disappearance** temperatures is at `https://rcastoragev2.blob.core.windows.net/fac5385c3077144472cce8760da0bd29/ao1c07072.PMC9134231.pdf` (leads doc, SNIPPET). WDT is not WAT or cloud point.
- **Folder / filename:** `data/raw/sle_papers/acs_omega_2022_quaternary_wdt.pdf`
- **Fallback:** Task 13 uses ThermoML SLE records only.

### M14. JPCRD 29, 1 (2000) n-alkane vapour pressure
- **Contains:** the Psat equation for n-alkanes up to C36. Doi 10.1063/1.556054, reprint `https://srd.nist.gov/jpcrdreprint/1.556054.pdf` (leads doc, SNIPPET; likely free). Use for C21-C24. Status conflict: AGENTS.md section 8.5 lists C21-C23 as ok and C20/C24 as `check`, while the leads doc calls C21-C24 `check`. Confirm in `data/processed/antoine_refit.csv` before relying on either.
- **Folder / filename:** `data/raw/psat_alkanes/1.556054.pdf`
- **Fallback:** C21-C24 keep their current refit status.

### M16. MISSING.md fuels with no source identified in the repo
No document is named for these; they are listed so none is orphaned. If you know a source, drop it in `data/raw/unsourced_fuels/` and tell me which fuels it covers.
- Gasoline (21): CoOptima_4comp, CoOptima_4comp_E10, CoOptima_4comp_E20, CoOptima_4comp_E30, CoOptima_4comp_M10, CoOptima_4comp_M20, NREL_BOB_01, PRF_E10, PRF_E20, PRF_E30, PRF_E40, PRF_E50, PRF_E60, PRF_E70, PRF_E80, PRF_E90, SURR_TRF_60, SURR_TRF_70, SURR_TRF_80, SURR_TRF_90, TIGAS_Green_Gasoline.
- Jet: POSF_10264, POSF_10289, POSF_11498, POSF_12344, POSF_12345, Dooley2010_1stGen, Dooley2012_2ndGen, KimVioli2014_UM1, KimVioli2014_UM2, Violi2002_JP8_6comp, NASA_LFPA3, NASA_Shale_JP5, BrunoHuber2010_3comp, Chen2016_JI, Honnet2009_Aachen, HyChem_NJFCP_A1, HyChem_NJFCP_A2, HyChem_NJFCP_A3 (only the A2 specification web page is named in the repo, see M1).
- Fallback: these stay dropped. The v4 claims cover only the fuels in `data/real_fuels/v4/`.

## 3. Cross-reference: every `data/MISSING.md` row maps to an entry or is marked out of scope

| MISSING.md row (date) | Checklist entry |
|---|---|
| NREL Cetane Compendium, 56 unmeasured isomers (2026-11-04) | Not in the repo as a raw file; out of scope for Tasks 9-19 (cetane/DCN is not a target). A processed copy is at `data/processed/v2/cetane_nrel_compendium.csv`. No entry |
| 8 unverified Antoine sets (2026-11-04) | No download exists (NIST WebBook tables do not exist for these compounds). Not blocking. Affects Task 9 only via M14 (C21-C24). No entry needed |
| RP-1 11-blend speciated compositions (2026-11-04) | M6 (confidential; not obtainable by download). Fallback in M6 |
| v4 Task 3 dropped: diesel (5 incomplete + 11 no clean) | M8 (all 16 named) |
| v4 Task 3 dropped: gasoline (55 no clean composition) | CRC_FACE_* : M9. The 21 others: M16 |
| v4 Task 3 dropped: jet incomplete (BrunoHuber2010, Chen2016, Honnet2009) | M5 note and M16 (no source identified) |
| v4 Task 3 dropped: jet no clean composition (18): POSF_10264/10289/10325/11498/12344/12345/13511/13516, HyChem A1-A3, Dooley 2010/2012, Kim-Violi 2014 UM1/UM2, NASA_LFPA3, NASA_Shale_JP5, Violi2002 | POSF_10325/13511/13516: M1, M2, M4. POSF_10264/10289/11498/12344/12345, Dooley, Kim-Violi, Violi, NASA_LFPA3, NASA_Shale_JP5: M16. HyChem A1-A3: M16 (A2 web page also in M1) |
| Task 6: JetA_API14avg | M3 |
| Task 6: UDRI_HTA_repfuel | M2 |
| Task 6: POSF-10325 | M4, M1, M7 |
| Task 6: POSF-13511 | M4, M2, M1 |
| Task 6: POSF-13516 | M4, M1 |
| Task 6: JP5/JP8 Potter-Simmons | M5 |
| Task 6: inventory_v4 RP-1 fuels | M6 |
| SD-02 diesel gap | M8 |
| Task 6 carry-forward: AJFTD id 2 vs POSF-10325 | M7 |
| Task 7 carry-forward: octane locked set thin | M9 |

## 4. ThermoML release check (scripted, 2026-11-10)

- **Method:** one GET with Python `requests` (20 s timeout) to the public NIST Midas record API: `https://data.nist.gov/rmm/records?@id=ark:/88434/mds2-2422` (read only; no download, nothing executed).
- **Result:** HTTP 200, JSON. Record title `ThermoML/Data Archive`, version `1.2.6`, `modified` 2021-03-05, `issued` 2021-08-05. Files listed: `ThermoML.v2020-09-30.tgz` (189,433,115 bytes, the file we already hold), its `.sha256`, `ThermoML.xsd` and `ThermoML.xsd.sha256`.
- **Conclusion:** no release newer than 2020-09-30 is listed on data.nist.gov for this record. This matches the leads doc finding. Limit of the check: it covers only this one record. A newer archive under a different record id, or on trc.nist.gov (not queried), would not be seen. If you know of one, tell me. Otherwise treat Task 13 as ThermoML-2020 only.

## 5. How to tell me a file has arrived

1. Save the file in the folder shown, with the exact filename shown (create the folder if needed). Keep the original download; do not edit it.
2. Say in chat: "M4 arrived", plus anything the file's title page does not make obvious (which POSF ids, which fuels). For M7 give the one-line identity answer.
3. I will then: hash it, record it in `extras/fuelblend_data/data/raw/MANIFEST.json`, open it and check the identifiers above against its title page (clearing the UNVERIFIED flags), transcribe the tables into a raw CSV under the same folder with page references, and write the ingester plus a test that uses a fixture mirroring the real file.
4. Only after that do I rerun the affected task and update the matching `data/MISSING.md` rows. If the document turns out not to hold the rows expected, I will say so rather than force a mapping.
5. Do not drop licensed or login-gated material into a folder that is committed, and never share `ajftd_session.json`. Tell me if a file's licence restricts redistribution.
