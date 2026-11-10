# PROJECT HANDOFF — UGRP Hydrocarbons / fuelblend

Date of last update: 2026-11-09
Written for: the next coding agent taking over this project.
Scope of this document: complete, factual record of what exists, how it got there, and how to operate it. It contains **no suggestions, no recommendations, and no future plans** — only state.

---

## 1. What the project is

Predict **properties of hydrocarbon fuel blends** (jet / diesel / rocket kerosene such as RP-1) from blend **composition**, with a model trustworthy on real measured blends, then use it for **inverse design** (find a blend meeting property targets).

Property priority order (as established at project start):

| Property | Role | Notes |
|---|---|---|
| Kinematic viscosity (cSt, any T) | core | Refutas blending is the physics baseline |
| Flash point (°C) | core | Le Chatelier + Raoult + Antoine baseline; one-sided safety error matters (over-predicting is the dangerous direction) |
| Pour point (°C) | core | PPBI baseline; weakest data and weakest physics |
| Density | auxiliary | needed to convert mass↔volume fractions and dynamic↔kinematic viscosity |
| Freezing/cloud point, distillation (D86), smoke point, DCN/MON, heat of combustion | extension | partly available |

Background: conference abstract (Shah et al., ChemCon 2026) + `Fuel_Blend_Property_Methods_final.docx` (both in `extras\fuelblend_data` era material) described a pure-component stage feeding a feed-forward network over ~40 component inputs trained on **synthetic** blends labelled by the blending rule itself, reporting ~0.4 % error on held-out synthetic blends and ~3–5 % on a real RP-1 check. The methodological conclusion reached in an early session: the 0.4 % number is **circular** (the net learns the labelling rule). A robust model needs (a) real measured blend data, (b) physics baseline + learned residual, (c) splits holding out whole blend *systems*, (d) calibrated uncertainty. The project moved to real data + physics baselines; the v2/v3 iterations (`docs\FULL_REPORT_v2.md`, `docs\FULL_REPORT_v3.md`) are the current result layers.

---

## 2. Machine and environment

- Windows, PowerShell. Project root: `C:\Users\kahaa\OneDrive\Desktop\UGRP_hydrocarbons` (**inside OneDrive** — deliberate; deep zips can hit Windows' 260-char limit, `git config core.longpaths true` is set in this repo).
- Python: root `.venv` (Python 3.12.x). All pipeline runs in this session used `.venv\Scripts\python.exe`.
- Pinned dependencies (`requirements.txt` at root): `requests>=2.31`, `pandas>=2.0`, `rdkit>=2023.9`, `lxml>=4.9`, `scipy>=1.10`, `beautifulsoup4>=4.12`, `playwright>=1.40`, `chemicals>=1.5`, `thermo>=0.6`, `openpyxl>=3.1`. The venv additionally carries `scikit-learn` (GBR/HistGradientBoosting, used by `api/predict_blend.py`, `mixture/t8_vs_hgb.py`), `torch` (used by `mixture/bilinear_pairwise.py`), `joblib`, `pytest` — **these are not listed in requirements.txt** even though the code imports them.
- pandas 3 gotcha encountered: assigning mixed-type lists via `.loc[k, [cols]] = [...]` fails → assign per column.
- Git 2.52; user identity configured globally: `Kahaan83 / kahaan83@github.com`.
- University CS-lab compute is available for model training (stated by user).

---

## 3. Full directory structure (current, exact)

```
C:\Users\kahaa\OneDrive\Desktop\UGRP_hydrocarbons\        <- git repo root (branch: main, origin: github.com/Kahaan83/UGRP_Hydrocarbons)
│   AGENTS.md                  <- previous handoff (v1 / data-collection era; still accurate for that era)
│   PROJECT_HANDOFF.md         <- this document
│   .gitignore                 <- root .gitignore (see §10)
│   compounds.csv              <- canonical compound registry (root level)
│   fuels.csv
│   Makefile
│   pytest.ini
│   requirements.txt
│   ajftd_data.csv             <- scraper-exported AJFTD metrics (untracked)
│   ajftd_data.json            <- scraper-exported AJFTD metrics (untracked)
│
├── src\fuelblend\             <- the library (47 .py incl. __init__.py files)
│   ├── __init__.py
│   ├── api\predict_blend.py   <- public predict(composition, T) wrapper; loads runs\..joblib bundles, lazily retrains mixture residuals if absent
│   ├── baselines\             <- flash.py (Le Chatelier), cold_flow.py (SLE/WAT), distillation.py (API D86), cetane.py, octane.py, heat_of_combustion.py (v1), heat_of_combustion_v2.py (Benson/Joback)
│   ├── data\                  <- model_table.py, validate.py (v3 tests import validate_measurement_df)
│   ├── design\                <- inverse_design_v2.py (hard-coded 60-compound registry, see F1)
│   ├── eval\                  <- consistency.py, metrics.py
│   ├── features\              <- pure_features.py (RDKit descriptors + cache), mixture_features.py
│   ├── mixture\               <- common.py (CompoundRegistry, Refutas ln(eta), Le Chatelier flash), cold_flow.py, fuel_hindcast.py (real-fuel composition vectors), t8_deepsets.py, t8_pairwise.py, t8_vs_hgb.py, bilinear_pairwise.py
│   ├── models\                <- dvpe_activity.py, octane_ghosh_lumps.py (v3 additions)
│   ├── splits\                <- split_generator.py, locked_sets.py
│   ├── thermo\                <- coutinho_wax.py, walther_viscosity.py, vapor_pressure_anchoring.py (v3 additions)
│   ├── train\__init__.py, uq\__init__.py
│
├── scripts\
│   ├── t02_leakage_audit.py   <- only pre-v2 script in the core
│   ├── v2\ (35 scripts)       <- Phase A (audits), B (real fuels), C (property fixes), D (model choice), E (designs); orchestrators reproduce_all.py, run_phase_c.py, run_phase_d.py, run_phase_e.py
│   └── v3\ (14 scripts)       <- ingest_phase1_expansion.py, audit_flash_anchoring.py, audit_thermoml_filters.py, analyze_hoc_propagation.py, build_viscosity_leaderboard.py, hindcast_real_fuels_flash.py, test_freeze_regimes.py, summarize_surrogate_check.py, run_phase2_octane.py, run_phase2_dvpe.py, run_phase2_repairs.py, run_v3_inverse_designs.py, verify_designs_independent.py, reproduce_all.py
│
├── tests\ (57 .py)            <- tests/test_v2_* (39), tests/test_v3_* (11), plus test_bilinear_pairwise, test_predict_blend, test_model_invariances, test_harness, test_validate_schema, test_result_regression, test_migration_integrity
│
├── configs\ (7 yaml)          <- specs_astm.yaml, specs_bis.yaml, specs_diesel.yaml, specs_gasoline.yaml, specs_jet.yaml, uncertainty_v2.yaml, uncertainty_v3.yaml
│
├── splits\ (6 json)           <- harder_splits.json, locked.json, real_fuels_holdout.json, real_fuels_nested_v2.json, s1s2s3.json, v3\real_fuels_v3_splits.json
│
├── data\
│   ├── MISSING.md             <- tracked documents which sources never arrived
│   ├── antoine_sources.csv    <- tracked
│   ├── real_fuels\            <- tracked entirely: inventory.csv, inventory_measured.csv, inventory_reconstructed.csv, compositions.csv, compositions_measured.csv, compositions_reconstructed.csv, blend_links.csv, candidates_v2.csv, v3\inventory_v3.csv, v3\compositions_v3.csv
│   ├── processed\             <- v2\ (9 tracked deliverable tables), v3\purchasable_compounds.csv (tracked), cetane_pure_curated.csv (tracked);
│   │                             UNTRACKED heavy/derived: model_table.csv (9.7 MB), thermoml_hc.csv (10.7 MB tracked, whitelist), measurements_migrated_core.csv, measurements_flash_expanded.csv, measurements_low_temp_flow.csv, measurements_real_fuels_heat_of_combustion.csv, crc_fuels.csv, dcn_hc.csv, distillation_ajftd.csv, excluded_by_R1.csv, fortin2012_rp1_blends.csv, gcv_MJ_kg.csv, rp1_characterization.csv, rp1_surrogates.csv, antoine_refit.csv, compounds_canonical.csv
│   ├── unified\               <- UNTRACKED: measurements.csv, measurements_dedup.csv (16.6 MB), UNIFIED_REPORT.md
│   ├── cache\pure_grid.npz    <- untracked
│   ├── runs\                  <- baseline_preds.csv etc. (untracked)
│   ├── raw\pour_point\        <- crc647_fuel_properties.csv, av2718_fuel_properties.csv (untracked)
│   └── interim\               <- download/parse intermediates (untracked)
│
├── docs\
│   ├── FULL_REPORT_v2.md, FULL_REPORT_v3.md   <- the v2 and v3 result reports (tracked)
│   ├── claims_register.md, model_decision.md, scope_decisions.md, data_sources.csv  <- tracked
│
├── results\
│   ├── t0\ (4 leakage CSVs), t6\inverse_designs_three_families.csv, t8\ (4 CSVs: ablation, metrics, hgb_comparison, deepsets predictions), eval\flash_loso_predictions.csv, hindcast\hindcast_predictions.csv, pairwise\extrapolation.csv   <- v1-era outputs, tracked
│   ├── v2\                    <- 40 tracked deliverables: A1..A9, A12, B6, C1..C8, D1..D6, E3..E6, improvement\ (7 files incl. flash_anchored_loso.csv, viscosity_walther_eval.csv, upgraded_designs.csv, paper_ablation_table.csv)
│   └── v3\                    <- 17 deliverables: S4\coverage_matrix.csv, R1\viscosity_leaderboard.csv, R2\flash_anchoring_audit.csv, R3\real_fuels_flash_hindcast.csv, R4\surrogate_check.csv, R5\freeze_regime_details.csv + freeze_regime_test.csv, R6\filter_drop_breakdown.csv, R6\property_distributions_before_after.csv, R7\hoc_propagation_error.csv, G2\octane_lump_model_eval.csv + blending_octane_values.csv, G6\dvpe_activity_eval.csv, G8\diesel_cetane_cold_flow_eval.csv, G10\thermoml_integration_eval.csv, U6\designs_v3.csv + design_robustness_audit.csv
│
├── runs\
│   ├── pure_models\           <- 20 .joblib bundles (per-property, per-seed) — untracked
│   └── mixture_models\        <- 10 .joblib bundles — untracked
│
├── archive\
│   ├── baseline_2026-11-03_code.zip
│   ├── baseline_manifest.json, v2_manifest.json   <- tracked (checksum manifests pinning baseline state)
│
└── extras\                    <- EVERYTHING non-core, moved here on 2026-11-07
    ├── fuelblend_data\        <- the FULL master repo (has its own .git with complete v1/v2/v3 history, remote origin = same GitHub URL; its HEAD was the pre-restructure main tip b23bbe3)
    │   ├── (its own src/, scripts/t*/v2/v3, tests/ (96 files incl. v1-era tests), configs, splits, data/, docs/ (FULL_REPORT.md v1, PROGRESS_REPORT.md, AUDIT.md, T8_T9_INTERFACE_SPEC.md, PROJECT_HANDOFF.md, DATA_LOG.md, DATA_ACQUISITION_TASKS.md, SOURCE_CHECK.md, CHANGES.md, SHA256SUMS.txt),
    │   ├── data\raw\ (ThermoML tgz 189 MB, Zenodo zips, git clones, PDFs, ajftd/, lfl/, rp1/, pour_point/, flashpoint_li/, surrogates/ etc.), data\interim\, data\processed\, data\unified\, results\ (incl. t0_repro\reference_backup\, rerun_outputs\), runs\, archive\, debug\
    │   └── root scripts: fetch.py, inspect_downloads.py, build_viscosity.py, thermoml_extract.py, webbook_antoine.py (+_extended), refit_antoine.py, build_unified.py, audit.py, ingest_ajftd.py, ingest_crc.py, ingest_dcn.py, ingest_rp1.py, ajftd_parse_html.py, ajftd_scraper.py, baselines.py, baselines_flash.py, common.py, fill_lfl.py, parse_zabetakis_a1.py, fortin_extract.py, build_* (costald, critical_constants, fortin2012_processed, freezing_point_hc, gcv, heat_of_combustion_real_fuels, low_temp_flow, octane_hc, sooting_real_fuels, surface_tension_hc, tsi_table), extract_* (ait, ajftd_properties, critical_constants, nist_density), fetch_* (pubchem, ysi), check_s2.py, + ajftd_session.json
    ├── fuelblend_university_bundle\ + fuelblend_university_bundle.zip   (56 MB, 80 files)
    ├── build_bundle.py        (built the university bundle)
    ├── pre_abstract\          (8.3 MB, 27 files — abstract-era material incl. the two source .docx files)
    ├── reports\               (ANALYSIS_PLAN.md, todo.md… actually reports\ holds WP0.md/plan.md-class notes, tasks\ holds plan.md, todo.md)
    ├── data\                  (4 files: YSI sooting database copies, old root data dir)
    ├── src\                   (3 stale .py copies of src/fuelblend files — verified older than repo versions before moving)
    ├── ajftd_session.json     (SENSITIVE: AJFTD login cookies — never commit or share)
    └── .pytest_cache\, __pycache__\
```

Key facts about the structure:

- The root repo (`main`, commit `759dcbf`, 268 tracked files, ~20.7 MB) tracks: all code (src/scripts/tests), configs, splits, docs, curated data whitelists, v2+v3 results, manifest JSONs. Untracked: heavy data (`data/processed/model_table.csv`, `data/unified/*`, `runs/`, `extras/`, `.venv/`). The GitHub copy therefore can't run the heavy pipeline alone but the local root can — the local untracked files are all physically present.
- `extras\fuelblend_data` is a separate, intact git repo with the full history (baseline tag 2026-11-03, v2, v2.1, v3 tags all point at `main`'s old tip `b23bbe3`). Nothing in `extras/` gets touched by root git operations because the root `.gitignore` excludes `extras/`.
- `ajftd_session.json` (sensitive login cookies) sits in `extras\` — git-ignored both at root and in the inner repo's `.gitignore`.

---

## 4. Chronological record of all work done

### Stage 1 — Data acquisition and automation (early sessions)

1. `fetch.py` written: automated downloads of every machine-obtainable source (git clones, Zenodo, PDFs, ThermoML archive). Features: idempotent/resumable, manifest + sha256 verification, `.part` resume, 429 back-off, per-file deadline, Ctrl-C summary, HTML-consent-page rejection (DTIC blocks scripts with HTTP 403 → lesson learned; two fake 0.0 MB DTIC "PDFs" remain flagged). Flags: `--skip zenodo | zenodo --list | --include REGEX | --max-gb N`.
2. `inspect_downloads.py`: unzips Zenodo `datasets.zip`/`pretrained_models.zip` (skips macOS `__MACOSX/._*` junk — hundreds of fake parse errors avoided), writes `data/interim/INVENTORY.md`, `read_errors.txt`.
3. Downloaded sources (all in `extras\fuelblend_data\data\raw\`):
   - Leenhouts et al. benchmark (Zenodo 19914760): curated pure + mixture flash/viscosity/DCN/MON. `flashpoint_pure_exp.csv` 585 rows; `flashpoint_mix_exp.csv` 1006 rows (123 binary, 3 ternary HC); `viscosity_pure_exp.csv` 3027; `viscosity_mix_exp.csv` 35,938; `dcn_mix_exp.csv` 484; `mon_mix_exp.csv` 292. Warning encoded: `DatasetsFuelProperties.xlsx` = 50,000 SYNTHETIC rows (never ground truth); `PureFuels.xlsx` is its source table.
   - CheMixHub (NeurIPS 2025, github.com/chemcognition-lab/chemixhub): cloned; `nist-logV` used (`processed_NISTlogV.csv`, 239,200 rows; label says "Log viscosity / cP" but the column is **ln(η / Pa·s)**).
   - Bilodeau et al. (Zenodo 8042966): raw NIST binary liquid viscosity inside `pretrained_models.zip`.
   - ThermoML archive (NIST TRC, doi 10.18434/mds2-2422, `ThermoML.v2020-09-30.tgz`, 189,433,115 B, sha256 `231161b5…58d7ec2`): downloaded and verified.
   - NIST WebBook: Antoine constants scraped from `cbook.cgi?ID=C<cas>&Mask=4` (bar, K).
   - AJFTD (Illinois/FAA ASCENT real jet-fuel database, login-gated): 81 fuels scraped via Playwright.
   - NIST RP-1 (NISTIR 6646 + Energy & Fuels 2012 doi 10.1021/ef3004009 + Fortin 2012 NIST copy).
   - CRC reports 647 / AV-27-18 / AV-11-09 (crcao.org PDFs).
   - Zabetakis Bulletin 627 (LFL), Li et al. 2020 (arXiv 2002.11315, 10,575 pure flash points — PDF obtained; figshare 9275210 dataset via `load_flashpoint_li()`), DTIC AFRL-TR-2020-0017 (two 403-blocked fake PDFs), DIPPR (not used, proprietary).

### Stage 2 — Per-source ingestion pipelines (early sessions)

4. `build_viscosity.py`: hydrocarbon-only NIST/TDE tables → `data/processed/viscosity_pure_hc.csv` (4,832 raw rows) and `viscosity_mix_hc.csv` (10,947 binary raw rows, 109 compounds); unit conversion documented (ln(Pa·s) − 6.908 = ln(mPa·s)).
5. `thermoml_extract.py` **rewritten from scratch** after verifying against the official XSD (ThermoPyL `thermoml_schema.py`) and a real record (`tests/data/je8006138.xml`):
   - bug 1: fixed T/P live under `Constraint > ConstraintID > ConstraintType`, not a direct `ConstraintType` child → 61 % of temperatures were being lost.
   - exact property-name whitelist (no "Excess viscosity"), `PropPhaseID/ePropPhase == Liquid` only (crystal vapour pressures dropped), composition variables with `VarPhaseID/ePhase != Liquid` ignored (VLE vapour fractions), ternary last component implied, pressure >1000 kPa dropped (flag `--max-pressure-kpa`), molality etc. rejected, P part of the dedupe key, old-format files (no `prop_phase`) refused.
   - → `data/processed/thermoml_hc.csv` (native ThermoML units Pa·s / m²·s / kg/m³ / kPa; `build_unified.py` converts).
   - A pre-fix run had produced a **contaminated** unified set (−0.063 excess viscosity, Psat up to 1.4×10⁶ kPa, density 0.52 kg/m³ at 873 K, only 39 % rows with T, 1,728 spurious conflicts, 6,451 skipped rows); that run was discarded and re-extracted after the fix.
6. `webbook_antoine.py`: converts WebBook Antoine (bar) to kPa via `A_kPa = A_bar + 2`; writes `antoine_best.csv` / `antoine_all_sets.csv`. Bug fixed: `page_name` returned the site banner; corrected to `<h1 id="Top">` → `<title>`. Result: 22/24 CAS-matched (methylcyclohexane, n-eicosane had no table); 20/22 sets do **not** cover the 40–100 °C window for C10+ alkanes (ranges start ~95–190 °C).
7. `refit_antoine.py`: refits Antoine from measured vapour pressures (ThermoML), kPa-native, per-compound status (`ok` | `check: rms > 0.03` | `skipped: <6 points or span <30 K` | `failed: …`), guard against negative-B starts. One run: 220 compounds → 143 ok, 43 check, 34 skipped; 94 ok-and-cover-window. Heavy odd alkanes missing: C13, C15, C17, C19 no data; C18 skipped (5 points); C11 only 393–459 K.
8. AJFTD extraction: `ajftd_scraper.py` (Playwright; patched: resume, retries, real fuel label from tombstone instead of site banner; pagination cross-check against "Showing x of N"), `ajftd_parse_html.py` (offline re-parse of saved `debug/detail_<id>.html` — 81 pages saved), `ingest_ajftd.py` → `ajftd_fuels.csv` (long numeric) + `ajftd_composition.csv` (wide `comp::<metric>`); API gravity → density at 15.56 °C via `ρ = 141.5/(API+131.5) · 999.016` kg/m³.
   Numeric coverage of the 81 fuels: API gravity 81; distillation T0…T100 ≈80 fuels (561 rows); smoke point 82; aromatics (% vol) 82 (8.6–22.6); diaromatics 80; existent gum 82; sulfur mercaptan 81; thermal stability 81. **Flash point: 1 fuel; kinematic viscosity: 1 fuel; density kg/m³: 1 fuel** (the rest "No Test").
9. LFL: `parse_zabetakis_a1.py` + hand-reading Zabetakis Bulletin 627 App. A-1 OCR → `data/raw/lfl/lfl_measured.csv` (12 hydrocarbons, ratio-gated vs 0.55·Cst); `fill_lfl.py` → 695 compounds (25 measured, rest estimated). Pure-flash physics became runnable on ~84–152/902–923 records; Antoine coverage remained the limiter.
10. `load_flashpoint_li()`: Sun/Krakauer/Morgan 2020 figshare 9275210 — 14,696 rows (flash in K), +550 HC pure flash records, duplicates collapsed to median; suspect values flagged not deleted (e.g. C25 paraffin at 120 °C).
11. `fortin_extract.py` + `build_fortin2012_processed.py`: NIST Fortin 2012 RP-1 blend paper → 252 distillation rows (T0.025…T85-90) for 11 orthogonal RP-1 blends + RP-1 4572 + RP-2 (loader `load_fortin()`). **The 11 blend compositions are not in that paper** (recipes confidential; NISTIR 6646 phase-1 covers one sample only).
12. RP-1 NISTIR 6646 ingestion (`ingest_rp1.py`): one characterized sample (33 components w/ SMILES + composition tiers), 14-component surrogate, 1480 property rows (density ~83–484 kg/m³ incl. high P/T, kinematic viscosity 1.13–7.67 cSt, thermal conductivity).
13. CRC ingest (`ingest_crc.py`): crc647 → 4,559 raw rows (57 aviation fuels: pour 57, freezing 57, flash 57, D86 grids, ±20/−40 °C kinematic viscosity ×171, heat, conductivity, speed of sound); av2718 → 869 rows (aromatics, D2887 distillation, viscosity+density, vapour pressure, gas solubility); av1109 → 3 rows. After `build_unified`: pour_point 57 / freezing 78 / viscosity_kin 361 / thermal_stability 53 in dedup (~70 k raw rows total).
14. Additional builders at repo root (output CSVs listed in §3 data/processed): `build_costald_parameters.py`, `build_critical_constants.py`, `build_gcv.py`, `build_heat_of_combustion_real_fuels.py`, `build_low_temp_flow.py`, `build_octane_hc.py`, `build_sooting_real_fuels.py`, `build_surface_tension_hc.py`, `build_tsi_table.py`, `fill_gcv_fallback.py`, `extract_ait.py`, `extract_critical_constants.py`, `extract_nist_density.py`, `fetch_pubchem.py`, `fetch_ysi.py`.

### Stage 3 — Unified dataset and audit (early sessions)

15. `build_unified.py` merges every source into the canonical schema (§7): `data/unified/measurements.csv` (every observation) and `measurements_dedup.csv` (one per condition). De-dup rules: (1) drop exact repeats (same source+ref+value — ~⅓ of raw NIST rows), (2) group by (property, system, composition, T, P), (3) median, (4) flag (never silently delete) outliers/conflicts. `system_id` = md5[:12] of the sorted component InChIKey set (`fuel:<id>` for real fuels) — the GroupKFold key.
16. Real run numbers (documented in AGENTS.md §8): of 4,341 pure viscosity points with ≥3 independent refs only 13 (0.3 %) deviate >25 % from the group median, but one reference reports benzene@298 K = 0.395 mPa·s vs ~70 refs ≈0.60. A real NIST run: 15,779 raw rows → 11,030 deduplicated records, 30 genuine conflicts.
17. Hydrocarbon-only open-data state at that stage: flash 244 pure + 126 mixtures usable; viscosity 692 pure + 2,752 binary (293–323 K); DCN/MON mixtures 346/106; **≥3-component HC data essentially none** (3 flash ternary rows); pour point none.
18. `audit.py` → AUDIT.md + `audit_bundle.zip` (no raw data inside). Last recorded bundle: FAIL 0 · MISSING 6 · WARN 8 · OK 49.
19. Physics-input coverage: 370 HC flash records, 248 components; Antoine available 52/248, in-range 29; LFL was 0/248 before the Zabetakis work.

### Stage 4 — First physics baselines + residual learning (v1.5, `baselines.py` at master-repo root)

20. Viscosity: physics baseline + residual GBM on 15,281 mixture records (GroupKFold on system_id): ln-linear MAPE 8.8 %, **Refutas 4.8 %** (binary winner), GBM residual helps only multi-component (4.85 % vs 6.03 %; over-corrects binaries — only 319 unique systems). Decision encoded: Refutas for binaries, GBM residual only for n_comp≥3. Report: `BASELINE_VISCOSITY_REPORT.md` (in extras\fuelblend_data).
21. Flash: domain filter (−90..250 °C), on 147 predicted records overall MAE 7.17 °C; binary 1.59 °C / 98 % within 5 °C; pure 10.7 °C (pubchem tail).

### Stage 5 — v2 iteration (scripts/v2; docs/FULL_REPORT_v2.md; tagged `v2-2026-11-05`, `v2.1-2026-11-05`)

v2 was built as a self-contained "paper core" (see `extras\fuelblend_data\FILE_MAP.md`, also in `docs\` of the core as `docs\data_sources.csv`, `scope_decisions.md`, `model_decision.md`). v2 keeps only 60 Python files + curated data; data-collection and v1 scripts were archived. Structure: Phase A (audits, CLM-01..10 claims) → B (real fuels) → C (property fixes) → D (model choice) → E (design).

Task → claim → output mapping (all outputs in `results\v2\`, all verified present with expected row counts):

| Script | Task | Claim | Output |
|---|---|---|---|
| audit_inverse_designs.py | A1 | CLM-01, CLM-02 | A1/design_audit.csv |
| eval_recalibration_loo.py | A2 | CLM-04 | A2/recalibration_loo.csv |
| check_cetane_basis.py | A3 | CLM-05 | A3/cetane_basis_rerun.csv |
| eval_gn_baselines.py | A4 | CLM-06 | A4/gn_baselines.csv |
| recheck_hoc.py | A5 | CLM-07 | A5/hoc_recheck.csv |
| refit_uncertainty_growth.py | A6 | CLM-08 | A6/uncertainty_growth_refit.csv |
| eval_flash_margin_semantics.py | A7 | CLM-09 | A7/flash_margin_semantics.csv |
| audit_flash_data.py | A8 | CLM-10 | A8/flash_usage_audit.csv |
| reproduce_hindcast.py | A9 | CLM-03 | A9/hindcast_*.csv |
| audit_data_redistribution.py | A12 | — | A12/redistribution_audit.csv |
| build_candidates_inventory.py | B1 | — | data/real_fuels/candidates_v2.csv |
| ingest_real_fuels_v2.py | B2-B6 | Sec. 3 | data/real_fuels/*, B6/ |
| assemble_distillation_real_fuels.py | B | — | data/processed/v2/distillation_real_fuels.csv |
| eval_hoc_v2.py | C1 | Sec 2.4 | C1/hoc_eval.csv |
| eval_flash_loso_v2.py | C2 | Sec 2.2 | C2/flash_loso_v2.csv |
| verify_c3_leakage.py | C3 | — | results/t0/leakage_*.csv |
| eval_cold_flow_v2.py | C4 | Sec 2.3 | C4/cold_flow_v2.csv |
| digitize_cetane_compounds.py | C5 | — | data/processed/*cetane* |
| eval_octane_real_fuels.py | C6 | Sec 2.6 | C6/octane_real_fuel.csv |
| eval_d86_prediction.py | C7 | Sec 2.5 | C7/d86_eval.csv |
| eval_dvpe_baseline.py | C8 | Sec 2.6 | C8/dvpe_eval.csv |
| eval_bilinear_hindcast.py | D1 | Sec 4.2 | D1/bilinear_hindcast.csv |
| eval_class_interactions.py | D2 | Sec 4.2 | D2/class_interactions.csv |
| eval_shrinkage_sweep.py | D3 | Sec 4.3 | D3/shrinkage_sweep.csv |
| eval_stacked_model.py | D4 | Sec 4.2 | D4/stacked_eval.csv |
| eval_error_vs_n_components.py | D5 | Sec 4.4 | D5/error_vs_n.csv |
| eval_truncation_sensitivity.py | D6 | Sec 4.5 | D6/truncation.csv |
| run_phase_e.py | E3-E6 | Sec 6 | E3..E6/*.csv |
| eval_coutinho_wax.py | v2.1 T2 | — | improvement/freezing_point_coutinho_eval.csv |
| eval_walther_viscosity.py | v2.1 T3 | — | improvement/viscosity_walther_eval.csv |
| run_inverse_design_upgraded.py | v2.1 T4 | — | improvement/upgraded_designs.csv |
| generate_paper_ablation_table.py | v2.1 T5 | — | improvement/paper_ablation_table.csv |

Concrete v2 numbers pinned by `tests/test_result_regression.py` (exact tolerances): A1 design audit — jet_claimed_text design: wat_C computed −15.2 °C (±0.2), kin_visc@−20C shows SOLID; jet_csv design freezes at 20.1 °C (±0.2). A2 recalibration: Physics (Raw Refutas) overall_rmse_ln = 0.4421; Recalibrated 2-param LOO = 0.1411; lab-only recalibration also present.

v2 result I/O: A8 found the flash dataset contains 5,339 rows of which pure-compound parameter-registry content is N=1 per compound (no mixing physics); D2 fits 15 symmetric class-pair interaction parameters on lab binary/ternary mixtures (N=15,507 rows of model_table) with Ridge L2, applied to the 57 CRC 647 aviation fuels under zero-shot and LOO; R1 (t02_leakage_audit) established **0 CRC rows in model_table S1 training rows** (Rule R1 requires composition_level L2/L3 — real fuels entirely excluded from `model_table.csv`).

Known unresolved inconsistencies carried in the v2 material (FILE_MAP §6, factual):

- F1: `inverse_design_v2.py` carries a hard-coded 60-compound table (density, mp, heat of fusion, bp, viscosities, DCN, RON/MON typed by hand) and does not read compounds.csv; every number in the three headline designs depends on those values.
- F2: the viscosity recalibration `a0=-2.40597, a1=0.70534` in the same file is the in-sample fit on all 57 CRC fuels; the reported 0.1411 is leave-one-fuel-out.
- F3: the report calls the bilinear model 40-parameter; `BilinearPairwiseInteraction`'s U matrix is in_features × rank (default rank 16) + lambdas; which object produced D1 is not recorded.
- F4: Bilinear and Stacked rows in the D1 table are identical to every digit.
- F5: spec limits exist in both `configs/specs_*.yaml` and hard-coded in the design code.
- F6: the master repo README still describes the data-collection stage with T8/T9 open.

### Stage 6 — v3 iteration (scripts/v3; docs/FULL_REPORT_v3.md 227 lines + docs/claims_register.md 22 rows; tagged `v3-2026-11-06`; `docs/model_decision.md` records the model choice)

v3 = "validation-gap pipeline". Orchestration: `scripts/v3/reproduce_all.py` runs 13 steps in order, then the full pytest suite, then deliverable verification, then writes `SHA256SUMS.txt` manifest. Flags: `--verify-only`, `--skip-tests`.

Pipeline steps in order:
1. Phase S: `ingest_phase1_expansion.py` (Task S1-S6, real-fuel expansion + oxygenates) → `data/real_fuels/v3/inventory_v3.csv` (153 rows), `compositions_v3.csv` (3,298 rows), `splits/v3/real_fuels_v3_splits.json` (6 entries), `results/v3/S4/coverage_matrix.csv` (20 rows).
2. `audit_flash_anchoring.py` (R2) → R2/flash_anchoring_audit.csv (5 rows).
3. `audit_thermoml_filters.py` (R6) → R6/filter_drop_breakdown.csv (7 rows) + R6/property_distributions_before_after.csv (reconstructs exact masks from `filter_thermoml_tidy.py`).
4. `analyze_hoc_propagation.py` (R7) → R7/hoc_propagation_error.csv (3 rows).
5. `build_viscosity_leaderboard.py` (R1; imports walther_viscosity, common, cold_flow, fuel_hindcast, bilinear_pairwise) → R1/viscosity_leaderboard.csv (32 rows).
6. `hindcast_real_fuels_flash.py` (R3) → R3/real_fuels_flash_hindcast.csv (70 rows) (+ flash_anchored_loso.csv, flash_loso_predictions.csv, hindcast_predictions.csv from earlier results dirs it reuses).
7. `test_freeze_regimes.py` (R5; imports coutinho_wax) → R5/freeze_regime_details.csv (113 rows) + freeze_regime_test.csv (3 rows).
8. `summarize_surrogate_check.py` (R4) → R4/surrogate_check.csv (14 rows).
9. `run_phase2_octane.py` (G2; imports octane_ghosh_lumps) → G2/octane_lump_model_eval.csv (3 rows) + G2/blending_octane_values.csv.
10. `run_phase2_dvpe.py` (G6; imports dvpe_activity) → G6/dvpe_activity_eval.csv (2 rows).
11. `run_phase2_repairs.py` (G8/G10; imports octane_ghosh_lumps, dvpe_activity, coutinho_wax) → G8/diesel_cetane_cold_flow_eval.csv (16 rows), G10/thermoml_integration_eval.csv (7 rows).
12. `run_v3_inverse_designs.py` (U1-U5) → U6/designs_v3.csv (6 rows; also consumes purchasable_compounds.csv, wax_dauphin_coutinho.csv, fuel_surrogates.csv).
13. `verify_designs_independent.py` (U7) → U6/design_robustness_audit.csv (30 rows); references `verify_designs_independent` results incl. `surrogate_comparison.csv`.

Verification state (run twice in this session, both at the old core dir and at the current root): **all 20 v3 deliverables PASS with expected row counts; all v2 deliverables PASS with expected row counts** (`python scripts/v3/reproduce_all.py --verify-only`, `python scripts/v2/reproduce_all.py --verify-only`; exit 0).

v3 data dependencies resolved by basename copy (all found outside reference_backup): thermoml_hc.csv, thermoml_clean_hc.csv, av2718_fuel_properties.csv, crc647_fuel_properties.csv, blending_octane_values.csv, flash_anchored_loso.csv, fuel_surrogates.csv, purchasable_compounds.csv, wax_dauphin_coutinho.csv, measurements_real_fuels_heat_of_combustion.csv, pure_grid.npz.

### Stage 7 — Core-directory build and restructuring (this session, 2026-11-07)

14. `build_core_dir.py` (received from user; located `C:\Users\kahaa\Downloads\files (3)\build_core_dir.py` + `FILE_MAP.md`): copies only paper-core files to a clean folder — 60-file CORE_PY list, EXTRA_DIRS (`configs`, `splits`, `data/real_fuels`, `data/processed/v2`, `results/v2`, `archive`), EXTRA_FILES (root CSVs, requirements.txt, Makefile, pytest.ini, MISSING.md, antoine_sources.csv, cetane_pure_curated.csv, docs), TEST_GLOBS (v2 tests), plus `src/**/__init__.py`, plus a **string-literal scanner**: every `.csv/.json/.yaml/.yml/.npz/.pt/.joblib` basename quoted in copied code is looked up by basename anywhere in the source repo (skipping `.git/.venv/__pycache__`, and skipping hits in `reference_backup`/`rerun_outputs`) and copied. It prints copied count, fixed-list misses, and unresolved data names. Idempotent (`t.exists()` → skip).
15. Updates applied to the script across this session, in order:
    - Added `"runs"` to EXTRA_DIRS (the 30 pre-trained `.joblib` bundles in `runs/pure_models` + `runs/mixture_models` were missed because `predict_blend.py` builds those paths programmatically; they were copied manually first, then the script was updated).
    - v3 support: CORE_PY += `src/fuelblend/data/model_table.py`, `data/validate.py`, `models/dvpe_activity.py`, `models/octane_ghosh_lumps.py`, `thermo/vapor_pressure_anchoring.py`; EXTRA_DIRS += `results/v3`, `data/processed/v3`; EXTRA_FILES += `docs/FULL_REPORT_v3.md`, `docs/claims_register.md`, `docs/model_decision.md`; TEST_GLOBS += `tests/test_v3_*.py`, `test_validate_schema.py`, `test_result_regression.py`, `test_migration_integrity.py`; copy loop also globs `scripts/v3/*.py` (reason: v3's reproduce_all runs the whole pytest suite, so its tests must ship; v3 tests import only `fuelblend.data.validate` + stdlib/pandas).
    - Bug fix 1: the data scanner previously rglob'd **all** `.py` under dst; once the master repo was moved inside dst, it raked third-party code in `extras` — hundreds of junk references and unintended copies (620-file run), plus…

    - Bug fix 2: crash at end-of-run printing a Cyrillic filename ("й", `\u0439`) on the cp1252 console (UnicodeEncodeError). Unresolved names now printed ASCII-safe (`{n!a}`).

    - Bug fix 3 (git stage): the `.gitignore` initially inherited `data/processed/<file>` re-includes directly under `data/processed/*` exclusion — git cannot re-include children of an excluded directory, so the v2/v3 whitelist files silently didn't stage. Rewritten to `!data/processed/` + `data/processed/*.csv` + `!data/processed/cetane_pure_curated.csv`.
16. Restructure: at the root, non-core items were **moved** into a new `extras\` folder — `fuelblend_data` (the master repo with its .git), `fuelblend_university_bundle` + zip, `build_bundle.py`, `pre_abstract`, `reports`, `tasks`, the stale root `data\` (4 YSI sooting files) and stale root `src\` (3 .py copies — mtime-verified **older** than the repo versions before moving), `.pytest_cache`, `__pycache__`, and hidden `ajftd_session.json` (sensitive AJFTD cookies, found with a hidden-file listing; moved out of root).
17. The core was then built **into the root itself**: `python build_core_dir.py --src ...\extras\fuelblend_data --dst ...\UGRP_hydrocarbons` → **337 files**, `fixed-list: none`, unresolved = only the literal `.json` (a glob suffix in code, not a real file). The earlier stray copy at `C:\Users\kahaa\fuelblend_core` was deleted after the root build verified.
18. Verification at the root with `.venv\Scripts\python.exe`: `scripts/v3/reproduce_all.py --verify-only` → all 20 deliverables PASS; `scripts/v2/reproduce_all.py --verify-only` → all PASS.

### Stage 8 — Git setup and push (2026-11-07/09)

19. Facts at that point: root was **not** a repo; `extras\fuelblend_data` was (and is) the repo with history; its remote is `https://github.com/Kahaan83/UGRP_Hydrocarbons.git` with `main` at `b23bbe3` and tags `baseline-2026-11-03`, `v2-2026-11-05`, `v2.1-2026-11-05`, `v3-2026-11-06`.
20. Chosen path (user decision: "Option B", then "Replace main (force-push)"): `git init` at root → `git config core.longpaths true` → root `.gitignore` written (§10) → `git add -A` audited (268 files, 20.66 MB, staged >1 MB only `thermoml_clean_hc.csv` 10.7 MB — explicitly whitelisted in the repo's own ignore policy — and `leakage_crc_training.csv` 5.4 MB v2 deliverable; leak check for extras/ajftd/model_table/unified/pure_grid/session: 0 hits) → commit `759dcbf` ("Paper core v2+v3: scripts, fuelblend library, curated data, results, tests", 268 files, 210,894 insertions) → `git push --force origin master:main` (remote main moved `b23bbe3 → 759dcbf`, forced update) → local branch renamed `master` → `main`, upstream set to `origin/main`, tree clean.
21. The old main tip remains reachable via the GitHub tags and via the intact repo in `extras\fuelblend_data`.

---

## 5. Canonical data schema (the contract)

### 5.1 Properties and units (`data/unified/measurements.csv`)

| property | unit | notes |
|---|---|---|
| `flash_point_C` | °C | closed vs open cup ~5–10 K; method kept in ref/AJFTD method column |
| `viscosity_dyn_mPas` | mPa·s | all conversions in the loaders |
| `viscosity_kin_cSt` | cSt = mm²/s | = dyn / density[g/cm³] at same T |
| `density_kg_m3` | kg/m³ | liquid only |
| `vapor_pressure_kPa` | kPa | pure liquids; feeds reffit of Antoine |
| `freezing_point_C`, `pour_point_C`, `cloud_point_C` | °C | |
| `dcn`, `mon` | – | with `build_unified.py --octane` |
| AJFTD extras | – | `distillation_T{0,5,10,20,…,95,100}_C`, `distillation_{loss,recovery,residue}_pct`, `smoke_point_mm`, `thermal_stability_spent_fuel_mL` |

### 5.2 `measurements.csv` / `measurements_dedup.csv` columns

`source` (`leenhouts|nist_tde|thermoml|ajftd`), `property`, `value` (canonical unit), `T_K`, `P_kPa`, `kind` (`pure|binary|multi|real_fuel`), `n_comp`, `is_hc` (True iff all components are hydrocarbons), `system_id` (= md5[:12] of sorted InChIKey set; `fuel:<id>` for real fuels; used as GroupKFold key), `comp_inchikeys` (JSON sorted list), `x_mole`/`w_mass` (JSON aligned lists, both always provided), `fuel_id`, `ref`, `raw_value`, `raw_unit`.

`measurements_dedup.csv` adds: `T_key` (0.5 K bins), `P_key` (≤150 kPa or missing → 101; else round to 100 kPa), `x_key` (3 decimals), `value` (median), `vmin`, `vmax`, `q25`, `q75`, `n_obs`, `n_refs`, `n_sources`, `sources`, `n_outlier` (>25 % from median, needs ≥3 obs), `conflict` (IQR/range exceeds tolerance: flash 2 °C, pour/cloud/DCN/MON 3 °C, freezing 1 °C; 5 % relative otherwise).

### 5.3 `compounds.csv`

`inchikey` (PK), `smiles`, `inchi`, `formula`, `mw`, `is_hc`, `cls` (aromatic | alkene_alkyne | cycloalkane | iso_paraffin | n_paraffin | other), `n_carbon`, `antoine_A/B/C`, `antoine_Tlo/Thi`, `antoine_source` (`refit_measured|webbook`), `antoine_covers_window` (True iff valid over 313.15–373.15 K), `antoine_unit`, `lfl_vol_pct`, `lfl_T_C`.

**Antoine convention in this repo: `log10(P/kPa) = A − B/(T/K + C)`.** WebBook bar-form: `A_kPa = A_bar + 2`.

### 5.4 Derived tables feeding v2/v3

`model_table.csv` (9.7 MB): the core training table; contains only composition_level ∈ {L2, L3} rows (Rule R1 excludes real fuels, which are `real_fuel` kind without speciated lists). Splits live in `splits/*.json` (train/test grouped, locked sets, holdout, nested v2, s1s2s3, v3).

---

## 6. Physics models (what each baseline is)

- **Refutas viscosity** (basis of T9/v2 model): `VBN = 14.534·ln(ln(ν+0.8)) + 10.975`; `VBN_mix = Σ w_i·VBN_i`; `ν_mix = exp(exp((VBN_mix−10.975)/14.534)) − 0.8`. Mass fractions. Valid ν ≳ 0.2 cSt.
- **Walther (ASTM D341)** `thermo/walther_viscosity.py`: log-log η vs T fit per composition, with cryogenic extrapolation; used in R1 leaderboard against recalibrated Refutas, and in upgraded inverse designs.
- **PPBI pour point**: `PPBI_i = T_i^12.5` (T in K), volume-weighted, `T_mix = PPBI_mix^0.08`. Poor for waxy/paraffinic blends.
- **Coutinho-style wax model** `thermo/coutinho_wax.py`: non-ideal multi-component wax appearance/freezing (fitted γ/α/β, 3 parameters); Schroeder-van Laar ideal SLE in `baselines/cold_flow.py` (WAT = max over components).
- **Le Chatelier flash**: closed-cup flash point = T where `H_T = Σ x_i γ_i P_i^sat(T) / (P_total·LFL_i) = 1` under Raoult (γ=1); root-finded with brentq; Zabetakis LFL temperature correction `LFL(T) = LFL_25·(1 − 0.000721(T−25))`. Composition basis: mole fraction from mass via MW.
- **Vapor-pressure anchoring** `thermo/vapor_pressure_anchoring.py`: anchors Psat curves into the useful window (feeding `data/processed/v2/compounds_anchored_vp.csv`); used by R2/R3 flash work.
- **Heat of combustion**: v1 `baselines/heat_of_combustion.py` (elemental fallback, kept for comparison) and v2 `baselines/heat_of_combustion_v2.py` (Benson/Joback group additivity via RDKit SMILES, GCV→NCV correction).
- **Distillation**: `baselines/distillation.py` — TBP cuts from composition → API D86 conversion `T = a·TBP^b` (T10/T50/T90); DVPE activity model `models/dvpe_activity.py` for gasoline DVPE.
- **Octane/cetane**: `baselines/octane.py` (RON/MON linear-by-volume + interaction blend "BON", sensitivity) + `models/octane_ghosh_lumps.py` (Ghosh-style lump model, G2); `baselines/cetane.py` (linear DCN blend on mass or volume basis; carries the v1 basis fix CLM-05).
- **Composition conversions** (`common.py`): `x_i = (w_i/M_i)/Σ(w_j/M_j)`; `v_i = (w_i/ρ_i)/Σ(w_j/ρ_j)` (all ρ at the SAME temperature); all equation temperatures absolute (K).

---

## 7. Data-source inventory (what exists where)

| Source | What | Where | State |
|---|---|---|---|
| Leenhouts benchmark | curated pure+mixture flash/visc/DCN/MON | Zenodo 19914760 → `extras\fuelblend_data\data\raw\` | downloaded, verified; code gitlab.kuleuven.be/creas/vermeiregroup/mixture_evaluation |
| CheMixHub | ~500k mixture points incl. NIST logV | cloned at raw\git\mixture_evaluation sibling path | uses nist-logV only |
| Bilodeau NIST viscosities | raw binary diffs | Zenodo 8042966 → pretrained_models.zip → data\interim\bilodeau_nist_viscosity | downloaded |
| ThermoML archive | 5-journal measured density/viscosity/VP | ThermoML.v2020-09-30.tgz sha256 verified | extraction re-run post-fix |
| NIST WebBook | Antoine constants (bar) | scraped once | C10+ ranges don't reach 40–100 °C; refits cover some |
| AJFTD | real jet fuels | login-gated; 81 fuels scraped | thin coverage on flash/visc/density (1 fuel each) |
| NIST RP-1 | 11 orthogonal blends + 1 characterized sample | NISTIR 6646 + Fortin 2012 | blend recipes NOT in papers; one sample characterized |
| CRC 647 / AV-27-18 / AV-11-09 | aviation fuel property tables | crcao.org PDFs | ingested (4,559+869+3 rows) |
| Zabetakis 627 | LFL bulletin | PDF | 12 HC transcribed + fill_lfl → 695 compounds |
| Li et al. 2020 | 10,575 pure flash points | figshare 9275210 | loaded (14,696 rows) |
| YSI/sooting | YSI database volumes | data\raw\sooting (extras) | builder present |
| DIPPR | proprietary | — | not used |

`data/MISSING.md` (tracked) records what never arrived. `SOURCE_CHECK.md`, `DATA_LOG.md`, `SHA256SUMS.txt` (extras\fuelblend_data root) carry the source-verification trail.

---

## 8. Operation guide (how to run what)

```powershell
cd C:\Users\kahaa\OneDrive\Desktop\UGRP_hydrocarbons
# venv already at .\.venv  (activate: .\.venv\Scripts\Activate.ps1)

# reproduce v2 results (Phases A+B direct, then run_phase_c.py, run_phase_d.py, run_phase_e.py)
.\.venv\Scripts\python.exe scripts\v2\reproduce_all.py            # add --verify-only to just check deliverables

# reproduce v3 results (13 steps + pytest + SHA256SUMS manifest)
.\.venv\Scripts\python.exe scripts\v3\reproduce_all.py            # --verify-only | --skip-tests

# run the test suite
.\.venv\Scripts\python.exe -m pytest tests -q

# public prediction API
from fuelblend.api.predict_blend import predict   # src on sys.path; loads runs\<pure|mixture>_models\*.joblib, lazily retrains mix residual models if absent

# rebuild the core folder elsewhere (the build script with all v3 edits):
python "C:\Users\kahaa\Downloads\files (3)\build_core_dir.py" --src C:\Users\kahaa\OneDrive\Desktop\UGRP_hydrocarbons\extras\fuelblend_data --dst <target>
# NOTE: the script will not overwrite existing files; wipe target core dirs to rebuild fresh.
```

Data-collection-era commands (live only in `extras\fuelblend_data`): `python fetch.py all`, `inspect_downloads.py`, `build_viscosity.py`, `thermoml_extract.py`, `webbook_antoine.py`, `refit_antoine.py`, `build_unified.py`, `ingest_ajftd.py`, `audit.py --bundle`, `fill_lfl.py`, `ingest_crc.py`, `ingest_rp1.py`, `ingest_dcn.py`, `build_octane_hc.py`, `build_gcv.py`, etc.

---

## 9. Test coverage in the core (57 files)

- `tests/test_v2_*.py` (39): one suite per v2 task (A1-A12, B, C1-C8, D1-D7, E, F-freeze, inverse-design-upgraded, new sources, thermoml tidy, vapor-pressure anchoring, walther, coutinho, ablation, licenses).
- `tests/test_v3_*.py` (11): phase0 regression, phase1 expansion, phase2 repairs, phase3 designs, R1-R7 (7 suites).
- Plus: `test_bilinear_pairwise.py`, `test_predict_blend.py`, `test_model_invariances.py`, `test_harness.py`, `test_validate_schema.py` (uses `fuelblend.data.validate.validate_measurement_df`), `test_result_regression.py` (pins v2 headline numbers within ±0.2 °C / ±0.005 ln), `test_migration_integrity.py` (stdlib-only imports of v1-era tests like test_offline/test_thermoml/test_unified are NOT in the core — those live with the data-collection scripts in the inner repo's `tests\`).
- `test_migration_integrity.py` requires `data/processed/measurements_migrated_core.csv` == core rows of `data/processed/model_table.csv` (sha256 of value-column rounding must match).

---

## 10. Git mechanics in force

- Root `.gitignore` (tracked): excludes `extras/`, `.venv/`, `debug/`, `__pycache__/`, `*.pyc`, `.pytest_cache/`, `runs/`, `*.pt`, `ajftd_session.json`, `ajftd_data.csv/json`; `data/*` with re-includes (`!data/MISSING.md`, `!data/antoine_sources.csv`, `!data/real_fuels/`, `!data/processed/` then `data/processed/*.csv` + `!data/processed/cetane_pure_curated.csv`); `archive/*` + `!archive/baseline_manifest.json`, `!archive/v2_manifest.json`; `results/t0_repro/reference_backup/`, `results/t0_repro/rerun_outputs/`.
- The git-can-only-reinclude-inside-reincluded-directory subtlety is documented inline in the `.gitignore`.
- 268 tracked files ≈ 20.7 MB; >1 MB tracked: `data/processed/v2/thermoml_clean_hc.csv` (10.7 MB), `results/t0/leakage_crc_training.csv` (5.4 MB) — both deliberate.
- Remote: `origin = https://github.com/Kahaan83/UGRP_Hydrocarbons.git`; local `main` ↔ `origin/main` in sync; `core.longpaths true` set.
- The inner repo `extras\fuelblend_data\.git` still holds the complete layered history (its remote points at the same URL; its branch main was left at b23bbe3 equivalent state).

---

## 11. Unit traps and bugs — fixed list (do not reintroduce)

1. NIST logV column: value is ln(η/Pa·s) despite the "Log viscosity / cP" label. Leenhouts viscosity: ln(mPa·s). Difference: 6.908. (Exception: thermoml_clean_hc.csv was re-derived by v2's tidy step, uses its own documented units.)
2. WebBook Antoine in bar → kPa: `A_kPa = A_bar + 2`. ThermoML Pa·s→mPa·s ×1e3; m²/s→cSt ×1e6. Flash unit string "Celcius" (sic).
3. ThermoML `Constraint` nesting bug (61 % temperature loss) — fixed in the rewritten `thermoml_extract.py`.
4. Antoine refit negative-B start crash → clipped starts, correlation guard, per-compound `failed:` status.
5. Zenodo macOS junk (`__MACOSX/._*`) → skipped everywhere by name parts.
6. compounds.csv had an unquoted comma-containing name (`1,2-...`).
7. WebBook `page_name` contained the site banner → fixed to real title.
8. AJFTD scraper named every fuel "Alternative Jet Fuels" → fixed to tombstone label; no resume → fixed; pagination not crossed-checked → added.
9. `fetch.py`: fake DTIC HTML-as-PDF validation, `.part` resume, 429 back-off, per-file deadline, Ctrl-C summary.
10. pandas 3: `.loc[k, [cols]] = [mixed types]` fails → assign per column.
11. `build_core_dir.py` (this session): scanner scope restricted to `src/scripts/tests` (extras-contained repos raked junk); unicode-print crash fixed with ASCII-escape presentation; gitignore re-include rule corrected as described in §4 step 15.
12. Synthetic-data circularity: v1's 0.4 % came from training on rule-generated labels; the entire real-data pipeline exists as a counter-move.

---

## 12. Known limitations / unresolved (state only, as of this handoff)

- `requirements.txt` does not pin `scikit-learn`, `torch`, `joblib`, `pytest`, though core code imports them; venv has them installed.
- F1–F6 (from §4 Stage 5) were still open when v2 was tagged.
- The two fake DTIC PDFs remain 0.0 MB placeholders.
- AJFTD pagination fit-check ("Showing x to y of N entries") never re-verified against N > 81 after the scraper fix.
- The old (pre-fix) ThermoML contamination story is historical; the current `thermoml_hc.csv` at core root is post-fix output.
- Antoine for heavy odd alkanes (C13/15/17/19) absent (no vapour-pressure data found); C11 only 393–459 K; C12/C20/C24 status `check` (rms 0.06–0.07 log10).
- Core-root `results/` also contains v1-era dirs (`t0`, `t6`, `t8`, `eval`, `hindcast`, `pairwise`) that were tracked by the inner repo's policy and kept.
- `archive/baseline_2026-11-03_code.zip` preserves the code state at the 2026-11-03 baseline tag; `baseline_manifest.json` / `v2_manifest.json` pin checksums.
- Open-data composition absolute limit recorded: essentially no ≥3-component hydrocarbon lab data exists anywhere (3 ternary flash rows) — real-fuel sources (AJFTD/RP-1/CRC) carry all multi-component validation.

---

## 13. Working conventions used throughout

- Hydrocarbon-only (`is_hc`) is the default filter; `--all-molecules` used deliberately.
- Training/evaluation always reads `measurements_dedup.csv` (or its v2 descendant `model_table.csv`), never raw rows; `conflict`, `n_outlier`, `n_refs` travel as quality weights.
- All unit conversions happen in the loaders (`build_unified.py` and friends); downstream code consumes canonical units only.
- Every new source got a loader to the canonical schema + a test with a fixture mirroring the **real** file structure + an audit check.
- Never commit `ajftd_session.json`, raw AJFTD pages, or licensed data; `audit_bundle.zip` contains no raw datasets by design.
- Reports distinguish what ran on real data vs fixtures.
- Grouped validation (`system_id` GroupKFold; leave-fuel-out on CRC 57; locked holdout splits in `splits/`) is how v2/v3 numbers were produced.
- Splits/locked sets and heavy derived CSVs are treated as data; code changes that alter them are logged in the claims register and manifests.

---

## 14. Pointers for the next agent (where everything lives)

- SUPERSEDED for numbers by docs/FULL_REPORT_v4.md (read it first). v2/v3 reports contain retracted claims: `docs\FULL_REPORT_v3.md`, `docs\FULL_REPORT_v2.md`, `docs\claims_register.md`, `docs\model_decision.md`, `docs\scope_decisions.md`, `docs\data_sources.csv`.
- Previous-era handoff: `AGENTS.md` at root (v1/data-collection era; its §12/§12.1 contain plans that are NOT part of current state).
- Inner-repo-era peripherals: `extras\fuelblend_data\{FILE_MAP.md, PROJECT_HANDOFF.md, T8_T9_INTERFACE_SPEC.md, FULL_REPORT.md, PROGRESS_REPORT.md, AUDIT.md, DATA_LOG.md, SOURCE_CHECK.md, CHANGES.md, SHA256SUMS.txt, BASELINE_VISCOSITY_REPORT.md, DATA_ACQUISITION_TASKS.md}`.
- The build script with all edits: `C:\Users\kahaa\Downloads\files (3)\build_core_dir.py`.
- Scratch notes (reports\tasks in `extras\`): `ANALYSIS_PLAN.md`, `WP0.md`, `plan.md`, `todo.md`.

---

## 15. v4 status (appended in v4 Task 19, branch `v4-integrity`)

* The v3 real-fuel corpus and design optimizer were formula-derived; they are quarantined (`data/quarantine/`, `results/v4/P1/`) and the dependent v3 claims are retracted (CLM-20..24) or revised (CLM-29) (`docs/claims_register.md`). v4 tables live in `data/real_fuels/v4/` (61 jet fuels; 57 CRC 647 with reconstructed compositions; 210 gasoline octane rows; diesel: none).
* Read `docs/FULL_REPORT_v4.md` first (every number cited to a `results/v4` file), then `docs/claims_register.md` (v4-updated rows CLM-09..28, new CLM-30..37), `docs/model_decision.md` sections 5-6.
* Reproduce: `python scripts/v4/reproduce_all.py` (use `.venv\Scripts\python.exe`; ~12 min plus pytest ~2.5 min). `--verify-only` checks row counts and rewrites `SHA256SUMS_v4.txt`; `--include-encoder` retrains the Task 16 encoder (slow). Steps that need `extras/` are skipped and reported when it is missing. Two consecutive runs gave an identical manifest.
* Tests: `pytest tests -q` -> 308 passed, 42 failed. All 42 failures are pre-existing "deliverable file missing" asserts for v2/v3 files (the baseline was 45 failed; three were removed by the quarantine move); list in `results/v4/baseline_pytest.txt`.
* Open items: ThermoML default-mode `phase_of()` reads only `ePhase` (CLM-37, not fixed, existing `thermoml_hc.csv` may be contaminated); mesitylene / n-hexylbenzene Psat inputs (CLM-36); density bias cause; no GPU run of the set encoder (`jobs/v4_gpu/` untested); diesel data; D1655 smoke limits and BIS spec values UNVERIFIED; DVPE units; IS 2796 DVPE minimum not citable.
* Uncommitted user work (not part of v4 tasks): `results/v3/R1/viscosity_leaderboard.csv`, `splits/locked.json`, `src/fuelblend/data/model_table.py`, `scripts/v4/compare_antoine_cox_gate.py`, `results/v4/antoine_cox_gate/`.
