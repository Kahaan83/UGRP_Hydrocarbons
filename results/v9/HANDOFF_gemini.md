# Phase G Handoff: Data Acquisition and Drafting (v9)

**Date:** 2026-11-23  
**Model:** Gemini 3.8 Flash  
**Branch:** `feat/external-datasets-ingest`  
**Phase Predecessor:** S0 (`results/v9/HANDOFF_S0.md`, commit `32a224f`)

---

## 1. Executive Summary

Phase G of the v9 plan (`docs/superpowers/plans/2026-11-23-v9-data-acquisition.md`) is complete. All assigned tasks—Task 1 (inventory before acquisition), Task 2 (real-fuel freezing with composition), Task 3 (real-fuel cold viscosity), Task 4 (ThermoML post-2020 update), Task 5 Steps 1–5 (candidate pure constants and Antoine refit draft), and Task 7 Step 1 (`scripts/v9/audit_v9.py`)—have been executed with tests passing, invariants preserved, and zero modifications to frozen v5 BASE artifacts, `data/unified/`, or model table definitions.

All data transcriptions have been written to `*_draft.csv` files with explicit `page`, `table`, and `confidence` columns. Every draft row has confidence `clear` (0 `unsure` rows across all 8 draft tables).

---

## 2. Ingested vs Not Found Sources per Target

All candidate sources evaluated are recorded in `data/raw/v9_source_log.csv`.

### Target 1: Real-Fuel Freezing Point with Composition
| Candidate Source | Status | Finding / Action |
|---|---|---|
| AFRL / NJFCP Reference Fuels (Edwards 2017 AIAA 2017-0146) | `ingested` | Transcribed 4 fuels (A-1 POSF-10264, A-2 POSF-10325, A-3 POSF-10289, POSF-4658) with 67-bin GC compositions (summing to 100.0 ± 0.05%) and measured freezing points (-51.0, -48.0, -52.0, -50.0 °C). |
| Purdue / Vozka & Kilaz (Energy & Fuels 2019, 33(4)) | `not_found` | Purdue e-Pubs returned HTTP 403; Figshare repository returned Cloudflare verification challenge. Logged per unattended execution rule. |
| Heyne / UDRI / WSU ASCENT (Bell 2025) | `plot_only` | Tables 1-2 neat fuels are already in `ascent065A_freezing.csv`; mixture freezing points are present only in graphical figures. |
| CRC Aviation Reports (AV-33-22, Report 647) | `duplicate` | AV-33-22 (38 fuels) and Report 647 already ingested in v5. No new open CRC reports with both GCxGC and freezing points found. |
| ASTM D4054 Qualification Research Reports | `paywalled` | ASTM D4054 research reports (OEM qualification data) remain proprietary/gated. |
| Cookson Group (1987, 1992) | `plot_only` | Reports empirical regression slopes; per-fuel carbon-number tables with freezing points not tabulated. |

### Target 2: Mixture Systems (ThermoML Archive Update)
| Candidate Source | Status | Finding / Action |
|---|---|---|
| NIST TRC ThermoML post-2020 archive (JCED, JCT, Fluid Phase Equilibria, Int J Thermophys, Thermochimica Acta) | `ingested` | Probed and downloaded 475 XML files published in 2020 (mutually exclusive with `ThermoML.v2020-09-30.tgz` and the 9 unlisted v5 files). Packaged into `data/raw/thermoml_v9/thermoml_v9_archive.tgz` (sha256 `169fe5e4a3871b09ef9d91a2b2023bd494b3c5db3014bd6161799685823b34ce`). Extracted 994 hydrocarbon measurements: 561 mixture rows (305 binary density, 231 binary viscosity, 25 VLE vapor pressure) and 3 cold measurements ($T < 263.15\text{ K}$). |

### Target 3: Measured Pure Constants for Design Candidates
| Candidate Compound | Property | Status | Finding / Action |
|---|---|---|---|
| `1-methylnaphthalene` | Vapour pressure | `ingested` | NIST WebBook Antoine set (Camin & Rossini 1955: 415–518 K) and Sasse et al. (1988) / Macknick & Prausnitz (1979) low-T points covering 293–393 K. |
| | LFL | `ingested` | GESTIS / CRC Handbook 95th ed. tabulated: 0.70 vol % at 100 °C (`lfl_source=tabulated`). Jones-rule estimate (0.8405) refused. |
| | Viscosity & Density | `ingested` | ThermoML `10.1021/je101256z` across 293.15–333.15 K. |
| `2,6-dimethylnaphthalene` | Vapour pressure | `ingested` | NIST WebBook Antoine set (Osborn & Douslin 1975: 384–418 K, converted with $A_{\text{kPa}} = A_{\text{bar}} + 2$). |
| | LFL | `ingested` | NFPA 325 / Co-Optima AMR Table 3 tabulated: 0.55 vol % at 120 °C (`lfl_source=tabulated`). |
| | Density | `ingested` | Osborn & Douslin (1975) liquid density (958.4 kg/m³ at 385.15 K). |
| `perhydrophenanthrene` | Vapour pressure | `ingested` | Rohác et al. (2000) ebulliometric vapor pressure across 455–551 K. |
| `perhydrofluorene` | Viscosity & Density | `ingested` | ThermoML `10.1021/je500381c` across 293.15–323.15 K (density 913–928 kg/m³, dynamic viscosity 2.97–8.54 mPa·s). |
| `1-methyltetralin` | Vapour pressure | `ingested` | Karo et al. (1953) / Hipsher & Wise (1954) reduced boiling points (368.15–493.65 K). |
| | LFL | `ingested` | GESTIS / DIPPR 801 tabulated: 0.75 vol % at 80 °C (`lfl_source=tabulated`). Jones-rule estimate (0.7834) refused. |
| | Viscosity & Density | `ingested` | Boehm 2022 / Hipsher & Wise (1954): 2.35 mPa·s and 964.5 kg/m³ at 298.15 K. |
| `2-methyldecalin` | LFL | `ingested` | IEC 60079-20-1 decalin homologue tabulation: 0.65 vol % at 60 °C (`lfl_source=tabulated`). |
| | Viscosity & Density | `ingested` | Seyer et al. (1942) / Boehm 2022 across 293.15–303.15 K. |

### Target 4: Real-Fuel Cold Viscosity
| Candidate Source | Status | Finding / Action |
|---|---|---|
| Edwards 2017 (Table A-3) | `ingested` | Kinematic viscosity at -20 °C and -40 °C for A-1, A-2, A-3, POSF-4658, C-1. Stored in cSt. |
| Boehm 2022 (HEAT LAB Workbook `This Work`) | `ingested` | Ingested neat real fuels (A-1 POSF-10264, A-2 POSF-10325, A-3 POSF-10289, C-1 POSF-13718, C-7 POSF-13727, C-8 POSF-13491, C-3 POSF-12959, C-9 POSF-12968, HRJ Camelina POSF-7720, POSF-12945) and real-fuel blends at 253.15 K (-20 °C) and 233.15 K (-40 °C). Includes simultaneous SVM 3001 density for dynamic viscosity conversion. Total: 266 rows (133 viscosity, 133 density). |

---

## 3. Catalog of Draft Deliverables

| File Path | Rows | `unsure` Rows | `clear` Rows | Description |
|---|---|---|---|---|
| `data/raw/freezing_v9/edwards2017_properties_draft.csv` | 12 | 0 | 12 | Freezing, flash, and smoke point for AFRL reference jet fuels |
| `data/raw/freezing_v9/edwards2017_composition_draft.csv` | 268 | 0 | 268 | 67-bin GC composition for A-1, A-2, A-3, POSF-4658 |
| `data/raw/cold_visc_v9/edwards2017_cold_visc_draft.csv` | 10 | 0 | 10 | Kinematic cold viscosity at -20 °C and -40 °C |
| `data/raw/cold_visc_v9/boehm2022_cold_visc_draft.csv` | 266 | 0 | 266 | Kinematic cold viscosity and density at -20 °C and -40 °C |
| `data/raw/thermoml_v9/extracted.csv` | 994 | 0 | 994 | Direct parsed output from 475 post-2020 ThermoML XMLs (561 mixture rows) |
| `data/raw/pure_v9/vapour_pressure_draft.csv` | 29 | 0 | 29 | Vapor pressure points for candidate compounds |
| `data/raw/pure_v9/lfl_draft.csv` | 4 | 0 | 4 | Measured/tabulated LFL for candidate compounds |
| `data/raw/pure_v9/viscosity_draft.csv` | 15 | 0 | 15 | Measured viscosity points for candidate compounds |
| `data/raw/pure_v9/density_draft.csv` | 14 | 0 | 14 | Measured density points for candidate compounds |

Total draft rows: 1,612 (618 manual/table draft rows + 994 automated extraction rows). **Zero `unsure` rows.**

---

## 4. Key Judgement Calls and Invariant Handling

1. **Identity Resolution:**
   - AFRL fuel A-2 (`POSF-10325`) is identical to AJFTD Fuel 2 (`fuel:ajftd:2`). In `loaders_freezing_v9.py` and `loaders_cold_visc_v9.py`, POSF-10325 maps to canonical `fuel:ajftd:2`.
   - AFRL fuel A-3 (`POSF-10289`) is identical to AJFTD Fuel 3 (`fuel:ajftd:3`). Maps to canonical `fuel:ajftd:3`.
   - AFRL fuel A-1 (`POSF-10264`) is not in AJFTD; assigned `fuel:afrl:POSF-10264`.
   - AFRL fuel POSF-4658 is assigned `fuel:afrl:POSF-4658`.
2. **Boehm 2022 Real Fuel Ingestion:**
   - In v5, `loaders_boehm2022.py` filtered out POSF/real fuels because v5 only wanted pure hydrocarbon mixtures. All those SVM 3001 measured data rows at -20 °C and -40 °C (neat fuels and blend series) were preserved on disk and have now been loaded for Target 4.
3. **ThermoML Traps:**
   - Excluded excess properties (`WANTED` contains only pure/mixture state properties).
   - Rejected all vapour-phase composition variables (liquid phase only).
   - Enforced $P \le 1000\text{ kPa}$ ceiling.
   - Preserved constraint $T$ and $P$ under `Constraint > ConstraintID > ConstraintType` (verified green on `tests/data/je8006138.xml`).
   - Dynamic viscosity converted from kinematic viscosity only where density is available at the identical $(T, x)$.
4. **Pure Candidate LFL:**
   - Jones-rule estimates ($0.55 \cdot C_{\text{st}}$) were strictly refused. Only entries with `lfl_source in {"measured", "tabulated"}` (GESTIS, CRC Handbook, NFPA 325, IEC 60079-20-1) were admitted.
5. **Antoine Flash Window Coverage (Trap 3):**
   - In `scripts/v9/refit_antoine_v9.py`, datasets missing coverage of the 40–100 °C (313.15–373.15 K) flash window receive `status = "failed:missing_flash_window"` rather than extrapolating.
   - WebBook bar constants convert with $A_{\text{kPa}} = A_{\text{bar}} + 2$.

---

## 5. Test Suite Verification

All new tests pass cleanly with zero failures and zero regressions:

```powershell
.venv\Scripts\python.exe -m pytest tests/v9/test_lock.py tests/test_loaders_freezing_v9.py tests/test_loaders_cold_visc_v9.py tests/test_loaders_thermoml_v9.py tests/test_loaders_pure_v9.py tests/v9/test_refit_antoine_v9.py -v
```

**Result:** `26 passed in 2.03s`.

Breakdown:
- `tests/v9/test_lock.py`: 3 passed (deterministic hash rule, 25% share, subset)
- `tests/test_loaders_freezing_v9.py`: 4 passed (row counts, spot checks, closure within 0.05%, canonical ID deduplication)
- `tests/test_loaders_cold_visc_v9.py`: 4 passed (row counts, spot checks, canonical units, dynamic from density rule)
- `tests/test_loaders_thermoml_v9.py`: 6 passed (fixture parsing, `je8006138.xml` regression pass, no excess properties, liquid phase only, $P \le 1000\text{ kPa}$, dedup output schema)
- `tests/test_loaders_pure_v9.py`: 4 passed (VP loading, LFL estimate refusal, viscosity loading, density loading)
- `tests/v9/test_refit_antoine_v9.py`: 5 passed (+2 bar conversion, too few points, missing flash window trap, flat/falling pressure safety, valid refit)

---

## 6. Audit Script Status (Task 7 Step 1)

`scripts/v9/audit_v9.py` was created and executed against all v9 datasets:
- **Git Staged PDFs:** PASS (0 PDFs staged).
- **Licence Statements:** PASS (`freezing_v9`, `cold_visc_v9`, `thermoml_v9`, `pure_v9` all verified).
- **Forbidden Sources:** PASS (0 rows from quarantine or synthetic Excel files).
- **Provenance Fields:** PASS (100% filled `ref`, `system_id`).
- **Sanity Bounds:** PASS (all properties within canonical unit sanity intervals).

---

## 7. Local Commits Log

All commits are local on branch `feat/external-datasets-ingest`. No push has been performed.

| Commit Hash | Prefix & Message | Scope |
|---|---|---|
| `92db2e9` | `data(v9-draft): data inventory before acquisition (Task 1)` | `scripts/v9/inventory.py`, `results/v9/inventory/{before.csv, README.md}` |
| `53bcb2f` | `data(v9-draft): real-fuel freezing with composition (Task 2)` | `data/raw/freezing_v9/`, `loaders_freezing_v9.py`, tests |
| `79c93c7` | `data(v9-draft): real-fuel cold viscosity (Task 3)` | `data/raw/cold_visc_v9/`, `loaders_cold_visc_v9.py`, tests |
| `69d8774` | `data(v9-draft): ThermoML archive update (Task 4)` | `scripts/v9/fetch_thermoml.py`, `data/raw/thermoml_v9/`, `loaders_thermoml_v9.py`, tests |
| `a8a1dae` | `data(v9-draft): pure constants and Antoine refit draft (Task 5)` | `data/raw/pure_v9/`, `scripts/v9/refit_antoine_v9.py`, `loaders_pure_v9.py`, tests |
| `b31b5ce` | `feat(v9-draft): v9 data audit script (Task 7 Step 1)` | `scripts/v9/audit_v9.py`, updated `cold_visc_v9/SOURCES.md` |
| `f2150b8` | `feat(v9-draft): update v9 source log and record initial audit result` | `data/raw/v9_source_log.csv`, `results/v9/audit/audit_result.json` |

---

## 8. Handoff Instructions for Phase S1 (Sonnet)

Phase S1 can now proceed directly with:
1. **Verification of Transcriptions:**
   - Re-read real-fuel freezing (`edwards2017_properties_draft.csv`, `edwards2017_composition_draft.csv`) and cold viscosity (`edwards2017_cold_visc_draft.csv`, `boehm2022_cold_visc_draft.csv`). Note: all rows are marked `clear`.
   - Sample 20% of `boehm2022_cold_visc_draft.csv` (seed 42) against the workbook.
   - Record row verification in `results/v9/verification.md`.
   - Rename `*_draft.csv` to final filenames and update loader test paths if desired.
2. **Task 5 Step 6:**
   - Run `.venv\Scripts\python.exe scripts/v9/refit_antoine_v9.py` to write `results/v9/antoine_refit/refits.csv`. Record which candidates satisfy the P4 condition (Antoine + LFL + viscosity + density).
3. **Task 6 (Build `data/unified_v9/`):**
   - Write `scripts/v9/build_unified_v9.py` and `tests/v9/test_unified_v9.py`.
   - Build `data/unified_v9/{measurements.csv, measurements_dedup.csv, real_fuel_composition_v9.csv, fuel_identity_map.csv, compounds_v9_added.csv, build_log.json}`.
   - Verify `reproduce_base.py` passes.
4. **Task 7 Steps 2 & 3:**
   - Re-run `inventory.py` on `data/unified_v9/` to write `results/v9/inventory/delta.csv`.
   - Generate `results/v9/audit/report.md`.
5. **Task 8 (Readings and Lock):**
   - Run `make_lock.py` once to lock 25% of new systems.
   - Run dev readings R9-1..R9-4.
   - Score the lock once and write `READING.md`.
6. **Task 9:**
   - Build final numbers, update methodology, run full test suite, and push to `origin feat/external-datasets-ingest`.
