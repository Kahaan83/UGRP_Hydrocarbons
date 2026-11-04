# fuelblend_data - data collection pipeline

Run on a machine with open internet (university network is fine). Python 3.10+.

```bash
pip install -r requirements.txt
python fetch.py all                 # git repos, Zenodo records, PDFs/reports, ThermoML archive
python inspect_downloads.py        # (skips macOS junk; summarises splits)
# unzip Zenodo files, write data/interim/INVENTORY.md
python build_viscosity.py           # hydrocarbon-only NIST/TDE viscosity (pure + binary mixtures)
python thermoml_extract.py          # hydrocarbon viscosity/density/vapour pressure from ThermoML
python webbook_antoine.py           # Antoine constants (bar->kPa, range check, Psat at 40 C)
# --- after the downloads: build, unify, refit ---
python thermoml_extract.py          # -> data/processed/thermoml_hc.csv   (a few minutes; prints progress)
python webbook_antoine.py           # WebBook Antoine sets (high-T only for C10+ alkanes!)
python ingest_ajftd.py --describe   # AJFTD: writes data/interim/ajftd_description.txt (send it), then: python ingest_ajftd.py
python build_unified.py             # ONE canonical table: data/unified/{compounds,measurements,measurements_dedup}.csv + UNIFIED_REPORT.md
python refit_antoine.py             # Antoine refit from measured vapour pressures (needs thermoml_hc.csv)
python build_unified.py             # rerun so compounds.csv carries the refit constants
python audit.py --bundle             # checks everything, writes AUDIT.md + audit_bundle.zip to send back
python tests/test_offline.py        # sanity tests (no network)
python tests/test_fetch_server.py   # downloader behaviour vs a local mock server
```

Useful flags: `python fetch.py zenodo --list` (see files/sizes, download nothing),
`--include REGEX` (only matching Zenodo files), `--max-gb N` (size guard), `--skip zenodo`.
Ctrl-C skips only the current group; the failure/manual summary is always printed.
If Zenodo still fails from your network, open https://zenodo.org/records/19914760 and
https://zenodo.org/records/8042966 in a browser and put the files in data/raw/zenodo/<name>/.

`fetch.py` is idempotent, resumes partial downloads, backs off on HTTP 429, times out stalled files, writes `data/raw/MANIFEST.json` (url, size, sha256, time),
and prints at the end what failed and what must be done by hand.

## Status (what was actually tested)
| Piece | Status |
|---|---|
| build_viscosity.py on CheMixHub NIST data | RAN on real data (outputs in data/processed/) |
| git clone (CheMixHub) | RAN |
| Zenodo / ThermoML / PDF downloads | NOT run against the real hosts (blocked where this was built). Downloader tested against a local mock server: 429 backoff, stalled-file timeout, resume, size guard, Ctrl-C. ThermoML direct URL is a best guess - manual fallback printed. |
| thermoml_extract.py | Tested on a synthetic ThermoML file only. Run `--debug` on the first real archive and check counts. |
| webbook_antoine.py | Parser + unit conversion tested on synthetic HTML (benzene check: 24.3 kPa at 40 C). Not run live. Verify `page_name` matches `name` and CAS numbers in compounds.csv. |

## Things that bit us - keep these in your cleaning code
* NIST logV column `value` is labelled "Log viscosity / cP" but is **ln(eta / Pa.s)** (hexane 298 K -> 0.296 mPa.s).
* NIST data is **dynamic** viscosity. Refutas needs **kinematic**: nu[cSt] = eta[mPa.s] / rho[g/cm3] at the same T.
* WebBook Antoine is log10(P/**bar**) with T in **K**; the methods doc wants kPa (x100).
* WebBook Antoine validity ranges are often narrow or at high T. `covers_window=False` => refit.
* Benchmarks store **mole** fractions; your viscosity blending uses **mass** fractions (build_viscosity.py converts).
* Flash point: record closed vs open cup (differ by ~5-10 K).

## Not automatable (see fetch.py MANUAL list)
AJFTD (registration + admin approval), RP-1 11-blend data (Energy & Fuels SI), pour/cloud-point
literature tables, Li et al. flash-point SI, DIPPR (proprietary), LFL tables (Zabetakis PDF is
downloaded; values must be transcribed).

## Corrections / known issues (after first real run)
* ThermoML filename is `ThermoML.v2020-09-30.tgz` (dash, not dot); fetch.py now checks size and sha256.
* CORRECTION of an earlier note: Zenodo record 8042966 (`pretrained_models.zip`) DOES contain the raw NIST viscosity
  table (`nist_dippr_source/NIST_Visc_Data.csv`, 239,200 rows, Visc in Pa.s). Keep it. (It is also what CheMixHub processes.)
* The Zenodo zips were made on a Mac: they contain `__MACOSX/._*.csv` binary shadow files that look like tables but
  are not. inspect_downloads.py now skips them (that was the source of the hundreds of parse errors).
* apps.dtic.mil blocks scripts (HTTP 403). fetch.py retries with browser-like headers; if that still fails
  download the two PDFs in a browser into data/raw/pdfs/ using the filenames printed in the FAILED list.
* (audit 2) fetch.py used to save an HTML consent page from apps.dtic.mil as a "successful" 1.4 KB PDF. It now
  checks the %PDF header, deletes bad files and tells you to download them in a browser. The old fetch_failures.txt
  is rewritten (or deleted) on every run.
* `~$Something.xlsx` files are Excel lock files (a workbook was open in Excel while unzipping) - ignored.
* `DatasetsFuelProperties.xlsx` holds 50,000 SYNTHETIC rows ('Synthetic value'): never use it as ground truth.

## Unified schema (data/unified/)
* `measurements_dedup.csv`: one row per (property, system, composition, T); `value` = median of independent observations;
  `conflict` / `n_outlier` / `n_refs` say how much to trust it. `system_id` = hash of the component SET -> GroupKFold key.
* Units: flash_point_C, viscosity_dyn_mPas, viscosity_kin_cSt, density_kg_m3, vapor_pressure_kPa, freezing/pour/cloud_point_C.
* Antoine in compounds.csv is always log10(P/kPa) = A - B/(T/K + C) (WebBook bar sets are converted: A+2).
* Raw NIST/TDE viscosity: ~1/3 of rows are exact repeats and a few references are clear outliers (benzene 25 C: ~70 refs give
  ~0.60 mPa.s, one gives 0.395). Always train on measurements_dedup.csv, never on raw rows.

---

## Progress (real-run results — 2026-11-02 session)

Pipeline now runs end-to-end on this machine (Windows PowerShell, venv `.venv`). Key files updated: `build_unified.py` (CRC loader), `fetch.py`, `audit.py`. Reports (`reports/T*`), scripts (`scripts/t*`), and tests (`tests/test_t1_leakage_audit.py`, `tests/test_predict_blend.py`) all green.

| Task | Deliverable | Status | Key number / note |
|---|---|---|---|
| T0 Coverage audit | `reports/T0_COVERAGE.md` | Done | 710 pure-flash rows; mixture-flash critical gap (123 binary / 3 multi before T7 fix) |
| T1 Family leakage audit | `src/fuelblend/features/family.py`, `reports/T1_LEAKAGE.md` | Done | Family-holdout: MAE 0.75 / RMSE 1.73 °C (honest claim; 0.56 is molecule-grouped optimistic) |
| T2 Wrapper/API | `src/fuelblend/api/predict_blend.py` (subagent: `google/antigravity-claude-sonnet-4-6`), `tests/test_predict_blend.py` (5/5 PASS) | Done | `predict(composition, T)` returns value/std/domain-warnings; uses lazy-train for mixture residuals |
| T3 Complexity errors | `reports/T3_COMPLEXITY.md` | Done | Binary viscosity RMSE 0.27 ln (median +3.8%); outlier system `a175950749e8` flagged |
| T4 Real-fuel scoring | `reports/T4_REAL_FUELS.md` | Done | CRC 647 (57 fuels) — flash MAE 53.6 °C bias (surrogate granularity failure, not math error); density MAPE 6.9%, viscosity MAPE 76% |
| T5 Mixture flash | `reports/T5_FLASH_MIXTURE.md` | Done | Bias calibrator +1.69 °C; S2 MAE 0.835 °C; 81.8% 90%-interval; 0% over-prediction >5 °C |
| T6 Extensions (GCV, MON/DCN, freezing/pour/cloud) | `reports/T6_EXTENSIONS.md` | Done | GCV excluded (0 pure); freezing/pour/cloud excluded (L0 only); MON/DCN null baselines: MAE 16.4 / 17.3 °C |
| T7 Ingest missing flash mixtures | `data/processed/t7_leenhouts_flash_mix.csv` (1006 rows; 681 binary / 325 multi) | Done | 813 new mixture-flash rows; loader `build_unified.py` updated |
| T8 DeepSets comparison | — | Open | Depends on T7-enlarged flash-mix data |
| T9 Inverse design | — | Open | Pending |

Repair notes added to this session:
- `build_unified.py` line 350: CRC `density_kg_m3` now maps `T_C=60` → `T_K=288.71`, avoiding `NaN`; added `shortcut:` upgrade path.
- `tests/test_unified.py`: schema verification runs; note some fixtures may differ.
- `t4_real_fuels.py`: negative flash result documented honestly; density loader `T_K` fixed.
- `t5_flash_mixture.py`: split-conformal quantiles calibrated on S1 mixture residuals only (n=34); S2 coverage 81.8%.
- `t6_extensions.py`: MON/DCN null baselines implemented; GCV/freezing/pour/cloud excluded with reasons.
- `t7_ingest_leenhouts_flash.py`: 1006 additional mixture-flash rows added from `Leenhouts` open data.
- `t1_leakage_audit.py`: family-level CV on pure flash (122 families, 710 rows); 5-fold leaves-family-out.
- `tests/test_predict_blend.py`: wrapper API determinism, domain warnings, unknown compound fallback verified.
