# UGRP Hydrocarbons — Project Handoff

> Purpose of this file: let a new agent/developer (opencode) continue **exactly** where the previous
> work stopped. It states the goal, the physics, the data schema, what was built, what was verified,
> what is broken or missing, and the ordered next steps. Read sections 1–3 first, then jump to section 12.
>
> Tip: save/copy this file as `AGENTS.md` (or reference it from your opencode config) so it is loaded as context.

---

## 1. Goal

Predict **properties of hydrocarbon fuel blends** (jet / diesel / rocket kerosene such as RP-1) from the
**composition of the blend**, with a model that is *trustworthy on real measured blends*, then use it for
**inverse design** (find a blend that meets property targets).

Properties in scope, in priority order:

| Property | Role | Notes |
|---|---|---|
| Kinematic viscosity (cSt, any T) | core | Refutas blending is the physics baseline |
| Flash point (°C) | core | Le Chatelier + Raoult + Antoine baseline; **one-sided safety error matters** (over-predicting is the dangerous direction) |
| Pour point (°C) | core | PPBI baseline; weakest data and weakest physics |
| Density | auxiliary / needed | needed to convert mass↔volume fractions and dynamic↔kinematic viscosity |
| Freezing / cloud point, distillation (D86), smoke point, DCN/MON, heat of combustion | extension | partly available (see §8) |

The existing work (conference abstract *Shah et al., ChemCon 2026*; `Abstract_Shah_et_al_Chemcon2026_final.docx`
and `Fuel_Blend_Property_Methods_final.docx`, both provided by the user) is, as understood in this
conversation: a pure-component property stage (Chew et al. style) feeding a feed-forward network over ~40
component inputs, trained on **synthetic** blends whose labels come from the blending rule itself
(Refutas), reporting ~0.4 % error on held-out synthetic blends and ~3–5 % on a real RP-1 check.
**Verify these figures against the abstract before quoting them.**

**Key methodological conclusion reached in this conversation:** the 0.4 % number is circular (the net learns
the labelling rule). A robust model needs (a) real measured blend data, (b) physics baseline + learned
residual, (c) splits that hold out whole blend *systems*, (d) calibrated uncertainty. The whole data
pipeline in this repo exists to build (a) and to make (b)–(d) possible.

### 1.1 What the user wants next
1. Finish assembling one clean, unified, de-duplicated dataset (this repo).
2. Compute the physics baselines for all three core properties on real data and see where they fail.
3. Train a proper model (university CS-lab compute is available) — residual learning on top of physics,
   permutation-invariant mixture encoder, uncertainty, system-grouped validation.
4. Extend to other properties; then inverse design.

---

## 2. Physics reference (from `Fuel_Blend_Property_Methods_final.docx`)

Three composition bases; **use the right one per property**:

| Symbol | Meaning | Used for |
|---|---|---|
| `w_i` | mass fraction | viscosity blending, input composition |
| `x_i` | mole fraction (liquid) | flash point |
| `v_i` | volume fraction | pour point |

```
x_i = (w_i/M_i) / Σ_j (w_j/M_j)          v_i = (w_i/ρ_i) / Σ_j (w_j/ρ_j)     (all ρ at the SAME temperature)
```
All equation temperatures are absolute (K): `T_K = T_°C + 273.15`.

**Viscosity — Refutas.** All ν are **kinematic, cSt (mm²/s), same temperature**. Valid for ν ≳ 0.2 cSt.
```
VBN_i  = 14.534 · ln(ln(ν_i + 0.8)) + 10.975
VBN_mix = Σ w_i · VBN_i
ν_mix  = exp( exp( (VBN_mix − 10.975) / 14.534 ) ) − 0.8
```
**Pour point — PPBI** (empirical; poor for waxy/paraffinic blends). `T_i` in K:
```
PPBI_i = T_i^12.5      PPBI_mix = Σ v_i · PPBI_i      T_mix = PPBI_mix^0.08   (then −273.15 for °C)
```
**Flash point — Le Chatelier + Raoult + Antoine.** Flash point is the T where the cumulative hazard ratio H_T = 1.
```
log10 P_i^sat = A_i − B_i / (T + C_i)              (units: see §6.3 — this repo uses P in kPa, T in K)
y_i = x_i · γ_i · P_i^sat(T) / P_total             P_total = 101.325 kPa, γ_i = 1 (ideal)
H_T = Σ_i y_i / LFL_i  = Σ_i x_i γ_i P_i^sat(T) / (P_total · LFL_i)      (LFL_i as a mole fraction = vol% / 100)
solve f(T) = H_T − 1 = 0   (Newton–Raphson / secant; bisection if it diverges)
optional: LFL(T) = LFL_25 · (1 − 0.000721 (T − 25))   (Zabetakis, T in °C)
```
Per-component data these methods need: composition + class (aromatic/alkane), MW, density at a stated T,
Antoine A/B/C with units and valid range, Psat at 40 °C (cross-check), pour point, **LFL with its measurement
temperature**, kinematic viscosity at the temperature(s) of interest.

---

## 3. Environment

* User machine: Windows, PowerShell, venv at
  `C:\Users\kahaa\OneDrive\Desktop\UGRP_hydrocarbons\.venv`; project folder `...\UGRP_hydrocarbons\fuelblend_data`.
* Python 3.12.5; pandas 3.0.6; numpy 2.5.3; rdkit 2026.03.6; requests 2.34.2; lxml 6.1.3; openpyxl 3.1.5; scipy.
* **Move the project out of OneDrive** (≈800 small files + ~600 MB of raw data sync badly; deep zip paths can hit
  Windows' 260-char limit). The audit script warns about this.
* Extra packages needed for AJFTD tools (not yet in `requirements.txt` — **add them**):
  `beautifulsoup4` (offline HTML parser), `playwright` + `python -m playwright install chromium` (scraper).
* pandas 3 gotchas already hit: assigning mixed-type lists via `.loc[k, [cols]] = [...]` fails → assign per column.
* Windows/PowerShell: use `python script.py`; long-running scripts print progress; Ctrl-C handling exists in `fetch.py`.
* University CS-lab compute is available for model training.

---

## 4. Repository layout (current, v4/r3; the v1 `fuelblend_data/` layout no longer exists)

```
src/fuelblend/     package: api, baselines, data (loaders, ThermoML extractor, provenance), design, eval, features,
                   mixture, models, splits, thermo, train, uq
scripts/           t02_leakage_audit.py; v2/, v3/ (historic, v3 outputs quarantined); v4/ (current pipeline,
                   entry point scripts/v4/reproduce_all.py)
tests/             pytest suites; tests/quarantine/ holds pinned v3 tests; tests/data/ fixtures
configs/           YAML (e.g. configs/uncertainty_v4.yaml)      splits/   split definitions (splits/v4/)
data/              raw/, interim/, processed/, unified/ (measurements, compounds.csv), real_fuels/ (v4/ = current
                   61 jet fuels + 210 gasoline octane rows), quarantine/ (dropped v3 data)
results/           v2/, v3/ (historic), v4/ (current evidence, one folder per task: P1, P2, B, D, D2, E, F, L, R2, R3, SLE, SP, ...)
docs/              FULL_REPORT_v4.md, claims_register.md, model_decision.md, PROJECT_HANDOFF.md, PROJECT_SUMMARY_v4.md,
                   data_sources.csv, abstracts/, superpowers/ (plans)
extras/            EXTERNAL raw-data staging (old fuelblend_data repo copy, raw downloads); gitignored, not part of the repo
```
Other root items: `fuels.csv` (read by v2/v3 scripts), `Makefile`, `pytest.ini`, `requirements.txt`, `SHA256SUMS_v4.txt`
(written by `reproduce_all.py`), `archive/` (baseline manifests), `jobs/v4_gpu/`, `runs/` (gitignored). The v1 module names
(`build_unified.py`, `fetch.py`, ...) live on only inside `extras/` and in git history; `common.py` helpers now sit under `src/fuelblend/`.

## 5. Pipeline — commands in order

```powershell
pip install -r requirements.txt      # + beautifulsoup4 playwright
python fetch.py all                  # git, Zenodo, PDFs, ThermoML.  Flags: --skip zenodo | zenodo --list | --include REGEX | --max-gb N
python inspect_downloads.py          # unzip + INVENTORY.md (skips __MACOSX / ._* / ~$* junk)
python build_viscosity.py            # → data/processed/viscosity_{pure,mix}_hc.csv, compounds_hc.csv
python thermoml_extract.py           # → data/processed/thermoml_hc.csv  (a few min; prints "rows NOT kept, by reason")
python webbook_antoine.py            # → data/processed/antoine_best.csv, antoine_all_sets.csv
python ajftd_parse_html.py --compare ajftd_data.csv   # offline re-parse of saved AJFTD pages (run where debug/ lives)
python ingest_ajftd.py               # → data/processed/ajftd_fuels.csv, ajftd_composition.csv
python build_unified.py              # → data/unified/*   (1st pass)
python refit_antoine.py              # → data/processed/antoine_refit.csv   (needs unified measurements incl. ThermoML vapour pressure)
python build_unified.py              # (2nd pass so compounds.csv carries refit constants)
python audit.py --bundle             # → AUDIT.md + audit_bundle.zip (small; safe to share)
```
`fetch.py` is idempotent, resumable, rate-limit aware, verifies PDF headers and the ThermoML SHA-256.

---

## 6. Data schema (the contract — everything downstream depends on this)

### 6.1 Canonical properties and units

| `property` | unit | notes |
|---|---|---|
| `flash_point_C` | °C | closed vs open cup differs ~5–10 K; method kept in `ref`/AJFTD `method` where known |
| `viscosity_dyn_mPas` | mPa·s | all sources converted (see §9 unit traps) |
| `viscosity_kin_cSt` | cSt = mm²/s | needed by Refutas; = `viscosity_dyn_mPas / density[g/cm³]` at same T |
| `density_kg_m3` | kg/m³ | liquid phase only |
| `vapor_pressure_kPa` | kPa | pure liquids; used to refit Antoine |
| `freezing_point_C`, `pour_point_C`, `cloud_point_C` | °C | AJFTD / manual |
| `dcn`, `mon` | – | only with `build_unified.py --octane` |
| AJFTD extras | – | `distillation_T{0,5,10,20,…,95,100}_C`, `distillation_{loss,recovery,residue}_pct`, `smoke_point_mm`, `thermal_stability_spent_fuel_mL` |

### 6.2 `data/unified/measurements.csv` (every observation) and `measurements_dedup.csv` (one per condition)

| column | meaning |
|---|---|
| `source` | `leenhouts` \| `nist_tde` \| `thermoml` \| `ajftd` |
| `property` | canonical name above |
| `value` | in canonical unit |
| `T_K`, `P_kPa` | conditions (NaN if not applicable/reported) |
| `kind` | `pure` \| `binary` \| `multi` (≥3) \| `real_fuel` (AJFTD; no molecule list) |
| `n_comp`, `is_hc` | number of components; True iff **all** components are hydrocarbons |
| `system_id` | md5[:12] of the sorted component **InChIKey set** (`fuel:<id>` for real fuels). **Use as the GroupKFold key** |
| `comp_inchikeys` | JSON list, sorted |
| `x_mole`, `w_mass` | JSON lists aligned to `comp_inchikeys` (both always provided) |
| `fuel_id`, `ref`, `raw_value`, `raw_unit` | provenance (raw = before unit conversion) |

Extra columns in `measurements_dedup.csv`:

| column | meaning |
|---|---|
| `T_key`, `P_key`, `x_key` | binning keys: T to 0.5 K; P ≤150 kPa or missing → 101, else rounded to 100 kPa; x to 3 decimals |
| `value` | **median** of independent observations |
| `vmin`, `vmax`, `q25`, `q75` | spread |
| `n_obs`, `n_refs`, `n_sources`, `sources` | support (`sources` e.g. `leenhouts+nist_tde`) |
| `n_outlier` | observations >25 % from the median (needs ≥3 obs) |
| `conflict` | IQR (or range if 2 obs) exceeds tolerance: abs 2 °C flash, 3 °C pour/cloud/DCN/MON, 1 °C freezing; 5 % relative otherwise |

De-dup rules: (1) drop exact repeats (same source+ref+value) — ~⅓ of raw NIST rows; (2) group by
(property, system, composition, T, P); (3) median; (4) flag, never silently delete, outliers/conflicts.

### 6.3 `data/unified/compounds.csv` (registry; only compounds actually used)

`inchikey` (primary key), `smiles`, `inchi`, `formula`, `mw`, `is_hc`, `cls` (§4 classes or `other`), `n_carbon`,
`antoine_A`, `antoine_B`, `antoine_C`, `antoine_Tlo`, `antoine_Thi`, `antoine_source` (`refit_measured` | `webbook`),
`antoine_covers_window` (True iff valid over 313.15–373.15 K), `antoine_unit`, `lfl_vol_pct`, `lfl_T_C`.

> **Antoine unit convention in this repo: `log10(P/kPa) = A − B/(T/K + C)`.** WebBook gives log10(P/**bar**):
> `A_kPa = A_bar + 2`. `build_unified.py` already converts WebBook sets; `refit_antoine.py` outputs kPa natively.

### 6.4 Intermediate / source schemas

* **Leenhouts benchmark** (`data/interim/leenhouts_benchmark/datasets/datasets/...`; Zenodo 19914760):
  `flashpoint_pure_exp.csv` (585 rows: `value` [unit string `Celcius`, sic], `component_inchi_0`),
  `flashpoint_mix_exp.csv` (1006; `fraction_type`=`molar`, `component_inchi_0..3`, `fraction_0..3`),
  `viscosity_pure_exp.csv` (3027) / `viscosity_mix_exp.csv` (35,938): `value` is **ln(mPa·s)** (`unit` = `ln(mPa*s)`), `temperature` in K.
  `dcn_mix_exp.csv` (484), `mon_mix_exp.csv` (292, up to ~34 components). The tables also carry ideal-mixing/excess
  columns (see `data/interim/INVENTORY.md`). **Never use `DatasetsFuelProperties.xlsx` as ground truth — 50,000 SYNTHETIC rows.**
  `PureFuels.xlsx` = pure fuel table behind those synthetic rows.
* **NIST/TDE viscosity** (CheMixHub `datasets/nist-logV/processed_data/processed_NISTlogV.csv`, 239,200 rows; raw
  `NIST_Visc_Data.csv` is inside Zenodo 8042966 `pretrained_models.zip`, 239,200 rows, Visc in Pa·s).
  **The `value` column is labelled "Log viscosity / cP" but is ln(η / Pa·s).**
* `viscosity_pure_hc.csv`: `compound_id, smiles, cls, n_carbon, mw, T_K, P_kPa, eta_mPas, ln_eta_Pas, Visc_Unc, Ref_ID`.
  `viscosity_mix_hc.csv`: `smiles_list, x_mole, w_mass, classes (JSON), n_comp, T_K, P_kPa, eta_mPas, ln_eta_Pas, Visc_Unc, Ref_ID`.
* `thermoml_hc.csv` (new extractor): `file, prop, prop_phase, value, T_K, P_kPa, n_comp, inchis (JSON), names (JSON),
  comp_type (pure|mole fraction|mass fraction|volume fraction), comp_values (JSON)`. Native ThermoML units
  (Pa·s, m²/s, kg/m³, kPa); `build_unified.py` converts. **An old-format file (no `prop_phase`) is refused.**
* `antoine_best.csv` / `antoine_all_sets.csv` (WebBook, **bar**): `T_lo, T_hi, A, B, C, reference, comment, name, cas, smiles,
  page_name, coverage, covers_window, psat_40C_kPa`.
* `antoine_refit.csv` (kPa): `inchikey, smiles, n_points, A, B, C, T_min, T_max, rms_log10P, max_abs_log10P, covers_window,
  anchored, webbook_maxdiff, status` (`ok` | `check: rms > 0.03` | `skipped: <6 points or span <30 K` | `failed: …`).
* `data/raw/lfl/lfl_template.csv` → fill and save as `lfl.csv`: `inchikey|smiles, lfl_vol_pct, T_C, source, note`.
* **AJFTD** — scraper CSV (`ajftd_data.csv` / `ajftd_all_metrics.csv`): `fuel_id, fuel_name, fuel_label, record_url, page_title, category,
  metric, table_id, measurement_index, Value, Unit, temperature_value, temperature_unit, volume_evaporated_value,
  volume_evaporated_unit, specification_limits, info, Test Date, Test Method`. Non-values: `No Test`, `No search`, `Not Found`, blank.
  → `ajftd_fuels.csv` (long, numeric only): `fuel_id, fuel_label, category, method, test_date, qualifier, info, src_file, src_col, property,
  value, unit_std, T_K, derived`. → `ajftd_composition.csv` (wide): `fuel_id, fuel_label, comp::<metric>[<unit>]…`.
  API gravity → density at 15.56 °C: `ρ = 141.5/(API+131.5) · 999.016` kg/m³ (`derived = "from API gravity"`).

---

## 7. Data sources (where everything comes from)

| Source | What | Access | Status |
|---|---|---|---|
| Leenhouts et al. benchmark | curated pure + mixture flash/viscosity/DCN/MON (CC BY/MIT) | Zenodo `19914760` (`datasets.zip`, `splits.zip`); code `gitlab.kuleuven.be/creas/vermeiregroup/mixture_evaluation`; paper arXiv 2605.29698 | downloaded, unzipped, row counts verified |
| CheMixHub (NeurIPS 2025) | ~500k mixture points incl. NIST viscosity | `github.com/chemcognition-lab/chemixhub` | cloned; `nist-logV` used |
| Bilodeau et al. | raw NIST viscosity of binary liquid mixtures | Zenodo `8042966` → `pretrained_models.zip` | downloaded (raw table inside) |
| ThermoML archive (NIST TRC) | measured density/viscosity/vapour pressure, 5 journals through 2019 | doi 10.18434/mds2-2422, file `ThermoML.v2020-09-30.tgz` (189,433,115 B, sha256 `231161b5…58d7ec2`) | downloaded & verified; **needs re-extraction with the fixed parser** |
| NIST WebBook | Antoine constants (bar, K) | scrape `cbook.cgi?ID=C<cas>&Mask=4` | run once; **valid ranges don't reach 40–100 °C for C10+ alkanes** |
| AJFTD (Illinois/FAA ASCENT) | real jet-fuel property database | login required (access granted) | 81 fuels scraped; thin (see §8.3) |
| NIST RP-1 (NISTIR 6646; Energy & Fuels 2012 doi 10.1021/ef3004009) | 11 orthogonal RP-1 blends — the best real multi-component validation | PDF downloaded; **blend data (SI) not yet obtained** | TODO |
| CRC reports 647, AV-27-18, AV-11-09 | real-fuel property tables (pour, freezing, viscosity, flash) | crcao.org PDFs | downloaded; **tables not yet transcribed** |
| Zabetakis Bulletin 627 | LFL data (paraffins, aromatics, T and MW effects) | UNT digital library PDF | downloaded; **LFL not yet transcribed** |
| Li et al. 2020 (arXiv 2002.11315) | 10,575 pure-compound flash points | SI | PDF downloaded; **dataset not obtained** |
| Chew et al. 2024 (doi 10.1186/s13321-024-00820-5) | pure viscosity vs T dataset (dynamic, cP) | SI PDF downloaded | dataset location unchecked |
| DTIC AFRL-TR-2020-0017, CRC Handbook | jet fuel properties | apps.dtic.mil blocks scripts (HTTP 403 → HTML) | optional; download in a browser; delete the fake 0.0 MB files |
| DIPPR | proprietary | – | not used |

---

## 8. Current state

### 8.0 v4/r3 state (tag `v4-r3-final`; latest commit "v4 r3: CRC 647 compositions from measured hydrocarbon types (Table G.1), re-scored")
Source of truth: `docs/FULL_REPORT_v4.md` (sections 0, 1, 11, 13) and `docs/claims_register.md`. Facts as reported there:
* The v3 corpus and design optimizer were formula-derived and are quarantined (`data/quarantine/`); 92 of 153 v3 fuels dropped, 61 kept (all jet).
* v4 evidence: 57 CRC 647 jet fuels (compositions reconstructed; r3 uses measured hydrocarbon types from Table G.1) + 4 other jet fuels,
  210 gasoline octane blend rows, no diesel. Only POSF-10325 is a headline GC fuel (n = 1).
* Designs: 0 of 13 certified (D2 re-run with known-composition lab margins); diesel unvalidated; freezing models show no skill over a constant.
* r2 reconstruction flipped CRC flash bias to the unsafe (over-prediction) side (+4.2 C raw, dev); causes not isolated.
* Open: ThermoML `phase_of()` default-mode bug (recorded, not fixed), thin locked sets, cold-viscosity lab margin is `inf`; full list in report section 11 and 13.4.
* Tests: pytest 308 passed / 42 failed at the v4 baseline (all pre-existing missing-deliverable asserts; `results/v4/baseline_pytest.txt`).

### 8.1 HISTORIC (v1-era, 2026-10-29 and earlier; numbers below predate v4 and several rest on v3-style or contaminated data; do not quote without checking `docs/claims_register.md`)

#### 8.1.0 Session 2026-10-29 additions
* **LFL measured rows**: `data/raw/lfl/lfl_measured.csv` (12 HC compounds hand-read from Zabetakis 627 App. A-1 OCR lines, ratio-gated vs 0.55*Cst); `fill_lfl.py` → 695 compounds, 25 measured / rest estimated. Pure-flash model now runs on ~84-152/902-923 records; Antoine coverage is the limiter.
* **flashpoint_li** (Sun/Krakauer/Morgan 2020, figshare 9275210, CC BY): 14,696 rows, flash in K, loader `load_flashpoint_li()`; +550 HC pure flash records (dedup). Duplicates collapse to median. Some suspect pure values remain (e.g. C25 paraffin at 120 C) - flagged not deleted.
* **fortin2012** (NIST copy, `data/raw/rp1/fortin2012_rp1_blends.pdf`): Tables 1-4 extracted → 252 rows: distillation T0.025..T85-90 for 11 orthogonal RP-1 + 1 RP-1 4572 + 2 RP-2 real fuels (loader `load_fortin()`). THE 11-blend compositions are NOT in this paper (recipes confidential; per NISTIR 6646 phase-1 covers one sample only).
* **RP-1 NISTIR 6646** (`data/raw/rp1/` sub-agent): one characterized RP-1 sample (33 components w/ SMILES + composition tiers), surrogate 14-comp, 1480 property rows (density ~83-484 kg/m3 incl high P/T, kinematic visc 1.13-7.67 cSt, thermal conductivity).
* **baselines.py**: viscosity physics baseline + residual GBM on 15,281 mixture records (GroupKFold system_id): loglin MAPE 8.8 %, **Refutas 4.8 %** (binary winner), GBM residual only helps multi-comp (4.85 % vs 6.03 %; over-corrects binaries because only 319 unique systems) → use Refutas for binaries, GBM residual only for n_comp>=3. See BASELINE_VISCOSITY_REPORT.md.
* **Flash baseline** (after exp) on 147 predicted records: overall MAE 7.17 C; **binary 1.59 C / 98 % within 5 C**; pure 10.7 C (pubchem tail). Domain filter for flash now (−90..250 C).
* **CRC ingest final**: crc647 4,559 raw rows (57 aviation fuels: pour 57, freezing 57, flash 57, D86 grids, ±20/−40 °C kinematic viscosity ×171, heat, conductivity, speed of sound), av2718 869 (aromatics, D2887 distillation, viscosity+density, vapour pressure, gas solubility), av1109 3. After build: **pour_point_C 57 / freezing_point_C 78 / viscosity_kin_cSt 361 / thermal_stability 53** in dedup (~70k raw rows total).
* CRC parse + Li/LiS1 + audit rerun: see the sub-agent outputs in data/raw/pour_point/ and data/raw/flashpoint_li/.

#### 8.1.1 Audit (last bundle): FAIL 0 · MISSING 6 · WARN 8 · OK 49
Missing: AJFTD export folder (since replaced by scraped data), RP-1 blends, pour/cloud tables, Li et al. SI, **LFL (0/611 compounds)**,
real-fuel records in unified (before AJFTD ingest). Warnings: OneDrive; 1 macOS junk file; 2 inspect read errors (2 empty files, 1 Excel
lock `~$…xlsx` — harmless); Antoine 22/24 found (methylcyclohexane, n-eicosane have no Antoine table); 20/22 sets don't cover 40–100 °C;
2 fake DTIC PDFs; WebBook `page_name` held the site banner (bug fixed in code — rerun `webbook_antoine.py`).

#### 8.1.2 Hydrocarbon-only open data (usable for training/validation)

| Property | Pure | Mixtures | Notes |
|---|---|---|---|
| Flash point | 244 of 585 (−136…204 °C) | **126** of 1006 (123 binary, 3 ternary) | binary alkane/ternary show real non-ideality (e.g. a decane/dodecane/tetradecane mix flashes ~9 °C below ideal) |
| Viscosity (Leenhouts, NIST-derived) | 692 | 2,752 binary | 293–323 K only |
| Viscosity (own NIST build) | 4,832 raw rows | 10,947 binary raw rows | 109 compounds; ≈1/3 raw rows are exact repeats |
| DCN / MON | – | 346 / 106 (MON up to 34 components) | not target properties |
| Pour point | none | none | no open curated dataset exists |
| **≥3-component hydrocarbon data** | | **essentially none** (3 flash rows) | → real-fuel sources (AJFTD, RP-1, CRC) are the critical path |

NIST/TDE noise: of 4,341 pure points with ≥3 independent refs only 13 (0.3 %) deviate >25 % from the group median, but one reference
reports benzene@298 K = 0.395 mPa·s while ~70 refs give ≈0.60. Dedup (median + flags) handles it.
Real run: 15,779 raw NIST rows → 11,030 de-duplicated records, 30 genuine conflicts.

#### 8.1.3 AJFTD (real fuels)
* 81 fuels only (ids 1 and 365–462), tests dated 2015-04-17 … 2016-11-29. **Check the site's "Showing x to y of N entries": if N > 81,
  pagination stopped early or filters were still active.**
* Numeric coverage of scraped metrics: API gravity 81; distillation T0…T100 ≈80 fuels (561 numeric rows); smoke point 82; aromatics (% vol)
  82 (8.6–22.6); diaromatics 80; existent gum 82; sulfur mercaptan 81; thermal stability 81.
  **Flash point 1 fuel, kinematic viscosity 1 fuel, density (kg/m³) 1 fuel** — the rest are "No Test". Other metrics on the page may exist
  but the scraper only exports aliases/composition/volume metrics → run `ajftd_parse_html.py` on the saved HTML to see everything.
* Fuel name was the site banner for every record (fixed: real label e.g. `(02239-070416)T108`).
* Possible use of the distillation data: estimate flash point from D86 T10 via a published correlation — an *estimate*, validation only.

#### 8.1.4 ThermoML (user's run, **before** the parser fix → contaminated)
Unified report from that data: viscosity records (HC) 18,322 with min value −0.063 (excess viscosity!), vapour pressure up to 1.4×10⁶ kPa,
density 0.52 kg/m³ at up to 873 K (gas phase), only 39 % of rows had T, 1,728 spurious conflicts, 6,451 rows skipped. **Discard that
`thermoml_hc.csv` and re-extract.**

#### 8.1.5 Antoine refit (also fitted on the contaminated data → redo)
220 compounds with vapour-pressure data: 143 ok, 43 `check`, 34 skipped; 94 ok-and-cover 40–100 °C. n-alkanes (kPa form):

| alkane | status | note |
|---|---|---|
| C5, C8, C9 | ok but T-range doesn't cover 40–100 °C window | |
| C6, C7, C10, C14, C16, C21–23 | ok, cover window | C10: P(25 °C)≈0.18 kPa ✓ literature; matches WebBook to 0.002 log units |
| C11 | ok, only 10 points at 393–459 K | doesn't cover window |
| C12, C20, C24 | `check` (rms 0.06–0.07 log10) | inspect |
| C13, C15, C17, C19 | **no data** | candidates for interpolation of A,B,C vs carbon number |
| C18 | skipped (5 points) | |

#### 8.1.6 Physics-input coverage (flash baseline), from the user's unified report (contaminated ThermoML)
370 HC flash records, 248 distinct components; Antoine available 52/248, in-range 29; **LFL 0/248 → 0/370 records can run the physics flash model.**
Viscosity: 17,020 of 18,322 HC records have all components with a measured pure density (needed for kinematic/Refutas).

---

## 9. Bugs found and fixed (do not reintroduce)

1. **ThermoML `Constraint` bug**: fixed T/P live in `Constraint > ConstraintID > ConstraintType`, not a direct `ConstraintType` child →
   61 % of temperatures lost. Also: exact property-name whitelist (no `Excess viscosity`), `PropPhaseID/ePropPhase` must be Liquid
   (crystal vapour pressures dropped), composition variables with `VarPhaseID/ePhase` ≠ Liquid (VLE vapour fractions) ignored, pressure
   >1000 kPa dropped (flag `--max-pressure-kpa`), `P` is part of the dedupe key, ternary last component implied, molality etc. rejected.
   Verified against the official XSD (ThermoPyL `thermoml_schema.py`) and a real record (`tests/data/je8006138.xml`).
2. **Antoine refit crash** (`Initial guess is outside of provided bounds`): falling/flat pressure data → negative B start → clipped starts,
   correlation guard, per-compound `failed:` status.
3. **Units**: NIST logV = ln(Pa·s) despite label; Leenhouts viscosity = ln(mPa·s); WebBook Antoine in bar (+2 → kPa); ThermoML
   Pa·s/m²/s need ×1e3/×1e6; flash unit string `Celcius`.
4. `fetch.py`: HTML consent page from DTIC saved as a "PDF" → now validated; `.part` resume; 429 back-off; per-file deadline; Ctrl-C summary.
5. Zenodo zips contain macOS `__MACOSX/._*` shadow files → hundreds of fake parse errors → skipped everywhere.
6. `compounds.csv` had an unquoted name containing commas (`1,2,4-trimethylbenzene`).
7. WebBook `page_name` returned the site banner; now `<h1 id="Top">` → `<title>`.
8. Scraper: every fuel named "Alternative Jet Fuels"; no resume; no pagination cross-check.
9. (Process lesson) A synthetic test built on a guessed XML structure passed while real data failed — keep the real-record regression test.

---

## 10. Known limitations / unverified

* New `thermoml_extract.py` verified on one real record + synthetic cases, **not yet on the full 189 MB archive**.
* `ajftd_parse_html.py` verified on fixtures built like the pages the scraper parses; **not on real saved pages** (it prints a warning if a page has no metric tables).
* `webbook_antoine.py` name fix not yet re-run live; CAS numbers in `compounds.csv` were typed from memory — check `page_name` vs `name`.
* The open data covers pure + binary only; no validation on realistic ≥3-component fuel blends yet.
* No LFL values anywhere yet; pour point has no data.
* Antoine for heavy odd alkanes absent; refit not yet repeated after the ThermoML fix.
* All development was done in a sandbox that could only reach GitHub; every network-dependent step was tested against mock servers or on the user's machine.

---

## 11. Tests (`python tests/<name>.py`, plain asserts, no pytest needed)

| suite | covers |
|---|---|
| `test_offline.py` | composition conversions, hydrocarbon filter/classes, Antoine (benzene ≈24.3 kPa @40 °C), WebBook parsing & name, Zenodo file list |
| `test_fetch_server.py` | local mock server: 429 back-off, stall timeout, resume, size guard, Ctrl-C, 403→browser headers, sha256, HTML-as-PDF rejection |
| `test_thermoml.py` | real record + synthetic: constraints, excess viscosity, VLE vapour composition, ternary completion, pressure/phase filters, end-to-end tarball |
| `test_unified.py` | schema, InChIKey ordering, units, mass/mole conversion, dedupe/conflict, pressure separation, physics join, LFL pickup, old-ThermoML refusal |
| `test_refit.py` | Antoine refit recovers synthetic truth to <2.5 % at 40/100 °C; bad data fails per-compound |
| `test_ajftd.py`, `test_ajftd_html.py`, `test_scraper_patch.py` | generic wide/long ingest; scraper-format ingest; offline HTML parsing; scraper name/resume helpers |

---

## 12. v4 status and v5 (read first)

v4 (branches `v4-integrity`, `v4-review-fixes`; final tag `v4-r3-final`) quarantined the formula-derived v3 corpus and designs and rebuilt the evidence on `data/real_fuels/v4/`. Entry points:
* `docs/FULL_REPORT_v4.md` (state, numbers, open items in section 11, review fixes in section 13)
* `docs/claims_register.md` (claim-level status, CLM-xx)
* `docs/model_decision.md` (sections 5-6) and `docs/PROJECT_HANDOFF.md` (section 15)
Reproduce with `python scripts/v4/reproduce_all.py` (venv python); pytest baseline in `results/v4/baseline_pytest.txt`.

v5: work continues on branch `v5` (current). A written v5 plan is not yet in the repo (`docs/superpowers/plans/` holds only the v4 plans); check there and `docs/data_source_leads_2026-11-09.md` for candidate new data before assuming scope.

## 12b. (v1-era plan, superseded by section 12 above) Exactly where we stopped

The last deliverable was `fuelblend_update.zip` (14 files: `thermoml_extract.py`, `build_unified.py`, `refit_antoine.py`, `ingest_ajftd.py`,
`ajftd_parse_html.py`, `ajftd_scraper.py`, `audit.py`, and tests). **Confirm those are applied**, then:

1. `python tests/test_*.py` (all 8 should pass).
2. `python thermoml_extract.py` → check the "rows NOT kept, by reason" table and that temperature coverage ≈100 %.
3. `python build_unified.py` → `python refit_antoine.py` → `python build_unified.py`. Re-read `UNIFIED_REPORT.md`.
4. `python ajftd_parse_html.py --compare ajftd_data.csv` (from the folder with `debug/`); move the old `ajftd_data.csv` out of
   `data/raw/ajftd/` so fuels aren't double counted; inspect `ajftd_metrics_inventory.csv` for what exists (freezing/pour/viscosity at T?).
5. Check AJFTD "of N entries" vs 81; clear filters and rerun the scraper if N is larger (it now resumes from saved HTML).
6. `python ingest_ajftd.py` → `python build_unified.py` → `python audit.py --bundle`.
7. **LFL**: fill `data/raw/lfl/lfl_template.csv` from Zabetakis Bulletin 627 (PDF already in `data/raw/pdfs/`; can be extracted programmatically), save as `lfl.csv`.
8. Get the RP-1 11-blend data, Li et al. flash-point SI, CRC 647 / AV-11-09 pour/freezing tables → `data/raw/rp1|flashpoint_li|pour_point/` and write ingesters in the same canonical schema.
9. For missing heavy n-alkane Antoine data (C13/15/17/19, C11 low-T) interpolate A,B,C vs carbon number (flag `antoine_source = interpolated`).

### 12.1 Then the modelling work (not started)
1. **Baselines** (new module, e.g. `baselines.py`): Refutas (needs kinematic ν → use ThermoML/pure density at same T), PPBI, Le Chatelier flash point
   (Antoine from `compounds.csv`, LFL, γ=1), plus ideal mixing. Evaluate on `measurements_dedup.csv`, report error by property, `kind`, class, T.
2. **Residual model**: `property = physics_baseline + learned_correction`; start with gradient boosting on blend features
   (mass/mole fractions, class fractions, mean/std of n_carbon, MW, ideal values), then a **DeepSets / attention mixture encoder**
   (component embedding from fingerprints/descriptors/GNN × fraction, pooled, shared trunk, per-property heads, masked multi-task loss).
   Learn an activity-coefficient-like correction for flash point; excess-viscosity correction for viscosity.
3. **Validation**: GroupKFold on `system_id`; leave-molecule-out; keep real RP-1 blends and a slice of AJFTD fuels as a never-tuned test set.
   Metrics: viscosity in log space (relative error), flash point absolute °C **plus one-sided over-prediction error**, pour point °C.
   Benchmark context (Leenhouts): flash RMSE ≈9–15 °C and the ideal baseline ≈13 °C on strict splits; ~10³ unique training mixtures are needed to beat ideal mixing.
4. **Uncertainty**: deep ensembles/MC-dropout + calibrated intervals (check 90 % coverage); propagate pure-component uncertainty.
5. **Real fuels** (AJFTD/RP-1/CRC) have no molecule lists: represent them by class × carbon-number distributions (GC×GC if available),
   use as validation and for a distribution-input model.
6. **Inverse design**: constrained optimisation on the composition simplex (Bayesian opt / CMA-ES / gradient), multi-property targets,
   cost/availability constraints, penalise high-uncertainty regions, verify experimentally or against an independent source.
7. Extension properties: density (ThermoML), freezing/cloud (AJFTD, CRC AV-11-09; wax thermodynamics), distillation, DCN/MON (Leenhouts).

---

## 13. Conventions for the next agent

* Keep **hydrocarbon-only** as the default filter (`is_hc`); use `--all-molecules` only deliberately.
* Always train/evaluate from `measurements_dedup.csv`, never raw rows; carry `conflict`, `n_outlier`, `n_refs` as quality weights/filters.
* Never mix units across sources; all conversions happen in the loaders (`build_unified.py`), nowhere else.
* Every new source gets: a loader to the §6 schema, a unit test with a fixture that mirrors the **real** file structure, and an audit check.
* Never commit `ajftd_session.json`, raw AJFTD pages, or licensed data; `audit_bundle.zip` contains no raw datasets by design.
* When a result looks too good, look for circular labelling (synthetic data, the blending rule used as the label).
* Report honestly what was run on real data versus fixtures.

---

## 14. Quick reference — key facts

* Real fuel property targets live at: flash 38 °C min (Jet A, D56) / viscosity limit 8 cSt max at −20 °C / density 775–840 kg/m³ (AJFTD spec limits).
* WebBook Antoine ranges for C10+ n-alkanes start ~95–190 °C → never use them unrefit at 40–100 °C.
* Leenhouts viscosity: ln(mPa·s); NIST logV: ln(Pa·s); difference 6.908.
* `system_id` = group key; `kind` ∈ {pure, binary, multi, real_fuel}.
* Files to read first in the repo: `build_unified.py` (schema + loaders), `thermoml_extract.py`, `ingest_ajftd.py`, `refit_antoine.py`, `audit.py`, `data/unified/UNIFIED_REPORT.md`.
