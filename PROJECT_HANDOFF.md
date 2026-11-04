# UGRP Hydrocarbons — Project Handoff

> Purpose of this file: let a new agent/developer (opencode) continue **exactly** where the previous
> work stopped. It states the goal, the physics, the data schema, what was built, what was verified,
> what is broken or missing, and the ordered next steps. Read sections 1–3 first, then jump to section 12.
>
> Tip: save/copy this file as `AGENTS.md` (or reference it from your opencode config) so it is loaded as context.
>
> **v2 (2026-10-29): scope extended to seven targets (viscosity, density, flash point, pour point, GCV, octane, smoke point).
> Read §17 first: it supersedes the property table in §1 and the decision list in §16.3-C. The data-gathering
> work is itemised in `DATA_ACQUISITION_TASKS.md` (task ids T00–T42).**

---

## 1. Goal

Predict **properties of hydrocarbon fuel blends** (jet / diesel / rocket kerosene such as RP-1) from the
**composition of the blend**, with a model that is *trustworthy on real measured blends*, then use it for
**inverse design** (find a blend that meets property targets).

Properties in scope (**v2, see §17**). Targets (seven): kinematic viscosity, liquid density, flash point, pour point,
GCV (gross calorific value), octane number (RON/MON), smoke point (ASTM D1322, via the threshold sooting index TSI).
Auxiliary / extension: vapour pressure (feeds flash point and COSTALD density), DCN, freezing/cloud point, distillation (D86).

| Property | Role | Notes |
|---|---|---|
| Kinematic viscosity (cSt, any T) | target | Refutas blending is the physics baseline |
| Liquid density (kg/m³, any T) | target | ideal volume additivity baseline; also needed to convert mass↔volume fractions and dynamic↔kinematic viscosity |
| Flash point (°C) | target | Le Chatelier + Raoult + Antoine baseline; **one-sided safety error matters** (over-predicting is the dangerous direction) |
| Pour point (°C) | target | PPBI baseline; weakest data and weakest physics |
| GCV (MJ/kg) | target | mass-weighted additivity is expected to be near-exact for liquids; the hard part is real fuels without composition and the measurement basis (§17.3) |
| Octane number (RON, MON) | target, **gasoline domain only** | volume-linear baseline plus learned non-linear residual; RON data is currently absent (§17) |
| Smoke point (mm, ASTM D1322) | target, **middle-distillate domain only** | modelled through TSI (blending is done on TSI, then converted back to mm with lamp constants taken from the source paper); pure-compound data are scarce, real-fuel smoke points exist (AJFTD, 82 fuels) but without composition (§17.3) |
| Vapour pressure, DCN, freezing/cloud, distillation | auxiliary / extension | partly available (see §8, §16) |

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

## 4. Repository layout (all in `fuelblend_data/`)

```
README.md  requirements.txt  compounds.csv  common.py
fetch.py                 download everything automatable (git, Zenodo, PDFs, ThermoML) with manifest/sha256
inspect_downloads.py     unzip Zenodo zips (skips macOS junk), write data/interim/INVENTORY.md
build_viscosity.py       hydrocarbon-only NIST/TDE viscosity tables from CheMixHub (pure + binary)
thermoml_extract.py      hydrocarbon liquid viscosity/density/vapour-pressure from the ThermoML archive
webbook_antoine.py       scrape Antoine sets from NIST WebBook (bar→kPa, range checks)
ingest_ajftd.py          AJFTD exports → canonical long table (+ composition table); `--describe` mode
ajftd_scraper.py         YOUR Playwright scraper (patched: resume, retries, real fuel name)
ajftd_parse_html.py      re-parse saved debug/detail_*.html offline for ALL metrics (no site access)
build_unified.py         ONE canonical schema: compounds / measurements / measurements_dedup / report
refit_antoine.py         refit Antoine constants from measured vapour pressures (valid at 40–100 °C)
audit.py                 checks everything, writes AUDIT.md, `--bundle` → audit_bundle.zip
tests/                   8 suites (see §11); tests/data/je8006138.xml = a REAL ThermoML record
data/raw/    git/ zenodo/ pdfs/ thermoml/ ajftd/ lfl/ rp1/ pour_point/ flashpoint_li/  (+ MANIFEST.json)
data/interim/ unzipped Zenodo, INVENTORY.md, read_errors.txt, ajftd_description.txt, ajftd_unmapped.txt
data/processed/ viscosity_*_hc.csv, compounds_hc.csv, thermoml_hc.csv, antoine_*.csv, ajftd_fuels.csv, ajftd_composition.csv
data/unified/ compounds.csv, measurements.csv, measurements_dedup.csv, UNIFIED_REPORT.md
```
Also present on the user's machine (not in the generated repo): `ajftd_scraper.py` original, `debug/detail_<id>.html`
(81 saved AJFTD pages), `ajftd_data.csv/.json`, `ajftd_session.json` (**login cookies — sensitive, never commit/share**).

`common.py`: `mol_from(smiles|inchi)`, `is_hydrocarbon(mol)` (only C/H, neutral, ≥1 C),
`hc_class(mol)` → `aromatic | alkene_alkyne | cycloalkane | iso_paraffin | n_paraffin`, `carbon_number`, `mol_weight`,
`mole_to_mass`, `mass_to_mole`, `mass_to_volume`.

---

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

## 8. Current state (measured numbers)

### 8.1 Audit (last bundle): FAIL 0 · MISSING 6 · WARN 8 · OK 49
Missing: AJFTD export folder (since replaced by scraped data), RP-1 blends, pour/cloud tables, Li et al. SI, **LFL (0/611 compounds)**,
real-fuel records in unified (before AJFTD ingest). Warnings: OneDrive; 1 macOS junk file; 2 inspect read errors (2 empty files, 1 Excel
lock `~$…xlsx` — harmless); Antoine 22/24 found (methylcyclohexane, n-eicosane have no Antoine table); 20/22 sets don't cover 40–100 °C;
2 fake DTIC PDFs; WebBook `page_name` held the site banner (bug fixed in code — rerun `webbook_antoine.py`).

### 8.2 Hydrocarbon-only open data (usable for training/validation)

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

### 8.3 AJFTD (real fuels)
* 81 fuels only (ids 1 and 365–462), tests dated 2015-04-17 … 2016-11-29. **Check the site's "Showing x to y of N entries": if N > 81,
  pagination stopped early or filters were still active.**
* Numeric coverage of scraped metrics: API gravity 81; distillation T0…T100 ≈80 fuels (561 numeric rows); smoke point 82; aromatics (% vol)
  82 (8.6–22.6); diaromatics 80; existent gum 82; sulfur mercaptan 81; thermal stability 81.
  **Flash point 1 fuel, kinematic viscosity 1 fuel, density (kg/m³) 1 fuel** — the rest are "No Test". Other metrics on the page may exist
  but the scraper only exports aliases/composition/volume metrics → run `ajftd_parse_html.py` on the saved HTML to see everything.
* Fuel name was the site banner for every record (fixed: real label e.g. `(02239-070416)T108`).
* Possible use of the distillation data: estimate flash point from D86 T10 via a published correlation — an *estimate*, validation only.

### 8.4 ThermoML (user's run, **before** the parser fix → contaminated)
Unified report from that data: viscosity records (HC) 18,322 with min value −0.063 (excess viscosity!), vapour pressure up to 1.4×10⁶ kPa,
density 0.52 kg/m³ at up to 873 K (gas phase), only 39 % of rows had T, 1,728 spurious conflicts, 6,451 rows skipped. **Discard that
`thermoml_hc.csv` and re-extract.**

### 8.5 Antoine refit (also fitted on the contaminated data → redo)
220 compounds with vapour-pressure data: 143 ok, 43 `check`, 34 skipped; 94 ok-and-cover 40–100 °C. n-alkanes (kPa form):

| alkane | status | note |
|---|---|---|
| C5, C8, C9 | ok but T-range doesn't cover 40–100 °C window | |
| C6, C7, C10, C14, C16, C21–23 | ok, cover window | C10: P(25 °C)≈0.18 kPa ✓ literature; matches WebBook to 0.002 log units |
| C11 | ok, only 10 points at 393–459 K | doesn't cover window |
| C12, C20, C24 | `check` (rms 0.06–0.07 log10) | inspect |
| C13, C15, C17, C19 | **no data** | candidates for interpolation of A,B,C vs carbon number |
| C18 | skipped (5 points) | |

### 8.6 Physics-input coverage (flash baseline), from the user's unified report (contaminated ThermoML)
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

## 12. Exactly where we stopped — do this next

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

---

## 15. Addendum — state after `audit_bundle_60f5ef` (supersedes conflicting numbers in §8)

**What the latest real run showed**
* ThermoML re-extraction with the fixed parser worked: **46,747** hydrocarbon rows, **100 %** with temperature (density ≈20.7k, viscosity ≈12.2k, vapour pressure ≈13.8k rows).
* Unified: 66,580 raw → **45,500** dedup records; conflicts fell from 1,728 to **253**; ≥3-component hydrocarbon records now **1,703** (density 910, viscosity 625, vapour pressure 165, flash 3).
* Antoine refit (no crash): 210 compounds → **158 ok**, 20 `check`, 32 skipped; **111 ok fits cover 40–100 °C** (was 94). n-alkanes good: C6, C7, C10, C14, C16; `check`: C12, C20; no data: C13, C15, C17, C18, C19; C11 only >393 K.
* Audit: FAIL 0 · MISSING 4 · WARN 8 · OK 63.
* **AJFTD: only 18 fuels were parsed** (240 rows) — the `debug/` folder holds only 18 saved pages, not 81. New metrics seen: `freezing_point` (−50 … −42.9 °C), fuel-details tables (900 cells), `doctor_test` (text POSITIVE/NEGATIVE — currently unmapped). Real-fuel flash point: 0.
* Unified ranges still contained non-fuel rows from real ThermoML/NIST data (negative viscosity, 25 kg/m³ "liquid" densities at 777 K, up to 1.4×10⁶ kPa "vapour pressures", plus 4 NIST rows at 8.8×10⁵ and ~4×10⁻⁴ mPa·s) → handled by the new domain filter.

**New/changed files since §4** (all tested, see `tests/`)
* `fill_lfl.py` (new): writes `data/raw/lfl/lfl.csv`. Priority: your `lfl_measured.csv` > Zabetakis Table 2 n-alkanes C1–C16 (OCR-read; nonane/decane cells reconstructed; **verify against PDF page 21**) > estimate `0.55·Cst` (Jones rule; reproduces C6–C16 within 10 %; ≈2 K flash-point effect per 15 % LFL error). `compounds.csv` gains `lfl_source`.
* `baselines_flash.py` (new): Le Chatelier + Raoult + Antoine on `measurements_dedup.csv`; writes `baseline_flash.csv` and `BASELINE_FLASH_REPORT.md` (errors by pure/mixture, Antoine validity, LFL quality, chain length, aromatics; **positive error = unsafe direction**).
* `build_unified.py`: fuel-domain filter (`DOMAIN`: density 400–1300 kg/m³ & 230–470 K; viscosity 0.05–1e5 mPa·s; vapour pressure ≤1000 kPa; …), dropped rows saved to `measurements_out_of_domain.csv` (`--no-domain-filter` to disable); `lfl_source` column.
* `audit.py`: LFL template no longer counts as data; reports measured vs estimated LFL, domain-filter count, baseline report presence.
* Sanity check (user's refit constants + Table 2 LFL): n-decane 47.4 °C, n-tetradecane 102 °C, n-hexadecane 129 °C flash — close to typical literature values (from memory; verify). Dodecane ≈81 °C (its Antoine fit is flagged `check`).

**Next steps, in order**
1. `python build_unified.py` → `python fill_lfl.py` → `python build_unified.py` → `python baselines_flash.py` → read `data/unified/BASELINE_FLASH_REPORT.md` (first real physics-baseline result; not yet seen).
2. Re-run `python webbook_antoine.py` (page-name fix not yet exercised live; last run still showed the site banner) and delete the two fake 0.0 MB DTIC PDFs.
3. AJFTD: check the site's "Showing x to y of N entries"; re-run `ajftd_scraper.py` (it resumes from `debug/`) until all fuels are saved; then `ajftd_parse_html.py` → `ingest_ajftd.py` → `build_unified.py`; inspect `ajftd_metrics_inventory.csv` and `ajftd_fuel_details.csv`.
4. Manual data: RP-1 11 blends, Li et al. flash-point SI (many more pure hydrocarbon flash points), CRC 647 / AV-11-09 tables, pour/cloud-point literature; add measured LFLs (aromatics/cycloalkanes) to `data/raw/lfl/lfl_measured.csv`.
5. Build `baselines_viscosity.py` (Refutas with kinematic ν from ThermoML densities) and a PPBI baseline; then the residual/mixture model (§12.1).

---

## 16. Addendum — state after `audit_bundle_9608b8`, gap analysis, and adding vapour pressure + octane

### 16.1 What this bundle showed (audit: FAIL 0 · MISSING 0 · WARN 8 · OK 121)
**The pipeline is now ahead of this document.** `UNIFIED_REPORT.md` lists sources whose loaders were NOT written by the assistant and were not seen:
`flashpoint_li` (14,696 rows), `fortin2012` (252), `crc` (5,431), `chew2024` (3,582), `saldana11` (631), `pure_freezing` (27), plus
`dcn_pure`, `mon_pure`, and a much richer AJFTD set (77 fuels). `build_unified.py` also gained a flash-point domain rule (−90…250 °C).
**Before editing `build_unified.py` again, read the user's current version — do not overwrite it with the version in an old zip.**

Hydrocarbon-only records after dedupe (71,129 raw → 47,531; 700 conflicts):

| property | pure | binary | ≥3 comp | real fuel (AJFTD) |
|---|---|---|---|---|
| flash point | 789 | 123 | 3 | 74 |
| vapour pressure (kPa, vs T) | 4,415 | 2,021 | 64 | 76 (D2879) + 76 (D6378) |
| viscosity dynamic (mPa·s) | 2,270 | 16,037 | 625 | – |
| viscosity kinematic (cSt) | 67 | 28 | 0 | 266 |
| density | 1,752 | 12,747 | 897 | 511 |
| DCN | 80 | 213 | 133 | – |
| MON | 267 | 30 | 76 (up to 34 comps) | – |
| **RON** | **0** | **0** | **0** | – |
| pour point | 0 | 0 | 0 | 57 (all −80…−60 °C → no signal) |
| freezing point | 25 | 0 | 0 | 78 |
Other AJFTD properties (real fuels only): distillation T0…T100 (~94–108 fuels), aromatics/saturates/olefins/naphthalenes vol %, sulfur classes,
net heat content, hydrogen wt %, surface tension, speed of sound, specific heat, dielectric constant, conductivity, lubricity, smoke point, thermal stability.
**No molecular/GC×GC composition for any real fuel.**

### 16.2 First real physics baseline (flash point; `BASELINE_FLASH_REPORT.md`)
* 915 records; only **147 predictable**: 439 miss Antoine, 329 miss Antoine **and** LFL. Antoine exists for 90 of 789 flash-point compounds; `lfl.csv`
  has 695 compounds (670 estimated, 25 measured) but the registry now has 1,038 → **stale** (rerun `fill_lfl.py` after `build_unified.py`, then `build_unified.py` again).
* Predicted subset: MAE 7.3 °C, RMSE 15.3, median |err| 2.4, 70.8 % within 5 °C, 40 % over-predicted (unsafe), 17 % over by >5 °C.
* **Works**: C6–C8 (MAE 2.1), C9–C12 (2.9), all binary/ternary mixtures (MAE 1.6, n=57 → ideal mixing is adequate here), records with all-measured LFL (MAE 2.2 vs 10.1 with estimated LFL).
* **Fails**: C13–C16 (MAE 14.9, 46 % over by >5 °C), >C16 (MAE 22.3, bias +13.6), pure aromatics (MAE 12.8, bias +8.6). Worst records mostly have `conflict=True`.
* Interpretation: the dominant problems are (a) missing/extrapolated vapour pressure for heavy and unusual compounds, (b) estimated LFL, (c) noisy/conflicting labels
  (`flashpoint_li+leenhouts+saldana11` 124 conflicts, `crc647` 285), not mixing non-ideality.
* Antoine-coverage "regression" (111 → 69 compounds covering 313–373 K): **not a fit-quality loss.** The 42 compounds are all heavy (C21–C42 alkanes, PAHs); the domain filter's
  vapour-pressure lower bound (1e-6 kPa) removed their 298 K points. Heavy compounds need Antoine at their *flash temperature* (≈400–520 K), not at 313–373 K.
  Fix: lower the bound to ~1e-9 kPa (or drop it), and judge coverage by "valid at predicted flash T ± 20 K" (what `baselines_flash.py` already reports).

### 16.3 What is still needed to train the model
**A. Data (ranked)**
1. **Vapour pressure for the ~700 compounds without it** — scrape NIST WebBook (Antoine, Tb, Tc, ω) for every registry compound; for the rest use a structure-based
   estimator (§16.4). This unlocks the flash baseline for 768 records *and* is itself a target.
2. **Realistic multi-component data** (the largest gap): ≥3-component flash = 3 rows; viscosity ternary 625; no pour point. Get: NIST RP-1 11 blends (+ composition),
   AFRL-TR-2020-0017 (jet-fuel properties incl. GC×GC; download `AD1093317` in a normal browser — DTIC blocks scripts), CRC AV-27-18/647 composition tables, any
   detailed-hydrocarbon-analysis (DHA/GC×GC) data tied to measured flash/viscosity/vapour pressure.
3. **Real-fuel representation**: with only class vol % and a distillation curve, real fuels cannot be fed to a molecule-based mixture model. Either obtain composition
   (item 2) or build a pseudo-component representation from D86 + density + class % (and validate flash/VP on the 74/76 AJFTD fuels).
4. **Pour point**: no usable data (AJFTD values are bounds). Sources: Lira thesis (d.lib.msu.edu/etd/3010), Daridon 2002, Coutinho, Joshi/IIP, CRC AV-11-09 (freezing). Needs n-paraffin content as a feature.
5. **Mixture viscosity at low T** (jet spec −20 °C): mixture data is 293–323 K only; ThermoML extends pure/binary to 230–470 K; AJFTD kinematic viscosity (266 rows) is the real-fuel check.
6. **Octane**: RON absent, MON small; see §16.5.
7. **Label governance**: record measurement method (open/closed cup, ASTM), use `conflict`, `n_refs`, `n_outlier` as weights; decide per-property tolerance; handle deuterated/odd compounds in the registry (e.g. `[2H]` pyrene).

**B. Code still to write**: baselines for viscosity (Refutas with ThermoML densities), density (ideal volume), vapour pressure (Raoult), pour (PPBI), octane (linear volume);
a `model_table` builder (features + masks + `system_id` folds + baseline columns); molecule/mixture encoder; trainer with ensembles; evaluation (log-viscosity error, one-sided flash error);
pseudo-component builder for real fuels; later inverse design.
**C. Decisions needed from the user** *(v2: answers and the revised list are in §17.2; the text below is kept for history)*
1. Fuel domain: **jet/diesel/RP-1** (flash, viscosity, pour, DCN) vs **gasoline** (RON/MON/RVP). Octane and RVP are gasoline properties and are dominated by **oxygenates (ethanol, MTBE/ETBE)** — the current `is_hc` (C/H only) filter would drop them. Recommended: element whitelist `{C,H,O}` configurable per property.
2. Meaning of "vapour pressure": (i) pure Psat(T), (ii) mixture bubble/total pressure at (T,x), (iii) Reid/ASTM vapour pressure of a fuel. All three are supportable; pick the primary.
3. Which properties are targets vs inputs (e.g. density as input to kinematic viscosity).
4. Held-out test sets that are never tuned on (RP-1 blends, a slice of AJFTD fuels, whole molecule families).

### 16.4 Vapour pressure as a target
* **Data in hand**: pure 4,415 (T 230–520 K, ≤1000 kPa), binary 2,021, ternary+ 64, AJFTD D2879/D6378 (0.19–101 kPa at 298–339 K, 76 fuels each).
* **Pure Psat(T)**: predict `log10 P = A − B/(T+C)` (or a Clausius–Clapeyron-form head) from molecular features; with only ~268 compounds having data, use physics-informed features (Tb, Tc, ω, group counts),
  augment with WebBook Antoine sets inside their validity ranges, and validate leave-compound-out. Report error in log10 P and, more importantly, in flash-point units.
* **Mixtures**: Raoult baseline `P = Σ x_i Psat_i(T)` + learned excess (activity-coefficient) correction using the 2,021 binary rows.
* **Fuel Reid vapour pressure**: baseline blending index `RVPBI_i = RVP_i^1.25`, `RVP_mix = (Σ v_i·RVPBI_i)^(1/1.25)` (verify before relying on it); no gasoline RVP blend data is in hand yet.
* The same Psat model feeds the flash-point physics baseline (§2).

### 16.5 Octane numbers (RON/MON) as targets
* **In hand**: MON (pure 267, binary 30, multi 76 up to 34 components) from the Leenhouts benchmark; DCN 426 records (diesel side). **No RON.**
* **Leads to acquire (not downloaded/verified)**: Rittig et al. 2022 *Graph ML for design of high-octane fuels* — data in GitLab
  `https://git.rwth-aachen.de/avt-svt/public/graph_ML_fuel_design` (RON, octane sensitivity); Jameel et al. 2018 (KAUST; SI: RON/MON of 128 pure hydrocarbons, 123 hydrocarbon–ethanol blends, 30 FACE gasoline–ethanol blends);
  *AI-driven design of fuel mixtures* (Commun Chem 2022, doi 10.1038/s42004-022-00722-3; pure species from 19 classes, surrogate fuels of 2–10 components, gasolines; RON/MON/YSI — check data availability);
  vom Lehn et al. (Mendeley doi 10.17632/gzh2z38pyh.2, models/features); Spencer UCT thesis (RON from molar composition).
* **Do NOT use** `glezdiazh/BIOFUELS.PTML` `SI03.xlsx` ("5,000 fuel blends"): inspected — random compositions with carbon dioxide, aniline, hexanoic acid and RON values from −25 to 140: synthetic.
* **Physics baseline**: volume-fraction linear blending `ON_mix = Σ v_i ON_i` (needs densities for mole/mass→volume); learn the non-linear blending residual (ethanol, aromatics, synergy/antagonism);
  model sensitivity `S = RON − MON` jointly. Octane data are volume-based; MON/RON are measured at fixed engine conditions (no T).
* New loader needs (same schema as §6): `property ∈ {ron, mon}`, `fraction` basis = volume (convert to mole/mass with densities), `is_hc` relaxed per §16.3-C1.

### 16.6 Recommended order from here
1. `build_unified.py` → `refit_antoine.py` → `fill_lfl.py` → `build_unified.py` → `baselines_flash.py`; relax the VP lower domain bound.
2. Download the two DTIC PDFs in a browser; rerun `webbook_antoine.py` (page-name fix unverified live) and extend it to all registry compounds.
3. Decide §16.3-C; then build the vapour-pressure model (it unblocks flash), then octane loaders + baselines, then the `model_table` builder and first residual model.


---

## 17. Addendum (2026-10-29) — scope v2: seven targets, new data plan, feasibility gates

### 17.1 What changed
* Targets are now **kinematic viscosity, liquid density, flash point, pour point, GCV, octane (RON/MON) and smoke point**
  (smoke point added at the user's request on 2026-10-29). Vapour pressure, DCN, freezing/cloud point and distillation stay as
  auxiliary properties (vapour pressure because the flash baseline and COSTALD density depend on it).
* The user's rule: **what is finally modelled depends on what data turns out to be available.** Section 17.4 turns that into
  explicit gates so the decision is mechanical rather than argued per property.
* New data work is itemised as tasks T00–T42 in `DATA_ACQUISITION_TASKS.md`. That file is the ledger; this section is the plan.
* New canonical properties and columns (`tier`, `method`, `derived`, `domain`, `source_table`) are defined in that file §0.

### 17.2 Decisions (replaces §16.3-C; confirm with the user if any is wrong)
1. **Two fuel domains**, stored as a `domain` column: `middle_distillate` (jet, diesel, RP-1) and `gasoline`.
   Octane is trained/evaluated only on `gasoline` (+ pure/blend hydrocarbon and oxygenate rows); pour point, flash point,
   smoke point and viscosity are reported mainly on `middle_distillate` (smoke point is a kerosene/jet specification test). Density and GCV span both. Use per-property loss masks, not one model that is
   forced to predict every property for every row.
2. **Element whitelist is per property**: `{C,H}` by default for middle distillates, `{C,H,O}` for octane and GCV rows involving oxygenates (ethanol, MTBE/ETBE).
   The old `is_hc` flag stays; add `elements_ok` per property.
3. **Vapour pressure** means primarily pure-component Psat(T) (it feeds the flash baseline). Mixture Raoult pressure and Reid vapour pressure are later extensions.
4. **Inputs vs targets**: density and vapour pressure are inputs to other baselines AND targets; train them first and freeze them
   (or use measured values) when computing the baselines of the other properties, to avoid leakage.
5. **Held-out sets that are never tuned on**: the 11 NIST RP-1 blends (once composition is obtained, T33/T37), a random 20 % slice of AJFTD fuels,
   all measurements of selected whole molecule families (e.g. one aromatic subclass), and (for octane) the FACE/real-gasoline records if obtained (T18).

### 17.3 Per-target plan

| target | physics baseline | data now (§16.1, dedup, hydrocarbon-only) | main gap | new sources (task ids) |
|---|---|---|---|---|
| viscosity (kin., cSt) | Refutas (§2); needs kinematic ν = η/ρ at the same T | dynamic 2,270 pure / 16,037 binary / 625 ≥3-comp; kinematic 67 pure, 28 binary, 266 AJFTD | low-T mixture data (mixtures 293–323 K only); real-fuel composition | ThermoML (done), NIST RP-1 blends 263–373 K (T33), AJFTD |
| density | ideal volume additivity `ρ_mix = 1/Σ(w_i/ρ_i)` at the same T; COSTALD optional (needs Tc, Pc, SRK-ω, V*) | 1,752 pure / 12,747 binary / 897 ≥3-comp / 511 AJFTD | density vs T for real fuels (AJFTD only 15.56 °C) | NIST RP-1 / Jet-A / JP-8 / S-8 papers (T33, T34), critical constants (T39, T40) |
| flash point | Le Chatelier + Raoult + Antoine + LFL (§2); baseline MAE 7.3 °C on the 147 predictable records (§16.2) | 789 pure / 123 binary / 3 ≥3-comp / 74 AJFTD | Antoine for heavy compounds, measured LFL, cup-type labels | T01–T05, T07–T10 |
| pour point | PPBI (§2); poor for waxy blends | **none usable** (AJFTD values are bounds) | any curated data; needs n-paraffin content as a feature | T28–T32 (literature transcription; expect low volume) |
| GCV | `GCV_mix = Σ w_i GCV_i` (mass-weighted; heats of combustion add for ideal liquids) | none in unified yet | pure-compound ΔcH (liquid), real-fuel GCV from net heat content + hydrogen wt% | NIST WebBook ΔcH°liquid (T04, T24), NBS/API-44 cross-checks (T25), chemicals combustion data (T26), AJFTD (T27) |
| octane (RON, MON) | volume-linear `ON_mix = Σ v_i ON_i` | MON 267 pure / 30 binary / 76 multi; **RON 0** | RON entirely; blends of known composition; oxygenates | vom Lehn 2021 (T13), Rittig 2022 (T14), Jameel 2018 (T15), Sarathy blends (T16), API-45 (T17) |
| smoke point (mm) | threshold sooting index (TSI): blend component TSIs (mixing basis as used in the TSI papers — **verify mole vs mass vs volume basis in Mensch 2010 / Boehm 2022 before coding**), convert to mm with the lamp constants given in the source paper | none in unified; AJFTD smoke point for 82 real fuels (no composition); literature: TSI for ~65 molecules / 124 points (Boehm 2022 training set), YSI ≥400 compounds (YSI is a different scale; correlation to smoke point is tier C) | pure-compound and blend smoke points on one scale; composition of the real fuels | Boehm 2022 (T19), Mensch 2010 (T20), YSI (T21), blends/jet fuels with composition (T41), AJFTD (T23) |

GCV definitions (state them in every GCV record): **gross** = water formed ends as liquid; **net** = water as vapour. Reference temperature
25 °C unless a source says otherwise (API-44 used 25 °C). Conversion used for real fuels: `GCV ≈ NCV + 0.218 × H_wt%` MJ/kg
(8.94 kg water per kg H, latent heat ≈ 2.44 MJ/kg); re-derive if a source uses a different reference temperature.
**Honest expectation for GCV:** because the baseline is nearly exact for known compositions, a learned residual may add nothing.
Report baseline error first; the learning problem is predicting GCV for real fuels from class/carbon-number/hydrogen-content features.

### 17.4 Feasibility gates (decide per property after the data tasks are done)
Compute, from `measurements_dedup.csv` (tier A/B only, `conflict=False`, one count per `system_id`):
`n_multi` = independent systems with ≥3 components; `n_binary_ternary` = unique binary + ternary mixtures; `n_real` = real-fuel/blend records with composition.

| gate | condition (proposed thresholds; adjust with the user) | what is delivered for that property |
|---|---|---|
| **G1** residual model | `n_binary_ternary` ≳ 1,000 (the Leenhouts benchmark found ~10³ unique training mixtures were needed to beat ideal mixing) **and** `n_multi` ≳ 100 **and** `n_real` ≥ 20 for a held-out test | physics baseline + learned residual + calibrated uncertainty, grouped-CV by `system_id` |
| **G2** calibrated baseline | `n_binary_ternary` 100–1,000 or `n_real` 5–19 | physics baseline with a low-capacity correction (few parameters per compound class), conformal intervals |
| **G3** baseline only | below G2 | physics baseline, reported with its error on whatever real data exists; **do not** train a network on synthetic labels |

Preliminary placement from §16.1 numbers (to be recomputed after the tasks):

| property | likely gate today | what would move it |
|---|---|---|
| viscosity | G1 for dynamic binary; G2/G3 for real-fuel kinematic until composition is obtained | T33, T36, T37 |
| density | G1 on binary/multi; real-fuel T-dependence limited | T33, T34 |
| flash point | G3 (3 ≥3-comp systems; 147 predictable records) | T01–T05, T07 (baseline coverage), T36 (real fuels) |
| pour point | below G3 (no data) | T28–T32, honest outcome may stay G3 |
| GCV | baseline only (additive), G3 by definition until real-fuel GCV with composition exists | T24, T27, T36 |
| octane | G3 (RON 0; MON multi 76) | T13–T16 |
| smoke point | G3 today (82 AJFTD fuels without composition; few blend points); G2 is realistic if T36/T41 give composition | T19–T21, T23, T36, T41 |

If a property ends at G3, the deliverable for it is the baseline plus a clear statement of why no model was trained.
**Never** replace missing real data with synthetic blends labelled by the blending rule (§1, §13).

### 17.5 What else is needed beyond §16 (checklist)
1. **Real-fuel composition (GC×GC / DHA)** — affects every property. Sources: AFRL-TR-2020-0017 (T36), NIST RP-1 blend SI (T33/T37), CRC reports (T38). Fallback: pseudo-component representation (§16.3-A3).
2. **Cup type and test method per measurement** (flash: open/closed, ASTM D56/D93/D3828; octane: D2699/D2700; pour: D97/D5853; GCV: D240/D4809/D3338). Without it, conflict flags mix systematic method differences with noise.
3. **Measurement tolerances per property** for the `conflict` flag. Existing: flash 2 °C, pour/cloud/DCN/MON 3 °C, freezing 1 °C, others 5 % relative. Proposed (verify against the ASTM precision statements before adopting): octane ±1 ON, GCV ±0.5 %, density ±0.3 %, smoke point to be taken from the D1322 precision statement (not yet looked up; compare in TSI space as well, since a small mm error at high smoke point is a larger TSI error).
4. **Critical constants and COSTALD inputs** (Tc, Pc, SRK-tuned ω, V*) only if the COSTALD route from the original method document is kept (T39, T40). The SRK ω must be fitted to Psat points, not copied from Pitzer tables.
5. **LFL with measurement temperature**, measured where possible (T07, T10); estimated LFL gave MAE 10.1 °C vs 2.2 °C with measured values (§16.2).
6. **n-Paraffin content / carbon-number distribution** as a model feature for pour point (and freezing/cloud point).
7. **Real gasoline set** (FACE/AVFL with PIONA and RON/MON) as the octane hold-out (T18).
8. **Licensed/paywalled material** the user must fetch: ASTM STP 225 (API-45), Yaws, CRC Handbook, DIPPR (not used), DTIC reports, journal SI behind paywalls.
9. A `domain` and `tier` column in the unified table and in the model-table builder.
10. **Smoke-point specifics.** (a) Store TSI as the primary quantity and keep mm as a derived view; (b) keep the lamp/method (`D1322`, instrument) and the TSI anchors used by each source, because different campaigns are on different scales until stitched (methylcyclohexane = 5, 1-methylnaphthalene = 100); (c) YSI values are not smoke points — only use the YSI→smoke-point correlation as tier C; (d) fuels with smoke point above the lamp range are censored (">50 mm"-type values): store as bounds with a `qualifier`, as for AJFTD pour points; (e) the real-fuel smoke points need class/aromatic/naphthalene content as features (AJFTD has aromatics and naphthalenes vol %).
11. **Verify-before-use step.** The harness must first run task T00 (check every URL, function name and API heading in `DATA_ACQUISITION_TASKS.md` and write `SOURCE_CHECK.md`) and report any source that is blocked, moved or different from what this document says, before any loader is written.

### 17.6 Source register (new in v2; status after this addendum: all TODO)
Full detail, URLs and acceptance tests are in `DATA_ACQUISITION_TASKS.md`. Verified to exist by web search on 2026-10-29 (contents of repositories and SI files **not yet inspected**):

| area | source | note |
|---|---|---|
| Antoine / Psat | Růžička & Majer 1994 JPCRD 23,1 (free NIST SRD PDF); Chickos & Lipkind 2008 JCED 53, 2432; `chemicals` and `thermo` data tables; NIST WebBook | Yaws Handbook is paywalled |
| flash / LFL | `chemicals.safety` (IEC 60079-20-1 (2010) preferred, NFPA 497 (2008) secondary); Li et al. 2020 (already loaded); PubChem experimental properties | Cheméo not used (mixes calculated values) |
| octane | vom Lehn et al. 2021 (615 molecules; some values ANN-estimated → filter); Rittig et al. 2022 GitLab; Jameel et al. 2018 SI (128 pure HC, 123 HC–ethanol, 30 FACE); Sarathy blend papers; API-45 / ASTM STP 225 (paywalled) | exclude `BIOFUELS.PTML SI03.xlsx` (synthetic) |
| smoke point | Boehm et al. 2022 (TSI of 65 molecules / 124 points); Mensch et al. 2010; YSI database (Harvard Dataverse); AJFTD | YSI-derived smoke points are tier C |
| GCV | NIST WebBook ΔcH°liquid; NBS/API-44 papers (free NIST PDFs); `chemicals.combustion` (~1,400 chemicals, tier C); AJFTD net heat content | gross vs net basis must be recorded |
| pour point | Krishna et al. 1989 (24 diesel blends); Coutinho & Daridon 2002; NASA NTRS 19860018789; Kirgina 2019; CRC AV-11-09 | expect low volume |
| density | ThermoML (done); NIST RP-1 (11 blends, 278–343 K), Jet-A, JP-8, S-8 papers | composition of RP-1 blends still needed |

Not verified (from background knowledge, treat as leads): PubChem PUG-View heading names, GESTIS terms of use, KDB/CHERIC, Stull 1947, Morgan & Kobayashi 1994, Calcote & Manos 1983 TSI formula.

### 17.7 Next steps, in order (replaces §12/§16.6 ordering for new work)
0. Run T00 (source check) and report. Then confirm the §17.2 decisions with the user (domains, element whitelist, held-out sets).
2. Run T03, T04, T07 → `build_unified.py` → `refit_antoine.py` (T05) → `build_unified.py` → `baselines_flash.py`; compare coverage with §16.2 (147/915 predictable).
3. T24/T27 (GCV) and `baselines_gcv.py`; T13–T15 (octane) and `baselines_octane.py`; T19–T21 and `baselines_smoke.py` (TSI blending); then `baselines_viscosity.py`, `baselines_density.py`, PPBI baseline.
4. T41 (blends and jet fuels with composition and sooting data); T33/T36/T37 real-fuel composition and validation sets.
5. Pour-point literature tasks T28–T32.
6. Recompute §17.4 gate table in `UNIFIED_REPORT.md`; decide per property: G1 model, G2 calibrated baseline, or G3 baseline-only.
7. Only then build `model_table`, the mixture encoder and the trainer (§12.1), for G1 properties only.
