# PROJECT HANDOFF: UGRP Hydrocarbons (fuel-blend property prediction)

State as of 2026-11-23 (v10 closed; paper-writing phase). Written so that a person or an agent who has never seen this project can pick it up without asking anyone.

**Read order for a newcomer:** this file (all of it), then `docs/README.md` (index of current and historical docs), `docs/GUIDE.md` (teaches every concept from zero), `docs/RESULTS_FOR_PAPER.md` (every paper number, generated), then `README.md`.

**What this file replaces.** The v4-era handoff of the same name and the v1-era handoff in `AGENTS.md` are superseded. Where older docs conflict with this file, `docs/RESULTS_FOR_PAPER.md` and the `docs/final_numbers_v*.md` files win.

---

## 1. One-paragraph summary

The project predicts six properties of a hydrocarbon fuel blend (kinematic viscosity, flash point, density, freezing point, net heat of combustion, smoke point) from its composition, then uses the predictions for inverse design (search for a blend that meets a Jet A or Jet A-1 specification). The backbone is classical physics (Refutas viscosity blending, ideal volume density, Le Chatelier + Raoult + Antoine flash point, Schroeder-van Laar + UNIFAC freezing). Two small LightGBM pair-correction models (viscosity and density excess volume) passed pre-declared gates on locked held-out systems and are used inside a stated domain. Every number carries a one-sided 90 % conformal margin computed from held-out lab systems. Where a model did not beat a trivial constant on real fuels (freezing, heat of combustion on reconstructed fuels, cold viscosity, smoke on the broad set) it is labelled "no skill shown". Of 16 designs (8 pools x 2 specs), 3 pass the certification test **under the model only**. No laboratory measurement of any design exists and none is possible (user decision 2026-11-14).

---

## 2. Who, what, and the hard constraints

| Item | Value |
|---|---|
| Project | UGRP (undergraduate research), Ahmedabad University |
| Authors, in order | Kahaan Shah (first), Namyaa Parmar (second); faculty: Arijit Ganguli, A. S. K. Sinha |
| Presentation | faculty presentation 2026-11-17; a conference deck (ChemCon-style, 16 slides) also exists |
| Repo root | `C:\Users\kahaa\OneDrive\Desktop\UGRP_hydrocarbons` |
| Working branch | `feat/external-datasets-ingest` (holds v5 to v10) |
| Full pre-simplification state | git tag `v5-full-2026-11-16` |
| Other branches / tags | branches `main`, `v4-integrity`, `v4-review-fixes`, `v5`; tags `v2.1-2026-11-05`, `v3-2026-11-06`, `v4`, `v4-r3-final`, `v5-full-2026-11-16` |
| Python | 3.12.5 in the repo `.venv` (Windows). The deck/PDF tooling uses a separate `py -3.11` |
| Source of truth for numbers | `docs/CURRENT_NUMBERS.md` (latest values changed since v5) and `docs/RESULTS_FOR_PAPER.md` (built from `docs/final_numbers_v5.md`, `_v7`, `_v8`, `_v9`, `_v10`) |

**Rules the user has set (keep following them):**

1. **Do not change numbers** to make a story nicer. Every figure must trace to `docs/final_numbers_v5.md` or a committed results file.
2. **Keep the qualifiers next to the claim**: no laboratory validation; the three certified designs are model-level only; the 70-fuel flash result is descriptive (those fuels were used to build the route); freezing, heat and cold viscosity show no real-fuel skill.
3. **Keep negative results.** Failed models stay in code and in `docs/v5_negative_results.md`. Never delete them.
4. **No lab test.** Certification stays model-level and must say so.
5. **Wording.** Use "calibrated", "scope", "validated range", "status: no skill shown". Do not write "honest" or "we say so" in presentation or paper text. Prose should be direct and academic, no em dashes in decks or manuscript.
6. **No attribution lines** in git commit messages or pull request descriptions.
7. **Push only when the user asks** (pushing this branch after the v10 and paper-prep commits was approved on 2026-11-23).
8. **Treat personal and project details as private** (do not paste them into external services).
9. **The manuscript is active** (`manuscript/main.tex`, `references.bib`, `figures/`). Verify every `references.bib` entry before citing.

---

## 3. Repository map (what each top-level thing is)

```
src/fuelblend/        the Python package (all science code)
  api/                predict_route.py (predict_fuel, the public function), predict_blend.py (predict, v4-style wrapper)
  baselines/          the physics: flash.py, flash_unifac.py, cold_flow.py, heat_of_combustion*.py,
                      sooting.py, distillation*.py, cetane.py, octane.py
  mixture/            common.py (ideal-volume density, Refutas, pure-property lookup, registry) + earlier mixture-ML experiments
  eval/               composition_route.py (the designated scoring route), cellkit.py, group_stats.py, metrics.py, cold_offset.py, ...
  ml/                 viscosity_pair.py, density_pair.py (kept, passed gates); set_models.py, ysi_model.py (failed, kept)
  uq/                 conformal.py (margins)
  design/             engine_v4.py (optimiser + judge), library_route.py, blendstock_v5.py, smoke_margin_v5.py
  data/               one loader per source, ThermoML extractor, provenance, CRC composition reconstruction
  features/ splits/ thermo/ models/ train/   descriptors, split generators, thermo helpers, older model code
scripts/v5/           v5 experiment runners; scripts/v5/base/ rebuilds the frozen scoreboard (reproduce_base.py)
                      and holds one delta_*.py per experiment; scripts/v5/final/ builds docs/final_numbers_v5.md
scripts/v4/           earlier evidence pipeline that v5 builds on (kept because src/ and tests import it)
tests/                pytest suites (test_base_* guard the frozen scoreboard); tests/data = fixtures
configs/              uncertainty_v5.yaml (every applied margin), specs_astm.yaml / specs_bis.yaml / specs_jet.yaml, ...
splits/               locked.json and other split definitions (which systems are held out)
data/                 raw/, interim/, processed/, unified/, real_fuels/, cache/pure_grid.npz (frozen), quarantine/
                      NOTE: most of data/ is gitignored (heavy, rebuildable, or licensed); see data/MISSING.md
models/v5/            trained ml_viscosity_pair, ml_density_pair, ml_ysi
results/v5/BASE/      the frozen scoreboard + one delta_* folder per experiment (hash-protected)
results/v5/final/     tables/figures behind docs/final_numbers_v5.md
results/v4, t*, pairwise, eval, hindcast   older results still read by code or tests
docs/                 GUIDE.md, final_numbers_v5.md, claims_register.md, model_decision.md, v5_negative_results.md,
                      decisions_2026-11-15.md, overnight_2026-11-15.md, FULL_REPORT_v4.md, scope_decisions.md, ...
SHA256SUMS_v5_base.txt, SHA256SUMS_v4.txt   hash manifests (detect silent changes to frozen results)
Fuel_Blend_Methodology.pdf                  the methodology document
```

**Deliverables folder (outside the repo):** `C:\Users\kahaa\OneDrive\Desktop\UGRP_deliverables\`

| File | What |
|---|---|
| `UGRP_Faculty_Presentation.pptx` | about 42 slides (about 30 main plus a Backup section with code slides and a flowchart appendix) |
| `UGRP_Conference_Presentation.pptx` | 16 slides |
| `UGRP_Faculty_Talking_Points.pdf` | 15 pages, per-slide SAY / EXPLAIN / IF ASKED / NUMBERS |
| `UGRP_Conference_Talking_Points.pdf` | 6 pages, the speaker notes as sentence bullets |
| `UGRP_Flowchart.pdf/.png`, `UGRP_Flowchart_Conference.pdf/.png` | faculty poster flowchart and the 6-stage conference flowchart |
| `GUIDE.md` | copy of `docs/GUIDE.md` (refresh it after editing the repo copy) |
| `PROJECT_HANDOFF.md` | copy of this file |
| `build_faculty.js`, `build_conference.js`, `lib.js`, `make_flowchart*.py`, `crop_flow.py`, `mkeq.py`, `assets/` | the sources that generate the decks and charts (`*_v1` / `_prefit` files are backups) |
| `manuscript/` | only `refs.bib` (entries marked [VERIFY] are NOT verified) and `figs/` |

---

## 4. How to set up and run (Windows, PowerShell)

```powershell
cd C:\Users\kahaa\OneDrive\Desktop\UGRP_hydrocarbons
.venv\Scripts\python.exe -m pip install -r requirements.txt     # only if the venv is new

# Rebuild the frozen scoreboard, run the BASE tests, verify every file hash (a few minutes; exit code 1 means drift)
.venv\Scripts\python.exe scripts/v5/base/reproduce_base.py

# Fast guard of the frozen scoreboard (47 passed, 1 skipped, about 2.5 minutes, run 2026-11-16)
.venv\Scripts\python.exe -m pytest tests/test_base_*.py -q

# Everything (long)
.venv\Scripts\python.exe -m pytest -q

# Predict one blend
$env:PYTHONPATH="src"
.venv\Scripts\python.exe -c "from fuelblend.api.predict_route import predict_fuel; r=predict_fuel([('CCCCCCCCCCCC',.45),('CCCCCCCCCCCCCC',.20),('Cc1ccccc1',.15),('CC1CCCCC1C',.20)],T_K=288.15,fuel_type='jet'); print(r['flash_point_C']['value'])"

# Inverse design on a pool
.venv\Scripts\python.exe scripts/v5/run_design_userpool_v5.py --spec astm:jet_a,bis:jet_a1 --pool v5pool --n-starts 30 --seed 42 --out results/v5/demo_design
```

`make reproduce` and `make test` do the first and third (`make PYTHON=.venv/Scripts/python.exe test`).

**Environment gotchas**

- The project folder is inside OneDrive. It works, but sync can be slow and deep paths can hit the 260-character limit (`git config core.longpaths true` is set).
- Do not `git checkout` BASE CSV files: the manifest hashes the CRLF working-tree files while `.gitattributes` stores LF, so a checkout creates false drift.
- pandas 3: assigning mixed-type lists through `.loc[k, [cols]] = [...]` fails; assign per column.
- Python heredocs with non-ASCII characters break under the Windows cp1252 code page: write a UTF-8 script file or set `PYTHONIOENCODING=utf-8`.
- Git Bash and PowerShell have different syntax; `&&` is not available in Windows PowerShell 5.1.
- AJFTD data is login-gated and not committed. `data/` is mostly gitignored, so a fresh clone cannot rerun data ingestion; it can rerun the frozen scoreboard only if the cached inputs it reads are present (see `data/MISSING.md`).

---

## 5. The science in the order a prediction uses it

Full detail and worked examples: `docs/GUIDE.md` parts A and 3. The designated route is `src/fuelblend/eval/composition_route.py`; the public entry point is `src/fuelblend/api/predict_route.py::predict_fuel`.

| Property | Designated method | Learned part | Applied margin |
|---|---|---|---|
| Density | ideal volume mixing, 1/rho = sum w_i/rho_i | `ml_pair_excess_volume`: a labelled variant for known binaries only (v5 settings failed ternaries, 3.92 against 3.64 kg/m3 for ideal; after v11b tuning ternaries are about equal, 3.64, so still a variant) | +4.59 / -3.25 kg/m3 (ML variant +/-4.42, unchanged to 4 s.f. by the retune) |
| Kinematic viscosity | Refutas (VBN = 14.534 ln(ln(nu+0.8)) + 10.975) | `ml_pair_residual` (LightGBM on pair features) is the designated value **inside its domain** (known molecules, 233 to 393 K, class pairs seen at least 5 times, no stand-ins) | ML: 0.120 ln (283 K and above; v11c hybrid, 0.134 before), 0.608 (263 to 283 K), 0.287 (below 263 K); Refutas 0.186 |
| Flash point | Le Chatelier + Raoult + Antoine (anchored to measured vapour pressure) + Zabetakis LFL correction, solved with `brentq` on 150 to 650 K; registry overlay adds bicyclohexyl and cis-decalin | none (nothing fitted to blends) | 2.207 C (35 lab systems; the yaml value 3.699 is superseded) |
| Freezing point | Schroeder-van Laar with UNIFAC activity coefficients (`svl_unifac`); ideal SvL is a reported variant | none | 13.84 C (variant 16.98 C); lab margin only |
| Net heat of combustion | group additivity, sum w_i dH_c,i | none | none; reported only |
| Smoke point | regression SP = 28.84 - 0.324 A - 0.791 D (A aromatics vol %, D diaromatics vol %) | none (ML YSI failed its gate) | 4.11 mm, in range only |

**Margin rule:** for n held-out systems, margin = k-th smallest worst dangerous-side system error, k = ceil((n+1)(1-alpha)), alpha = 0.10; infinity if k > n (so at least 9 systems are needed). Implemented in `src/fuelblend/uq/conformal.py`; applied values in `configs/uncertainty_v5.yaml`.

**Certified** = finite margin AND blend `in_domain` AND spec still met after the margin. It is a model-level label, not an approval.

---

## 6. What the evidence says (headline numbers; all traceable to `docs/final_numbers_v5.md`)

| Property | Result | Verdict |
|---|---|---|
| Viscosity (lab, ML pair correction) | Refutas 0.082 to 0.060 ln dev (264 systems); 0.054 to 0.031 locked (70 systems); ternary+ 0.098 to 0.042; cold 0.099 to 0.061; gates R-MLV-1..4 all passed | beats Refutas |
| Viscosity (real fuels, Refutas) | 263 K and above beats trivial; below 263 K MAE 0.527 ln vs trivial 0.253 | cold: worse than trivial, margin infinity |
| Density (lab) | MAE 1.35 kg/m3 (about 0.17 %); real fuels T2 beat trivial; CRC reconstructed compositions 14.3 vs 13.2 | CRC: no skill (a composition-reconstruction limit) |
| Flash (real fuels, T2, 70 groups) | MAE 2.196 C vs trivial 5.156, bias +0.06 C | beats trivial, **descriptive** |
| Flash (AJFTD GCxGC, 31 fuels) | bias -0.152, MAE 2.995 | beats trivial after the overlay |
| Flash (D56) | rests on one GC fuel (POSF-10325) | no margin |
| Freezing (lab liquidus, 12 systems, 73 points) | UNIFAC MAE 3.65 C, bias +2.18, margin 13.84; ideal 4.19 / +0.88 / 16.98; constant 15.63 (v11c xylene correction; was 3.93 and 4.47) | lab skill only |
| Freezing (real fuels) | pooled n 50: 8.425 vs trivial 6.413; locked n 8: 9.31 vs 5.176 | no skill / worse than trivial; bias +4 to +6 C (safe side) |
| Smoke | 0.904 mm MAE on the 100-fuel set (beats trivial); broad set 1.679 vs about 1.65 | no skill on broad set |
| Heat of combustion | reconstructed CRC fuels worse than trivial (0.376 vs 0.114); AJFTD bulk models beat trivial | no margin issued |
| Designs | 16 designs: 3 certified under the model (P0y Jet A, P0y-b Jet A, P0y-b Jet A-1) | none measured |

---

## 7. Data: what exists and where

| Source | Use | Where / note |
|---|---|---|
| Leenhouts benchmark (Zenodo 19914760) | lab mixture flash, viscosity, DCN, MON | primary T1 evidence |
| NIST/TDE viscosity (CheMixHub, Bilodeau), ThermoML | viscosity, density, vapour pressure | re-extracted after parser fixes |
| NIST WebBook, Zabetakis Bulletin 627 | Antoine constants, LFL | anchored/refit because ranges miss 40 to 100 C |
| AJFTD | 2,601 fuels; 116 with GCxGC, 113 mappable | `data/real_fuels/v5/`; login-gated source |
| CRC 647 | 57 jet fuels | compositions **reconstructed** (r3: from measured hydrocarbon types, Table G.1), `data/real_fuels/v4/` |
| Co-Optima, NISTIR 6646, Fortin 2012 | octane, RP-1 | research or characterisation only |
| Frozen pure grid | `data/cache/pure_grid.npz`, SHA-256 checked at every load (`GRID_SHA256`) | 1,042 compounds (114 with a full grid) |

**Refused:** `DatasetsFuelProperties.xlsx` (50,000 synthetic rows) and the v3 corpus + optimiser (formula-derived; quarantined). Unit traps already fixed: NIST "log viscosity / cP" is really ln(Pa.s); Leenhouts viscosity is ln(mPa.s); WebBook Antoine is in bar (add 2 for kPa); ThermoML uses Pa.s and m2/s.

---

## 8. Decisions on record (why things are the way they are)

User decisions, applied 2026-11-15 (`docs/decisions_2026-11-15.md`):

| # | Decision | Effect |
|---|---|---|
| 1 | Freezing: adopt UNIFAC-SvL, ideal SvL kept as a variant | margin 13.84 C (CLM-73) |
| 2 | D86: Engler only; Rayleigh is a reported check | label `d86_model_dependent` |
| 3 | Flash-registry overlay on by default, re-score | `FUELBLEND_FLASH_REGISTRY_ADDITIONS=0` restores the old registry; adopted as the frozen BASE (U4, CLM-74) |
| 4 | Bicyclohexyl dH_fus = fusion + transition, 13.94 kJ/mol | over-predicts freezing (safe side) by up to 33.6 C on that system |
| 5 | Coutinho lookup: exact aromatic match | 112 misassigned aromatics fixed (substring match kept as default for v4 callers) |
| 6 | `strict_constraints=True` by default | engine raises `ModelUnavailable` instead of silently dropping a constraint |
| 7 | cis-decalin MW verified (138.25 g/mol) | |
| U6 | no dicycloparaffin bound added | the 37.8 wt % caveat on the certified Jet A-1 stays documented |
| U7 | `design_pool(alt_standin=True)` default | CLM-67 |

Other standing decisions: pour point scoped out (`docs/scope_decisions.md` SD-01: 0 speciated mixtures with pour point); no lab test; WS6 cold-viscosity offsets not adopted; D86 margins adopted only after limits are verified; ML density only a variant; ML YSI and ML DCN not exposed. Full table with reasons: `docs/GUIDE.md` section 10.

---

## 9. Frozen evidence and how to change it safely

- `results/v5/BASE/` is the frozen scoreboard. `SHA256SUMS_v5_base.txt` lists a hash for every file. `reproduce_base.py` recomputes and compares; any difference is reported as drift and exits with code 1.
- To make a deliberate change: create a new `results/v5/BASE/delta_<name>/` folder with a README that **declares the question and the pass rule before any result exists**, run it, commit, and only then (if the user agrees) re-write the manifest with `reproduce_base.py --update-manifest`. The 2026-11-15 re-score was adopted this way (commit `3b75437`); the pre-adoption BASE is recoverable at `3ab033d`.
- Claims are tracked as CLM-nn in `docs/claims_register.md`; update it when a claim changes.
- `SHA256SUMS_v4.txt` pins the v4 artefacts; do not touch lines other than the one deliberately rehashed (`docs/claims_register.md`, done on 2026-11-15).

---

## 10. Known issues, caveats and open items (all disclosed; none hidden)

**Scientific caveats**

1. No laboratory validation of any design; certification is model-level. The three certified designs sit on the margin boundary (minimum slack 0.003 to 0.013): optimiser edge, not robustness. They carry labels `d86_model_dependent`, `freezing_route_dependent`, `overlay_dependent`.
2. The certified Jet A-1 design has 37.8 wt % dicycloparaffin, above the AJFTD GCxGC maximum (30 %).
3. Freezing and cold-viscosity margins come from mostly **binary** lab systems and are applied to 6-component designs: a transfer, not a validation. No multicomponent freezing margin exists (WS2 failed: scorable systems alkane-only 2, alkane+naphthene 1, with-aromatic 0 of 8 needed).
4. D56 flash rests on one GC fuel; thin locked sets (3 and 4 cold viscosity systems; CRC fold_b 12 fuels).
5. T2 real-fuel cells are descriptive (those fuels informed the route), not held-out claims.
6. CRC 647 compositions are reconstructed, not measured.
7. Jet A-1 limits and the new D86 limits come from secondary sources; the user intends to read D7566, D1655 and IS 1571 and upgrade the `verification:` flags in the spec yamls. Until then nothing counts as certified against verified limits.
8. The o-/m-xylene eutectic composition was transposed in `ingest_calorimetry_primary.py` (EUT); corrected in v11c (checked by ideal Schroeder-van Laar consistency, not against the page image). Lab liquidus S_F MAE is now 4.187 ideal / 3.648 UNIFAC (was 4.474 / 3.935); margins unchanged. `docs/final_numbers_v5.md` still shows the old S_F values.
9. AJFTD:325 has a freezing prediction of -144.7 C against a measured -81 C under the designated constants (error -63.7 C); the reason it is in the pooled n of 50 and not 49 is documented in `delta_freeze_t2_unifac`.
10. The engine's freezing constraint label still reads `ideal_sle_mole_fraction` in `designs.csv` even though the value and margin are UNIFAC (label issue only).
11. Data in the repo but unscored: lab freezing of hexadecane mixtures (Oakley 2017) and smoke of 23 coal/JP-8 fuels (Yang 2006). The ML density model has not been tested on real fuels. Diesel is unvalidated.
12. Flash: the "wrong routing" bias of -21.7 C in the old `api_as_shipped` path was fixed by moving `predict()` onto the BASE route; the cause of the earlier CRC reconstruction over-prediction (+4.2 C raw) was never isolated.

**Engineering state**

13. Tests: last full run 2026-11-23: 909 passed, 2 skipped, 2 failed. `test_committed_splits_match_rebuild` (input hash of the AJFTD long/wide files differs) is a known failure, deliberately left unchanged. `test_fit_ensemble_logger_none_bitwise` is load-sensitive and passes when run alone. `tests/test_base_*.py`: 48 passed; `reproduce_base.py` reports 0 changed, 0 missing, 0 new of 23 hashed files.
14. Git: branch `feat/external-datasets-ingest`; `results/v5/BASE/delta_*` shows as modified in `git status` only because of line endings (`git diff --ignore-cr-at-eol` is empty and the manifest hashes match). Do not stage, checkout or normalise them.
15. `data/unified/` is the v5 production set (76,328 dedup rows). The v9 additions are in `data/unified_v9/` (gitignored, rebuild with `scripts/v9/build_unified_v9.py`); v5 and earlier results are not retrained on them.
16. The simplified-branch history is kept: the full pre-simplification tree is at tag `v5-full-2026-11-16`.

**Deliverable caveats**

17. `manuscript/refs.bib` contains [VERIFY] entries that were not checked against the real publications; do not cite them until verified.
18. The decks, talking points and flowcharts were checked against `final_numbers_v5.md` and the guide on 2026-11-16; the unverified items list given with the delivery should be re-read before presenting.

---

## 11. Open items (none started)

1. Paper: write from `docs/RESULTS_FOR_PAPER.md`; verify `manuscript/references.bib`; qualifiers stay beside every claim.
2. Spec limits: verify against ASTM D1655, D7566, IS 1571, Def Stan 91-091 and upgrade the `verification:` flags in `configs/specs_*.yaml`.
3. Data gap: real-fuel freezing and cold viscosity need verified fuel-level data (Edwards 2017, AIAA 2017-0146, and the true AFRL-TR-2020-0017 / DTIC AD1093317 were not obtained; the supplied AFRL file is a different paper). Multicomponent flash, freezing and viscosity margins need 8 or more systems per class.
4. v7 fresh flash lock: not scorable with the v5 route (all 61 systems lack Antoine or LFL constants; HANNA is binary only). Needs constants for the lock compounds first.
5. Confirm the xylene eutectic correction against the Pitzer and Scott 1943 page image; supply verified Antoine and LFL for the 34 oxygenates that block the v7 flash lock (`results/v11/lock_missing_constants.csv`).

## 12. Troubleshooting

| Symptom | Likely cause and fix |
|---|---|
| `reproduce_base.py` exits 1 with drift | a BASE file changed (often a CRLF/LF issue after a git checkout, or a real code change). Compare the reported file with the manifest; if the change is deliberate follow section 9 |
| `predict_fuel` raises about the grid hash | `data/cache/pure_grid.npz` differs from `GRID_SHA256`; restore the frozen file, do not regenerate silently |
| `ModelUnavailable` from the design engine | `strict_constraints=True`: a constraint has no model (for example D86 limits or margins missing). Provide the model/margin or pass an explicit opt-out for reproduction runs only |
| A property comes back with `in_domain: False` | a molecule has no own data (stand-in), temperature outside 233 to 393 K, or a class pair seen fewer than 5 times; the main value falls back to Refutas with a warning and cannot be certified |
| Margin is `inf` | fewer than 9 independent held-out systems for that property/bin; this is intentional |
| Pytest failures in older suites | the v4 baseline had 42 failing missing-deliverable asserts (`results/v4/baseline_pytest.txt`), later resolved; re-check on the current branch before assuming a regression |
| `.pptx` will not open for writing | stale PowerPoint process (item 15 above) |
| Non-ASCII characters crash a Python one-liner | set `PYTHONIOENCODING=utf-8` or use a script file |

---

## 13. Rebuilding the presentation deliverables

The decks use `pptxgenjs` (Node) with the helper `lib.js`. From the deliverables folder:

```powershell
$env:NODE_PATH="C:\Users\kahaa\AppData\Local\Temp\sp\npm\node_modules"   # where pptxgenjs was installed; reinstall with npm if missing
node build_faculty.js
node build_conference.js
py -3.11 make_flowchart.py        # faculty poster flowchart (GROW=2.3, FS_CAP=30)
py -3.11 make_flowchart_conf.py   # conference flowchart
py -3.11 crop_flow.py             # cuts the poster into bands A..E used on slides 41 and 42
```

Validate a deck with the pptx skill's `scripts/office/validate.py`, and render slides through PowerPoint (COM export) to check layout. The talking-point PDFs were generated with reportlab (DejaVu fonts) from the decks' speaker notes plus the expanded SAY/EXPLAIN/IF ASKED text; the generator scripts were deleted after use, so regenerate by reading the notes from the .pptx XML if the decks change.

Deck rules carried over from the brief: authors Kahaan Shah then Namyaa Parmar; faculty names on the faculty title slide; qualifiers stay on the same slide as the claim; vocabulary as in section 2.

---

## 14. Glossary pointer

Every term (Refutas, VBN, Antoine, LFL, SvL, UNIFAC, conformal margin, locked set, gate, verdict, tiers T1/T2/T3, CRC647_r3, AJFTD, GCxGC, BASE, delta, CLM) is defined in plain words in `docs/GUIDE.md` Part A and section 2.

---

## 15. Since v5 (v6 to v10), one line each

| Round | Outcome | Source |
|---|---|---|
| v6 | freezing holdout scored once; pure-compound GNN and mixture retries; failures kept | `docs/v6_summary.md`, `docs/v6_negative_results.md` |
| v7 | pure-compound models (random forest or D-MPNN per property) deployed in `predict_pure`; no blend route replaced the v5 designated route; fresh flash lock not scorable | `docs/final_numbers_v7.md`, `docs/methodology_v7.md` |
| v8 | no ML-constant substitution passed G1; real-fuel freezing: status: no skill shown | `docs/final_numbers_v8.md` |
| v9 | 539 new dedup groups (ThermoML, Boehm 2022) built in `data/unified_v9/`; 14 new systems, 4 locked; no freezing fuel added; Antoine refits mostly failed the flash-window gate | `docs/final_numbers_v9.md` |
| v10 | 3 PDFs supplied, 25 perhydrophenanthrene vapour-pressure points verified against page images, 0 new fuel rows; closed as a negative result | `docs/final_numbers_v10.md` |
| v11c | minor optimisations: xylene eutectic fixed (S_F MAE 3.935 to 3.648 UNIFAC); viscosity hybrid adopted post hoc (tuned warm model for T >= 283.15 K, v5 below; warm margin 0.134 to 0.120 ln; designs unaffected, cold only); seed averaging, pure-compound ensembles, tuned random forest and lock constants: not adopted or not done | `docs/final_numbers_v11c.md` |
| v11, v11b | nested-CV tuning of the viscosity and density pair models. Viscosity: not adopted (cold-bin margin guard failed). Density: adopted (loss l1, 296 trees, lr 0.03, 11 leaves); dev ML MAE 0.958 to 0.908 kg/m3, locked second look 0.897 to 0.839; still a labelled variant, ternaries about equal to ideal | `docs/final_numbers_v11.md`, `docs/v11b_predeclaration.md` |

## 16. Handoff checklist

- [ ] Read this file, then `docs/GUIDE.md` Part A and sections 0, 6, 8.
- [ ] `git status` and `git log --oneline -5` on branch `feat/external-datasets-ingest`; confirm tag `v5-full-2026-11-16` exists.
- [ ] `.venv\Scripts\python.exe -m pytest tests/test_base_*.py -q` shows 47 passed, 1 skipped.
- [ ] `scripts/v5/base/reproduce_base.py` reports no drift.
- [ ] Run the full test suite before editing code (item 13); expect one known failure.
- [ ] Re-read the rules in section 2 (no number changes, qualifiers stay, negative results stay, no attribution lines, push only when asked).
