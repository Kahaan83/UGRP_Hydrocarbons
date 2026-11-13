# v5 Certification Path and ML Roadmap

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking. Subagents and reviewers: Opus at most, never Fable (user rule).

**Goal:** Get inverse designs to certify under fixed, pre-declared gates (or report each one as blocked with a named reason), then improve the models with ML where the data supports it and fill the pure-component data gaps.

**Architecture:** Workstream 1 (WS1) clears the smoke-point blocker, the one thing stopping all 8 current designs. It does this by keeping recipes in the model's domain (aromatics minimum), adding the D1655 alternative limit (18 mm with naphthalenes ≤ 3 vol %), shrinking the smoke margin with real-fuel data, and letting the optimiser see recipe-dependent margins. WS2–WS7 are a roadmap. Each one becomes its own detailed plan when it starts, and each one follows the delta pattern: README with declared readings first, results appended later, user decisions quoted.

**Tech Stack:** Python 3.12, pandas, numpy, scipy, LightGBM (existing pair models), PyTorch + PyG (new, WS4 only), RDKit, pytest. Venv `.venv`.

**Spec:** This document; evidence it rests on: `results/v5/BASE/delta_design_ml_visc/README.md`, `results/v5/BASE/delta_smoke_margin/README.md`, `results/v5/BASE/SCOREBOARD.md`, `docs/model_decision.md` ("v5 ML components"), `configs/specs_astm.yaml`.

## Global Constraints

- Success means "certifies under fixed declared gates, or is blocked with the reason". It never means "certifies at any cost".
- Every margin change, spec-verification flag change, or new design constraint is either a declared reading (README written before any result) or a user decision quoted verbatim. Gates are never moved after results are seen.
- An unverified spec limit is never certified. Only the user can upgrade `verification` in `configs/specs_*.yaml`, after reading the primary standard text.
- Train and evaluate only from de-duplicated real measurements. Synthetic or rule-labelled rows never serve as ground truth (AGENTS.md §1, §13). Group splits use `system_id` (fuels use `fuel:<id>`).
- Flash: over-prediction is the unsafe side. Cold viscosity: under-prediction is the unsafe side. Every margin states its direction.
- `engine_v4.py` changes need a test that shows v4/v5 results are unchanged when the new feature is off.
- Do not commit or push unless the user asks (2026-11-12 standing approval covers pushing v5 after commits the user asked for).

## Where things stand (2026-11-13, commit 399a454)

| Constraint (Jet A / Jet A-1, 8 designs) | Status | Why |
|---|---|---|
| density, flash, viscosity −20 °C, NHOC | pass | margins 4.59 kg/m³, 2.21 °C (set by one lab system), 0.287 ln (ML), 0.194 MJ/kg |
| freezing | pass, slack 0.02 | margin 16.98 °C measured on binaries, transferred to multicomponent recipes |
| aromatics | by construction, 0 vol % | the optimiser drives aromatics to 0 |
| **smoke point** | **uncertified in all 8** | 0 vol % aromatics is outside the bulk-regression domain (8–23 vol %, trusted 19–25 mm); no margin in the engine; spec limit `secondary_sources_only` / `UNVERIFIED` |

The smoke arithmetic, which sets what WS1 must achieve:
- The bulk regression is trusted only for SP 19–25 mm, so the 25 mm primary limit can never pass with any margin above 0.
- The only route through this model is the D1655 alternative: SP ≥ 18 mm **and** naphthalenes ≤ 3.0 vol %. With the current 7.8 mm margin, that needs SP ≥ 25.8 mm, which is outside the trusted range, so it still fails.
- **Certification needs, together:**
  - (a) the alternative limit implemented (`engine_v4` ignores `alt_min` / `alt_requires` today; grep finds no use);
  - (b) aromatics ≥ 8 vol % so recipes are in domain (D7566 Table 1 minimum for synthetic blends; secondary source D7566-12a, verify the current edition);
  - (c) an in-domain margin small enough that predicted SP − margin ≥ 18 mm. `smoke_point_from_bulk` gives 24.97 / 24.33 / 23.68 / 22.71 / 21.09 mm at 8 / 10 / 12 / 15 / 20 vol % aromatics, so the largest margin that can certify is **6.97 / 6.33 / 5.68 / 4.71 / 3.09 mm**. The current 7.8 mm misses even at 8 vol %;
  - (d) the user verifying the smoke/aromatics limits from the primary standard text.
- Freezing is the next blocker: forcing aromatics in changes recipes, and freezing already failed once when recipes moved (delta_smoke_margin R-SM-4).

## Review Focus

1. A recipe with aromatics exactly at 8.0 vol % (domain edge). The smoke model must not warn or flip status because of float noise at the boundary. Test it in Task 1.3.
2. A recipe whose naphthalenes are unknown, because the pool has a component with no class or ring data. The alternative limit must be `uncertified`, never `pass`. Test it in Task 1.2.
3. Pools with no aromatic component at all (e.g. `ext_n8` variants). The aromatics minimum makes them infeasible. The design must report `blocked: aromatics_min infeasible for pool` rather than return a 0 % recipe as category B. Test it in Task 1.3.
4. A recipe-dependent margin in the optimiser must give the same final judgement as the post-hoc rejudge for the same recipe. The two paths must never disagree. Test it in Task 1.5.
5. Refitting the smoke margin on real AJFTD fuels must exclude fuels whose smoke value is censored or out of range ("> 25", "No Test"). Censored rows inflate or deflate the margin silently. Test it in Task 1.4.
6. Class-mix transfer. The designs are 68–82 % butylcyclohexane, which petroleum jet fuels never are, yet the bulk model checks only aromatics vol %. A real-fuel margin applies only when the recipe's cycloalkane (and n-/iso-alkane) mass fractions lie inside the range of the AJFTD GCxGC fuels. Outside that range, smoke is `uncertified` with reason `class_mix_outside_real_fuel_range`. Test it in Task 1.4 (`test_class_mix_gate`).

---

## WS1 — Smoke/aromatics path to certification (detailed)

Delta folder: `results/v5/BASE/delta_smoke_cert/`. Declared readings (write them in its README before Task 1.1):
- **R-SC-1:** number of AJFTD fuels with smoke point, aromatics and naphthalenes all present (descriptive).
- **R-SC-2:** smoke margin `alpha=0.10`, direction `over`, refit on in-domain real jet fuels (8–23 vol % aromatics, uncensored smoke).
  - Grouping is chosen now, knowing the levels. `source_group` has only 10 levels on fuels with aromatics (JAA 672, JA 452, JFK 348, LAB_REPORT 271, NUMERIC_ID 190, OTHER 170, PAREN_BATCH 65, HR 49, CRC 2, JETSCREEN 2), so it cannot carry a ≥ 30-group gate.
  - Declared: (i) split conformal over individual fuels (unit `fuel_id`), with leave-one-`source_group`-out fits, so no fuel is scored by a model that saw its group; (ii) the margin is computed per held-out `source_group` too. The **certifying margin = max(fuel-level margin, largest per-group margin over groups with ≥ 40 fuels)**.
  - Gate: at least 200 in-domain fuels in at least 5 groups with ≥ 40 fuels, else the margin stays `inf`.
- **R-SC-2b (user decision, not automatic):** if R-SC-2 passes, the user decides whether the real-fuel margin **replaces** the 7.8 mm lab-mixture margin for in-domain jet recipes. `_combined_margin` takes the max over keys, so adding a key alongside the 7.8 mm one changes nothing. The 7.8 mm margin was the user's own choice ("Yes, add it"), so replacing it needs a new quoted answer. Without that answer the smaller margin is not used.
- **R-SC-3:** the 4 design runs (seed 42, same pools) plus the 2 aromatic pools of Task 1.6, re-run with aromatics_min 8 vol % and the alternative smoke limit. Per design: aromatics vol %, smoke pred, margin, the maximum margin that certifies (pred − 18), smoke status, freezing slack, category. This is a reading, not a unit test (about 1.5 h).
- **R-SC-4:** with every new option off, the new code gives identical constraint results to the pre-change code on the same in-memory weights (the stored CSV weights are rounded to 4 decimals, so CSV comparison is not bit-exact).
- Stated now: certification also needs the user's verification flag. Without it, the best possible outcome is "would certify; blocked only by `secondary_sources_only`".

### Task 1.0: User action — verify spec limits (blocking for certification, not for Tasks 1.1–1.5)

- [ ] Ask the user to read ASTM D1655 (current edition) and D7566 Table 1 for: smoke point 25 mm / 18 mm + naphthalenes 3.0 vol %, and aromatics min 8 vol % (D1319) / 8.4 (D6379), max 25. University library or ASTM Compass access.
- [ ] Also ask for BIS IS 1571 for the Jet A-1 limits.
- [ ] The user edits `verification:` in `configs/specs_astm.yaml` / `specs_bis.yaml` and their decision is quoted in the delta README. The agent never edits this flag.

### Task 1.1: Scrape AJFTD smoke point and naphthalenes for the 2,601 fuels

**Files:**
- Modify: the v5 AJFTD scraper/export (find with `grep -rn "ajftd_fuels_long" scripts/v5`) so it exports `smoke_point_mm` and `naphthalenes_vol_pct`
- Create: `scripts/v5/ingest/ajftd_smoke_naph.py`, output `data/real_fuels/v5/ajftd_smoke_naph.csv`
- Test: `tests/test_ajftd_smoke_naph.py`

**Interfaces:**
- Produces: `load_ajftd_smoke_naph() -> pd.DataFrame` in `src/fuelblend/data/loaders_ajftd_v5.py` (create the file if absent), with columns `fuel_id, source_group, smoke_point_mm, smoke_qualifier, naphthalenes_vol_pct, aromatics_vol_pct, method_smoke`

`data/real_fuels/v5/ajftd_fuels_long.csv` has `aromatics_vol_pct` for 2,290 fuels but **no smoke point rows** (the v1 scrape had 82). Re-parse the saved HTML first (`ajftd_parse_html.py` pattern). Re-scrape only if the pages were not saved; site access is granted, so re-scraping is allowed.

- [ ] **Step 1: Write the failing tests** `test_censored_smoke_kept_as_qualifier` (an HTML fixture value `>25` gives NaN `smoke_point_mm` and `smoke_qualifier == ">"`), `test_no_test_dropped` ("No Test" gives no row), and `test_columns`. Build the fixture from a real saved page structure (AGENTS.md §9 lesson 9).
- [ ] **Step 2:** `pytest tests/test_ajftd_smoke_naph.py -v`. Expected: FAIL (module missing).
- [ ] **Step 3:** Implement the parser and loader.
- [ ] **Step 4:** Tests pass. Run the script and record the counts (R-SC-1) in the README result section.
- [ ] **Step 5:** Commit (only if the user asked for commits this session).

### Task 1.2: Implement the alternative smoke limit (18 mm with naphthalenes ≤ 3.0 vol %)

**Files:**
- Modify: `src/fuelblend/design/engine_v4.py` (`evaluate_constraint`, around lines 267–372), `src/fuelblend/design/library_route.py` (a naphthalenes model)
- Test: `tests/test_v5_smoke_alt_limit.py`

**Interfaces:**
- Produces: `PropertyModel` entry `props["naphthalenes_vol_pct"]` (class volume fraction of naphthalene-type components; `by_construction=True` when every pool component has a known class, else `ModelUnavailable`)
- Produces: `evaluate_constraint` honours `c["alt_min"]` and `c["alt_requires"]`. The result `status == pass` if the primary limit passes, or if `alt_min` passes with margin **and** every `alt_requires` constraint passes. `reason` names the limit that was used (`primary` | `alt_min`).
- Feature flag: `spec["use_alt_limits"]` (default False keeps v4/v5 behaviour)

- [ ] **Step 1: Failing tests:**
  - `test_alt_limit_off_unchanged`: the flag is off, so for the v5pool recipes (weights read from `designs.csv`, held in memory) `evaluate_blend` returns the same statuses, predictions and slacks as a frozen copy of the pre-change results generated by the same in-memory call (fixture written before the change).
  - `test_alt_limit_used`: smoke pred 22.0, margin 3.0, naphthalenes 0.0 gives `pass`, reason contains `alt_min`.
  - `test_alt_limit_naph_unknown`: the naphthalenes model is unavailable, so the result is `uncertified`.
  - `test_alt_limit_naph_over`: naphthalenes 3.5 vol %, pred 22 gives `fail`.
- [ ] **Step 2:** Run. Expected: FAIL.
- [ ] **Step 3:** Implement. The `alt_requires` keys map to constraint keys (`naphthalenes_vol_pct_max` gives an upper limit on `naphthalenes_vol_pct`).
- [ ] **Step 4:** Tests pass. Run `pytest tests/test_v5_design_forward.py tests/test_predict_route.py -q`. Expected: no new failures.
- [ ] **Step 5:** Commit.

### Task 1.3: Aromatics minimum constraint (8 vol %) as a declared design option

**Files:**
- Modify: `configs/specs_astm.yaml` and `configs/specs_bis.yaml`: add `aromatics_vol_pct: {min: 8.0, max: 25.0, source: "D7566 Table 1 (secondary, D7566-12a)", verification: secondary_sources_only, applies_to: synthetic_blend}`, used only when `--synthetic-blend` is passed
- Modify: `scripts/v5/run_design_userpool_v5.py` (add flags `--synthetic-blend`, `--alt-limits`)
- Test: `tests/test_v5_aromatics_min.py`

**Interfaces:**
- Produces: `design_routed(..., synthetic_blend: bool = False, alt_limits: bool = False)`. When a pool has no aromatic member and `synthetic_blend` is set, the result is `category = "blocked_aromatics_min_infeasible"` with no optimisation.

- [ ] **Step 1: Failing tests:**
  - `test_aromatics_edge_no_warning`: a recipe at 8.0 vol % (to 1e-9) raises no `OutOfDomainWarning` for aromatics.
  - `test_pool_without_aromatics_blocked`: a pool of n-alkanes only gives the blocked category.
  - `test_flag_off_identical`: R-SC-4.
- [ ] **Step 2:** Run. Expected: FAIL.
- [ ] **Step 3:** Implement. Clamp the domain check with a tolerance of 1e-6 vol %, documented in `smoke_point_from_bulk`.
- [ ] **Step 4:** Tests pass.
- [ ] **Step 5:** Commit.

### Task 1.4: Refit the smoke margin on real in-domain jet fuels

**Files:**
- Create: `scripts/v5/base/delta_smoke_cert.py`, outputs `results/v5/BASE/delta_smoke_cert/margin_fit.csv` and `readings.json`
- Modify: `configs/uncertainty_v5.yaml` (new key `smoke_point_mm.ajftd_v5_in_domain.bulk_regression`), **only after** R-SC-2 passes its gate **and** the user answers R-SC-2b. The new key then *replaces* `smoke_point_mm.lab_mixtures.bulk_regression` for in-domain jet recipes (in the design smoke model and in `predict_route._smoke`). Adding it alongside would leave the max at 7.8 mm.
- Test: `tests/test_delta_smoke_cert.py`

**Interfaces:**
- Consumes: `load_ajftd_smoke_naph()` (Task 1.1), `smoke_point_from_bulk` (unchanged coefficients if the 100-fuel training set is disjoint from the scoring set; otherwise refit with leave-one-group-out folds and score out-of-fold)
- Produces: `fit_smoke_margin(df: pd.DataFrame, alpha: float = 0.10, group_col: str = "source_group") -> dict` with keys `margin, n_groups, n_rows, direction, coverage_oof`

- [ ] **Step 1: Failing tests:**
  - `test_censored_excluded`: rows with `smoke_qualifier` not empty are dropped and counted.
  - `test_gate_too_few`: 150 fuels, or fewer than 5 groups with ≥ 40 fuels, give `margin == inf`.
  - `test_direction_over`: the margin is the conformal (1−alpha) quantile of (pred − meas) over out-of-fold residuals (leave-one-group-out fits).
  - `test_certifying_margin_is_max`: certifying margin == max(fuel-level, largest eligible per-group margin).
  - `test_class_mix_gate`: a 75 % butylcyclohexane recipe gets `uncertified` with reason `class_mix_outside_real_fuel_range`, even at 10 vol % aromatics.
- [ ] **Step 2:** Run. Expected: FAIL.
- [ ] **Step 3:** Implement.
- [ ] **Step 4:** Tests pass. Run the script and append R-SC-1/R-SC-2 to the README. Ask the user before adding the yaml key (quote the answer).
- [ ] **Step 5:** Commit.

### Task 1.5: Recipe-dependent margins inside the optimiser

**Files:**
- Modify: `src/fuelblend/design/engine_v4.py` (`PropertyModel` gains `margin_fn: Optional[Callable[[list, np.ndarray, Optional[float]], Optional[float]]] = None`. When set, it overrides the yaml margin for that recipe; `None` from it means uncertified.)
- Modify: `src/fuelblend/design/library_route.py`. The viscosity (ML vs Refutas by `visc_ml_status`) and smoke (in-domain margin, else None) models use `margin_fn`. `_rejudge_smoke` and `_rejudge_viscosity` become consistency checks.
- Test: `tests/test_v5_margin_fn.py`

- [ ] **Step 1: Failing tests:**
  - `test_margin_fn_matches_rejudge`: for the 8 stored recipes, the engine status, slack and margin with `margin_fn` equal the post-hoc rejudge to 1e-9.
  - `test_margin_fn_absent_unchanged`: R-SC-4.
  - `test_optimiser_sees_smoke_margin`: a 3-component pool with toluene and a **stub** `margin_fn` returning a fixed 3.0 mm in domain (independent of Task 1.4) makes the optimiser return a recipe whose smoke passes the 18 mm alternative limit (seed 42, `n_starts=5`).
- [ ] **Step 2:** Run. Expected: FAIL.
- [ ] **Step 3:** Implement. `_Sides.vector` calls `margin_fn` per evaluation. Keep the memoisation pattern from commit d1b2b60.
- [ ] **Step 4:** Tests pass, and the full `pytest -q` matches the last full-run baseline (581 passed at delta_smoke_margin, plus later additions) with no new failures.
- [ ] **Step 5:** Commit.

### Task 1.6: Re-run the designs (R-SC-3) with new pools that contain aromatics

- [ ] Run the 4 existing pools plus 2 new pools containing an aromatic blendstock (e.g. v5pool + toluene/1,2,4-TMB, and v5pool + an alkylbenzene C10–C12 with verified constants), using `--synthetic-blend --alt-limits`, seed 42, into `results/v5/BASE/delta_smoke_cert/designs/`. Use `run_in_background` (each run takes 20+ minutes).
- [ ] Append the R-SC-3 table: per design, the status of every constraint, the freezing slack, and the category. Expected outcome, stated now: smoke passes only where the replacing margin (R-SC-2b) is below pred − 18 at the recipe's aromatics (6.97 mm at 8 vol %, 4.71 mm at 15 vol %), and only where the class-mix gate holds. The current butylcyclohexane-heavy recipes are likely to fail that gate, which pushes the optimiser toward alkane-richer recipes and moves freezing. Freezing may fail or block, and that feeds WS2. Certification also needs Task 1.0.
- [ ] Update `docs/claims_register.md` (new CLM) and `docs/model_decision.md` (v5 ML components / design side).

---

## WS2 — Freezing margin valid for multicomponent recipes (roadmap)

Why: the freezing slack is 0.02 on a margin measured on binaries. Recipes with aromatics will move it.
- Data: multicomponent liquidus / wax-disappearance data (Coutinho quaternary n-C11/C14/C16/C18, ACS Omega 2022; Coutinho 2000 jet/diesel wax model paper, Energy & Fuels; Bell/Boehm/Heyne fuel freezing set, Energy & Fuels, WSU). Already listed in `data/raw/freezing_lit` and `sle_constants`, so check what is ingested before fetching anything new.
- Pure constants: Tm/ΔHfus from the Acree–Chickos compendia (JPCRD 2016 C1–C10, 2017 C11–C192, 2022 update) for every pool component and every GCxGC representative. Verified-constant coverage is the gate for the freezing margin.
- Models compared on the same groups: ideal SvL (designated), Coutinho predictive UNIQUAC solid phase, and SvL + an aromatic–alkane liquid activity term (UNIFAC, then GH-GNN from WS5).
- Declared readings: margin by class mix (alkane-only vs alkane+naphthene vs +aromatic); the margin applies only to the class mixes that were validated. Gate: at least 8 independent multicomponent systems per class mix, else the margin stays transferred, which is stated as such.

## WS3 — Pure-component gap filling (ML used for gaps; highest leverage)

Each filled value carries `source = ml_estimated`, its out-of-fold error, and a propagated effect on the mixture prediction (vary by ±1 OOF RMSE and report the change, as in R-MLY-5). A filled value never replaces a measured one.
- **Tm/ΔHfus** (freezing): GNN or group contribution trained on Acree–Chickos plus the existing verified CSV. Gate: OOF MAE in Tm ≤ 5 K on hydrocarbons C8–C20.
- **Antoine/Psat for heavy n-alkanes** (C11 low-T, C13/15/17/19, C18): carbon-number interpolation of the refit (A,B,C) (AGENTS.md §12b item 9), checked against ThermoML 2020 vapour pressures.
- **YSI/TSI for saturated C13–C18 representatives**: delta_ml_ysi failed because the YSI→smoke map under-predicts mixtures by about 7.6 mm. Next attempt: recalibrate `smoke_point_from_ysi` on the mixture data (Mensch binaries, Rainville surrogates, the JP-8 drum), with grouped CV, using the TSI relation TSI = 4.07·MW/SP − 4.8 that the Mensch data obey exactly. This gives a smoke model that works at 0–8 vol % aromatics and above 25 mm, which is the second route to smoke certification if WS1's margin gate fails.
- **LFL**: 25 measured / about 670 estimated; add Zabetakis-style estimates with an error bar. Low priority (flash already passes).

## WS4 — Mixture ML: DeepSets / set transformer / GNN versus the LightGBM pair residual (declared comparison)

Expected outcome stated now: **a loss or a tie is likely.** Viscosity ternaries are 18 systems, density ternaries failed, DCN has 31 systems, and Leenhouts reports roughly 10³ unique mixtures are needed to beat ideal mixing. The value is the evidence either way.
- Baseline to beat: `ml_pair_residual` (viscosity) and `ml_pair_excess_volume` (density) on the same GroupKFold `system_id` folds, the same locked set, and the R-MLV/R-MLD readings.
- Candidates, all trained as residuals on top of the physics value (never raw targets):
  1. **Pairwise DeepSets**: `ln η = ln η_Refutas + Σ_{i<j} x_i x_j φ(h_i, h_j, T)` with a shared MLP φ on component embeddings h (RDKit descriptors + class + pure properties at T). This is the neural version of the current pair term and keeps its physics form.
  2. **MolSets-style GNN + attention pooling** (Zhang et al. 2023, arXiv 2312.16473, code github.com/Henrium/MolSets): molecular graphs, attention weights by fraction.
  3. **Set transformer** (multi-task heads for viscosity, density and flash residual, masked loss).
  Use deep ensembles (5 seeds) for uncertainty and split conformal by group for margins.
- Gate: the paired bootstrap CI of the MAE gain over the pair residual excludes 0 on dev **and** locked, and ternaries do not get worse. Only then does a variant become eligible for designation, which is a user decision.
- Compute: university CS lab GPU (`jobs/v4_gpu/` has the job pattern). The CPU smoke test runs locally first.
- Report the Leenhouts benchmark splits too, so the results are comparable to the published numbers (`delta_leenhouts_splits`).

## WS5 — Non-ideality for flash and SLE: GH-GNN (test, not a fix)

UNIFAC made flash worse (dMAE +0.45 °C, v4 F). GH-GNN (Sanchez Medina et al., Digital Discovery 2023, code github.com/edgarsmdn/GH-GNN) predicts ln γ∞(T). Use it through a two-parameter Margules/Wilson fit per pair, then modified Raoult in Le Chatelier and SVL.
- Declared readings: flash on the lab mixtures (31 systems), CRC 647 r3, and T2 measured; freezing on the WS2 sets. Check the applicability domain (Tanimoto > 0.35, ≥ 25 training systems per class pair) for every pair and report it.
- Gate: improvement with CI excluding 0 and no rise in the unsafe-side rate. Also check whether it removes the single mesitylene system that sets the 2.21 °C flash margin.

## WS6 — Real fuels from composition distributions (GCxGC)

- 113 mappable AJFTD GCxGC fuels + CRC 647 + AV-33-22. Models: physics on representatives (current), + PLS / GP / small NN residual on the class × carbon-number matrix, grouped by fuel family.
- Compare with the published numbers. Illinois deep learning on GCxGC: 106 petroleum + SAF samples, density, distillation, flash, viscosity. Dayton GCxGC-VUV, ten properties with UQ. DLR iso-alkane subgrouping (Lüdtke et al. 2025; Pütz et al. 2025, Energy & Fuels). Bell/Boehm/Heyne freezing (MAE 4.4 °C).
- Target: the cold-viscosity margin (`inf` at −40 °C for reconstructed fuels) and freezing on real fuels, where models now show no skill over a constant.

## WS7 — Design realism and the final check

- Pools of real blendstocks: Jet A base (measured), HEFA-SPK, ATJ-SPK, SAK, each represented by its GCxGC composition, blended with pure components. This matches how SAF blends are certified.
- Add D86 distillation constraints (T10 ≤ 205 °C, FBP ≤ 300 °C; D7566 also has T50−T10 ≥ 15 °C and T90−T10 ≥ 40 °C). Use the existing distillation baseline, with a margin from CRC D86 grids.
- Optimiser: CMA-ES or multi-start SLSQP with an uncertainty penalty and the domain kept as hard constraints.
- **Final check: lab-test 1–2 designs.** Ask the user whether the university lab can measure D56 flash, D445 viscosity at −20 °C, D4052 density, D5972 freezing, and D1322 smoke. One measured designed blend is worth more than any margin.

---

## New data and method leads (web search 2026-11-13; de-duplicated against `docs/data_sources.csv`, `data/raw/*/SOURCES.md`)

| Lead | Use | Access | Status |
|---|---|---|---|
| AJFTD smoke point + naphthalenes for 2,601 fuels | WS1 smoke margin on real fuels | login (granted) | **not exported in v5 scrape**; highest priority |
| ASTM D7566 Table 1 (aromatics min 8 vol %, smoke 25 or 18 + naph 3.0) | WS1 limits | ASTM (paid; free iTeh previews may not include Table 1) | secondary only (D7566-12a); user must verify |
| Prak, Vozka, Kilaz et al., JP-5 composition + surrogate density/viscosity, J. Chem. Eng. Data 2019 (doi 10.1021/acs.jced.8b01233) | WS4 multi-component density/viscosity | NIST ThermoML page (open) | **not in local ThermoML 2020 folder** (9 files under 10.1021); fetch |
| Vozka et al., alternative blending components in Jet A, Energy & Fuels 2019 (ACS figshare 7844741) | WS6 blend properties vs composition | figshare SI (open?) | unverified |
| Illinois, "Deep learning algorithms for assessing sustainable jet fuels from GCxGC" (106 samples) | WS6 benchmark / data | journal; SI unchecked | unverified |
| Dayton GCxGC-VUV property prediction with UQ (1,226-molecule reference database) | WS6 benchmark | open report (rosap.ntl.bts.gov dot/64804, 64805) | unverified |
| DLR Lüdtke et al. 2025, Pütz et al. 2025 (iso-alkane GCxGC subgroups) | WS6 features | elib.dlr.de (open PDFs) | unverified |
| Mendeley p8p6b6tjr5 (ATR-FTIR spectra + properties, 50 SAF blends, listed 10 Nov 2026) | WS6 extra real blends (no GCxGC) | open | very recent; verify content |
| Acree & Chickos JPCRD 2016/2017/2022 fusion enthalpies | WS2/WS3 Tm, ΔHfus | journal (UNT copy of Part 2 open) | Acree is in `sle_constants/SOURCES.md`; check how much is ingested |
| MolSets code (github.com/Henrium/MolSets) | WS4 | open (PyTorch/PyG) | new |
| GH-GNN code (github.com/edgarsmdn/GH-GNN) | WS5 | open | new |
| Kuzhagaliyeva et al. 2022, Commun. Chem. 5:111, AI fuel-mixture design (CC BY) | WS4/WS7 method reference (mixing operator, search) | open (PMC9814251) | already in `data_sources.csv`; method only |
| Hunt 1953, Ind. Eng. Chem. 45:602, pure smoke points | WS3 TSI route | paywalled | already in `smoke_lit/SOURCES.md` |

## Order and dependencies

1. WS1 Tasks 1.1–1.5 (code); Task 1.0 runs in parallel (user).
2. WS1 Task 1.6, then decide: if the smoke margin gate fails, go to the WS3 TSI route; if freezing blocks, go to WS2.
3. WS3 (Tm/ΔHfus, Antoine) feeds WS2 and WS5.
4. WS4 and WS5 are independent of WS1 and can run as separate sessions (GPU).
5. WS6 and then WS7 last; the WS7 lab test is the true final check.
