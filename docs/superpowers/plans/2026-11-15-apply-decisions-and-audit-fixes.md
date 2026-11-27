# Apply the 2026-11-15 decisions and close the audit tail: Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** apply the seven user decisions in `docs/decisions_2026-11-15.md` as declared deltas. Then re-derive the freezing margins, re-run the design judgements and re-score BASE. Also close the open items of `docs/audit_2026-11-15.md`.

**Architecture:**
- Constants change first: bicyclohexyl dHfus, cis-decalin MW and the Coutinho matching.
- Both lab-liquidus freezing margins are then **re-derived by script** on the new constants, before any route switches. Bicyclohexyl is one of the 9 S_F systems, so decision 4 changes the inputs to both margins, and 13.84 C is stale until it is recomputed.
- The designated freezing route (API and design judge) then moves to UNIFAC-SvL, and the flash-registry overlay becomes the default.
- Designs and BASE are re-scored as declared deltas. Manifests are rewritten only on the user's yes.

**Tech Stack:**
- Python 3.12 (`.venv/Scripts/python.exe`): pandas, numpy, scipy (`brentq`), RDKit, `thermo.unifac`, pytest.

**Spec:**
- `docs/decisions_2026-11-15.md` (the decisions).
- `docs/audit_2026-11-15.md` (the audit).
- `docs/overnight_2026-11-15.md` (the readings both of them rest on).

## Global Constraints

- **Gate.** Nothing in Tasks 2-10 starts until the other session's new data has landed **and** the user says go. The user's words: "dont start any reruns yet, another session will add new data soon". Task 1 is docs only and may run earlier.
- **Declarations.** Declare readings before computing them. Every re-run is a delta folder under `results/v5/BASE/delta_*/`, and its README declaration is committed **before** the compute step. Do not quote the "Expected consequence" lines in the decisions doc as results.
- **Gates and margins.** Never move a gate, a spec limit or a spec `verification:` flag. Never hand-type a margin. Margins enter `configs/uncertainty_v5.yaml` only through a builder function that reads a `summary.csv`, following the `scripts/v5/base/build_uncertainty_v5_q2.py::add_freezing` pattern. The existing key `freezing_point_C.lab_liquidus.svl_verified` is never edited or removed.
- **Frozen v4 artefacts.** v4 artefacts (tag `v4-r3-final`, `SHA256SUMS_v4.txt`) are frozen. No change here may alter a v4 script's output.
- **Pinned fixtures.** `tests/data/margin_fn_baseline.json`, `tests/data/smoke_alt_baseline.json` and `tests/data/ws7_feature_off_fixture.json` are not regenerated. Tests that pin the old freezing behaviour are pointed at the legacy model explicitly.
- **Negative results.** Keep them (`docs/v5_negative_results.md`). Do not run lab tests.
- **Subagent models.** Never use the fable model; opus is the maximum.
- **Commits.** Commit steps name explicit paths. `docs/data_audit_2026-11-11.md` (the user's edit) and the untracked files `Fuel_Blend_Methodology.pdf`, `scripts/v5/base/ingest_calorimetry_primary.py` and `tests/test_calorimetry_primary.py` are never staged.
- **Never committed.** PDFs, licensed text, `ajftd_session.json` and raw AJFTD pages.
- **Pushes.** No push without the user's approval.
- **Long runs.** Runs of more than about 2 minutes (design runs, `reproduce_base.py --full`, the full pytest suite of about 22 min) run in the background.

## Review Focus

1. **A recipe member with no UNIFAC groups** (UNIFAC fallback). Freezing must be `uncertified` with warning `unifac_fallback`. It must not be silently scored by ideal SvL under the UNIFAC margin. Test in Task 6.
2. **The overlay explicitly disabled** (`FUELBLEND_FLASH_REGISTRY_ADDITIONS=0`). The registry must equal the base v2 build frame-for-frame. Test in Task 7.
3. **Stale caches across a switch.** `_default_registry`, `predict_route._tmap` and `_tmap_freezing` are `lru_cache`d. A test or script that flips the env var or the freezing model in-process must `cache_clear()`, or it reads the old state. Tests in Tasks 6 and 7 assert the switch takes effect after `cache_clear()`.
4. **A recipe whose only would-be crystallisers are alkylbenzenes** (for example toluene + ethylbenzene) once Coutinho matching is exact. `predict_fuel` must return freezing `value=None` with reason `non_finite`, not raise. Test in Task 4.
5. **A zero-weight pool member under UNIFAC** (the optimiser visits w_i = 0). The freezing value must equal the value with that member removed, to 1e-6 C. Test in Task 6.

---

### Task 0: Gate and data intake (no code)

**Files:** none changed.

- [ ] **Step 1:** Confirm with the user that the other session's data has landed and that applying the decisions may start. If not, stop.
- [ ] **Step 2:** Find what landed. Run `git log --oneline d9309d3..HEAD` and `git status --short`, then read any new README or SOURCES files.
- [ ] **Step 3:** Write a 5-line intake note at the top of `docs/decisions_2026-11-15.md` under `## Applying (date)`, recording which tasks the new data touches:
  - new SLE / liquidus rows: Task 5;
  - new flash or Antoine data: Task 7;
  - new real fuels or GCxGC data: Task 9;
  - anything else: say so.
- [ ] **Step 4:** If the new data changes `data/raw/sle_constants/tm_dhfus_verified.csv`, `data/processed/v5/pure_props_v5.csv` or the flash registry, re-read the affected rows before Tasks 2-4.

### Task 1: Audit doc tail (docs only; may run before the gate)

**Files:**
- Modify: `AGENTS.md` (§5 header and opening line, §11 header and table intro, §12b first line, §12 new paragraph).

The audit (§7) found that §5, §11 and §12b still present the v1 commands (`fetch.py`, `build_unified.py`, ...) and "plain asserts" tests as current.

- [ ] **Step 1:** In §5, change the header to `## 5. Pipeline (HISTORIC v1 commands; current entry points below)`. Add one line under it:
  - `Current: python scripts/v4/reproduce_all.py (v4), python scripts/v5/base/reproduce_base.py [--full] (v5 BASE), .venv/Scripts/python.exe -m pytest -q (tests). The v1 scripts below survive only in extras/ and git history.`
- [ ] **Step 2:** In §11, change the header to `## 11. Tests (HISTORIC v1 list; the suite now runs with pytest)`. Add one line:
  - `Current: .venv/Scripts/python.exe -m pytest -q (741 passed, 11 skipped on 2026-11-15, docs/audit_2026-11-15.md §1).`
- [ ] **Step 3:** In §12, append a paragraph:
  - `v5 overnight 2026-11-15: docs/overnight_2026-11-15.md (D86 Rayleigh opt-in, UNIFAC-SvL freezing, P0y-b conditional, strict_constraints guard); user decisions in docs/decisions_2026-11-15.md; their application plan docs/superpowers/plans/2026-11-15-apply-decisions-and-audit-fixes.md.`
  - Leave the "plan not yet committed" sentence until the user decides item U1 (see the User-gated items section).
- [ ] **Step 4:** Commit.

```bash
git add AGENTS.md
git commit -m "docs(AGENTS): mark v1 pipeline/test sections historic, point to current entry points; add overnight 2026-11-15 pointer"
```

### Task 2: Decision 6, `strict_constraints=True` by default

**Files:**
- Modify: `src/fuelblend/design/blendstock_v5.py:534` (signature default) and the docstring line at :538.
- Create: `scripts/v5/base/check_sides_excluded.py` (pre-check, no optimisation).
- Test: `tests/test_delta_pool_bicyclohexyl.py` (add one test).

**Interfaces:**
- Produces: `design_pool(..., strict_constraints: bool = True)`.
- `check_sides_excluded.main() -> int` prints `pool, spec, constraint, reason` for every `_Sides.excluded` entry whose reason does not start with `"no model"`. It returns the number of such entries.

- [ ] **Step 1: Pre-check.** Build the models and `E._Sides(spec, pool, use, margins, graded=..., alt_point_standin=...)` exactly as `blendstock_v5.design_pool` does, with no optimisation. Do this for every pool behind a committed design run, on `astm:jet_a` and `bis:jet_a1`:
  - P0, P0x-a, P0x-b and P0y-b: pool definitions in `scripts/v5/base/delta_smoke_branch_opt.py` and `scripts/v5/base/delta_pool_bicyclohexyl.py`, with the overlay env as each committed run used it.
  - P1, P2 and P3: `scripts/v5/base/delta_ws7_design_realism.py`.
- [ ] **Step 2:** Run `.venv/Scripts/python.exe scripts/v5/base/check_sides_excluded.py`. Expected: only P0y lists entries (flash, nhoc; `decalin_cis.MW status=unverified`).
  - If any other committed pool lists one, stop and report it to the user: flipping the default would make that frozen run irreproducible.
- [ ] **Step 3: Write the failing test.**

```python
def test_strict_constraints_is_default():
    import inspect
    from fuelblend.design import blendstock_v5 as B
    assert inspect.signature(B.design_pool).parameters["strict_constraints"].default is True
```

- [ ] **Step 4:** Run `.venv/Scripts/python.exe -m pytest tests/test_delta_pool_bicyclohexyl.py -q`. Expected: the new test FAILS.
- [ ] **Step 5:** Flip the default and update the docstring: "default True by user decision 2026-11-15 (decisions doc #6)".
- [ ] **Step 6:** In `scripts/v5/base/delta_pool_bicyclohexyl.py:109`, pass `strict_constraints=False` explicitly, but only when `--pool P0y` is used, so the committed R-BC-2 run stays reproducible. Add a comment naming R-BC-2.
- [ ] **Step 7:** Run the tests. Expected: all pass.

```bash
.venv/Scripts/python.exe -m pytest tests/test_delta_pool_bicyclohexyl.py tests/test_v5_sak_blendstock.py tests/test_v5_ws7_design_realism.py tests/test_v5_design_forward.py tests/test_delta_smoke_branch_opt.py -q
```

- [ ] **Step 8:** Commit `src/fuelblend/design/blendstock_v5.py`, `scripts/v5/base/check_sides_excluded.py`, `scripts/v5/base/delta_pool_bicyclohexyl.py` and `tests/test_delta_pool_bicyclohexyl.py`.

### Task 3: Decisions 7 and 4, constants (cis-decalin MW; bicyclohexyl dHfus)

**Files:**
- Modify: `scripts/v5/audit_pure_constants_v5.py` (`build()`, add a third loop).
- Regenerate: `data/processed/v5/pure_props_v5.csv` and `results/v5/P1/pure_props_changes_v5.csv`.
- Modify: `data/raw/sle_constants/tm_dhfus_verified.csv`, bicyclohexyl row (`WVIIMZNLDWSIRH-UHFFFAOYSA-N`), plus its `SOURCES.md` if present in that folder.
- Create: `tests/test_decisions_2026_11_15_constants.py`.

**Interfaces:**
- Produces: `audit_pure_constants_v5.FORMULA_MW_VERIFIED: dict[str, str] = {"decalin_cis": "C10H18"}`.
- The verified tmap entry for bicyclohexyl has `dh_fus_J_mol == 13940.0`.

- [ ] **Step 1: Write the failing tests.**

```python
import subprocess, io
import pandas as pd
from fuelblend.eval.composition_route import build_tmap_verified

def test_cis_decalin_mw_verified_formula():
    pp = pd.read_csv("data/processed/v5/pure_props_v5.csv")
    r = pp[(pp.compound == "decalin_cis") & (pp.prop == "MW")].iloc[0]
    assert r.status == "verified" and abs(r.ref_value - 138.25) < 0.01 and "C10H18" in r.ref_source

def test_pure_props_only_cis_decalin_mw_changed():
    old = pd.read_csv(io.StringIO(subprocess.check_output(["git", "show", "d9309d3:data/processed/v5/pure_props_v5.csv"], text=True)))
    new = pd.read_csv("data/processed/v5/pure_props_v5.csv")
    diff = old.merge(new, how="outer", indicator=True).query("_merge != 'both'")
    assert set(zip(diff.compound, diff.prop)) == {("decalin_cis", "MW")}

def test_bicyclohexyl_svl_dhfus_total():
    t = build_tmap_verified()["WVIIMZNLDWSIRH-UHFFFAOYSA-N"]
    assert t["dh_fus_J_mol"] == 13940.0 and t["tm_K"] == 276.8
```

- [ ] **Step 2:** Run `.venv/Scripts/python.exe -m pytest tests/test_decisions_2026_11_15_constants.py -q`. Expected: 2 FAIL (cis-decalin MW, bicyclohexyl dH); `only_changed` FAILS or passes trivially.
- [ ] **Step 3: cis-decalin MW.** In `build()`, after the listings loop, add a loop over `FORMULA_MW_VERIFIED`:
  - For each compound's `MW` row: `ref_value = Descriptors.MolWt(MolFromSmiles(smiles))` (RDKit, 138.2499) and `ref_source = "formula C10H18 (RDKit MolWt); user decision 2026-11-15, docs/decisions_2026-11-15.md #7"`.
  - Then set `delta` and `status` through `_status("MW", listed, ref)` and `log(i, "MW reference -> formula (user decision)")`.
  - Probable cause of the gap, to record in the comment: the chemicals lookup that gave `decalin_trans` its `ref_value` missed the cis stereo SMILES.
  - Regenerate with `.venv/Scripts/python.exe scripts/v5/audit_pure_constants_v5.py`.
- [ ] **Step 4: bicyclohexyl dHfus.** In `tm_dhfus_verified.csv`, bicyclohexyl row: set `svl_dhfus_choice` from `fusion` to `total` and `svl_dhfus_kJ_mol` from 6.86 to 13.94. Append to `note`: `SvL uses fusion + II->I transition (13.94) by user decision 2026-11-15 (#4)`. No other cell changes.
- [ ] **Step 5:** Run the Step 2 command. Expected: 3 PASS. If `only_changed` lists more rows, the generator's other inputs moved since d9309d3 (possibly the new data). Stop and report the rows; do not commit.
- [ ] **Step 6:** Run `.venv/Scripts/python.exe -m pytest tests/test_delta_freeze_margin.py tests/test_predict_route.py -q`. If a test pins a bicyclohexyl-dependent freezing value, it now fails. Record which ones; they are fixed in Task 6 (legacy pointer), not by editing values here.
- [ ] **Step 7:** Commit the script, the two regenerated csvs, `tm_dhfus_verified.csv` and the test.

### Task 4: Decision 5, exact aromatic matching in the Coutinho table (v5 callers only)

**Files:**
- Modify: `src/fuelblend/thermo/coutinho_wax.py:26-61`.
- Modify: `src/fuelblend/eval/composition_route.py:150-158` (`build_tmap`).
- Test: `tests/test_coutinho_exact_match.py` (create).

**Interfaces:**
- Produces: `get_pure_thermal_properties(comp_df=None, aromatic_match: str = "substring")`, with `"substring"` (legacy, v4 callers) and `"exact"`.
- Produces: `build_tmap(reg, idx, aromatic_match: str = "exact")`.
- `"exact"` compares the RDKit canonical SMILES of the row to the canonical SMILES of the four special molecules:
  - durene `Cc1cc(C)c(C)cc1C`;
  - naphthalene `c1ccc2ccccc2c1`;
  - p-xylene `Cc1ccc(C)cc1`;
  - benzene `c1ccccc1`.
- Every other aromatic gets `nan, nan`. Non-aromatic branches are unchanged.

- [ ] **Step 1: Write the failing tests.** Use a 6-row `comp_df` (benzene, p-xylene, naphthalene, toluene, fluorene `c1ccc2c(c1)Cc1ccccc1-2`, phenanthrene `c1ccc2c(c1)ccc1ccccc12`; `cls="aromatic"`, `n_carbon` set):

```python
def test_exact_keeps_specials_and_drops_others(df6):
    t = get_pure_thermal_properties(df6, aromatic_match="exact")
    assert t["BENZ"]["tm_K"] == 278.7 and t["PXYL"]["tm_K"] == 286.4 and t["NAPH"]["tm_K"] == 353.4
    assert all(np.isnan(t[k]["tm_K"]) for k in ("TOL", "FLUO", "PHEN"))

def test_substring_default_unchanged(df6):
    t = get_pure_thermal_properties(df6)
    assert t["TOL"]["tm_K"] == 278.7 and t["FLUO"]["tm_K"] == 278.7      # legacy defect kept for v4 reproduction

def test_no_misassigned_under_exact():
    # rows of results/v5/BASE/delta_coutinho_lookup_audit/misassigned.csv: none may keep its assigned_tm_K under "exact"
    ...  # build comp_df from the csv (smiles, cls="aromatic"); assert every tm_K is nan

def test_predict_fuel_alkylbenzene_only_freezing_none():           # Review Focus 4
    out = predict_fuel({"toluene": 0.5, "ethylbenzene": 0.5}, basis="mass")
    fp = out["freezing_point_C"]
    assert fp["value"] is None or fp["status"] != "ok" or fp["uncertified"]   # no exception; verified toluene may still give a value
```

Pin the last assert to what the code actually does once it runs. It must not raise, and `uncertified` must be truthy unless every component is verified.
- [ ] **Step 2:** Run `.venv/Scripts/python.exe -m pytest tests/test_coutinho_exact_match.py -q`. Expected: FAIL (unknown keyword).
- [ ] **Step 3:** Implement the keyword in `get_pure_thermal_properties`, and pass `aromatic_match="exact"` from `build_tmap`. Leave every `scripts/v4/*` call and `eval/freezing_f1.py:60` as they are. Before changing `freezing_f1`, grep its callers: it is shared with v4.
- [ ] **Step 4:** Run the new test file plus the suite below. Expected: all pass.

```bash
.venv/Scripts/python.exe -m pytest tests/test_coutinho_exact_match.py tests/test_predict_route.py tests/test_delta_freeze_margin.py tests/test_v5_design_forward.py -q
```

  If a v5 test pins a Coutinho-fallback point value (the P2/P3 printed freezing values), record it for Task 9 rather than editing it.
- [ ] **Step 5:** Run `.venv/Scripts/python.exe scripts/v4/reproduce_all.py --verify-only`. Expected: the only drift is the already-known `docs/claims_register.md` hash (U3). Any other v4 drift means a v4 path was touched; revert it.
- [ ] **Step 6:** Commit the two src files and the test.

### Task 5: Re-derive both lab-liquidus freezing margins (declared delta `delta_freeze_adopt`)

**Files:**
- Create: `results/v5/BASE/delta_freeze_adopt/README.md` (declaration first).
- Create: `scripts/v5/base/delta_freeze_adopt.py`.
- Modify: `scripts/v5/base/build_uncertainty_v5_q2.py` (add `add_freezing_unifac(doc)`).
- Modify: `configs/uncertainty_v5.yaml` (new key only, written by the builder).

**Interfaces:**
- Consumes: Task 3 constants (`build_tmap_verified()`).
- Consumes: `FMG.scope_rows(tmap)` and `FMG.SLE_CSV` from `scripts/v5/base/delta_freeze_margin.py`.
- Consumes: `FM.svl_activity_C(x, consts, gamma_fn)` and `make_unifac_gamma_fn(smiles) -> (gamma_fn, fallback_idx)`, the same calls as `scripts/v5/base/delta_freeze_unifac.py:40-58`.
- Produces: `results/v5/BASE/delta_freeze_adopt/summary.csv` (same columns as `delta_freeze_unifac/summary.csv`), `rows.csv` and `readings.json` with keys `margin_ideal_C`, `margin_unifac_C`, `mae_ideal_C`, `mae_unifac_C`, `n_systems`, `n_rows`, `unifac_fallbacks`, `rule_prefers_unifac`, `stop_reason`.
- Produces: yaml key `freezing_point_C.lab_liquidus.svl_unifac`.

- [ ] **Step 1: Write and commit the declaration** before any compute. The README declares:
  - **Scope S_F.** As in `delta_freeze_unifac`: ThermoML liquidus, C8-C16 HCs with verified constants, T <= 273.15 K. One row per liquidus point; the 10 duplicate rows found overnight are dropped (67 rows on the old data). Plus any new SLE systems from Task 0 that meet the same scope rule, listed by system_id.
  - **Constants.** `tm_dhfus_verified.csv` after Task 3, including bicyclohexyl 13.94.
  - **R-FA-1.** Ideal and UNIFAC margins (conformal, one-sided under, alpha 0.1, k = ceil((n+1)(1-alpha)) over systems) and MAE.
  - **Preference rule (unchanged from `delta_freeze_unifac`).** UNIFAC is designated only if its margin < the ideal margin AND its MAE <= the ideal MAE, with 0 UNIFAC fallbacks.
  - **Stop rule.** If the rule does not prefer UNIFAC, or n_systems < 9, write `stop_reason`, write no yaml key, and stop the plan at Task 6 to ask the user. The decision was made on 13.84 vs 16.98, so a reversal goes back to them.
  - **R-FA-2 (descriptive).** The 40-system set, as R-FU-2.
- [ ] **Step 2:** Commit the README alone with the message `v5 delta_freeze_adopt: declaration (before compute)`.
- [ ] **Step 3:** Implement `delta_freeze_adopt.py`. Reuse `delta_freeze_unifac.main()`'s row loop on the new tmap. Drop duplicate rows the same way `delta_freeze_unifac` did to get 67 rows from 77; read its row handling and reuse it, don't invent a key. Assert 67 rows on the pre-new-data SLE file. Compute both models and write the outputs.
- [ ] **Step 4:** Run `.venv/Scripts/python.exe scripts/v5/base/delta_freeze_adopt.py`. Expected: `readings.json` written. Report the margin and MAE numbers to the user as readings.
- [ ] **Step 5:** Write `add_freezing_unifac(doc)`, mirroring `add_freezing`:
  - Raise unless `readings["rule_prefers_unifac"]`.
  - Assert the key is absent.
  - Fill `value`, `applied_margin`, `n_groups`, `n_rows` and `residual_file` from `summary.csv` row `S_F, svl_unifac`.
  - Set `model`: "UNIFAC (thermo, Dortmund groups as flash_unifac.make_unifac_gamma_fn) Schroeder-van Laar, verified Tm / dHfus".
  - Copy the scope and caveats of the ideal key, and add: "exploratory: ideal residuals were seen before the delta_freeze_unifac declaration".
  - Set `user_decision`: "2026-11-15 'Adopt UNIFAC-SvL (Recommended)' (docs/decisions_2026-11-15.md #1)".
  - Run it with the existing builder entry point.
- [ ] **Step 6:** Run `.venv/Scripts/python.exe -m pytest tests/test_delta_freeze_margin.py -q`, then `git diff configs/uncertainty_v5.yaml`. Expected: one added block; the `svl_verified` block is byte-identical.
- [ ] **Step 7:** Commit the script, the results folder, the builder and the yaml.

### Task 6: Decision 1, designated freezing route becomes UNIFAC-SvL (API and design judge)

**Files:**
- Modify: `src/fuelblend/eval/composition_route.py`: add `FREEZING_MODELS`, `DESIGNATED_FREEZING_MODEL` and `freezing_C`.
- Modify: `src/fuelblend/api/predict_route.py:729-790` (`FREEZING_MARGIN_KEY` becomes a lookup; `_freezing`).
- Modify: `src/fuelblend/design/library_route.py:208-240` (`FREEZING_MARGIN_KEY`; `_freezing_api_fn`).
- Modify (legacy pointer only): `tests/test_delta_freeze_margin.py`, `tests/test_predict_route.py`, `tests/test_v5_design_forward.py`, and any test listed in Task 3 Step 6. Each pins the legacy model with a fixture that sets `composition_route.DESIGNATED_FREEZING_MODEL = "ideal_svl"` and clears the caches.
- Create: `tests/test_freezing_unifac_route.py`.

**Interfaces:**
- `FREEZING_MODELS = {"ideal_svl": "freezing_point_C.lab_liquidus.svl_verified", "svl_unifac": "freezing_point_C.lab_liquidus.svl_unifac"}`: model name to margin key.
- `DESIGNATED_FREEZING_MODEL: str = "svl_unifac"`.
- `freezing_C(x, iks, smiles, tmap, model: str | None = None) -> tuple[float, dict]` returns `(value_C, info)`, where `info = {"model": str, "unifac_fallback": list[int]}`.
  - `"ideal_svl"` is `ideal_sle_C(x, iks, tmap)`.
  - `"svl_unifac"` is `svl_activity_C(x, [(tm, dh) or None per ik], gamma_fn)`, with `gamma_fn` from `make_unifac_gamma_fn(smiles)`.
- `predict_route.freezing_margin_key() -> str` and `library_route.freezing_margin_key() -> str` both return `FREEZING_MODELS[DESIGNATED_FREEZING_MODEL]`. The two old module constants are removed so they cannot drift apart.
- The API output for `freezing_point_C` adds:
  - `method` = `"svl_unifac_mole_fraction"`;
  - `variants["ideal_svl"] = {"value", "safe_bound", "label"}` (the ideal value with its 16.98 C margin, reported only);
  - the existing `variants["coutinho_constants"]`, unchanged.
- The design library caches `gamma_fn` per `tuple(pool)` on `L.__dict__["_gamma"]`, because `make_unifac_gamma_fn` is too slow to call per optimiser evaluation.

- [ ] **Step 1: Write the failing tests.**

```python
def test_designated_model_is_unifac(): assert CR.DESIGNATED_FREEZING_MODEL == "svl_unifac"

def test_margin_key_single_source():
    assert predict_route.freezing_margin_key() == library_route.freezing_margin_key() == "freezing_point_C.lab_liquidus.svl_unifac"

def test_api_freezing_unifac_with_ideal_variant():
    out = predict_fuel({"n-decane": 0.5, "n-dodecane": 0.3, "trans-decalin": 0.2}, basis="mass")   # all verified
    fp = out["freezing_point_C"]
    assert fp["method"] == "svl_unifac_mole_fraction" and fp["uncertified"] is False
    assert abs(fp["safe_bound"] - fp["value"] - MARGIN_UNIFAC) < 1e-9      # MARGIN_UNIFAC read from the yaml
    assert fp["variants"]["ideal_svl"]["value"] is not None

def test_unifac_fallback_is_uncertified(monkeypatch):                    # Review Focus 1
    # monkeypatch flash_unifac.unifac_subgroups to return None for one member
    ...; assert fp["uncertified"] and "unifac_fallback" in " ".join(fp["warnings"])

def test_zero_weight_member_matches_removed():                            # Review Focus 5
    v1, _ = CR.freezing_C(x=[0.6, 0.4, 0.0], iks=..., smiles=..., tmap=tmap)
    v2, _ = CR.freezing_C(x=[0.6, 0.4], iks=..., smiles=..., tmap=tmap)
    assert abs(v1 - v2) < 1e-6

def test_switch_after_cache_clear(monkeypatch):                           # Review Focus 3
    monkeypatch.setattr(CR, "DESIGNATED_FREEZING_MODEL", "ideal_svl"); predict_route._tmap_freezing.cache_clear()
    assert predict_fuel(BLEND)["freezing_point_C"]["method"] == "ideal_svl_mole_fraction"
```

- [ ] **Step 2:** Run `.venv/Scripts/python.exe -m pytest tests/test_freezing_unifac_route.py -q`. Expected: FAIL.
- [ ] **Step 3:** Implement `freezing_C` and the two `freezing_margin_key()` functions. Switch `_freezing` and `_freezing_api_fn` to them. The UNIFAC fallback rule:
  - A non-empty `info["unifac_fallback"]` among members with x > 1e-4 makes `_freezing_margin_v5` return False.
  - Add warning `unifac_fallback: <names>`.
  - In the design route, record `L.nonverified_used.add(f"{n}.UNIFAC=fallback")`, so that the judge leaves freezing uncertified.
- [ ] **Step 4:** Add the legacy-pointer fixture to the tests listed in Files.
- [ ] **Step 5:** Run the new tests, the legacy-pointed tests and the design suites. Expected: all pass, and the pinned JSON fixtures still match.

```bash
.venv/Scripts/python.exe -m pytest tests/test_freezing_unifac_route.py tests/test_delta_freeze_margin.py tests/test_predict_route.py tests/test_v5_design_forward.py tests/test_v5_ws7_design_realism.py tests/test_delta_smoke_branch_opt.py tests/test_delta_api_route.py -q
```

- [ ] **Step 6: Timing check.** Time one `models["props"]["freezing"].fn(pool, w, None)` call on the P0y-b Jet A recipe, under the ideal model and under UNIFAC, after warm-up. Record both in the Task 9 README. If UNIFAC is more than 20x slower, report it before Task 9: a P0y-b spec took about 25 min under ideal.
- [ ] **Step 7:** Commit the src files and tests.

### Task 7: Decision 3, flash-registry overlay on by default

**Files:**
- Modify: `src/fuelblend/baselines/flash_unifac.py:117-134`.
- Test: `tests/test_delta_pool_bicyclohexyl.py` (invert the two overlay tests).

**Interfaces:**
- `_default_registry()` merges `REGISTRY_ADDITIONS_CSV` unless `os.environ.get(REGISTRY_ADDITIONS_ENV) == "0"`.
- Docstring: "default on by user decision 2026-11-15 (#3); =0 restores the v2 registry".
- `data/processed/v2/compounds_anchored_vp.csv` is not edited.

- [ ] **Step 1: Rewrite the tests.**
  - `test_overlay_on_by_default`: env unset, `cache_clear()`, then bicyclohexyl and cis-decalin rows carry `antoine_source == "anchored_thermo_ambient_v5add"`.
  - `test_overlay_off_equals_base` (Review Focus 2): env `"0"`, `cache_clear()`, then `reg.equals(base.drop_duplicates("can").set_index("can"))`. This is the existing frame-equality assert, moved under `=0`.
  - Both tests clear the cache in teardown.
- [ ] **Step 2:** Run `.venv/Scripts/python.exe -m pytest tests/test_delta_pool_bicyclohexyl.py -q`. Expected: `on_by_default` FAILS.
- [ ] **Step 3:** Implement the default flip.
- [ ] **Step 4:** Run the 109-test API/flash group from `docs/overnight_2026-11-15.md` (Tests run tonight) plus `tests/test_delta_pool_bicyclohexyl.py`. Expected: all pass. A failure here is a value the overlay moves; list it for Task 9 or Task 10, do not edit it.
- [ ] **Step 5:** Commit `flash_unifac.py` and the test.

### Task 8: Claims and decision status (after Tasks 2-7, before re-runs)

**Files:**
- Modify: `docs/claims_register.md` (CRLF line endings, as CLM-69..72).
- Modify: `docs/decisions_2026-11-15.md` (add an `applied in` column: commit per decision).

- [ ] **Step 1:** Add CLM-73, "designated freezing route UNIFAC-SvL", with margin = the Task 5 reading, `user_decision` and status `adopted`. Add a one-line "superseded as designated by CLM-73" note to CLM-70.
- [ ] **Step 2:** Fill the `applied in` column for #1 and #3-#7 with commit hashes; #2 reads `no change`.
- [ ] **Step 3:** Commit both docs.

### Task 9: Design re-runs under the adopted defaults (declared delta `delta_adopt_designs`)

**Files:**
- Create: `results/v5/BASE/delta_adopt_designs/README.md` (declaration first).
- Create: `scripts/v5/base/delta_adopt_designs.py`.

**Interfaces:**
- Consumes: `blendstock_v5.design_pool(spec, bp, models, W.MARGINS, n_starts=10, seed=42, d86=True, lam=B.LAMBDA_UNC, alt_standin=..., strict_constraints=True)` with `alt_standin` as in `delta_smoke_branch_opt`.
- Pools: P0, P0x-a, P0x-b, P0y (now runnable: cis-decalin MW verified), P0y-b, P1, P2 and P3, built as in their committed scripts.

- [ ] **Step 1: Write and commit the declaration** before any run.
  - **Readings, per pool x spec:** judge status, min slack, recipe, flash, freezing (UNIFAC value, safe bound and the `ideal_svl` variant), Engler and Rayleigh D86, and the dicycloparaffin mass fraction (reported only; no gate).
  - **Conditions table**, as in `delta_pool_bicyclohexyl`. A design certified only under Engler is labelled `d86_model_dependent` (decision 2 keeps Engler designated and Rayleigh as a reported check).
  - **Compare:** each design against its committed reading, with the delta per field.
  - **Values to confirm or refute** (from the decisions doc, not quoted as results): P0y-b Jet A stays judge-certified; P0y-b Jet A-1 fails freezing.
- [ ] **Step 2:** Commit the README alone.
- [ ] **Step 3:** Run in the background, one process per spec, writing `designs.csv` and `log_<pool>_<spec>.txt`:

```bash
.venv/Scripts/python.exe scripts/v5/base/delta_adopt_designs.py --spec astm:jet_a
.venv/Scripts/python.exe scripts/v5/base/delta_adopt_designs.py --spec bis:jet_a1
```

  Expected: 16 rows and no `ModelUnavailable` from the strict guard. If the guard raises for a pool, that is a reading: record it and continue with the other pools.
- [ ] **Step 4:** Write the results section and the conditions table; keep every failure.
- [ ] **Step 5:** Commit the script, `designs.csv`, the README and the per-design csvs. Do not commit `log_*.txt` files over 5 MB.

### Task 10: BASE re-score (declared delta `delta_adopt_rescore`)

**Files:**
- Create: `results/v5/BASE/delta_adopt_rescore/README.md` (declaration first) and `drift.csv`.

Decisions 3 (overlay) and 5 (Coutinho) can move frozen T2 / WS7 / P2 / P3 numbers. Decision 1 moves every freezing cell that uses the designated route.

- [ ] **Step 1: Write and commit the declaration.**
  - Readings: every `results/v5/BASE/**` file that `reproduce_base.py --full` regenerates differently from `SHA256SUMS_v5_base.txt`. For each changed scoreboard cell, the old value, the new value and the decision that caused it (attributed by re-running with only that switch reverted where the cell is ambiguous).
  - Rule: the frozen BASE files and the manifest are **not** overwritten in this task.
- [ ] **Step 2:** Copy `results/v5/BASE` (excluding `delta_*`) to the scratchpad as `BASE_frozen`.
- [ ] **Step 3:** Run `.venv/Scripts/python.exe scripts/v5/base/reproduce_base.py --full` in the background (hours). Expected: exit 1 with drift.
- [ ] **Step 4:** Write `drift.csv` (`path, cell, old, new, cause`) and the README results.
- [ ] **Step 5:** Restore the frozen files from the scratchpad copy. Confirm with `.venv/Scripts/python.exe scripts/v5/base/reproduce_base.py`, which must show no drift.
- [ ] **Step 6:** Commit the delta folder. Ask the user item U4.

### Task 11: Final verification and report

- [ ] **Step 1:** Run `.venv/Scripts/python.exe -m pytest -q` in the background (about 22 min). Expected: 0 failed. If the count differs from 741 passed / 11 skipped, list each new pass or failure.
- [ ] **Step 2:** Update the memory file `overnight-2026-11-15-findings.md` (decisions applied, with commits) and its `MEMORY.md` line.
- [ ] **Step 3:** Report to the user:
  - the Task 5 margins;
  - the Task 9 table;
  - the Task 10 drift summary;
  - the user-gated items below.
  - Do not push.

---

## User-gated items (from the audit; listed, not done)

| id | item | source | what an agent does on "yes" |
|---|---|---|---|
| U1 | commit `docs/superpowers/plans/2026-11-12-v5-api-route-and-design.md` and `2026-11-12-v5-full-baseline.md` | audit §6.3, §7 | `git add` those two paths, commit, update the AGENTS §12 "not yet committed" sentence |
| U2 | untrack `data/raw/sak/dot_64804/dot_64804_DS1.pdf` (2.9 MB; the rule is never commit PDFs) | audit §7 | `git rm --cached` it, add it to `.gitignore`, add a `SOURCES.md` pointer in that folder (the file stays on disk and in history) |
| U3 | rehash `docs/claims_register.md` in `SHA256SUMS_v4.txt` (drift confirmed: committed ce613ab6..., current 82273fbb...) | audit §7.1 | `reproduce_all.py --update-manifest` only after checking that no other line drifts |
| U4 | adopt the Task 10 re-score as the new frozen BASE | Task 10 | re-run `--full`, `reproduce_base.py --update-manifest`, commit |
| U5 | push the commits of this plan to `origin/v5` | push rule | `git push -u origin v5` |
| U6 | question, not a fix: add a dicycloparaffin bound (0-30 %, the 113-fuel range) to the smoke class-mix domain rule (`smoke_margin_v5.class_mix_in_real_fuel_range`)? This would be a new gate, so it is the user's call; the overnight Jet A-1 optimum had 37.5 % | overnight observation | declare it as its own delta, opt-in flag, re-judge |
| U7 | question: turn on `alt_standin` (audit §3 fix, CLM-67) by default in design runs? | audit §3 | flip the default with a "flag off = unchanged" test |

Out of scope: the audit's §8 open research items (cold-viscosity margin `inf`, a freezing margin from binaries only, diesel unvalidated, the flash margin from one lab system), unless the new data from Task 0 bears on one of them.
