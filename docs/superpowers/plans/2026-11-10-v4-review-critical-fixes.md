# v4 Review Critical Fixes Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Fix the three critical findings of the 2026-11-10 review: unrealistic CRC 647 reconstructed compositions, the CRC-fitted viscosity recalibration applied to known-composition blends, and design margins calibrated on the wrong population.

**Architecture:** Add a second reconstruction ("r2") beside the existing one, never overwriting v4 data, and re-run the existing evaluators against it through a new `--compositions` flag. Add a known-composition lab evaluation that produces new conformal margin keys. Point the design engine and API at raw Refutas and the lab margins for known-composition blends.

**Tech Stack:** Python 3.12 (project venv `.venv/Scripts/python.exe`), pandas 3, numpy, rdkit, pytest, PyYAML.

**Spec:** No separate spec file. The findings this plan implements are summarised in "Background" below; they come from the review of `docs/PROJECT_SUMMARY_v4.md` on 2026-11-10.

## Background (the findings)

- **F1/F2.** `extras/fuelblend_data/scripts/build_real_fuels_dataset.py:199` gives every CRC fuel the same saturate split (0.45 n / 0.35 iso / 0.20 cyclo), and `fuel_hindcast.Reps.get` takes the first registry compound in each (class, carbon number) bin. Result: 47 % of mass is n-paraffin; cycloalkanes are JP-10 and diamantane; aromatics include biphenyl, fluorene, anthracene, pyrene; 11.5 % is perdeuterated alkane. The density bias (-11 kg/m3) and the +18 C freezing bias follow from this, not from the mixing rules.
- **F3.** `src/fuelblend/design/models_v4.py` (`visc_recal`) and `src/fuelblend/api/predict_blend.py:583` multiply Refutas by a factor fitted on those reconstructed fuels. Pure n-decane at -20 C comes out 3.26 cSt instead of 2.30, and the API result jumps by about 40 % between one component and two.
- **F4.** Density (22.25 kg/m3) and viscosity margins are calibrated on reconstructed-fuel residuals but applied to designed blends whose composition is exact.

**Known consequence, decided up front:** lab mixture viscosity below 263 K covers only 7 systems (40 records in `data/unified/measurements_dedup.csv`). The lab viscosity margin for that temperature bin will be `inf`, so the jet -20 C viscosity constraint will become `uncertified` for designs. That is the honest result. Do not substitute a number or borrow another bin's margin.

## Global Constraints

- Run everything with `.venv/Scripts/python.exe`; system Python lacks rdkit.
- Never modify `data/real_fuels/v4/compositions_v4.csv`, `splits/v4/real_fuels_v4_splits.json`, or any existing file under `results/v4/{B,F,V,D,U,E,O,SP,SLE,P1,P2}`.
- New outputs go to `data/real_fuels/v4/compositions_crc647_r2.csv`, `results/v4/R2/`, `results/v4/L/`, `results/v4/D2/`, `splits/v4/lab_systems_splits.json`.
- Existing keys in `configs/uncertainty_v4.yaml` must keep their exact values; this plan only adds keys.
- The reconstruction may read only descriptor rows (aromatics, diaromatics, olefins, D86). It must never read a target property (density, viscosity, flash, freezing, heat of combustion).
- No choice in this plan may be made by looking at locked fuels or locked lab systems. Locked slices are scored once, at the end of the task that owns them.
- A margin of `inf` is never replaced by a number. Consumers read `applied_margin` only.
- Hydrocarbon-only; residual convention `predicted - measured`.
- pandas 3: assign mixed-type values per column, not via `.loc[k, [cols]] = [...]`.

## Review Focus

1. A (class, carbon number) bin with no realistic representative: the mass must be reported as fallback mass, never silently given to a cage or polycyclic compound. (Test in Task 1.)
2. A fuel with fewer than 4 D86 points or no aromatics value: skipped with a logged reason, no crash. (Test in Task 1.)
3. A lab system containing a component with no measured pure density: excluded and counted in the README, not filled from a class regression. (Test in Task 3.)
4. A temperature bin whose margin is `inf`: the engine returns `uncertified`, and does not use a neighbouring bin. (Test in Task 4.)
5. A temperature exactly on a bin edge (263.15 K, 283.15 K): belongs to the warmer bin, deterministically. (Test in Task 3 and Task 4.)

---

### Task 1: Realistic representatives and the r2 reconstruction

**Files:**
- Create: `src/fuelblend/data/reconstruction_r2.py`
- Create: `scripts/v4/rebuild_crc_reconstruction.py`
- Test: `tests/test_v4_recon_r2.py`
- Output: `data/real_fuels/v4/compositions_crc647_r2.csv`, `results/v4/R2/reconstruction_report.md`

**Interfaces:**
- Consumes: `fuelblend.mixture.fuel_hindcast.carbon_distribution`, `RHO_CLASS`; `fuelblend.mixture.common.CompoundRegistry`, `PurePropertyProvider`; `fuelblend.mixture.cold_flow.compound_meta`.
- Produces:
  - `is_realistic(smiles: str, kind: str) -> bool`, `kind` in `{"n_paraffin","iso_paraffin","cycloalkane","aromatic_mono","aromatic_di","alkene_alkyne"}`
  - `class RealisticReps` with `get(kind: str, n_carbon: int, max_gap: int = 3) -> int | None` (registry index)
  - `nominal_sat_split(jet_class_csv: Path) -> tuple[float, float, float]` (n, iso, cyclo; sums to 1)
  - `reconstruct(row: Mapping, sat_split: tuple[float, float, float], reps: RealisticReps) -> tuple[pd.DataFrame, float]` returning rows `inchikey, smiles, fraction, kind` (mass basis, sums to 1 minus fallback) and the fallback mass fraction
  - CSV with the same columns as `compositions_v4.csv`, `provenance = "reconstructed_r2"`, plus one `LUMP_NO_REALISTIC_REP` row per fuel carrying the fallback mass.

Rules for `is_realistic` (RDKit; all kinds reject any atom with a non-zero isotope):
- `n_paraffin`: acyclic, saturated, every carbon degree <= 2.
- `iso_paraffin`: acyclic, saturated, 1 to 3 tertiary carbons, no quaternary carbon.
- `cycloalkane`: saturated, 1 or 2 rings, every ring size 5 or 6, `CalcNumBridgeheadAtoms == 0`, `CalcNumSpiroAtoms == 0`.
- `aromatic_mono`: exactly one aromatic ring; any second ring must be saturated (tetralin and indane allowed).
- `aromatic_di`: exactly two aromatic rings sharing a bond (naphthalenes).
- `alkene_alkyne`: acyclic, exactly one C=C, no triple bond.

`RealisticReps` considers only registry compounds that are usable, have valid pure eta(T)/rho(T) in the provider, and pass `is_realistic`. Within a bin, pick the compound with the lowest InChIKey (deterministic).

`reconstruct` algorithm: class masses from vol % with `RHO_CLASS` exactly as `fuel_hindcast.bulk_items` does; carbon distribution from `carbon_distribution(row)`; saturates split by `sat_split`; monoaromatic mass = aromatic mass x (1 - D/A) spread over the carbon distribution; diaromatic mass to the `aromatic_di` representative nearest C11. Any (kind, n) with `reps.get(...) is None` adds to fallback mass. Return an empty frame and fallback 1.0 when the row has fewer than 4 D86 points or no finite aromatics value.

`nominal_sat_split`: mean over the GC jet fuels returned by `fuelblend.data.loaders_jet.included_jet_fuel_ids` of their n / iso / cyclo mass shares renormalised to the saturates. Assert each fuel's class shares sum to 1 +/- 0.02 after the loader's own cleaning; skip and log any that do not. These fuels are not CRC targets, so this is a prior, not a fit.

- [ ] **Step 1: Confirm the descriptor names available per CRC fuel**

Run: `.venv/Scripts/python.exe -c "import pandas as pd; d=pd.read_csv('data/real_fuels/v4/jet_properties_v4.csv'); print(d[d.fuel_id.str.startswith('CRC647')&(d.role=='descriptor')].property.value_counts())"`
Expected: a list including an aromatics vol % property and distillation points. Map them to the `row` keys `aromatics_vol_pct`, `diaromatics_vol_pct`, `olefins_vol_pct`, `d86_ibp`, `d86_t5` ... `d86_fbp` in `rebuild_crc_reconstruction.py`. If CRC fuels have no descriptor rows here, read the same fields from the `fuels.csv` that `extras/fuelblend_data/scripts/build_real_fuels_dataset.py` uses (`FUELS_CSV_PATH`) and record that in the report.

- [ ] **Step 2: Write the failing tests**

```python
from fuelblend.data.reconstruction_r2 import is_realistic

def test_cage_and_polycyclic_rejected():
    assert not is_realistic("C1C[C@@H]2[C@H]3CC[C@H](C3)[C@@H]2C1", "cycloalkane")   # JP-10
    assert not is_realistic("C1C2CC3C4CC5CC(C14)C(C2)C3C5", "cycloalkane")          # diamantane
    assert not is_realistic("c1ccc2cc3ccccc3cc2c1", "aromatic_mono")                # anthracene
    assert not is_realistic("c1ccc(-c2ccccc2)cc1", "aromatic_mono")                 # biphenyl
    assert not is_realistic("CC(CC(C)(C)C)CC(C)(C)C", "iso_paraffin")               # quaternary carbons

def test_realistic_accepted():
    assert is_realistic("CCCCC1CCCCC1", "cycloalkane")
    assert is_realistic("C1CCC2CCCCC2C1", "cycloalkane")                            # decalin
    assert is_realistic("CCCCCc1ccccc1", "aromatic_mono")
    assert is_realistic("c1ccc2c(c1)CCCC2", "aromatic_mono")                        # tetralin
    assert is_realistic("Cc1cccc2ccccc12", "aromatic_di")
    assert is_realistic("CCCCCCCC(C)C", "iso_paraffin")

def test_deuterated_rejected_everywhere():
    d = "[2H]C([2H])([2H])C([2H])([2H])C([2H])([2H])[2H]"
    assert not is_realistic(d, "n_paraffin") and not is_realistic(d, "iso_paraffin")
```

Also write, using a stub reps object whose `get` returns `None` for cycloalkanes:
- `test_missing_bin_goes_to_fallback`: fallback mass equals the cycloalkane mass share within 1e-9, and no returned row has `kind == "cycloalkane"`.
- `test_short_d86_is_skipped`: a row with 3 D86 points returns an empty frame and fallback 1.0.
- `test_fractions_sum`: for a normal row, `frame.fraction.sum() + fallback == 1` within 1e-9.
- `test_reconstruction_reads_no_target`: `inspect.getsource(reconstruction_r2)` contains none of `density`, `flash`, `freez`, `visc`, `heat_of_combustion` outside the `RHO_CLASS` import line.

- [ ] **Step 3: Run the tests and confirm they fail**

Run: `.venv/Scripts/python.exe -m pytest tests/test_v4_recon_r2.py -q`
Expected: collection error, module not found.

- [ ] **Step 4: Implement `reconstruction_r2.py` to the interfaces above**

- [ ] **Step 5: Run the tests and confirm they pass**

Run: `.venv/Scripts/python.exe -m pytest tests/test_v4_recon_r2.py -q`
Expected: all pass.

- [ ] **Step 6: Implement and run `rebuild_crc_reconstruction.py`**

It writes `compositions_crc647_r2.csv` for the 57 CRC fuels at the nominal split, and `results/v4/R2/reconstruction_report.md` containing: the nominal split and the GC fuels it came from; per-kind mean mass share; the representative chosen for every bin; mean and max fallback mass.

Run: `.venv/Scripts/python.exe scripts/v4/rebuild_crc_reconstruction.py`
Expected: 57 fuels written. **Stop condition:** if mean fallback mass exceeds 0.20, stop and report which bins have no realistic representative. Do not relax the rules; the fix then needs pure-property estimates for missing compounds and is a separate decision for the user.

- [ ] **Step 7: Commit**

```bash
git add src/fuelblend/data/reconstruction_r2.py scripts/v4/rebuild_crc_reconstruction.py tests/test_v4_recon_r2.py data/real_fuels/v4/compositions_crc647_r2.csv results/v4/R2/reconstruction_report.md
git commit -m "v4 review fix 1: realistic representatives and r2 CRC reconstruction"
```

---

### Task 2: Re-score density, freezing, flash and viscosity on r2, with split sensitivity

**Files:**
- Modify: `scripts/v4/eval_bulk_props_v4.py`, `scripts/v4/eval_flash_v4.py`, `scripts/v4/build_viscosity_leaderboard_v4.py` (argument parsing and the `COMP` / `OUT` paths only)
- Create: `scripts/v4/recon_sensitivity.py`
- Test: `tests/test_v4_recon_r2.py` (add)
- Output: `results/v4/R2/{B,F,V}/`, `results/v4/R2/sat_split_sensitivity.csv`, `results/v4/R2/README.md`

**Interfaces:**
- Consumes: `compositions_crc647_r2.csv`; `reconstruct`, `RealisticReps`, `nominal_sat_split` from Task 1.
- Produces: each of the three scripts accepts `--compositions PATH` (default: the current v4 file) and `--out DIR` (default: the current output folder). With defaults, outputs are byte-identical to today.
- `recon_sensitivity.py` writes one row per (fuel, draw, property) with columns `fuel_id, split_role, draw, sat_n, sat_iso, sat_cyclo, property, T_C, predicted, measured`.

Sensitivity design: 32 Dirichlet draws, `np.random.default_rng(20261111).dirichlet(np.array(nominal) * 50)`, dev fuels only. Properties: density at 15.56 C (ideal volume), ideal-SLE freezing with mole fractions, raw Refutas at -20 C, ideal Le Chatelier flash.

Pre-declared reading, written into the README before looking at locked rows:
- Density: if dev bias at the nominal split is within +/- 3 kg/m3, the v4 bias is attributed to the old reconstruction. Otherwise report the remaining bias and the spread across draws (5th to 95th percentile of the per-fuel prediction).
- Freezing: report the ideal-SLE bias at the nominal split next to v4's +14.69 C (mole-fraction variant) and whether any model now beats the dev-mean constant.
- Viscosity: report raw Refutas dev bias next to v4's -0.301 ln. If |bias| drops below 0.10 ln, the recalibration was correcting the reconstruction.
- Flash: report MAE with a dev-fit constant offset (leave-one-fuel-out on dev, dev-fit applied to locked) next to the dev-mean constant.

- [ ] **Step 1: Write the failing test**

`test_default_paths_unchanged(tmp_path)`: run `eval_bulk_props_v4.py --out <tmp_path>` with the default `--compositions`, then assert `filecmp.cmp(tmp_path / "density_summary_v4.csv", "results/v4/B/density_summary_v4.csv", shallow=False)`.

- [ ] **Step 2: Run it and confirm it fails** (`unrecognized arguments: --out`)

Run: `.venv/Scripts/python.exe -m pytest tests/test_v4_recon_r2.py::test_default_paths_unchanged -q`

- [ ] **Step 3: Add `--compositions` and `--out` to the three scripts**

In `eval_bulk_props_v4.py`, `ALLOWED_COMPOSITION_PROVENANCE` gains `"reconstructed_r2"`. Drop `LUMP_NO_REALISTIC_REP` rows the same way `is_lumped` rows are dropped today.

- [ ] **Step 4: Run the test and confirm it passes**

- [ ] **Step 5: Run the three evaluators on r2**

```bash
.venv/Scripts/python.exe scripts/v4/eval_bulk_props_v4.py --compositions data/real_fuels/v4/compositions_crc647_r2.csv --out results/v4/R2/B
.venv/Scripts/python.exe scripts/v4/eval_flash_v4.py --compositions data/real_fuels/v4/compositions_crc647_r2.csv --out results/v4/R2/F
.venv/Scripts/python.exe scripts/v4/build_viscosity_leaderboard_v4.py --compositions data/real_fuels/v4/compositions_crc647_r2.csv --out results/v4/R2/V
```

Expected: each finishes and writes the same file set as the v4 folder.

- [ ] **Step 6: Implement and run `recon_sensitivity.py`; write `results/v4/R2/README.md`**

The README has one v4-versus-r2 table per property (dev and locked in separate rows), the pre-declared readings above with their outcome, and this disclosure sentence: "The 12 locked fuels were already scored under the v4 reconstruction; r2 design choices were fixed on dev fuels and the GC-fuel prior only."

- [ ] **Step 7: Commit**

```bash
git add scripts/v4 tests/test_v4_recon_r2.py results/v4/R2
git commit -m "v4 review fix 2: re-score CRC properties on r2 reconstruction with saturate-split sensitivity"
```

---

### Task 3: Known-composition lab residuals and margins

**Files:**
- Create: `scripts/v4/eval_lab_known_composition.py`
- Modify: `scripts/v4/build_uncertainty_v4.py` (add `lab_known()` and call it from `main()`)
- Test: `tests/test_v4_lab_known.py`
- Output: `splits/v4/lab_systems_splits.json`, `results/v4/L/{density_lab_residuals.csv, viscosity_lab_residuals.csv, README.md}`, new keys in `configs/uncertainty_v4.yaml`

**Interfaces:**
- Consumes: `eval_bulk_props_v4.PureDensity` (tier 1 only); `build_viscosity_leaderboard_v4.build_lab` (its `resid` column is ln measured minus ln raw Refutas; negate to get `predicted - measured`); `fuelblend.uq.conformal`.
- Produces:
  - `t_bin(T_K: float) -> str` in the new script: `"T_lt_263K"` if `T_K < 263.15`, `"T_263_283K"` if `263.15 <= T_K < 283.15`, else `"T_ge_283K"`.
  - Residual CSV columns: `system_id, split_role, T_K, t_bin, n_comp, predicted, measured, resid` (`resid` in kg/m3 for density, ln-ratio for viscosity).
  - yaml keys: `density_kg_m3.lab_mixtures.ideal_volume.over` and `.under` (one-sided), and `viscosity_kin_cSt.lab_mixtures.raw_refutas.T_lt_263K`, `.T_263_283K`, `.T_ge_283K` (two-sided, unit `ln-ratio`).

Data:
- Density: `data/processed/v2/thermoml_clean_hc.csv`, rows with `n_comp >= 2`, `prop` containing `Mass density`, `prop_phase == "Liquid"`, `P_kPa` missing or <= 110, every component a hydrocarbon with a tier-1 pure density valid at that T. System id = md5[:12] of the sorted InChIKey set, as in `build_unified.py`. Convert `comp_values` to mass fractions according to `comp_type`.
- Viscosity: what `build_lab` already selects, without its `LAB_T_RANGE` filter.

Lab split: all system ids from both datasets, sorted; 20 % locked, `np.random.default_rng(20261111)`; written to `splits/v4/lab_systems_splits.json` **before** any residual is computed. Margins use dev systems only; locked systems give a counts-only coverage check, like the existing `locked_check` entries.

- [ ] **Step 1: Write the failing tests**

```python
def test_t_bin_edges():          # t_bin loaded by path with importlib, as other v4 tests load scripts
    assert t_bin(263.14) == "T_lt_263K"
    assert t_bin(263.15) == "T_263_283K"
    assert t_bin(283.15) == "T_ge_283K"

def test_cold_viscosity_bin_has_no_substitute():
    m = yaml.safe_load(open("configs/uncertainty_v4.yaml"))["margins"]["viscosity_kin_cSt.lab_mixtures.raw_refutas.T_lt_263K"]
    if m["n_groups"] < 9:
        assert m["value"] == "inf" and "applied_margin" not in m
```

Also write:
- `test_existing_margin_keys_unchanged`: load the yaml from `git show HEAD:configs/uncertainty_v4.yaml` and the regenerated file; every key present at HEAD has an identical dict now.
- `test_system_without_pure_density_is_excluded`: a two-system input where one system has a component absent from `PureDensity.t1` yields rows for one system only and `n_excluded_no_pure_density == 1`.
- `test_lab_locked_systems_not_in_margin_rows`: for each new yaml key, the residual rows selected by its `recompute.where` filter contain no system with `split_role == "locked"`.

- [ ] **Step 2: Run and confirm they fail**

Run: `.venv/Scripts/python.exe -m pytest tests/test_v4_lab_known.py -q`

- [ ] **Step 3: Implement `eval_lab_known_composition.py`, write the split, run it**

Run: `.venv/Scripts/python.exe scripts/v4/eval_lab_known_composition.py`
Expected: both residual CSVs written; the README reports, per dataset and bin, records, systems, MAE, bias, and the count excluded for missing pure density.

- [ ] **Step 4: Add `lab_known()` to `build_uncertainty_v4.py` using the existing `add(...)` helper, then regenerate**

Run: `.venv/Scripts/python.exe scripts/v4/build_uncertainty_v4.py`
Expected: five new keys; existing keys unchanged.

- [ ] **Step 5: Run the tests and confirm they pass**

Run: `.venv/Scripts/python.exe -m pytest tests/test_v4_lab_known.py tests/test_v4_conformal.py -q`

- [ ] **Step 6: Commit**

```bash
git add scripts/v4/eval_lab_known_composition.py scripts/v4/build_uncertainty_v4.py tests/test_v4_lab_known.py splits/v4/lab_systems_splits.json results/v4/L results/v4/U configs/uncertainty_v4.yaml
git commit -m "v4 review fix 3: known-composition lab residuals and conformal margins"
```

---

### Task 4: Design engine and API use raw Refutas and lab margins for known compositions

**Files:**
- Modify: `src/fuelblend/design/engine_v4.py` (`PropertyModel`, `_combined_margin`, `evaluate_constraint`)
- Modify: `src/fuelblend/design/models_v4.py` (`build_models`: density and viscosity entries)
- Modify: `src/fuelblend/api/predict_blend.py` (`_predict_viscosity`, margin key table, `predict` signature)
- Modify: `scripts/v4/run_designs_v4.py`, `scripts/v4/verify_designs_v4.py` (add `--out`, default unchanged)
- Test: `tests/test_v4_design.py`, add `tests/test_v4_api_known_composition.py`
- Output: `results/v4/D2/`

**Interfaces:**
- Consumes: the five yaml keys from Task 3.
- Produces:
  - `viscosity_t_bin(T_C: float) -> str` in `engine_v4.py`, same bins as Task 3's `t_bin` (the engine must not import from `scripts/`).
  - `PropertyModel.margin_keys_by_side: Optional[dict[str, Sequence[str]]]` with keys `"over"` and `"under"`. When set, the lower limit uses the `"over"` margin and the upper limit uses the `"under"` margin; `margin_keys` is ignored.
  - `PropertyModel.margin_keys_fn: Optional[Callable[[Optional[float]], Sequence[str]]]`, called with the constraint temperature; when set it replaces `margin_keys`.
  - `predict(..., composition_kind: str = "known")` in the API; allowed values `"known"` and `"reconstructed"`.

Decisions:
- Design viscosity primary becomes `refutas_raw` (`f["visc_raw"]`) with `margin_keys_fn=lambda T_C: [f"viscosity_kin_cSt.lab_mixtures.raw_refutas.{viscosity_t_bin(T_C)}"]`, `T_valid_C=None`, `certified_kinds=None` (lab mixtures are not jet-specific; the class-extrapolation check still applies). `refutas_recalibrated` and the ridge stay as alternates in the robustness table only.
- Design density uses `margin_keys_by_side={"over": ["density_kg_m3.lab_mixtures.ideal_volume.over"], "under": ["density_kg_m3.lab_mixtures.ideal_volume.under"]}`, `certified_kinds=None`.
- API: `composition_kind="known"` always uses raw Refutas, for any number of components, with the lab margin for the temperature bin. `composition_kind="reconstructed"` keeps today's recalibrated path and CRC margin, and still requires `fuel_type="jet"`. The single-component special case at `predict_blend.py:583` is removed.
- If a side's margin is missing or `inf`, the constraint is `uncertified` with the yaml reason, exactly as for `margin_keys` today.

- [ ] **Step 1: Write the failing tests**

```python
def test_binary_is_continuous_with_pure():
    pure = predict(["CCCCCCCCCC"], [1.0], T_K=253.15, composition_kind="known")["viscosity_dyn_mPas"]["value"]
    near = predict(["CCCCCCCCCC", "CCCCCCCCCCCC"], [0.999, 0.001], T_K=253.15, composition_kind="known")["viscosity_dyn_mPas"]["value"]
    assert abs(near / pure - 1) < 0.01

def test_known_composition_never_uses_recalibration():
    r = predict(["CCCCCCCCCC", "CCCCCCCCCCCC"], [0.5, 0.5], T_K=253.15, composition_kind="known")
    assert r["viscosity_dyn_mPas"]["method"] == "refutas_raw"

def test_bin_edge_in_engine():
    assert viscosity_t_bin(-20.0) == "T_lt_263K"
    assert viscosity_t_bin(-10.0) == "T_263_283K"
    assert viscosity_t_bin(10.0) == "T_ge_283K"
```

Match the `predict` call signature that existing API tests use; adjust argument names here to it. Also write:
- `test_inf_bin_is_uncertified_not_borrowed`: with a margins dict where `...T_lt_263K` has no `applied_margin` and `...T_ge_283K` has 0.05, `evaluate_constraint("kin_visc_m20C_cSt", {"max": 8.0}, ...)` returns `status == "uncertified"`, `applied_margin is None`, and a reason containing `T_lt_263K`.
- `test_side_specific_density_margin`: with over = 1.0 and under = 5.0, prediction 780 and limits 775..840, the slack equals `min((780 - 1.0 - 775) / 65, (840 - 785.0) / 65)`.

- [ ] **Step 2: Run and confirm they fail**

Run: `.venv/Scripts/python.exe -m pytest tests/test_v4_design.py tests/test_v4_api_known_composition.py -q`

- [ ] **Step 3: Implement the engine changes, then `models_v4.py`, then the API**

`tests/test_v4_design.py` has an AST test that restricts numeric literals in the engine. Implement `viscosity_t_bin` by comparing `T_C` with `-10.0` and `10.0`, and add those two literals to that test's allow-list with a comment naming the bins.

- [ ] **Step 4: Run the v4 tests**

Run: `.venv/Scripts/python.exe -m pytest tests -q -k "v4"`
Expected: all pass. Any existing test that pinned `refutas_recalibrated` as the design primary is updated to the new primary, with the old expectation kept for the alternate.

- [ ] **Step 5: Re-run the designs into a new folder**

Run: `.venv/Scripts/python.exe scripts/v4/run_designs_v4.py --out results/v4/D2`, then `.venv/Scripts/python.exe scripts/v4/verify_designs_v4.py --out results/v4/D2`.
Expected: `designs_v4.csv` with 13 rows. Write `results/v4/D2/README.md` with the category counts next to the v4 counts (5 / 0 / 0 / 8) and, for every constraint whose status changed, the old and new status and reason. Report the certified count as it comes out. Do not adjust pools, margins or rules to change it.

- [ ] **Step 6: Commit**

```bash
git add src/fuelblend/design src/fuelblend/api scripts/v4/run_designs_v4.py scripts/v4/verify_designs_v4.py tests results/v4/D2
git commit -m "v4 review fix 4: known-composition viscosity and density path in design engine and API"
```

---

### Task 5: Record the results

**Files:**
- Modify: `docs/claims_register.md` (append rows), `docs/FULL_REPORT_v4.md` (append section 13), `docs/PROJECT_SUMMARY_v4.md` (section 7 items 3, 4, 7), `scripts/v4/reproduce_all.py` (new steps and manifest paths), `SHA256SUMS_v4.txt` (regenerated)

**Interfaces:**
- Consumes: `results/v4/R2/README.md`, `results/v4/L/README.md`, `results/v4/D2/README.md`.

- [ ] **Step 1: Append claims rows, each citing its CSV**
  - `CLM-38` r2 reconstruction: representative rules, nominal split, fallback mass, and the density / freezing / viscosity / flash outcomes against the pre-declared readings of Task 2.
  - `CLM-39` known-composition lab margins: the five keys, group counts, which are `inf`.
  - `CLM-40` design re-run: D2 category counts versus v4, with the statement that -20 C viscosity is uncertified where the cold bin is `inf`.
  - Update `CLM-25` (density bias "cause unknown") and `CLM-33` to point at CLM-38 and CLM-40. Old text stays; add a dated note, as the register's policy requires.

- [ ] **Step 2: Add the new scripts to `reproduce_all.py` in dependency order** (rebuild reconstruction, three r2 evaluations, sensitivity, lab evaluation, uncertainty build, D2 designs) and add the new output folders to the manifest globs.

- [ ] **Step 3: Run the pipeline twice and compare manifests**

Run: `.venv/Scripts/python.exe scripts/v4/reproduce_all.py --skip-tests --manifest-out SHA256SUMS_v4.txt`, copy the file aside, run again, and `diff` the two.
Expected: no differences.

- [ ] **Step 4: Commit**

```bash
git add docs scripts/v4/reproduce_all.py SHA256SUMS_v4.txt
git commit -m "v4 review fixes: claims register, report section 13, reproduction steps"
```

---

## Out of scope (separate plans)

- Flash inputs: the four `webbook` Antoine sets (mesitylene and three alkylbenzenes) and the Liaw formulation.
- The in-sample Coutinho margin (`freezing_point_C.crc647.coutinho_library`), joint-coverage correction, and bulk-descriptor baselines.
- New data: measured-composition fuels, low-temperature lab viscosity (needed to make the cold viscosity bin finite), and extending ideal SLE beyond 8 systems.
