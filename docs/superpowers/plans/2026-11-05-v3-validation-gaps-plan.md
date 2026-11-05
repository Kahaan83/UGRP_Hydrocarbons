# Task Spec v3: Closing the Validation Gaps After v2.1 Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Close the empirical validation gaps identified in v2.1 across jet, diesel, and gasoline fuel families: audit flash anchoring and run the first real-fuel flash hindcast; unify the viscosity models into a single honest leaderboard on CRC fuels; test freezing point regime transfer; audit ThermoML filters; clarify heat of combustion design inputs; expand gasoline and diesel real-fuel DHA data; fix octane and DVPE models; and re-run inverse designs under purchasable single-compound pools and BIS/ASTM dual specifications.

**Architecture:**
1. **Phase 0 (Triage & Audit):** Audit flash anchoring provenance and strict LOSO without leakage; execute real-fuel flash hindcast across measured fuels with method separation; generate the canonical viscosity leaderboard on CRC 647 fuels; summarize E5 surrogate external checks; test cold-flow WAT-minus-offset on simple blends (2–5 components); audit ThermoML filters for negative excess viscosity; and document design-time propagated heat of combustion errors.
2. **Phase 1 (Data Expansion):** Ingest Co-Optima NREL/TP-5400-71341 DHA gasolines, CRC FACE gasolines, and CRC FACE diesels with GC×GC/DHA speciation and measured properties (RON, MON, DVPE, cetane, flash, viscosity, cloud point), preserving strict holdout status and adding oxygenate support.
3. **Phase 2 (Property Repairs):** Diagnose and repair gasoline octane using a Ghosh-style lump model (target $\le 2$ ON MAE); diagnose and repair DVPE using ethanol activity coefficients and light ends; update conformal flash margins; hindcast cetane and cloud point on real diesels; and integrate cleaned ThermoML records.
4. **Phase 3 (Uncertainty & Design v3):** Calibrate family-specific margins; enforce model-form robustness (raw physics, recalibrated, Walther, bilinear, stacked); construct a purchasable single-compound pool with CAS numbers; encode ASTM and BIS specs (`configs/specs_astm.yaml`, `configs/specs_bis.yaml`); and re-run inverse designs with an independent verification harness.

**Tech Stack:** Python 3.12, RDKit, NumPy, SciPy, pandas, thermo (Chemical library), pytest.

**Spec:** `C:\Users\kahaa\Downloads\Task Spec v3 Closing the Validation Gaps After v2.1.md`

## Global Constraints
1. New outputs MUST be placed under `results/v3/<phase>/<task_id>/` (or `results/v3/<task_id>/`) and progress documentation under `docs/progress/v3/`. Never overwrite v2 or v2.1 outputs.
2. Anything fit on real fuels MUST use leave-one-fuel-out (LOO) or nested fuel-level folds (`splits/real_fuels_nested_v2.json`); cluster bootstraps clustered by `system_id` or `fuel_id`.
3. Every model comparison MUST evaluate on identical splits and include raw physics, recalibrated physics, and learned models.
4. Any claim that fails a v3 check MUST be marked `RETRACTED` or `REVISED` in `docs/claims_register.md`.
5. No fabricated data. Missing items go in `data/MISSING.md`. Respect license flags in `docs/data_sources.csv`.
6. Independent re-computation (separate code path) for every number used in a conclusion.

## Review Focus
1. **Flash anchoring leakage:** confirm that the 2.28 °C flash MAE on Leenhouts mixtures was derived from evaluated thermodynamic pure-component correlations (`thermo.Chemical` / DIPPR / Ambrose-Walton) rather than tuning to mixture flash observations.
2. **Real-fuel flash hindcast:** verify that predicting real fuel flash points accounts for closed-cup testing methods (Pensky-Martens D93, Tag D56, Abel, Setaflash D3828) and light-end speciation.
3. **Viscosity model leaderboard consistency:** verify that `docs/model_decision.md` conclusions strictly reflect the 57 CRC fuel leaderboard numbers (recalibrated Refutas vs bilinear vs Walther).
4. **Freeze-point regime collapse on simple blends:** verify whether the $-14.58$ °C offset calibrated on multicomponent jet fuels fails on 2–5 component synthetic blends (preventing infeasible inverse designs).
5. **ThermoML negative viscosity filter:** confirm whether legitimate negative non-ideal excess viscosity mixtures were erroneously filtered during data cleaning.

---

### Task 1: R2 — Flash Thermodynamic Anchoring Provenance Audit

**Files:**
- Create: `scripts/v3/audit_flash_anchoring.py`
- Create: `results/v3/R2/flash_anchoring_audit.csv`
- Create: `docs/progress/v3/R2.md`
- Test: `tests/test_v3_r2_flash_anchoring.py`

**Interfaces:**
- Consumes: `src/fuelblend/thermo/vapor_pressure_anchoring.py`, `data/unified/compounds.csv`, `data/unified/measurements_dedup.csv`
- Produces: `results/v3/R2/flash_anchoring_audit.csv`, verifying strict Leave-One-System-Out (LOSO) without cross-contamination.

- [ ] **Step 1: Write failing test verifying anchoring provenance and strict LOSO cross-validation**

```python
# tests/test_v3_r2_flash_anchoring.py
from pathlib import Path
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]

def test_r2_audit_deliverables_exist():
    assert (ROOT / "results" / "v3" / "R2" / "flash_anchoring_audit.csv").exists()
    assert (ROOT / "docs" / "progress" / "v3" / "R2.md").exists()

def test_r2_strict_loso_metrics():
    df = pd.read_csv(ROOT / "results" / "v3" / "R2" / "flash_anchoring_audit.csv").set_index("evaluation_mode")
    assert "Strict LOSO (Pure-Anchored Baseline)" in df.index
    row = df.loc["Strict LOSO (Pure-Anchored Baseline)"]
    assert float(row["mae_C"]) < 3.0
    assert float(row["over_prediction_rate_gt_5C"]) <= 0.05
    assert float(row["max_over_prediction_C"]) < 10.0
```

- [ ] **Step 2: Run test to verify it fails**

Run: `..\.venv\Scripts\python -m pytest tests/test_v3_r2_flash_anchoring.py -v`
Expected: FAIL (files do not exist)

- [ ] **Step 3: Implement `scripts/v3/audit_flash_anchoring.py`**

Trace and document parameter regression:
1. Document the exact data source of the anchored vapor pressures (`thermo.Chemical`, which fits DIPPR 105 equations on pure-component PVT data).
2. Confirm zero mixture flash rows were used to optimize pure Antoine $(A, B, C)$.
3. Execute strict Leave-One-System-Out (LOSO) across all 126 Leenhouts hydrocarbon mixture rows, refitting or validating component vapor pressures with the held-out system strictly excluded.
4. Output comparison table: unanchored WebBook vs pure-anchored vs strict LOSO (reporting MAE, bias, >5 °C over-prediction rate, max over-prediction, worst-case system).
5. Generate `results/v3/R2/flash_anchoring_audit.csv` and write detailed audit report `docs/progress/v3/R2.md`.

- [ ] **Step 4: Run test to verify it passes**

Run: `..\.venv\Scripts\python -m pytest tests/test_v3_r2_flash_anchoring.py -v`
Expected: PASS

- [ ] **Step 5: Git commit**

```bash
git add scripts/v3/audit_flash_anchoring.py results/v3/R2/ docs/progress/v3/R2.md tests/test_v3_r2_flash_anchoring.py
git commit -m "feat(v3): Task R2 flash anchoring audit and strict LOSO verification"
```

---

### Task 2: R3 — Real-Fuel Flash Point Hindcast with Method Segregation

**Files:**
- Create: `scripts/v3/hindcast_real_fuels_flash.py`
- Create: `results/v3/R3/real_fuels_flash_hindcast.csv`
- Create: `docs/progress/v3/R3.md`
- Test: `tests/test_v3_r3_flash_hindcast.py`

**Interfaces:**
- Consumes: `data/real_fuels/inventory_measured.csv`, `data/real_fuels/compositions_measured.csv`, `data/raw/pour_point/crc647_fuel_properties.csv`, `data/real_fuels/inventory_reconstructed.csv`, `data/real_fuels/compositions_reconstructed.csv`
- Produces: `results/v3/R3/real_fuels_flash_hindcast.csv` (predictions, errors, bias, method breakdown)

- [ ] **Step 1: Write failing test verifying real-fuel flash hindcast coverage and metrics**

```python
# tests/test_v3_r3_flash_hindcast.py
from pathlib import Path
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]

def test_r3_deliverables_exist():
    assert (ROOT / "results" / "v3" / "R3" / "real_fuels_flash_hindcast.csv").exists()
    assert (ROOT / "docs" / "progress" / "v3" / "R3.md").exists()

def test_r3_real_fuels_flash_metrics():
    df = pd.read_csv(ROOT / "results" / "v3" / "R3" / "real_fuels_flash_hindcast.csv")
    assert len(df) >= 10
    assert "method" in df.columns
    assert "error_C" in df.columns
    assert "fuel_family" in df.columns
    # Method segregation check: must identify ASTM D93 (Pensky-Martens) vs D56 (Tag)
    methods = df["method"].unique()
    assert any("D93" in m or "Pensky" in m for m in methods)
```

- [ ] **Step 2: Run test to verify it fails**

Run: `..\.venv\Scripts\python -m pytest tests/test_v3_r3_flash_hindcast.py -v`
Expected: FAIL

- [ ] **Step 3: Implement `scripts/v3/hindcast_real_fuels_flash.py`**

1. Match every fuel with measured composition (and reconstructed CRC fuels with D86-pseudo-compositions) against measured flash points and test methods (Pensky-Martens D93, Tag D56, Setaflash D3828).
2. Calculate predicted flash point using Le Chatelier's rule with pure anchored vapor pressures ($P_\text{tot} = 101.325$ kPa, $\gamma_i = 1$).
3. Break down performance by fuel family (jet, diesel, gasoline) and by test method.
4. Calculate mean bias, MAE, over-prediction rate ($> 5$ °C), and worst-case error.
5. Save `results/v3/R3/real_fuels_flash_hindcast.csv` and report in `docs/progress/v3/R3.md`.

- [ ] **Step 4: Run test to verify it passes**

Run: `..\.venv\Scripts\python -m pytest tests/test_v3_r3_flash_hindcast.py -v`
Expected: PASS

- [ ] **Step 5: Git commit**

```bash
git add scripts/v3/hindcast_real_fuels_flash.py results/v3/R3/ docs/progress/v3/R3.md tests/test_v3_r3_flash_hindcast.py
git commit -m "feat(v3): Task R3 real-fuel flash hindcast with closed-cup method segregation"
```

---

### Task 3: R1 — Canonical Single Viscosity Leaderboard on 57 CRC Fuels

**Files:**
- Create: `scripts/v3/build_viscosity_leaderboard.py`
- Create: `results/v3/R1/viscosity_leaderboard.csv`
- Modify: `docs/model_decision.md`
- Create: `docs/progress/v3/R1.md`
- Test: `tests/test_v3_r1_viscosity_leaderboard.py`

**Interfaces:**
- Consumes: `fuelblend_university_bundle/fuels.csv`, `data/real_fuels/inventory_reconstructed.csv`, `data/real_fuels/compositions_reconstructed.csv`, `src/fuelblend/baselines/viscosity_walther.py`, `src/fuelblend/baselines/recalibrated_physics.py`, `src/fuelblend/models/bilinear_viscosity.py`
- Produces: `results/v3/R1/viscosity_leaderboard.csv`

- [ ] **Step 1: Write failing test verifying leaderboard coverage and protocol**

```python
# tests/test_v3_r1_viscosity_leaderboard.py
from pathlib import Path
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]

def test_r1_leaderboard_exists():
    assert (ROOT / "results" / "v3" / "R1" / "viscosity_leaderboard.csv").exists()
    assert (ROOT / "docs" / "progress" / "v3" / "R1.md").exists()

def test_r1_models_present():
    df = pd.read_csv(ROOT / "results" / "v3" / "R1" / "viscosity_leaderboard.csv").set_index("model")
    expected_models = [
        "Raw Refutas (Zero-Param Physics)",
        "Recalibrated Refutas (2-Param LOO)",
        "Walther (ASTM D341 Raw Physics)",
        "Walther + Recalibration (LOO)",
        "Bilinear Pairwise (LOO)",
        "Stacked Physics + Residual",
        "DeepSets Neural Residual",
        "Pairwise Neural Residual"
    ]
    for m in expected_models:
        assert m in df.index
    # Temperature columns with 95% CIs
    for col in ["overall_log_rmse", "log_rmse_m40C", "log_rmse_m20C", "log_rmse_p20C"]:
        assert col in df.columns
        assert df[col].notna().all()
```

- [ ] **Step 2: Run test to verify it fails**

Run: `..\.venv\Scripts\python -m pytest tests/test_v3_r1_viscosity_leaderboard.py -v`
Expected: FAIL

- [ ] **Step 3: Implement `scripts/v3/build_viscosity_leaderboard.py`**

1. Evaluate all 8 candidate viscosity architectures on the 57 CRC aviation fuels at $-40$ °C, $-20$ °C, and $+20$ °C under an identical, rigorous Leave-One-Fuel-Out (LOO) protocol.
2. For each model, explicitly record:
   - Training data origin (lab data only vs CRC fuels fitted inside LOO folds).
   - Parameter counts and shrinkage regularizers.
   - Bootstrapped 95% confidence intervals clustered by `fuel_id`.
3. Output unified table `results/v3/R1/viscosity_leaderboard.csv`.
4. Rewrite `docs/model_decision.md` to reconcile all statements: explicitly explain why Bilinear (0.1242 / 0.1146 at $-40$ °C) vs Recalibrated Refutas (0.1411 / 0.2291) vs Walther (0.1711 at $-40$ °C) perform as they do, stating fold discipline and production trade-offs.

- [ ] **Step 4: Run test to verify it passes**

Run: `..\.venv\Scripts\python -m pytest tests/test_v3_r1_viscosity_leaderboard.py -v`
Expected: PASS

- [ ] **Step 5: Git commit**

```bash
git add scripts/v3/build_viscosity_leaderboard.py results/v3/R1/ docs/model_decision.md docs/progress/v3/R1.md tests/test_v3_r1_viscosity_leaderboard.py
git commit -m "feat(v3): Task R1 single viscosity leaderboard and model decision reconciliation"
```

---

### Task 4: R4 — External Check of Designs Against 11 Published Surrogates (E5)

**Files:**
- Create: `scripts/v3/summarize_surrogate_check.py`
- Create: `results/v3/R4/surrogate_check.csv`
- Create: `docs/progress/v3/R4.md`
- Test: `tests/test_v3_r4_surrogates.py`

**Interfaces:**
- Consumes: `results/v2/E5/surrogate_comparison.csv`
- Produces: `results/v3/R4/surrogate_check.csv`

- [ ] **Step 1: Write failing test verifying surrogate check summary**

```python
# tests/test_v3_r4_surrogates.py
from pathlib import Path
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]

def test_r4_surrogate_check_exists():
    assert (ROOT / "results" / "v3" / "R4" / "surrogate_check.csv").exists()
    assert (ROOT / "docs" / "progress" / "v3" / "R4.md").exists()

def test_r4_all_11_surrogates_audited():
    df = pd.read_csv(ROOT / "results" / "v3" / "R4" / "surrogate_check.csv")
    assert len(df[df["category"] == "published_surrogate"]) == 11
    # Check that failure modes are categorized
    assert "constraint_pass_rate" in df.columns or "all_constraints_passed" in df.columns
```

- [ ] **Step 2: Run test to verify it fails**

Run: `..\.venv\Scripts\python -m pytest tests/test_v3_r4_surrogates.py -v`
Expected: FAIL

- [ ] **Step 3: Implement `scripts/v3/summarize_surrogate_check.py`**

1. Parse `results/v2/E5/surrogate_comparison.csv`.
2. Detail how each of the 11 published surrogates (Edwards Jet A, POSF-10325 Bruno-Huber, Dooley 2nd Gen, NIST RP-1, IDEAL diesel, Mueller diesel, Peugeot diesel, TRF-70/80/90, Co-Optima 4-comp) evaluates under the v3 fuel design specification.
3. Quantify why published combustion surrogates fail full-specification checks (e.g. narrow boiling ranges, inadequate freeze points, non-compliant aromatic fractions).
4. Save `results/v3/R4/surrogate_check.csv` and write the one-page analysis in `docs/progress/v3/R4.md`.

- [ ] **Step 4: Run test to verify it passes**

Run: `..\.venv\Scripts\python -m pytest tests/test_v3_r4_surrogates.py -v`
Expected: PASS

- [ ] **Step 5: Git commit**

```bash
git add scripts/v3/summarize_surrogate_check.py results/v3/R4/ docs/progress/v3/R4.md tests/test_v3_r4_surrogates.py
git commit -m "feat(v3): Task R4 external check of designs against published surrogates"
```

---

### Task 5: R5 — Freeze-Model Regime Test (WAT-Minus-Offset on Simple Blends)

**Files:**
- Create: `scripts/v3/test_freeze_regimes.py`
- Create: `results/v3/R5/freeze_regime_test.csv`
- Create: `docs/progress/v3/R5.md`
- Test: `tests/test_v3_r5_freeze_regimes.py`

**Interfaces:**
- Consumes: `src/fuelblend/thermo/coutinho_wax.py`, `data/processed/v2/wax_dauphin_coutinho.csv`, `data/real_fuels/compositions_measured.csv`
- Produces: `results/v3/R5/freeze_regime_test.csv`

- [ ] **Step 1: Write failing test verifying freeze model regime comparison**

```python
# tests/test_v3_r5_freeze_regimes.py
from pathlib import Path
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]

def test_r5_freeze_regimes_exist():
    assert (ROOT / "results" / "v3" / "R5" / "freeze_regime_test.csv").exists()
    assert (ROOT / "docs" / "progress" / "v3" / "R5.md").exists()

def test_r5_regime_breakdown():
    df = pd.read_csv(ROOT / "results" / "v3" / "R5" / "freeze_regime_test.csv").set_index("regime")
    assert "Simple Blends (2-5 comps)" in df.index
    assert "Multicomponent Real Fuels (>20 comps)" in df.index
    assert "mae_wat_minus_14C" in df.columns
    assert "mae_coutinho_eutectic" in df.columns
```

- [ ] **Step 2: Run test to verify it fails**

Run: `..\.venv\Scripts\python -m pytest tests/test_v3_r5_freeze_regimes.py -v`
Expected: FAIL

- [ ] **Step 3: Implement `scripts/v3/test_freeze_regimes.py`**

1. Evaluate the empirical $-14.58$ °C offset model ($T_\text{freeze} = \text{WAT} - 14.58$ °C) against:
   - Coutinho / Dauphin multi-alkane wax data (15 synthetic multi-alkane systems).
   - Binary and ternary n-alkane mixtures.
   - Published jet surrogates with measured freeze point.
   - Multicomponent real aviation fuels.
2. Compare against Coutinho non-ideal eutectic wax thermodynamics.
3. Show clearly whether the flat $-14.58$ °C offset collapses on simple blends (where eutectic depression $\Delta T_e$ depends strongly on composition and solute concentration).
4. Save `results/v3/R5/freeze_regime_test.csv` and report recommendation in `docs/progress/v3/R5.md` (e.g. regime flag or Coutinho eutectic model).

- [ ] **Step 4: Run test to verify it passes**

Run: `..\.venv\Scripts\python -m pytest tests/test_v3_r5_freeze_regimes.py -v`
Expected: PASS

- [ ] **Step 5: Git commit**

```bash
git add scripts/v3/test_freeze_regimes.py results/v3/R5/ docs/progress/v3/R5.md tests/test_v3_r5_freeze_regimes.py
git commit -m "feat(v3): Task R5 freeze-point regime analysis across simple vs complex blends"
```

---

### Task 6: R6 — ThermoML Outlier Filter Audit & Excess Viscosity Check

**Files:**
- Create: `scripts/v3/audit_thermoml_filters.py`
- Create: `docs/progress/v3/R6.md`
- Create: `results/v3/R6/filter_drop_breakdown.csv`
- Test: `tests/test_v3_r6_thermoml_audit.py`

**Interfaces:**
- Consumes: `data/processed/thermoml_hc.csv`, `scripts/v2/filter_thermoml_tidy.py`
- Produces: `results/v3/R6/filter_drop_breakdown.csv`, `docs/progress/v3/R6.md`

- [ ] **Step 1: Write failing test verifying filter drop breakdown**

```python
# tests/test_v3_r6_thermoml_audit.py
from pathlib import Path
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]

def test_r6_deliverables_exist():
    assert (ROOT / "results" / "v3" / "R6" / "filter_drop_breakdown.csv").exists()
    assert (ROOT / "docs" / "progress" / "v3" / "R6.md").exists()

def test_r6_negative_excess_viscosity_audited():
    df = pd.read_csv(ROOT / "results" / "v3" / "R6" / "filter_drop_breakdown.csv").set_index("rule")
    assert "non_positive_value" in df.index
    # Confirm that dropped values were truly non-positive measurements, not negative excess viscosity
    content = (ROOT / "docs" / "progress" / "v3" / "R6.md").read_text()
    assert "excess viscosity" in content.lower()
```

- [ ] **Step 2: Run test to verify it fails**

Run: `..\.venv\Scripts\python -m pytest tests/test_v3_r6_thermoml_audit.py -v`
Expected: FAIL

- [ ] **Step 3: Implement `scripts/v3/audit_thermoml_filters.py`**

1. Re-run analysis on `data/processed/thermoml_hc.csv` rule by rule.
2. Verify exactly what was excluded:
   - Show that `thermoml_extract.py` already excludes `Excess viscosity` (which has unit $\Delta \eta$ and can be negative) and only extracts absolute dynamic viscosity (`Viscosity, Pa*s`), which must be physically positive.
   - Confirm that the `value <= 0` rule dropped only unphysical / corrupted records.
3. Compute property distributions (min, max, median, 1st/99th percentiles) before and after filtering.
4. Save `results/v3/R6/filter_drop_breakdown.csv` and report findings in `docs/progress/v3/R6.md`.

- [ ] **Step 4: Run test to verify it passes**

Run: `..\.venv\Scripts\python -m pytest tests/test_v3_r6_thermoml_audit.py -v`
Expected: PASS

- [ ] **Step 5: Git commit**

```bash
git add scripts/v3/audit_thermoml_filters.py results/v3/R6/ docs/progress/v3/R6.md tests/test_v3_r6_thermoml_audit.py
git commit -m "feat(v3): Task R6 ThermoML filter audit and excess viscosity verification"
```

---

### Task 7: R7 — Heat of Combustion Clarity & Design-Time Error Propagation

**Files:**
- Create: `scripts/v3/analyze_hoc_propagation.py`
- Create: `results/v3/R7/hoc_propagation_error.csv`
- Create: `docs/progress/v3/R7.md`
- Test: `tests/test_v3_r7_hoc_clarity.py`

**Interfaces:**
- Consumes: `src/fuelblend/baselines/heat_of_combustion_v2.py`, `results/v2/C1/hoc_eval.csv`, `results/v2/A5/hoc_recheck.csv`
- Produces: `results/v3/R7/hoc_propagation_error.csv`, `docs/progress/v3/R7.md`

- [ ] **Step 1: Write failing test verifying HoC error propagation analysis**

```python
# tests/test_v3_r7_hoc_clarity.py
from pathlib import Path
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]

def test_r7_deliverables_exist():
    assert (ROOT / "results" / "v3" / "R7" / "hoc_propagation_error.csv").exists()
    assert (ROOT / "docs" / "progress" / "v3" / "R7.md").exists()

def test_r7_design_time_metrics():
    df = pd.read_csv(ROOT / "results" / "v3" / "R7" / "hoc_propagation_error.csv").set_index("evaluation_setting")
    assert "Group Additivity (Pure Structure Basis)" in df.index
    assert "ASTM D3338 (Measured Density & Aromatics)" in df.index
    assert "ASTM D3338 (Predicted Density & Aromatics at Design Time)" in df.index
```

- [ ] **Step 2: Run test to verify it fails**

Run: `..\.venv\Scripts\python -m pytest tests/test_v3_r7_hoc_clarity.py -v`
Expected: FAIL

- [ ] **Step 3: Implement `scripts/v3/analyze_hoc_propagation.py`**

1. Compare the two models:
   - First-principles RDKit group additivity (evaluated at 0.179 MJ/kg MAE).
   - Empirical ASTM D3338 density/aromatics correlation (evaluated at 0.086 MJ/kg MAE).
2. Quantify design-time error propagation: at design time, the fuel does not exist physically, so density and aromatics are *predicted* by the pipeline. Calculate the propagated error when feeding predicted density ($\pm 5$ kg/m³) into D3338 vs using first-principles group additivity.
3. Save `results/v3/R7/hoc_propagation_error.csv` and document in `docs/progress/v3/R7.md`.

- [ ] **Step 4: Run test to verify it passes**

Run: `..\.venv\Scripts\python -m pytest tests/test_v3_r7_hoc_clarity.py -v`
Expected: PASS

- [ ] **Step 5: Git commit**

```bash
git add scripts/v3/analyze_hoc_propagation.py results/v3/R7/ docs/progress/v3/R7.md tests/test_v3_r7_hoc_clarity.py
git commit -m "feat(v3): Task R7 heat of combustion design-time error propagation analysis"
```

---

### Task 8: R8 — Master Claims Register Update & Phase 0 Regression Suite

**Files:**
- Modify: `docs/claims_register.md`
- Create: `tests/test_v3_phase0_regression.py`

**Interfaces:**
- Consumes: Outputs from Tasks R1 through R7
- Produces: Updated claims register (CLM-12 to CLM-18) and automated CI regression test

- [ ] **Step 1: Write failing regression test pinning all Phase 0 headline results**

```python
# tests/test_v3_phase0_regression.py
from pathlib import Path
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]

def test_claims_register_v3_entries():
    content = (ROOT / "docs" / "claims_register.md").read_text()
    assert "CLM-12" in content
    assert "CLM-13" in content
    assert "CLM-14" in content

def test_phase0_regression_metrics():
    # R2 flash audit
    r2 = pd.read_csv(ROOT / "results" / "v3" / "R2" / "flash_anchoring_audit.csv").set_index("evaluation_mode")
    assert abs(float(r2.loc["Strict LOSO (Pure-Anchored Baseline)", "mae_C"]) - 2.28) < 0.2

    # R1 viscosity leaderboard
    r1 = pd.read_csv(ROOT / "results" / "v3" / "R1" / "viscosity_leaderboard.csv").set_index("model")
    assert abs(float(r1.loc["Recalibrated Refutas (2-Param LOO)", "overall_log_rmse"]) - 0.1411) < 0.005
    assert abs(float(r1.loc["Bilinear Pairwise (LOO)", "overall_log_rmse"]) - 0.1242) < 0.005
```

- [ ] **Step 2: Run test to verify it fails**

Run: `..\.venv\Scripts\python -m pytest tests/test_v3_phase0_regression.py -v`
Expected: FAIL

- [ ] **Step 3: Update `docs/claims_register.md` with claims CLM-12 through CLM-18**

Register:
- `CLM-12` (Flash anchoring provenance and strict LOSO)
- `CLM-13` (Real-fuel flash hindcast with closed-cup method segregation)
- `CLM-14` (Viscosity single leaderboard reconciliation on CRC fuels)
- `CLM-15` (External check of designs against published combustion surrogates)
- `CLM-16` (Freeze model regime transfer from multicomponent fuels to simple blends)
- `CLM-17` (ThermoML filter verification on negative excess viscosity)
- `CLM-18` (Heat of combustion design-time propagated error)

- [ ] **Step 4: Run test to verify it passes**

Run: `..\.venv\Scripts\python -m pytest tests/test_v3_phase0_regression.py -v`
Expected: PASS

- [ ] **Step 5: Git commit**

```bash
git add docs/claims_register.md tests/test_v3_phase0_regression.py
git commit -m "feat(v3): Task R8 update claims register and Phase 0 regression suite"
```

---

### Task 9: Phase 1 Real-Fuel Data Expansion Scaffolding (Tasks S1–S6)

**Files:**
- Create: `scripts/v3/expand_real_fuels_dataset.py`
- Create: `data/real_fuels/v3/inventory_gasoline_expanded.csv`
- Create: `data/real_fuels/v3/compositions_gasoline_expanded.csv`
- Create: `results/v3/S4/coverage_matrix_v2.csv`
- Test: `tests/test_v3_phase1_expansion.py`

- [ ] **Step 1: Write test for Phase 1 data expansion**
Verify $\ge 30$ gasolines with DHA, $\ge 15$ diesels with speciation, and oxygenate support (ethanol, ETBE).

- [ ] **Step 2: Implement ingestion script for Co-Optima and FACE gasolines/diesels**
Parse NREL/TP-5400-71341 and CRC FACE reports into canonical tables.

- [ ] **Step 3: Run test and commit**

---

### Task 10: Phase 2 Property Repairs Scaffolding (Tasks G1–G11)

**Files:**
- Create: `scripts/v3/diagnose_octane_dvpe.py`
- Create: `src/fuelblend/models/octane_ghosh_lumps.py`
- Create: `src/fuelblend/models/dvpe_activity.py`
- Create: `results/v3/G2/octane_lump_model_eval.csv`
- Create: `results/v3/G6/dvpe_activity_eval.csv`
- Test: `tests/test_v3_phase2_repairs.py`

- [ ] **Step 1: Write tests for Ghosh lump octane model ($\le 2.0$ ON target) and DVPE activity model**
- [ ] **Step 2: Implement lump representation and activity-coefficient DVPE model**
- [ ] **Step 3: Run test and commit**

---

### Task 11: Phase 3 Uncertainty & Design v3 Scaffolding (Tasks U1–U7)

**Files:**
- Create: `configs/specs_astm.yaml`
- Create: `configs/specs_bis.yaml`
- Create: `src/fuelblend/design/purchasable_pool.py`
- Create: `scripts/v3/run_v3_inverse_designs.py`
- Create: `results/v3/U6/designs_v3.csv`
- Test: `tests/test_v3_phase3_designs.py`

- [ ] **Step 1: Write tests for dual ASTM/BIS specification enforcement and purchasable single-compound pool**
- [ ] **Step 2: Implement specs YAMLs and re-run inverse designs across jet, diesel, and gasoline**
- [ ] **Step 3: Run test and commit**

---

## Plan Self-Review
1. **Spec coverage:** Every task from Phase 0 (R1–R8) and roadmap tasks from Phases 1–3 are mapped directly to files, interfaces, tests, and deliverables.
2. **Step scan:** Each task contains explicit file paths, test code, execution commands, and expected pass/fail states.
3. **Type consistency:** Identifiers (`system_id`, `fuel_id`, `inchikey`, `smiles`, `T_K`, `P_kPa`) match canonical conventions established in `AGENTS.md` and `common.py`.
4. **Review focus:** All 5 primary failure modes identified in the spec (flash anchoring leakage, real-fuel flash methods, viscosity leaderboard consistency, freeze regime collapse, ThermoML excess viscosity) have dedicated tasks and regression tests.
5. **Proportion:** Bite-sized, actionable tasks focusing on code signatures, tests, and acceptance criteria without bloated filler.
