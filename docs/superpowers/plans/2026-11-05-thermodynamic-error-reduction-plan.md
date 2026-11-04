# Thermodynamic Error Reduction & Physical Accuracy Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Drive down prediction errors across the primary fuel properties (Flash Point, Freezing Point, Viscosity) by replacing empirical/extrapolated baselines with rigorous thermodynamic physics (low-T vapor pressure anchoring, Coutinho multi-component wax crystallization, ASTM D341 Walther viscosity curvature), then re-run inverse blend designs and generate the publication-ready ablation table.

**Architecture:** 
1. **Flash Point:** Anchor low-temperature vapor pressures (280–373 K) for C10–C20 alkanes and iso-alkanes using ThermoML evaluated experimental data and Cox/Růžička-Majer parameters, eliminating the +42 °C extrapolation bug in Le Chatelier predictions.
2. **Freezing Point:** Replace the flat -14.58 °C offset with a composition-dependent non-ideal solid-solution wax crystallization model (Coutinho SP-Wax formulation) fitted on multi-component n-alkane systems (`wax_dauphin_coutinho.csv`), predicting true eutectic depression.
3. **Viscosity:** Upgrade the cryogenic real-fuel model from linear $1000/T$ to ASTM D341 / Walther double-logarithmic curvature ($\log \log(\nu + 0.7) = A - B \log T$) combined with low-rank bilinear class interactions, eliminating the residual error at -40 °C.
4. **Inverse Design & Ablation:** Re-run Jet A and Diesel inverse designs using the upgraded thermodynamic engine, tracking errors version-by-version as a structured ablation story for the paper.

**Tech Stack:** Python 3.12, RDKit, NumPy, SciPy (optimization & root-finding), pandas, pytest.

**Spec:** User instructions following v2 master technical report (`docs/FULL_REPORT_v2.md`).

## Global Constraints
1. No overwriting v1 or existing v2 baseline files. All new deliverables go to `results/v2/improvement/` and progress documentation to `docs/progress/v2/thermodynamic_upgrades.md`.
2. Evaluate all model changes with strict cross-validation: Leave-One-System-Out (LOSO) for lab mixtures and Leave-One-Fuel-Out (LOO) for real fuels.
3. Every step includes automated verification tests asserting error reduction without regressions.
4. Tag new compound parameter sources explicitly in registry (`antoine_source = "thermoml_cox_anchored"`).

## Review Focus
- **C12+ Iso-alkane ambient vapor pressure:** ensure anchored Antoine/Cox curves do not diverge at sub-ambient temperatures (250–300 K).
- **Eutectic depression consistency:** ensure the Coutinho multi-solid wax model predicts $T_\text{freeze} \le \text{WAT}_\text{pure}$ for all mixtures (freezing point must never exceed the highest-melting pure component).
- **Walther viscosity low-temperature limit:** verify that Walther parameters maintain monotonicity ($\frac{d\nu}{dT} < 0$) down to -40 °C without mathematical singularity ($Z = \nu + 0.7 > 0$).
- **Inverse design constraint feasibility:** confirm that newly designed Jet A blends maintain $T_\text{freeze} \le -40$ °C and $\nu_{-20^\circ\text{C}} \le 8.0$ cSt with calibrated conformal safety margins.

---

### Task 1: Low-Temperature Vapor Pressure Anchoring for C10–C20 Alkanes & Iso-alkanes (Flash Point)

**Files:**
- Create: `src/fuelblend/thermo/vapor_pressure_anchoring.py`
- Create: `scripts/v2/anchor_vapor_pressures.py`
- Test: `tests/test_v2_vapor_pressure_anchoring.py`
- Deliverable: `results/v2/improvement/flash_anchored_loso.csv`

**Interfaces:**
- Consumes: `data/processed/thermoml_hc.csv`, `data/unified/compounds.csv`, `data/raw/ruzicka_majer/t01_transcription.md`
- Produces: `data/processed/v2/compounds_anchored_vp.csv`, `fit_anchored_antoine(inchikey, T_min, T_max)` returning $(A, B, C)$ valid over 273–373 K.

- [ ] **Step 1: Write failing unit test for low-T vapor pressure anchoring and isododecane ambient pressure**
```python
def test_anchored_vapor_pressure_isododecane():
    from fuelblend.thermo.vapor_pressure_anchoring import get_anchored_vapor_pressure_kpa
    # Isododecane (2,2,4,6,6-pentamethylheptane) at 25 C (298.15 K)
    # Old WebBook extrapolated: ~0.005 kPa (wrong, boiling point 177 C gives flash 45 C -> Psat must be ~0.2-0.3 kPa)
    p_sat = get_anchored_vapor_pressure_kpa("GWHJZXXIDMPWGX-UHFFFAOYSA-N", 298.15)
    assert 0.10 <= p_sat <= 0.40
```

- [ ] **Step 2: Run test to verify it fails**
Run: `pytest tests/test_v2_vapor_pressure_anchoring.py::test_anchored_vapor_pressure_isododecane -v`
Expected: FAIL with "ModuleNotFoundError"

- [ ] **Step 3: Implement `vapor_pressure_anchoring.py` and `scripts/v2/anchor_vapor_pressures.py`**
Extract low-temperature vapor pressure points (273–373 K) from ThermoML liquid VP records and Cox parameters; refit Antoine $(A, B, C)$ with explicit low-temperature weighting so $P_\text{sat}(298.15\text{ K})$ matches experimental values within 5%.

- [ ] **Step 4: Re-evaluate Flash Point baseline on all 126 Leenhouts hydrocarbon mixture points under LOSO CV**
Compute MAE and over-prediction rate > 5 °C on C12+ mixtures (previously 8.98 °C MAE and +42.4 °C max over-prediction). Verify MAE drops below 3.0 °C and max over-prediction drops below 15 °C.

- [ ] **Step 5: Run tests, commit, push**
Run: `pytest tests/test_v2_vapor_pressure_anchoring.py -v`
Expected: PASS

---

### Task 2: Coutinho-Style Multi-Component Non-Ideal Wax Crystallization (Freezing Point / Cold Flow)

**Files:**
- Create: `src/fuelblend/thermo/coutinho_wax.py`
- Create: `scripts/v2/eval_coutinho_wax.py`
- Test: `tests/test_v2_coutinho_wax.py`
- Deliverable: `results/v2/improvement/freezing_point_coutinho_eval.csv`

**Interfaces:**
- Consumes: `data/processed/v2/wax_dauphin_coutinho.csv`, `data/real_fuels/inventory_reconstructed.csv`, `data/real_fuels/compositions_reconstructed.csv`
- Produces: `predict_freezing_point_coutinho(mole_fractions, pure_tm_K, pure_dhm_J_mol)` returning $(T_\text{freeze\_C}, \text{solid\_phase\_speciation})$.

- [ ] **Step 1: Write failing test for Coutinho multi-component wax crystallization**
```python
def test_coutinho_eutectic_depression():
    from fuelblend.thermo.coutinho_wax import predict_freezing_point_coutinho
    # Binary n-C10 + n-C16 mixture (x_C16 = 0.10)
    # Pure n-C16 Tm = 18.2 C. Ideal single-solute WAT = ~3 C. Eutectic depressed WAT = ~-2 C.
    t_freeze = predict_freezing_point_coutinho([0.90, 0.10], ["CCCCCCCCCC", "CCCCCCCCCCCCCCCC"])
    assert -10.0 <= t_freeze <= 0.0
```

- [ ] **Step 2: Run test to verify it fails**
Run: `pytest tests/test_v2_coutinho_wax.py::test_coutinho_eutectic_depression -v`
Expected: FAIL

- [ ] **Step 3: Implement `src/fuelblend/thermo/coutinho_wax.py`**
Implement the predictive local composition solid-solution model (UNIQUAC solid-state activity coefficient) calibrated on Dauphin binary and multi-component n-alkane crystallization curves, accounting for carbon-number polydispersity and solid-solid orthorhombic transitions.

- [ ] **Step 4: Re-evaluate on 57 CRC aviation real fuels under strict LOO CV**
Compare against the v2 flat offset model ($T_\text{freeze} = \text{WAT} - 14.58$ °C, MAE 4.05 °C, $r=0.443$). Verify correlation $r$ improves from 0.44 to $>0.65$ and MAE drops below 3.2 °C.

- [ ] **Step 5: Run tests, commit, push**
Run: `pytest tests/test_v2_coutinho_wax.py -v`
Expected: PASS

---

### Task 3: ASTM D341 / Walther Temperature-Viscosity Curvature Model (Viscosity)

**Files:**
- Create: `src/fuelblend/thermo/walther_viscosity.py`
- Create: `scripts/v2/eval_walther_viscosity.py`
- Test: `tests/test_v2_walther_viscosity.py`
- Deliverable: `results/v2/improvement/viscosity_walther_eval.csv`

**Interfaces:**
- Consumes: `data/processed/v2/distillation_real_fuels.csv`, `results/v2/A9/hindcast_predictions.csv`
- Produces: `walther_kinematic_viscosity(T_C, A, B)` returning $\nu$ in cSt satisfying $\log_{10} \log_{10}(\nu + 0.7) = A - B \log_{10}(T + 273.15)$.

- [ ] **Step 1: Write failing test for ASTM D341 Walther viscosity curvature**
```python
def test_walther_cryogenic_curvature():
    from fuelblend.thermo.walther_viscosity import fit_and_predict_walther_loo
    # Verify Walther curvature beats linear 1000/T at -40 C on real fuels
    res = fit_and_predict_walther_loo()
    assert res["rmse_minus40C"] < 0.150  # Lower than 2-param Refutas (0.229)
```

- [ ] **Step 2: Run test to verify it fails**
Run: `pytest tests/test_v2_walther_viscosity.py::test_walther_cryogenic_curvature -v`
Expected: FAIL

- [ ] **Step 3: Implement `src/fuelblend/thermo/walther_viscosity.py`**
Implement the ASTM D341 double-logarithmic viscosity equation combined with Refutas VBN mixture blending:
$$\log_{10} \log_{10}(\nu_\text{mix}(T) + 0.7) = A_\text{mix} - B_\text{mix} \log_{10}(T)$$
Evaluate under Leave-One-Fuel-Out CV across all 57 CRC fuels at $-40$ °C, $-20$ °C, and $+20$ °C.

- [ ] **Step 4: Verify error reduction**
Verify overall log-RMSE drops from 0.1411 (linear slope) to $< 0.1100$, and $-40$ °C log-RMSE drops from 0.2291 to $< 0.1250$.

- [ ] **Step 5: Run tests, commit, push**
Run: `pytest tests/test_v2_walther_viscosity.py -v`
Expected: PASS

---

### Task 4: Re-running Inverse Blend Designs (T9 Hardening for Jet A & Diesel)

**Files:**
- Modify: `src/fuelblend/design/inverse_design_v2.py`
- Create: `scripts/v2/run_inverse_design_upgraded.py`
- Test: `tests/test_v2_inverse_design_upgraded.py`
- Deliverable: `results/v2/improvement/upgraded_designs.csv`

**Interfaces:**
- Consumes: Upgraded flash model (Task 1), upgraded Coutinho wax model (Task 2), upgraded Walther model (Task 3)
- Produces: Optimal Jet A and Diesel No. 2 blend formulations evaluated under the upgraded physics engine.

- [ ] **Step 1: Write test for upgraded inverse design evaluation**
```python
def test_upgraded_inverse_design_feasibility():
    from fuelblend.design.inverse_design_v2 import evaluate_fuel_design_upgraded
    # Test Jet A 4-component blend with Coutinho wax model and Walther viscosity
    res = evaluate_fuel_design_upgraded("jet", {"isododecane": 0.58, "n-dodecane": 0.12, "decalin": 0.15, "tetralin": 0.15})
    assert res["all_constraints_passed"] is True
    assert res["freezing_point_coutinho_C"] <= -40.0
    assert res["viscosity_walther_m20C"] <= 8.0
```

- [ ] **Step 2: Run test to verify it fails**
Run: `pytest tests/test_v2_inverse_design_upgraded.py::test_upgraded_inverse_design_feasibility -v`
Expected: FAIL

- [ ] **Step 3: Update `inverse_design_v2.py` and run optimizer**
Incorporate the upgraded thermodynamic evaluators into the simplex objective and constraint functions. Re-run SLSQP multi-start optimization.

- [ ] **Step 4: Save `results/v2/improvement/upgraded_designs.csv`**
Document all constraints, safety margins, and compliance flags.

- [ ] **Step 5: Run tests, commit, push**
Run: `pytest tests/test_v2_inverse_design_upgraded.py -v`
Expected: PASS

---

### Task 5: Master Version-by-Version Ablation & Validation Audit (The Paper Ablation Story)

**Files:**
- Create: `scripts/v2/generate_paper_ablation_table.py`
- Deliverable: `results/v2/improvement/paper_ablation_table.csv`
- Document: `docs/progress/v2/THERMODYNAMIC_IMPROVEMENTS_SUMMARY.md`
- Test: `tests/test_v2_ablation_table.py`

**Interfaces:**
- Consumes: All baseline and upgrade deliverables (A2, A4, A5, C1, C2, C4, D1, Task 1, Task 2, Task 3)
- Produces: Unified version-by-version error comparison table tracing accuracy progression from v1 $\rightarrow$ v2 $\rightarrow$ v2-Thermodynamic.

- [ ] **Step 1: Write test asserting monotonic error reduction across ablation stages**
```python
def test_ablation_monotonic_error_reduction():
    df = pd.read_csv("results/v2/improvement/paper_ablation_table.csv").set_index("version")
    # Verify flash error: v1 baseline > v2 baseline > v2 anchored
    assert df.loc["v1_raw_physics", "flash_mae_C"] > df.loc["v2_raoult_lfl", "flash_mae_C"] > df.loc["v2_anchored_thermo", "flash_mae_C"]
    # Verify freezing error: v1 ideal WAT > v2 offset > v2 Coutinho
    assert df.loc["v1_ideal_wat", "freeze_mae_C"] > df.loc["v2_offset", "freeze_mae_C"] > df.loc["v2_coutinho", "freeze_mae_C"]
```

- [ ] **Step 2: Run test to verify it fails**
Run: `pytest tests/test_v2_ablation_table.py::test_ablation_monotonic_error_reduction -v`
Expected: FAIL

- [ ] **Step 3: Implement `scripts/v2/generate_paper_ablation_table.py`**
Aggregate cross-validated test metrics across all versions:
1. `v1_raw_physics` (uncalibrated Refutas, raw Le Chatelier, Dulong NHOC, raw WAT)
2. `v2_recalibrated` (2-param Refutas, group-additivity NHOC, flat WAT offset)
3. `v2_thermodynamic_anchored` (ASTM D341 Walther, anchored low-T VP, Coutinho multi-solid wax SLE)

- [ ] **Step 4: Save `results/v2/improvement/paper_ablation_table.csv` and write report**
Write `docs/progress/v2/THERMODYNAMIC_IMPROVEMENTS_SUMMARY.md` presenting the complete progression for the publication manuscript.

- [ ] **Step 5: Run tests, commit, push**
Run: `pytest tests/test_v2_ablation_table.py -v`
Expected: PASS
