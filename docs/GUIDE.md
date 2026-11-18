# UGRP Hydrocarbons — the complete guide

*Predicting jet-fuel blend properties from composition: physics first, machine learning only where it was proven to help, a measured safety margin on every number, and an honest label wherever the model has no skill.*

Authors: Namyaa Parmar, Kahaan Shah, Arijit Ganguli, A. S. K. Sinha — Ahmedabad University.
This guide explains **every term, every decision, the current status, every important code file and how to run it**. Numbers come from `docs/final_numbers_v5.md` and `Fuel_Blend_Methodology.pdf`; if a number here disagrees with those, they win.

---

## 0. The 90-second version (say this first)

1. **Problem.** A fuel is a blend of hundreds of molecules. Specs (Jet A: flash point ≥ 38 °C, freezing point ≤ −40 °C, density 775–840 kg/m³, viscosity ≤ 8 cSt at −20 °C) are checked by lab tests. We want to predict these from the *composition* so we can design blends on a computer.
2. **Approach.** Every property is computed with a classical physics/chemistry rule fed with measured properties of the pure molecules. A small learned correction is added only where it beat the physics on mixtures the model never saw.
3. **Honesty.** Each model is tested on whole mixture *systems* held out from training, against a "just guess the average" baseline. A model that does not beat the guess is labelled *no skill shown* and is not used to certify.
4. **Safety margin.** From the held-out errors we compute a one-sided 90 % conformal margin. A blend is "certified" only if it still meets the spec after the margin is applied.
5. **Result.** Viscosity error 8 % → 3 % with the learned pair correction; flash point within ≈ 2.2 °C on 70 real-fuel groups; density within ≈ 0.15 %. Freezing, smoke point and heat of combustion have **no proven skill on real fuels** and we say so. Inverse design finds 3 of 16 blend/spec combinations that pass **under the model only** — nothing was measured in a lab.
6. **Why it is credible.** Many negative results are kept (neural mixture models, GNN activity coefficients, ML sooting, ML cetane): they failed pre-declared tests and we left them in the repo on purpose.

---

## 1. The problem and the goal

**Goal.** Predict viscosity, density, flash point, freezing point, heat of combustion and smoke point of a hydrocarbon fuel blend (jet / diesel / RP-1 style) from its composition, with calibrated uncertainty, and then run **inverse design**: find a blend of available molecules that meets the specification.

**Why a plain neural network is wrong here.** Earlier work (the ChemCon abstract) trained a network on *synthetic* blends whose labels came from the blending rule itself (Refutas). It "achieved 0.4 % error" — but only because it learned to copy the rule. That number says nothing about real fuels. The whole of this repository exists to avoid that circularity: **real measured data, physics baselines, whole-system hold-out, and a trivial-baseline comparison.**

**What "done" means for each property** (priority order): viscosity, flash point, density (core, certifiable); freezing point (certifiable only for known molecules with verified constants, margin from lab data); smoke point (regression, in range only); heat of combustion (reported, no margin); octane (gasoline research result, not in the API).

---

## 2. Glossary — every term you may be asked about

| Term | Meaning (plain words) |
|---|---|
| **Composition bases** | Mass fraction *w*, mole fraction *x*, volume fraction *v*. Different rules need different ones. x_i = (w_i/M_i)/Σ(w_j/M_j); v_i = (w_i/ρ_i)/Σ(w_j/ρ_j). Temperatures inside equations are always kelvin: T_K = T_°C + 273.15. |
| **SMILES / InChIKey** | Text codes that identify a molecule exactly (SMILES = structure string, InChIKey = fixed-length hash). The API accepts either. |
| **cSt (centistokes)** | Unit of *kinematic* viscosity = mm²/s. Water at 20 °C ≈ 1 cSt. ν = 1000·η/ρ (η in Pa·s, ρ in kg/m³). |
| **Dynamic vs kinematic viscosity** | η (mPa·s) is resistance to shear; ν = η/ρ is how fast the fluid flows under gravity. Specs use kinematic. |
| **Refutas** | Viscosity blending rule: map each ν to a Viscosity Blending Number (VBN) that averages linearly by mass. |
| **VBN** | VBN = 14.534·ln(ln(ν+0.8)) + 10.975. |
| **Grunberg–Nissan** | Classic rule: ln η_mix = Σ x_i ln η_i + Σ_{i<j} x_i x_j G_ij. Our ML learns the pair term. |
| **Ideal volume mixing** | Volumes add: 1/ρ_mix = Σ w_i/ρ_i. Hydrocarbon mixtures deviate by well under 1 %. |
| **Excess volume** | The (small) volume change on mixing. ML option learns it per pair of molecules. |
| **Antoine equation** | Vapour pressure of a pure liquid: log10 P_sat = A − B/(T + C). Here P in kPa, T in K. |
| **Raoult's law** | Ideal mixture: partial pressure = x_i · P_sat,i. |
| **LFL** | Lower flammability limit: smallest fuel-vapour fraction in air that can burn (mole fraction). |
| **Le Chatelier's rule** | A vapour mixture is flammable when Σ y_i/LFL_i = 1. |
| **Flash point** | Lowest temperature at which the vapour above the liquid ignites. Spec: ≥ 38 °C for Jet A. Over-predicting is the *dangerous* direction. |
| **Zabetakis correction** | LFL falls with temperature: LFL(T) = LFL₂₅ (1 − 0.000721 (T − 25)). |
| **Anchored Antoine** | Antoine constants re-anchored to measured vapour pressure in the relevant window (0–100 °C) because literature ranges for C10+ alkanes start at 95–190 °C. |
| **Registry overlay** | Extra Antoine + LFL entries (bicyclohexyl, cis-decalin) so real-fuel compositions mapped to those molecules can be flash-predicted. On by default; `FUELBLEND_FLASH_REGISTRY_ADDITIONS=0` turns it off. |
| **Freezing point (jet)** | Temperature at which the first wax crystals appear. Spec: ≤ −40 °C (Jet A), ≤ −47 °C (Jet A-1). Under-predicting is dangerous. |
| **Schroeder–van Laar (SvL)** | Solid–liquid equilibrium: 1/T_i = 1/T_m,i − R ln(x_i γ_i)/ΔH_fus,i; the blend freezes at max_i T_i. |
| **UNIFAC** | Group-contribution method that estimates activity coefficient γ_i from chemical groups. Main freezing route; ideal (γ = 1) is the variant. |
| **Activity coefficient γ** | How far a molecule departs from ideal behaviour in a liquid mixture. |
| **ΔH_fus, T_m** | Heat of fusion and melting point of the pure solid. |
| **Smoke point** | Tallest flame (mm) that burns without smoking; aromatics lower it. Regression SP = 28.84 − 0.324 A − 0.791 D (A aromatics vol %, D diaromatics vol %). |
| **YSI** | Yield sooting index, a measured soot-tendency scale per molecule. Our YSI route failed its gate. |
| **Heat of combustion (net)** | Energy per kg released burning the fuel (MJ/kg). Group additivity: Σ w_i ΔH_c,i. Spec ≥ 42.8. |
| **D86 / Engler** | Distillation curve (T10, T50, T90…). Engler = our default TBP-offset model; Rayleigh = opt-in alternative that failed the adoption rule. |
| **DCN / cetane** | Diesel ignition quality. ML version "not feasible" (only 31 usable systems). |
| **RON / MON / Co-Optima** | Gasoline octane measures; Co-Optima = DOE dataset of 65 blends. Research only. |
| **Simplex** | The set of non-negative fractions summing to 1; the search space of inverse design. |
| **System** | One set of molecules (e.g. decane + toluene) measured at any ratio/temperature. All its points stay together in train **or** test. |
| **GroupKFold** | Cross-validation that never splits a system across train/test. |
| **Locked (held-out) set** | Systems set aside *before* any model was built and scored exactly once. |
| **Gate (R-xxx-n)** | A pre-declared pass/fail rule, committed to git before results exist (e.g. R-MLV-1). |
| **Trivial comparator** | "Predict the average of the other groups" (leave-one-group-out constant). A model has skill only if it beats this. |
| **Verdict** | `beats_trivial`, `no_skill_shown`, `worse_than_trivial` (group-bootstrap comparison vs trivial). |
| **Bias / MAE** | Bias = mean signed error (pred − meas). MAE = mean absolute error. |
| **ln-error for viscosity** | |ln η_pred − ln η_meas|; for small values ≈ relative error (0.031 ≈ 3 %). |
| **Conformal margin** | Distribution-free allowance from ranked held-out errors (see §6). |
| **Certified** | Margin exists ∧ blend in tested domain ∧ spec met after margin. *Model-level only.* |
| **Domain / `in_domain`** | Known molecules with their own data, 233–393 K, class pairs seen ≥ 5 times in training; no stand-ins. |
| **T1 / T2 / T3 / CRC647_r3** | Evidence tiers: T1 lab mixtures; T2 real fuels with measured composition; T3 AJFTD bulk-property regressions; CRC647_r3 real fuels with *reconstructed* composition. |
| **AJFTD** | Alternative Jet Fuels Test Database (2,601 fuels scraped; 116 with GC×GC composition, 113 mappable). |
| **GC×GC** | Two-dimensional gas chromatography: composition by molecule class × carbon number. |
| **CRC 647** | Coordinating Research Council report: 57 aviation fuels with property tables. |
| **ThermoML** | NIST archive of measured thermophysical data (density, viscosity, vapour pressure…). |
| **Leenhouts benchmark** | Curated pure + mixture flash/viscosity/DCN/MON benchmark (Zenodo 19914760). |
| **LightGBM** | Gradient-boosted decision trees; the learner for pair corrections. |
| **DeepSets / GNN / set transformer** | Neural architectures for mixtures; all failed our gates. |
| **BASE / delta_*** | `results/v5/BASE/` is the frozen scoreboard; each `delta_*` folder is one experiment measured against it. |
| **CLM-nn** | Claim number in `docs/claims_register.md`. |

---

## 3. The physics, property by property

All of this is implemented once, in the **designated scoring route** (`src/fuelblend/eval/composition_route.py`) and exposed through `src/fuelblend/api/predict_route.py::predict_fuel`.

### 3.1 Density — ideal volume mixing
`1/ρ_mix = Σ w_i/ρ_i` (all ρ at the same T). Code: `mixture/common.py::density_idvol`.
Evidence: lab 1.33 kg/m³ MAE (≈ 0.17 %); real fuels beat the trivial guess. **Reconstructed-composition CRC fuels do not** (14.3 vs 13.2 kg/m³) — that is a composition-reconstruction limit, not a density-law failure.
**ML option** (`ml/density_pair.py`): `V = Σ x_i M_i/ρ_i + Σ_{i<j} x_i x_j g(d_i,d_j,T)`, `ρ = Σ x_i M_i / V`. 1.20 → 0.90 kg/m³ on 41 held-out binary systems, margin ±4.42. **Fails on ternaries** (3.64 → 3.92), so it is a labelled variant for known binaries only; ideal volume stays the main value.

### 3.2 Viscosity — Refutas, plus learned pair correction
Convert each pure ν_i to VBN_i = 14.534 ln(ln(ν_i+0.8)) + 10.975; average by *mass*: VBN_mix = Σ w_i VBN_i; invert: ν_mix = exp(exp((VBN_mix − 10.975)/14.534)) − 0.8. Valid ν ≳ 0.2 cSt. Code: `mixture/common.py::refutas_ln_eta`, `eval/composition_route.py::refutas_nu_safe`.
**Why double-log:** viscosities of blend components differ by orders of magnitude; a mass-weighted average of ν would be dominated by the thick part. VBN puts them on a scale where mass-averaging works.
**Pure-component data:** a *frozen* table (`data/cache/pure_grid.npz`, SHA-256 pinned and checked at every load) of ln η versus T for 1,042 compounds (114 with a full grid), starting at 230 K. A close relative may stand in if a molecule has no data — flagged, never certified.
**ML pair correction** (`ml/viscosity_pair.py`, models in `models/v5/ml_viscosity_pair`):
`ln η_mix = ln η_Refutas + Σ_{i<j} x_i x_j g(d_i, d_j, T)`.
LightGBM (400 trees, 15 leaves, lr 0.05, seed 42, deterministic) on *pair features* — mean and |difference| of carbon number, molar mass, ring counts, CH₃ count, quaternary C, ln η, molar volume, the class pair and 1/T. Symmetric in the pair, zero for a pure compound, works for any number of components because only pairs are learned. Trained on 264 binary systems, graded once on 70 locked systems.
Gate results (Refutas → ML, MAE in ln η): R-MLV-1 dev 0.082 → 0.060; R-MLV-2 locked 0.054 → 0.031; R-MLV-3 ternary+ 0.098 → 0.042 (per-system); R-MLV-4 cold (<283 K) 0.099 → 0.061. All passed.
**Designation rule (user decision):** ML is the main viscosity *inside its domain* (known molecules, ≥ 283 K gate for scoring, class pairs with ≥ 5 training systems). Elsewhere the main value is Refutas with a warning.
**Honest limit:** on real fuels below 263 K Refutas is *worse than a constant* (MAE 0.53 ln vs 0.25), so no real-fuel cold margin exists (margin = ∞).

### 3.3 Flash point — Le Chatelier + Raoult + Antoine (nothing fitted to blends)
1. Pure vapour pressure: `log10 P_i^sat = A_i − B_i/(T + C_i)` (kPa, K), constants anchored to measured vapour pressure in the 0–100 °C window.
2. Raoult (ideal): `y_i = x_i P_i^sat(T)/P_total`, P_total = 101.325 kPa.
3. Flammable when `H(T) = Σ x_i P_i^sat(T)/(P_total LFL_i) = 1`; solve by `brentq` on 150–650 K. LFL corrected to temperature (Zabetakis).
Code: `baselines/flash.py` (`psat_antoine`, `lfl_temperature_correction`, `hazard_function`, `predict_flash_le_chatelier`).
**Why this is the fair test:** nothing is fitted to blend flash points, so testing on blends is honest.
Evidence: 70 real-fuel groups MAE 2.20 °C (trivial 5.16), bias +0.06 °C; AJFTD GC×GC fuels (31) bias −0.15, MAE 3.00, beats trivial; lab mixtures margin 2.207 °C (35 systems).
Applied flash margin is **2.207 °C** (the yaml value 3.699 and the scoreboard `margin90` 3.699 are the superseded pre-Q2 figures).
Limits: D56 real-fuel check rests on **one** GC fuel (POSF-10325) → no margin; ideal-mixture assumption fails for some systems (mesitylene outlier excluded, disclosed).

### 3.4 Freezing point — Schroeder–van Laar with UNIFAC
`1/T_i = 1/T_m,i − R ln(x_i γ_i)/ΔH_fus,i`; `T_freeze = max_i T_i`. Main route: γ from original UNIFAC (`baselines/flash_unifac.py` and the SLE code in `eval/composition_route.py::freezing_C`, `DESIGNATED_FREEZING_MODEL = "svl_unifac"`). Variant: ideal (γ=1). Melting points / ΔH_fus come from a *verified* constants table (`data/raw/sle_constants/tm_dhfus_verified.csv`); bicyclohexyl uses ΔH_fus plus its solid–solid transition (13.94 kJ/mol).
Lab liquidus set (73 points, 12 systems): UNIFAC MAE 3.93 °C, bias +2.47, **margin 13.84 °C**; ideal 4.47 / +1.17 / 16.98; constant guess 15.63.
**Real fuels: no skill.** Pooled 50 fuels MAE 8.42 vs constant 6.41; locked 8 fuels 9.31 vs 5.18 (`worse_than_trivial`). Bias is +4 to +6 °C (safe side), so the 13.84 margin "covers" 98 % trivially. We use freezing *only* through the lab margin, for known molecules with verified constants.
Known data issue: the o-/m-xylene eutectic composition appears transposed in the ingest script (not fixed; margins unaffected).

### 3.5 Smoke point and heat of combustion
Smoke: `SP = 28.84 − 0.324 A − 0.791 D` (fitted on 100 jet fuels; trusted for 8 ≤ A ≤ 23 vol % and 19–25 mm). 0.90 mm MAE on the 100-fuel set (beats trivial), no skill shown on the broader set (1.68 vs 1.65). YSI route failed gate (median rel. error 0.229 vs 0.20).
Heat: group additivity `Σ w_i ΔH_c,i`; reported, **no margin**; worse than trivial on reconstructed CRC fuels.

### 3.6 Octane (research only)
Blending octane value: `RON = (1−v) RON_base + v BOV`. Gaussian process / GBM on blendstock descriptors, 65 Co-Optima blends, 18 dev blendstocks + 4 held out: RON 2.29, MON 1.15 (CV). Too few held out to certify; not in the API.

---

## 4. The machine-learning story

**Principle:** ML learns only what the ideal rule misses — pairwise molecular interactions — and is used only where a pre-declared gate passed on locked systems.

| Model | Base | Held-out result | Status |
|---|---|---|---|
| Viscosity pair correction (LightGBM) | Refutas | 0.054 → 0.031 ln (70 systems) | **main value in domain** |
| Density excess volume (LightGBM) | ideal volume | 1.20 → 0.90 kg/m³ (41 systems); ternaries fail | option (known binaries) |
| Octane GP | BOV rule | RON 2.29 / MON 1.15 (CV) | research only |
| Neural mixture encoders (DeepSets, GNN, set transformer) | Refutas / ideal volume | none pass; GNN worse (0.031 → 0.040) | **failed, kept** |
| GH-GNN activity coefficients in flash | ideal Le Chatelier | +0.70 °C worse lab, +1.09 °C worse real fuels | **failed, kept** |
| GC×GC real-fuel models | physics + offset | worse in 10 of 12 leave-one-source-out cases | **does not transfer** |
| ML YSI for smoke | class mean | rel. error 0.229 > 0.20 | **failed** |
| ML DCN | – | only 31 usable systems (~40 needed) | **not feasible** |
| Antoine interpolation, Tm/ΔH_fus estimators, multicomponent freezing | – | fail declared gates | failed, flagged unused |

**Process every model followed:** (1) pass rules committed to git first; (2) development by cross-validation over whole systems; (3) locked systems scored once; (4) bootstrap CI over systems; (5) compare to trivial constant; (6) failed models stay in code and in `docs/v5_negative_results.md`.

---

## 5. Data — what we used and what we refused to use

| Source | What | Notes |
|---|---|---|
| Leenhouts benchmark (Zenodo 19914760) | curated pure + mixture flash, viscosity (ln mPa·s), DCN, MON | primary lab mixture evidence |
| NIST/TDE viscosity (CheMixHub, Bilodeau) | ~239k binary viscosity rows | value column is ln(η/Pa·s) despite label |
| ThermoML archive | measured density, viscosity, vapour pressure | re-extracted after parser fixes (phase, constraint, units) |
| NIST WebBook | Antoine constants (bar) → kPa (+2) | ranges don't reach 40–100 °C for C10+ alkanes → anchored/refit |
| Zabetakis Bulletin 627 | measured LFL values | hand-read, ratio-gated |
| Li et al. (Sun/Krakauer/Morgan 2020) | pure flash points | loader `load_flashpoint_li` |
| AJFTD | 2,601 real fuels; 116 GC×GC | login-gated, not committed |
| CRC 647 | 57 jet fuels | compositions *reconstructed* from measured hydrocarbon types (Table G.1) |
| NISTIR 6646 / Fortin 2012 | RP-1 characterisation, distillation | RP-1 11-blend recipes confidential |
| Co-Optima | 65 octane blends | research |

**Refused:** `DatasetsFuelProperties.xlsx` (50,000 *synthetic* rows) and the v3 fuel corpus + optimiser (formula-derived; 92 of 153 fuels dropped, quarantined in git history). **Rule:** when a result looks too good, look for circular labelling.
De-duplication: drop exact repeats, group by (property, system, composition, T, P), take the median, flag (never silently delete) outliers and conflicts; the group key for CV is `system_id`.

---

## 6. Uncertainty — the conformal margin and "certified"

For each held-out system *g* take the **worst dangerous-side error** `s_g = max_{points in g}(error in the dangerous direction)`. With *n* systems and α = 0.10:
`k = ⌈(n+1)(1−α)⌉`, `margin = s_(k)` (the k-th smallest score). If *k > n* the margin is **∞** — not enough systems (needs ≥ 9 at 90 %). Code: `src/fuelblend/uq/conformal.py` (`conformal_rank`, `group_scores`, `one_sided_margin`, `two_sided_margin`, `coverage`). Applied margins live in `configs/uncertainty_v5.yaml`.

| Property | Dangerous side | Margin | Basis |
|---|---|---|---|
| Flash | too high | **2.207 °C** | 35 lab systems |
| Freezing (UNIFAC) | too low | **13.84 °C** | 12 lab systems |
| Freezing (ideal, variant) | too low | 16.98 °C | 12 |
| Density (ideal) | both | +4.59 / −3.25 kg/m³ | 167 |
| Density (ML, binaries) | both | ±4.42 kg/m³ | 41 locked |
| Viscosity ML ≥ 283 K | both | ±0.134 ln (≈ ±14 %) | 264 |
| Viscosity ML 263–283 K | both | ±0.608 ln | 16 |
| Viscosity ML < 263 K | both | ±0.287 ln | 20 |
| Viscosity Refutas ≥ 283 K | both | ±0.186 ln | 255 |
| Smoke (jet class mix) | too high | 4.11 mm | 7 source groups |
| Heat | – | none | – |

**Worked example.** Jet A needs flash ≥ 38 °C; margin 2.207 → certified only if predicted flash ≥ 40.207 °C. "Predicted 38.1 °C" is *not* certified.
**Certified =** (a) a finite margin exists, (b) the blend is `in_domain`, (c) the limit still holds after the margin. Otherwise `meets_at_point_estimate` or `fails`.
**Important:** the freezing and cold-viscosity margins come from mostly **binary** lab systems; applying them to 6-component designs is a transfer, not a validation.

---

## 7. Honest evaluation — the rules that make the numbers trustworthy

1. **Real data only**; synthetic labels banned.
2. **Hold out whole systems** (GroupKFold on `system_id`; leave-one-source-out for real fuels).
3. **Pre-declare** every test (README with the pass rule committed *before* the result).
4. **Beat a trivial guess** (leave-one-group-out constant) or the verdict is `no_skill_shown`.
5. **Bootstrap over groups** for CIs.
6. **Keep negative results**; never delete a failed model.
7. **Descriptive vs gated:** T2 real-fuel cells are *descriptive* (every T2 fuel has been seen), not held-out claims.
8. **Frozen BASE:** `results/v5/BASE/` is hash-protected by `SHA256SUMS_v5_base.txt`; any change is a visible drift.

---

## 8. Inverse design

`design/engine_v4.py` searches the simplex (non-negative fractions summing to 1) with multistart SLSQP-style optimisation (10 starts, seed 42), maximising the minimum slack across constraints, then **re-predicts the result with the same API** and judges limit by limit — so design and prediction cannot disagree. Constraints come from config only (specs and margins). `strict_constraints=True` (default) makes the engine stop with `ModelUnavailable` instead of silently dropping a constraint.
Pools: P0 = 7 pure molecules (n-decane, n-undecane, n-dodecane, butylcyclohexane, trans-decalin, 1,2,4-trimethylbenzene, isooctane); P0x-a (+tetralin, n-tridecane); P0x-b (+hexylbenzene); P0y-b (+bicyclohexyl); P0y (+bicyclohexyl, cis-decalin); P1–P3 mix real blendstocks (GC×GC Jet A base, HEFA, ATJ, SAK) with pure molecules.
**Result (8 pools × Jet A / Jet A-1 = 16 designs):** certified under the model — P0y Jet A, P0y-b Jet A, P0y-b Jet A-1. All three carry labels `d86_model_dependent`, `freezing_route_dependent`, `overlay_dependent`: they pass on Engler D86 but fail the Rayleigh check; they pass freezing on UNIFAC+13.84 but would exceed the limit on ideal SvL+16.98; they need the bicyclohexyl/cis-decalin registry additions. The certified Jet A-1 has 37.8 wt % dicycloparaffin, above the AJFTD GC×GC maximum (30 %). Minimum slacks (0.003–0.013) show the optimiser sitting on the margin boundary, **not robustness**. No design has been measured; a lab test was ruled out (decision 2026-11-14).

---

## 9. Where we are (status, 17 Nov 2026)

* **Working & evidenced:** physics routes for all properties; viscosity and density ML corrections (lab-gated); conformal margins; 16-design table; frozen, hash-checked scoreboard; API accepting SMILES/InChIKeys.
* **No proven skill on real fuels:** freezing, heat of combustion (reconstructed), pour point, cold viscosity (<263 K, Refutas), CRC reconstructed density.
* **Not done / not claimed:** any laboratory validation; diesel validation; multicomponent freezing margin; real-fuel ML density testing; Jet A-1 limits verified against Def Stan 91-091; lab freezing of hexadecane mixtures (Oakley 2017) and smoke of 23 coal/JP-8 fuels (Yang 2006) are in the data but unscored.
* **Tests:** last full run 785 passed / 0 failed (before the simplification); after simplification the BASE suite and the touched suites were re-run (see the repository's final commit message).

---

## 10. Decision log — what we decided and why

| Decision | Why |
|---|---|
| Physics-first, ML as residual | blending rules are right to ~1 % for density and a few % for viscosity; a net learning the rule proves nothing |
| Hold out systems, not rows | rows of one system are near-copies; row-level CV hides over-fitting |
| One-sided conformal margin | specs are one-sided; direction of the error matters for safety |
| Margin ∞ when n too small | refuse to invent a number the data cannot support |
| ML viscosity designated in domain | passed all 4 gates; Refutas stays as `refutas_raw` outside the domain |
| ML density is only a variant | failed the ternary gate |
| UNIFAC-SvL designated for freezing | smaller margin (13.84 vs 16.98), smaller error, no case it could not compute (declared rule) |
| Flash-registry overlay on by default | enables bicyclohexyl/cis-decalin; `FUELBLEND_FLASH_REGISTRY_ADDITIONS=0` restores old |
| Engler D86 kept; Rayleigh opt-in | Rayleigh T10 bias +11.6 °C failed the adoption rule |
| Smoke = bulk regression, not YSI | YSI had 0 evaluable fuels / failed gate (user: "pick the less-error one") |
| No lab test | ruled out by the user; certification stays model-level and says so |
| Keep failed models | reproducibility of negative results is part of the science |
| Frozen BASE + hash manifest | any silent change to evidence is detected |
| v3 corpus quarantined | its fuel compositions and optimiser were derived from the formulas being tested |

History in one line: v1 data pipeline → v2/v3 (circular, quarantined) → v4 integrity rebuild on real fuels → v5 full baseline, ML components, adoption of user decisions, final numbers pack. The full pre-simplification state is the git tag `v5-full-2026-11-16`.

---

## 11. Code map — which file does what

**Entry points**
| File | Role |
|---|---|
| `src/fuelblend/api/predict_route.py` | `predict_fuel(composition, T_K, fuel_type)` — the public prediction function; loads the frozen grid (`GRID_SHA256` check), resolves names/SMILES/InChIKeys, returns each property with `value`, `method`, `status`, `in_domain`, `lo/hi`, `variants`, `warnings` |
| `src/fuelblend/api/predict_blend.py` | `predict(...)` — v4-compatible wrapper on the same route; reads margins from the yaml (`_load_margins`, `_applied_margin`) |
| `scripts/v5/base/reproduce_base.py` | rebuilds the scoreboard, runs the BASE tests, checks every file in `results/v5/BASE` against `SHA256SUMS_v5_base.txt` |
| `scripts/v5/run_design_userpool_v5.py`, `run_jet_designs_v5.py` | inverse design on a user pool / the jet pools |

**Physics**
| File | Role |
|---|---|
| `baselines/flash.py` | Antoine, LFL temperature correction, hazard function, `brentq` flash solver |
| `baselines/flash_unifac.py` | UNIFAC activity coefficients (variant / freezing) |
| `baselines/cold_flow.py`, `mixture/cold_flow.py` | freezing/cloud helpers |
| `baselines/heat_of_combustion*.py`, `sooting.py`, `distillation*.py`, `cetane.py`, `octane.py` | heat, smoke/YSI, D86 (Engler/Rayleigh), DCN, octane |
| `mixture/common.py` | `density_idvol`, `refutas_ln_eta`, torch Newton flash, `CompoundRegistry`, `PurePropertyProvider` |
| `eval/composition_route.py` | the **designated scoring route**: `refutas_nu_safe`, `freezing_C` (`ideal_svl` / `svl_unifac`), `RouteCtx`, `predict_properties` |

**Machine learning**
| File | Role |
|---|---|
| `ml/viscosity_pair.py` | LightGBM pair correction; `FEATURES`, `LGB_PARAMS`, `ln_correction` (= 0.5·xᵀGx), `domain_status` |
| `ml/density_pair.py` | same design for excess volume |
| `ml/set_models.py` | DeepSets / GNN / set transformer (failed, kept) |
| `ml/ysi_model.py` | ML YSI (failed, kept) |
| `mixture/t8_deepsets.py`, `t8_pairwise.py`, `set_encoder_v4.py`, `structured_residual.py`, `bilinear_pairwise.py` | earlier mixture-ML experiments |
| `models/v5/` | trained model files |

**Uncertainty & evaluation**
| File | Role |
|---|---|
| `uq/conformal.py` | conformal rank/margins/coverage |
| `eval/cellkit.py`, `group_stats.py`, `metrics.py` | scoreboard cells, group bootstrap, verdict logic |
| `splits/locked_sets.py`, `splits/` (repo root) | locked-set definitions |
| `configs/uncertainty_v5.yaml` | every applied margin and its provenance |

**Design**
| File | Role |
|---|---|
| `design/engine_v4.py` | optimiser + judge (statuses, `_Sides`, strict constraints) |
| `design/models_v4.py`, `library_route.py` | models per constraint; route-backed molecule library |
| `design/blendstock_v5.py`, `smoke_margin_v5.py` | real blendstocks, smoke margin |

**Data**
| Folder / file | Role |
|---|---|
| `src/fuelblend/data/` | one loader per source (`loaders_*.py`), ThermoML extractor, provenance ledger, reconstruction of CRC compositions (`reconstruction_r3.py`) |
| `data/unified/` | `measurements.csv`, `measurements_dedup.csv`, `compounds.csv` (schema in `AGENTS.md` §6) |
| `data/real_fuels/v4/`, `v5/` | AJFTD / CRC real-fuel tables |
| `data/cache/pure_grid.npz` | frozen pure-property grid |

**Evidence & docs**
| Path | Role |
|---|---|
| `results/v5/BASE/` | frozen scoreboard + one `delta_*` folder per experiment (each has a README declaring its readings before results) |
| `results/v5/final/` | tables/figures behind `docs/final_numbers_v5.md` |
| `docs/claims_register.md` | every claim CLM-xx and its status |
| `docs/model_decision.md` | model choices and ML component decisions |
| `docs/v5_negative_results.md` | the failed attempts |

---

## 12. How to run things

```powershell
# one-time
.venv\Scripts\python.exe -m pip install -r requirements.txt

# 1. Rebuild the frozen scoreboard and verify hashes (minutes; exit code 1 on drift)
.venv\Scripts\python.exe scripts/v5/base/reproduce_base.py

# 2. Predict a blend
$env:PYTHONPATH="src"
.venv\Scripts\python.exe -c "from fuelblend.api.predict_route import predict_fuel; r=predict_fuel([('CCCCCCCCCCCC',.45),('CCCCCCCCCCCCCC',.20),('Cc1ccccc1',.15),('CC1CCCCC1C',.20)],T_K=288.15,fuel_type='jet'); print(r['flash_point_C']['value'], r['density_kg_m3']['value'])"

# 3. Inverse design on a pool
.venv\Scripts\python.exe scripts/v5/run_design_userpool_v5.py --spec astm:jet_a,bis:jet_a1 --pool v5pool --n-starts 30 --seed 42 --out results/v5/demo_design

# 4. Tests
.venv\Scripts\python.exe -m pytest tests/test_base_*.py -q     # fast guard of the frozen scoreboard
.venv\Scripts\python.exe -m pytest -q                            # everything (long)
```
Worked example output (4-component jet-like blend at 15 °C): flash 17.7 °C (`ideal_le_chatelier_anchored`, in domain), density 779.7 kg/m³ (`dens_idvol`), kinematic viscosity 1.65 cSt (`refutas_raw` — outside the ML domain), freezing −20.1 °C (`svl_unifac_mole_fraction`), net heat 43.7 MJ/kg, smoke 23.2 mm (not certified). Always read `in_domain`, `status` and `warnings` before trusting a number.

---

## 13. Likely questions from the faculty (and honest answers)

**Why not just train a neural network on lots of blends?** The only abundant "blend" data were synthetic labels from the blending rule; a net learns the rule (0.4 % error) which proves nothing. Real measured blends are few (hundreds of systems), so we use physics for the bulk and ML only for the pair correction.

**What exactly is learned?** One symmetric pair function g(i,j,T) (LightGBM, 400 trees) added to the Refutas log-viscosity weighted by x_i x_j. Features are descriptors of the two molecules, not the identities.

**How do you know it isn't over-fit?** Whole systems are held out; 70 locked systems were scored once; bootstrap CI over systems excludes zero (gain 0.0233, CI [0.0163, 0.0302]); the model is compared with a constant-offset comparator (0.043) too.

**What does 'certified' mean? Is it safe to fly?** No. It means: finite margin, in domain, and the spec is met after the margin — *under the model*. Nothing was measured; real certification needs ASTM tests.

**Why is freezing 'no skill' but you still use it?** Only through its lab-derived margin, for known molecules; on real fuels it's a safe-side screening estimate. We state that it does not beat a constant.

**Why did flash get so much better (−21.7 → +0.1 bias)?** The old API path used wrong routing (`api_as_shipped`); re-routing `predict()` onto the BASE route fixed it. The new numbers are descriptive on fuels used to develop the route.

**Why is your margin 2.2 °C but the yaml says 3.7?** 3.699 was derived on 31 groups before mesitylene/new systems; 2.207 (35 systems, mesitylene excluded) is applied. Mesitylene is disclosed as an outlier.

**Why is the Jet A-1 design above the AJFTD dicyclo range?** The optimiser is free to leave the observed real-fuel range; we flag it (37.8 % vs 30 % max) and do not claim realism.

**Why is the ML density only an option?** It failed the ternary gate (3.64 → 3.92 kg/m³).

**Where are the weaknesses?** Thin locked sets (3 and 4 cold viscosity systems), single GC fuel for D56, binary-only margins applied to multicomponent blends, no lab validation, reconstructed CRC compositions, Jet A-1 limits from secondary sources.

**What would you do next?** Obtain lab flash/freezing/viscosity on the three certified blends; collect ternary+ lab data for margins; verify Def Stan 91-091 limits; score the Oakley/Yang sets; test ML density on real fuels.

**Can you explain a file?** Use §11: start with `api/predict_route.py` (entry), `eval/composition_route.py` (what is computed), `baselines/flash.py` (145 lines, easiest to explain), `ml/viscosity_pair.py` (173 lines), `uq/conformal.py` (90 lines).

---

## 14. One-page equation sheet

```
x_i = (w_i/M_i)/Σ(w_j/M_j)            v_i = (w_i/ρ_i)/Σ(w_j/ρ_j)
1/ρ_mix = Σ w_i/ρ_i                   ρ_mix(ML) = Σ x_i M_i / [Σ x_i M_i/ρ_i + Σ_{i<j} x_i x_j g]
VBN_i = 14.534 ln(ln(ν_i+0.8)) + 10.975;  VBN_mix = Σ w_i VBN_i;  ν_mix = exp(exp((VBN_mix−10.975)/14.534)) − 0.8
ln η_mix = ln η_Refutas + Σ_{i<j} x_i x_j g(d_i,d_j,T);   ν = 1000 η/ρ
log10 P_sat = A − B/(T+C);   H(T) = Σ x_i P_sat,i(T)/(P_tot LFL_i) = 1  → T_flash
1/T_i = 1/T_m,i − R ln(x_i γ_i)/ΔH_fus,i;  T_freeze = max_i T_i
SP = 28.84 − 0.324 A − 0.791 D;   ΔH_c = Σ w_i ΔH_c,i
s_g = max_{p∈g}(error, dangerous side);  k = ⌈(n+1)(1−α)⌉;  margin = s_(k)  (∞ if k>n)
certified ⇔ margin finite ∧ in_domain ∧ (flash − margin ≥ 38 °C, freezing + margin ≤ limit, …)
```
