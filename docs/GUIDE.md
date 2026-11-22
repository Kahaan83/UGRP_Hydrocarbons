# UGRP Hydrocarbons: the complete guide

*Predicting jet-fuel blend properties from composition: physics first, machine learning only where it was proven to help, a measured safety margin on every number, and a clear label wherever the model has no skill.*

Authors: Kahaan Shah, Namyaa Parmar, Arijit Ganguli, A. S. K. Sinha, Ahmedabad University.

**Who this guide is for.** Anyone, including a reader who has never studied chemistry, fuels, statistics or machine learning. Part A builds the background from nothing. Parts 0 to 14 then explain every term, every decision, the current status, every important code file and how to run it. If you only have five minutes, read Part A section A9 and then section 0.

**Source of truth for numbers.** `docs/final_numbers_v5.md` (blends), `docs/final_numbers_v7.md` (pure compounds, best method per property) and `Fuel_Blend_Methodology.pdf`. If a number in this guide ever disagrees with those, they win. Worked examples in this guide that use made-up round numbers are marked **(illustration)** and are not results.

**Contents**

- Part A. Background from zero (A1 to A9)
- 0. The 90-second version
- 1. The problem and the goal
- 2. Glossary
- 3. The physics, property by property
- 4. The machine-learning story
- 5. Data: what we used and what we refused to use
- 6. Uncertainty: the conformal margin and "certified"
- 7. Honest evaluation rules
- 8. Inverse design
- 9. Where we are (status)
- 10. Decision log
- 11. Code map
- 12. How to run things
- 13. Likely questions and honest answers
- 14. Equation sheet
- 15. Reading the results folders (how to find any number)
- 16. What a prediction looks like end to end

---

# Part A. Background from zero

## A1. What is a fuel, and what is a hydrocarbon?

A **hydrocarbon** is a molecule built only from carbon (C) and hydrogen (H) atoms. Methane (natural gas) is one carbon with four hydrogens. The molecules in jet fuel are bigger: chains or rings of roughly 8 to 16 carbon atoms.

The kinds of hydrocarbon that matter here (called **classes**):

| Class | Shape | Example | Why it matters |
|---|---|---|---|
| **n-alkane** (straight-chain paraffin) | a straight chain, every bond single | n-decane (10 carbons), n-dodecane (12) | burns cleanly, high energy per kilogram, but long chains solidify (freeze) at fairly high temperature |
| **iso-alkane** (branched paraffin) | a chain with side branches | isooctane | branches stop the molecules stacking, so they freeze at lower temperature |
| **cycloalkane** (naphthene) | one or more saturated rings | butylcyclohexane, decalin | dense, good energy per litre, low freezing point |
| **aromatic** | contains a benzene ring | toluene, 1,2,4-trimethylbenzene, tetralin | dense and good solvents, but they make soot (low smoke point) |

**A real fuel is not one molecule.** Jet A is a blend of hundreds of different hydrocarbons produced by refining crude oil (or by newer routes such as hydroprocessed esters and fatty acids, HEFA, or alcohol-to-jet, ATJ). Its properties come from what is in the blend and in what proportions. The proportions are the **composition**.

## A2. Why predict properties from composition?

Fuels must satisfy a **specification** (spec): a list of limits written in standards such as ASTM D1655 for Jet A. A few of the limits we use:

| Property | Jet A limit | Jet A-1 limit | Why the limit exists |
|---|---|---|---|
| Flash point | at least 38 C | at least 38 C | below this the fuel gives off enough vapour to ignite near a spark, which is a fire risk in storage and handling |
| Freezing point | at most -40 C | at most -47 C | at high altitude the fuel is very cold; if wax crystals form they can block filters and fuel lines |
| Density at 15 C | 775 to 840 kg/m3 | 775 to 840 kg/m3 | energy content and fuel metering |
| Viscosity at -20 C | at most 8 cSt | at most 8 cSt | the fuel must still flow and atomise in the cold |
| Net heat of combustion | at least 42.8 MJ/kg | at least 42.8 MJ/kg | energy delivered per kilogram |
| Smoke point | at least 25 mm | at least 25 mm | how sooty the flame is |

Each of these is normally checked by a **laboratory test** on a physical sample. Lab tests cost time and money and need the fuel to exist already.

If we could **predict** the properties from the composition on a computer, we could (1) check a recipe before making it, and (2) run the process backwards (**inverse design**): ask the computer to search the space of recipes for one that meets the spec. That is the goal of this project.

## A3. The properties in plain words

- **Flash point.** Warm a liquid fuel slowly in a closed cup and hold a small flame near it. The flash point is the lowest temperature at which the vapour above the liquid ignites with a flash. Low flash point means more volatile and more dangerous. A prediction that is too high says "safe" when the fuel is not, so **over-prediction is the dangerous direction**.
- **Freezing point (jet fuel).** Cool the fuel until wax crystals first appear; then warm it until they disappear. That temperature is the freezing point. **Under-prediction is the dangerous direction**: the model says "stays liquid to -45 C" but wax appears at -40 C.
- **Density.** Mass per volume, in kg/m3. Hydrocarbon fuels are about 0.75 to 0.85 times the density of water.
- **Viscosity.** Resistance to flow (honey is more viscous than water). **Dynamic viscosity** eta (unit mPa.s) is the resistance to shear; **kinematic viscosity** nu = eta / density (unit cSt, centistokes, which equals mm2/s) is how quickly the fluid flows under gravity. Specs use kinematic. Water at 20 C is about 1 cSt. Viscosity rises steeply as the fuel gets colder, which is why the spec is stated at -20 C.
- **Heat of combustion.** Energy released by burning one kilogram, in MJ/kg. The "net" value assumes the water produced stays as vapour.
- **Smoke point.** The tallest flame, in millimetres, that burns without smoking. Aromatic molecules make soot, so more aromatics means a lower smoke point.
- **Distillation (D86).** Heat a sample and record the temperature at which 10 %, 50 %, 90 % has boiled off. It describes how spread out the boiling range is.
- **Octane and cetane.** Gasoline and diesel ignition-quality scales. Octane (RON/MON) is a research-only part of this project; cetane (DCN) was found not feasible.

## A4. Why a blend's property is not just the average

If you mix equal masses of two liquids with viscosities 1 cSt and 4 cSt, the blend is **not** 2.5 cSt. Viscosity mixing is strongly non-linear. The classical rule used here (Refutas, section 3.2) gives about **1.81 cSt** for that case **(illustration)**. Densities are close to a plain volume average (equal masses of 780 and 880 kg/m3 give about 827 kg/m3 **(illustration)**, not 830). Flash point is even less linear: a small amount of a volatile component can pull the flash point of the whole blend down a lot, because the vapour above the liquid is dominated by the most volatile molecules.

So each property needs its own **mixing rule**: a formula, derived from physics or chemistry, that takes the pure-component properties and the composition and returns the blend property.

## A5. What "machine learning" means here, in plain words

A **model** is a function that takes inputs and returns a prediction. **Machine learning** means the function has adjustable numbers (parameters) that are tuned automatically to fit a table of past examples (the **training data**). Two ideas matter for everything that follows:

1. **Generalisation.** A model is only useful if it works on cases it has not seen. Testing it on the same examples it was tuned on is meaningless, because a model can simply memorise them (**overfitting**). So we keep some data aside, never tune on it, and measure the error there (**held-out** or **test** data).
2. **Leakage and circularity.** If the test data are secretly related to the training data, the test is too easy. In this project the key danger is that several measurements belong to the same chemical mixture (the same two molecules at many ratios and temperatures). Those rows are near-copies. So we hold out **whole systems** (whole sets of molecules) rather than random rows.

The cautionary tale that started the project: an earlier approach trained a neural network on **synthetic** blends whose "measured" values were produced by the Refutas formula itself. It reported 0.4 % error, but only because it had learned to copy the formula. That says nothing about real fuels. The rule going forward: **real measurements, whole-system hold-out, comparison against a trivial baseline.**

## A6. The idea of the project in one paragraph

Use the physics and chemistry formulas as the backbone (they are tested, they extrapolate, and nothing in them is fitted to the blend data we later test on). Where measured blends show the formula leaves a systematic error that a small model can learn from molecule descriptors (for example, viscosity of two particular molecules mixing a little differently than ideal), add a **learned correction**, but accept it only if it beat the plain formula on locked held-out systems under a rule written down in advance. For every property, compute a **safety margin** from the held-out errors, so a prediction comes with "plus or minus this much, with 90 % confidence, in the dangerous direction". Finally, run inverse design with those margins and report which recipes pass. Be explicit about what is *not* proven.

## A7. The vocabulary of evidence

- **Trivial comparator (baseline).** The dumbest sensible predictor: "guess the average of the other fuels." Any model must beat this to claim it has skill. We use a leave-one-group-out version: for each test group, the guess is the average over all the other groups.
- **Verdict.** One of `beats_trivial`, `no_skill_shown`, `worse_than_trivial`, decided by comparing the model's error to the trivial one with a bootstrap over groups.
- **Gate.** A pass/fail rule committed to the repository **before** the result is computed, so it cannot be bent to fit the result. Named like R-MLV-1.
- **Locked set.** Systems set aside before modelling began and scored exactly once.
- **Evidence tiers.** T1 = lab mixtures of known pure molecules (strongest). T2 = real fuels whose composition was measured. T3 = bulk-property regressions on a large fuel database. CRC647_r3 = real fuels whose composition was reconstructed from reported hydrocarbon types (weaker).

## A8. How to read the status labels in this guide

- **Certified** never means "approved for flight". It means: the model has a finite margin for that property, the blend is inside the tested domain, and the spec is still met after the margin is applied. It is a statement about the model only. No design in this project has been measured in a laboratory.
- **No skill shown** is a result, not a failure to report. It means the model did not beat the trivial guess on real fuels, so we do not use it to claim anything.
- **Descriptive** means the fuels were used while building the route, so the number describes fit and is not a held-out claim.

## A9. The five-minute path

1. Read section 0 (90-second version).
2. Read A4 and A6 above if the logic is unclear.
3. Skim section 6 (margins) and section 8 (designs).
4. For questions, section 13.
5. For "where is that number", section 15.

---

# 0. The 90-second version (say this first)

1. **Problem.** A fuel is a blend of hundreds of molecules. Specs (Jet A: flash point at least 38 C, freezing point at most -40 C, density 775 to 840 kg/m3, viscosity at most 8 cSt at -20 C) are checked by lab tests. We want to predict these from the *composition* so we can design blends on a computer.
2. **Approach.** Every property is computed with a classical physics/chemistry rule fed with measured properties of the pure molecules. A small learned correction is added only where it beat the physics on mixtures the model never saw.
3. **Honesty.** Each model is tested on whole mixture *systems* held out from training, against a "just guess the average" baseline. A model that does not beat the guess is labelled *no skill shown* and is not used to certify.
4. **Safety margin.** From the held-out errors we compute a one-sided 90 % conformal margin. A blend is "certified" only if it still meets the spec after the margin is applied.
5. **Result.** Viscosity error 8 % to 3 % with the learned pair correction; flash point within about 2.2 C on 70 real-fuel groups; density within about 0.15 %. Freezing, smoke point and heat of combustion have **no proven skill on real fuels** and we say so. Inverse design finds 3 of 16 blend/spec combinations that pass **under the model only**; nothing was measured in a lab.
6. **Why it is credible.** Many negative results are kept (neural mixture models, GNN activity coefficients, ML sooting, ML cetane): they failed pre-declared tests and we left them in the repository on purpose.

---

# 1. The problem and the goal

**Goal.** Predict viscosity, density, flash point, freezing point, heat of combustion and smoke point of a hydrocarbon fuel blend (jet / diesel / RP-1 style) from its composition, with calibrated uncertainty, and then run **inverse design**: find a blend of available molecules that meets the specification.

**Why a plain neural network is wrong here.** Earlier work (the ChemCon abstract) trained a network on *synthetic* blends whose labels came from the blending rule itself (Refutas). It "achieved 0.4 % error", but only because it learned to copy the rule. That number says nothing about real fuels. The whole of this repository exists to avoid that circularity: **real measured data, physics baselines, whole-system hold-out, and a trivial-baseline comparison.**

**What "done" means for each property** (priority order): viscosity, flash point, density (core, certifiable); freezing point (certifiable only for known molecules with verified constants, margin from lab data); smoke point (regression, in range only); heat of combustion (reported, no margin); octane (gasoline research result, not in the API).

**What the project is *not*.** It is not a replacement for ASTM laboratory testing, it does not model engine performance, and it does not claim any blend is safe to use.

---

# 2. Glossary: every term you may be asked about

| Term | Meaning (plain words) |
|---|---|
| **Composition bases** | Mass fraction *w*, mole fraction *x*, volume fraction *v*. Different rules need different ones. x_i = (w_i/M_i)/sum(w_j/M_j); v_i = (w_i/rho_i)/sum(w_j/rho_j). A **mole** is a fixed count of molecules (6.022e23), so mole fraction is "what share of the molecules". M_i is molar mass in g/mol. Temperatures inside equations are always kelvin: T_K = T_C + 273.15. |
| **SMILES / InChIKey** | Text codes that identify a molecule exactly (SMILES = structure string such as `CCCCCCCCCC` for n-decane; InChIKey = fixed-length hash). The API accepts either. |
| **cSt (centistokes)** | Unit of *kinematic* viscosity = mm2/s. Water at 20 C is about 1 cSt. nu = 1000 eta / rho (eta in Pa.s, rho in kg/m3). |
| **Dynamic vs kinematic viscosity** | eta (mPa.s) is resistance to shear; nu = eta/rho is how fast the fluid flows under gravity. Specs use kinematic. |
| **Refutas** | Viscosity blending rule: map each nu to a Viscosity Blending Number (VBN) that averages linearly by mass. |
| **VBN** | VBN = 14.534 ln(ln(nu+0.8)) + 10.975. |
| **Grunberg-Nissan** | Classic rule: ln eta_mix = sum x_i ln eta_i + sum_{i<j} x_i x_j G_ij. Our ML learns the pair term. |
| **Ideal volume mixing** | Volumes add: 1/rho_mix = sum w_i/rho_i. Hydrocarbon mixtures deviate by well under 1 %. |
| **Excess volume** | The (small) volume change on mixing. ML option learns it per pair of molecules. |
| **Antoine equation** | Vapour pressure of a pure liquid: log10 P_sat = A - B/(T + C). Here P in kPa, T in K. A, B, C are constants fitted per molecule. |
| **Vapour pressure** | The pressure of vapour in equilibrium with the liquid; high means volatile. |
| **Raoult's law** | Ideal mixture: partial pressure = x_i * P_sat,i. |
| **LFL** | Lower flammability limit: smallest fuel-vapour fraction in air that can burn (mole fraction). |
| **Le Chatelier's rule** | A vapour mixture is flammable when sum y_i/LFL_i = 1. |
| **Flash point** | Lowest temperature at which the vapour above the liquid ignites. Spec: at least 38 C for Jet A. Over-predicting is the *dangerous* direction. |
| **Zabetakis correction** | LFL falls with temperature: LFL(T) = LFL_25 (1 - 0.000721 (T - 25)), T in C. |
| **Anchored Antoine** | Antoine constants re-anchored to measured vapour pressure in the relevant window (0 to 100 C) because literature ranges for C10+ alkanes start at 95 to 190 C. |
| **Registry overlay** | Extra Antoine + LFL entries (bicyclohexyl, cis-decalin) so real-fuel compositions mapped to those molecules can be flash-predicted. On by default; `FUELBLEND_FLASH_REGISTRY_ADDITIONS=0` turns it off. |
| **Freezing point (jet)** | Temperature at which the first wax crystals appear. Spec: at most -40 C (Jet A), at most -47 C (Jet A-1). Under-predicting is dangerous. |
| **Schroeder-van Laar (SvL)** | Solid-liquid equilibrium: 1/T_i = 1/T_m,i - R ln(x_i gamma_i)/dH_fus,i; the blend freezes at max_i T_i. |
| **UNIFAC** | Group-contribution method that estimates activity coefficient gamma_i from chemical groups. Main freezing route; ideal (gamma = 1) is the variant. |
| **Activity coefficient gamma** | How far a molecule departs from ideal behaviour in a liquid mixture (gamma = 1 means ideal). |
| **dH_fus, T_m** | Heat of fusion and melting point of the pure solid. |
| **Smoke point** | Tallest flame (mm) that burns without smoking; aromatics lower it. Regression SP = 28.84 - 0.324 A - 0.791 D (A aromatics vol %, D diaromatics vol %). |
| **YSI** | Yield sooting index, a measured soot-tendency scale per molecule. Our YSI route failed its gate. |
| **Heat of combustion (net)** | Energy per kg released burning the fuel (MJ/kg). Group additivity: sum w_i dH_c,i. Spec at least 42.8. |
| **Group additivity** | Estimate a molecule's property by adding contributions of its structural pieces (CH3, CH2, ring, ...). |
| **D86 / Engler** | Distillation curve (T10, T50, T90...). Engler = our default TBP-offset model; Rayleigh = opt-in alternative that failed the adoption rule. |
| **DCN / cetane** | Diesel ignition quality. ML version "not feasible" (only 31 usable systems). |
| **RON / MON / Co-Optima** | Gasoline octane measures; Co-Optima = DOE dataset of 65 blends. Research only. |
| **Simplex** | The set of non-negative fractions summing to 1; the search space of inverse design. |
| **System** | One set of molecules (for example decane + toluene) measured at any ratio/temperature. All its points stay together in train **or** test. |
| **GroupKFold** | Cross-validation that never splits a system across train/test. Cross-validation = split the data into folds, hold each out in turn. |
| **Locked (held-out) set** | Systems set aside *before* any model was built and scored exactly once. |
| **Gate (R-xxx-n)** | A pre-declared pass/fail rule, committed to git before results exist (for example R-MLV-1). |
| **Trivial comparator** | "Predict the average of the other groups" (leave-one-group-out constant). A model has skill only if it beats this. |
| **Verdict** | `beats_trivial`, `no_skill_shown`, `worse_than_trivial` (group-bootstrap comparison vs trivial). |
| **Bootstrap** | Resample the groups with replacement many times to see how much a statistic wobbles; gives a confidence interval (CI). |
| **Bias / MAE** | Bias = mean signed error (pred - meas). MAE = mean absolute error. |
| **ln-error for viscosity** | abs(ln eta_pred - ln eta_meas); for small values about the relative error (0.031 is about 3 %). |
| **Conformal margin** | Distribution-free allowance from ranked held-out errors (see section 6). |
| **Certified** | Margin exists AND blend in tested domain AND spec met after margin. *Model-level only.* |
| **Domain / `in_domain`** | Known molecules with their own data, 233 to 393 K, class pairs seen at least 5 times in training; no stand-ins. |
| **Stand-in** | A close relative molecule used when the exact one has no data; flagged, never certified. |
| **T1 / T2 / T3 / CRC647_r3** | Evidence tiers: T1 lab mixtures; T2 real fuels with measured composition; T3 AJFTD bulk-property regressions; CRC647_r3 real fuels with *reconstructed* composition. |
| **AJFTD** | Alternative Jet Fuels Test Database (2,601 fuels scraped; GCxGC composition count being updated by the AJFTD rebuild: 132 with bins, 126 complete). |
| **GCxGC** | Two-dimensional gas chromatography: composition by molecule class x carbon number. |
| **CRC 647** | Coordinating Research Council report: 57 aviation fuels with property tables. |
| **ThermoML** | NIST archive of measured thermophysical data (density, viscosity, vapour pressure...). |
| **Leenhouts benchmark** | Curated pure + mixture flash/viscosity/DCN/MON benchmark (Zenodo 19914760). |
| **LightGBM** | Gradient-boosted decision trees (many small decision trees added together); the learner for pair corrections. |
| **DeepSets / GNN / set transformer** | Neural architectures for mixtures (a GNN works on the molecular graph); all failed our gates. |
| **BASE / delta_*** | `results/v5/BASE/` is the frozen scoreboard; each `delta_*` folder is one experiment measured against it. |
| **CLM-nn** | Claim number in `docs/claims_register.md`. |
| **Hash / SHA-256** | A fingerprint of a file's contents; if one byte changes the fingerprint changes. Used to detect silent changes to frozen results. |

---

# 3. The physics, property by property

All of this is implemented once, in the **designated scoring route** (`src/fuelblend/eval/composition_route.py`) and exposed through `src/fuelblend/api/predict_route.py::predict_fuel`.

## 3.0 The general recipe for every property

1. Look up the property of each **pure** molecule at the temperature of interest (from the frozen table, section 5).
2. Convert the composition to the basis the rule needs (mass, mole or volume fractions).
3. Apply the mixing rule.
4. Optionally add the learned pair correction (viscosity, density) if it passed its gates and the blend is in its domain.
5. Attach the margin from section 6 and the status flags.

## 3.1 Density: ideal volume mixing

`1/rho_mix = sum w_i/rho_i` (all rho at the same T). Code: `mixture/common.py::density_idvol`.

*In plain words.* Take 1 kg of blend. The mass of component i is w_i kg; its volume is w_i / rho_i. Add the volumes, and divide the total mass (1 kg) by the total volume. This assumes mixing neither shrinks nor expands the liquid; for hydrocarbons that is true to well under 1 %.

*Worked example (illustration).* 50 % by mass of a 780 kg/m3 liquid and 50 % of an 880 kg/m3 liquid: 1/rho = 0.5/780 + 0.5/880 = 0.0012093, so rho = 827.0 kg/m3.

Evidence: lab 1.33 kg/m3 MAE (about 0.17 %); real fuels beat the trivial guess. **Reconstructed-composition CRC fuels do not** (14.3 vs 13.2 kg/m3); that is a composition-reconstruction limit, not a density-law failure.

**ML option** (`ml/density_pair.py`): `V = sum x_i M_i/rho_i + sum_{i<j} x_i x_j g(d_i,d_j,T)`, `rho = sum x_i M_i / V`. Here V is the molar volume and g is a learned pair term for the excess volume, with d_i the descriptors of molecule i. 1.20 to 0.90 kg/m3 on 41 held-out binary systems, margin +/-4.42. **Fails on ternaries** (3.64 to 3.92), so it is a labelled variant for known binaries only; ideal volume stays the main value.

## 3.2 Viscosity: Refutas, plus learned pair correction

Convert each pure nu_i to VBN_i = 14.534 ln(ln(nu_i+0.8)) + 10.975; average by *mass*: VBN_mix = sum w_i VBN_i; invert: nu_mix = exp(exp((VBN_mix - 10.975)/14.534)) - 0.8. Valid for nu of about 0.2 cSt and above. Code: `mixture/common.py::refutas_ln_eta`, `eval/composition_route.py::refutas_nu_safe`.

**Why the double logarithm.** Viscosities of blend components differ by orders of magnitude; a mass-weighted average of nu would be dominated by the thick component. The double log compresses the scale so that mass-averaging works (this is an empirical rule that works well for petroleum fractions).

*Worked example (illustration).* 50/50 by mass of a 1.0 cSt and a 4.0 cSt liquid: VBN = 3.252 and 17.518; mix = 10.385; nu_mix = exp(exp((10.385 - 10.975)/14.534)) - 0.8 = 1.81 cSt. The plain average would be 2.5 cSt.

**Pure-component data:** a *frozen* table (`data/cache/pure_grid.npz`, SHA-256 pinned and checked at every load) of ln eta versus T for 1,042 compounds (114 with a full grid), starting at 230 K. A close relative may stand in if a molecule has no data, flagged and never certified.

**ML pair correction** (`ml/viscosity_pair.py`, models in `models/v5/ml_viscosity_pair`):
`ln eta_mix = ln eta_Refutas + sum_{i<j} x_i x_j g(d_i, d_j, T)`.

*In plain words.* The physics rule is right on average but misses how particular pairs of molecules interact. For every pair of molecules in the blend, a small model predicts a correction from *descriptions* of the two molecules (not their names). The correction is weighted by x_i x_j, so it vanishes for a pure compound and grows when both are present in quantity. Because only pairs are learned, the same model handles blends of any number of components.

LightGBM (400 trees, 15 leaves, learning rate 0.05, seed 42, deterministic) on *pair features*: mean and absolute difference of carbon number, molar mass, ring counts, CH3 count, quaternary carbon count, ln eta, molar volume, the class pair and 1/T. Symmetric in the pair, zero for a pure compound. Trained on 264 binary systems, graded once on 70 locked systems.

Gate results (Refutas to ML, MAE in ln eta): R-MLV-1 dev 0.082 to 0.060; R-MLV-2 locked 0.054 to 0.031; R-MLV-3 ternary+ 0.098 to 0.042 (per-system); R-MLV-4 cold (below 283 K) 0.099 to 0.061. All passed.

**Designation rule (user decision):** ML is the main viscosity *inside its domain* (known molecules, 283 K and above for scoring, class pairs with at least 5 training systems). Elsewhere the main value is Refutas with a warning.

**Honest limit:** on real fuels below 263 K Refutas is *worse than a constant* (MAE 0.53 ln vs 0.25), so no real-fuel cold margin exists (margin = infinity).

## 3.3 Flash point: Le Chatelier + Raoult + Antoine (nothing fitted to blends)

*The physical picture.* A closed container of liquid fuel fills its headspace with vapour. If the vapour/air mixture contains enough fuel (above the lower flammability limit, LFL) a spark ignites it. As the liquid warms, more vapour forms, so there is a temperature at which the vapour just reaches the LFL: the flash point.

1. Pure vapour pressure: `log10 P_i^sat = A_i - B_i/(T + C_i)` (kPa, K), constants anchored to measured vapour pressure in the 0 to 100 C window.
2. Raoult (ideal): the partial pressure of component i above the blend is `x_i P_i^sat(T)`; its mole fraction in the vapour is `y_i = x_i P_i^sat(T)/P_total`, P_total = 101.325 kPa.
3. Le Chatelier: the vapour is flammable when `H(T) = sum y_i/LFL_i = sum x_i P_i^sat(T)/(P_total LFL_i) = 1`. H rises with T, so solve H(T) = 1 for T with the root finder `brentq` on 150 to 650 K. LFL is corrected to temperature (Zabetakis).

Code: `baselines/flash.py` (`psat_antoine`, `lfl_temperature_correction`, `hazard_function`, `predict_flash_le_chatelier`).

**Why this is the fair test:** nothing is fitted to blend flash points, so testing on blends is honest.

Evidence: 70 real-fuel groups MAE 2.20 C (trivial 5.16), bias +0.06 C; AJFTD GCxGC fuels (31) bias -0.15, MAE 3.00, beats trivial; lab mixtures margin 2.207 C (35 systems).

Applied flash margin is **2.207 C** (the yaml value 3.699 and the scoreboard `margin90` 3.699 are the superseded pre-Q2 figures).

Limits: the D56 real-fuel check rests on **one** GC fuel (POSF-10325), so no margin; the ideal-mixture assumption fails for some systems (a mesitylene mixture is an excluded outlier, disclosed).

## 3.4 Freezing point: Schroeder-van Laar with UNIFAC

*The physical picture.* Each component, if it were alone, would start to crystallise out of the liquid at a temperature set by its melting point, its heat of fusion and how much of it is present (its mole fraction, and how "happy" it is in the mixture, gamma). The blend begins to freeze at the *highest* of these temperatures.

`1/T_i = 1/T_m,i - R ln(x_i gamma_i)/dH_fus,i`; `T_freeze = max_i T_i`. (R = gas constant.) Main route: gamma from original UNIFAC (`baselines/flash_unifac.py` and the SLE code in `eval/composition_route.py::freezing_C`, `DESIGNATED_FREEZING_MODEL = "svl_unifac"`). Variant: ideal (gamma=1). Melting points / dH_fus come from a *verified* constants table (`data/raw/sle_constants/tm_dhfus_verified.csv`); bicyclohexyl uses dH_fus plus its solid-solid transition (13.94 kJ/mol).

Lab liquidus set (73 points, 12 systems): UNIFAC MAE 3.93 C, bias +2.47, **margin 13.84 C**; ideal 4.47 / +1.17 / 16.98; constant guess 15.63.

**Real fuels: no skill.** Pooled 50 fuels MAE 8.42 vs constant 6.41; locked 8 fuels 9.31 vs 5.18 (`worse_than_trivial`). Bias is +4 to +6 C (safe side), so the 13.84 margin "covers" 98 % trivially. We use freezing *only* through the lab margin, for known molecules with verified constants.

Known data issue: the o-/m-xylene eutectic composition appears transposed in the ingest script (not fixed; margins unaffected).

## 3.5 Smoke point and heat of combustion

Smoke: `SP = 28.84 - 0.324 A - 0.791 D` (fitted on 100 jet fuels; trusted for 8 to 23 vol % aromatics and 19 to 25 mm). A and D are the bulk aromatic and diaromatic volume percentages (a bulk measurement, not a molecule list). 0.90 mm MAE on the 100-fuel set (beats trivial), no skill shown on the broader set (1.68 vs 1.65). YSI route failed gate (median relative error 0.229 vs 0.20).

Heat: group additivity `sum w_i dH_c,i`; reported, **no margin**; worse than trivial on reconstructed CRC fuels.

## 3.6 Octane (research only)

Blending octane value: `RON = (1-v) RON_base + v BOV` (v volume fraction, BOV the blending octane value of the added stock). Gaussian process / GBM on blendstock descriptors, 65 Co-Optima blends, 18 dev blendstocks + 4 held out: RON 2.29, MON 1.15 (CV). Too few held out to certify; not in the API.

---

# 4. The machine-learning story

**Principle:** ML learns only what the ideal rule misses (pairwise molecular interactions) and is used only where a pre-declared gate passed on locked systems.

| Model | Base | Held-out result | Status |
|---|---|---|---|
| Viscosity pair correction (LightGBM) | Refutas | 0.054 to 0.031 ln (70 systems) | **main value in domain** |
| Density excess volume (LightGBM) | ideal volume | 1.20 to 0.90 kg/m3 (41 systems); ternaries fail | option (known binaries) |
| Octane GP | BOV rule | RON 2.29 / MON 1.15 (CV) | research only |
| Neural mixture encoders (DeepSets, GNN, set transformer) | Refutas / ideal volume | none pass; GNN worse (0.031 to 0.040) | **failed, kept** |
| GH-GNN activity coefficients in flash | ideal Le Chatelier | +0.70 C worse lab, +1.09 C worse real fuels | **failed, kept** |
| GCxGC real-fuel models | physics + offset | worse in 10 of 12 leave-one-source-out cases | **does not transfer** |
| ML YSI for smoke | class mean | rel. error 0.229 > 0.20 | **failed** |
| ML DCN | none | only 31 usable systems (about 40 needed) | **not feasible** |
| Antoine interpolation, Tm/dH_fus estimators, multicomponent freezing | none | fail declared gates | failed, flagged unused |
| Pure-compound property models (v7): random forest on RDKit descriptors, D-MPNN | constant | MAE per property in `docs/final_numbers_v7.md` (dev rows, optimistic) | **deployed in `predict_pure`** for 10 pure-compound properties |
| v7 blend candidates: HANNA gamma, parameter-head Antoine, MolPool, set models | v5 route | none beats v5 with a CI above 0 | **not adopted, kept** |

**Process every model followed:** (1) pass rules committed to git first; (2) development by cross-validation over whole systems; (3) locked systems scored once; (4) bootstrap CI over systems; (5) compare to trivial constant; (6) failed models stay in code and in `docs/v5_negative_results.md`.

**Why keep failures?** A reader can then see what was tried and that the accepted models were chosen by rule rather than cherry-picked, and nobody wastes time repeating a dead end.

---

# 5. Data: what we used and what we refused to use

| Source | What | Notes |
|---|---|---|
| Leenhouts benchmark (Zenodo 19914760) | curated pure + mixture flash, viscosity (ln mPa.s), DCN, MON | primary lab mixture evidence |
| NIST/TDE viscosity (CheMixHub, Bilodeau) | about 239k binary viscosity rows | value column is ln(eta/Pa.s) despite label |
| ThermoML archive | measured density, viscosity, vapour pressure | re-extracted after parser fixes (phase, constraint, units) |
| NIST WebBook | Antoine constants (bar), converted to kPa (+2 in log10) | ranges do not reach 40 to 100 C for C10+ alkanes, so anchored/refit |
| Zabetakis Bulletin 627 | measured LFL values | hand-read, ratio-gated |
| Li et al. (Sun/Krakauer/Morgan 2020) | pure flash points | loader `load_flashpoint_li` |
| AJFTD | 2,601 real fuels; GCxGC count being updated by the AJFTD rebuild (132 with bins, 126 complete) | login-gated, not committed |
| CRC 647 | 57 jet fuels | compositions *reconstructed* from measured hydrocarbon types (Table G.1) |
| NISTIR 6646 / Fortin 2012 | RP-1 characterisation, distillation | RP-1 11-blend recipes confidential |
| Co-Optima | 65 octane blends | research |

**Refused:** `DatasetsFuelProperties.xlsx` (50,000 *synthetic* rows) and the v3 fuel corpus + optimiser (formula-derived; 92 of 153 fuels dropped, quarantined in git history). **Rule:** when a result looks too good, look for circular labelling.

**De-duplication:** drop exact repeats, group by (property, system, composition, T, P), take the median, flag (never silently delete) outliers and conflicts; the group key for cross-validation is `system_id` (a hash of the sorted set of the components' InChIKeys, so the same molecules always form the same group).

**Unit traps that were fixed (do not reintroduce):** NIST "log viscosity / cP" is really ln(Pa.s); Leenhouts viscosity is ln(mPa.s) (difference 6.908); WebBook Antoine is in bar; ThermoML uses Pa.s and m2/s; the flash unit string is spelled `Celcius` in the source.

---

# 6. Uncertainty: the conformal margin and "certified"

*The idea in plain words.* Suppose we have n held-out test systems and, for each, the model's worst mistake in the dangerous direction. If we sort those n mistakes from smallest to largest and take one near the top, then a new system drawn the same way will have a smaller mistake than that value about 90 % of the time. That value is the **margin**. It needs no assumption about the shape of the error distribution, which is why it is called distribution-free (conformal).

For each held-out system *g* take the **worst dangerous-side error** `s_g = max_{points in g}(error in the dangerous direction)`. With *n* systems and alpha = 0.10:
`k = ceil((n+1)(1-alpha))`, `margin = s_(k)` (the k-th smallest score). If *k > n* the margin is **infinity**: not enough systems (needs at least 9 at 90 %). Code: `src/fuelblend/uq/conformal.py` (`conformal_rank`, `group_scores`, `one_sided_margin`, `two_sided_margin`, `coverage`). Applied margins live in `configs/uncertainty_v5.yaml`.

*Worked rank examples (these follow directly from the formula).* n = 8: k = ceil(9 x 0.9) = 9 > 8, so margin = infinity. n = 9: k = 9, the largest score. n = 12 (the freezing set): k = ceil(13 x 0.9) = 12, the largest of the 12 scores. n = 35 (flash): k = 33.

| Property | Dangerous side | Margin | Basis |
|---|---|---|---|
| Flash | too high | **2.207 C** | 35 lab systems |
| Freezing (UNIFAC) | too low | **13.84 C** | 12 lab systems |
| Freezing (ideal, variant) | too low | 16.98 C | 12 |
| Density (ideal) | both | +4.59 / -3.25 kg/m3 | 167 |
| Density (ML, binaries) | both | +/-4.42 kg/m3 | 41 locked |
| Viscosity ML 283 K and above | both | +/-0.134 ln (about +/-14 %) | 264 |
| Viscosity ML 263 to 283 K | both | +/-0.608 ln | 16 |
| Viscosity ML below 263 K | both | +/-0.287 ln | 20 |
| Viscosity Refutas 283 K and above | both | +/-0.186 ln | 255 |
| Smoke (jet class mix) | too high | 4.11 mm | 7 source groups |
| Heat | n/a | none | n/a |

**Worked example.** Jet A needs flash at least 38 C; margin 2.207 means certified only if predicted flash is at least 40.207 C. "Predicted 38.1 C" is *not* certified.

**Reading a ln margin.** A viscosity margin of 0.134 in ln means the true value lies within a factor exp(0.134) = 1.14 either way, about plus or minus 14 %.

**Certified =** (a) a finite margin exists, (b) the blend is `in_domain`, (c) the limit still holds after the margin. Otherwise `meets_at_point_estimate` or `fails`.

**Important:** the freezing and cold-viscosity margins come from mostly **binary** lab systems; applying them to 6-component designs is a transfer, not a validation.

---

# 7. Honest evaluation: the rules that make the numbers trustworthy

1. **Real data only**; synthetic labels banned.
2. **Hold out whole systems** (GroupKFold on `system_id`; leave-one-source-out for real fuels).
3. **Pre-declare** every test (README with the pass rule committed *before* the result).
4. **Beat a trivial guess** (leave-one-group-out constant) or the verdict is `no_skill_shown`.
5. **Bootstrap over groups** for CIs.
6. **Keep negative results**; never delete a failed model.
7. **Descriptive vs gated:** T2 real-fuel cells are *descriptive* (every T2 fuel has been seen), not held-out claims.
8. **Frozen BASE:** `results/v5/BASE/` is hash-protected by `SHA256SUMS_v5_base.txt`; any change is a visible drift.

---

# 8. Inverse design

*Plain words.* Instead of asking "what are the properties of this recipe?", ask "which recipe has the properties I want?" The recipe is a list of fractions, one per available molecule. The set of all valid recipes (fractions at least zero, sum one) is the **simplex**. A search algorithm tries many recipes, calling the same prediction code each time, and keeps the best.

`design/engine_v4.py` searches the simplex with multistart SLSQP-style optimisation (10 starts, seed 42), maximising the minimum slack across constraints (slack = how far inside the limit the margin-adjusted prediction sits), then **re-predicts the result with the same API** and judges limit by limit, so design and prediction cannot disagree. Constraints come from config only (specs and margins). `strict_constraints=True` (default) makes the engine stop with `ModelUnavailable` instead of silently dropping a constraint.

Pools: P0 = 7 pure molecules (n-decane, n-undecane, n-dodecane, butylcyclohexane, trans-decalin, 1,2,4-trimethylbenzene, isooctane); P0x-a (+tetralin, n-tridecane); P0x-b (+hexylbenzene); P0y-b (+bicyclohexyl); P0y (+bicyclohexyl, cis-decalin); P1 to P3 mix real blendstocks (GCxGC Jet A base, HEFA, ATJ, SAK) with pure molecules.

**Result (8 pools x Jet A / Jet A-1 = 16 designs):** certified under the model: P0y Jet A, P0y-b Jet A, P0y-b Jet A-1. All three carry labels `d86_model_dependent`, `freezing_route_dependent`, `overlay_dependent`: they pass on Engler D86 but fail the Rayleigh check; they pass freezing on UNIFAC+13.84 but would exceed the limit on ideal SvL+16.98; they need the bicyclohexyl/cis-decalin registry additions. The certified Jet A-1 has 37.8 wt % dicycloparaffin, above the AJFTD GCxGC maximum (30 %). Minimum slacks (0.003 to 0.013) show the optimiser sitting on the margin boundary, **not robustness**. No design has been measured; a lab test was ruled out (decision 2026-11-14).

---

# 9. Where we are (status, 22 Nov 2026)

* **Pure-compound layer (v7):** `predict_pure` gives the best-evidenced method for ten pure-compound properties (random forest for flash point, cetane, heat of fusion, log viscosity, melting point, YSI; D-MPNN for derived cetane, density, MON, RON). Significant wins: flash point, heat of fusion, log viscosity. The rest are tied or nominal. Errors are dev-row CV and optimistic; no conformal margin. Blend routes are unchanged because no v7 candidate beat v5.
* **Working and evidenced:** physics routes for all properties; viscosity and density ML corrections (lab-gated); conformal margins; 16-design table; frozen, hash-checked scoreboard; API accepting SMILES/InChIKeys.
* **No proven skill on real fuels:** freezing, heat of combustion (reconstructed), pour point, cold viscosity (below 263 K, Refutas), CRC reconstructed density.
* **Not done / not claimed:** any laboratory validation; diesel validation; multicomponent freezing margin; real-fuel ML density testing; Jet A-1 limits verified against Def Stan 91-091; lab freezing of hexadecane mixtures (Oakley 2017) and smoke of 23 coal/JP-8 fuels (Yang 2006) are in the data but unscored.
* **Tests:** last full run 785 passed / 0 failed (before the simplification); after simplification the BASE suite and the touched suites were re-run (see the repository's final commit message). The full suite has not been re-run on the simplified branch since.
* **Manuscript:** paused by the user; only `manuscript/refs.bib` (with unverified entries) and figures exist in the deliverables folder.

---

# 10. Decision log: what we decided and why

| Decision | Why |
|---|---|
| Physics-first, ML as residual | blending rules are right to about 1 % for density and a few % for viscosity; a net learning the rule proves nothing |
| Hold out systems, not rows | rows of one system are near-copies; row-level CV hides over-fitting |
| One-sided conformal margin | specs are one-sided; direction of the error matters for safety |
| Margin infinity when n too small | refuse to invent a number the data cannot support |
| ML viscosity designated in domain | passed all 4 gates; Refutas stays as `refutas_raw` outside the domain |
| ML density is only a variant | failed the ternary gate |
| UNIFAC-SvL designated for freezing | smaller margin (13.84 vs 16.98), smaller error, no case it could not compute (declared rule) |
| Flash-registry overlay on by default | enables bicyclohexyl/cis-decalin; `FUELBLEND_FLASH_REGISTRY_ADDITIONS=0` restores old |
| Engler D86 kept; Rayleigh opt-in | Rayleigh T10 bias +11.6 C failed the adoption rule |
| Smoke = bulk regression, not YSI | YSI had 0 evaluable fuels / failed gate (user: "pick the less-error one") |
| No lab test | ruled out by the user; certification stays model-level and says so |
| Keep failed models | reproducibility of negative results is part of the science |
| Frozen BASE + hash manifest | any silent change to evidence is detected |
| v3 corpus quarantined | its fuel compositions and optimiser were derived from the formulas being tested |

History in one line: v1 data pipeline, then v2/v3 (circular, quarantined), then v4 integrity rebuild on real fuels, then v5 full baseline, ML components, adoption of user decisions, final numbers pack. The full pre-simplification state is the git tag `v5-full-2026-11-16`.

---

# 11. Code map: which file does what

**Entry points**

| File | Role |
|---|---|
| `src/fuelblend/api/predict_route.py` | `predict_fuel(composition, T_K, fuel_type)`: the public prediction function; loads the frozen grid (`GRID_SHA256` check), resolves names/SMILES/InChIKeys, returns each property with `value`, `method`, `status`, `in_domain`, `lo/hi`, `variants`, `warnings` |
| `src/fuelblend/api/best_methods.py` | registry: best method per property for pure compounds and blends, with the committed CSV and column behind each error; `best_method(prop, scope)` |
| `src/fuelblend/api/predict_pure.py` | `predict_pure(smiles, props)`: best-method predictions for pure compounds; RF fits on first use into `data/cache/pure_models/`, D-MPNN weights come from `scripts/v7/final/fit_pure_models.py` |
| `src/fuelblend/api/predict_blend.py` | `predict(...)`: v4-compatible wrapper on the same route; reads margins from the yaml (`_load_margins`, `_applied_margin`) |
| `scripts/v5/base/reproduce_base.py` | rebuilds the scoreboard, runs the BASE tests, checks every file in `results/v5/BASE` against `SHA256SUMS_v5_base.txt` |
| `scripts/v5/run_design_userpool_v5.py`, `run_jet_designs_v5.py` | inverse design on a user pool / the jet pools |

**Physics**

| File | Role |
|---|---|
| `baselines/flash.py` | Antoine, LFL temperature correction, hazard function, `brentq` flash solver (145 lines, the easiest file to read) |
| `baselines/flash_unifac.py` | UNIFAC activity coefficients (variant / freezing) |
| `baselines/cold_flow.py`, `mixture/cold_flow.py` | freezing/cloud helpers |
| `baselines/heat_of_combustion*.py`, `sooting.py`, `distillation*.py`, `cetane.py`, `octane.py` | heat, smoke/YSI, D86 (Engler/Rayleigh), DCN, octane |
| `mixture/common.py` | `density_idvol`, `refutas_ln_eta`, torch Newton flash, `CompoundRegistry`, `PurePropertyProvider` |
| `eval/composition_route.py` | the **designated scoring route**: `refutas_nu_safe`, `freezing_C` (`ideal_svl` / `svl_unifac`), `RouteCtx`, `predict_properties` |

**Machine learning**

| File | Role |
|---|---|
| `ml/viscosity_pair.py` | LightGBM pair correction; `FEATURES`, `LGB_PARAMS`, `ln_correction` (= 0.5 x^T G x), `domain_status` |
| `ml/density_pair.py` | same design for excess volume |
| `ml/set_models.py` | DeepSets / GNN / set transformer (failed, kept) |
| `ml/ysi_model.py` | ML YSI (failed, kept) |
| `mixture/t8_deepsets.py`, `t8_pairwise.py`, `set_encoder_v4.py`, `structured_residual.py`, `bilinear_pairwise.py` | earlier mixture-ML experiments |
| `models/v5/` | trained model files |

**Uncertainty and evaluation**

| File | Role |
|---|---|
| `uq/conformal.py` | conformal rank/margins/coverage (90 lines) |
| `eval/cellkit.py`, `group_stats.py`, `metrics.py` | scoreboard cells, group bootstrap, verdict logic |
| `src/fuelblend/splits/locked_sets.py`, `splits/` (repo root) | locked-set definitions |
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
| `data/unified/` | `measurements.csv`, `measurements_dedup.csv`, `compounds.csv` (schema in `AGENTS.md` section 6) |
| `data/real_fuels/v4/`, `v5/` | AJFTD / CRC real-fuel tables |
| `data/cache/pure_grid.npz` | frozen pure-property grid |

**Evidence and docs**

| Path | Role |
|---|---|
| `results/v5/BASE/` | frozen scoreboard + one `delta_*` folder per experiment (each has a README declaring its readings before results) |
| `results/v5/final/` | tables/figures behind `docs/final_numbers_v5.md` |
| `docs/claims_register.md` | every claim CLM-xx and its status |
| `docs/model_decision.md` | model choices and ML component decisions |
| `docs/v5_negative_results.md` | the failed attempts |
| `PROJECT_HANDOFF.md` (repo root) | state of the project for whoever continues |

---

# 12. How to run things

```powershell
# one-time
.venv\Scripts\python.exe -m pip install -r requirements.txt

# 1. Rebuild the frozen scoreboard and verify hashes (minutes; exit code 1 on drift)
.venv\Scripts\python.exe scripts/v5/base/reproduce_base.py

# 2. Predict a blend
$env:PYTHONPATH="src"
.venv\Scripts\python.exe -c "from fuelblend.api.predict_route import predict_fuel; r=predict_fuel([('CCCCCCCCCCCC',.45),('CCCCCCCCCCCCCC',.20),('Cc1ccccc1',.15),('CC1CCCCC1C',.20)],T_K=288.15,fuel_type='jet'); print(r['flash_point_C']['value'], r['density_kg_m3']['value'])"

# 2b. Predict pure-compound properties with the best method (random forests fit on first use)
.venv\Scripts\python.exe -c "from fuelblend.api.predict_pure import predict_pure; print(predict_pure('CCCCCCCCCCCC', props=['flash_C','tm_K']))"
# fit the D-MPNN properties (derived cetane, density, MON, RON) in the GPU environment
$env:FUELBLEND_DEVICE='cuda'; C:uelblend-gpu\.venv\Scripts\python.exe scripts/v7/final/fit_pure_models.py

# 3. Inverse design on a pool
.venv\Scripts\python.exe scripts/v5/run_design_userpool_v5.py --spec astm:jet_a,bis:jet_a1 --pool v5pool --n-starts 30 --seed 42 --out results/v5/demo_design

# 4. Tests
.venv\Scripts\python.exe -m pytest tests/test_base_*.py -q     # fast guard of the frozen scoreboard
.venv\Scripts\python.exe -m pytest -q                            # everything (long)
```

Worked example output (4-component jet-like blend at 15 C: 45 % n-dodecane, 20 % n-tetradecane, 15 % toluene, 20 % 1,2-dimethylcyclohexane, by mass): flash 17.7 C (`ideal_le_chatelier_anchored`, in domain), density 779.7 kg/m3 (`dens_idvol`), kinematic viscosity 1.65 cSt (`refutas_raw`, outside the ML domain), freezing -20.1 C (`svl_unifac_mole_fraction`), net heat 43.7 MJ/kg, smoke 23.2 mm (not certified). Always read `in_domain`, `status` and `warnings` before trusting a number.

---

# 13. Likely questions from the faculty (and honest answers)

**Why not just train a neural network on lots of blends?** The only abundant "blend" data were synthetic labels from the blending rule; a net learns the rule (0.4 % error) which proves nothing. Real measured blends are few (hundreds of systems), so we use physics for the bulk and ML only for the pair correction.

**What exactly is learned?** One symmetric pair function g(i,j,T) (LightGBM, 400 trees) added to the Refutas log-viscosity weighted by x_i x_j. Features are descriptors of the two molecules, not the identities.

**How do you know it is not over-fit?** Whole systems are held out; 70 locked systems were scored once; bootstrap CI over systems excludes zero (gain 0.0233, CI [0.0163, 0.0302]); the model is compared with a constant-offset comparator (0.043) too.

**What does 'certified' mean? Is it safe to fly?** No. It means: finite margin, in domain, and the spec is met after the margin, *under the model*. Nothing was measured; real certification needs ASTM tests.

**Why is freezing 'no skill' but you still use it?** Only through its lab-derived margin, for known molecules; on real fuels it is a safe-side screening estimate. We state that it does not beat a constant.

**Why did flash get so much better (-21.7 to +0.1 bias)?** The old API path used wrong routing (`api_as_shipped`); re-routing `predict()` onto the BASE route fixed it. The new numbers are descriptive on fuels used to develop the route.

**Why is your margin 2.2 C but the yaml says 3.7?** 3.699 was derived on 31 groups before mesitylene/new systems; 2.207 (35 systems, mesitylene excluded) is applied. Mesitylene is disclosed as an outlier.

**Why is the Jet A-1 design above the AJFTD dicyclo range?** The optimiser is free to leave the observed real-fuel range; we flag it (37.8 % vs 30 % max) and do not claim realism.

**Why is the ML density only an option?** It failed the ternary gate (3.64 to 3.92 kg/m3).

**Why Refutas and not something fancier?** It is a standard petroleum-industry rule, it needs only each component's viscosity, it has no fitted parameters of our own, and it scored well against real fuels above 263 K. The ML part is a small correction on top of it.

**What is a conformal margin, in one sentence?** A distribution-free allowance read off the sorted worst errors on held-out systems so that about 90 % of new systems stay inside it in the dangerous direction.

**Where are the weaknesses?** Thin locked sets (3 and 4 cold viscosity systems), single GC fuel for D56, binary-only margins applied to multicomponent blends, no lab validation, reconstructed CRC compositions, Jet A-1 limits from secondary sources.

**What would you do next?** Obtain lab flash/freezing/viscosity on the three certified blends; collect ternary+ lab data for margins; verify Def Stan 91-091 limits; score the Oakley/Yang sets; test ML density on real fuels.

**Can you explain a file?** Use section 11: start with `api/predict_route.py` (entry), `eval/composition_route.py` (what is computed), `baselines/flash.py` (145 lines, easiest to explain), `ml/viscosity_pair.py` (173 lines), `uq/conformal.py` (90 lines).

---

# 14. One-page equation sheet

```
x_i = (w_i/M_i)/sum(w_j/M_j)            v_i = (w_i/rho_i)/sum(w_j/rho_j)
1/rho_mix = sum w_i/rho_i               rho_mix(ML) = sum x_i M_i / [sum x_i M_i/rho_i + sum_{i<j} x_i x_j g]
VBN_i = 14.534 ln(ln(nu_i+0.8)) + 10.975;  VBN_mix = sum w_i VBN_i;  nu_mix = exp(exp((VBN_mix-10.975)/14.534)) - 0.8
ln eta_mix = ln eta_Refutas + sum_{i<j} x_i x_j g(d_i,d_j,T);   nu = 1000 eta / rho
log10 P_sat = A - B/(T+C);   H(T) = sum x_i P_sat,i(T)/(P_tot LFL_i) = 1  ->  T_flash
1/T_i = 1/T_m,i - R ln(x_i gamma_i)/dH_fus,i;  T_freeze = max_i T_i
SP = 28.84 - 0.324 A - 0.791 D;   dH_c = sum w_i dH_c,i
s_g = max_{p in g}(error, dangerous side);  k = ceil((n+1)(1-alpha));  margin = s_(k)  (infinity if k>n)
certified <=> margin finite AND in_domain AND (flash - margin >= 38 C, freezing + margin <= limit, ...)
```

---

# 15. Reading the results folders: how to find any number

1. Start at `docs/final_numbers_v5.md`. Sections A to H hold every headline number in tables; the CSV behind each is in `results/v5/final/` (table_A_headline.csv and so on; `results/v5/final/README.md` maps each table to its sources).
2. For one experiment, open `results/v5/BASE/delta_<name>/README.md`. It states the question, the pass rule (written before the result) and the reading. The numbers are in `readings.json`, `summary.csv` or `rows.csv` next to it.
3. The frozen scoreboard is `results/v5/BASE/scoreboard.csv` (human summary `SCOREBOARD.md`); one row per property x tier x set x model.
4. A claim such as CLM-73 is in `docs/claims_register.md` with its status and evidence path.
5. Margins as applied: `configs/uncertainty_v5.yaml`. Spec limits: `configs/specs_astm.yaml`, `specs_bis.yaml`, `specs_jet.yaml`.
6. Failed models: `docs/v5_negative_results.md`, `docs/model_decision.md`.
7. Anything changed since the baseline is a visible hash difference: run `scripts/v5/base/reproduce_base.py`.

**Column names in the scoreboard.** `n_groups` = number of independent systems/fuels; `bias` = mean(pred - meas); `mae` = mean absolute error; `mae_trivial` = the leave-one-group-out constant's MAE; `verdict` as in A7; `margin90` = the 90 % one-sided margin (infinity when too few groups); `label` = `descriptive`, `dev` or `already_read`.

---

# 16. What a prediction looks like end to end

1. You pass a list of (molecule, mass fraction) pairs, a temperature in kelvin and a fuel type. Molecules can be names, SMILES or InChIKeys.
2. The API resolves every molecule against the compound registry and the frozen pure-property table (`data/cache/pure_grid.npz`, whose SHA-256 is checked at load; a mismatch stops the run).
3. It converts mass fractions to the mole and volume fractions the rules need.
4. For each property it runs the designated route: ideal volume for density; Refutas (plus the ML pair correction if the blend is in the ML domain) for viscosity; Le Chatelier with anchored Antoine constants and the LFL overlay for flash point; UNIFAC Schroeder-van Laar for freezing; group additivity for heat; the aromatics regression for smoke.
5. It attaches, per property, the `value`, the `method` that produced it, the `status`, `in_domain`, a `lo/hi` interval from the applicable margin (or none), the alternative `variants`, and `warnings` (for example "stand-in used", "outside ML domain").
6. In design mode the optimiser calls this same function, then the judge checks each limit with the margin applied and labels the outcome `certified`, `meets_at_point_estimate`, `fails` or uncertified.

Reading the output: trust a number only after checking `in_domain`, `status` and `warnings`; a value without a margin is a point estimate, not a guarantee.
