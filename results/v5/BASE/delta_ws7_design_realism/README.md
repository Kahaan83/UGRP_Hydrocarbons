# delta_ws7_design_realism: real blendstocks, D86 constraints, a blendstock optimiser, and the lab-test request (WS7)

Written 2026-11-14 BEFORE any declared reading was computed. Plan: `docs/superpowers/plans/2026-11-13-v5-certification-and-ml-roadmap.md`, section WS7.
Base commit 50f4333 (v5). The smoke-pool result (6dc66f8 on v5; pool C, 0/14 certified, counterfactual True for both specs) is read with
`git show` and is the comparison for the optimiser (R-W7-5). Spec `verification:` flags are never edited. Nothing here is certified unless the
judge says so, and the judge never certifies an unverified limit.

## Pre-declaration orientation (no reading computed; what shaped the declarations)

* Blendstock compositions present in the repo, mapped with the T2 route (`composition_route.kind_bins` -> `map_bins` with `reconstruction_r3.SubReps`,
  whose candidates all have a pure grid and a default-registry Antoine/LFL):
  * Jet A base: AJFTD GCxGC fuel_id 2 (label 10325, POSF-10325 Jet A), 42 representatives, unmapped 2.47 % (diaromatics without rep).
  * HEFA-SPK: FuelLib `HEFA-CAME`, `HEFA-MFAT`, `HEFA-TALL` (category SATF), each maps 100 %.
  * ATJ-SPK: FuelLib `C-1 POSF11498 ATJ`, unmapped 0.44 %.
  * SAK: **no named SAK composition exists in the repo** -> SAK is **blocked** (not invented). AJFTD GCxGC fuels 47 and 139 are ~97 % aromatic and
    several AJFTD fuels are ~98 % iso-paraffin, but they carry no pathway label; they are recorded as "unlabelled aromatic / iso-paraffinic stream"
    and are not used as SAK / ATJ.
* The existing distillation baseline (`baselines/distillation.py`, TBP on volume fractions + fixed Engler offsets +17.24 / +2.61 / -15.06 C) was
  never scored: every T2 real-fuel row was `tb_coverage_lt_0.90` because Tb came from the `compounds.csv` Antoine. The default flash registry has an
  Antoine for every SubReps representative. Tb from it: n-C10 172.2, n-C12 210.9, butylcyclohexane 179.0, isooctane 99.2, 1,2,4-TMB 169.1,
  trans-decalin 180.4 C (literature ~174.1 / 216.3 / 180.9 / 99.2 / 169.4 / 187.3: 0-6 C low; Antoine extrapolated to Tb).
* Pushing pool C's stored Jet A recipe through the baseline gives D86 T10 175.8, T50 177.5, **T90 170.8 C** (non-monotone, T90 < T10): the fixed
  offsets shrink T90-T10 by 32 C, so a narrow-boiling pure-compound recipe gets an impossible curve. A domain rule is needed (below).
* Provenance of the Engler offsets: the docstring says "ASTM D86 vs D2887 on CRC aviation fuels" (no file or fit script found). They are NOT fitted
  here; whether they were fitted on AV-33-22 cannot be verified, so the margin fit below may be partly in-sample for the offsets (stated, not fixed).
* `cma` is not installed: the optimiser is multi-start SLSQP.

## Declared data rules

* D86 margin data: measured-GCxGC real fuels only. AV-33-22 (`crc_av3322_properties.csv`, method `D86` rows only, never `distillation_simdis_*`)
  and AJFTD GCxGC fuels (`ajftd_fuels_long.csv` distillation rows whose `method` contains `D86`; rows without a method are excluded and counted).
  **CRC 647 is excluded** (its compositions were reconstructed from D86, ruling R10: circular). FuelLib fuels are not used for the margin (4 fuels,
  of which HEFA/ATJ are the design blendstocks themselves).
* One fuel per double-count group: `results/v5/BASE/t2_fuels.csv` `chosen_for_group`; a fuel absent from that table is its own group.
* Mapped coverage >= 0.90 of mass (`COVERAGE_MIN`), else excluded and counted. Measured values: mean of repeat tests per fuel and quantity.
* FBP = `distillation_T100_C` (D86 end point).

## Declared models (fixed now)

* **Tb**: inversion of the default flash-registry Antoine at 101.325 kPa (`log10 P/kPa = A - B/(T+C)`); representatives and pure components alike.
* **D86 T10 / T50 / T90**: `predict_d86_from_composition` unchanged (volume fractions from the pure grid rho at 288.15 K, as T2 `d86_preds`).
* **FBP model (new, declared before any residual)**: TBP cut interpolated at 99 vol % on the same cumulative curve, no offset.
* **T50-T10, T90-T10**: differences of the D86 cuts above.
* **Distillation domain** (hard in the judge, graded distance in the optimiser): the D86 curve must be monotone (T10 <= T50 <= T90 <= FBP) and the
  recipe's TBP T90-T10 must lie inside [min, max] of the TBP T90-T10 of the margin-fit fuels. Outside: `uncertified`, reason `outside_distillation_domain`.

## Declared readings and gates (never moved after results)

* **R-W7-1 (descriptive)**: blendstock table: name, source, provenance, n representatives, unmapped fraction, class mass fractions, pathway label
  or "blocked"/"unlabelled".
* **R-W7-2 (D86 margins)**: per quantity q in {T10, T50, T90, FBP, T50-T10, T90-T10}: residual r = meas - pred for the maximum limits (T10, FBP:
  under-prediction unsafe, direction `under`) and r = pred - meas for the minimum limits (T50-T10, T90-T10: over-prediction unsafe, direction
  `over`); T50, T90 are reported with both sides (no limit). Fit on each constraint's OWN residual (a difference is never built from two point margins).
  Split conformal over fuels, alpha = 0.10: margin = the ceil((n+1)(1-alpha))-th smallest r, floored at 0. Per source too; **certifying margin =
  max(pooled margin, per-source margin of every source with >= 10 fuels)**. **Gate: >= 30 independent fuels with both AV-33-22 and AJFTD contributing
  >= 10 each, else the margin is `inf`** (constraint stays `uncertified`). Also reported: bias, MAE, MAE of the constant (mean) predictor.
  Margins are written to `d86_margins.json` in this folder; adopting them into `configs/uncertainty_v5.yaml` is a **user decision** (not done).
* **R-W7-3 (Tb sanity)**: the Tb values above vs literature, reported as a table (descriptive).
* **R-W7-4 (feature off unchanged, unit test)**: with D86 off (`load_spec(..., d86=False)`, the default) every existing spec loads exactly the
  same constraints, and the pool C stored recipes give identical in-memory statuses, predictions, slacks and margins through `evaluate_recipe`
  before and after the change. With a single blendstock at weight 1 the blendstock evaluation equals the rep-level `evaluate_recipe` exactly.
* **R-W7-5 (optimiser comparison, seed 42)**: the new blendstock optimiser (multi-start SLSQP in top-level blendstock fractions, graded objective,
  domain hard constraints) run on pool C itself (every member a pure compound, D86 off, uncertainty penalty off), compared with the 6dc66f8 engine
  result: min normalised slack (Jet A 0.0187, Jet A-1 0.0155) and statuses. Reading: equal to 1e-3 = "reproduces"; higher = better optimum found;
  lower = worse (reported plainly). Also: the same with the uncertainty penalty on (lambda declared below), reported as the change in min slack and recipe.
* **R-W7-6 (designs)**: specs `astm:jet_a` and `bis:jet_a1`, `synthetic_blend`, `alt_limits`, `recipe_margins`, D86 on, seed 42. Pools:
  * **P0 control = pool C** (7 verified pure compounds) with D86 on;
  * **P1 = Jet A base + pool C pure compounds** (8 top-level members);
  * **P2 SAF = Jet A base + HEFA-SPK (FuelLib HEFA-CAME, first alphabetically; all three HEFA map 100 %) + ATJ-SPK (POSF11498) + 1,2,4-TMB**
    (verified aromatic), with the D7566 blending rule as a hard constraint: HEFA + ATJ <= 50 % (mass, used as the proxy for vol %; stated).
  Per design and constraint: status, slack, margin, direction, prediction, reason; certified or not; binding constraints (slack within 0.003 of min);
  `would_certify_if_specs_verified` (COUNTERFACTUAL, rule of delta_smoke_pool ported unchanged: only the `secondary_sources_only` / `UNVERIFIED` clauses
  and the alt-met encoding are removed; every other reason blocks).
* **Represented-blendstock rule (declared)**: a constraint of a recipe with > 0 mass from a GCxGC-represented blendstock is judged by margins
  validated on known compositions (lab mixtures) and real-fuel margins were not adopted (WS6, user decision). Every model-predicted (not
  by-construction) constraint of such a recipe is therefore at most `uncertified`, reason `represented_blendstock_margin_not_validated`, and that
  reason blocks the counterfactual.
* **Uncertainty penalty (declared)**: objective = min margin-adjusted normalised slack - lambda * U, lambda = 0.25, U = max over constraints whose
  property has alternate model forms (`models["alternates"]`) of |alternate - primary| / spec scale (largest alternate). Optimiser only; the judge never sees it.
* **Domain hard constraints in the optimiser**: (i) extrapolation (absent-class wt fraction <= 0.6 per model), (ii) the synthetic-blendstock cap
  in P2, (iii) the distillation domain (D86 on). The remaining smoke/class-mix and ML-viscosity domains enter as before through the graded objective.
* **Expected outcomes, stated now**: 0 designs certified (aromatics vol % and smoke limits unverified; the new D86 limits are also
  `secondary_sources_only`). P0 (pure compounds) fails or is uncertified on the D7566 slope limits (narrow boiling). P1/P2 are uncertified for constant
  provenance (representatives are not in the verified table: route constants for density/freezing/flash) and for the represented-blendstock rule, so their
  counterfactual is False with those blockers listed. That is a valid negative: realism costs certifiability until real-fuel margins exist.
* **R-W7-7 (final check)**: `LAB_TEST_REQUEST.md` lists 1-2 candidate designs and the tests (D56 flash, D445 viscosity at -20 C, D4052 density,
  D5972 freezing, D1322 smoke; plus D86 because of the new constraints) as **needs user decision**. Nobody is contacted.

## Results

### Implementation deviations (mechanical, stated before the design readings were read)

* `build_route_models` normally refuses a pool member without Tm/dHfus. Many GCxGC representatives have none (iso-paraffins: 8 of 8 ATJ reps,
  19-20 compounds in P1/P2), so `blendstock_v5.build_models` requires grid + flash data only; a member without Tm/dHfus makes freezing
  `model unavailable` (uncertified, and the constraint drops out of the optimiser) instead of refusing the pool. Listed in `models['sle_gaps']`.
* A Jet A representative shares its InChIKey with pool C's isooctane under the table name `2,2,4-trimethylpentane` (no verified MW/rho15); a pure
  member now keeps its own table name and shared InChIKeys resolve to the most-verified table name (commit bd4647a).
* Design runs use 10 starts (P0/P1/P2, about 43-50 compounds per pool, 0.06-0.14 s per objective evaluation plus 0.07-0.2 s for the
  uncertainty term); the optimiser comparison on pool C uses 30 starts, like the 6dc66f8 engine runs.

### R-W7-1 blendstocks (`blendstocks.csv`)

| name | source | reps | unmapped | class wt (aro / cyclo / iso / n) | status |
|---|---|---|---|---|---|
| jetA_posf10325 | AJFTD GCxGC fuel 2 (POSF-10325 Jet A) | 42 | 0.0247 | 0.166 / 0.324 / 0.305 / 0.205 | used (Jet A base) |
| hefa_came | FuelLib HEFA-CAME | 25 | 0.0 | 0.000 / 0.028 / 0.901 / 0.071 | used (HEFA-SPK) |
| atj_posf11498 | FuelLib C-1 POSF11498 ATJ | 8 | 0.0044 | 0 / 0 / 1.000 / 0 | used (ATJ-SPK) |
| SAK | - | - | - | - | **blocked**: no named SAK composition in the repo |
| ajftd 47, 139 | AJFTD GCxGC | - | - | ~97 % aromatic | not used: no pathway label |

### R-W7-2 D86 margins (`d86_fit_rows.csv`, `d86_margins.json`)

65 fuels used (AV-33-22 35, AJFTD GCxGC 30); 9 dropped as double-count duplicates; 0 below coverage; every AJFTD distillation row had a D86 method.
**Gate passed for every quantity.** Distillation domain (TBP T90-T10 of the fit fuels): 30.1-163.6 C.

| quantity | direction | pooled | AV-33-22 | AJFTD | **certifying margin (C)** | bias pred-meas | MAE | MAE constant |
|---|---|---|---|---|---|---|---|---|
| T10 (max 205) | under | 6.37 | 4.00 | 7.93 | **7.93** | -0.69 | 4.52 | 8.19 |
| FBP (max 300) | under | 0.82 | 0.82 | 3.26 | **3.26** | +15.11 | 15.75 | 16.52 |
| T50-T10 (min 15) | over | 11.02 | 3.60 | 18.68 | **18.68** | -1.58 | 7.67 | 9.32 |
| T90-T10 (min 40) | over | 8.79 | 0.00 | 18.25 | **18.25** | -8.98 | 12.56 | 16.53 |
| T50 (no limit) | under | 12.46 | 14.96 | 8.95 | 14.96 | -2.27 | 8.75 | 11.45 |
| T90 (no limit) | under | 20.54 | 17.56 | 26.66 | 26.66 | -9.68 | 13.06 | 15.31 |

Reading: the baseline beats the constant predictor on every quantity (T10 MAE 4.5 vs 8.2 C) but the slope quantities need margins of about
18 C because of the AJFTD fuels, so a design needs predicted T50-T10 >= 33.7 C and T90-T10 >= 58.2 C to clear the D7566 minimums with the
margin. The FBP model (TBP at 99 vol %) over-predicts by 15 C on average, so its under-prediction margin is small (3.3 C); it is conservative,
not accurate. The Engler offsets' provenance is unverified (possible in-sample use of AV-33-22). The margins live in `d86_margins.json` only;
adopting them into `configs/uncertainty_v5.yaml` is a user decision.

### R-W7-3 Tb sanity (`tb_sanity.csv`)

Registry-Antoine Tb minus literature: n-C10 -1.9, n-C11 -3.4, n-C12 -5.4, butylcyclohexane -1.9, isooctane 0.0, 1,2,4-TMB -0.3, trans-decalin
-6.9 C. Tb is biased low by up to 7 C for C11+ and bicyclics (Antoine extrapolated beyond its fitted range); the D86 margins absorb this
only to the extent the real fuels share it.

### R-W7-4 feature off unchanged (unit tests, `tests/test_v5_ws7_design_realism.py`, 13 passed)

With `d86=False` (default) all 24 spec loads equal the pre-change fixture (`tests/data/ws7_feature_off_fixture.json`, written before the
engine/yaml change) and the pool C stored recipes give identical statuses, reasons, predictions, slacks and margins. One blendstock at weight 1
(ATJ) equals the rep-level `evaluate_recipe` exactly. `tests/test_v4_design.py`, `test_v5_aromatics_min.py`, `test_v5_smoke_alt_limit.py`: 43 passed.

### R-W7-5 optimiser comparison on pool C (seed 42, 30 starts, D86 off) (`compare/`, `r_designs.csv`)

| spec | engine 6dc66f8 min slack | new optimiser, lambda 0 | new optimiser, lambda 0.25 | reading |
|---|---|---|---|---|
| Jet A | 0.018733 | **0.018733** (same recipe to 4 decimals, 17/30 converged) | 0.018222 (BCH 0.579, TMB 0.093; U lower) | **reproduces**; penalty costs 0.0005 slack |
| Jet A-1 | 0.015545 | **0.015545** (same recipe) | 0.014701 (TMB 0.100, aromatics 9.05 vol %) | **reproduces**; penalty costs 0.0008 slack |

The uncertainty penalty moves the recipe slightly away from the alternate-model disagreement (freezing slack 0.019 -> 0.042) at the cost of a
small loss in the binding slacks. All four keep `would_certify_if_specs_verified = True` (blocked only by the aromatics vol % and smoke spec flags).

### R-W7-6 designs with D86 on (seed 42, 10 starts, lambda 0.25) (`designs/P*/`, `r_designs.csv`)

| design | top-level recipe (wt) | n pass / fail / uncertified | failing constraints (point estimate or with margin) | counterfactual |
|---|---|---|---|---|
| P0 Jet A (pool C) | n-C10 0.18, n-C11 0.21, n-C12 0.06, BCH 0.05, decalin 0.18, TMB 0.26, isooctane 0.06 | 3 / 5 / 4 | flash 32.3 C, freezing -49.4 (margin), aromatics wt 25.5, smoke 20.1, **T90-T10 49.9 < 40+18.2** | False |
| P0 Jet A-1 | similar | 3 / 4 / 4 | flash, freezing, smoke, T90-T10 49.5 | False |
| P1 Jet A (Jet A base + pure) | Jet A base 0.85, n-C12 0.09, isooctane 0.03, n-C11 0.02 | 0 / 4 / 7 | flash 35.8 < 38 (nominal), freezing -39.8 > -40 (nominal), T50-T10 32.8 and T90-T10 55.9 (with margin) | False |
| P1 Jet A-1 | Jet A base 0.91, n-C12 0.06, isooctane 0.03 | 0 / 3 / 8 | flash 35.3, freezing -43.6 > -47, T50-T10 32.6 | False |
| P2 Jet A (SAF) | Jet A base 0.54, HEFA 0.46, ATJ 0, TMB 0.001 | 0 / 2 / 9 | density 774.9 < 775 and aromatics 7.99 < 8 vol % (both nominal, normalised slack -0.0008) | False |
| P2 Jet A-1 (SAF) | Jet A base 0.49, HEFA 0.10, ATJ 0.40, TMB 0.01 (synthetic cap 0.50 active) | 0 / 2 / 9 | freezing -46.96 > -47 and aromatics 7.99 < 8 (nominal, slack -0.0008) | False |

Readings:
* **0 of 10 designs certified** (and 0 of the 6 R-W7-6 designs). As declared, the 4 pool C comparisons would certify only if the two spec flags
  were verified; every design with D86 on or with a real blendstock has a non-spec blocker, so their counterfactual is False.
* **D86 is the binding problem for pure-compound pools (P0)**: no pool C recipe reaches T90-T10 >= 58.2 C (40 + the 18.2 C margin); widening
  the curve with isooctane drives flash below 38 C. This was the expected outcome.
* **A real Jet A base widens the curve (P1/P2 TBP spread inside the domain; P2 T90-T10 86-116 C passes with margin, T10 and FBP pass with
  margin)** but brings its own limits: POSF-10325 as represented here predicts flash ~35-36 C and freezing near -40 C, so P1 fails flash and
  freezing at the point estimate. The SAF pool P2 meets flash, freezing (Jet A), viscosity, NHOC and all four D86 limits at the point estimate
  and ends at the very edge of density (775) and aromatics (8 vol %): the max-min optimum sits at normalised slack -0.0008, i.e. it cannot meet
  both the D7566 aromatics minimum and the density minimum with the declared penalty on. Not a pass; reported as fail.
* Every P1/P2 model-predicted constraint is uncertified for two reasons that are not spec flags: (a) non-verified constants of the GCxGC
  representatives (route_grid rho15, route_anchored flash, Coutinho Tm/dHfus; 19-20 representatives lack Tm/dHfus altogether), so density,
  flash, freezing and aromatics vol % have no margin and the optimiser steered on point estimates; (b) the declared represented-blendstock rule.
  Realism therefore costs certifiability until real-fuel margins (WS6) and verified representative constants exist.
* n_converged: P0 2/10, P2 10/10 (other values in `r_designs.csv`); a flag, not proof of a global optimum.

### R-W7-7 final check

`LAB_TEST_REQUEST.md` (needs user decision; nobody contacted).

## SAK pool P3 (declared 2026-11-14 BEFORE any run)

Now that an open SAK composition exists (delta_sak_blendstock, Feldhausen 2022), P3 = P2 + `sak_virent_feldhausen2022` (Jet A base POSF-10325, HEFA-CAME, ATJ POSF11498, SAK, 1,2,4-trimethylbenzene; synthetic cap 0.50 as in P2). Motivation: P2 missed density and aromatics vol % by 0.0008 normalised slack; SAK is an aromatic (about 62 % C9 alkylbenzene) synthetic stream. Same protocol as P0 to P2: D86 on, 10 starts, seed 42, uncertainty penalty lambda as in P2, margins `configs/uncertainty_v5.yaml`, run on the merged v5 code with the user-verified smoke and aromatics flags and the 7 upgraded freezing constants (toluene, n-hexane, n-octane, tetralin, o-xylene, mesitylene, 2-methylpentane; none is in P3 except as pure members, so they should not matter here).
* **Reading, not a unit test.** Per spec: every constraint status, slack, margin and reason; binding constraints; certified or not; counterfactual column unchanged.
* **Stated expectation:** uncertified, because the four D86 limits are still `secondary_sources_only` and the GCxGC fuel constants for the Jet A base, HEFA, ATJ and SAK are not verified (non-flag blockers 2(a)/(b) above). A design cannot certify under either. Pass or fail of the point-estimate constraints is the result.
* No gate is moved; SAK composition is one 2022 development batch, not the D7566 annex product (see delta_sak_blendstock/README.md).

### SAK pool P3 result (2026-11-14, 10 starts, seed 42, lam 0.25, D86 on, synthetic cap)
Pool = Jet A POSF-10325, HEFA-CAME, ATJ POSF11498, SAK (Feldhausen 2022), 1,2,4-trimethylbenzene (43 compounds). Logs `logs/P3_run.log`, designs in `designs/P3/`.
* ASTM Jet A: SAK weight 0.0011 (Jet A 0.54, HEFA 0.46). Density fails by 0.0008 slack (774.9) and aromatics vol% by the same 0.0008 (7.986 vs 8); flash, freezing, viscosity, NHOC, D86 uncertified.
* BIS Jet A-1: SAK weight 0.0 (Jet A 0.49, ATJ 0.40, HEFA 0.10). Freezing fails by 0.0008 and aromatics vol% by 0.0008.
* Reading: SAK is not selected. It adds nothing the HEFA/ATJ blendstocks do not already give, so the pool-P2 conclusion (0 of 6 certified) is unchanged. Matches the declared expectation (uncertified). Reasons for `uncertified` remain: D86 limits `secondary_sources_only`, GCxGC-represented blendstock margins not validated, route-grid constants for several pool members. No gate moved.
