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

(appended after the runs)
