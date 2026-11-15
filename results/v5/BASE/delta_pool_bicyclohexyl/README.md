# delta_pool_bicyclohexyl: flash-registry entries for bicyclohexyl and cis-decalin (opt-in overlay); pure pool P0y

Declared 2026-11-15 BEFORE any declared reading in this folder was computed. Base: v5 after commit e92d5c6 (delta_d86_rayleigh).
Origin: `docs/audit_2026-11-15.md` section 5.1 ("Getting a flash-registry entry (Antoine + LFL) for bicyclohexyl would add a 239 C cycloalkane").
Gates, margins, spec limits, spec `verification:` flags, the judge and `configs/uncertainty_v5.yaml` are not changed.

## Orientation (no declared reading)

* Pure-table compounds with verified freezing constants but no default-registry Antoine: **bicyclohexyl** (Tb 239 C) and **decalin_cis**
  (Tb 196 C) (also ethanol and 1-methylnaphthalene, not used: not a hydrocarbon jet component / counts as a naphthalene). Both have a registry row
  with an LFL (Jones estimate, the same tier as most registry LFLs) and a blank Antoine.
* The registry tier `anchored_thermo_ambient` (399 compounds) is an Antoine fitted to `thermo`'s default `VaporPressure` method over
  273.15-373.15 K (checked on 6 entries: ratio to thermo 0.99-1.01 across the window). For bicyclohexyl thermo's default is `LANDOLT`
  (Landolt-Boernstein Antoine, fitted range 424-577 K, so the 273-373 K window is an extrapolation); P(372.5 K) = 1.26 kPa against the
  WebBook reduced-pressure boiling point 0.01 bar at 372-373 K (Buckingham and Donaghy 1982): 25 % high, i.e. more volatile (safe side for flash).
* Measured pure flash points of bicyclohexyl in the unified table: 74.0 C (x5: Li 2020 set and Saldana 2011) and 91.85-92.2 C (x4). Two clusters.
* Bicyclohexyl has a 7.08 kJ/mol solid-solid transition 3.3 K below its melt; `tm_dhfus_verified.csv` governs SvL with fusion-only 6.86
  kJ/mol (`svl_dhfus_choice = fusion`), which is the less conservative choice. Not changed here; reported as a sensitivity (R-BC-3).

## Why an overlay

`data/processed/v2/compounds_anchored_vp.csv` also feeds `SubReps` (GCxGC representative choice for real fuels) through `_lookup`. Writing the
new constants into it could move representatives and so every frozen T2 / WS7 / P2 / P3 number. The new rows therefore live in
`data/processed/v5/flash_registry_additions.csv` and are merged into the default registry **only** when the environment variable
`FUELBLEND_FLASH_REGISTRY_ADDITIONS=1` is set before the first registry read. Off (default): bit-identical lookups.

## Declared readings

* **R-BC-1 (data).** For bicyclohexyl and cis-decalin: Antoine (kPa, K) by least squares in log10 P on thermo's default method at 21 equally
  spaced temperatures 273.15-373.15 K (the `anchored_thermo_ambient` procedure), `antoine_source = anchored_thermo_ambient_v5add`,
  LFL = the existing registry value. Reported, not gates: method used, max |log10 ratio| on the window, Tb from the Antoine vs literature
  (bicyclohexyl 239, cis-decalin 195.8 C), pure Le Chatelier flash vs measured flash in `data/unified/measurements.csv`.
  **Gate for using an entry in R-BC-2:** Antoine-extrapolated Tb within 10 C of the literature Tb (the D86 route inverts it), else the
  compound is dropped from the pool and recorded.
  Test: with the variable unset, `_lookup` returns None for both and the registry frame equals the base frame.
* **R-BC-2 (pool P0y).** P0y = P0x-b (pool C + tetralin + n-tridecane + hexylbenzene) + the R-BC-1 compounds that pass the gate.
  Protocol of delta_smoke_branch_opt unchanged: ASTM Jet A and BIS Jet A-1, D86 on (designated Engler route and WS7 margins), 10 starts,
  seed 42, lambda 0.25, `alt_standin` on, margins `configs/uncertainty_v5.yaml`. Reported: every constraint (status, prediction, slack,
  reason), binding constraints, certified or not. Certified only if the unchanged judge certifies.
  **Also reported for each optimum (checks, not gates):** the Rayleigh D86 quantities (delta_d86_rayleigh model) with the Rayleigh margins;
  if the judge certifies but Rayleigh T90-T10 < 40 + 5.31 C or T50-T10 < 15 + 1.24 C, the design is flagged `d86_model_dependent`.
* **R-BC-3 (freezing sensitivity).** For each optimum, freezing recomputed with bicyclohexyl's SvL dHfus = 13.94 kJ/mol (fusion + the
  7.08 kJ/mol transition). Reported; a design whose freezing slack turns negative is flagged `freezing_dH_choice_dependent`.

Stated expectation (before running): the optimiser uses bicyclohexyl to lift T90 and flash; whether T90-T10 >= 58.25 C (Engler + margin) can
be met together with flash >= 40.2 C, smoke and density is unknown. A smoke `uncertified` from the class-mix domain rule
(`smoke_margin_v5`) is possible for a recipe rich in dicycloparaffins.

## Results R-BC-1 and R-BC-2 (appended 2026-11-15 after the runs; nothing above was edited)

### R-BC-1 (`r_bc_1.csv`, `data/processed/v5/flash_registry_additions.csv`)
| compound | thermo method | A / B / C (kPa, K) | Tb Antoine vs literature (C) | pure flash pred vs measured (C) | gate |
|---|---|---|---|---|---|
| bicyclohexyl | LANDOLT (Clausius extrapolation below its 424 K range) | 7.3398 / 2697.29 / 0.000 | 232.5 vs 239.0 (-6.5) | 85.5 vs 74.0 (x5) / 91.9-92.2 (x4) | pass |
| cis-decalin | ANTOINE_POLING | 7.2854 / 2406.36 / -5.07 | 187.7 vs 195.8 (-8.1) | 57.7 vs 57.5 (median of 6) | pass |
Overlay tests: `tests/test_delta_pool_bicyclohexyl.py` (off by default: no entry; on: entry) pass.

### R-BC-2 (pool P0y): not certified; the optimum is invalid
Both specs converged 10/10 to about **37-39 wt % isooctane**, 24 % butylcyclohexane, 13 % bicyclohexyl, 13 % cis-decalin, 9.5 % TMB.
`flash_point_C` and `nhoc_MJ_kg` are `uncertified`, reason `model unavailable: decalin_cis.MW status=unverified` (the pure table holds
cis-decalin's MW at an unverified tier). Every other constraint passes. Because the flash model raised `ModelUnavailable` on every recipe
containing cis-decalin, the optimiser received no flash signal and drove isooctane up (it buys T10, density and smoke slack). The real flash
point of such a recipe is far below 38 C. The judge did not certify it: **0 of 2 certified**.

Finding (optimiser, not the judge): a constraint whose model is unavailable for a pool member is silently dropped from the optimiser's
objective. Pool members must therefore carry every tier the judge needs, or the optimum can be meaningless. Not fixed here (engine change);
the follow-up below removes the member instead.

## Declared follow-up R-BC-4 (declared after R-BC-2, before running it)
* Pool **P0y-b** = P0x-b + bicyclohexyl (cis-decalin dropped: its MW tier makes flash and NHOC unavailable). Protocol, checks (Rayleigh D86,
  `d86_model_dependent`) and R-BC-3 sensitivity exactly as R-BC-2. Certified only if the unchanged judge certifies. A constraint reported
  `uncertified` with reason `model unavailable` makes the run invalid as a design (reported as such).

## Results R-BC-4 and R-BC-3 (appended 2026-11-15 after the runs; nothing above was edited)

`designs/P0yb/`, `r_designs.csv`, `r_bc_checks.csv`, logs in `logs/P0yb_*.log` (10 starts, seed 42, about 25 min each).

### R-BC-4 (pool P0y-b = P0x-b + bicyclohexyl): the unchanged judge certifies both optima, with the overlay ON
| spec | judge | min slack | binding | recipe (wt %) |
|---|---|---|---|---|
| ASTM Jet A | **certified** (10 pass, 0 fail, 0 uncertified) | 0.0115 | flash, freezing, NHOC, smoke, T50-T10 | butylcyclohexane 27.7, TMB 17.8, bicyclohexyl 17.3, n-C10 15.1, n-C11 13.6, isooctane 3.0, n-C12 2.9, n-C13 2.6 |
| BIS Jet A-1 | **certified** (10 pass) | 0.0058 | flash, freezing, NHOC, smoke, T50-T10 | bicyclohexyl 35.9, butylcyclohexane 24.7, TMB 14.4, n-C10 9.3, n-C11 7.9, isooctane 3.5, trans-decalin 1.6, n-C12 1.4, n-C13 1.3 |

Predictions (Jet A / Jet A-1): flash 40.64 / 40.43 C (limit 38 + margin 2.21), freezing -57.4 / -64.2 C (+16.98 margin vs -40 / -47),
density 801.8 / 823.7, -20 C viscosity 3.34 / 4.23 cSt, smoke 22.3 / 23.2 mm (alt limit with naphthalenes 0), Engler T90-T10 70.4 / 69.2 C.
No constraint is `uncertified ... model unavailable`, so the runs are valid as designs under the R-BC-4 rule.

### Declared checks (`r_bc_checks.csv`)
| design | Rayleigh T50-T10 / T90-T10 (C) | needs (limit + Rayleigh margin) | freezing with bicyclohexyl dH 13.94 kJ/mol | flags |
|---|---|---|---|---|
| P0y-b Jet A | 9.5 / 35.4 | >= 16.24 / >= 45.31 | -57.4 (bicyclohexyl not first to freeze; unchanged) | **d86_model_dependent** |
| P0y-b Jet A-1 | 14.6 / 44.0 | >= 16.24 / >= 45.31 | -40.7 (+16.98 = -23.7 > -47: fails) | **d86_model_dependent, freezing_dH_choice_dependent** |

Reading. These are the first judge-certified designs of v5, and both rest on choices that are not designated:
1. the flash-registry overlay (bicyclohexyl Antoine extrapolated from LANDOLT data at 424-577 K), opt-in and not adopted;
2. the Engler D86 route: under the zero-parameter Rayleigh model both recipes miss the D7566 slope requirements (the R-DR-3 pattern:
   pure-compound pools have flat distillation curves);
3. for Jet A-1, the fusion-only dHfus of bicyclohexyl.
They are therefore reported as **certified by the judge, conditional**, not as certified designs. 0 unconditional certifications.

## Declared check R-BC-5 (declared after R-BC-4, before running it)
* **Flash data sensitivity.** The overlay's pure bicyclohexyl flash is 85.5 C. The measured values form two clusters, 74.0 C (x5) and
  91.85-92.2 C (x4), so the overlay is on the unsafe side of the lower cluster. For each P0y-b optimum the mixture flash is recomputed with the
  same flash model and bicyclohexyl's Antoine A shifted (B, C unchanged) so that its pure Le Chatelier flash equals 74.0 C. If prediction
  minus the 2.21 C margin falls below 38 C, the design is flagged `flash_data_dependent`. Reported; it does not change the judge's result.

## Results R-BC-5 (appended 2026-11-15 after the run; nothing above was edited)

`python scripts/v5/base/delta_pool_bicyclohexyl.py flash5 --variant overlay|shift74` -> `r_bc_5.csv` (shifted overlay in `r_bc_5_additions_shift74.csv`, never read by default).

| variant | bicyclohexyl A | Tb from Antoine (C) | pure flash, Le Chatelier / model (C) | mixture flash Jet A / Jet A-1 (C) | bound = pred - 2.21 (limit 38) | flag |
|---|---|---|---|---|---|---|
| overlay (as certified) | 7.3398 | 232.5 | 85.5 / 84.6 | 40.65 / 40.43 (reproduces the judge) | 38.44 / 38.22 | - |
| shift to 74.0 C | 7.5889 | 210.0 | 74.0 / 73.3 | 40.45 / 40.02 | 38.25 / **37.82** | Jet A-1 **flash_data_dependent** |

Reading. As expected the mixture flash is set mostly by the light members (n-decane, TMB, butylcyclohexane, isooctane), so the shift moves
it by 0.2 C (Jet A, 17 % bicyclohexyl) and 0.4 C (Jet A-1, 36 %). Jet A keeps a 0.25 C flash slack; Jet A-1 loses its slack.
Data-consistency finding: with the registry (Jones) LFL, a pure flash of 74 C needs an Antoine whose Tb is 210 C, against the literature
239 C (the overlay's 232.5 C already reads 6.5 C low). The 74 C measured cluster, the literature Tb and the Jones LFL cannot all hold; the
92 C cluster is the one consistent with Tb. Not a gate; recorded.

## Summary of the conditions (for the user decision)
| design | judge | overlay needed | Rayleigh D86 slopes | bicyclohexyl dH with transition | flash at the 74 C cluster | UNIFAC-SvL + 13.84 (R-FU-3, counterfactual) |
|---|---|---|---|---|---|---|
| P0y-b Jet A | certified | yes | miss (9.5 / 35.4 vs 16.2 / 45.3) | passes (bicyclohexyl not first to freeze) | passes (38.25) | -43.1 vs -40, passes |
| P0y-b Jet A-1 | certified | yes | miss (14.6 / 44.0) | **fails** (-23.7 vs -47) | **fails** (37.82) | -49.7 vs -47, passes |

The Jet A design is the stronger one: its only open conditions are the overlay and the D86 model. Thermodynamically the fusion-only dHfus
is the less defensible choice for a liquidus far below bicyclohexyl's solid-solid transition (about 3 K below its melt), and the Jet A-1
liquidus (-64 C) is about 64 K below it (Tm 276.8 K, transition about 273.5 K), so the Jet A-1 design should be read as failing freezing unless the user decides otherwise.
