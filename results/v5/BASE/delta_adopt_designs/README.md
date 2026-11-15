# delta_adopt_designs: design re-runs under the adopted defaults (decisions 2026-11-15)

Declared 2026-11-15 BEFORE any design run in this folder (no smoke test, no timing, no full run has been made for it). Base: v5 after
commit `1b196e3`. Plan: `docs/superpowers/plans/2026-11-15-apply-decisions-and-audit-fixes.md`, Task 9. Spec of the decisions:
`docs/decisions_2026-11-15.md`.
Gates, margins, spec limits, spec `verification:` flags, the judge and `configs/uncertainty_v5.yaml` are not changed. Runner:
`scripts/v5/base/delta_adopt_designs.py`.

## What changed since the committed design readings (the adopted defaults, all live by default)

| decision | default now | commit |
|---|---|---|
| 1 | designated freezing route = UNIFAC Schroeder-van Laar (`svl_unifac`), margin `freezing_point_C.lab_liquidus.svl_unifac` = 13.8427 C; ideal SvL (`ideal_svl`, margin 16.9755 C) is a reported variant | `bcab8b4`, `9fbff24` |
| 2 | D86: Engler stays designated; Rayleigh is a reported check (`d86_model_dependent`) | no change |
| 3 | flash-registry additions overlay ON (`FUELBLEND_FLASH_REGISTRY_ADDITIONS=0` restores the v2 registry) | `5986a6b` |
| 4 | bicyclohexyl SvL dHfus = fusion + transition, 13.94 kJ/mol | `3ec5aa1` |
| 5 | exact aromatic matching in the Coutinho table (`build_tmap`) | `0b248e6` |
| 6 | `strict_constraints=True` in `design_pool` | `e93dc22` |
| 7 | cis-decalin MW verified (P0y now runnable under the strict guard) | `3ec5aa1` |

## Protocol (declared; identical for every row)

* Pools x specs: **P0, P0x-a, P0x-b, P0y, P0y-b, P1, P2, P3** x **`astm:jet_a`, `bis:jet_a1`** = **16 rows**.
  Folder / design-id names drop the hyphen: P0xa, P0xb, P0y, P0yb.
  * P0 = pool C (n-decane, n-undecane, n-dodecane, butylcyclohexane, trans-decalin, 1,2,4-trimethylbenzene, isooctane), as `delta_ws7_design_realism.make_pool`.
  * P0x-a = pool C + tetralin + n-tridecane; P0x-b = P0x-a + hexylbenzene (`delta_smoke_branch_opt.POOL_X`).
  * P0y = P0x-b + bicyclohexyl + cis-decalin (both passed the R-BC-1 gate in `delta_pool_bicyclohexyl/r_bc_1.csv`); P0y-b = P0x-b + bicyclohexyl
    (`delta_pool_bicyclohexyl.pool_members`).
  * P1, P2, P3 as `delta_ws7_design_realism.make_pool` (P1 = POSF-10325 + pool C; P2 = POSF-10325 + HEFA-CAME + ATJ POSF-11498 + TMB; P3 = P2 + SAK).
* Design call, every row: `design_pool(spec, bp, models, W.MARGINS, n_starts=10, seed=42, d86=True, lam=B.LAMBDA_UNC, alt_standin=True,
  strict_constraints=True, synthetic_cap=...)` with `models = B.add_d86_models(B.build_models(bp))` (Engler D86) and
  `W.MARGINS = configs/uncertainty_v5.yaml`.
* **`synthetic_cap`**: `B.SYNTHETIC_CAP` (0.50, the D7566 synthetic-blendstock cap) for P2 and P3, None otherwise, exactly as the committed
  P2 / P3 runs (`delta_ws7_design_realism`, `delta_smoke_branch_opt`). It is part of the protocol: dropping it would remove a spec constraint.
* **n_starts = 10, seed = 42. They will not be reduced.** The runner refuses to write into this folder with any other value; smaller
  values are only for smoke tests written elsewhere.
* The runner does not set `FUELBLEND_FLASH_REGISTRY_ADDITIONS` (the overlay is ON by default) and refuses to run if it is `0`.
  It records the git HEAD, `DESIGNATED_FREEZING_MODEL` and the overlay state in every row and log.
* **A pool for which the strict guard (or the model build) raises `ModelUnavailable` is a reading**: the row is recorded with status
  `model_unavailable` and the message, and the run continues with the next pool. Expected (from the plan, Task 9): 16 rows and no `ModelUnavailable`.

## Declared readings, per pool x spec

1. **Judge status**: certified (`passed`) or not; counts of pass / fail / uncertified; failing and uncertified constraint keys; every
   constraint (status, prediction, applied margin, slack, reason) in the per-design `__constraints.csv`.
2. **Min slack** and the binding constraints (slack <= min + 0.003, the rule of the committed tables).
3. **Recipe**: top-level weights and compound weights.
4. **Flash**: prediction, applied margin, bound = prediction - margin, limit.
5. **Freezing**: the judge's UNIFAC-SvL value (designated), its applied margin and **safe bound = value + margin** (13.8427 C) vs the limit
   (-40 C Jet A / -47 C Jet A-1); and the **`ideal_svl` variant**: `composition_route.freezing_C(..., model="ideal_svl")` on the same
   designated freezing constants and RDKit-MW mole fractions, with its bound = value + 16.9755 C. The variant is reported, never judged.
6. **D86, Engler** (designated): the judge's T10, T50, T90, FBP, T50-T10, T90-T10 predictions and slacks.
   **D86, Rayleigh** (check): `B.rayleigh_quantities_smiles` with rho at 15 C, margins from `delta_d86_rayleigh/d86_rayleigh_margins.json`
   (T50-T10 1.2427 C, T90-T10 5.3057 C). If a member lacks an Antoine or rho the Rayleigh check is recorded `unavailable` (not a failure).
7. **Dicycloparaffin mass fraction** (reported only; **no gate**): sum of compound mass fractions whose molecule is non-aromatic, has only
   single bonds and exactly 2 rings (RDKit SSSR). Bicyclohexyl and the decalins count; tetralin (aromatic) does not. Context only:
   the 113 mapped AJFTD GCxGC fuels have max 30.0 % and 95th percentile 17.3 % (`delta_pool_bicyclohexyl`, undeclared observation).

## Conditions table (filled after the runs, one row per design, as in `delta_pool_bicyclohexyl`)

| design | judge | overlay needed | Engler D86 | Rayleigh D86 slopes | freezing UNIFAC bound vs limit | ideal_svl bound vs limit | dicyclo wt % | label |
|---|---|---|---|---|---|---|---|---|

Labels (declared):
* **`d86_model_dependent`**: the judge certifies the design (Engler, designated by decision 2) but Rayleigh T90-T10 < 40 + 5.3057 C or
  T50-T10 < 15 + 1.2427 C. Such a design is reported as *certified by the judge, conditional on the Engler D86 route*, not as an
  unconditional certification.
* **`freezing_route_dependent`**: the judge certifies (UNIFAC) but the `ideal_svl` bound exceeds the limit. Reported, not a gate.
* **`overlay_dependent`**: a certified design that contains a member whose flash-registry entry comes only from the overlay
  (bicyclohexyl, cis-decalin).
* R-BC-5 (the 74 C bicyclohexyl flash cluster) is **not re-run** here; its committed reading (`delta_pool_bicyclohexyl/r_bc_5.csv`)
  stands for the committed recipes only and is cited, not transferred, to any new recipe.

## Comparison to the committed readings (filled after the runs)

Each design is compared to its committed reading, field by field: passed, min slack, binding set, per-constraint status / prediction /
slack (delta = new - committed), recipe (per-member weight delta and the L1 norm of the top-level weight difference).

| pool | committed reading | design ids | protocol differences beyond the adopted defaults |
|---|---|---|---|
| P0 | `delta_ws7_design_realism/designs/P0/` | `<spec>__P0` | committed run had `alt_standin` off; the delta mixes the adopted defaults and the stand-in |
| P1 | `delta_ws7_design_realism/designs/P1/` | `<spec>__P1` | same: committed run had `alt_standin` off |
| P0x-a, P0x-b | `delta_smoke_branch_opt/designs/P0xa`, `P0xb` | `<spec>__P0xa_standin`, `<spec>__P0xb_standin` | none |
| P2, P3 | `delta_smoke_branch_opt/designs/P2`, `P3` (not the WS7 runs, which had `alt_standin` off) | `<spec>__P2_standin`, `<spec>__P3_standin` | none |
| P0y | `delta_pool_bicyclohexyl/designs/P0y` | `<spec>__P0y_standin` | committed run had `strict_constraints=False` and was declared invalid as a design (flash / NHOC model unavailable) |
| P0y-b | `delta_pool_bicyclohexyl/designs/P0yb` | `<spec>__P0yb_standin` | none (the overlay was on there by env var; now default) |

Committed headline readings (for reference; quoted from the committed `r_designs.csv` files):

| pool | Jet A: passed / min slack | Jet A-1: passed / min slack |
|---|---|---|
| P0 | False / -0.2162 | False / -0.2244 |
| P1 | False / -0.0590 | False / -0.0717 |
| P0x-a | False / -0.0981 | False / -0.1509 |
| P0x-b | False / -0.0741 | False / -0.1063 |
| P2 (standin) | False / 0.0155 (11 uncertified) | False / 0.0156 (11 uncertified) |
| P3 (standin) | False / 0.0156 (11 uncertified) | False / 0.0156 (11 uncertified) |
| P0y | False / 0.1220 (invalid) | False / 0.0762 (invalid) |
| P0y-b | True / 0.0115 | True / 0.0058 |

## Values to confirm or refute (expectations from `docs/decisions_2026-11-15.md`, NOT results)

* P0y-b Jet A stays judge-certified (conditional only on the Engler D86 route).
* P0y-b Jet A-1 fails freezing under the 13.94 kJ/mol bicyclohexyl dHfus.

Both are expectations written before any run. The re-optimised recipes can differ from the committed ones, so either may be refuted;
the results section will say which, with the numbers.

## Outputs

* `designs.csv`: one row per pool x spec (rebuilt from the per-design files after every pool, so an interrupted run keeps finished rows).
* `designs/<pool>/<design_id>__designs.csv`, `__constraints.csv`, `__readings.csv` (the readings 4-7 above).
* `log_<pool>_<spec>.txt` (not committed if over 5 MB).

## Results

(Appended after the runs; nothing above will be edited.)

## Results (2026-11-15)

Run: 16 designs (8 pools x 2 specs), n_starts = 10, seed = 42, Engler D86, alt_standin on, strict_constraints on, synthetic cap 0.50 for P2/P3 only, registry overlay ON,
designated freezing route `svl_unifac`. All 16 rows have status `designed` (no `ModelUnavailable` from the strict guard, as expected). The runs were 8 processes (pool x spec)
on a clean tree (`code_dirty` False). Rows carry git HEAD `fc618b5` (9 rows) or `adb5972` (7 rows); the two commits differ only by `results/v5/BASE/delta_adopt_rescore/README.md`
(no code change). Files: `designs.csv`, `compare.csv`, `results_table.csv` (one tidy row per pool x spec with every field below), `designs/<pool>/*`, `log_*.txt`.

**Headline.** Judge-certified (every constraint `pass` or `by_construction` under the declared margins): **P0y Jet A, P0y-b Jet A and P0y-b Jet A-1** (3 of 16).
The other 13 are kept as readings: 8 fail at least one constraint (P0, P0x-a and P0x-b on both specs; P0y Jet A-1; P1 Jet A-1) and the other 5 (P1 Jet A; P2 and P3 on both specs) are
uncertified by construction (represented-blendstock margins not validated, non-verified constants). All three certified designs carry the labels `d86_model_dependent`,
`freezing_route_dependent` and `overlay_dependent` (conditions table), so none is an unconditional certification. Their min slacks (0.0027 to 0.0132) are the optimiser sitting on the
margin boundary, so a small positive slack is the expected outcome and is not a robustness measure.

### Readings 1-2: judge status, min slack, binding and failing constraints

P/F/U = pass / fail / uncertified counts; binding = slack <= min + 0.003; "uncertified" is a count of keys (full lists in `designs.csv`).

| pool | spec | certified | min slack | P/F/U | binding | failing | uncertified |
|---|---|---|---|---|---|---|---|
| P0 | astm:jet_a | N | -0.2047 | 6/5/0 | flash_point_C, distillation_T90_minus_T10_C | flash_point_C, freezing_point_C, aromatics_wt_pct, smoke_point_mm, distillation_T90_minus_T10_C | - |
| P0 | bis:jet_a1 | N | -0.2134 | 6/3/1 | flash_point_C, freezing_point_C, distillation_T90_minus_T10_C | flash_point_C, freezing_point_C, distillation_T90_minus_T10_C | 1 keys |
| P0x-a | astm:jet_a | N | -0.0569 | 3/3/4 | flash_point_C, freezing_point_C, smoke_point_mm | flash_point_C, freezing_point_C, smoke_point_mm | 4 keys |
| P0x-a | bis:jet_a1 | N | -0.1363 | 6/4/0 | flash_point_C, freezing_point_C, smoke_point_mm, distillation_T90_minus_T10_C | flash_point_C, freezing_point_C, smoke_point_mm, distillation_T90_minus_T10_C | - |
| P0x-b | astm:jet_a | N | -0.0567 | 3/3/4 | flash_point_C, freezing_point_C, smoke_point_mm | flash_point_C, freezing_point_C, smoke_point_mm | 4 keys |
| P0x-b | bis:jet_a1 | N | -0.0980 | 3/4/3 | flash_point_C, freezing_point_C, smoke_point_mm | flash_point_C, freezing_point_C, smoke_point_mm, distillation_T90_minus_T10_C | 3 keys |
| P0y | astm:jet_a | Y | +0.0132 | 10/0/0 | flash_point_C, freezing_point_C, nhoc_MJ_kg, smoke_point_mm, distillation_T50_minus_T10_C | - | - |
| P0y | bis:jet_a1 | N | -0.0807 | 9/1/0 | distillation_T50_minus_T10_C | distillation_T50_minus_T10_C | - |
| P0y-b | astm:jet_a | Y | +0.0132 | 10/0/0 | flash_point_C, freezing_point_C, nhoc_MJ_kg, smoke_point_mm, distillation_T50_minus_T10_C | - | - |
| P0y-b | bis:jet_a1 | Y | +0.0027 | 10/0/0 | flash_point_C, freezing_point_C, nhoc_MJ_kg, smoke_point_mm, distillation_T50_minus_T10_C | - | - |
| P1 | astm:jet_a | N | +0.0152 | 0/0/11 | flash_point_C, nhoc_MJ_kg, distillation_T50_minus_T10_C, distillation_T90_minus_T10_C | - | 11 keys |
| P1 | bis:jet_a1 | N | -0.0519 | 0/4/7 | flash_point_C, freezing_point_C, distillation_T50_minus_T10_C | flash_point_C, freezing_point_C, distillation_T50_minus_T10_C, distillation_T90_minus_T10_C | 7 keys |
| P2 | astm:jet_a | N | +0.0155 | 0/0/11 | nhoc_MJ_kg, aromatics_wt_pct, distillation_FBP_C, distillation_T50_minus_T10_C | - | 11 keys |
| P2 | bis:jet_a1 | N | +0.0155 | 0/0/11 | nhoc_MJ_kg, distillation_FBP_C, distillation_T50_minus_T10_C | - | 11 keys |
| P3 | astm:jet_a | N | +0.0156 | 0/0/11 | nhoc_MJ_kg, aromatics_wt_pct, distillation_FBP_C, distillation_T50_minus_T10_C | - | 11 keys |
| P3 | bis:jet_a1 | N | +0.0155 | 0/0/11 | nhoc_MJ_kg, distillation_FBP_C, distillation_T50_minus_T10_C | - | 11 keys |

### Readings 4-5: flash and freezing

Flash: prediction / applied margin / bound (= prediction - margin) (slack), limit 38 C. Freezing UNIFAC (designated): value / safe bound (= value + 13.8427 C) vs limit (slack).
`ideal_svl` variant (reported, never judged): value / bound (= value + 16.9755 C) (bound - limit; positive = exceeds the limit, i.e. would fail on that route). `n/a` = no declared margin
for the real-fuel pools (their slack is the judge's uncertified slack).

| pool | spec | flash pred / margin / bound (slack) | freezing UNIFAC value / bound vs limit (slack) | freezing ideal_svl value / bound (bound - limit) |
|---|---|---|---|---|
| P0 | astm:jet_a | 32.43 / 2.21 / 30.22 (-0.2047) | -48.40 / -34.56 vs -40 (-0.1361) | -48.60 / -31.62 (8.38) |
| P0 | bis:jet_a1 | 32.10 / 2.21 / 29.89 (-0.2134) | -50.82 / -36.98 vs -47 (-0.2132) | -51.03 / -34.05 (12.95) |
| P0x-a | astm:jet_a | 38.13 / 2.21 / 35.92 (-0.0548) | -51.57 / -37.72 vs -40 (-0.0569) | -51.54 / -34.57 (5.43) |
| P0x-a | bis:jet_a1 | 35.03 / 2.21 / 32.82 (-0.1363) | -54.44 / -40.60 vs -47 (-0.1363) | -54.54 / -37.56 (9.44) |
| P0x-b | astm:jet_a | 38.13 / 2.21 / 35.92 (-0.0547) | -51.57 / -37.73 vs -40 (-0.0567) | -51.55 / -34.57 (5.43) |
| P0x-b | bis:jet_a1 | 36.51 / 2.21 / 34.31 (-0.0972) | -56.24 / -42.40 vs -47 (-0.0980) | -56.23 / -39.26 (7.74) |
| P0y | astm:jet_a | 40.71 / 2.21 / 38.50 (0.0132) | -54.37 / -40.53 vs -40 (0.0132) | -54.30 / -37.32 (2.68) |
| P0y | bis:jet_a1 | 40.50 / 2.21 / 38.29 (0.0077) | -61.21 / -47.36 vs -47 (0.0077) | -61.20 / -44.22 (2.78) |
| P0y-b | astm:jet_a | 40.71 / 2.21 / 38.50 (0.0132) | -54.37 / -40.53 vs -40 (0.0132) | -54.30 / -37.33 (2.67) |
| P0y-b | bis:jet_a1 | 40.31 / 2.21 / 38.10 (0.0027) | -60.97 / -47.13 vs -47 (0.0027) | -61.06 / -44.08 (2.92) |
| P1 | astm:jet_a | 38.58 / n/a / n/a (0.0152) | -43.56 / n/a vs -40 (0.0889) | -43.05 / -26.08 (13.92) |
| P1 | bis:jet_a1 | 36.03 / n/a / n/a (-0.0519) | -44.56 / n/a vs -47 (-0.0519) | -43.99 / -27.02 (19.98) |
| P2 | astm:jet_a | 43.82 / n/a / n/a (0.1532) | -48.27 / n/a vs -40 (0.2068) | -47.73 / -30.76 (9.24) |
| P2 | bis:jet_a1 | 44.29 / n/a / n/a (0.1656) | -48.96 / n/a vs -47 (0.0416) | -48.45 / -31.47 (15.53) |
| P3 | astm:jet_a | 44.68 / n/a / n/a (0.1759) | -47.80 / n/a vs -40 (0.1949) | -47.25 / -30.28 (9.72) |
| P3 | bis:jet_a1 | 44.29 / n/a / n/a (0.1656) | -48.96 / n/a vs -47 (0.0416) | -48.45 / -31.47 (15.53) |

### Reading 6: D86

Engler (designated): T10 / T50 / T90 / FBP in C (T50 and T90 = T10 + the slope predictions), then the judged slopes with slack. Rayleigh (check): slopes with slack against
15 + 1.2427 (T50-T10) and 40 + 5.3057 (T90-T10); negative = below the limit after the Rayleigh margin. `rayleigh_status` is `ok` in all 16 rows.

| pool | spec | Engler T10 / T50 / T90 / FBP | Engler T50-T10 (slack) | Engler T90-T10 (slack) | Rayleigh T50-T10 (slack) | Rayleigh T90-T10 (slack) |
|---|---|---|---|---|---|---|
| P0 | astm:jet_a | 126.0 / 175.6 / 176.0 / 208.2 | 49.63 (1.0634) | 50.09 (-0.2040) | 10.13 (-6.12) | 17.59 (-27.72) |
| P0 | bis:jet_a1 | 125.4 / 174.9 / 175.2 / 207.6 | 49.47 (1.0525) | 49.71 (-0.2134) | 9.86 (-6.38) | 16.54 (-28.76) |
| P0x-a | astm:jet_a | 136.3 / 180.2 / 177.4 / 224.6 | 43.85 (1.9235) | 41.07 (0.0267) | 7.27 (-8.97) | 16.88 (-28.43) |
| P0x-a | bis:jet_a1 | 134.1 / 180.2 / 186.9 / 223.6 | 46.07 (0.8261) | 52.80 (-0.1363) | 9.53 (-6.71) | 18.92 (-26.39) |
| P0x-b | astm:jet_a | 136.4 / 179.9 / 177.4 / 226.9 | 43.56 (1.9040) | 41.04 (0.0260) | 7.27 (-8.98) | 16.88 (-28.42) |
| P0x-b | bis:jet_a1 | 138.0 / 180.3 / 177.4 / 226.2 | 42.25 (1.8169) | 39.39 (-0.0153) | 7.96 (-8.28) | 17.51 (-27.79) |
| P0y | astm:jet_a | 144.6 / 178.5 / 214.1 / 232.2 | 33.88 (0.0132) | 69.52 (0.2819) | 9.15 (-7.09) | 32.84 (-12.47) |
| P0y | bis:jet_a1 | 149.2 / 181.7 / 214.9 / 232.3 | 32.47 (-0.0807) | 65.72 (0.1869) | 9.76 (-6.49) | 32.85 (-12.45) |
| P0y-b | astm:jet_a | 144.6 / 178.5 / 214.3 / 232.2 | 33.88 (0.0132) | 69.75 (0.2875) | 9.37 (-6.87) | 33.89 (-11.41) |
| P0y-b | bis:jet_a1 | 148.0 / 181.7 / 215.0 / 232.3 | 33.72 (0.0027) | 67.00 (0.2189) | 9.26 (-6.98) | 33.81 (-11.50) |
| P1 | astm:jet_a | 170.0 / 203.9 / 228.9 / 286.8 | 33.91 (0.0152) | 58.85 (0.0152) | 21.61 (5.37) | 55.89 (10.59) |
| P1 | bis:jet_a1 | 168.8 / 201.7 / 226.5 / 285.7 | 32.90 (-0.0519) | 57.74 (-0.0126) | 22.12 (5.88) | 55.99 (10.68) |
| P2 | astm:jet_a | 170.2 / 204.1 / 276.3 / 292.1 | 33.91 (0.0155) | 106.14 (1.1974) | 21.92 (5.67) | 94.79 (49.48) |
| P2 | bis:jet_a1 | 173.4 / 207.3 / 273.8 / 292.1 | 33.91 (0.0155) | 100.43 (1.0545) | 20.70 (4.46) | 92.76 (47.46) |
| P3 | astm:jet_a | 172.6 / 206.6 / 265.2 / 292.1 | 33.91 (0.0156) | 92.55 (0.8577) | 20.97 (4.73) | 88.48 (43.18) |
| P3 | bis:jet_a1 | 173.4 / 207.3 / 273.8 / 292.1 | 33.91 (0.0155) | 100.43 (1.0545) | 20.70 (4.46) | 92.76 (47.46) |

### Readings 3 and 7: recipe (compound wt %, entries >= 0.5 wt %) and dicycloparaffin mass fraction

Dicyclo = non-aromatic, saturated, exactly 2 rings; reported only, no gate. AJFTD reference: max 30.0 %, 95th percentile 17.3 %. For P1-P3 the compound list is the GCxGC
representation of the real-fuel blendstocks plus pool members (top-level weights are in `designs.csv`, `top_weights`).

| pool | spec | recipe (wt %) | dicyclo wt % |
|---|---|---|---|
| P0 | astm:jet_a | n-undecane 26.1; 1,2,4-trimethylbenzene 25.8; n-decane 17.6; butylcyclohexane 11.6; n-dodecane 6.6; decalin_trans 6.2; isooctane 6.1 | 6.2 |
| P0 | bis:jet_a1 | 1,2,4-trimethylbenzene 25.6; n-undecane 22.1; n-decane 18.4; decalin_trans 11.2; butylcyclohexane 11.0; isooctane 6.2; n-dodecane 5.4 | 11.2 |
| P0x-a | astm:jet_a | butylcyclohexane 32.1; 1,2,4-trimethylbenzene 21.8; n-undecane 20.6; decalin_trans 11.9; n-dodecane 4.9; n-tridecane 4.5; isooctane 3.8 | 11.9 |
| P0x-a | bis:jet_a1 | butylcyclohexane 32.7; 1,2,4-trimethylbenzene 17.0; n-undecane 16.6; decalin_trans 14.0; tetralin 7.3; isooctane 5.2; n-dodecane 3.8; n-tridecane 3.5 | 14.0 |
| P0x-b | astm:jet_a | butylcyclohexane 34.2; 1,2,4-trimethylbenzene 21.8; n-undecane 20.6; decalin_trans 10.0; n-dodecane 4.9; n-tridecane 4.5; isooctane 3.8 | 10.0 |
| P0x-b | bis:jet_a1 | butylcyclohexane 28.6; decalin_trans 19.8; 1,2,4-trimethylbenzene 17.5; n-undecane 14.6; n-decane 5.4; isooctane 4.3; hexylbenzene 3.7; n-dodecane 3.2; n-tridecane 2.9 | 19.8 |
| P0y | astm:jet_a | butylcyclohexane 23.8; 1,2,4-trimethylbenzene 17.8; n-undecane 16.5; n-decane 16.2; bicyclohexyl 14.2; n-dodecane 3.8; n-tridecane 3.4; isooctane 3.0; decalin_cis 1.3 | 15.5 |
| P0y | bis:jet_a1 | decalin_cis 20.2; butylcyclohexane 19.7; bicyclohexyl 17.9; 1,2,4-trimethylbenzene 13.7; n-decane 11.5; n-undecane 9.8; isooctane 3.4; n-dodecane 1.9; n-tridecane 1.8 | 38.1 |
| P0y-b | astm:jet_a | butylcyclohexane 24.4; 1,2,4-trimethylbenzene 17.8; n-undecane 16.5; n-decane 16.0; bicyclohexyl 15.1; n-dodecane 3.7; n-tridecane 3.4; isooctane 3.0 | 15.1 |
| P0y-b | bis:jet_a1 | butylcyclohexane 30.8; decalin_trans 19.6; bicyclohexyl 18.2; 1,2,4-trimethylbenzene 14.5; n-undecane 9.9; isooctane 3.2; n-dodecane 2.0; n-tridecane 1.8 | 37.8 |
| P1 | astm:jet_a | CCCCCCCCC(C)CC 13.9; CCCCCCCCC(C)C 5.7; Cc1ccc(C(C)C)cc1 5.5; CCCCCC1CCCCC1 5.4; n-dodecane 4.8; n-undecane 4.5; CCCCCC(C)CCC 4.2; CCCCCCCC1CCCCC1 4.1; butylcyclohexane 4.1; bicyclohexyl 3.8; CCCCCCC1CCCCC1 3.7; CCCCCCCCCCCCCC(C)C 3.5; n-decane 3.2; tetralin 3.2; mesitylene 2.9; n-tridecane 2.8; decalin_cis 2.6; CCCC1CCCCC1 2.3; isooctane 2.3; CCCCCCCCC1CCCCC1 2.2; hexylbenzene 2.2; n-tetradecane 2.0; CCCCCCCCC 1.4; CCCCCCCCCC1CCCCC1 1.3; o-xylene 1.1; CCCCC(C)CCC 1.0; CCCCCCCCCCCCCCC 1.0; CCCCCCCCc1ccccc1 0.8; CCC1CCCCC1 0.8; CCCCCCCCCCC1CCCCC1 0.6; CCCCCCCCCC(C)CCCCCC 0.5; n-octane 0.5 | 6.4 |
| P1 | bis:jet_a1 | CCCCCCCCC(C)CC 12.6; n-dodecane 8.3; CCCCCCCCC(C)C 5.1; Cc1ccc(C(C)C)cc1 5.0; CCCCCC1CCCCC1 4.9; 1,2,4-trimethylbenzene 4.8; n-undecane 3.8; CCCCCC(C)CCC 3.8; CCCCCCCC1CCCCC1 3.7; butylcyclohexane 3.7; bicyclohexyl 3.4; CCCCCCC1CCCCC1 3.3; isooctane 3.2; CCCCCCCCCCCCCC(C)C 3.1; n-decane 2.9; tetralin 2.9; mesitylene 2.7; n-tridecane 2.5; decalin_cis 2.4; CCCC1CCCCC1 2.1; CCCCCCCCC1CCCCC1 2.0; hexylbenzene 1.9; n-tetradecane 1.8; CCCCCCCCC 1.3; CCCCCCCCCC1CCCCC1 1.2; o-xylene 1.0; CCCCC(C)CCC 0.9; CCCCCCCCCCCCCCC 0.9; CCCCCCCCc1ccccc1 0.7; CCC1CCCCC1 0.7; CCCCCCCCCCC1CCCCC1 0.5 | 5.8 |
| P2 | astm:jet_a | CCCCCCCCC(C)CC 22.0; 1,2,4-trimethylbenzene 19.6; CCCCCCCCCC(C)CCCCCC 10.7; CCCCCCCCCCCCCC(C)C 5.8; CCCCCC(C)CCC 5.7; CCCCCCCCC(C)C 5.6; CCCCC(C)CCC 4.7; Cc1ccc(C(C)C)cc1 1.7; CCCCCC1CCCCC1 1.7; n-undecane 1.7; n-decane 1.6; n-dodecane 1.4; butylcyclohexane 1.3; CCCCCCCC1CCCCC1 1.3; CCCCCCCCC 1.2; bicyclohexyl 1.2; CCCCCCC1CCCCC1 1.2; n-tridecane 1.1; tetralin 1.0; CCC1CCCCC1 1.0; mesitylene 0.9; CCCC1CCCCC1 0.9; decalin_cis 0.8; CCC(C)CC(C)C 0.8; n-tetradecane 0.7; CCCCCCCCC1CCCCC1 0.7; hexylbenzene 0.7; CCCCCCCCCCCCCCC 0.5 | 2.0 |
| P2 | bis:jet_a1 | CCCCCCCCC(C)CC 23.6; 1,2,4-trimethylbenzene 22.2; CCCCCCCCCC(C)CCCCCC 9.8; CCCCCCCCCCCCCC(C)C 6.0; CCCCCC(C)CCC 5.2; CCCCCCCCC(C)C 5.2; CCCCC(C)CCC 4.3; Cc1ccc(C(C)C)cc1 1.6; CCCCCC1CCCCC1 1.6; n-undecane 1.5; n-decane 1.4; n-dodecane 1.3; butylcyclohexane 1.2; CCCCCCCC1CCCCC1 1.2; CCCCCCCCC 1.1; bicyclohexyl 1.1; CCCCCCC1CCCCC1 1.1; n-tridecane 1.0; tetralin 0.9; CCC1CCCCC1 0.9; mesitylene 0.9; CCCC1CCCCC1 0.8; decalin_cis 0.8; CCC(C)CC(C)C 0.8; n-tetradecane 0.7; CCCCCCCCC1CCCCC1 0.6; hexylbenzene 0.6 | 1.9 |
| P3 | astm:jet_a | CCCCCCCCC(C)CC 21.9; 1,2,4-trimethylbenzene 12.6; CCCCCCCCCC(C)CCCCCC 9.1; CCCCCCCCCCCCCC(C)C 5.7; CCCCCCCCC(C)C 5.4; CCCCCC(C)CCC 5.3; mesitylene 4.9; CCCCC(C)CCC 4.1; Cc1ccc(C(C)C)cc1 3.8; CCCCCC1CCCCC1 2.1; n-undecane 1.9; n-decane 1.7; n-dodecane 1.6; butylcyclohexane 1.6; CCCCCCCC1CCCCC1 1.6; bicyclohexyl 1.5; CCCCCCC1CCCCC1 1.4; tetralin 1.3; n-tridecane 1.2; CCCCCCCCC 1.2; CCCC1CCCCC1 1.1; decalin_cis 1.0; CCC1CCCCC1 0.9; n-tetradecane 0.8; CCCCCCCCC1CCCCC1 0.8; hexylbenzene 0.8; CCC(C)CC(C)C 0.7; o-xylene 0.6; CCCCCCCCCCCCCCC 0.5; CCCCCCCCCC1CCCCC1 0.5 | 2.5 |
| P3 | bis:jet_a1 | CCCCCCCCC(C)CC 23.6; 1,2,4-trimethylbenzene 22.2; CCCCCCCCCC(C)CCCCCC 9.8; CCCCCCCCCCCCCC(C)C 6.0; CCCCCC(C)CCC 5.2; CCCCCCCCC(C)C 5.2; CCCCC(C)CCC 4.3; Cc1ccc(C(C)C)cc1 1.6; CCCCCC1CCCCC1 1.6; n-undecane 1.5; n-decane 1.4; n-dodecane 1.3; butylcyclohexane 1.2; CCCCCCCC1CCCCC1 1.2; CCCCCCCCC 1.1; bicyclohexyl 1.1; CCCCCCC1CCCCC1 1.1; n-tridecane 1.0; tetralin 0.9; CCC1CCCCC1 0.9; mesitylene 0.9; CCCC1CCCCC1 0.8; decalin_cis 0.8; CCC(C)CC(C)C 0.8; n-tetradecane 0.7; CCCCCCCCC1CCCCC1 0.6; hexylbenzene 0.6 | 1.9 |

### Conditions table (one row per design, as in `delta_pool_bicyclohexyl`)

Overlay column: the recipe contains bicyclohexyl or cis-decalin (>= 0.1 wt %), the only members whose flash-registry entry comes from the overlay. For P1-P3 these two compounds appear
only through the GCxGC representation of the real-fuel blendstock, so the column is not applicable. Engler column: status of the four Engler D86 constraints. Label rules as declared above.

| design | judge | overlay needed | Engler D86 | Rayleigh D86 slopes | freezing UNIFAC bound vs limit | ideal_svl bound vs limit | dicyclo wt % | label |
|---|---|---|---|---|---|---|---|---|
| P0 astm:jet_a | N (fail) | no | fail | T50-T10 10.1 (need 16.2), T90-T10 17.6 (need 45.3) | -34.56 vs -40 | -31.62 vs -40 | 6.2 | - |
| P0 bis:jet_a1 | N (fail) | no | fail | T50-T10 9.9 (need 16.2), T90-T10 16.5 (need 45.3) | -36.98 vs -47 | -34.05 vs -47 | 11.2 | - |
| P0x-a astm:jet_a | N (fail) | no | pass | T50-T10 7.3 (need 16.2), T90-T10 16.9 (need 45.3) | -37.72 vs -40 | -34.57 vs -40 | 11.9 | - |
| P0x-a bis:jet_a1 | N (fail) | no | fail | T50-T10 9.5 (need 16.2), T90-T10 18.9 (need 45.3) | -40.60 vs -47 | -37.56 vs -47 | 14.0 | - |
| P0x-b astm:jet_a | N (fail) | no | pass | T50-T10 7.3 (need 16.2), T90-T10 16.9 (need 45.3) | -37.73 vs -40 | -34.57 vs -40 | 10.0 | - |
| P0x-b bis:jet_a1 | N (fail) | no | fail | T50-T10 8.0 (need 16.2), T90-T10 17.5 (need 45.3) | -42.40 vs -47 | -39.26 vs -47 | 19.8 | - |
| P0y astm:jet_a | Y | yes | pass | T50-T10 9.2 (need 16.2), T90-T10 32.8 (need 45.3) | -40.53 vs -40 | -37.32 vs -40 | 15.5 | d86_model_dependent,freezing_route_dependent,overlay_dependent |
| P0y bis:jet_a1 | N (fail) | yes | fail | T50-T10 9.8 (need 16.2), T90-T10 32.9 (need 45.3) | -47.36 vs -47 | -44.22 vs -47 | 38.1 | - |
| P0y-b astm:jet_a | Y | yes | pass | T50-T10 9.4 (need 16.2), T90-T10 33.9 (need 45.3) | -40.53 vs -40 | -37.33 vs -40 | 15.1 | d86_model_dependent,freezing_route_dependent,overlay_dependent |
| P0y-b bis:jet_a1 | Y | yes | pass | T50-T10 9.3 (need 16.2), T90-T10 33.8 (need 45.3) | -47.13 vs -47 | -44.08 vs -47 | 37.8 | d86_model_dependent,freezing_route_dependent,overlay_dependent |
| P1 astm:jet_a | N (uncertified only) | yes | pass | T50-T10 21.6 (need 16.2), T90-T10 55.9 (need 45.3) | n/a vs -40 | -26.08 vs -40 | 6.4 | - |
| P1 bis:jet_a1 | N (fail) | yes | fail | T50-T10 22.1 (need 16.2), T90-T10 56.0 (need 45.3) | n/a vs -47 | -27.02 vs -47 | 5.8 | - |
| P2 astm:jet_a | N (uncertified only) | yes | pass | T50-T10 21.9 (need 16.2), T90-T10 94.8 (need 45.3) | n/a vs -40 | -30.76 vs -40 | 2.0 | - |
| P2 bis:jet_a1 | N (uncertified only) | yes | pass | T50-T10 20.7 (need 16.2), T90-T10 92.8 (need 45.3) | n/a vs -47 | -31.47 vs -47 | 1.9 | - |
| P3 astm:jet_a | N (uncertified only) | yes | pass | T50-T10 21.0 (need 16.2), T90-T10 88.5 (need 45.3) | n/a vs -40 | -30.28 vs -40 | 2.5 | - |
| P3 bis:jet_a1 | N (uncertified only) | yes | pass | T50-T10 20.7 (need 16.2), T90-T10 92.8 (need 45.3) | n/a vs -47 | -31.47 vs -47 | 1.9 | - |

Reading of the labels: every certified design is `d86_model_dependent` (Rayleigh T50-T10 is 9.2-9.4 C against >= 16.24 C and T90-T10 32.8-33.9 C against >= 45.31 C; far below, not borderline),
`freezing_route_dependent` (the `ideal_svl` bound exceeds the limit by 2.68, 2.67 and 2.92 C) and `overlay_dependent`. No certified design is certified only under Engler: each also depends on
the UNIFAC freezing route and on the overlay, so **none is certified under every reported route**. The Rayleigh slopes are below the thresholds for all ten P0-family designs (certified or not);
only the six P1-P3 real-fuel rows exceed them.

### Comparison to the committed readings (delta = new - committed)

Committed values come from the per-design files named in the declared comparison table; per-constraint deltas (status, prediction, slack, applied margin, model) are in `compare.csv`.
Recipe L1 = L1 distance of the compound-weight vectors (0 = same recipe, 2 = disjoint).

| pool | spec | passed (committed -> new) | min slack (committed -> new, delta) | recipe L1 |
|---|---|---|---|---|
| P0 | astm:jet_a | False -> False | -0.2162 -> -0.2047 (+0.0115) | 0.279 |
| P0 | bis:jet_a1 | False -> False | -0.2244 -> -0.2134 (+0.0110) | 0.304 |
| P0x-a | astm:jet_a | False -> False | -0.0981 -> -0.0569 (+0.0411) | 0.275 |
| P0x-a | bis:jet_a1 | False -> False | -0.1509 -> -0.1363 (+0.0146) | 0.188 |
| P0x-b | astm:jet_a | False -> False | -0.0741 -> -0.0567 (+0.0174) | 0.123 |
| P0x-b | bis:jet_a1 | False -> False | -0.1063 -> -0.0980 (+0.0083) | 0.216 |
| P0y | astm:jet_a | False -> True | +0.1220 -> +0.0132 (-0.1088) | 0.988 |
| P0y | bis:jet_a1 | False -> False | +0.0762 -> -0.0807 (-0.1569) | 0.837 |
| P0y-b | astm:jet_a | True -> True | +0.0115 -> +0.0132 (+0.0017) | 0.111 |
| P0y-b | bis:jet_a1 | True -> True | +0.0058 -> +0.0027 (-0.0030) | 0.545 |
| P1 | astm:jet_a | False -> False | -0.0590 -> +0.0152 (+0.0742) | 0.226 |
| P1 | bis:jet_a1 | False -> False | -0.0717 -> -0.0519 (+0.0197) | 0.092 |
| P2 | astm:jet_a | False -> False | +0.0155 -> +0.0155 (-0.0000) | 0.065 |
| P2 | bis:jet_a1 | False -> False | +0.0156 -> +0.0155 (-0.0000) | 0.049 |
| P3 | astm:jet_a | False -> False | +0.0156 -> +0.0156 (-0.0000) | 0.087 |
| P3 | bis:jet_a1 | False -> False | +0.0156 -> +0.0155 (-0.0000) | 0.049 |

What moved (from `compare.csv`):
* The applied freezing margin changed from 16.9755 to 13.8427 C (-3.133 C) in every row that has a margin, and every pool was re-optimised (recipe L1 0.09-0.30 for P0, P0x-a, P0x-b and P1; P2/P3 <= 0.09).
  The freezing prediction at the new optimum is 1.5-3.3 C higher (less cold) than at the committed optimum for P0 to P0y-b; recipes and route changed together, so this is not decomposed into a pure route effect.
* P0y: the committed run was `strict_constraints=False` and invalid (committed freezing prediction -100.4 / -100.6 C, flash and NHOC models unavailable), so its committed slack (+0.122 / +0.076)
  is not a valid reference; the new run is the first valid P0y reading (recipe L1 0.99 / 0.84).
* P0y-b is the only pool with a valid committed pass: Jet A 0.0115 -> 0.0132 (L1 0.11); Jet A-1 0.0058 -> 0.0027 (L1 0.55, a materially different recipe).
* P2 / P3: slack within 2e-5 of the committed value and L1 <= 0.09, as expected (nothing is certifiable there).

### Values to confirm or refute

| expectation (decisions doc, written before any run) | outcome | numbers |
|---|---|---|
| P0y-b Jet A stays judge-certified (conditional only on the Engler D86 route) | **Confirmed for the judge; "only on Engler" is too narrow** | `passed` True, min slack +0.0132 (committed +0.0115), 10 pass / 0 fail / 0 uncertified. It is also `freezing_route_dependent` (ideal_svl bound -37.33 C vs the -40 C limit, exceeds by 2.67 C) and `overlay_dependent` (bicyclohexyl 15.1 wt %), besides `d86_model_dependent` (Rayleigh T50-T10 9.4 C, T90-T10 33.9 C). |
| P0y-b Jet A-1 fails freezing under the 13.94 kJ/mol bicyclohexyl dHfus | **Refuted for the designated route (UNIFAC); true only on the ideal_svl variant** | Under the designated route the re-optimised recipe passes freezing: UNIFAC value -60.97 C, + 13.84 C margin = bound -47.13 C vs limit -47 C, slack +0.0027 (the design is certified, min slack +0.0027). On the `ideal_svl` variant the same recipe has value -61.06 C and bound -61.06 + 16.98 = -44.08 C vs -47 C, so it would fail by 2.92 C. The expectation therefore matches the ideal SvL + 16.98 C reading, not the adopted route. The recipe is re-optimised (L1 0.55 from the committed one) and carries 37.8 wt % dicycloparaffin (bicyclohexyl 18.2 + decalin_trans 19.6), above the AJFTD maximum of 30.0 %. |

Freezing numbers per route for the three certified designs (value / margin / bound / limit, C):
* P0y Jet A: UNIFAC -54.37 / 13.84 / -40.53 / -40; ideal_svl -54.30 / 16.98 / -37.32 / -40 (exceeds by 2.68).
* P0y-b Jet A: UNIFAC -54.37 / 13.84 / -40.53 / -40; ideal_svl -54.30 / 16.98 / -37.33 / -40 (exceeds by 2.67).
* P0y-b Jet A-1: UNIFAC -60.97 / 13.84 / -47.13 / -47; ideal_svl -61.06 / 16.98 / -44.08 / -47 (exceeds by 2.92).
On the UNIFAC route the freezing constraint is in the binding set of all three.

Other findings (not declared expectations): P0y Jet A is certified (its committed reading was invalid) but P0y Jet A-1 fails, on Engler T50-T10 only (32.47 C, slack -0.0807); that recipe contains
20.2 wt % cis-decalin and 38.1 wt % dicycloparaffin. Every failure and every uncertified row is kept in the tables above and in `designs.csv`.

### Caveats

* P0 and P1 comparisons mix the adopted defaults with `alt_standin` on, while the committed P0/P1 runs had it off, so their deltas are not a pure test of the adopted defaults. P0y comparisons are against an invalid committed run.
* The engine's freezing constraint still carries the label `ideal_sle_mole_fraction` in `__constraints.csv` (`model` column) although the value and the applied margin (13.8427 C) are the UNIFAC SvL ones;
  the route is recorded in the `freezing_route` and `freezing_model_*` columns of `designs.csv`. A label issue only; no number is affected.
* Certification means: passes under the model and the declared margins (no lab test, per the 2026-11-15 decision). The flash margin has reason `lab_margin_set_by_one_system`.
* R-BC-5 (the 74 C bicyclohexyl flash cluster) was not re-run and is not transferred to the new recipes (bicyclohexyl 15.1 / 18.2 wt %). The certified Jet A-1 design is at 37.8 wt % dicycloparaffin, outside the range of the 113 AJFTD GCxGC fuels (max 30.0 %).
* A one-start smoke run is not a reading; only the 16 rows with n_starts = 10 are. The Rayleigh check is reported, not a gate.
* P1-P3 rows are uncertified by construction (represented-blendstock margins not validated; `route_grid`, `route_anchored` and `coutinho_table` constants non-verified); their positive min slacks say nothing about certification.
