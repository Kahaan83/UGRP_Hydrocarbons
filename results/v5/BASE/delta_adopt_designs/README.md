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
