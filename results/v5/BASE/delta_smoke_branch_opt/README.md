# delta_smoke_branch_opt: optimiser sees the smoke alternative branch the judge reports; wider-boiling pure pools (P0x)

Declared 2026-11-15 BEFORE any run in this folder. Base commit 9f17534 (v5). Origin: audit `docs/audit_2026-11-15.md` section 3 (diagnostic, not a
declared reading) and section 5.1. User decisions 2026-11-15, quoted from the reply to the audit: "Yes delete", "No lab test is happening",
"Yes verified" (hexylbenzene Tm/dHfus; values to be confirmed, see P0x-b), "YES" (declared optimiser smoke fix), "YES" (commit).

Gates, margins and spec flags are not changed. The judge (`blendstock_v5.judge` / `engine_v4.evaluate_constraint`) is not changed. Nothing here can
certify a recipe that the unchanged judge would not certify.

## What changes (optimiser only)

The audit found that `engine_v4._Sides._build_alt` returns `None` (alternative branch dead) when an `alt_requires` requirement has no margin, e.g. the
naphthalenes of a GCxGC-represented blendstock. The optimiser then scores smoke on the 25 mm primary limit only, which the bulk regression cannot reach
at >= 8 vol % aromatics (24.97 mm at 8 vol %). The max-min optimum then pushes aromatics to 7.986 vol % (slack -0.00082), and the judge reports
aromatics `fail`. The judge itself computes the alternative-branch slack at the point estimate for such a requirement and reports `uncertified`.

New flag `alt_point_standin` (default **False** = existing behaviour, bit-identical):
* In `_Sides`, with the flag on, a requirement whose margin is missing (yaml margin None, or `margin_fn` returning None) uses margin 0 in the
  optimiser's alternative-branch slack instead of killing the branch. This is the same slack the judge reports for that requirement.
* `blendstock_v5.design_pool(..., alt_standin=False)` passes it through. Off by default everywhere; only the runs below turn it on.

## Declared readings

* **R-SB-1 (feature off unchanged; unit tests).** With the flag off, `_Sides.vector` is identical (bit for bit) to the current code on pool C and on the
  recorded P2 ASTM recipe; the existing WS7 feature-off tests still pass. With the flag on, on the recorded P2 ASTM recipe, the optimiser's
  smoke slack equals the slack the judge reports for `smoke_point_mm` (to 1e-9). A flag-on pool with no alternative requirement at all is unchanged.
* **R-SB-2 (P2/P3 re-run, flag on).** Pools P2 and P3 x ASTM Jet A and BIS Jet A-1, the WS7 protocol unchanged otherwise: D86 on, 10 starts, seed 42,
  lambda 0.25, synthetic cap 0.50, margins `configs/uncertainty_v5.yaml`. Reported: every constraint (status, prediction, slack, reason), binding
  constraints, certified or not. **Stated expectation:** aromatics vol % >= 8 and density (ASTM) / freezing (BIS) meet their limits at the point
  estimate, so no `fail` row remains; every model-predicted constraint stays `uncertified` (represented-blendstock rule, non-verified representative
  constants); **0 of 4 certified**. A `fail` that remains is reported as such.
* **R-SB-3 (P0x-a, pure pool, needs no new verification).** Pool C + `tetralin` + `n-tridecane` (both have verified freezing constants in
  `tm_dhfus_verified.csv`, a pure grid and flash-registry data). D86 on, 10 starts, seed 42, lambda 0.25, flag on (no effect expected: naphthalenes are
  by construction in a pure pool). **Stated expectation:** T90-T10 still below 40 + 18.25 C because no verified-freezing component boils above about
  235 C in more than about 2 mol %; reported pass or fail as computed. Certified only if the unchanged judge certifies.
* **R-SB-4 (P0x-b, pending).** P0x-a + `hexylbenzene`, run ONLY after the user's verified hexylbenzene Tm and dHfus values (and source) are
  written to `data/raw/sle_constants/tm_dhfus_verified.csv` with the user's words quoted. Same protocol as R-SB-3. Not run before that.

## Results

(appended 2026-11-15 after the runs; nothing above was edited)

Run details: `scripts/v5/base/delta_smoke_branch_opt.py`. Logs are in `logs/`, tables in `designs/<pool>/` and `r_designs.csv`. Optimiser
convergence (`opt_n_converged`, out of 10) is a flag, not proof of a global optimum: P2/P3 10/10; P0x 0 or 1/10.

### R-SB-1 (unit tests, `tests/test_delta_smoke_branch_opt.py`): pass
* The default is off and gives a vector identical to an explicit off.
* With the flag on, on the recorded P2 ASTM recipe, exactly one optimiser entry changes (`smoke_point_mm_min`, -0.000806 to +0.387769). It
  equals the judge's reported smoke slack (0.3877691) to 1e-9. The judge still says `uncertified`.
* Pool C is identical with the flag on and off.
* The WS7 feature-off tests still pass (16 passed together).

### R-SB-2 (P2/P3, flag on): matches the declared expectation
| design | top-level recipe (wt) | fail | uncertified | min slack | notes |
|---|---|---|---|---|---|
| P2 Jet A | see `r_designs.csv` | 0 | 11 | +0.0155 | aromatics 22.0 vol %, density 788.9, freezing -48.0, smoke 20.4 mm (alternative branch) |
| P2 Jet A-1 | | 0 | 11 | +0.0156 | aromatics 23.9 vol %, density 790.8, freezing -48.6, smoke 19.8 mm |
| P3 Jet A | Jet A 0.364, HEFA 0.307, ATJ 0.144, SAK 0.049, TMB 0.136 | 0 | 11 | +0.0156 | SAK now gets weight (0.0011 before) |
| P3 Jet A-1 | Jet A 0.281, HEFA 0.338, ATJ 0.162, SAK 0, TMB 0.219 | 0 | 11 | +0.0156 | |

Reading:
* The four WS7 `fail` rows at -0.00082 (density / freezing / aromatics) are gone. Every limit is met at the point estimate. Every
  model-predicted constraint stays `uncertified` (represented-blendstock rule, non-verified representative constants). **0 of 4 certified**, as declared.
* The recipes moved to high aromatics (22 to 24 vol %, near the 25 vol % maximum) and rely on the smoke alternative branch (18 / 19 mm with
  naphthalenes <= 3 vol %). The binding constraints are NHOC, FBP, T50-T10 and aromatics wt (+0.0155).
* This is a better design, not a certified one.

### R-SB-3 (P0x-a = pool C + tetralin + n-tridecane): fails, as expected
| spec | fails (slack) | T90-T10 (needs >= 58.25 with the margin) |
|---|---|---|
| Jet A | flash 36.5 C, freezing -53.1 C, smoke 20.3 mm, T90-T10 (all -0.098) | 54.3 C |
| Jet A-1 | flash 34.5, freezing -56.9, smoke 20.2, T90-T10 (all -0.151), NHOC 42.97 (-0.0005) | 52.2 C |

### R-SB-4 (P0x-b = P0x-a + hexylbenzene, user-verified 209.7 K / 20.877 kJ/mol): fails
| spec | fails (slack) | T90-T10 |
|---|---|---|
| Jet A | flash 37.4 C, freezing -54.0, smoke 21.0, T90-T10 (-0.063 to -0.074) | 55.7 C (hexylbenzene 0.02 wt) |
| Jet A-1 | flash 36.2, freezing -59.0, smoke 21.1 (-0.106); D86 `uncertified` (TBP outside the D86 fit domain) | 40.3 C |

Reading (R-SB-3/4):
* Pure-compound pools with D86 on stay structurally infeasible under the declared margins.
* Widening the boiling range needs heavy components. The verified-freezing ones are aromatics (tetralin, hexylbenzene) or n-alkanes capped at
  about 2 mol %. Heavy aromatics lower smoke below the alternative limit plus its 4.11 mm margin and raise freezing against the 16.98 C margin.
  A lighter front end to widen the curve drops flash below 38 + 2.21 C.
* The max-min optimum equalises flash, freezing, smoke and T90-T10 at a negative slack.
* **0 of 4 certified.** No gate was moved.

### Overall: 0 of 8 designs certified (P2, P3, P0x-a, P0x-b x Jet A / Jet A-1)
The remaining routes to a model-level certificate are the ones in `docs/audit_2026-11-15.md` section 5. The user ruled out a lab test on 2026-11-15
("No lab test is happening"), so the SAF pools need validated margins for represented blendstocks and verified representative constants. A pure
pool would need smaller freezing / smoke / D86 margins (new data, not a gate change) or heavy cycloalkanes with flash-registry data (e.g. bicyclohexyl).

Side finding (recorded in the hexylbenzene row note): the Coutinho table lookup returned the **benzene** entry (278.7 K, 9.87 kJ/mol) for
hexylbenzene (SMILES-substring defect, as for toluene in CLM-65). Every earlier freezing prediction with hexylbenzene (it is a GCxGC
representative in the Jet A base POSF-10325) used benzene's constants. From 2026-11-15 the verified values replace them. Earlier committed
results are not re-run here.
* Constant state per run:
  * The R-SB-2 and R-SB-3 runs used the declared base state. Hexylbenzene was not yet in the verified csv: the row was reverted while they
    ran, and each process caches the table at start-up.
  * R-SB-4 ran with the row in place. Its Jet A-1 run was concurrent with the P0x-a Jet A-1 run, which had loaded the table before the row
    was re-applied.
  * So the R-SB-2 P2/P3 freezing values still carry hexylbenzene's benzene-entry constants. All of those rows are `uncertified` either way.
