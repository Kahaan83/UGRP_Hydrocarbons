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

(appended after the runs; nothing above is edited afterwards)
