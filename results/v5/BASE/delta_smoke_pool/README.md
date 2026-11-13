# delta_smoke_pool: a pool whose constants are verified, to see what else blocks smoke certification (WS1 follow-up)

Written 2026-11-14 BEFORE any run. Follows `results/v5/BASE/delta_smoke_cert/README.md` (CLM-49): pools A/B used toluene / hexylbenzene whose Tm/dHfus are not verified, so freezing was `uncertified` for a reason that is not about the design.

## Verified-constants check (done before choosing, reading tables only, no design run)

* Freezing constants come from `data/raw/sle_constants/tm_dhfus_verified.csv`; any pool member absent from it uses the Coutinho table and makes freezing `uncertified` (`library_route.build_route_models`, `margin_missing_reason`). Aromatic rows in that csv: **1,2,4-trimethylbenzene (list A, Tm 229.33 K, dHfus 13.19 kJ/mol)**, benzene, naphthalene and methyl/dimethylnaphthalenes (list B). Toluene, hexylbenzene, mesitylene, xylenes, ethylbenzene, tetralin are NOT in it.
* 1,2,4-trimethylbenzene is already in v5pool; it is the verified aromatic. Naphthalenes are not used (they would feed the naphthalenes <= 3.0 vol % alternative-limit requirement). If one wants a different aromatic, none other qualifies; so the best available aromatic IS the verified one. Its `pure_props_v5` flash status is `disputed`; the registry flash route (LFL overlay) is what the model uses and `build_route_models` shows no non-verified constant for the pool (checked: `margin_missing_reason` empty for every constraint model except the structurally excluded ron/mon/dcn/dvpe).
* Iso-paraffin: isooctane (2,2,4-trimethylpentane) has verified rho15, MW, Tm, dHfus (in the csv, list B) and registry flash; isododecane (rho15 disputed) and 2,2,4,4,6,8,8-heptamethylnonane (Tm, dHfus disputed) are not usable with verified constants.

## Declared pool and run

* **Pool C = v5pool + isooctane**: n-decane, n-undecane, n-dodecane, butylcyclohexane, trans-decalin, 1,2,4-trimethylbenzene, isooctane (`designs/poolC_ids.txt`). No toluene, no hexylbenzene (they are exactly what blocked A/B). Only one new pool: no second verified aromatic exists, and adding unverified compounds would repeat the A/B confound.
* Run (identical to the R-SC-3 protocol): `python scripts/v5/run_design_userpool_v5.py --pool results/v5/BASE/delta_smoke_pool/designs/poolC_ids.txt --spec astm:jet_a,bis:jet_a1 --synthetic-blend --alt-limits --recipe-margins --seed 42 --out results/v5/BASE/delta_smoke_pool/designs/poolC_isooctane_tmb` (n_starts default 30 for 7 components; the n_starts 8 variant applies only to the 20+-component ext pool, which is not re-run here because it has no verified iso-paraffin and was shown not to move).
* Margins `configs/uncertainty_v5.yaml` (real-fuel smoke key 4.107 mm, user decision R-SC-2b). Spec flags are not edited.

## Declared readings and gates (not moved after results)

* **Reading, not a unit test.** Per design: every constraint status, slack, margin and reason; binding constraint(s) = those with slack within 0.003 of the minimum; certified or not; the exact reason when not.
* **Expected (stated now): `uncertified`, not `pass`,** because the smoke and aromatics vol % spec limits are `secondary_sources_only` (nothing is certified under an unverified limit). Also stated now: that outcome alone would not tell us whether anything else blocks; hence the counterfactual.
* **Counterfactual column `would_certify_if_specs_verified` (NOT a certification).** A design gets True iff, for every constraint, status is `pass` / `by_construction`, or status is `uncertified` AND its slack > 0 AND, after removing from its reason the clauses that only say the spec limit is `secondary_sources_only` (and the naphthalenes `by_construction` note inside the alt-met clause), nothing else remains (no class-mix, no non-verified constant, no extrapolation, no margin gap, no 'limit violated after applying the margin'). Any `fail`, or any other remaining reason, gives False and the remaining reason is listed. The smoke primary limit failing with margin counts as satisfied only through the D1655 alternative branch exactly as the judge already encodes it (alt-met clause present and slack > 0).
* Each design is also reported with category from the runner (category B means alternative branch met at the point estimate, margins ignored; not a pass).
* Pass rule: `passed` stays what the judge says; the agent never edits `verification:` flags. Any design with `passed = True` would be reported as a surprise and re-checked.
* Table builder: `scripts/v5/base/delta_smoke_pool.py` -> `results/v5/BASE/delta_smoke_pool/r_pool.csv`.

## Amendment to the counterfactual rule (declared before the counterfactual column was read; made while fixing the table builder)

* The unverified-limit clause appears in two wordings: `secondary_sources_only` (aromatics vol %, both specs; smoke in Jet A) and `UNVERIFIED` (smoke in Jet A-1). Both are the same blocker class (a spec limit that is not verified); the builder removes both.
* The judge encodes a smoke recipe whose primary 25 mm limit fails with margin but whose D1655 alternative branch (>= 18 mm, naphthalenes <= 3.0 vol %) holds as `limit violated after applying the margin; alt_min 18 met but not certified`. The builder treats "limit violated after applying the margin" as satisfied only when `alt_min 18 met` is in the same reason (as declared: carried by the alternative branch exactly as the judge encodes it). Nothing else is removed.

## Results (run 2026-11-14; 7-component pool, 30 starts, seed 42, ~17-18 min per spec on a loaded machine)

Pool C = n-decane, n-undecane, n-dodecane, butylcyclohexane, trans-decalin, 1,2,4-trimethylbenzene, isooctane. `build_route_models` shows no non-verified constant for any design constraint model (freezing, flash, density, viscosity, smoke), so no `uncertified` here comes from pool constants. Table `r_pool.csv`, builder `scripts/v5/base/delta_smoke_pool.py`, per-design files in `designs/poolC_isooctane_tmb/`.

| | Jet A (astm) | Jet A-1 (bis) |
|---|---|---|
| recipe (wt) | BCH 0.557, n-C10 0.160, n-C11 0.136, TMB 0.093, n-C12 0.028, isooctane 0.026, decalin 0 | BCH 0.695, n-C10 0.094, n-C11 0.080, TMB 0.092, isooctane 0.026, n-C12 0.014, decalin 0 |
| density 15 C | pass, 784.0, slack 0.068 | pass, 792.8, slack 0.203 |
| flash | pass, 40.92 C (margin 2.21), slack 0.0187 | pass, 40.80 C, slack 0.0155 |
| freezing | pass, -57.7 C (margin 16.98), slack 0.0187 | pass, -64.7 C, slack 0.0155 |
| viscosity -20 C | pass, 2.90 cSt, slack 0.517 | pass, 2.91 cSt, slack 0.515 |
| NHOC | pass, 43.80 MJ/kg, slack 0.0187 | pass, 43.66, slack 0.0155 |
| aromatics wt % | by_construction (9.34) | no constraint row in Jet A-1 |
| aromatics vol % (min 8) | uncertified, 8.32, slack 0.0187, reason: spec `secondary_sources_only` | uncertified, 8.26, slack 0.0155, same reason |
| smoke | uncertified, pred 24.87 mm, margin 4.107 (class-mix gate inside), slack 0.154; primary 25 mm fails with margin, alt branch (>= 18, naphthalenes 0 by construction) holds; reason: spec `secondary_sources_only` | same, pred 24.89, spec `UNVERIFIED` |
| binding constraints (slack within 0.003 of min) | flash, freezing, NHOC, aromatics vol % | same |
| category / passed | B_meets_both_not_certifiable / False | same |
| status: pass vs uncertified | no `pass` overall; 0 `fail`; 2 `uncertified` rows (aromatics vol %, smoke) | same |
| **would certify if specs verified (COUNTERFACTUAL, not a certification)** | True | True |

Reading:
* The only thing blocking certification of either design is the verification flag of two spec limits (aromatics vol % min 8, smoke point). The class-mix gate holds (isooctane 2.6 wt %, cycloalkane 0.557 <= 0.709, n-paraffin 0.32 <= 0.508), freezing/flash/density/viscosity/NHOC pass with their margins, no constant is non-verified. That answers the question WS1 left open: with a verified-constants pool, nothing other than the spec flag blocks.
* Caveats of the counterfactual (not to be read as certification): (1) the smoke result needs the D1655 alternative branch (>= 18 mm), which is itself on an unverified limit; the primary 25 mm limit fails with the 4.107 mm margin (pred 24.87 < 25). (2) Slacks are thin: 0.016-0.019 at flash, freezing, NHOC and aromatics (an equalised max-min, so they are tied by construction, not a safety reserve). (3) Aromatics sit at 8.3 vol %, at the very edge of the 8 vol % minimum and the 8-23 vol % smoke regression domain; flash margin 2.21 C is the lab margin set by one system (`lab_margin_set_by_one_system`), not a real-fuel margin. (4) Freezing margin 16.98 C is measured on binaries with verified constants. (5) The optimum is from 30 starts of a graded SLSQP objective; `n_converged` is 15/30 (Jet A) and 10/30 (Jet A-1) and is a flag, not proof of a global optimum.
* Spec flags untouched. Needs user decision: whether to upgrade the `verification` flag of the aromatics-vol-% and smoke limits (Task 1.0); until then 0 of 14 designs (12 earlier + these 2) are certified.
* No other verified aromatic exists in the constants table (toluene, hexylbenzene, mesitylene, xylenes, ethylbenzene and tetralin are Coutinho-only); 1,2,4-trimethylbenzene is both the verified aromatic and the one the earlier v5 pool already contained. The "iso-paraffin" role is isooctane only (isododecane rho15 disputed, heptamethylnonane Tm/dHfus disputed).
