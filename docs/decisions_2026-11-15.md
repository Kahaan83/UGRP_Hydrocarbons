# User decisions 2026-11-15 (on the overnight run, docs/overnight_2026-11-15.md)

Recorded 2026-11-15 from the user's answers. **Applied 2026-11-15** (after the other session's new data landed): see the `applied in` column and the plan
`docs/superpowers/plans/2026-11-15-apply-decisions-and-audit-fixes.md`. Originally the user asked that no re-runs start until new data was added.

## Applying (2026-11-15)

- Gate: user gave the go in chat on 2026-11-15; other session's data landed as commits 25fafd9 and 1f9b89d (primary-paper calorimetry).
- New verified Tm/dHfus rows (n-nonane, 2,3-dimethylbutane, m-xylene, p-xylene) and 3 xylene eutectics (Pitzer & Scott 1943) in thermoml_sle_hc.csv -> Tasks 3 and 5.
- New vapour-pressure rows (Pitzer & Scott 1943) -> read for Task 7 only, not merged.
- New real fuels / GCxGC data: none -> Task 9 not touched by data.
- Application plan: docs/superpowers/plans/2026-11-15-apply-decisions-and-audit-fixes.md.

| # | decision | user's choice | what applying it needs (later) | applied in |
|---|---|---|---|---|
| 1 | freezing model | **Adopt UNIFAC-SvL** (`svl_unifac`, margin 13.84 C); ideal SvL kept as a reported variant | `configs/uncertainty_v5.yaml` key + `predict_fuel` + judge freezing route; re-derive the margin if the new data adds SLE systems; re-run design judgements | `bcab8b4` (route); margin re-derived in `16070c2` (declaration) and `9fbff24` (results): 13.8427 C UNIFAC vs 16.9755 C ideal; CLM-73 |
| 2 | D86 route | **Engler only (current)**; Rayleigh stays a reported check with `d86_model_dependent` | none | no change |
| 3 | flash-registry overlay | **Adopt fully and re-score** | merge `data/processed/v5/flash_registry_additions.csv` into the default registry (or make the env default on); regenerate frozen T2 / WS7 / P2 / P3 numbers that depend on GCxGC representative choice | `5986a6b` (env `FUELBLEND_FLASH_REGISTRY_ADDITIONS=0` restores the v2 registry); T2 re-score `a5afab1` (delta_adopt_rescore), adopted as the frozen BASE in `3b75437` (U4, CLM-74) |
| 4 | bicyclohexyl dHfus for SvL | **Fusion + transition, 13.94 kJ/mol** | `tm_dhfus_verified.csv` `svl_dhfus_choice` for bicyclohexyl -> total; P0y-b Jet A-1 then fails freezing (R-BC-3) | `3ec5aa1` (constant), `daa3764` (SOURCES note) |
| 5 | Coutinho lookup table | **Fix and re-score** | correct SMILES matching (112 misassigned aromatics, `results/v5/BASE/delta_coutinho_lookup_audit/misassigned.csv`); regenerate results it feeds | `0b248e6` (exact aromatic match; substring kept as default for v4 callers) |
| 6 | optimiser guard | **`strict_constraints=True` by default** | flip the default in `blendstock_v5.design_pool`; check tests | `e93dc22` (P0y keeps an explicit opt-out to reproduce R-BC-2) |
| 7 | cis-decalin MW | **Verify** (C10H18, 138.25 g/mol, formula-derived) | pure-table tier -> verified; re-run pool P0y as a declared follow-up | `3ec5aa1` (MW verified 138.254, RDKit C10H18) |
| 8 | push | **Pushed as is** (e92d5c6..4a7b33d to origin/v5, 2026-11-15) | done | done (e92d5c6..4a7b33d pushed) |

U4 (adopt the delta_adopt_rescore re-score as the new frozen BASE; user answer 2026-11-15 "Yea adopt"): **applied in `3b75437`** (re-run `--full`,
`reproduce_base.py --update-manifest`, pins and claims; CLM-74). Pre-adoption BASE recoverable at `3ab033d`.

Expected consequence (from the committed readings, before any re-run): with decisions 1 + 3 + 4, P0y-b Jet A stays judge-certified
(conditional only on the Engler D86 route, which decision 2 keeps designated); P0y-b Jet A-1 fails freezing under the 13.94 kJ/mol dH
(UNIFAC + 13.84 not yet computed with that dH). To be confirmed by the re-run, not quoted as a result.

Reading from the application (`results/v5/BASE/delta_freeze_adopt/`): decision 4 makes bicyclohexyl + butylcyclohexane over-predict freezing point by up to 33.6 C
(safe side), system MAE 14.7 C; margins unchanged (13.8427 UNIFAC / 16.9755 ideal). The o-xylene + m-xylene eutectic composition in the Pitzer & Scott ingest
appears transposed (MAE would be 4.19 ideal / 3.65 UNIFAC with the swap; margins unchanged); a data-owner check is open.

## Answers to the open items (user)

- U1: the two untracked 2026-11-12 plan files stay untracked (not committed).
- U2: `data/raw/sak/dot_64804/dot_64804_DS1.pdf` untracked (`git rm --cached`); the file stays on disk and is covered by the `data/*` ignore rule.
- U3: the `docs/claims_register.md` line in `SHA256SUMS_v4.txt` rehashed (only that line); no other v4 artefact touched.
- U5: pushed.
- U6: dicycloparaffin bound (0-30 %) in the smoke class-mix rule: **not added, no re-run.** The user's condition was to re-run only if the bound would improve results. It cannot raise the
  certified count: it only restricts. Of the 3 certified designs, P0y Jet A (15.5 wt %) and P0y-b Jet A (15.1 wt %) are already inside 30 %; only P0y-b Jet A-1 (37.8 wt %) is outside, so the
  best case is 3 of 16 unchanged and the likely case is 2 of 16 (freezing pushes the Jet A-1 optimum towards dicycloparaffins: P0x-b Jet A-1 at 19.8 wt % misses the limit). The 37.8 wt % caveat stays documented.
- U7: `design_pool(alt_standin=...)` now defaults to True (CLM-67). The 16 designs in `delta_adopt_designs` already ran with it on, so their results are unchanged. The WS7 runner passes
  `alt_standin=False` explicitly to keep reproducing its committed protocol. Test: `tests/test_alt_standin_default.py`.
