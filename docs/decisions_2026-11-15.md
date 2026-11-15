# User decisions 2026-11-15 (on the overnight run, docs/overnight_2026-11-15.md)

Recorded 2026-11-15 from the user's answers. **Not applied yet**: the user asked that no re-runs start until another session has added new
data. Apply after that data lands, as declared deltas, then re-score.

## Applying (2026-11-15)

- Gate: user gave the go in chat on 2026-11-15; other session's data landed as commits 25fafd9 and 1f9b89d (primary-paper calorimetry).
- New verified Tm/dHfus rows (n-nonane, 2,3-dimethylbutane, m-xylene, p-xylene) and 3 xylene eutectics (Pitzer & Scott 1943) in thermoml_sle_hc.csv -> Tasks 3 and 5.
- New vapour-pressure rows (Pitzer & Scott 1943) -> read for Task 7 only, not merged.
- New real fuels / GCxGC data: none -> Task 9 not touched by data.
- Application plan: docs/superpowers/plans/2026-11-15-apply-decisions-and-audit-fixes.md.

| # | decision | user's choice | what applying it needs (later) |
|---|---|---|---|
| 1 | freezing model | **Adopt UNIFAC-SvL** (`svl_unifac`, margin 13.84 C); ideal SvL kept as a reported variant | `configs/uncertainty_v5.yaml` key + `predict_fuel` + judge freezing route; re-derive the margin if the new data adds SLE systems; re-run design judgements |
| 2 | D86 route | **Engler only (current)**; Rayleigh stays a reported check with `d86_model_dependent` | none |
| 3 | flash-registry overlay | **Adopt fully and re-score** | merge `data/processed/v5/flash_registry_additions.csv` into the default registry (or make the env default on); regenerate frozen T2 / WS7 / P2 / P3 numbers that depend on GCxGC representative choice |
| 4 | bicyclohexyl dHfus for SvL | **Fusion + transition, 13.94 kJ/mol** | `tm_dhfus_verified.csv` `svl_dhfus_choice` for bicyclohexyl -> total; P0y-b Jet A-1 then fails freezing (R-BC-3) |
| 5 | Coutinho lookup table | **Fix and re-score** | correct SMILES matching (112 misassigned aromatics, `results/v5/BASE/delta_coutinho_lookup_audit/misassigned.csv`); regenerate results it feeds |
| 6 | optimiser guard | **`strict_constraints=True` by default** | flip the default in `blendstock_v5.design_pool`; check tests |
| 7 | cis-decalin MW | **Verify** (C10H18, 138.25 g/mol, formula-derived) | pure-table tier -> verified; re-run pool P0y as a declared follow-up |
| 8 | push | **Pushed as is** (e92d5c6..4a7b33d to origin/v5, 2026-11-15) | done |

Expected consequence (from the committed readings, before any re-run): with decisions 1 + 3 + 4, P0y-b Jet A stays judge-certified
(conditional only on the Engler D86 route, which decision 2 keeps designated); P0y-b Jet A-1 fails freezing under the 13.94 kJ/mol dH
(UNIFAC + 13.84 not yet computed with that dH). To be confirmed by the re-run, not quoted as a result.
