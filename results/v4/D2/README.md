# results/v4/D2 - design engine with raw Refutas and lab margins for known compositions (review fix 4)

Re-run of `scripts/v4/run_designs_v4.py --out results/v4/D2` and `scripts/v4/verify_designs_v4.py --out results/v4/D2`. `results/v4/D` is untouched. Pools, margins and rules were NOT adjusted.

What changed in the engine: design viscosity is now raw Refutas judged by `viscosity_kin_cSt.lab_mixtures.raw_refutas.<T bin>` (bins T_lt_263K | T_263_283K | T_ge_283K; edges -10 C and 10 C belong to the warmer bin); design density is ideal volume mixing judged by `density_kg_m3.lab_mixtures.ideal_volume.over` (lower limit) and `.under` (upper limit). Both are hydrocarbon-only (a component with w > 0 that is not a hydrocarbon makes the model unavailable) and not jet specific. The recalibrated Refutas (CRC margin) and the ridge remain as alternates in `robustness_v4.csv`. Both cold viscosity bins have no valid margin in the yaml (`inf`), so every constraint below 10 C is `uncertified`; no neighbouring bin or number is substituted.

## Category counts (A fails at point estimate / B meets at point estimate but not certifiable / C certified / D undetermined)

| run | A | B | C | D |
|---|---|---|---|---|
| v4 (`results/v4/D`) | 5 | 0 | 0 | 8 |
| review fix 4 (`results/v4/D2`) | 7 | 0 | 0 | 6 |

Certified designs: 0 of 13.

## Design-level category moves

| design | v4 | D2 |
|---|---|---|
| gasoline_bis_is2796_e10__v3pool | D | A |
| gasoline_bis_is2796_e10__v4pool | D | A |

## Constraints whose status changed

Caution: the optimiser is re-run, so the recipe (weights) can differ between D and D2 for the same design id; a status change is a change of the whole design evaluation, not only of the margin.

| design | constraint | old status | new status | old reason | new reason |
|---|---|---|---|---|---|
| jet_astm_d1655__v3pool | kin_visc_m20C_cSt | pass (refutas_recalibrated) | uncertified (refutas_raw) | - | viscosity_kin_cSt.lab_mixtures.raw_refutas.T_lt_263K: no valid margin (applied_margin absent: too few independent groups, see yaml entry) |
| jet_astm_d1655__v4pool | kin_visc_m20C_cSt | pass (refutas_recalibrated) | uncertified (refutas_raw) | - | viscosity_kin_cSt.lab_mixtures.raw_refutas.T_lt_263K: no valid margin (applied_margin absent: too few independent groups, see yaml entry) |
| jet_bis_is1571__v3pool | kin_visc_m20C_cSt | pass (refutas_recalibrated) | uncertified (refutas_raw) | - | viscosity_kin_cSt.lab_mixtures.raw_refutas.T_lt_263K: no valid margin (applied_margin absent: too few independent groups, see yaml entry) |
| jet_bis_is1571__v4pool | kin_visc_m20C_cSt | pass (refutas_recalibrated) | uncertified (refutas_raw) | - | viscosity_kin_cSt.lab_mixtures.raw_refutas.T_lt_263K: no valid margin (applied_margin absent: too few independent groups, see yaml entry) |
| jet_bis_is1571__v4pool | smoke_point_mm | fail (ajftd_bulk_regression) | uncertified (ajftd_bulk_regression) | smoke_point_mm.ajftd100_bulk.selected_in_fold: no valid margin (applied_margin absent: too few independent groups, see yaml entry); spec limit verification: UNVERIFIED (an unverified limit is never certified); limit violated | model emitted OutOfDomainWarning: aromatics 4.494589525768724 vol % outside training domain 8-23; predicted smoke point 26.1 mm outside trusted range 19-25 mm; smoke_point_mm.ajftd100_bulk.selected_in_fold: no valid margin (applie |
| diesel_astm_d975__v3pool | kin_visc_40C_cSt | fail (refutas_recalibrated) | pass (refutas_raw) | T=40 C outside the margin's validated range (-40.0, 20.0); limit violated after applying the margin | - |
| diesel_bis_is1460__v3pool | kin_visc_40C_cSt | fail (refutas_recalibrated) | pass (refutas_raw) | T=40 C outside the margin's validated range (-40.0, 20.0); limit violated after applying the margin | - |
| gasoline_astm_d4814__v3pool | density_15C_kg_m3 | pass_jet_derived_margin_only (ideal_volume_rho15) | pass (ideal_volume_rho15) | margin validated on ['jet'] fuels only; this is a gasoline spec, so the pass is not a certification | - |
| gasoline_astm_d4814__v4pool | density_15C_kg_m3 | pass_jet_derived_margin_only (ideal_volume_rho15) | pass (ideal_volume_rho15) | margin validated on ['jet'] fuels only; this is a gasoline spec, so the pass is not a certification | - |
| gasoline_astm_d4814__v4pool_dvpe | density_15C_kg_m3 | pass_jet_derived_margin_only (ideal_volume_rho15) | pass (ideal_volume_rho15) | margin validated on ['jet'] fuels only; this is a gasoline spec, so the pass is not a certification | - |
| gasoline_bis_is2796_e10__v3pool | density_15C_kg_m3 | pass_jet_derived_margin_only (ideal_volume_rho15) | fail (ideal_volume_rho15) | margin validated on ['jet'] fuels only; this is a gasoline spec, so the pass is not a certification | limit violated after applying the margin |
| gasoline_bis_is2796_e10__v4pool | density_15C_kg_m3 | pass_jet_derived_margin_only (ideal_volume_rho15) | fail (ideal_volume_rho15) | margin validated on ['jet'] fuels only; this is a gasoline spec, so the pass is not a certification | limit violated after applying the margin |

## Reading the table

* `kin_visc_m20C_cSt` (jet, 4 designs): pass -> uncertified. The -20 C constraint falls in the T_lt_263K bin, whose lab margin is `inf` (7 independent groups). The old pass relied on the CRC recalibration margin, which belongs to reconstructed compositions.
* `kin_visc_40C_cSt` (diesel v3 pools): fail -> pass. The old model was refused outside its -40..20 C window and the recalibration over-predicted; raw Refutas with the T_ge_283K margin (0.186 ln) passes.
* `density_15C_kg_m3` (gasoline ASTM): pass_jet_derived_margin_only -> pass, because the lab margin is not jet specific. This is a margin-scope change, not a better design.
* `density_15C_kg_m3` (E10, 2 designs): the density model is unavailable whenever ethanol has w > 0, so at the optimiser start the constraint is excluded, the optimiser drops ethanol to 0, and the resulting hydrocarbon-only recipe then violates the density minimum (fail). The gasoline E10 'designs' are therefore ethanol-free; this is the documented consequence of the hydrocarbon-only ruling.
* `smoke_point_mm` (jet BIS v4pool): fail -> uncertified is not caused by the viscosity or density change (the smoke-point model and margin are untouched); the re-run optimiser returned a different recipe (the new reason cites 4.5 vol % aromatics, outside the smoke model's 8-23 domain).

