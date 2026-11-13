# delta_ws2_freeze_multicomp (WS2: freezing margin valid for multicomponent recipes)

Declared 2026-11-14 before any script of this delta exists; this commit contains no result file. Frozen reference: branch v5 commit
`50f4333`. Plan: `docs/superpowers/plans/2026-11-14-v5-ws2-ws3-freezing.md`. The declared text below is never edited after results; results
are appended under "Result".

## Why
The freezing margin in `configs/uncertainty_v5.yaml` (`freezing_point_C.lab_liquidus.svl_verified` = 16.98 C, direction under) was measured
on 9 ThermoML liquidus systems that are all binaries (`delta_freeze_margin`) and is applied to multicomponent design recipes by transfer.
This delta asks whether multicomponent liquidus data exist to measure it per class mix.

## Data (inventory found before any scoring)
* `data/processed/v5/thermoml_sle_hc.csv`: liquidus rows with n_comp >= 3: two systems (1-methylnaphthalene + decalin + hexatriacontane;
  1,4-dimethylbenzene + hexane + hexadecane). No system with >= 4 components.
* `data/raw/freezing_lit/oakley2017_liquidus.csv`: one ternary (hexane + p-xylene + hexadecane, 260.80 K) = the same paper (DOI
  10.1016/j.jct.2016.12.030) and point as the ThermoML ternary; it is removed as a duplicate (same system_id, |dT| < 0.5 K, |dx| < 0.005).
* `data/raw/freezing_lit/ascent065A_freezing.csv` (Bell, Boehm, Heyne; ASCENT 065A): `kind == multi` rows: nC16/nC13/iso-octane (2
  compositions), nC16/nC13/nC8, nC16/bicyclohexyl/iso-octane, the 8-component n-alkane fuel (C9-C16). Basis: mole (as recorded there).
* NEW: Shariatrad, Javanmardi, Rasoolzadeh, Mohammadi, ACS Omega 7 (2022) 16928, doi 10.1021/acsomega.1c07072 (open access, PMC9134231):
  Table 1, 30 wax disappearance temperatures of n-C11 + n-C14 + n-C16 + n-C18 (weight fractions; converted to mole fractions with RDKit
  MW), visual method, step heating, U(T) = 0.1 K. One system (one component set). Ingested to `data/raw/wdt_quaternary_shariatrad2022/`.
* Searched, nothing new: Bell et al. 2025 Energy Fuels full text (WSU repository) tabulates only pure compounds and the n-alkane fuel already
  in ASCENT 065A (mixture data in figures only; not digitised). Affens et al. 1981 (NRL model jet fuels): no digital full text.
* Vozka 2018 Jet A / HEFA blends are real fuels with reconstructed composition, not defined multicomponent systems: excluded from M.

## Declared choices
* **Set M.** Liquidus / wax-disappearance rows of systems with n_comp >= 3, defined composition, all hydrocarbons, from the four sources
  above, duplicates removed as stated. Eutectic rows excluded. Group = `system_id` (md5 of the sorted InChIKey set). Independent system =
  one `system_id` (30 Shariatrad compositions = 1 system).
* **Class mix of a system** (`class_mix()`), from the components present (any fraction > 0): `alkane_only` (only n-/iso-alkanes),
  `alkane_naphthene` (alkanes and cycloalkanes, no aromatic), `with_aromatic` (any aromatic ring present).
* **Constants (tier rule, fixed now).** Per component: `verified` (`tm_dhfus_verified.csv`, `svl_dhfus_kJ_mol`) else `thermoml_measured`
  (`data/processed/v5/thermoml_fusion_pure_hc.csv`, delta_ws3_pure_gaps); `ml_estimated` is never used. A row with a component lacking
  both is `no_constant`, counted, not dropped. Sensitivity: verified-only constants.
* **Models**, all on the same rows and constants:
  1. `svl_ideal` (designated): ideal Schroeder-van Laar, pure solids, mole fraction, `fuelblend.eval.composition_route.ideal_sle_C`.
  2. `svl_unifac`: per component solve ln(x_i gamma_i(T, x)) = -dH_i/R (1/T - 1/Tm_i) for T (liquid gamma from original UNIFAC,
     `fuelblend.baselines.flash_unifac.make_unifac_gamma_fn`, at the stored liquid composition), maximum over components. A component whose
     UNIFAC groups cannot be assigned falls back to gamma = 1 and the row is flagged.
  3. `ideal_solid_solution_nalkanes`: a labelled stand-in, NOT predictive UNIQUAC: the n-alkanes form one ideal solid solution (liquidus T
     where sum_i x_i / x_i_sat(T) = 1 over the n-alkanes), other components crystallise pure (SvL); maximum of the two.
  4. `coutinho_empirical_depression`: the existing `fuelblend.thermo.coutinho_wax.predict_freezing_point_coutinho` (ideal WAT minus a
     hand-set eutectic depression). It is NOT Coutinho's predictive UNIQUAC and is labelled so. Called with the same constants (thermal map
     built from the tier rule) and mole fractions.
  * **Predictive UNIQUAC solid phase (Coutinho): not implemented.** A correct version needs solid-solution interaction energies from
    sublimation enthalpies, rotator-phase handling and a multiphase solid split; a quick version would give an unverified number. It cannot
    change the outcome of the gate below, which is limited by the number of systems, not by the model.
* **Margin**, only where the gate passes: `fuelblend.uq.conformal.one_sided_margin(resid, groups, alpha=0.10, direction="under")`,
  resid = predicted - measured (C), group = `system_id`, designated model only.

## Declared readings
* **R-W2-1:** number of independent systems in M per class mix (after dedup and the constant rule).
* **R-W2-2 (gate, per class mix): >= 8 independent multicomponent systems.** Where it passes, a class-mix margin is computed and PROPOSED
  in this README as "needs user decision" (never written to the yaml here). Where it fails, the margin for that class mix stays the binary
  margin, applied by transfer, and is stated as "transferred, not validated on multicomponent data".
  **Expected outcome stated now:** alkane_only about 4 systems, alkane_naphthene 1, with_aromatic 2; all three fail; no new margin.
* **R-W2-3 (descriptive):** per model and class mix: n rows, n systems, MAE, bias, fraction under-predicted, largest under-prediction (the
  unsafe side against a maximum limit), per-system MAE. A model comparison is reported, no model is designated or re-designated.
* **R-W2-4 (descriptive):** does the transferred 16.98 C margin cover every row of M (pred + 16.98 >= measured)? Reported as counts; it is
  not a validation (too few systems).

## Declared addendum (2026-11-14, after the first README commit, before any script or result)
A fifth already-ingested source was found after the README above was committed: `data/processed/v2/wax_dauphin_coutinho.csv` (Coutinho
SoftX 2018 database, github.com/ElsevierSoftwareX/SOFTX_2018_207, Apache 2.0; listed in `docs/data_sources.csv`). It holds 4 multicomponent
systems (Dauphin et al. 1999 BIM mixtures: n-decane solvent + n-C18...C36 waxes, 7-20 components, one wax appearance temperature each) and
53 binary n-heptane + heavy-wax points (binaries: not in M). Decision fixed now: the 4 multicomponent systems ARE included in M, class mix
`alkane_only`, with method label `wax_appearance` (appearance on cooling, not the equilibrium liquidus; `delta_freezing` found ideal SvL
under-predicts these). The gate counts them like any other system; summary rows are also given without them (`excluding_wax_appearance`)
so a reader can see how much of a count rests on C18-C36 waxes far above the jet range. Expected count changes to alkane_only about 8
(if every Dauphin component has a verified or thermoml_measured constant; otherwise fewer). If alkane_only reaches the gate only through these
waxy systems, the proposed margin is reported with that caveat and remains "needs user decision".

## Outputs
`scripts/v5/base/delta_ws2_freeze_multicomp.py` -> `rows.csv` (one per scored row and model), `summary.csv`, `systems.csv`, `readings.json`.
Tests `tests/test_ws2_freeze_multicomp.py`.
