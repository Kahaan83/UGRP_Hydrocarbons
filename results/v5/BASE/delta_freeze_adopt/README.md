# delta_freeze_adopt: re-derive both lab-liquidus freezing margins (ideal SvL and UNIFAC SvL) on the new constants and data

Declared 2026-11-15 BEFORE any reading in this folder was computed. Base: v5 at commit 0b248e6 (Tasks 1-4 of
`docs/superpowers/plans/2026-11-15-apply-decisions-and-audit-fixes.md` applied). Script: `scripts/v5/base/delta_freeze_adopt.py`.
Nothing in `freezing_point_C.lab_liquidus.svl_verified`, `delta_freeze_margin/`, `delta_freeze_unifac/`, any gate, spec limit or spec
verification flag is changed. Routes are not touched (that is Task 6).

## Why

The user chose UNIFAC-SvL (docs/decisions_2026-11-15.md #1) on the reading 13.84 C (UNIFAC) vs 16.98 C (ideal) from `delta_freeze_unifac`
R-FU-1. Since then the constants changed (bicyclohexyl SvL dHfus 6.86 -> 13.94 kJ/mol, total incl. the II->I transition, decision 4;
four new verified rows: n-nonane, 2,3-dimethylbutane, m-xylene, p-xylene from Finke 1954 / Douslin 1946 / Pitzer & Scott 1943) and three
measured xylene eutectics (Pitzer & Scott 1943 Table X) were added to `data/processed/v5/thermoml_sle_hc.csv` (commit 1f9b89d).
Both margins are re-derived before the user's choice is written to the yaml.

## Orientation done before this declaration (scope only, no prediction or residual was computed)

* Old SLE file (`git show 1f9b89d~1:data/processed/v5/thermoml_sle_hc.csv`): 849 rows, 838 liquidus + 11 eutectic. New file: 852 rows,
  838 liquidus + 14 eutectic. The three new rows are `point_kind = eutectic`, so `delta_freeze_margin.scope_rows` (liquidus only) would drop
  them silently.
* Liquidus rows in S_F: 67 with the old tmap (d9309d3) and 67 with the new tmap. The new verified constants pull in no new liquidus rows
  (no in-scope SLE system contains n-nonane or 2,3-dimethylbutane).
* Eutectic rows meeting the scope rule (C8-C16, all HC, all verified, T <= 273.15 K) with the new tmap: 7. Four are ThermoML eutectics of
  systems already in S_F (cbc85653bf3b naphthalene + trans-decalin 241.98 K; 26a166195d71 naphthalene + cis-decalin 225.74 K;
  898f827ffb96 dodecane + tetradecane 260.55 K; 3cb8853b2842 tridecane + pentadecane 267.45 K), three are the new Pitzer & Scott rows.
  The 3cb8853b2842 eutectic is identical (same system_id, T_K and x_liquid) to a liquidus row already in S_F.

## Scope S_F (definition, unchanged rule) and the three sub-readings

Scope rule (as `delta_freeze_unifac` / `delta_freeze_margin.scope_rows`): all components hydrocarbons, all C8-C16, every component has a
verified Tm and dHfus in `build_tmap_verified()`, T <= 273.15 K. One row per point from `scope_rows` (no `add_old` merge, which made the 10
duplicate rows found overnight; 67 rows on the old data). Scored with the unchanged `delta_newdata_freeze_smoke.summarize`.

**Eutectic points enter S_F in (C).** A eutectic is the lowest point of the liquidus: at the eutectic composition the Schroeder-van Laar
prediction (maximum over components of the solubility temperature) is the liquidus temperature there, so a measured eutectic is scored
exactly like a liquidus point. The three Pitzer & Scott eutectics meet the scope rule (C8 aromatics, 209-220 K, all three xylenes
verified, the last two only since 25fafd9), so excluding them would be a scope change, not keeping one. The rule is applied uniformly:
every eutectic row of the file that meets the scope rule enters (C), the ThermoML ones as well as the primary-paper ones, except a
eutectic row identical (system_id, T_K, x_liquid) to a liquidus row already in scope (one row per point; drops the 3cb8853b2842 row).
Mechanism: the script passes a temporary copy of the SLE csv, with eutectic rows relabelled `liquidus` and the original kind kept in
`orig_kind`, through the unchanged `FMG.scope_rows` (`FMG.SLE_CSV` patched per call). The Pitzer & Scott rows come from a primary paper,
not ThermoML; the scope label "ThermoML liquidus" in (C) therefore reads "ThermoML liquidus and eutectic points plus Pitzer & Scott 1943
eutectics".

* **(A) attribution check, old data, old constants.** Old SLE file (before 1f9b89d), liquidus only, tmap built from
  `git show d9309d3:data/raw/sle_constants/tm_dhfus_verified.csv` (bicyclohexyl 6.86, no n-nonane / 2,3-DMB / m- / p-xylene).
  Asserted: 67 rows, 9 systems. Must reproduce `delta_freeze_unifac` R-FU-1: ideal margin 16.9755, UNIFAC 13.8427, MAE 3.3271 / 2.6238
  (to 1e-3). If it does not, the run stops and B vs C is not reported as attribution.
* **(B) same 67 rows, new constants.** The (A) row identities (the in-scope mask of (A)), predicted with the current
  `build_tmap_verified()`. Isolates the constant change (bicyclohexyl 13.94 is the only one that touches these rows).
* **(C) full new scope, new constants (HEADLINE).** Current SLE file, liquidus + eutectic as above, current tmap. Expected from the
  orientation: 73 rows, 12 systems; the new systems are a44fe03d18dd (m-xylene + p-xylene), 74aede3c6579 (o-xylene + m-xylene),
  d7d3558e3a17 (o-xylene + m-xylene + p-xylene, the first ternary in S_F). n_rows and n_systems are reported as computed.
* **(C-ps) descriptive.** (B) rows + the three Pitzer & Scott rows only (no ThermoML eutectics), to separate the effect of the new paper
  data from the uniform eutectic rule. Not a gate.

## Declared readings

* **R-FA-1 (on A, B, C).** Ideal (`svl_verified`, `ideal_sle_C`) and UNIFAC (`svl_unifac`: `FM.svl_activity_C(x, consts,
  make_unifac_gamma_fn(smiles))`, the calls of `delta_freeze_unifac.py:40-58`) margins: conformal, one-sided `under`, alpha 0.1,
  k = ceil((n+1)(1-alpha)) over systems; MAE; bias; UNIFAC fallbacks.
  Note: for 9 <= n <= 18 systems k = n, so the margin is the worst single system. Adding systems can keep or raise a margin, never lower it
  below the worst old system.
* **Preference rule (unchanged from `delta_freeze_unifac`), on (C).** UNIFAC is designated only if its margin < the ideal margin AND its
  MAE <= the ideal MAE AND 0 UNIFAC fallbacks in S_F.
* **Stop rule.** If the rule does not prefer UNIFAC on (C), or n_systems < 9, `readings.json` gets a `stop_reason`, no yaml key is written,
  and the plan stops at Task 6 for the user (the decision was made on 13.84 vs 16.98; a reversal is theirs).
* **Yaml (only if the rule prefers UNIFAC).** New key `freezing_point_C.lab_liquidus.svl_unifac`, written only by the builder function
  `add_freezing_unifac(doc)` in `scripts/v5/base/build_uncertainty_v5_q2.py` from `summary.csv` row `S_F, svl_unifac` (no hand-typed
  number). The existing `svl_verified` key stays byte-identical (its 16.98 value is not re-derived into the yaml here; the (C) ideal margin
  is reported beside it).
* **R-FA-2 (descriptive).** Both models on `thermoml_all_verified_any_T_any_nC` (all-HC, all-verified rows, any T and carbon number, as
  R-FU-2, 40 systems on the old data), liquidus + eutectic with the same rule, current file and constants. Not a gate.

## Outputs

`summary.csv` (columns as `delta_freeze_unifac/summary.csv`; scopes `S_F` = (C), `S_F_A_old`, `S_F_B_newconst`, `S_F_Cps`,
`thermoml_all_verified_any_T_any_nC`), `rows.csv` (C rows and the R-FA-2 rows, both predictions, fallback flag, `orig_kind`),
`rows_A.csv`, `rows_B.csv`, `per_system_S_F.csv`, `readings.json` (keys `margin_ideal_C`, `margin_unifac_C`, `mae_ideal_C`, `mae_unifac_C`,
`n_systems`, `n_rows`, `unifac_fallbacks`, `rule_prefers_unifac`, `stop_reason` for (C), plus the A / B / C-ps sub-readings).

Status: **exploratory / dev** (the ideal residuals were seen before the `delta_freeze_unifac` declaration; no locked SLE set exists).
