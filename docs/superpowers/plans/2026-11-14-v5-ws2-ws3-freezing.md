# v5 WS3 (Tm / dHfus and heavy n-alkane Antoine gap filling) then WS2 (freezing margin for multicomponent recipes)

> Detailed plan for WS3 (Tm/dHfus + Antoine parts only; YSI/TSI and LFL parts are out of scope) and WS2 of
> `docs/superpowers/plans/2026-11-13-v5-certification-and-ml-roadmap.md`. Written 2026-11-14 before any result. Executed inline
> (autonomous session, user asleep). Global Constraints of the roadmap apply unchanged. Subagents: none.

**Goal.** (1) Fill missing pure-component inputs of the freezing (Tm, dHfus) and flash (heavy n-alkane Antoine) models with
measured data first and a gated estimate second, every filled value flagged and never replacing a measured one. (2) Test whether
the freezing margin, today measured on binaries and transferred to multicomponent recipes, can be measured on real multicomponent
liquidus data per class mix.

**Not in scope / never touched.** `configs/uncertainty_v5.yaml` and `configs/specs_*.yaml` (no margin or flag edits; any proposal is
"needs user decision"), `src/fuelblend/api/predict_route.py`, `src/fuelblend/design/library_route.py`, `engine_v4.py`,
`data/raw/sle_constants/tm_dhfus_verified.csv` (curated, read-only), the v2 flash registry `data/processed/v2/compounds_anchored_vp.csv`.
Nothing in the API or the design path reads any file written here.

## Delta folders (README with declared readings committed before any result)
* `results/v5/BASE/delta_ws3_pure_gaps/` - WS3 readings R-W3-*.
* `results/v5/BASE/delta_ws2_freeze_multicomp/` - WS2 readings R-W2-*.

## Data (state found 2026-11-14, before any result)
* Verified Tm/dHfus: `data/raw/sle_constants/tm_dhfus_verified.csv`, 39 compounds (WebBook calorimetry rows). Acree-Chickos compendia are
  NOT ingested (cited once, for 1,8-DMN). Not fetched in this session (very large PDF tables; ThermoML gives machine-readable measured rows).
* ThermoML v2020-09-30 archive was no longer on disk; re-downloaded to `data/raw/thermoml_archive_2020/` (189,433,115 B, sha256
  231161b5...58d7ec2, identical to AGENTS.md; public domain). Not committed.
* Multicomponent liquidus already ingested: ThermoML ternaries (2 systems), Oakley 2017 ternary (same data as the ThermoML
  1,4-dimethylbenzene + hexane + hexadecane system, DOI 10.1016/j.jct.2016.12.030), ASCENT 065A (Bell/Boehm/Heyne) mixtures.
* New open data fetched: Shariatrad et al., ACS Omega 7 (2022) 16928, n-C11 + n-C14 + n-C16 + n-C18 wax disappearance temperatures,
  30 mixtures (Table 1), PMC9134231. Bell et al. 2025 full text (WSU repository, CC BY-NC-ND) has no new tabulated mixture data beyond
  what ASCENT 065A already gave (figures only). Affens et al. 1981 (NRL, model jet fuels): no digital full text (OSTI 6034336).

## Tasks
### WS3
- [ ] T1 README `delta_ws3_pure_gaps` (readings, gates, family definition, tier rules, expected outcomes). Commit.
- [ ] T2 Antoine heavy n-alkanes: `src/fuelblend/thermo/antoine_interp.py` (pure functions: grid-log10P interpolation + refit; literal
  per-parameter interpolation comparator), `scripts/v5/base/delta_ws3_antoine.py`; tests `tests/test_ws3_antoine_interp.py`
  (synthetic homologous series recovered; knots never overwritten; flag `interpolated`). Output
  `results/v5/BASE/delta_ws3_pure_gaps/antoine_*.csv`, `data/processed/v5/antoine_interpolated_nalkanes.csv` (flag `antoine_source=interpolated`).
- [ ] T3 ThermoML pure fusion harvest: `src/fuelblend/data/thermoml_fusion.py` (pure-component `Normal melting temperature`,
  `Triple point temperature`, `Molar enthalpy of transition or fusion` with phases exactly {one Crystal, Liquid}); test with a real-record-shaped
  fixture `tests/test_thermoml_fusion.py`. Output `data/processed/v5/thermoml_fusion_pure_hc.csv` (per-row) and per-compound medians, tier
  `thermoml_measured`.
- [ ] T4 Gap table: every pool member (v5pool, ext, pool A, pool B of `delta_smoke_cert`) and every SubRepsV5 representative used for the
  AJFTD GCxGC fuels, marked `verified` / `thermoml_measured` / `needs_estimate`.
- [ ] T5 Tm / dSfus estimator (`src/fuelblend/thermo/fusion_estimate.py`): LightGBM on RDKit descriptors + symmetry + group counts; Joback
  comparator; leave-family-out GroupKFold; gate R-W3-3. Write `data/processed/v5/fusion_ml_estimated.csv` (`source=ml_estimated`, OOF
  errors, gate status) for the gap-table compounds only. Test `tests/test_ws3_fusion_estimate.py` (never replaces measured; family function).
- [ ] T6 Propagation: +-1 OOF RMSE in Tm (and dHfus) on the pool A / pool B design recipes of `delta_smoke_cert` and the components that
  needed an estimate; report the change in the SvL freezing value.
### WS2
- [ ] T7 README `delta_ws2_freeze_multicomp` (readings, class mixes, gate, tier rule, model list incl. what is not implemented). Commit.
- [ ] T8 Ingest Shariatrad 2022 Table 1 (`data/raw/wdt_quaternary_shariatrad2022/shariatrad2022_wdt.csv` + SOURCES.md; loader in
  `loaders_freezing_lit.py`-style; test checks 30 rows, weight fractions sum to 1, values vs PDF spot checks).
- [ ] T9 Multicomponent set builder + scorer `scripts/v5/base/delta_ws2_freeze_multicomp.py`: n_comp >= 3 liquidus rows from ThermoML,
  ASCENT 065A, Oakley (deduplicated against ThermoML), Shariatrad; class mix per system; models: ideal SvL (designated), SvL + UNIFAC liquid
  activity (`make_unifac_gamma_fn`), the existing empirical `predict_freezing_point_coutinho` depression (labelled as such, not p.UNIQUAC),
  an ideal n-alkane solid-solution stand-in (labelled). p.UNIQUAC solid phase: not implemented (reason in README).
  Tests `tests/test_ws2_freeze_multicomp.py` (class-mix function, SvL-UNIFAC reduces to ideal when gamma = 1, dedup).
- [ ] T10 Append results to both READMEs; claims (new CLM numbers after CLM-49; collision with parallel WS agents is possible, flagged in
  the report) and a `docs/model_decision.md` note. One full pytest at the end (baseline: 1 failed
  `tests/test_v4_sle.py::test_default_mode_unchanged`). Commit on the worktree branch, no push.

## Review focus
1. Antoine knots: a refit is a knot only if its (A,B,C) reproduces the compound's own measured ThermoML points (rms <= 0.03 log10); the
   C13/C15/C17/C18/C19 refits in `antoine_refit.csv` have no measured points in `measurements_dedup.csv` and are NOT knots.
2. A filled Tm/dHfus never replaces a verified or measured one; ml_estimated values are not read by the API or designs.
3. Leave-family-out is strict; the family function is fixed in the README before the fit.
4. Multicomponent means n_comp >= 3; duplicates across sources removed by (system_id, composition, T) before counting systems.
5. Margin wording: if the per-class-mix gate fails, the existing 16.98 C margin stays "transferred from binaries", stated as such.
