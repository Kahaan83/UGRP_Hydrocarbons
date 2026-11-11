# UGRP Hydrocarbons: Project Summary (v4 state, for external review)

Purpose of this file: give a reviewer one place to understand what the project is, what was done in v4, what the evidence says, and where it is weak. Numbers here are copied from `docs/FULL_REPORT_v4.md` and `results/v4/*`, which are the authoritative sources. Where this file and those disagree, those win. Please check the claims against them.

Branch: `v4-integrity`. Open PR: Kahaan83/UGRP_Hydrocarbons#5 (Tasks 17-19). Earlier tasks merged in PRs #1-#4.

---

## 1. Goal

Predict properties of hydrocarbon fuel blends (jet, diesel, gasoline, rocket kerosene) from composition, with a model that is trustworthy on real measured blends, then use it for inverse design (find a blend meeting property targets).

Properties: kinematic viscosity (Refutas physics baseline), flash point (Le Chatelier + Raoult + Antoine, ideal gamma baseline; one-sided safety error matters because over-prediction is the dangerous direction), pour point, density, freezing point, heat of combustion, smoke point, octane (RON/MON), vapour pressure (DVPE), distillation.

Methodological position taken in v4: a low error on synthetic blends whose labels come from the blending rule is circular. A robust model needs real measured data, a physics baseline plus a learned residual, validation that holds out whole fuels or blend systems, and calibrated uncertainty.

## 2. What went wrong before v4 (why v4 exists)

- The v3 fuel corpus and the v3 design optimiser were formula-derived or otherwise fabricated: compositions and labels came from the same formulas the models then "predicted". v3 headline numbers (including a retracted 3.10 ON octane figure and the v3 designs) are not valid evidence.
- Provenance of every input was not tracked, so circular or placeholder values could reach headline numbers.

v4 response: quarantine the fabricated data (`data/quarantine/`, `tests/quarantine/`, excluded by `pytest.ini`), add provenance classes to every input (`raw_file`, `cited_typed`, `formula_derived`, `target_derived_composition`, `placeholder_constant`, `unknown`), and re-run everything on the cleaned corpus.

## 3. Pipeline and layout (v4 parts)

```
src/fuelblend/
  baselines/        physics baselines (Refutas, Le Chatelier/UNIFAC flash, sooting/YSI, heat of combustion, SLE)
  mixture/          structured residual (ridge/GP), set encoder (attention/sum)
  models/           octane (bRON blend, Ghosh lump), density, etc.
  uq/conformal.py   split-conformal margins
  design/           engine_v4.py, models_v4.py (design on library physics)
  api/predict_blend.py   prediction API (no lazy retrain; yaml margins)
scripts/v4/         build/eval scripts per task, reproduce_all.py, verify_designs_v4.py
configs/            uncertainty_v4.yaml (GENERATED margins), specs_astm.yaml, specs_bis.yaml
results/v4/<X>/     per-task outputs + READMEs (V viscosity, B bulk props, SLE, SP smoke, O octane, U uncertainty, E encoder, D designs, F flash)
docs/               FULL_REPORT_v4.md, claims_register.md, model_decision.md, scope_decisions.md
SHA256SUMS_v4.txt   83-file manifest; reproduce_all.py --verify-only
```

Run tests with the project venv (`.venv/Scripts/python.exe`); the system Python lacks rdkit.

## 4. Data: what is real and what is not

- Real-fuel corpus: CRC Report 647 (57 aviation fuels with measured pour, freezing, flash, D86, viscosity) plus AJFTD (81 fuels, thin), RP-1 data (one characterised sample, NISTIR 6646), and pure/binary data from Leenhouts benchmark, NIST/TDE viscosity, ThermoML, Li et al. flash points, Zabetakis LFL.
- **CRC647 compositions are reconstructed, not measured** (formula-derived). Real-fuel errors therefore mix reconstruction error and physics error.
- Locked (never-tuned) test set is thin: 12 CRC fuels, and 4 octane blendstock (CAS) groups. Diesel has 0 qualifying fuels, so nothing is validated for diesel.
- Open data covers pure and binary mixtures well; realistic 3+ component hydrocarbon blends with measured properties are almost absent.
- Pure-component constants (`data/processed/v4/pure_props_v4.csv`, 411 rows) carry status counts 216 / 78 / 65 / 26 / 26; some constants are disputed or anchored to a reference (see mesitylene below).

## 5. Methods used in v4

- Splits: GroupKFold by system, leave-one-fuel-out, outer CV with selection inside training folds only (ruling R6), cluster bootstrap CIs, dev set separate from locked set. bRON octane grouped by blendstock CAS (R19). FACE rows are never targets (R17).
- Models: physics baselines; structured residuals linear in 1000/T for viscosity (bilinear, ridge, GP); group-additivity for heat of combustion; ideal Schroeder-van Laar and Coutinho for freezing/SLE; YSI mixing and a bulk regression for smoke point; blend-RON models for octane; a set encoder (attention/sum) tried and not promoted.
- Uncertainty: split-conformal margins generated into `configs/uncertainty_v4.yaml`. Consumers must read `applied_margin` only; `inf` means no valid margin and the property is uncertified. A zero-floored margin must not be used for novel blends.
- Design engine: optimises on library physics, certifies only where a valid margin exists; unverified spec limits block certification.

## 6. Results (headline, with caveats)

| Area | Result | Caveat |
|---|---|---|
| Viscosity | Dev leave-one-fuel-out log-RMSE: raw Refutas 0.427, recalibrated 0.135, bilinear 0.090, ridge 0.083, GP 0.089. Learned models beat recalibrated Refutas (paired CI excludes 0). | No pairwise ranking among bilinear/ridge/GP was tested. GP intervals under-cover (0.77 dev / 0.69 locked). A lab zero-shot model does not fix -40 C. |
| Set encoder | Lost to the ridge: 0.0862 vs 0.0833, CI includes 0. Not promoted. | Locked fuel-level coverage 8/12. |
| Flash | Lab conformal margin 3.70 C (31 systems); UNIFAC baseline checked. | Margin depends on the mesitylene system: without it 0.31 C (CLM-36). Mesitylene anchored Psat is 0.52x its reference (9 lab records); n-hexylbenzene 1.20x (1 system). D93 margin floored to 0 and must not be used alone. |
| Density | Dev MAE 11.81 kg/m3, bias -11.1 (locked -7.9, all 57 fuels -10.4). Margin 22.25 / 23.61 kg/m3. | Cause of the systematic low bias is unknown. |
| Heat of combustion | Group additivity, margins about 0.19-0.22 MJ/kg. | |
| Freezing | Dev-only offset +18.04 C; locked MAE 3.61 vs constant-mean 3.32. | No skill over a constant mean. SLE: ideal SvL defined only on 8 n-alkane systems / 103 points (MAE 2.02 K); Coutinho margin 38.27 K over 44 groups; the two freezing models disagree on 2 jet designs. |
| Smoke point | Conservative CV 1.29 vs 1.47 mm (grouped by inferred tank code). | No skill on the broad set; D1655 limits come from secondary sources only; BIS spec values UNVERIFIED. No valid margin. |
| Octane (bRON) | Linear blending is worse than leaving the base fuel unchanged. Pooled-BOV beats the constant mean in both views; GBM/GP beat it on stored folds only (leave-one-blendstock-out CIs include 0). | 4 locked CAS groups. No valid margin (rank 18 of 18). |
| DVPE | Evaluated on 51 BOB+oxygenate rows. | Units unconfirmed (psi inferred); not used for claims. |
| Design | 0 of 13 designs certified (5 fail at point estimate, 8 undetermined); v3 recipes re-scored 0 of 6. Independent verifier: 17 agree / 8 skipped / 0 disagree. | Certification is jet-only; diesel/gasoline passes are `pass_jet_derived_margin_only` and not counted. Viscosity ridge unavailable in robustness rows because the envelope guard refuses pure n-decane-like extrapolations. |
| API | No lazy retrain; yaml `applied_margin` intervals; only `fuel_type="jet"` can be in domain. | `fuel_type` is trusted from the caller (composition not checked). Viscosity certified only -40..20 C. Residual GBM opt-in, no margin. |

## 7. Known weak spots and open items (the reviewer should pressure these)

1. Reconstructed CRC647 compositions: how much of the "real fuel" error is reconstruction rather than physics? No measured GC x GC compositions for these fuels.
2. Thin locked sets (12 fuels, 4 CAS groups): confidence intervals are wide; negative or null results may be power problems, not model problems.
3. Flash margin hinges on one compound (mesitylene). Is the Antoine/LFL input for mesitylene (and n-hexylbenzene) correct? Conformal margin from 31 systems is fragile. **Update 2026-11-10:** unchanged by the review fixes (flash inputs are out of scope), but on the r2 reconstruction the raw ideal flash bias flips from -9.44 to +4.22 C, to the unsafe over-prediction side (CLM-38).
4. Density bias of about -11 kg/m3 with unknown cause (possible reconstruction, reference temperature, or pure-density inputs). **Update 2026-11-10:** on the r2 reconstruction the dev bias is -4.69 kg/m3 (was -11.11), still outside +/-3; about 4.7 of the bias remains and is not attributed to the old reconstruction (CLM-38). Known-composition lab density margins: over 4.587 / under 3.250 kg/m3 (CLM-39).
5. Freezing and smoke have no skill over trivial baselines, and no valid conformal margin. Is the problem the data, the physics, or the validation?
6. The Gaussian-process viscosity model under-covers; no pairwise test among bilinear/ridge/GP.
7. Design engine certifies nothing; margins were validated on jet fuels only, so extrapolation to other fuel types is unsupported. Design vs verifier agreement (17/8/0) rests on the same physics libraries. **Update 2026-11-10:** D2 re-run with lab margins: categories 7/0/0/6 (v4 5/0/0/8), still 0 of 13 certified; jet -20 C viscosity constraints are uncertified because the cold lab bin is `inf`; E10 designs are now ethanol-free and fail density (CLM-40).
8. DVPE units unconfirmed; D1655 smoke limits secondary-source only; BIS spec limits UNVERIFIED; no citable IS 2796 DVPE minimum (ruling R4 not applied for BIS E10).
9. Latent ThermoML bug: default-mode `phase_of()` reads only `ePhase`, so VLE vapour compositions may exist in the existing `thermoml_hc.csv`. Recorded as CLM-37, not fixed, no re-extraction run. v4 density tier-1 inputs read from the same extraction (impact unverified, probably small).
10. `reproduce_all.py` does not re-run encoder training, post-hoc encoder fixes, the Antoine/Cox gate comparison, or ThermoML re-extraction. `--verify-only` rewrites `SHA256SUMS_v4.txt` rather than diffing it (cannot detect drift).
11. Test suite: 308 passed, 42 failed; all 42 failures are assertions on v2/v3 deliverable files that no longer exist (pre-existing, but they make the suite look red).
12. Deferred minor review items from Tasks 0-18 are listed in `.superpowers/sdd/2026-11-09-v4-integrity-and-upgrade/progress.md` (local, git-ignored).
13. Uncommitted user changes exist (`splits/locked.json`, `src/fuelblend/data/model_table.py`, `results/v3/R1/viscosity_leaderboard.csv`, the Cox-gate script). Their effect on results is unassessed; v4 code does not read `splits/locked.json`.
14. No validation on realistic 3+ component hydrocarbon blends with measured properties beyond reconstructed CRC fuels; pour point has little data and weak physics.

## 8. Questions for the reviewer

- Are the claims in `docs/claims_register.md` (especially CLM-12, 13, 16, 18, 19, 25, 26, 28, 30-37) supported by the cited results, with no overstatement?
- Is the split design (dev vs locked, grouped CV, selection inside folds) leak-free? Look for any use of locked fuels in offsets, margins or hyperparameters.
- Is split-conformal calibration appropriate with 31-45 groups, and is the "larger of lab and D93 margin" rule sound?
- Is it reasonable to report null results (freezing, smoke, octane) as is, or is there a better baseline or model family worth trying?
- What is the highest-value next data to acquire (measured compositions, more locked fuels, diesel, RP-1 blends, real 3+ component blends)?
- What would make the design engine produce certifiable designs without weakening certification rules?
- Any cross-document contradictions among `FULL_REPORT_v4.md`, `claims_register.md`, `model_decision.md`, `docs/PROJECT_HANDOFF.md`, `AGENTS.md`?

## 9. Where to look first

- `docs/FULL_REPORT_v4.md` (per-property results, all numbers cited)
- `docs/claims_register.md` (claim status vs evidence)
- `docs/model_decision.md` (section 5/6: v4 routing and decisions)
- `results/v4/*/README.md` (per-task detail and limits)
- `configs/uncertainty_v4.yaml` (generated margins; `applied_margin`, `inf`)
- `scripts/v4/reproduce_all.py` and `SHA256SUMS_v4.txt` (reproducibility)
- `AGENTS.md` section 12 (project handoff) and `docs/PROJECT_HANDOFF.md`
