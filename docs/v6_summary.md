# v6 summary

Status: 2026-11-21. v5 BASE stays frozen (`results/v6` and `scripts/v6` are new; the BASE manifest hashes are unchanged and
`scripts/v5/base/reproduce_base.py` passes). Every number below is read from a committed file named in the last column. All results are
model-level only: no lab measurements exist, nothing here is exposed in `predict_fuel`, and nothing is certified. Gates were declared and
locked before fitting; failures are kept in `docs/v6_negative_results.md`. Plan: `docs/v6_plan_2026-11-21.md`; literature:
`docs/literature_review_v6.md`.

## 1. Outcome at a glance
| Track | Question | Result | Status |
|---|---|---|---|
| A freezing point | Can anything beat the constant on real fuels? | Holdout (12 fuels, scored once): best model MAE 9.30 vs constant 6.86; ridge 6.06 (CI [-1.90, 4.07]) | **no skill shown** |
| B pure-compound D-MPNN | Does a multi-task GNN beat constants and existing estimators? | Beats the constant on all 10 targets (CV); Tm loses to the v5 estimator; **YSI holdout passes** | positive for YSI, weak elsewhere |
| C mixture ML | Does a pooled GNN head beat physics and a descriptor GBM on mixtures? | 1 of 5 heads passed dev and locked, **did not replicate** on 51 fresh rows | **no skill shown** |
| Flash, oxygenate arm | Does MolPool help flash point of oxygenate blends? | Dev gain CI [0.16, 1.04], locked [-0.47, 1.38] | not established |

## 2. Data work
- **External ingest** (`scripts/external_ingest/`): KAUST Table 4/5, FACE surrogate and AI4fuel fixes with tests; superseded Boehm, CRC FACE-1
  and ternary-flash staging removed (duplicates of production). Known defect: AI4fuel heat-of-combustion values span roughly -22,100 to
  -199 under an "MJ/kg" label (sign or unit error), so LHV was not used.
- **AJFTD v6 tables** (`data/real_fuels/v6/`, `scripts/v6/data/build_ajftd_v6.py`, 12 tests): fuel groups (`dup_group`, used as the CV unit),
  long and wide property tables, GCxGC composition. Audit entry: `docs/data_audit_2026-11-11.md` section 7.

## 3. Track A: freezing point (`results/v6/freezing/`)
- Data: 64 groups (v5 T2 entities plus AJFTD 20, 21, 22, 64, 154, 155). Dev 52 groups (leave-one-group-out), holdout 12 groups, locked
  before models (`holdout_lock.json`). Reproduction gate: baseline (b) reproduces v5 (MAE 8.4246 on the 50 headline rows).
- Baselines: constant, v5 UNIFAC-SvL, mass-weighted Tm average, ridge on coarse composition, single-species control curve, isomer-averaged
  tiles; then a monotone LightGBM residual.
- Dev LOGO: the residual reached MAE 3.31 vs constant 6.64 (CI of gain [1.56, 5.25]). **On the locked holdout it did not replicate**
  (9.30 vs 6.86, CI [-6.31, 1.31]). UNIFAC-SvL holdout MAE 13.48; the fixed Tm average fails (23.99, bias +24 C). Reading in
  `docs/v6_negative_results.md`: the dev gain was driven by the AV-33-22 fuels (65 % of dev), and the 90 % margin covered 8 of 12 holdout
  fuels (0.67), so it is not adopted. Live freezing margin stays 13.84 C.
- Post-hoc step 5 (dev only, holdout untouched): filling the 1,042-compound registry's missing isomer Tm and dHfus with Track B or
  fusion_estimate values changes 31 tiles (4 before). Dev MAE 7.88 (D-MPNN fill) and 8.09 (fusion Tm fill) vs v5 8.83 vs constant 6.64. The
  entire gain is fuel 64; on the other 51 fuels the result equals v5 (7.31). No skill shown.

## 4. Track B: pure-compound D-MPNN (`results/v6/gnn/`)
- Model: directed MPNN in plain torch (`src/fuelblend/ml/dmpnn.py`), multi-task masked loss, 5-model ensemble, early stopping on an
  inner group-held-out split, two configurations chosen per outer fold. Pool: 1,208 compounds, 10 targets (Tm, dHfus, flash, density,
  log viscosity at 293 K, YSI, DCN, cetane, RON, MON); oxygenates carry octane labels only; Tb and LHV not trained. CV: GroupKFold(5) on
  generic Murcko scaffold.
- Grouped CV MAE (GNN vs pooled-median constant): flash 19.1 vs 55.2 C, density 21.9 vs 107.3 kg/m3, Tm 30.4 vs 71.8 K, YSI 69.4 vs 290.7,
  RON 9.0 vs 16.6, MON 8.2 vs 14.5, cetane 9.0 vs 18.6, DCN 5.6 vs 18.4, dHfus 12.8 vs 19.1 kJ/mol, log10 viscosity 0.225 vs 0.457. Every
  90 % bootstrap gain interval excludes zero. The constant is a weak bar: absolute errors are too large to replace physics or measurements
  for flash (mean one-sided over-prediction 9.4 C) or density.
- Existing estimators on the same rows: Tm 30.4 K vs v5 fusion_estimate 28.4 K (**no improvement**, so Tm is not used to fill Track A
  constants); YSI 69.4 vs v4 ridge 106.5.
- **YSI locked holdout** (`results/v6/gnn/ysi_holdout/`, declared and locked before fitting, scored once): 56 compounds, 17 scaffold
  groups. MAE GNN 59.6, ridge 133.6, constant 366.0; gain over ridge [32.6, 154.0], over constant [223.9, 435.2]. **Gate passes.** The 90 %
  absolute margin (190.1) covered 91 % of the holdout; it is a rough margin, not a clean conformal calibration (the CV folds it came from
  had seen the holdout compounds).
- Exploratory v2 (`results/v6/gnn_v2/`, wider config grid, longer patience): 0.5 to 8 % lower MAE per target (YSI 63.9, flash 18.1, density
  20.9, Tm 29.6). Chosen after seeing v1, no holdout, differences within plausible CV noise; v1 stays the declared result.
- A CUDA path exists (`FUELBLEND_DEVICE=cuda`, CUDA torch installed outside the venv); about 2x faster in a benchmark; the reported v1 and v2
  CV numbers are CPU, the YSI holdout fit ran on the GPU.

## 5. Track C: mixture ML (`results/v6/mixture/`)
- Model: frozen Track B encoder, per-component embeddings pooled by mole fraction (MolPool), bilinear pair term, residual over the physics
  baseline, config chosen by inner grouped CV, 5-seed ensemble (`src/fuelblend/ml/molpool.py`). Pretraining on all pure labels is fair
  (baselines also get measured pure values; the pool has no mixture labels).
- Heads: KAUST RON/MON (120 blends, 8 families), Leenhouts MON (241 rows, 52 systems; 51 low-coverage rows set aside), Boehm viscosity and
  density (23 mixtures at 5 temperatures). Co-Optima blends are not in the repo (excluded). Comparators: linear molar, corrected-pure (our
  definition), descriptor GBM; Boehm: physics, per-T bias, descriptor GBM. Gate: dev OOF CI above 0, locked CI above 0, n_comp >= 3 no worse.
- Result: KAUST and Boehm fail (Boehm viscosity and density significant dev losses). Leenhouts MON passed dev (gain CI [0.04, 1.33]) and
  locked ([0.07, 0.98]). **Replication on 51 unscored system-disjoint rows** (15 systems, pure MON gaps filled by the Track B head): MAE
  MolPool 5.28 vs descriptor GBM 4.97, gain CI [-1.24, 0.58]; not replicated. Post-hoc head tuning (v2) did not help. Status: no skill shown
  over the strongest comparator; the earlier pass is best read as chance.

## 6. Flash arm, oxygenate mixtures (`results/v6/flash_oxy/`)
Run in a separate session, summarised from its README. MolPool arm B vs descriptor GBM: dev (57 systems) MAE 3.44 vs 4.08, gain CI [0.16,
1.04]; locked (21 systems, once) 3.88 vs 4.48, CI [-0.47, 1.38]. Primary gate fails (point estimates favour MolPool on both sets, so
"not established", not "worse"). Hydrocarbon non-inferiority passes on 14 locked systems (upper bound 0.65) but fails on dev (3.29), so
it is not established; the oxygenate model is not to be used for hydrocarbon fuels. Le Chatelier is not evaluable on oxygenates.

## 7. What v6 adds, and what it does not
- Adds: a working multi-task GNN with a passing YSI holdout; a documented set of negative results for freezing point and mixture ML; AJFTD
  v6 tables; a reusable MolPool head and device option.
- Does not add: any property model that replaces a v5 route. Freezing point stays "no skill shown"; the live margin and designated routes
  are unchanged. Mixture ML does not beat a descriptor GBM on new systems.
- Limits: tiny mixture sets (5 KAUST dev families, 52 Leenhouts and 23 Boehm systems), one head of five passing without multiplicity
  correction before it failed replication, weak constant baselines for Track B, no lab data.

## 8. Reproduce and verify
- `python scripts/v5/base/reproduce_base.py` (BASE manifest and 48 tests; passed, 0 files changed).
- `pytest tests/test_v6_ajftd_tables.py tests/test_v6_freezing.py tests/test_v6_gnn.py tests/test_v6_mixture.py tests/external_ingest -q`
  (38 passed). Legacy v2/S4 test failures are accepted and out of scope (plan, "Known test failures").
- Track C needs `results/v6/mixture/encoder.pt` (gitignored); rebuild with `scripts/v6/mixture/pretrain_encoder.py` (about 6 min CPU).
- Declared-protocol order is visible in git history: each lock or README declaration commit precedes its scoring commit.
