# ANALYSIS PLAN — Fuel-Blend ML (draft for approval)

Status: DRAFT — needs explicit user approval before any locked-set access or model tuning.
Date: 2026-10-30 · split seed 20261030 · model seeds [11, 23, 37, 41, 59]

## 0. Hard rules (non-negotiable)
- Train/evaluate only from `measurements_dedup.csv` (via `data/processed/model_table.csv`); never raw rows, never synthetic labels.
- Never invent numbers; every metric comes from a run with a logged command + seed.
- Locked sets (`splits/locked.json`: 24 real fuel_ids + 12 multi system_ids under `splits/s1s2s3.json` exclusions) are accessed exactly once at the end; every access is appended to `reports/LOCKED_ACCESS.log`.
- Conflict-flagged rows (`conflict=True`) excluded from headline metrics; sensitivity run reported separately.
- Hydrocarbon-only (`is_hc`) is the default filter.

## 1. Evaluation modes
- **Mode A (known-pure)**: every blend component appears in the pure-compound training pool.
- **Mode B (novel-molecule)**: S2/S3 molecule holdouts — measure generalization to unseen molecules (S2: partial nan-features allowed; S3: full holdout).

## 2. Splits
- S1 (interpolation, 28,729), S2 (molecule out, 8,050), S3 (all-new molecules, 826) — generator `src/fuelblend/splits/split_generator.py`, manifests in `splits/s1s2s3.json`.
- Locked test: 20 % per source bucket, never tuned on. Real_fuel rows excluded from S1/S2 by construction.
- Secondary stress splits (built later from same generator): class-out, order-extrapolation (carbon number), T-extrapolation, source-grouped.
- All splits deterministic under seed 20261030; hashes logged in DATA_CARD.

## 3. Metrics (src/fuelblend/eval/metrics.py)
- Per property: MAE, RMSE, median APE, fraction within tolerance (viscosity/density: ±5 % rel; flash: ±5 °C abs; cold-flow: ±3 °C).
- Flash extra: **one-sided over-prediction rate** (predicted − true > 5 °C) — the safety-relevant error direction; separate display, never averaged away.
- Uncertainty: cluster bootstrap over `system_id`, 2,000 resamples, percentile 95 % CI.
- Success gate: model beats physics baseline by ≥ the measured noise floor (NOISE_FLOOR.md) with a non-overlapping CI, and conformal 90 % interval coverage within 2 pp of nominal.

## 4. Noise floor (reference values, reports/NOISE_FLOOR.md)
- viscosity_dyn_mPas: median spread 0.016, p90 0.12 (824 rows ≥3 obs)
- density_kg_m3: 0.73 / 5.44 (448)
- flash_point_C: 6 / 29.4 (289, conflict rate 24 %) — treat flash noise floor as per-method; ASTM statements to be inserted (T00).
- vapor_pressure_kPa: 0.41 / 4.14 (380)
- viscosity_kin_cSt, freezing, mon, dcn, smoke, net heat: n < 3 → no floor estimated; report as "insufficient replication".

## 5. Baselines (Phase 3, before any learned model)
Refutas (nu log-blend) → Grunberg-Nissan → Kay mixing → ideal-volume density → Le Chatelier/Raoult flash → PPBI pour → volume-linear octane → mass-additive GCV → TSI smoke. Each evaluated identically on S1/S2/S3, per property × mode × split.

## 6. Models (after baselines are locked in)
- Pure-compound models per property (features from molecular descriptors/class/n_carbon).
- Mixture residual model: physics baseline + learned correction (GBM first, then permutation-invariant DeepSets encoder with per-property heads); flash model trained with asymmetric loss on over-prediction, weights {1, 3, 5, 10} ablation.
- Conformal (Mondrian by property × kind) for calibrated intervals.
- Gated properties → only what GATES.md classifies as G2 gets a learned model; G3 = baseline-only report.

## 7. Ablations (planned)
- quality_weight on/off · source-tier filter (A-only vs A+B vs A+B+C) · mode A vs B · residual on/off · asymmetric-loss weight sweep · class-out robustness.

## 8. Open questions for approval
1. Approve locks as generated (locked.json hashes) → sign-off unlocks WP3 repro work.
2. Provide/approve ASTM precision quotes for T00 SOURCE_CHECK.
3. Confirm G2 threshold set (viscosity_dyn, density, smoke learnable now; others baseline-only).
4. OneDrive → local move + `git init` still pending; requested before heavy training runs (file-loss risk observed).
5. Timing: with NOISE_FLOOR gaps for most properties, retry counts below floor are sample-limited — accept reporting "no floor" rather than imputing?

---

## 9. Round 2 Data Addendum & Composition Rule R1 (Approved 2026-10-30)

### 9.1 Rule R1 and Composition Hierarchy
**Rule R1**: Use a fuel blend only if its composition is available; otherwise use pure compounds.
Composition tiers:
- **L3 (Molecule-resolved)**: Gravimetric/molar fractions with identified molecular structures (SMILES/InChIKey). Eligible for training and testing across all model stages. Primary benchmark set.
- **L2 (Group-resolved)**: Class by carbon number (e.g. GC×GC 106-group distributions or n-alkane distributions with lumped remainder). **Eligible with restrictions**: Used strictly for the real-fuel bridge stage (training & evaluation of real-fuel bridge/surrogates). Reported separately from L3 in every table; never mixed into an aggregate metric. NIST RP-1 is test-only.
- **L1 (Class totals only)**: PIONA totals, ASTM D2789 mass specs, bulk NMR groups without carbon-number resolution. **Excluded** from model training and test sets. Routed to `data/processed/excluded_by_R1.csv`.
- **L0 (No composition)**: Bulk distillation curves (D86), bulk properties without chemical breakdown. **Excluded** from model training and test sets. Routed to `data/processed/excluded_by_R1.csv`.

### 9.2 Impact on Work Package 7 (Real-Fuel Bridge)
The previous bulk class-fraction / distillation stacking model (WP7) trained on L1/L0 aviation fuels (AJFTD, CRC bulk without GC×GC) is discontinued as a blend-learning model under R1. Real-fuel bridge modeling is re-anchored exclusively on L2 fuels with quantitative group breakdowns (Vozka 2019, Caltech FM99-4, and AJFTD GC×GC fuels) plus L3 surrogates. L1/L0 records are preserved in `excluded_by_R1.csv` for provenance and correlation checks.

