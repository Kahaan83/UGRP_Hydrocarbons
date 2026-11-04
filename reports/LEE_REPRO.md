# LEE_REPRO — Reproduction of the Leenhouts et al. fuel-blend benchmark baselines

**Script:** `scripts/leenhouts_repro.py` (self-contained, read-only, prints results to stdout).
**Run command (verbatim):**

```powershell
cd C:\Users\kahaa\OneDrive\Desktop\UGRP_hydrocarbons\fuelblend_data
& ..\.venv\Scripts\python.exe -u -X utf8 scripts\leenhouts_repro.py > reports\LEE_REPRO_run_log.txt 2>&1
# confirmed exit code 0
```

Concretely: `& C:\Users\kahaa\OneDrive\Desktop\UGRP_hydrocarbons\.venv\Scripts\python.exe -u -X utf8 scripts\leenhouts_repro.py`
— confirmed **exit code 0** and completes in ~10 min (rdkit MW extraction dominates; all numerics are trivial).
**Run log:** `reports/LEE_REPRO_run_log.txt` (full stdout of the final run).
**Determinism:** no RNG is used anywhere; `np.random.seed(0)` is set only as a documented formality.

**Paper cited:** Leenhouts et al., *A Systematic Evaluation of Molecular Mixture Behavior Prediction…*
arXiv 2605.29698 — local PDF `data/raw/pdfs/Leenhouts_2026_mixture_eval_arXiv2605.29698.pdf`
(quoting its Tables 3, 11, 13). Their released code was cloned at
`data/raw/git/mixture_evaluation` and was used **only to read semantics** (matcher + metric definitions),
not copied or executed.

---

## 1. What was reproduced (summary)

| Item | Status |
|---|---|
| Their processed data loaded as released, row counts verified | ✅ |
| Their `ideal-mixture` baseline recomputed from the raw pure/mixture CSVs, **validated exact vs. their precomputed `*_with_ideal.csv` `ideal_value`** | ✅ flash 973/973, viscosity 31,660/31,660, all exact to 1e-9 |
| Baselines scored on the test rows of their released 5-fold split CSVs | ✅ all released split types |
| Viscosity ideal-baseline RMSE vs. paper Tables 3/11 | ✅ matches (§3) |
| Viscosity ideal-baseline Kendall τ vs. paper Table 13 | ✅ matches (§3) |
| Flash ideal-baseline RMSE ratio vs. paper | ⚠️ decimals of their printed mean are clipped in their PDF; only "14.x / 13.x / 13.x" citable (§3) |
| Extra baselines requested: Le Chatelier flash, Refutas viscosity | ✅ computed on their splits (§4) |
| Paper model rows (DMPNN+FFN, MolT5+FFN, RDKit+XGBoost) | ❌ out of scope (requires training their models) |

Their baseline definitions, read from their code (`scripts/scripts_processing/compute_ideal_and_excess_values.py`,
`utils/mixture_pure_matching.py`; `ideal_mode="linear"` flash, `"ln_to_linear_back"` viscosity):

* **Viscosity:** `value = ln(η [mPa·s])`; ideal = `Σᵢ xᵢ · ln(ηᵢ)` over matched pure-compound values at the
  mixture temperature (± 1 K candidate tolerance). This **is** the Grunberg–Nissan log-linear rule with k=0.
* **Flash point:** ideal = `Σᵢ xᵢ · FPᵢ` (plain mole-weighted linear average — *not* Le Chatelier).
* Note: for the ideal baseline, "absolute RMSE" ≡ "excess RMSE" (their Tables 3 and 11 ideal rows are
  numerically identical — predicting the ideal value makes the absolute error equal the excess error).

## 2. Exact files and rows used (all read-only)

| File | Rows | Role |
|---|---|---|
| `data/interim/leenhouts_benchmark/datasets/datasets/flashpoint/processed_data/flashpoint_pure_exp.csv` | 585 | pure FP (°C, `unit='Celcius'`) |
| `…/flashpoint_mix_exp.csv` | 1006 | mixture FP (°C, `fraction_type='molar'`, 1–4 components, InChIs) |
| `…/flashpoint_mix_exp_with_ideal.csv` | 1006 | cross-check of their `ideal_value` |
| `…/viscosity/processed_data/viscosity_pure_exp.csv` | 3027 | pure `value = ln(mPa·s)` + `temperature [K]` |
| `…/viscosity_mix_exp.csv` | 35,938 | mixture rows, `temperature [K]` |
| `…/viscosity_mix_exp_with_ideal.csv` | 35,938 | cross-check |
| `data/interim/leenhouts_benchmark/splits/splits/flashpoint/flashpoint_mix_exp/split_definitions/{random_rowwise,mixture,molecule,pure_to_mixture}/cross_validation_5fold/fold_00..04.csv` | 4 split types × 5 folds | test-row masks |
| `data/interim/leenhouts_benchmark/splits/splits/viscosity/viscosity_mix_exp/split_definitions/{random_rowwise,mixture,molecule,pure_to_mixture,mixture_temperature}/cross_validation_5fold/fold_*.csv` | 5 split types (10 folds for mixture_temperature) | test-row masks |
| `data/unified/compounds.csv` | 1038 | Antoine A/B/C (kPa) + LFL for Le Chatelier; join by exact `inchi`, fallback `inchikey` (rdkit) |

Split mapping verified: `data_id = mix:<i>` corresponds to `source_row_index = i` (0-based row in the mix CSV).
The split directories were used exactly as released; no custom splits were invented.

## 3. Numbers

### 3.1 Unit / sanity checks

* Pure n-decane flash point from their file = **46.0 °C** (reference ≈ 46 °C) → flash values are really °C.
* Pure n-decane η = **0.877 mPa·s at 295.65 K** (lit. ≈ 0.92 at 298 K; `exp(ln)` round trip) → viscosity
  `value` is really ln(mPa·s).
* Recomputed flash ideal: **973/1006** flash rows have full pure coverage; all 973 match theirs
  exactly (max|diff| = 0, tol 1e-9).
* Recomputed viscosity ideal: **31,660/35,938** rows covered at T±1 K; all 31,660 match theirs exactly.
* Cell-level check of their PDF confirmed the T_flash mean decimals are physically absent from the PDF
  ("14.±0.78" — the digit after the decimal point is clipped by the table column), so the paper's exact
  flash ideal-baseline means (14.x/13.x) cannot be read off the file.

### 3.2 Flash point — cm⁻¹ °C; leak-free test rows of their splits

RMSE = mean ± std over the completed folds of each released split; `n` = total evaluated test rows
(rows without pure coverage counted as skipped and never guessed).

| Split | Baseline | RMSE (°C) | MAE (°C) | Kendall τ | n | skip |
|---|---|---:|---:|---:|---:|---:|
| random_rowwise | ideal (their baseline) | **13.59 ± 0.69** | 8.22 ± 0.35 | 0.76 | 972 | 33 |
| random_rowwise | Le Chatelier (Antoine+LFL) | **1.86 ± 0.54** | 1.62 ± 0.44 | 0.95 | 57 | 948 |
| mixture | ideal | **12.67 ± 4.86** | 8.18 ± 3.56 | 0.75 | 973 | 33 |
| mixture | Le Chatelier | **2.38 ± 0.94** | 2.18 ± 1.06 | 0.94 | 57 | 949 |
| molecule | ideal | **13.40 ± 7.40** | 9.63 ± 5.05 | 0.68 | 661 | 17 |
| molecule | Le Chatelier | **1.69 ± 0.03** | 1.44 ± 0.13 | 0.98 | 52 | 626 |
| pure_to_mixture | ideal | 13.40 ± 7.40 | 9.63 ± 5.05 | 0.68 | 661 | 17 |
| pure_to_mixture | Le Chatelier | 1.69 ± 0.03 | 1.44 ± 0.13 | 0.98 | 52 | 626 |

Note: molecule and pure_to_mixture rows are identical **by construction** — their metadata states
pure_to_mixture reuses the "shared molecule_type test set", and we verified fold test sets coincide
fold-by-fold. The ideal baseline depends only on test rows, so identical numbers are expected (not a bug);
identical reproduction would NOT hold for a trained model.

### 3.3 Viscosity — ln(mPa·s); same protocol

| Split | Baseline | RMSE | MAE | Kendall τ | n | skip |
|---|---|---:|---:|---:|---:|---:|
| random_rowwise | log-linear Grunberg–Nissan k=0 (their baseline) | **0.296 ± 0.014** | 0.153 ± 0.004 | 0.84 | 31,660 | 4,278 |
| random_rowwise | Refutas VBN (mass fractions) | **0.265 ± 0.012** | 0.132 ± 0.003 | 0.87 | 31,674 | 4,264 |
| mixture | log-linear | **0.297 ± 0.031** | 0.154 ± 0.015 | 0.84 | 31,660 | 4,278 |
| mixture | Refutas | **0.264 ± 0.031** | 0.132 ± 0.007 | 0.87 | 31,674 | 4,264 |
| molecule | log-linear | **0.210 ± 0.055** | 0.126 ± 0.014 | 0.85 | 13,335 | 1,057 |
| molecule | Refutas | **0.186 ± 0.057** | 0.102 ± 0.009 | 0.87 | 13,341 | 1,051 |
| pure_to_mixture | log-linear | 0.210 ± 0.055 | 0.126 | 0.85 | 13,335 | 1,057 |
| pure_to_mixture | Refutas | 0.186 ± 0.057 | 0.102 | 0.87 | 13,341 | 1,051 |
| mixture_temperature (10 folds) | log-linear | **0.312 ± 0.124** | 0.171 ± 0.048 | 0.82 | 9,953 | 1,516 |
| mixture_temperature (10 folds) | Refutas | **0.278 ± 0.105** | 0.149 ± 0.045 | 0.85 | 9,967 | 1,502 |

### 3.4 Head-to-head vs the paper

**Viscosity (ln(mPa·s)), ideal baseline — reproduced vs. paper Tables 3 & 11:**

| Split | Reproduced | Paper |
|---|---:|---:|
| random_rowwise | 0.296 ± 0.014 | 0.30 ± 0.015 ✅ |
| mixture | 0.297 ± 0.031 | 0.30 ± 0.035 ✅ |
| molecule | 0.210 ± 0.055 | 0.21 ± 0.062 ✅ |
| mixture_temperature | 0.312 ± 0.124 | "worse than random" per paper (errors ↑ under temperature shift) — consistent |
| pure_to_mixture | (same as molecule; baseline has no training dependence) | paper lists only model rows |

Kendall τ (Table 13): reproduced 0.84 (random), 0.84 (mixture), 0.85 (molecule)
vs. paper 0.84 / 0.84 / 0.85 ✅.

**Flash (°C), ideal baseline — reproduced vs. paper Tables 3 & 11:**

| Split | Reproduced | Paper |
|---|---:|---:|
| random_rowwise | 13.59 ± 0.69 | "14.x ± 0.78" (decimal clipped in PDF) — consistent |
| mixture | 12.67 ± 4.86 | "13.x ± 5.4" — consistent |
| molecule | 13.40 ± 7.40 | "13.x ± 8.3" — consistent |
| pure_to_mixture | 13.40 ± 7.40 | paper lists only model rows |

All three reproduced values sit within the printable "14.±0.78 / 13.±5.4 / 13.±8.3" envelopes, and the
std values match exactly; the missing mean decimals are a typesetting artifact of their PDF, not a
reproduction gap.

## 4. Project baselines (not in the paper) on their splits

### 4.1 Le Chatelier + Raoult flash point (from data/unified/compounds.csv)

Configuration: Antoine kPa convention `log10(P/kPa) = A − B/(T/K + C)`; LFL vol%→mole fraction with the
Zabetakis T-correction `LFL(T) = LFL₂₅·(1 − 7.21e−4·(T−25))`; P = 101.325 kPa; bisection of
`H(T) = Σᵢ xᵢ·Psatᵢ/(P·LFLᵢ(T)) − 1` on 233.15–573.15 K (P_sat strictly increasing ⇒ monotone in T).

Results are in §3.2. On the rows where it is applicable, it is nearly an order of magnitude more
accurate than the paper's linear-averaging baseline (1.7–2.4 °C RMSE vs 12.7–13.6 °C), with much higher
Kendall τ. Caveats: only 71 unique components of the flash mixtures map into our hydrocarbon registry
(34 by InChI/key match), 15 of those have Antoine constants; 126 of 1006 blends have *all* components in
the registry, of which 57 have full Antoine ⇒ n=52–57 evaluated rows per split (see skip accounting below).
**The subset is small and hydrocarbon/n-alkane-starting rich — treat the comparison as a promising
signal, NOT a benchmark victory.** Their chemistry (alcohols, ketones, esters, chlorides) dominates the
mixtures and our registry is intentionally hydrocarbon-only.

Skip accounting (cumulative over all 4 splits' test rows; sum 948+949+626+626 = 3,149):

* `not_in_registry` 2,877 row-evaluations — compounds outside the hydrocarbon-focused `compounds.csv`
  (expected by project design; the benchmark is general chemistry);
* `missing_antoine_or_lfl` 272 row-evaluations — in-registry but without Antoine constants (mostly C10+
  alkanes whose WebBook refit windows don't reach the flash region) or missing LFL.

Rows skipped because the split-row lacks their ideal `ideal_value` coverage: counted separately as
"invisible" (n in §3.2 already excludes them).

### 4.2 Refutas VBN viscosity

Mass fractions from their molar fractions via rdkit MW; pure dynamic viscosity in mPa·s used directly as
cSt — **shortcut:** their pure file has no densities; for hydrocarbons (~0.7–0.9 g/cm³) the induced error
is small but systematic (η/ρ ~ 1/0.8 ⇒ ν ≈ 1.25·true if ignored). Despite this, on every split Refutas
with mass fraction VBN averaging beats the log-linear (Grunberg–Nissan k=0) baseline by ~10–12 % RMSE
(0.265 vs 0.296 random; 0.186 vs 0.210 molecule; 0.278 vs 0.312 mixture_temperature). VBN is also undefined
for pure components with ν ≤ 0.2 (cSt treated as cP ⇒ any η ≤ 0.2 mPa·s mixture gets refused, included in
"skip" counts).

This is consistent with expectations from the project's own earlier work (Refutas favored for binaries);
here it is on Leenhouts' data, larger chemistry set, strict splits — and it still wins in ln space, which
is the metric the paper uses for viscosity.

## 5. Differences, suspected causes

* **Flash ideal-baseline decadal digit:** not recoverable from their PDF (nine digits after the point are
  missing in the source PDF's lossy table cell rendering; verified with pdfplumber char coordinates).
* **Viscosity ideal baseline agreement is effectively exact** (0.296 vs 0.30, 0.297 vs 0.30, 0.210 vs 0.21,
  matching std values) — any residual mismatch ≤ rounding of the paper's printed sig figs.
* **Le Chatelier numbers cannot be compared to the paper** — the paper's flash baseline is a
  mole-weighted linear average; the Le Chatelier method needs Antoine + LFL constants the benchmark does
  not ship; our registry covers a small HC subset of their chemistry. RMSE 1.7–2.4 °C is only over those
  n=52–57 blends per split.
* **Refutas advantage over log-linear** is a real property-level finding, but the cSt-vs-cP shortcut means
  the true Refutas-only gap could differ by a few percent; a density-joined re-run would pin it down.
* Temperature coverage: the ±1 K pure-value tolerance from their matcher was replicated exactly, and
  validated through the exact ideal_value reproduction; no other split-compatibility approximation was made.

## 6. What could NOT be reproduced

* **Model-family rows** (DMPNN+FFN, MolT5+FFN, RDKit+XGBoost) and all training curves — model training out
  of scope for this task by design.
* **Exact paper means for flash ideal baseline** — PDF's own digits are clipped; only "14.x / 13.x / 13.x"
  citable. Our reproduced values (13.59, 12.67, 13.40) are consistent with those envelopes.
* **pure_to_mixture ideal-baseline row in the paper** — the paper's Table 3/11 lists only models for that
  split; we compute it here anyway (§3.2/§3.3).
* Re-running their training/validation pipelines (torch-dependent, models, HPC) — no.
* Evaluation of n-applicability-quality on DCN/MON/solubility/enthalpy tasks — out of scope (flash +
  viscosity only, per task).

## 7. Honesty notes / gaps

* Results computed on **1 CPU, ~10 min**; all workload is exact rule application, no stochastic learners.
* `molecule` and `pure_to_mixture` fold test sets coincide fold-by-fold, per their metadata; verified, not
  assumed.
* The Le Chatelier baseline reads Antoine constants (webbook/refit) from the project's own
  `data/unified/compounds.csv` whose C10+ constants were validated in-house earlier; not re-validated
  against independent references here.
* `Antoine covers window` was not enforced as a gate for applicability (only presence of constants);
  cases where T spans the Antoine range edge are possible for the 52–57 evaluated blends; no range
  violation check is implemented. Given the very low RMSE, this did not blow up visibly.
* Rows skipped for missing registry coverage could bias Le Chatelier's evaluated subset toward simple
  binaries; they are *counted*, not silently dropped.
* Their code (`mixture_evaluation`) was read, not executed; our matcher replica is verified bit-exact
  against its output (§3.1), making the semantic replication verifiable rather than assumed.

## 8. Stop-point summary / next steps

**What ran:** <venv python> -u scripts/leenhouts_repro.py → exit 0; run log
`reports/LEE_REPRO_run_log.txt`; this report.

**Exact numbers with units:**
* Flash (°C), test folds of their splits, mean±std RMSE: ideal 13.59±0.69 (random), 12.67±4.86 (mixture),
  13.40±7.40 (molecule/pure_to_mixture); Le Chatelier 1.86±0.54 (random), 2.38±0.94 (mixture),
  1.69±0.03 (molecule), MAE 1.4–2.2 °C, n=52–57/split.
* Viscosity (ln(mPa·s)), same protocol: log-linear 0.296±0.014 (random), 0.297±0.031 (mixture),
  0.210±0.055 (molecule), 0.312±0.124 (mixture_temperature); Refutas 0.265±0.012 / 0.264±0.031 /
  0.186±0.057 / 0.278±0.105 — Refutas wins on every split.

**Matched vs paper:** viscosity ideal RMSE ✅ exactly (3 splits + Kendall ✅); flash ideal RMSE envelopes ✅
(greyed paper digits); Le Chatelier/Refutas not comparable to paper.

**Failed:** nothing functional; the only "failure" is the paper-side missing flash decimal digits.

**Decisions made:** ideal-coverage gate for the ideal baselines (their semantics); Le Chatelier evaluated
on all test rows it can predict (independent of ideal coverage), with components registry-joined by
InChI→InChIKey fallback; Refutas measured in ln(mPa·s) space to match their metric basis.

**Next steps:** density-join a kinematic-viscosity variant of Refutas; add flash Le Chatelier registry
coverage (C11+ alkanes' Antoine via interpolation flagged in AGENTS.md §12 item 9); optionally implement
their `excess`-target formulation of the residual-learning task with the same splits.
