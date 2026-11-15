# delta_freeze_unifac: UNIFAC liquid activity in Schroeder-van Laar on the scope that sets the 16.98 C freezing margin

Declared 2026-11-15 BEFORE any declared reading in this folder was computed. Base: v5 after commit 4fab6dd.
Nothing in `configs/`, the judge, `predict_fuel`, `tm_dhfus_verified.csv` or any spec flag is changed. Designating a different freezing model or
margin is a **user decision**.

## Orientation (no declared reading)

The live freezing margin `freezing_point_C.lab_liquidus.svl_verified` = 16.98 C is the worst of 9 binary systems (scope S_F,
`results/v5/BASE/delta_freeze_margin/rows.csv`). Per-system worst under-prediction (read before this declaration):
naphthalene + trans-decalin -16.98, naphthalene + cis-decalin -14.59, 1-methylnaphthalene + 2,6-dimethylnaphthalene -13.84,
bicyclohexyl + butylcyclohexane -9.08, then -4.69 ... +0.99. The three worst are aromatic + naphthene or substituted-naphthalene systems, where
ideal mixing is known to be a poor approximation. WS2 implemented `svl_unifac` (original UNIFAC liquid gamma, pure solids,
`fuelblend.eval.freeze_multicomp.svl_activity_C`) but scored it only on multicomponent class mixes, never on S_F.
Because I have seen the ideal residuals, this reading is labelled **exploratory / dev** (no locked set exists for SLE).

## Declared readings

* **R-FU-1 (S_F, same rows, same groups).** Model `svl_unifac`: for each S_F row, verified constants (`build_tmap_verified`, the governing
  `svl_dhfus_kJ_mol`), mole fractions, `make_unifac_gamma_fn(smiles)` (gamma = 1 fallback flagged), maximum over components of the T solving
  ln(x_i gamma_i(T, x)) = -dH_i/R (1/T - 1/Tm_i). Scored with the unchanged `delta_newdata_freeze_smoke.summarize` (direction `under`,
  alpha 0.10, conformal over systems). Reported: MAE, bias, margin, per-system worst, unsafe fraction, beside `svl_verified`.
  **Preference rule (for the user decision):** `svl_unifac` is preferred only if its margin is smaller than 16.98 C AND its MAE is not larger
  than `svl_verified`'s (3.01 C) AND no system has a UNIFAC fallback. Otherwise `svl_verified` stays.
* **R-FU-2 (wider check, descriptive).** The same comparison on `thermoml_all_verified_any_T_any_nC` (36 systems, any T and carbon number)
  and on the ASCENT mixtures and Vozka fuels if `svl_activity_C` can be applied with the existing row inputs. Reported, not a gate.
* **R-FU-3 (designs, counterfactual only).** If R-FU-1 prefers UNIFAC: freezing of the stored pure-pool optimum recipes (P0x-a, P0x-b,
  P0y, pool C) under `svl_unifac` with its own R-FU-1 margin, against the spec limit. Reported as "would pass / fail if adopted"; no design run.

## Results (appended 2026-11-15 after the run; nothing above was edited)

`scripts/v5/base/delta_freeze_unifac.py` -> `rows.csv`, `summary.csv`, `per_system_S_F.csv`, `readings.json`.

### Row count correction (found while running, affects a committed reading's MAE only)
The committed `delta_freeze_margin/rows.csv` has 77 S_F rows; `scope_rows` gives 67. The 10 extra rows are duplicates made by
`delta_freeze_margin.add_old` (a merge on system_id, T_K, sample against an old table that has repeated rows; tridecane + pentadecane 19 vs 11,
dodecane + tetradecane 18 vs 16). The **16.98 C margin is unchanged** (it is the per-system maximum), but the S_F MAE of the ideal model is
**3.33 C on the 67 distinct rows**, not 3.01 C. Both models below are scored on the 67 rows.

### R-FU-1 (S_F, 9 systems, 67 rows): UNIFAC preferred by the declared rule
| model | MAE (C) | bias (C) | margin alpha 0.10 (C) | unsafe fraction | UNIFAC fallback |
|---|---|---|---|---|---|
| svl_verified (live) | 3.33 | -1.21 | 16.98 | 0.45 | - |
| **svl_unifac** | **2.62** | +0.13 | **13.84** | 0.37 | none |

Per system (worst under-prediction, ideal -> UNIFAC): naphthalene + trans-decalin -16.98 -> **+0.06**; naphthalene + cis-decalin
-14.59 -> **+2.11**; 1-methylnaphthalene + 2,6-dimethylnaphthalene -13.84 -> -13.84 (aromatic-aromatic, gamma about 1: now the margin);
bicyclohexyl + butylcyclohexane -9.08 -> -8.82; the n-alkane and other dimethylnaphthalene systems change by < 0.7 C.

Reading: the three worst ideal systems were two aromatic + naphthene binaries that ideal mixing under-predicts by 15-17 C; UNIFAC's
activity coefficient (aromatic in a saturated ring solvent, gamma > 1) fixes both. The margin is now set by an aromatic-aromatic system
UNIFAC cannot help. **Exploratory / dev** (the ideal residuals were seen before this reading). Designating `svl_unifac` with 13.84 C is a
user decision.

### R-FU-2 (all verified ThermoML, any T and carbon number, 40 systems, 643 rows; descriptive)
svl_verified MAE 4.24, margin 17.81; svl_unifac MAE 3.92, margin 17.23. UNIFAC is a small improvement on the wider set too.
