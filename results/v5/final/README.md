# v5 final numbers pack: tables and figures

Built by `scripts/v5/final/build_final_numbers_v5.py` from files committed at `a3c5967` (`git show`, never the working tree).
Master document: `docs/final_numbers_v5.md`. Nothing here is recomputed; values are selected, joined and plotted.

| file | content | sources |
|---|---|---|
| `table_A_headline.csv` | headline designated rows (current defaults = re-score, plus frozen BASE side by side; join on the full key incl. cells_file/status/evidence_path); plus the UNIFAC SvL designated freezing rows on T2 pooled / pooled_locked (not scoreboard cells; `scoreboard_variant` = `(not a scoreboard cell)`) | `results/v5/BASE/delta_adopt_rescore/rerun/scoreboard.csv`, `results/v5/BASE/scoreboard.csv`, `results/v5/BASE/delta_freeze_t2_unifac/summary.csv` |
| `table_A_freezing_unifac_lab.csv` | UNIFAC vs ideal SvL on lab liquidus S_F | `results/v5/BASE/delta_freeze_adopt/summary.csv` |
| `table_A2_t2_per_source.csv` | T2 per-source designated rows | same as table A |
| `table_B_drift_by_cause.csv` | changed scoreboard cells by cause | `results/v5/BASE/delta_adopt_rescore/drift.csv` |
| `table_B2_drift_other_files.csv` | non-cell drift rows (rows_* summaries) by file and cause | `results/v5/BASE/delta_adopt_rescore/drift.csv` |
| `table_C_margins.csv` | applied margins + spec limits | `configs/uncertainty_v5.yaml`, `configs/specs_astm.yaml`, `configs/specs_bis.yaml`, `results/v5/BASE/delta_ml_viscosity/readings.json`, `results/v5/BASE/delta_ml_density/readings.json`, `results/v5/BASE/delta_adopt_rescore/rerun/scoreboard.csv` |
| `table_D_ml_gates.csv` | ML gate readings (R-MLV-1..5, R-MLD-1..3) | `results/v5/BASE/delta_ml_viscosity/readings.json`, `results/v5/BASE/delta_ml_viscosity/readings_R-MLV-5.json`, `results/v5/BASE/delta_ml_density/readings.json` |
| `table_D2_ml_margins.csv` | ML margins and locked coverage | `results/v5/BASE/delta_ml_viscosity/readings.json` (R-MLV-6), `results/v5/BASE/delta_ml_density/readings.json` (R-MLD-4) |
| `table_E_designs.csv` | 16 designs | `results/v5/BASE/delta_adopt_designs/results_table.csv` |
| `table_F_negative_results.csv` | negative results | `docs/v5_negative_results.md`, `docs/model_decision.md`, `results/v5/BASE/delta_ml_ysi/readings.json`, `results/v5/BASE/delta_ml_density/readings.json`, `results/v5/BASE/delta_adopt_rescore/rerun/scoreboard.csv`, `results/v5/BASE/delta_freeze_t2_unifac/summary.csv`, `results/v5/BASE/delta_freeze_t2_unifac/rows.csv` |
| `table_G_limitations.csv` | limitations | `configs/uncertainty_v5.yaml`, `results/v5/BASE/delta_freeze_adopt/per_system_S_F.csv`, `results/v5/BASE/delta_freeze_adopt/rows.csv`, `results/v5/BASE/delta_freeze_adopt/README.md`, `results/v5/BASE/delta_adopt_designs/results_table.csv`, `results/v5/BASE/delta_adopt_designs/designs.csv` |
| `table_H_commits.csv` | commits for reproduction | git log |

Headline selection rule: designated rows of the headline properties, no temperature split, T2 sets starting with `pooled`
(or `measured_pooled`), CRC647_r3 labels dev and fold_b, T1_lab label dev, T3_ajftd label dev excluding the trivial models.

## Figures (`figures/`, 200 dpi, Okabe-Ito palette)

| figure | content | source |
|---|---|---|
| `a_flash_parity_t2.png` | flash parity, T2 real fuels (70 chosen rows), designated `pred_designated`; one-sided over-prediction band of the applied margin (2.21 C). The script asserts the rows reproduce the scoreboard per source. | `results/v5/BASE/delta_adopt_rescore/rerun/rows_t2_flash.csv` |
| `b_viscosity_refutas_vs_ml.png` | lab |ln resid| by T bin, Refutas vs ML; real-fuel R-MLV-5 MAE per source x T | `results/v5/BASE/delta_ml_viscosity/rows_scored.csv`, `results/v5/BASE/delta_ml_viscosity/realfuel_transfer_R-MLV-5.csv` |
| `c_density_ideal_vs_ml.png` | lab |error| ideal vs ML by binary/ternary and dev/locked | `results/v5/BASE/delta_ml_density/rows_scored.csv` |
| `d_margins_vs_spec.png` | margins per property with the Jet A / Jet A-1 limits in the panel title | table C |
| `e_ml_before_after.png` | ML gate MAEs, baseline vs ML, with PASS/FAIL | table D |
| `f_design_headroom.png` | flash and freezing headroom of the 10 designs with declared margins | table E |
| `g_rescore_drift_by_cause.png` | changed scoreboard cells by cause | table B |

Deviations from the brief:
* (a) the brief asked for a +-margin band; the flash margin is one-sided (direction `over`), so a one-sided band is drawn.
* (b), (c) the brief asked for parity plots; `rows_scored.csv` of both ML deltas stores residuals only (no measured or predicted values), so residual distributions are plotted instead. No values were rebuilt by joining other files.
* (f) drawn as headroom bars (bound vs limit) for the 10 designs that have declared margins, instead of a heat map.
