# v7 summary (2026-11-22)

Dev-only and exploratory; every number is read from a committed file under `results/v7/`. No v5 or v6 verdict is reopened and no model is adopted. Protocol: `docs/v7_protocol.md`.

| Task | Result | Source |
|---|---|---|
| 10 Track B D-MPNN retrain | beats v1 on tm_K only (+3.20 K, Holm p 0.015); beats trees on dcn only; trees beat the GNN on flash_C, dhfus, logvisc | `PHASE3_READINGS.md`, `track_b/retrain/READING.md` |
| 11 MolPool mixture and flash oxygenates | none beats its best comparator; flash oxy +0.47 vs desc_gbm, Holm p 0.162 | `mixture/READING.md`, `flash_oxy/READING.md` |
| 12 set models | no cell passes the gate | `set_models/READING.md` |
| 13 tree retrains | only ysi_lgbm beats the refit original (+35.9, CI 11.9 to 64.2) | `PHASE3_READINGS.md` |
| 14 freezing limiting species | worse than constant | `freezing/baseline_d/READING.md` |
| 15 freezing GP and novelty flag | fails the skill rule; novelty flag rests on 4 to 5 fuels | `freezing/gp/READING.md` |
| 16 freezing joint MolPool | LOGO MAE 3.00 vs constant 6.64 and ridge 5.83; joint vs fuel-only CI includes 0 | `freezing/molpool_joint/READING.md` |
| 19 conformal | global margin under-covers some families (YSI poly_arom 0.28); Mondrian restores 0.92 to 0.96 at large width | `uq/READING.md` |
| 20 fresh lock | 61 unseen systems (658 rows) locked | `lock/LOCK.md`, `lock/overlap.json` |
| 17 HANNA vs UNIFAC vs ideal gamma (flash, 92 dev rows, 16 systems) | HANNA worse than ideal (MAE gain -0.48, CI -0.75 to -0.24); no difference from UNIFAC; fresh lock left unspent | `flash/hanna/READING.md` |
| 18 Antoine and Andrade heads | beat pooled curve, not the descriptor ridge on log P (-0.34); flash recompute over-predicts 16 degC on average | `track_b/param_heads/READING.md` |
