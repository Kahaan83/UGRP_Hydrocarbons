# Task 17 reading: HANNA vs UNIFAC vs ideal gamma in the v5 flash route (dev; exploratory)

Rows: 92 Leenhouts hydrocarbon binary flash rows, 16 systems (locked systems excluded), 0 dropped. Numbers: `summary.json`, `predictions.csv`. Script: `scripts/v7/flash/hanna_vs_unifac.py`.

- MAE (degC): ideal 3.20, UNIFAC 3.81, HANNA 3.66.
- Primary, HANNA vs ideal: MAE gain -0.48 (90 % CI -0.75 to -0.24): HANNA is significantly worse than the production ideal route. Gate fails.
- Secondary, HANNA vs UNIFAC: MAE gain -0.21 (CI -0.54 to 0.10): no difference. Gate fails.
- Over-prediction (the dangerous side) is negligible in all arms (mean positive part 0.09 ideal, 0.02 UNIFAC, 0.02 HANNA); the errors are under-prediction. The over-prediction criterion is met but the MAE criterion is not.
- Caveats: 16 systems only; HANNA training data may include some of these binaries; binary only (no use for multicomponent fuels).

Status: no skill shown for HANNA in the flash route; production stays ideal. Recommendation: do not spend the fresh lock (`results/v7/lock/`) on this comparison, since the dev reading gives no reason to expect a gain; the lock stays unspent unless the user directs otherwise.
