# v7 fresh flash lock (user decision 2026-11-22: "whichever is recommended")

Spent: the 61 Leenhouts flash-mixture systems (658 rows) unseen by `measurements_dedup.csv` and every system-id lock (`lock.json`, ids listed there).
Primary comparison (single): HANNA gamma vs UNIFAC gamma in the v5 flash route (Task 17), per-system MAE gain, UNIFAC minus HANNA, group bootstrap 90 % CI over systems.
Gate: MAE gain CI above 0 and upper CI bound of the increase in mean over-prediction below 0.5 degC. Holm rule: only this one comparison is scored here.
Rules: nothing is fitted, tuned or inspected on these systems before scoring; they are scored exactly once after Task 17 is final and committed. If HANNA cannot be obtained, the lock stays unspent.
Caveat: 35 staging systems overlap `flash_oxy_lock` and 31 are in the dedup data; these are excluded from the lock.
