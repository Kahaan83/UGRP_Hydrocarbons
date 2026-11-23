# v7 fresh flash lock: not scorable

`scripts/v7/flash/score_fresh_lock.py` was run once (user authorised spending the lock, 2026-11-23). The v5 flash route (`flash_unifac._prepare`) reports missing Antoine or LFL constants for every one of the 658 rows (61 systems), so no ideal, UNIFAC or HANNA flash point was produced and no prediction was compared with a measurement. Additionally HANNA is binary-only, and 215 of the 658 rows are ternary or quaternary.

Outcome: primary comparison (HANNA vs UNIFAC) not evaluated, status: no skill shown. The gate was neither passed nor failed. Scoring would need Antoine and LFL constants for the lock compounds; that is outside v7 scope and is not done. `predictions.csv` and `summary.json` here are empty or NaN by construction.
