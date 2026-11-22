# Task 5 reading (dev rows only)

Predeclared rule: longer training is not the limit if median `best_epoch_` < 75 and val loss is flat or rising afterwards.
Median best epoch: v1-cap 36.5 (cap 150, patience 25), long-cap 136 (cap 1000, patience 100). 136 is above 75, so the rule is not met: the v1 patience of 25 stopped models early.
Outer MAE, v1-cap vs long-cap: lower with long-cap on 8 of 10 targets (ysi 92.7 vs 75.6, tm_K 31.8 vs 27.8, ron 9.17 vs 8.24, mon 8.33 vs 7.64, dhfus 13.3 vs 10.7, flash_C 19.7 vs 18.0, logvisc 0.239 vs 0.218, density 22.5 vs 22.2); equal on cetane (8.87 vs 8.87); slightly higher on dcn (5.37 vs 5.54).
Seconds per epoch (median, GPU, no concurrent GPU job): 0.20 (h200 d4), 0.17 (h300 d3). Task 5 wall time about 54 min.
Caveat: the long-cap gain is not yet bootstrapped; it is a descriptive comparison, not a tested claim.
