# Debugging reading (flash_C, dev rows)

MAE: Chemprop 16.83; in-house Task 7 (constant lr, 150 epochs, patience 25) 17.18.
Variants vs Chemprop (gain = Chemprop minus variant, 90 % CI): cosine 150: MAE 17.65, [-1.19, -0.32] (worse); plateau 150: MAE 17.41, [-0.93, -0.28] (worse); constant lr, 1000 epochs, patience 100: MAE 17.15, [-1.62, 0.49] (CI contains 0).
Reading (predeclared rule): the long-cap variant closes the gap, so the difference is attributed to the epoch budget (H2), not the schedule (H1: both decaying schedules were worse than constant lr within 150 epochs). H3 (featurisation) is untested and not excluded; the long-cap variant is still 0.32 degC above Chemprop (CI contains 0, not shown equal).
Consequence: the stop rule on Task 10 is lifted on condition that Task 10 uses epochs 1000, patience 100 (as the dmpnn search space already specifies).
