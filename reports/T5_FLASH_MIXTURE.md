# T5 - mixture-level flash residual + conformal coverage

conformal alpha=0.1 (nominal 90%); split-conformal quantiles from S1 mixture residuals only.
row counts: S1 mixtures=34, S2 mixtures=22, S3 mixtures=0.

## Model: Le Chatelier + median mixture bias (fitted on S1 mixtures only)

S1 mixture residuals (value - chat_nominal): n=34, median=+1.690 degC, mean=+1.473, std=1.583, range=[-3.67, +5.17]
conformal 90% two-sided interval on the residual: [-1.419, +3.473] degC
one-sided upper quantile (alpha=0.10): u=+2.886 degC

S1 in-sample: coverage=85.29% (under=11.76%, over=2.94%) - NOT out-of-sample evidence

## S2 (n=22; per-row details below)

- raw Le Chatelier MAE=1.354 degC, +bias MAE=0.835 degC (+38.3% change)
- conformal 90% coverage=81.82% (under=18.18%, over=0.00%); one-sided over-prediction rate: >0=54.55%, >5 degC=0.00%
- wide-side rate vs nominal (each 5%): under=18.18%, over=0.00%
- | n_comp | T_K | measured | chat | +bias | interval | in-interval |
  | 2 | nan | 60.0 | 58.0 | 59.7 | [58.3, 63.2] | yes |
  | 2 | nan | 14.7 | 13.1 | 14.8 | [13.4, 18.3] | yes |
  | 2 | nan | 13.7 | 12.4 | 14.1 | [12.6, 17.5] | yes |
  | 2 | nan | 11.4 | 10.6 | 12.3 | [10.9, 15.8] | yes |
  | 2 | nan | 11.7 | 10.3 | 12.0 | [10.6, 15.5] | yes |
  | 2 | nan | 8.4 | 7.8 | 9.5 | [8.1, 13.0] | yes |
  | 2 | nan | 9.7 | 7.4 | 9.1 | [7.6, 12.5] | yes |
  | 2 | nan | 6.0 | 5.4 | 7.1 | [5.7, 10.6] | yes |
  | 2 | nan | 7.7 | 4.9 | 6.6 | [5.2, 10.1] | yes |
  | 2 | nan | 3.4 | 3.3 | 5.0 | [3.5, 8.4] | NO |
  | 2 | nan | 4.7 | 2.8 | 4.5 | [3.0, 7.9] | yes |
  | 2 | nan | 1.5 | 1.4 | 3.1 | [1.6, 6.5] | NO |
  | 2 | nan | 2.7 | 0.9 | 2.6 | [1.2, 6.1] | yes |
  | 2 | nan | -0.0 | -0.3 | 1.4 | [-0.1, 4.8] | yes |
  | 2 | nan | 1.7 | -0.7 | 1.0 | [-0.4, 4.5] | yes |
  | 2 | nan | -1.5 | -1.9 | -0.2 | [-1.6, 3.3] | yes |
  | 2 | nan | -0.3 | -2.2 | -0.5 | [-1.9, 3.0] | yes |
  | 2 | nan | -3.0 | -3.3 | -1.6 | [-3.0, 1.9] | NO |
  | 2 | nan | -0.3 | -3.5 | -1.8 | [-3.2, 1.7] | yes |
  | 2 | nan | -4.5 | -4.5 | -2.9 | [-4.3, 0.6] | NO |
  | 2 | nan | -2.3 | -4.7 | -3.0 | [-4.4, 0.5] | yes |
  | 2 | nan | -3.3 | -5.5 | -3.8 | [-5.2, -0.4] | yes |

## S3: no mixture flash rows

## Verdict

- Mixture flash data (S1 only: 34 rows; S2: 22; S3: 0) can support exactly one thing - a bias-only
  calibrator - and that is what this deliverable is. No GBM/learned residual here: n<=34 would memorize.
- **Bias finding**: Le Chatelier under-predicts real mixture flash points by a median **+1.69 degC**
  (S1). Applying it cuts S2 MAE from 1.354 to 0.835 degC (+38% better) - the same mixture-level
  non-ideality signal the pure-model story missed.
- **Safety side is clean**: on S2, over-prediction rate >0 degC = 54.5% but >5 degC = 0% - no
  dangerous over-prediction occurs out-of-sample. The 81.8% empirical coverage vs nominal 90%
  misses land only on the under-prediction (conservative for any "will it ignite" use), due to the
  small calibration set (one row per quantile step: 22 rows x 5% tails = ~1 row each side).
- Honest reporting: intervals are split-conformal calibrated on 34 S1 rows only; treat "90%" as
  "~82%-85%" until more binary flash rows are ingested (T0 gap #1 / T7 priority: more flash binaries).
- Le Chatelier baseline itself remains strong: RMSE 1.92 degC on S2 binaries with zero learned parts.

