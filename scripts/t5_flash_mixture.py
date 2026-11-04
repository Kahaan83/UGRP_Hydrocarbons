#!/usr/bin/env python
"""T5: mixture-level flash residual model + conformal coverage.

    python scripts/t5_flash_mixture.py

Data is intentionally thin (mixture flash: ~54 binary + 3 multi applicable rows, see
reports/T0_COVERAGE.md), so the "model" here is deliberately non-learned: Le Chatelier
(baseline `chat_nominal`) plus an additive median bias fitted ONLY on S1 mixture rows,
with split-conformal intervals taken from the SAME S1 mixture residuals. Coverage is
then checked once on S2 (and reported for S3, never tuned on).

Writes reports/T5_FLASH_MIXTURE.md.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

BASE = Path(__file__).resolve().parents[1]
PRED = BASE / "data/runs/baseline_preds.csv"
OUT = BASE / "reports/T5_FLASH_MIXTURE.md"
ALPHA = 0.10  # nominal 90% coverage, mirrors WP6


def flash_rows(df, split):
    d = df[(df["property"] == "flash_point_C") & (df["split"] == split)
           & (df["conflict"] == False) & df["chat_nominal"].notna()  # noqa: E712
           & df["value"].notna() & (df["kind"] != "real_fuel")]
    return d.copy()


def coverage(y, pred, lo, hi):
    return float(np.mean((y >= lo) & (y <= hi))), float(np.mean(y < lo)), float(np.mean(y > hi))


def main() -> None:
    df = pd.read_csv(PRED, low_memory=False)
    s1, s2, s3 = flash_rows(df, "S1"), flash_rows(df, "S2"), flash_rows(df, "S3")
    mix_s1 = s1[s1["n_comp"] >= 2]
    mix_s2, mix_s3 = s2[s2["n_comp"] >= 2], s3[s3["n_comp"] >= 2]

    lines = ["# T5 - mixture-level flash residual + conformal coverage", "",
             f"conformal alpha={ALPHA} (nominal 90%); split-conformal quantiles from S1 mixture residuals only.",
             f"row counts: S1 mixtures={len(mix_s1)}, S2 mixtures={len(mix_s2)}, S3 mixtures={len(mix_s3)}.", ""]

    if not len(mix_s1):
        OUT.write_text("\n".join(lines + ["No S1 mixture flash rows - cannot calibrate anything."]) + "\n", encoding="utf-8")
        print("no S1 mixture rows"); return

    r = (mix_s1["value"].to_numpy(float) - mix_s1["chat_nominal"].to_numpy(float))
    bias = float(np.median(r))
    q_lo, q_hi = np.quantile(r, [ALPHA / 2, 1 - ALPHA / 2])
    q_up = np.quantile(r, 1 - ALPHA)  # safety side: pred must not exceed true+q_up

    pred_s1 = mix_s1["chat_nominal"].to_numpy(float)
    lines += ["## Model: Le Chatelier + median mixture bias (fitted on S1 mixtures only)", ""]
    lines.append(f"S1 mixture residuals (value - chat_nominal): n={len(r)}, median={bias:+.3f} degC, "
                 f"mean={r.mean():+.3f}, std={r.std(ddof=1):.3f}, range=[{r.min():+.2f}, {r.max():+.2f}]")
    lines.append(f"conformal 90% two-sided interval on the residual: [{q_lo:+.3f}, {q_hi:+.3f}] degC")
    lines.append(f"one-sided upper quantile (alpha=0.10): u={q_up:+.3f} degC")
    lines.append("")

    # train-side in-sample diagnostic (honest note: not out-of-sample)
    cov, under, over = coverage(mix_s1["value"].to_numpy(float), pred_s1 + bias, pred_s1 + bias + q_lo, pred_s1 + bias + q_hi)
    lines.append(f"S1 in-sample: coverage={cov:.2%} (under={under:.2%}, over={over:.2%}) - NOT out-of-sample evidence")
    lines.append("")

    for name, d in (("S2", mix_s2), ("S3", mix_s3)):
        if not len(d):
            lines.append(f"## {name}: no mixture flash rows"); lines.append("")
            continue
        y = d["value"].to_numpy(float)
        p_raw = d["chat_nominal"].to_numpy(float)
        p_bias = p_raw + bias
        lo, hi = p_bias + q_lo, p_bias + q_hi
        cov, under, over = coverage(y, p_bias, lo, hi)
        mae_raw = float(np.mean(np.abs(y - p_raw)))
        mae_bias = float(np.mean(np.abs(y - p_bias)))
        over_rate = float(np.mean(p_bias > y))
        over5 = float(np.mean(p_bias - y > 5))
        lines += [f"## {name} (n={len(d)}; per-row details below)", ""]
        lines.append(f"- raw Le Chatelier MAE={mae_raw:.3f} degC, +bias MAE={mae_bias:.3f} degC ({(1 - mae_bias/mae_raw)*100:+.1f}% change)")
        lines.append(f"- conformal 90% coverage={cov:.2%} (under={under:.2%}, over={over:.2%}); "
                     f"one-sided over-prediction rate: >0={over_rate:.2%}, >5 degC={over5:.2%}")
        lines.append(f"- wide-side rate vs nominal (each {(ALPHA/2):.0%}): under={under:.2%}, over={over:.2%}")
        lines.append("- | n_comp | T_K | measured | chat | +bias | interval | in-interval |")
        for i, row in d.reset_index(drop=True).iterrows():
            inside = "yes" if lo[i] <= y[i] <= hi[i] else "NO"
            lines.append(f"  | {int(row.n_comp)} | {row.T_K:.2f} | {y[i]:.1f} | {p_raw[i]:.1f} | {p_bias[i]:.1f} | [{lo[i]:.1f}, {hi[i]:.1f}] | {inside} |")
        lines.append("")

    lines += ["## Verdict template (fill after run)", "",
              "The bias + conformal wrapper is the only defensible mixture-flash model at this n;",
              "a GBM residual on <=31 mixture rows would memorize. Quote coverage per split above.",
              ""]
    OUT.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
