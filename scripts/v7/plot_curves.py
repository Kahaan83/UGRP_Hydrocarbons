"""plot_curves.py <curves.csv> --out <png>: median and IQR train/val loss across seeds per fold, best epoch marked.

Works on a partial file from a run that is still going."""
import argparse

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("csv")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    d = pd.read_csv(a.csv)
    if "fold" not in d.columns or d.fold.isna().all():
        d["fold"] = 0
    folds = sorted(d.fold.unique())
    fig, axs = plt.subplots(1, len(folds), figsize=(4 * len(folds), 3.2), squeeze=False, sharey=False)
    for ax, f in zip(axs[0], folds):
        g = d[d.fold == f]
        for col, c in (("train_loss", "tab:blue"), ("val_loss", "tab:red")):
            if col not in g or g[col].isna().all():
                continue
            q = g.groupby("epoch")[col].quantile([.25, .5, .75]).unstack()
            ax.plot(q.index, q[.5], color=c, label=col)
            ax.fill_between(q.index, q[.25], q[.75], color=c, alpha=.2)
        if "val_loss" in g and g.val_loss.notna().any():
            m = g.groupby("epoch").val_loss.median()
            ax.axvline(m.idxmin(), color="k", ls=":", lw=1)
        ax.set_title(f"fold {f}"); ax.set_xlabel("epoch")
    axs[0][0].legend(); fig.tight_layout(); fig.savefig(a.out, dpi=120)


if __name__ == "__main__":
    main()
