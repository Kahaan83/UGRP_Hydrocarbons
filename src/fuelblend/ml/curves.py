"""Per-epoch curve logger (CSV, optional TensorBoard) and a LightGBM callback.

Rows are appended and flushed per epoch so a partial file always parses. The header is fixed by the first row written;
later unknown metric keys are ignored.
"""
from __future__ import annotations

import csv
import os
import time
from pathlib import Path

BASE_COLS = ["run_id", "model", "fold", "seed", "epoch", "train_loss", "val_loss", "lr", "sec", "gpu_mem_mb"]


def _gpu_mb() -> float:
    try:
        import torch
        if torch.cuda.is_available():
            return float(torch.cuda.max_memory_allocated() / 2 ** 20)
    except Exception:
        pass
    return float("nan")


class CurveLogger:
    def __init__(self, csv_path, *, run: dict | None = None, tb_dir=None, _shared=None):
        self.run = {k: str(v) for k, v in (run or {}).items()}
        self.run.setdefault("run_id", "run")
        if _shared is not None:
            self._s = _shared
            return
        p = Path(csv_path)
        p.parent.mkdir(parents=True, exist_ok=True)
        tb = None
        if tb_dir is not None:
            try:
                from torch.utils.tensorboard import SummaryWriter
                tb = SummaryWriter(str(tb_dir))
            except Exception:
                tb = None
        self._s = {"path": p, "fh": None, "writer": None, "cols": None, "tb": tb, "t0": time.time()}

    def child(self, **run_extra) -> "CurveLogger":
        return CurveLogger(None, run={**self.run, **{k: str(v) for k, v in run_extra.items()}}, _shared=self._s)

    def log(self, epoch: int, **metrics: float) -> None:
        s = self._s
        row = {**{c: "" for c in BASE_COLS}, **self.run, "epoch": epoch}
        row["sec"] = metrics.pop("sec", time.time() - s["t0"])
        row["gpu_mem_mb"] = metrics.pop("gpu_mem_mb", _gpu_mb())
        row.update(metrics)
        if s["writer"] is None:
            exists = s["path"].exists() and s["path"].stat().st_size > 0
            if exists:
                with open(s["path"], newline="") as f:
                    s["cols"] = next(csv.reader(f))
            else:
                s["cols"] = BASE_COLS + sorted(k for k in row if k not in BASE_COLS and k not in self.run)
            s["fh"] = open(s["path"], "a", newline="")
            s["writer"] = csv.DictWriter(s["fh"], fieldnames=s["cols"], extrasaction="ignore")
            if not exists:
                s["writer"].writeheader()
        s["writer"].writerow(row)
        s["fh"].flush()
        if s["tb"] is not None:
            tag = "/".join(self.run.get(k, "") for k in ("model", "fold", "seed"))
            for k, v in metrics.items():
                try:
                    s["tb"].add_scalar(f"{tag}/{k}", float(v), epoch)
                except (TypeError, ValueError):
                    pass

    def close(self) -> None:
        s = self._s
        if s["fh"] is not None:
            s["fh"].close()
            s["fh"] = None
            s["writer"] = None
        if s["tb"] is not None:
            s["tb"].close()
            s["tb"] = None


def lgb_curve_callback(logger: CurveLogger):
    """Callback for lgb.train: logs the first validation metric per boosting round as val_loss (and the training metric)."""
    def _cb(env):
        ev = env.evaluation_result_list or []
        m = {}
        for name, metric, val, *_ in ev:
            key = "train_loss" if name in ("training", "train") else "val_loss"
            m.setdefault(key, float(val))
        logger.log(env.iteration + 1, **m)
    _cb.order = 30
    return _cb
