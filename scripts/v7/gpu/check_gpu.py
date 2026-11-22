"""Print CUDA status and a 1-epoch D-MPNN timing on CPU and (if available) GPU. Run in the GPU venv with PYTHONPATH=src."""
import os
import sys
import time
from pathlib import Path

import numpy as np
import torch

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "src"))


def one_epoch(dev: str) -> float:
    os.environ["FUELBLEND_DEVICE"] = dev
    import importlib
    from fuelblend.ml import dmpnn as D
    importlib.reload(D)
    smi = ["CCCCCCCC", "Cc1ccccc1", "C1CCCCC1", "CC(C)CCCC", "CCCCCCCCCC", "c1ccc2ccccc2c1"] * 10
    f = [D.featurize(s) for s in smi]
    Y = np.random.default_rng(0).normal(size=(len(f), 2)); M = np.ones_like(Y, bool)
    t = time.time()
    D.fit_ensemble(f, Y, M, f[:6], Y[:6], M[:6], n_models=1, epochs=1)
    if dev == "cuda":
        torch.cuda.synchronize()
    return time.time() - t


ok = torch.cuda.is_available()
print("cuda", ok, torch.cuda.get_device_name(0) if ok else "none", "torch", torch.__version__, torch.version.cuda)
print("cpu 1 epoch s", round(one_epoch("cpu"), 3))
if ok:
    one_epoch("cuda")  # warm-up
    print("gpu 1 epoch s", round(one_epoch("cuda"), 3))
