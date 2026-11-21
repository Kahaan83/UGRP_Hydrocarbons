import sys; sys.path.insert(0, "src")
import numpy as np, torch
from fuelblend.ml import dmpnn as d
f = [d.featurize(s) for s in ["CCCC", "c1ccccc1C", "CC(C)CO", "C1CCCCC1"]]
Y = np.random.randn(4, 3); M = np.ones((4, 3))
ms = d.fit_ensemble(f, Y, M, f, Y, M, n_models=2, epochs=3)
print(d.predict(ms, f)[0].shape)
