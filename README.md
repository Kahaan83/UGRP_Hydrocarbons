# UGRP Hydrocarbons — physics-first fuel-blend property prediction

Predict flash point, viscosity, density, freezing point, heat of combustion and smoke point of a hydrocarbon blend from its
composition, with a measured error margin on every number and an honest "no skill" label wherever a model does not beat a constant.
Full explanation: `docs/GUIDE.md`. Methodology: `Fuel_Blend_Methodology.pdf`. All numbers: `docs/final_numbers_v5.md` (blends), `docs/final_numbers_v7.md` (pure compounds). Current methods: `docs/methodology_v7.md`.

## Quick start (Windows, venv python)

```powershell
.venv\Scripts\python.exe -m pip install -r requirements.txt
.venv\Scripts\python.exe scripts/v5/base/reproduce_base.py        # rebuild scoreboard + BASE tests + hash check (minutes)
.venv\Scripts\python.exe -m pytest -q                              # full test suite (long)
```

Predict a blend (SMILES or InChIKey, mass fractions):

```python
from fuelblend.api.predict_route import predict_fuel          # run with PYTHONPATH=src
r = predict_fuel([("CCCCCCCCCCCC", .45), ("CCCCCCCCCCCCCC", .20), ("Cc1ccccc1", .15), ("CC1CCCCC1C", .20)], T_K=288.15, fuel_type="jet")
r["flash_point_C"]["value"], r["flash_point_C"]["method"], r["flash_point_C"]["in_domain"]
```

## Where things are

| Folder | What is in it |
|---|---|
| `src/fuelblend/baselines/` | the physics: Le Chatelier flash (`flash.py`), freezing, heat of combustion, smoke, D86, cetane, octane |
| `src/fuelblend/mixture/common.py` | ideal-volume density, Refutas viscosity, pure-property lookup |
| `src/fuelblend/eval/composition_route.py` | the **designated scoring route**: how a composition becomes every property |
| `src/fuelblend/ml/` | the learned corrections: pair viscosity, pair excess volume (kept); ML smoke/YSI (failed gate) |
| `src/fuelblend/uq/conformal.py` | the error margins (group-level conformal, one-sided) |
| `src/fuelblend/api/` | `predict_route.py` (`predict_fuel`) and `predict_blend.py` (`predict`) — the public entry points |
| `src/fuelblend/design/` | inverse design on the composition simplex (`engine_v4.py`, `library_route.py`) |
| `src/fuelblend/data/`, `splits/`, `features/` | loaders for each source, split definitions, descriptors |
| `scripts/v5/base/` | reproduces the frozen scoreboard (`reproduce_base.py`) and every `delta_*` experiment |
| `scripts/v5/run_design_userpool_v5.py`, `run_jet_designs_v5.py` | inverse design runs |
| `scripts/v4/` | earlier evidence pipeline that v5 builds on (kept because `src/` and tests use it) |
| `results/v5/BASE/` | frozen v5 evidence, one folder per experiment, protected by `SHA256SUMS_v5_base.txt` |
| `results/v4`, `results/t*`, `results/pairwise` | earlier results still read by code/tests |
| `data/`, `models/v5/` | unified measurements, real-fuel tables, trained pair models |
| `docs/` | `GUIDE.md` (start here), `final_numbers_v5.md`, `claims_register.md`, `model_decision.md`, `v5_negative_results.md`, `FULL_REPORT_v4.md` |
| `tests/` | pytest suites (`test_base_*` guard the frozen scoreboard) |

## Status in one paragraph

Physics baselines are the backbone; two learned corrections (viscosity pair model, excess-volume pair model) passed pre-declared gates
and are used inside their domain. Freezing, smoke and diesel models show no skill over a constant on real fuels and are labelled so.
Inverse design returns 16 blends; 3 pass the certification test **under the model only** (no laboratory validation was done).
Failed approaches (neural mixture encoders, GNN flash activity, GCxGC models, ML YSI, ML DCN) are kept on purpose; see `docs/v5_negative_results.md`.

The full pre-simplification repository is the git tag `v5-full-2026-11-16`.
