# GH-GNN code and weights (WS5, 2026-11-14)

| item | value |
|---|---|
| source | https://github.com/edgarsmdn/GH-GNN (shallow clone to `repo/`, not committed) |
| commit | 0bbd2b317d86388c2ec7f531983d9d17f28342e0 (2023-08-31) |
| licence | MIT, Copyright (c) 2022 Edgar Sanchez (file `repo/LICENSE`); the README asks users to cite Sanchez Medina, Linke, Stoll, Sundmacher, Digital Discovery 2023, 2, 781-798, doi 10.1039/D2DD00142J |
| downloaded | 2026-11-14 |
| weights used | `repo/notebooks/GHGNN/GHGNN.pth` sha256 `091f70f4fb1f2b39678863eb5dae10c4fbd1390bd27fca45c09dc64edf956c14` |
| AD fingerprints | `repo/notebooks/GHGNN/training_fp.pickle` sha256 `d9fba3112981c5500bae9aa3f389b5d2b43474254049c9f2c98b16d5ee48de7e`; inspected with `pickletools` before loading: only the global `rdkit.DataStructs.cDataStructs.ExplicitBitVect`. The other three pickles (`training_classes`, `training_solutes`, `training_solvents`) contain no globals (strings, ints, lists) |
| training-data reference | `repo/data/raw/Brouwer_2021.csv` sha256 `7150705dd7a638f1bb67ed3a3e5e44ae80ee862446b2152afaefe166decca10c` (used for temperature-range reporting only) |
| pip dependencies | installed with `pip install --no-deps --target data/raw/ghgnn_code/site`: torch_geometric 2.8.1, mordredcommunity 2.0.7, tqdm 4.70.1, psutil 7.2.2, xxhash 4.1.0. The user's venv is unchanged. `torch` 2.15.0+cpu, `rdkit`, `networkx`, `numpy` come from the venv |
| torch_scatter | no wheel for torch 2.15 CPU on Windows; replaced by a pure-torch stand-in in `scripts/v5/ws5/shims/torch_scatter/` (scatter_add, scatter_mean over dim 0 only, which is all the GH-GNN architecture uses) |
| code read before running | `notebooks/GHGNN/ghgnn.py`, `GHGNN_architecture.py`, `mol2graph.py`, `scr/models/GHGNN_pred.py`, `scr/models/utilities/mol2graph.py` |
| reproduction check | `results/v5/BASE/delta_ws5_ghgnn/validation_repo_predictions.json` |
