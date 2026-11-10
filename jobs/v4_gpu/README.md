# v4 Task 16 - GPU job pack (set encoder)

Purpose: re-run `scripts/v4/train_set_encoder.py` with a much longer schedule than the CPU run recorded in `results/v4/E/`.
**No GPU was available when this pack was written; nothing here has been run on a GPU and no GPU result is claimed.**
The CPU result in `results/v4/E/README.md` (40 epochs, 5 seeds, 3 threads, hand-fixed hyper-parameters) is the only measured result.

## Requirements
Python 3.12 venv of this repo plus a CUDA build of PyTorch; the same data the CPU run needs
(`data/unified/compounds.csv`, `data/unified/measurements_dedup.csv`, `data/cache/pure_grid.npz`, `data/real_fuels/v4/*`,
`results/v4/F/flash_*`, `results/v4/V/viscosity_predictions_v4.csv`, `splits/v4/real_fuels_v4_splits.json`).
Copy the repo (incl. those untracked/ignored data folders) to the GPU machine, or mount it.

## Run
    bash jobs/v4_gpu/run_all.sh cuda        # or: cpu

Outputs: `results/v4/E_gpu_attention/`, `results/v4/E_gpu_sum/` (same CSV schema as `results/v4/E/`), checkpoints in `runs/v4/gpu/`
(`runs/` is gitignored; do not commit `.pt`).

## Rules the job must keep
* The 12 locked CRC 647 fuels are never in a training, selection or calibration index (`split_audit.json` records the checks; the unit test enforces it).
* Hyper-parameters are fixed in code (hidden 64, AdamW lr 3e-3, wd 1e-3, cosine schedule, no early stopping). If you tune epochs/width/lr on the GPU,
  do it by an INNER grouped split inside the dev fuels only, record the grid, and never look at locked rows.
* Promotion to `model_decision` still requires `promotion_decision.json: PROMOTE_TO_MODEL_DECISION = true`
  (paired fuel-cluster bootstrap CI vs Task 11's structured ridge excludes 0 under both protocols AND dev cross-conformal row coverage >= 0.85).
* Interpret a longer schedule with care: the dev set is 45 fuels, so more epochs mainly add over-fitting risk, not information.
