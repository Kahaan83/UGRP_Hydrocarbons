# v5 WS4 GPU job pack (mixture ML: pds / gnn / st versus the LightGBM pair residuals)

Purpose: re-run `scripts/v5/base/delta_ws4_mixture_ml.py` with the declared longer schedule of 300 epochs, versus 60 in the CPU reading.
Declaration: `results/v5/BASE/delta_ws4_mixture_ml/README.md`.
**No GPU was available when this pack was written. Nothing here has been run on a GPU, and no GPU result is claimed.**
The CPU reading in `results/v5/BASE/delta_ws4_mixture_ml/<cand>/` is the only measured result.

## Requirements
- The Python 3.12 venv of this repo with a CUDA build of PyTorch (plain torch; PyG is not needed). RDKit, LightGBM, scikit-learn, scipy.
- Inputs. Several are gitignored, so copy them or mount the repo:
  - `data/cache/pure_grid.npz` (sha256 b931e969...; asserted by the baseline script);
  - `data/unified/measurements_dedup.csv`, `data/unified/compounds.csv`;
  - `data/real_fuels/v5/boehm2022_viscosity.csv`, `data/processed/v2/thermoml_clean_hc.csv`;
  - `results/v5/S1/viscosity_lab_residuals_v5.csv`, `results/v4/L/density_lab_residuals.csv`;
  - `results/v5/BASE/delta_ml_{viscosity,density}/rows_scored.csv`. These are the baseline errors, and they are read, not refitted.

## Run
    bash jobs/v5_gpu/run_ws4.sh cuda        # or: cpu ; EPOCHS=300 by default

Outputs go to `results/v5/BASE/delta_ws4_mixture_ml_gpu/<cand>/`, with the same files as the CPU reading. No checkpoints are written.

## Rules the job must keep
- The hyper-parameters are fixed in code and in the declaration: width 64, AdamW lr 2e-3, wd 1e-4, batch 512, cosine schedule, 5 seeds.
  Only the epoch count differs, and it is declared now as 300.
- Any tuning (width, lr, epochs) is done only by an INNER grouped split inside the dev systems, the grid is recorded, and the locked rows are
  never looked at for it.
- The locked systems were already scored once by the CPU reading. A GPU locked score reuses that set, so it is a second, pre-declared
  reading, not a fresh test. A GPU-only pass would need confirmation on a fresh holdout (e.g. new ThermoML data) before any designation
  is offered.
- Nothing is exposed in `predict_fuel`; designation is a user decision.

## What a GPU run adds, and what it does not
- **Adds:**
  - longer schedules (300+ epochs) and wider or deeper encoders (a full MolSets-width GNN, more SAB blocks);
  - more ensemble seeds, for better-calibrated ensemble spread;
  - an inner-dev hyper-parameter grid at a cost the shared 6-core CPU cannot carry;
  - the option to pre-train the GNN on pure-compound data.
- **Does not add data.** The binding constraint is the number of distinct systems: 264 / 167 dev binary systems, and 18 viscosity and
  21 density ternary systems (2 locked) for the extrapolation test. It is not compute. Leenhouts et al. report that about 10^3 unique mixtures
  are needed to beat ideal mixing.
- A longer schedule on this data mainly adds over-fitting risk.
- The roadmap's flash head and the Leenhouts benchmark splits (`splits.zip`, Zenodo 19914760, not on disk) are not part of this pack.
