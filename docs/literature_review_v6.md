# Literature review for v6 (corrected)

Scope: sources that motivate the v6 tracks (freezing point, pure-compound GNN, mixture ML retry). This replaces the offline
draft `master_literature_synthesis_and_ml_expansion.md`, which is not kept in the repo. Citations were checked on
2026-11-21 (web search and Consensus). Repo numbers come only from `docs/final_numbers_v5.md`. No v6 targets are set here;
targets are defined by the gates in `docs/v6_plan_2026-11-21.md`.

## 1. Repo reference point (v5)

- Freezing point on real fuels: UNIFAC SvL (designated route) reads `no_skill_shown`: n 50, MAE 8.425 C vs constant
  6.413 C, bias +5.878 C (`docs/final_numbers_v5.md`, section on freezing). Locked set: n 8, MAE 9.31 C vs 5.176 C,
  `worse_than_trivial`.
- Cold viscosity: status "no skill shown" (same file).
- v5 mixture ML (WS4) is treated as a weak implementation, not a final verdict. Its results stay documented as the v5 attempt.

## 2. Verified literature

### 2.1 Freezing point and GC x GC composition
- **Bell et al. (2025), Energy & Fuels.** "Freezing Point of Hydrocarbon Fuels from Single Species Concentrations."
  Control curves built from single-species concentrations; reported MAE 4.4 C. Basis for baseline (d) in Track A.
  Volume and page numbers were not recorded at verification time; add before citing in a manuscript.
- **Shi et al. (2017), Fuel.** Modified weighted average (MWA) over a class x carbon-number matrix; freezing MAE 0.82 C on
  17 fuels. Small n, so the figure is not comparable with a 50-fuel grouped CV. Baseline (c).
- **Watanasiri et al. (2026), Fuel.** Revised MWA; freezing is the best-performing property. Baseline (c) variant.
- **Hall, Creton, Rauch, Bauder, Aigner (2022), Energy & Fuels 36(1), 463-479.** Probabilistic mean QSPR (M-QSPR) from
  GC x GC: each tile is represented by a probability-weighted average over its isomers, not one stand-in compound.
  Reports that freezing point and viscosity need real fuels in training. Basis for isomer-averaged tiles.
- **Boehm, Coburn, Yang, Wanstall, Heyne (2022), Energy & Fuels 36(17), 9980-9991.** Freeze-point blend model for jet-range
  hydrocarbons. Its cold-viscosity data are already in production (`loaders_boehm2022.py`).
- **Shao et al. (2024), Chem. Eng. Sci.** ML on 104 jet fuels. **Dong et al. (2026)**, SISSO on 131 fuels with GC x GC.
  Both are comparison points only; full bibliographic details to be added.
- **Coutinho (1998), Ind. Eng. Chem. Res. 37(12), 4870-4875.** Predictive UNIQUAC for solid-liquid equilibria; n-alkanes
  form non-ideal solid solutions. Context for the SvL baseline.

### 2.2 Mixture ML and pooling
- **Leenhouts, Larsson, Verhelst, Vermeire (2025), Fuel 381, 133218**, pooled GNNs for fuel mixtures (MolPool); and
  **Leenhouts, Morgan, Ibrahim, Green, Vermeire (2025), Chem. Eng. J. 513, 162232**, MolPool for solvent mixtures.
  Mole-fraction-weighted pooling is the Track C encoder design.
- **Schweidtmann et al. (2023), Comput. Chem. Eng. 172, 108202.** Physical pooling functions. Use only the supported
  point: intensive properties need fraction-normalised pooling, extensive ones need sum pooling. No "theorem" wording.
- **Zhang et al. (2024), PRX Energy 3(2), 023006** (MolSets). Permutation-invariant mixture model with set attention.
- **Dobbelaere et al. (2022), Ind. Eng. Chem. Res. 61(25), 8581-8594.** ML for hydrocarbon-mixture properties; linear mixing
  holds for density, not for viscosity, vapour pressure or flash point.
- **Larsson, Vermeire, Verhelst (2023), SAE 2023-01-0337.** Multi-task and transfer learning from pure-compound data.
  Basis for Track B (multi-task heads) and the transfer step into Track C.
- **Schweidtmann et al. (2020), Energy & Fuels 34(9), 11395-11407.** GNNs for ignition quality (DCN, RON, MON).
- **Xin et al. (2026), Front. Chem. Sci. Eng.**, DOI 10.1007/s11705-026-2698-7. Multimodal representation with latent-space
  mixing for RON of pure and mixed fuels. Octane head design in Track C.
- **Kuzhagaliyeva et al. (2022), Commun. Chem. 5, 111.** Linear mixing in latent space with multi-property decoders.
- **Correa Gonzalez et al. (2021), Energy & Fuels 35(11), 9332-9340.** Octane blending on molar basis is near-linear;
  volumetric basis and polar/aromatic species need non-linear terms. Linear molar is a required baseline.
- **Oh, Oldani, Solecki, Lee (2024), Fuel 356, 129508.** Bayesian neural-network ensemble for flash point from GC x GC.
  This replaces the unverifiable Radaideh and Violi 2025 citation.
- **ChemFlow (arXiv 2603.02810).** Preprint only, not peer reviewed; cite as such.

### 2.3 Thermodynamics and fuel context (unchanged from v5 basis)
Gmehling and Rasmussen (1982) and Liaw and Chiu (2003) for flash point; Quinones-Cisneros et al. (2000) and Allal et al.
(2001) for viscosity; Das et al. (2018) for YSI; Colket et al. (2017) and Flora et al. (2019) for NJFCP reference fuels and
blend optimisation; Bruno (2006) for distillation curves. Page and DOI details are in the offline draft and have not been
re-verified; re-check before manuscript use.

## 3. Removed from the offline draft

- Citations not found: LG-Transformer (Nat. Commun. "Vol 17 Art 1234"), Radaideh and Violi 2025 (ECM:X), Abraham et al. 2025.
- Wrong DOI fixed: Xin et al. (-2698-7, not -2702-7).
- False repo claims: the 48 % "sum pooling" story (WS4 never used plain sum pooling; 48 % was v2/v3 at -40 C); bicyclohexyl
  and P0y-b (the 33.6 C over-prediction is a lab binary, and P0y-b Jet A-1 is certified); the -0.50 ln bias (pooled value
  below 263 K, not a general figure); the +0.70/+1.09 C flash figures (GH-GNN only; mesitylene belongs to the ideal flash route).
- Unsourced: the +2.3 C D93 offset (fitting it on our own errors would leak) and invented v6 targets (MAE below 2.5 C,
  8 of 13 certified).
- Overreach: "Schweidtmann theorem".

## 4. v6 ideas carried forward

| Idea | Source | Where |
|---|---|---|
| Isomer-averaged tiles | Hall 2022 | Track A baseline (e) |
| Single-species freezing control curves | Bell 2025 | Track A baseline (d) |
| Class x carbon-number MWA matrix | Shi 2017, Watanasiri 2026 | Track A baseline (c) |
| Mole-fraction pooling | Leenhouts 2025 (MolPool) | Track C |
| Multi-task and pure-compound transfer | Larsson 2023 | Track B, C |
| Latent mixing head for octane | Kuzhagaliyeva 2022, Xin 2026 | Track C octane |
| Joint conformal regions for design | none verified; method idea, to be tested | after Track C |

Gates for every item: GroupKFold on `system_id`, locked holdout, beat-trivial, out-of-fold conformal margins.

Note (v7): the years given for Leenhouts et al. and Hall et al. are online dates; print years may differ by one. Pending a DOI check, do not quote them as final.
