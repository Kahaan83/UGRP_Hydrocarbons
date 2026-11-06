# Missing Data Log (`data/MISSING.md`)

Per the project ground rules:
> "Do not fabricate data or constants. If a value cannot be found in a credible source, leave it missing and log it in `data/MISSING.md`."

| Date Logged | Item / Source | Expected Content | Current Status & Mitigation |
|---|---|---|---|
| 2026-11-04 | NREL Cetane Compendium (Yanowitz et al. 2017, NREL/TP-5400-67585) | Digitized prioritized hydrocarbons for real-fuel components and candidate fuels | Updated `cetane_pure_curated.csv` to 88 pure hydrocarbons (covering 37/93 real-fuel components, 39.8%). Remaining unmeasured specialty isomers (56 components: C, C1=CCCCCCCCC1, C1C2CC3C4CC5CC(C14)C(C2)C3C5, C1CCCCCCCCCCC1, C1C[C@@H]2[C@H]3CC[C@H](C3)[C@@H]2C1, C=CCCCCCCCC=CC...) logged as unmeasured. |
| 2026-11-04 | 8 Unverified Antoine Sets (Nov 3 expansion) | NIST WebBook saturation vapor pressure tables for hexylbenzene, isocetane, octylbenzene, dodecylbenzene, 2,3,3-trimethylpentane, trans/cis-1,4-dimethylcyclohexane, cis-1,2-dimethylcyclohexane | Primary NIST WebBook tables do not exist for these 8 compounds. Flagged in `data/antoine_sources.csv` as `needs-removal-candidate`. |
| 2026-11-04 | RP-1 11-Blend Speciated Compositions (Fortin et al. 2012) | Speciated molecular formulations for the 11 orthogonal RP-1 blends | Confidential / proprietary blend recipes not published in open paper (only distillation curves published in Table 1-4). Reconstructed via POSF surrogates in hindcasts. |
