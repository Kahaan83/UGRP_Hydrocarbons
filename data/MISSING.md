# Missing Data Log (`data/MISSING.md`)

Per the project ground rules:
> "Do not fabricate data or constants. If a value cannot be found in a credible source, leave it missing and log it in `data/MISSING.md`."

| Date Logged | Item / Source | Expected Content | Current Status & Mitigation |
|---|---|---|---|
| 2026-11-04 | NREL Cetane Compendium (Yanowitz et al. 2017, NREL/TP-5400-67585) | Digitized prioritized hydrocarbons for real-fuel components and candidate fuels | Updated `cetane_pure_curated.csv` to 88 pure hydrocarbons (covering 37/93 real-fuel components, 39.8%). Remaining unmeasured specialty isomers (56 components: C, C1=CCCCCCCCC1, C1C2CC3C4CC5CC(C14)C(C2)C3C5, C1CCCCCCCCCCC1, C1C[C@@H]2[C@H]3CC[C@H](C3)[C@@H]2C1, C=CCCCCCCCC=CC...) logged as unmeasured. |
| 2026-11-04 | 8 Unverified Antoine Sets (Nov 3 expansion) | NIST WebBook saturation vapor pressure tables for hexylbenzene, isocetane, octylbenzene, dodecylbenzene, 2,3,3-trimethylpentane, trans/cis-1,4-dimethylcyclohexane, cis-1,2-dimethylcyclohexane | Primary NIST WebBook tables do not exist for these 8 compounds. Flagged in `data/antoine_sources.csv` as `needs-removal-candidate`. |
| 2026-11-04 | RP-1 11-Blend Speciated Compositions (Fortin et al. 2012) | Speciated molecular formulations for the 11 orthogonal RP-1 blends | Confidential / proprietary blend recipes not published in open paper (only distillation curves published in Table 1-4). Reconstructed via POSF surrogates in hindcasts. |

<!-- v4-task3-dropped:begin -->

### v4 Task 3: fuels dropped from `data/real_fuels/v4/` (2026-11-10)

Missing: a raw-file, cited or reconstructed (R10) molecule list summing to >= 0.95 for each fuel. Rows are in `data/quarantine/`; full list with reasons in `data/quarantine/dropped_fuels.csv`; cause in `results/v4/P1/AUDIT.md`.

- **incomplete composition after quarantine (kept fraction < 0.95)** (8 fuels):
  - diesel (5): ASTM_D975_ULSD_Reference [0.540], SURR_EN590_European_Diesel [0.580], SURR_IDEAL_4comp [0.500], SURR_Mueller_4comp [0.520], SURR_Peugeot_4comp [0.500]
  - jet (3): BrunoHuber2010_3comp [0.690], Chen2016_JI [0.341], Honnet2009_Aachen [0.800]
- **no clean composition** (84 fuels):
  - diesel (11): CRC_FACE_Diesel_1, CRC_FACE_Diesel_2, CRC_FACE_Diesel_3, CRC_FACE_Diesel_4, CRC_FACE_Diesel_5, CRC_FACE_Diesel_6, CRC_FACE_Diesel_7, CRC_FACE_Diesel_8, CRC_FACE_Diesel_9, NASA_LFP3, NASA_LFP5
  - gasoline (55): CRC_FACE_A, CRC_FACE_A_E10, CRC_FACE_A_E25, CRC_FACE_A_E40, CRC_FACE_A_E60, CRC_FACE_B, CRC_FACE_C, CRC_FACE_C_E10, CRC_FACE_C_E25, CRC_FACE_C_E40, CRC_FACE_C_E60, CRC_FACE_D, CRC_FACE_E, CRC_FACE_F, CRC_FACE_F_E10, CRC_FACE_F_E25, CRC_FACE_F_E40, CRC_FACE_F_E60, CRC_FACE_G, CRC_FACE_G_E10, CRC_FACE_G_E25, CRC_FACE_G_E40, CRC_FACE_G_E60, CRC_FACE_H, CRC_FACE_I, CRC_FACE_I_E10, CRC_FACE_I_E25, CRC_FACE_I_E40, CRC_FACE_I_E60, CRC_FACE_J, CRC_FACE_J_E10, CRC_FACE_J_E25, CRC_FACE_J_E40, CRC_FACE_J_E60, CoOptima_4comp, CoOptima_4comp_E10, CoOptima_4comp_E20, CoOptima_4comp_E30, CoOptima_4comp_M10, CoOptima_4comp_M20, NREL_BOB_01, PRF_E10, PRF_E20, PRF_E30, PRF_E40, PRF_E50, PRF_E60, PRF_E70, PRF_E80, PRF_E90, SURR_TRF_60, SURR_TRF_70, SURR_TRF_80, SURR_TRF_90, TIGAS_Green_Gasoline
  - jet (18): Dooley2010_1stGen, Dooley2012_2ndGen, HyChem_NJFCP_A1, HyChem_NJFCP_A2, HyChem_NJFCP_A3, KimVioli2014_UM1, KimVioli2014_UM2, NASA_LFPA3, NASA_Shale_JP5, POSF_10264, POSF_10289, POSF_10325, POSF_11498, POSF_12344, POSF_12345, POSF_13511, POSF_13516, Violi2002_JP8_6comp

<!-- v4-task3-dropped:end -->
