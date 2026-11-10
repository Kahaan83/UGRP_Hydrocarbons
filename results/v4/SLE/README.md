# v4 Task 13: ThermoML solid-liquid equilibrium (liquidus), ideal Schroeder-van Laar vs library Coutinho

All numbers below are recomputed from `sle_eval.csv`, `sle_summary.csv`, `sle_margin.csv` (written by `scripts/v4/eval_sle.py`).
Errors are predicted minus observed liquidus temperature, in K (negative = model predicts a COLDER liquidus than measured, the unsafe side for a freezing limit).

## Data
- Source: ThermoML archive `ThermoML.v2020-09-30.tgz` (sha256 `231161b5...58d7ec2`, matches `extras/fuelblend_data/data/raw/MANIFEST.json`). No newer release exists (docs/v4_manual_downloads.md section 4).
- Extractor: `thermoml_extract.py --sle` (patch kept in `thermoml_extract_sle.patch`, see Caveats on the repo location). Property `Solid-liquid equilibrium temperature, K`, LIQUID composition only, all-hydrocarbon systems, unsupported composition scales (solvent ratios) rejected. Default mode and its whitelist are unchanged.
- 11,923 XML files read; 1,058 hydrocarbon SLE rows from 23 files (1,051 binary, 7 ternary), saved as `data/processed/v4/thermoml_sle_hc.csv`. All rows are `Direct value, X` (absolute T, not a depression).
- Evaluation set: rows with P <= 150 kPa (drops 216 high-pressure rows: ethene + eicosane etc.) = 842 rows; 64 pure-component end points (x = 0 or 1) are not scored. **Scored: 778 points, 59 systems (57 binary, 2 ternary) from 21 files.** GroupKFold(5) by `system_id` (folds hold 154-158 points each).
- Pure constants: `data/processed/v4/pure_props_v4.csv`, `kind == pure` and `status == verified` only (R16); no substitution.

## Results (all scored points where the model is defined; pred - obs, K)
| model | systems | points | MAE | bias | RMSE | P90 abs |
|---|---|---|---|---|---|---|
| ideal SvL, mole fractions, verified constants | 8 | 103 | 2.02 | -1.10 | 3.10 | 4.81 |
| ideal SvL + fold-fitted offset | 8 | 103 | 3.29 | -1.63 | 4.11 | 6.40 |
| library Coutinho (mass fractions, as the repo calls it) | 44 | 525 | 9.97 | -7.16 | 18.18 | 28.79 |
| library Coutinho + fold-fitted offset | 44 | 525 | 10.39 | 0.00 | 16.77 | 22.27 |
| library Coutinho fed mole fractions | 44 | 525 | 13.06 | -11.03 | 19.27 | 29.24 |

The two model families are defined on different rows (ideal SvL needs verified Tm and dHfus for every component, 8 systems; Coutinho needs the library table to cover the crystallising component, 44 systems), so compare them only on the common subset (8 systems, 100 points):

| model (common 100 pts) | MAE | bias | RMSE | P90 abs |
|---|---|---|---|---|
| ideal SvL | 2.06 | -1.16 | 3.14 | 4.93 |
| ideal SvL + offset | 3.33 | -1.74 | 4.15 | 6.55 |
| Coutinho (mass) | 3.85 | -3.22 | 4.48 | 6.42 |
| Coutinho (mass) + offset | 3.74 | +3.40 | 4.76 | 8.38 |
| Coutinho (mole) | 6.19 | -6.09 | 7.25 | 10.99 |

By system class on the common subset: n-paraffin + n-paraffin (7 systems, 93 points) ideal MAE 2.11 / bias -1.34, Coutinho 4.08 / -3.48; cycloalkane + n-paraffin (1 system, 7 points) ideal 1.30 / +1.30, Coutinho 0.78 / +0.19 (one system, not a ranking). Coutinho's wider coverage shows up in the all-points table: aromatic + aromatic MAE 35.8 K (40 pts), aromatic + cycloalkane 24.2 K, aromatic + n-paraffin 10.6 K (188 pts), iso-paraffin + n-paraffin 1.3 K (67 pts), alkene + n-paraffin 3.6 K (33 pts); the full per-class table is `sle_summary.csv`.

Reading: on simple n-alkane blends ideal SvL is better than the library Coutinho (about 2 K vs 4 K MAE) and the fitted offset does not help (it makes the 8-system result worse: 3.3 K). Coutinho's depression term subtracts several K that the measured n-alkane liquidus does not show, which pushes it to the unsafe (too cold) side. Outside n-alkane/cycloalkane blends Coutinho errors are large (aromatic systems tens of K), consistent with its eutectic term being a fuel-wax calibration, not an SLE model.

## Caveats (honest limits)
1. **Coverage of verified constants is thin.** Only 8 of 59 systems (103 of 778 points) can be predicted by ideal SvL: n-alkane/n-alkane (7) and cyclohexane + hexadecane (1). Aromatic, naphthalene-type, benzene, C20+ alkanes and branched species have no verified Tm and dHfus in the v4 table, so nothing is claimed for them. The ideal SvL result is therefore an n-alkane (C6-C18) result.
2. **Coutinho is calibrated on CRC/Dauphin data** (library constants COUTINHO_ALPHA/BETA/GAMMA), and uses its own Tm/dHfus table, so it is not independent of the freezing work in Task 12; ThermoML SLE data are an out-of-source test of its eutectic term but its n-alkane polydispersity term is not what these binaries probe. Library code was used unmodified. The library silently skips components without constants; rows where the crystallising component (per the ThermoML PhaseID) is not in its table are excluded (69 of the 594 rows where it returned a number).
3. Library Coutinho takes a mass-fraction vector inside the SvL log (the repo convention); the mole-fraction variant is worse (6.19 vs 3.85 K MAE on the common subset). The mass-fraction log is not a physical SvL (SvL is defined on mole fraction); that the mass-fraction input scores better here is an empirical observation, not a validation of the convention.
4. Ideal SvL here ignores dCp and solid-solution formation; hexadecane/tetradecane/octadecane mixtures form solid solutions. The largest single miss is hexane + tridecane at x(C6) = 0.9659, T_obs 224.18 K vs ideal 211.84 K (+12.3 K under-prediction in the dilute-solute limit, crystal = tridecane).
5. 8 ideal-SvL systems is too few for a clean ranking by sub-class or for a reliable fitted offset; GroupKFold with 5 folds puts only 1-2 of them in a test fold. Same-system points from different papers are pooled in one group.
6. Constants matched by InChIKey connectivity (first 14 characters), so stereo isomers share constants; no SLE system here involves decalin or other stereoisomer pairs.
7. The extractor change lives in `extras/fuelblend_data/` which is a separate, git-ignored nested repository; it is not in the main-repo commit (patch saved here). `tests/test_v4_sle.py` imports it from that path.

## Proposal for the Task 10 freezing margin (proposal only, not adopted)
For a simple blend (all components n-alkane or cycloalkane with verified Tm and dHfus, composition known) the liquidus from ideal SvL, with NO fitted offset and NO Coutinho depression, can be used with a safety margin added to the predicted freezing temperature. Unsafe-side residuals (observed minus ideal-SvL prediction) over the 100 common points / 8 systems: median +0.83 K, 90th percentile +4.05 K, 95th percentile +7.36 K, maximum +12.34 K; the worst per-system mean is +3.41 K. Because no parameter is fitted, these residuals are out-of-fold by construction.
- Proposed: margin = 7 K (about the 95th point percentile) for such simple blends; 4 K if a 90 % level is acceptable.
- Coutinho on the same points is on the unsafe side more often (median +3.37 K, 95th percentile +8.72 K), so a margin derived this way would not be conservative-by-construction for Coutinho outputs.
- Not applicable to real fuels or to blends containing aromatics, naphthenes beyond cyclohexane, or heavy waxes; for those no ThermoML-derived margin exists and the Task 12 finding (no skill vs constant mean on CRC fuels) stands.

## Claims register (not edited here; for Task 19)
CLM-25 and CLM-28 (freezing/WAT models: add SLE evidence that ideal SvL beats library Coutinho on n-alkane binaries but only for 8 systems, Coutinho worse outside n-alkanes); CLM-18 and CLM-16 if they cite the Task 12 freezing numbers.
