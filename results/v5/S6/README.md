# S6 diesel: Day-9 gate and scoring (v5)

Status: SCORED 2026-11-11 (readings declared before scoring; see Results). Declared 2026-11-11. Do not edit the readings after scoring; add a dated note.

## Sources (downloaded 2026-11-11 with the user's go-ahead, raw dir gitignored: `data/raw/diesel/`)

* NIST TN 1949, https://nvlpubs.nist.gov/nistpubs/TechnicalNotes/NIST.TN.1949.pdf (2.7 MB): Table 1 exact mole fractions of the four CRC AVFL-18a surrogates V0a, V0b, V1, V2; Table 11 measured kinematic viscosity at 40 C.
* CRC FACE-1 (2010), https://crcao.org/wp-content/uploads/2019/06/FACE-1-Chem-and-Phys-Props-of-FACE-Research-Diesel-Fuels.pdf (7.5 MB, 242 pages): nine FACE diesels, 2-D GC group types (Table 5.2), PIONA/SOAP (7.1, 8.1-8.3), DHA (10.1), ASTM properties.

## Questions

A. (Gate, plan Day 9) Do these sources give at least 5 real diesel fuels with raw-file composition plus at least one measured target?
B. If yes, how well do the existing models predict diesel density, 40 C viscosity, flash and cetane/DCN from the measured composition?
C. How well do the four known-composition surrogates (T1 evidence) reproduce their measured 40 C viscosity?

## Readings (fixed before scoring)

R-S6-1 (gate): a fuel QUALIFIES if its composition comes from a table in the PDF (not reconstructed from distillation or bulk descriptors), is at least class-level (paraffin/iso/cyclo/aromatic or finer), sums to 95-105 percent after transcription, and has at least one measured target (density, viscosity 40 C, flash, cloud, cetane/DCN) in the same document. The gate is met if the count is 5 or more. The count and each fuel's reason for in/out are reported. No reconstructed diesel composition is created or used under any outcome.
R-S6-2 (only if the gate is met): with measured composition mapped to representative compounds, report bias, MAE and n for density at 15 C, ln viscosity at 40 C (Refutas), flash (ideal Le Chatelier) and cetane/DCN where the repo has a model. For viscosity also report the fraction of fuels inside the lab-mixture 40 C margin already in `configs/uncertainty_v4.yaml` (0.186 ln). Evidence label: T2-diesel, n stated, no claim stronger than n supports (n < 8 is "descriptive only").
R-S6-3 (surrogates): for V0a, V0b, V1, V2 compute Refutas 40 C viscosity from the Table 1 compositions and compare with Table 11; report ln error each, and whether each is inside 0.186 ln. n = 4: descriptive only.
R-S6-4: if the gate is not met, diesel is reported as known-composition evidence only (T1 covers C10-C20 molecules plus R-S6-3) and the two diesel designs are NOT re-run.
R-S6-5 (transcription check): for every table transcribed, at least 10 values (or all, if fewer) are compared against the PDF page image; any mismatch is corrected and listed. The check is recorded in the results folder.
R-S6-6: nothing is written to `configs/`. Diesel design re-runs (if the gate is met) are a separate step with their own declared reading and are not part of this scoring.

A reading that cannot be evaluated is stated "not evaluable" with the reason.

## Results

Scored 2026-11-11 with `scripts/v5/s6_diesel.py` (`build`, `gate`, `score`, `check`); tests `tests/test_s6_diesel.py` (8 pass). The declared readings above were not changed. Files: `gate_decision.csv`, `transcription_check.csv`, `composition_internal_checks.csv`, `surrogate_viscosity_R-S6-3.csv`, `mapping_summary.csv`, `mapped_representatives.csv`, `predictions_R-S6-2.csv`, `summary_R-S6-2.csv`, `summary_R-S6-2_sensitivity_all9.csv`; data in `data/real_fuels/v5/diesel_*.csv`.

**Transcription (R-S6-5).** 359 values compared by eye against rendered page images (NIST Table 1: 25, Table 11: 8 (all); FACE Tables 5.2: 18, 2.2/3.5 properties: about 106, Appendix B.1-B.9: 12-13 each, Appendix C.1-C.9: 10-12 each). 0 mismatches, so nothing was corrected. A separate internal-consistency check (printed row/column totals vs summed cells, 99 checks) found only 3 rounding differences. Note: the PDF labels FD5A as "Face5" in Table 5.2 (typo; treated as FD5A).

**R-S6-1 gate: MET, count 9 of 9 FACE diesels** (FD1A-FD9A). Each has a table-sourced composition by hydrocarbon type x carbon number (GC-FIMS C8-C21 for the >200 C part, PIONA C3-C11 for the <200 C part), sums 99.88-100.03 percent, and measured density, 40 C viscosity, flash, cloud and cetane (D613 and D6890 DCN) in the same PDF. The four NIST surrogates (V0a, V0b, V1, V2) are not real fuels and are not counted (they are scored under R-S6-3). No reconstructed composition was used. Composition is per carbon number (not class-only), and Table 5.2 gives a class-level cross-check.

**R-S6-3 surrogates (n = 4, descriptive only).** Refutas ln(pred/exp) at 40 C: V0a +0.136 (2.810 vs 2.452 cSt), V0b +0.108, V1 +0.115; all three are inside 0.186 ln, all over-predict. V2 is not evaluable: 28.8 percent of its moles (n-eicosane, 2-methylheptadecane, 1,3,5-triisopropylbenzene, perhydrophenanthrene) lack valid pure eta(T) data in the registry (a covered-only partial blend gives +0.045 but is not a valid evaluation). For comparison, the NIST REFPROP model in Table 11 is within -0.084 to +0.008 ln of experiment.

**R-S6-2 (T2-diesel, n = 4, descriptive only).** Eligibility (my choice, same 0.10 cut-off as G1, set in the script before scoring and not in the declared readings): unmapped mass of at most 10 percent. Eligible: FD1A, FD2A, FD3A, FD5A. Ineligible: FD4A 11 percent, FD6A 22 percent, FD7A 14 percent, FD8A 38 percent, FD9A 15 percent (mostly diaromatics and tricycloparaffins with no representative). The four eligible fuels are all light, high-PIONA, low-cetane fuels, so they are not representative of the nine. Primary lab is ORNL/SwRI (other labs are in the summary file):
* Density at 15.56 C (SG x 999.016, ideal volume): bias +5.8, MAE 15.0 kg/m3.
* ln viscosity at 40 C (Refutas): bias -0.41, MAE 0.41, 0 of 4 inside 0.186 ln (NCUT values: bias -0.31, 1 of 4). Under-prediction of 18-53 percent, which is much larger than the surrogate result, so the registry representatives (for example mesitylene for all C9 aromatics) and the light fractions of these fuels dominate the error, not Refutas itself. Inter-lab spread between ORNL and NCUT is up to about 0.12 ln.
* Flash (ideal Le Chatelier): bias -12.8 C, MAE 12.8 C, n = 4; all four predictions are below the measurement, so 0 of 4 over-predict by more than 5 C (unsafe direction not seen, but the sign is a model bias and not a safety margin).
* Cetane / DCN: not evaluable. The repo's only model is linear-by-volume blending of verified pure DCN values; verified values cover only 5-37 percent of the mass of the nine fuels (below the 90 percent coverage floor), so no prediction was made.

Sensitivity, all nine fuels with unmapped mass dropped and renormalised (not the declared reading; `summary_R-S6-2_sensitivity_all9.csv`): density bias -11.8, MAE 21.2 kg/m3 (n = 9, worst FD8A -61); viscosity bias -0.40 ln, 0 of 9 inside 0.186 ln; flash bias -14.5 C, MAE 14.5 C, 0 of 9 over-predicted by more than 5 C. n = 9 is only just above the n = 8 descriptive limit and the fuels are not independent (one FACE program, nine blends designed to span properties).

**R-S6-4 / R-S6-6.** The gate was met, so R-S6-4 does not apply. Nothing was written to `configs/`. Diesel design re-runs were not done (separate step). Caveats: diaromatics, benzodicycloalkanes, tricycloparaffins, fluorenes and aromatic sulfur have no registry representative and are dropped (reported per fuel in `mapping_summary.csv`); the 0.186 ln margin was derived from lab mixtures and is applied here to real fuels.
