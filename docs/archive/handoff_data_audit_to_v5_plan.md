# Handoff: data audit results for the v5 plan (`i-want-to-move-joyful-platypus.md`)

Written 2026-11-11 by the data-audit session, for the agent executing the v5 plan. Full detail: `docs/data_audit_2026-11-11.md`.
Everything below is on branch `v5`, in the working tree, **uncommitted**. Full suite after these changes: 383 passed, 0 failed.

## 1. What changes the plan

| plan item | what the audit found | where |
|---|---|---|
| D1 ("inspect the Composition Metrics rows for hydrocarbon-type data") | The answer is in a different file. `extras/fuelblend_data/data/raw/ajftd/ajftd_fuel_details.csv` holds AFRL GCxGC weight % by class and carbon number for 116 fuels; 113 close to 98-102 %. Labels are POSF numbers. `docs/ajftd_v5_note.md` said no composition exists; it now carries a correction. | loader `src/fuelblend/data/loaders_ajftd_gcxgc.py`; `data/real_fuels/v5/ajftd_gcxgc_{composition,fuels,mapped}.csv` |
| Tier T2 size | T2 was about 42 fuels (FuelLib 4 + AV-33-22 38). 39 more AJFTD fuels now have measured composition and at least one measured target (flash 37, freezing 31 scored, viscosity 32, density 30). | `results/v5/G1/` |
| S3 (bin to representative-compound mapping, needed for AV-33-22) | A mapping function exists and is tested: `map_fuel()` in `scripts/v5/eval_ajftd_gcxgc_v5.py`, built on `reconstruction_r3.SubReps`. It takes rows of (class, carbon number, mass %). AV-33-22's 104 bins need only a class-name map. No diaromatic or tricycloparaffin representative exists in the registry, so those bins are unmapped mass (at most 5 % on AJFTD). | `tests/test_v5_data_audit.py` |
| S2 part C (flash on measured GCxGC) | `results/v5/S2/README.md` is still declared-but-unscored. An independent answer now exists on AJFTD: ideal flash bias -0.43 C, MAE 3.54 C, 11 % over by more than 5 C (n = 37); UNIFAC is no better. S2's own readings on FuelLib and AV-33-22 are still owed. | `results/v5/G1/RESULTS.md` |
| S2 (pure inputs) | Besides Antoine ranges, LFL is an input issue: 82 of 98 pure inputs in the v4 flash run use the Jones estimate. An overlay with tabulated values (IEC table, Co-Optima) was built and tested in memory: MAE moves under 0.1 C, CRC 647 dev over-predictions fall from 6 to 2 of 45. Declared reading "safe to adopt" holds. **Adoption as default is waiting on the user; do not regenerate v4 flash numbers or margins from it until they decide.** | `data/processed/v5/lfl_v5.csv`, `results/v5/L1/` |
| S5 (freezing) | Ideal Schroeder-van Laar on measured GCxGC has no skill: MAE 15.2 C against 7.3 C for a leave-one-out constant (n = 31). The AJFTD GCxGC table gives n-alkane mass by carbon number per fuel, which is the input the Bell-Boehm-Heyne first-species test needs. | `results/v5/G1/` |
| S5 / design certification (freezing margin `inf`, 8 systems) | 51 of 59 ThermoML SLE systems were skipped only for lack of audited fusion constants. With wider constants: 11 systems give a 12.3 K margin, 38 systems give 24.3 K. Finite, but far wider than the 0.5-0.7 C by which the D4 designs clear the limit. `configs/` untouched. | `scripts/v5/eval_sle_extended_v5.py`, `results/v5/SLE2/` |
| S1 (cold viscosity) | On real fuels with measured composition Refutas is fine at 40 C (bias -0.2 %) and under-predicts cold: -26 % at -20 C (n = 32), -49 % at -40 C (n = 12). | `results/v5/G1/` |
| S3 density question | Ideal-volume density is still biased low with measured composition: -6.1 kg/m3 at 15 C (n = 30). This used the pure-property provider, not the v4 ThermoML tier model, so compare with care. | `results/v5/G1/` |
| D9 / S6 (diesel, "recalled, not verified") | Now partly verified. NIST TN 1949 (doi 10.6028/NIST.TN.1949) Table 1 gives exact mole fractions of the four CRC AVFL-18a surrogates V0a, V0b, V1, V2; Table 11 gives measured 40 C viscosity. CRC FACE-1 (2010) has GC group-type, PIONA and DHA tables for the nine FACE diesels; whether composition is per carbon number is not yet checked. Neither PDF is saved in the repo (downloads need the user's go-ahead). Mueller et al. 2016 (the measured-property source for the surrogates) is not opened. | `docs/data_audit_2026-11-11.md` section 5 |

## 2. Constraints to respect

* **Locked AJFTD slice.** G1 read 11 fuels with S4 role `locked` (nothing was fitted). The plan says locked sets are scored once: treat these 11 as already read for composition-based flash, viscosity, density and freezing.
* **Double counting.** AJFTD labels 4658, 10289, 10325 are also in FuelLib; 10289, 10325, 11769 are in CRC AV-27-18; 10325 is the v4 headline fuel. They are flagged `shared_with_other_source` in G1 and need one shared group before any pooled T2 statistic. Overlap between AJFTD GCxGC fuels and the 38 AV-33-22 fuels is not checked.
* **Readings.** `results/v5/G1/README.md` and `results/v5/L1/README.md` were written before scoring; do not edit them, add dated notes.
* **AJFTD is login-gated.** The new composition CSVs sit beside the existing AJFTD v5 tables; raw rows must not be published.
* Four AJFTD fuels over-predict flash by 8 to 17 C (labels 7695, 9324, 8092, 8438). Unexplained; do not drop them.

## 3. Unused data checked and ready to wire (not done)

* NISTIR 6659: 121 compressed-liquid density rows for POSF 4658, which has a composition (FuelLib, AJFTD) and no property rows (`extras/.../processed/nist_density_compressed_liquid.csv`).
* NISTIR 6646 RP-1: 1,465 property rows from the same sample as the composition; v4 treated them as different fuels. Composition lists only 80.6 wt %.
* `data/processed/v2/wax_dauphin_coutinho.csv`: 57 measured wax-appearance temperatures, never read by the SLE evaluation.
* 74 AJFTD fuels have a complete composition and no property row in the AJFTD bulk table; POSF labels may link to other property sources.

## 4. Still owed from this work

* `docs/data_sources.csv`: no rows yet for the AJFTD GCxGC table, the LFL overlay, or the extended SLE constants.
* `docs/claims_register.md`: no claims added for G1, L1 or SLE2.
* None of the new scripts is in `reproduce_all.py`.
* Phase 0 item still open: `extras/ajftd_session.json` and `extras/fuelblend_data/ajftd_session.json` (login cookies) are still on disk.

## 5. New files

`src/fuelblend/data/loaders_ajftd_gcxgc.py`; `scripts/v5/{eval_ajftd_gcxgc_v5, eval_flash_lfl_overlay_v5, eval_sle_extended_v5, build_lfl_overlay_v5}.py`;
`tests/test_v5_data_audit.py`; `results/v5/{G1, L1, SLE2}/`; `data/real_fuels/v5/ajftd_gcxgc_*`; `data/processed/v5/lfl_v5.csv`;
`docs/data_audit_2026-11-11.md`; edit to `docs/ajftd_v5_note.md`.
