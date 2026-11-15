# Primary-paper calorimetry data (2026-11-15) — data only, nothing adopted

Transcribed by `scripts/v5/base/ingest_calorimetry_primary.py` from page images of four papers (the OCR text layer is unusable for tables).
The PDFs are ACS copyright and are not in the repo; numbers + citation only. Tests: `tests/test_calorimetry_primary.py`.

| key | paper | used |
|---|---|---|
| F54 | Finke, Gross, Waddington, Huffman, J. Am. Chem. Soc. 76, 333 (1954) — n-C8..n-C16 | Table I (triple pt, impurity), IV (dHfus), V (solid-solid transitions) |
| DH46 | Douslin, Huffman, J. Am. Chem. Soc. 68, 1704 (1946) — five hexanes | Tables I, IV, V, VI |
| PS43 | Pitzer, Scott, J. Am. Chem. Soc. 65, 803 (1943) — benzene/toluene/xylenes | Tables I, VI, VII, X; text (dHvap, eqs 1-5) |
| S62 | Scott et al., J. Phys. Chem. 66, 911 (1962) — toluene | Table V + text (Ttp, dHfus), Table VII + eq 2 (dHvap) |

## Files
* `phase_change_constants.csv` — Tm/triple point, impurity (mol %), dHfus, up to two solid-solid transitions, cal and kJ (x 4.184). One row per compound per paper; every compound has a single paper here, so `selected` is True throughout.
* `verification_vs_tm_dhfus_verified.csv` — each row compared with the existing `data/raw/sle_constants/tm_dhfus_verified.csv`.
* `vapor_pressure_pitzer_scott1943.csv` (+ `_equations_`) — toluene and o/m/p-xylene, mmHg and kPa; the paper's fits (1)-(5) reproduce its own table (test).
* `dhvap.csv` — dHvap at 298.16 K (PS43, 4 compounds) and toluene at 4 T/P points + eq (2) (S62).
* `xylene_eutectics_pitzer_scott1943.csv` — o/m/p-xylene binary and ternary eutectics, schema aligned with `data/processed/v5/thermoml_sle_hc.csv` (`system_id`, `inchikeys`, `x_liquid`, `T_K`, `crystal_inchikeys`, `point_kind`).

## Verification result (existing Finke/other rows were NOT edited)
* All 9 Finke fusion values, all 3 transition heats (C11, C13, C15) and every Tm agree with the registry to the last printed digit (cal x 4.184 reproduces 28.715, 22.179, 36.836, 28.501, 6.858, 7.661, ...). The convention in `SOURCES.md` is confirmed from the primary paper: for C9/C11 Table IV footnote c, the printed fusion heat is fusion only (transition + fusion were measured together and the transition subtracted), so `dhfus_total` = fusion + transition is a sum of the paper's own numbers.
* n-hexane 13.079 vs 13.080 (rounding). Toluene, o-xylene, 2,2-dimethylbutane agree.
* **2-methylpentane**: registry 6.268, paper 6.266 kJ/mol (1497.6 +/- 1.6 cal). Inside the paper's uncertainty, but not an exact match. Unchanged.
* 2,2-dimethylbutane has two solid-solid transitions (126.81 K 5.409 kJ/mol; 140.79 K 0.285 kJ/mol) and 2,3-dimethylbutane one (136.07 K 6.494 kJ/mol); the registry carries none for 2,2-DMB. Flag only (hexanes are outside the jet range).
* Toluene: Pitzer & Scott 1943 did not measure toluene fusion, so there is one primary source (S62), not two.

## Integration into pipeline files (2026-11-15, user-approved)
* `data/raw/sle_constants/tm_dhfus_verified.csv`: +4 rows (list `C`): n-nonane (`svl_dhfus_choice` = total 21.748, odd-alkane convention as C11/C13/C15), 2,3-dimethylbutane (fusion only; its 136.07 K transition is in the `ss` columns), m-xylene, p-xylene. 47 -> 51 rows; `build_tmap_verified()` 45 -> 49 entries (`tests/test_delta_freeze_margin.py` pin updated). **Frozen v5 results under `results/` were not regenerated**; reruns will pick the new constants up (any blend containing these compounds moves from the Coutinho table to verified constants).
* `data/unified/measurements.csv` / `measurements_dedup.csv`: +21 liquid vapour-pressure observations (source `calorimetry_primary`, ref `Pitzer&Scott1943 Table VII`); no collision with existing keys. The two solid p-xylene points are not included (liquid-only schema). The Antoine refit was not rerun; `refit_antoine.py` will use them next time.
* `data/processed/v5/thermoml_sle_hc.csv`: +3 measured eutectic rows (`point_kind` eutectic, file `calorimetry_primary:PS43_TableX`). All scoring scripts filter `liquidus`, so no current score changes.
* Those three files are gitignored/generated: after regenerating them, rerun `python scripts/v5/base/ingest_calorimetry_primary.py --integrate` (idempotent) to re-add the rows.
* Not integrated (no canonical property in the schema): heats of vaporization, impurities, transition heats beyond the verified table, calculated/Nakatsuchi eutectics.

## Caveats
* 3-methylpentane could not be crystallised (glass); Tm and dHfus are blank by design.
* m-xylene: Tm 225.27 K and 2765 cal/mol are for pure m-xylene after the paper's impurity correction (sample 2.34 mol % impurity; impure sample melted at 224.38 K, 2784 cal/mol).
* Temperatures are as published (0 degC = 273.16 K); not shifted. Vapour-pressure `T_K_published_basis` = t + 273.16.
* Xylene eutectics: composition basis is **assumed mole fraction** (the table does not say; the paper uses Raoult's law). The "calculated" rows are the paper's ideal-solution computation (DCp of fusion = 3R), not measurements; measured eutectic temperatures are 220.2 +/- 0.2 K (m-p, 12.9 % p), 211.4 +/- 1 K (o-m), 209 +/- 1 K (ternary); the o-p eutectic has no measurement in this work (237.6 K calculated, 237.7 +/- 0.4 K Nakatsuchi). Experimental compositions are deduced, not analysed. This is new validation data for the xylene system; it has not been scored.
* Vapour-pressure rows were not used to refit or alter any Antoine constants.

## Skipped (stated, not forgotten)
Raw and smoothed heat-capacity tables (F54 II/III, DH46 II/III, PS43 II-V, S62 IV/VI; Cp only matters for DCp corrections and the OCR/underline phase marks make the smoothed tables ambiguous), the melting-curve series (T vs fraction melted), entropies, ideal-gas functions, vibrational assignments, benzene (not a fuel component here), and the benzene/toluene/xylene disproportionation equilibria.
