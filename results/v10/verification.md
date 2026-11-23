# v10 verification (Task 5)

- `Vapor_Pressure_and_Liquid_Heat_Capacity.pdf` (J. Chem. Eng. Data 45(6), 2000, p. 1207), Table 2: all 25 perhydrophenanthrene (T, p) pairs read from the rendered page image and compared with the draft: 25/25 match (`verification_rows.csv`). Draft promoted to `data/raw/pure_v10/vapour_pressure.csv` (gitignored, not committed).
- `AFRL-TR-2020-0017 (DTIC AD1093317).pdf` is NOT the Edwards report. Page 1 is "Jet engine combustor parametric analyzes" (Prisacariu et al., INCAS Bulletin). It holds no fuel property tables. The user must re-download AD1093317.
- `je500381c_si_001.pdf` is adamantane-derivative + alcohol binary excess properties (not hydrocarbon fuels, not a candidate). Nothing transcribed.
- AIAA 2017-0146 not provided.
- Fuel drafts: 0 rows, so nothing to verify.

# Antoine refit (Task 6)
Perhydrophenanthrene: 25 points, 454.684 to 551.61 K, `failed:missing_flash_window`, `in_flash_window False` (Trap 3). No usable refit for flash. R10-4: 0 candidates with the P4 condition met.
