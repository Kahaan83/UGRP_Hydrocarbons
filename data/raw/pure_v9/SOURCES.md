# Pure Constants for Candidate Compounds (v9 Target 3)

## Overview
Tabulated and measured physical constants for the seven candidate design compounds identified in `results/v5/BASE/delta_representatives/candidate_gaps.csv`:
1. `1-methylnaphthalene` (`Cc1cccc2ccccc12`, CAS 90-12-0)
2. `2,6-dimethylnaphthalene` (`Cc1ccc2cc(C)ccc2c1`, CAS 581-42-0)
3. `perhydrophenanthrene` (`C1CCC2C(C1)CCC1CCCCC12`, CAS 5743-97-5)
4. `perhydrofluorene` (`C1CCC2C(C1)CC1CCCCC12`, CAS 5744-01-4)
5. `1-methyltetralin` (`CC1CCCc2ccccc21`, CAS 1559-81-5)
6. `n-butyltetralin` (`CCCCC1CCCc2ccccc21`, CAS 28729-54-6)
7. `2-methyldecalin` (`CC1CCC2CCCCC2C1`, CAS 2958-76-3)

## Data Sources
1. **Vapour Pressure:**
   - NIST Chemistry WebBook (SRD 69), Camin & Rossini (1955), Osborn & Douslin (1975), Rohác et al. (2000), Sasse et al. (1988). WebBook Antoine constants in bar: $A_{\text{kPa}} = A_{\text{bar}} + 2$.
   - ThermoML / J. Chem. Eng. Data / J. Chem. Thermodyn. measured liquid vapor pressures.
2. **Lower Flammability Limit (LFL):**
   - IS IEC 60079-20-1 (2010), "Explosive atmospheres - Part 20-1: Material characteristics for gas and vapour classification".
   - US Bureau of Mines Bulletin 627 (Zabetakis, 1965), "Flammability characteristics of combustible gases and vapors".
   - CRC Handbook of Chemistry and Physics, 95th ed. (2014-2015), Section 16: Health and Safety Information.
   - GESTIS Substance Database (IFA).
   - Only tabulated or measured values are admitted (`lfl_source in {"measured", "tabulated"}`); Jones-rule estimates ($0.55 \cdot C_{\text{st}}$) are refused.
3. **Viscosity & Density:**
   - ThermoML archives: `10.1021/je101256z` (1-methylnaphthalene), `10.1021/je500381c` (perhydrofluorene).
   - Boehm et al. (2022) / HEAT Lab pure measurements.
   - Seyer et al. (1942), J. Am. Chem. Soc. 64: 2136-2139 (methyldecalin isomers).

## Licences & Attributions
* NIST Chemistry WebBook: U.S. Secretary of Commerce on behalf of the United States of America. Public domain / open standard reference data.
* US Bureau of Mines Bulletin 627: U.S. Government work, public domain.
* IEC 60079-20-1: Standardized reference property tables.
* ThermoML: NIST TRC open public domain archive.

## Transcription Verification
All transcriptions are written to `*_draft.csv` files with explicit `page`, `table`, and `confidence` columns. Every value was cross-referenced with primary literature references.
