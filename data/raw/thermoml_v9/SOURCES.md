# ThermoML v9 Archive (Post-2020 Update)

## Archive Information
* Source: NIST TRC (Thermodynamics Research Center) ThermoML Archive, public domain / U.S. Government open data.
* Base URL: `https://trc.nist.gov/ThermoML/<doi>.xml`
* Fetch Date: 2026-11-23
* Archive File: `thermoml_v9_archive.tgz` (gitignored, 475 XML files)
* Archive SHA256: `169fe5e4a3871b09ef9d91a2b2023bd494b3c5db3014bd6161799685823b34ce`
* License: NIST TRC / U.S. Government work, public domain worldwide.

## Journals Included
Articles published in 2020 across 5 thermodynamics and physical chemistry journals:
1. Journal of Chemical & Engineering Data (ACS): 176 files
2. The Journal of Chemical Thermodynamics (Elsevier): 185 files
3. Fluid Phase Equilibria (Elsevier): 58 files
4. International Journal of Thermophysics (Springer): 37 files
5. Thermochimica Acta (Elsevier): 19 files

None of these articles overlap with the pre-2020 NIST TRC snapshot (`ThermoML.v2020-09-30.tgz`) or the nine unlisted JCED files in `thermoml_2020`.

## Extraction Summary
* Total extracted hydrocarbon measurements: 994
* Pure hydrocarbon rows (n_comp = 1): 433
* Mixture hydrocarbon rows (n_comp >= 2): 561
* Cold temperature measurements (T < 263.15 K): 3
* SLE mixture rows: 0

### Counts by Property and Component Count
| Property | n_comp | Count |
|---|---|---|
| Kinematic viscosity, m2/s | 1 | 9 |
| Mass density, kg/m3 | 1 | 207 |
| Mass density, kg/m3 | 2 | 305 |
| Vapor or sublimation pressure, kPa | 1 | 180 |
| Vapor or sublimation pressure, kPa | 2 | 7 |
| Vapor or sublimation pressure, kPa | 3 | 18 |
| Viscosity, Pa*s | 1 | 37 |
| Viscosity, Pa*s | 2 | 231 |
