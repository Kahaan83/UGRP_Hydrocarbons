# cold_visc_v9: real-fuel low-temperature kinematic viscosity (v9 acquisition, 2026-11-23)

Loader: `src/fuelblend/data/loaders_cold_visc_v9.py`; test: `tests/test_loaders_cold_visc_v9.py`.

## edwards2017_cold_visc_draft.csv
* Edwards, "Reference Jet Fuels for Combustion Testing", 55th AIAA Aerospace Sciences Meeting (AIAA 2017-0146), doi 10.2514/6.2017-0146.
* Table A-3 (p14): kinematic viscosity (mm^2/s = cSt) measured at -20 C (253.15 K) and -40 C (233.15 K).
* Fuels:
  - POSF-10264 (A-1 JP-8): -20 C: 3.5 cSt; -40 C: 6.6 cSt
  - POSF-10289 (A-3 JP-5): -20 C: 6.5 cSt; -40 C: 14.1 cSt
  - POSF-10325 (A-2 Jet A): -20 C: 4.5 cSt; -40 C: 9.2 cSt
  - POSF-4658 (Jet A): -20 C: 4.4 cSt; -40 C: 8.8 cSt
  - POSF-11498 (C-1 ATJ): -20 C: 5.0 cSt; -40 C: 9.89 cSt
* Confidence: `clear` for all rows.

## boehm2022_cold_visc_draft.csv
* Boehm et al., "Low-Temperature Viscosity of Aviation Turbine Fuels and Model Mixtures", Front. Energy Res. 2022, doi 10.3389/fenrg.2022.1074699. CC BY 4.0.
* Raw workbook: `data/raw/boehm2022/HEAT_LAB_ViscosityData_2022.xlsx`, sheet `This Work`.
* Measured kinematic viscosity (cSt) and density (g/cc -> kg/m3) at -20 C (253.15 K) and -40 C (233.15 K).
* Neat fuels:
  - A-1 (POSF 10264): -20 C: 3.4392 cSt; -40 C: 6.5151 cSt
  - A-2 (POSF 10325): -20 C: 4.4669 cSt; -40 C: 9.1359 cSt
  - A-3 (POSF 10289): -20 C: 6.0979 cSt; -40 C: 14.103 cSt
  - C-1 (POSF 13718): -20 C: 5.3575 cSt; -40 C: 10.577 cSt
  - C-7 (POSF 13727): -20 C: 6.2241 cSt; -40 C: 14.136 cSt
  - C-8 (POSF 13491): -20 C: 4.9402 cSt; -40 C: 10.512 cSt
  - C-3 (POSF 12959): -20 C: 9.7551 cSt; -40 C: 25.206 cSt
  - C-9 (POSF 12968): -20 C: 5.6587 cSt
  - HRJ Camelina (POSF 7720): -20 C: 6.0998 cSt
  - POSF 12945: -20 C: 5.0796 cSt
* Also includes measured blends with real fuels at -20 C and -40 C.
* Values stored in cSt. Dynamic viscosity conversion uses the simultaneously measured density at the same T.
* Confidence: `clear` for all rows read from the digital spreadsheet.
