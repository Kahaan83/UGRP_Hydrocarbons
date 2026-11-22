# freezing_v9: real-fuel freezing points and fuel compositions (v9 acquisition, 2026-11-23)

PDFs and source extracts staged in `extras/staging/lit/v9/` (gitignored, not committed).
Loader: `src/fuelblend/data/loaders_freezing_v9.py`; test: `tests/test_loaders_freezing_v9.py`.

## edwards2017_properties_draft.csv, edwards2017_composition_draft.csv
* Edwards, "Reference Jet Fuels for Combustion Testing", 55th AIAA Aerospace Sciences Meeting (AIAA 2017-0146), doi 10.2514/6.2017-0146. Also reported in AFRL-TR-2020-0017 (AD1093317).
* Numeric facts transcribed from Table A-3 (properties) and GC composition tables for the NJFCP reference fuels (A-1 POSF-10264, A-2 POSF-10325, A-3 POSF-10289, and POSF-4658).
* Licences: Government / public research data transcribed with citation; numeric facts are not copyrightable.
* Confidence: `clear` for all rows transcribed from Table A-3.
* Freezing points:
  - POSF-10264 (A-1 JP-8): -52.0 C (Table A-3, p14)
  - POSF-10289 (A-3 JP-5): -50.0 C (Table A-3, p14)
  - POSF-10325 (A-2 Jet A): -51.0 C (Table A-3, p14)
  - POSF-4658 (Jet A): -48.0 C (Table A-3, p14)
* Compositions:
  - 67 GC bins per fuel with class and carbon number, summing to 99.95-100.02 wt % (within 0.5 wt % closure).
  - Classes mapped to canonical schema: `n_paraffin`, `iso_paraffin`, `monocycloparaffin`, `dicycloparaffin`, `alkylbenzene`, `cycloaromatic`, `diaromatic`.
* Duplicate identity resolution (Task 6 rule):
  - POSF-10325 is canonical AJFTD fuel 2 (`fuel:ajftd:2`).
  - POSF-10289 is canonical AJFTD fuel 39 (`fuel:ajftd:39`).
  - POSF-4658 is canonical AJFTD fuel 23 (`fuel:ajftd:23`).
  - POSF-10264 is `fuel:afrl:posf10264` (new fuel not in AJFTD).
