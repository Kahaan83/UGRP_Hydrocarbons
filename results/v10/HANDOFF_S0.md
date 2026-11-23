# v10 Phase S0 handoff

- v9 closed: e46bb3c is an ancestor of origin/feat/external-datasets-ingest.
- Register (`results/v10/pdf_register.csv`): 3 files ok (all start `%PDF-`, all have a text layer), 1 not provided.
  - `AFRL-TR-2020-0017 (DTIC AD1093317).pdf`: 10 pages (short for a full technical report; Gemini must check coverage in the page index).
  - `je500381c_si_001.pdf`: 23 pages.
  - `Vapor_Pressure_and_Liquid_Heat_Capacity.pdf`: 6 pages.
  - AIAA 2017-0146: not provided.
- `sources.txt` was not created, so titles are blank.
- Tests: `pytest tests/v10` 5 passed.
- Next: phase G (Gemini), Tasks 1 to 4.
