# Real-Fuel Drafts: AFRL-TR-2020-0017 and Edwards (AIAA 2017-0146) (v10)

## Sources
1. **AFRL-TR-2020-0017 (DTIC AD1093317.pdf):**
   - The file provided under `data/raw/v10_pdfs/` is an INCAS Bulletin article (Prisacariu, Andrei, Circiu, 2025, *Jet engine combustor parametric analyzes*, INCAS Bulletin 17(3): 71-82).
   - Content inspection: The document presents 1D GasTurb parametric numerical analysis of combustor operation. Table 1 (page 2) lists average formula and molecular weight (MW) for reference fuels A-1, A-2, A-3, C-1 to C-6, citing Edwards [2]. Table 2 lists LHV and HHV for general fuel types (Hydrogen, Gasoline, Methanol, etc.). Tables 3 to 9 list GasTurb engine cycle simulation input parameters.
   - Status: No measured fuel properties (freezing point, cold viscosity, density, flash point) or GCxGC compositions are tabulated in this document. All tables are cataloged under target `other` in `results/v10/page_index.csv`.
   - Result: `afrl2020_properties_draft.csv` and `afrl2020_composition_draft.csv` contain header-only draft schemas with 0 extracted property/composition rows.

2. **AIAA 2017-0146 (Edwards 2017):**
   - Reference Jet Fuels for Combustion Testing (55th AIAA Aerospace Sciences Meeting).
   - Status: `not provided` by the user in `data/raw/v10_pdfs/` (recorded in `results/v10/pdf_register.csv`).
   - Result: Logged as `not provided` per unattended execution rules. Draft files are initialized with schemas.

## Licences & Notes
* All files strictly follow the project schema.
* No numbers were invented or estimated.
