# Pure Physical Constants Drafts (v10)

## Candidate Compounds Covered
Design candidates targeted for physical constants:
1. `perhydrophenanthrene`: 25 primary vapor pressure points extracted from `Vapor_Pressure_and_Liquid_Heat_Capacity.pdf` Table 2 (J. Chem. Eng. Data 2000, 45, 1205-1210).
   - T range: 454.684 K to 551.610 K (181.534 to 278.460 degC)
   - P range: 8.292 kPa to 99.533 kPa
   - Method: Comparative ebulliometry with water reference
   - Extraction: Direct text layer (200 dpi rendered page verified)
2. Other candidate compounds (1-methylnaphthalene, 2,6-dimethylnaphthalene, perhydrofluorene, 1-methyltetralin, n-butyltetralin, 2-methyldecalin):
   - No primary viscosity, density, or LFL data present in the registered PDFs.
   - `je500381c_si_001.pdf` contains only binary alcohol-adamantane excess properties (excluded).
   - `AFRL-TR-2020-0017 (DTIC AD1093317).pdf` contains no pure property tables.

## Rules & Quality Controls
- No estimated LFL values (Jones rule, homologue extrapolations) are admitted.
- Unit conversions: pressures in mmHg, bar, or Pa are normalized to kPa in loaders.
- Any scanned image-extracted rows are forced to `unsure` confidence.
- InChIKeys are referenced from `data/unified/compounds.csv` or left blank if not registered.
