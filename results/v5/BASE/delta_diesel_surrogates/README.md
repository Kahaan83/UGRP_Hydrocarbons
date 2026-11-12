# delta_diesel_surrogates: Mueller et al. 2016 diesel surrogates through predict()

Why: a surrogate is a recipe of pure compounds with measured properties, so it is the cleanest real test of the any-size API on a known composition. The diesel cell of the scoreboard is `no_evaluation` today.

## Source (numbers transcribed by hand; no PDF, no quoted text in the repo)

Mueller, Cannella, Bays, Bruno, DeFabio, Dettman, Gieleciak, Huber, Kweon, McConnell, Pitz, Ratcliff, Energy & Fuels 30(2):1445-1461 (2016), doi 10.1021/acs.energyfuels.5b02879 (open access, ACS AuthorChoice; page 1 checked against the DOI). Only numbers are kept, each with table (or figure) and page.

* Recipes: Table 3, page 1449, mole % columns (mass % kept only in `tests/test_delta_diesel_surrogates.py` as a consistency check). Four surrogates of one target diesel (CF): V0a (4 components), V0b (5), V1 (8), V2 (9). Recipes sum to 100.0 mol %. Appended to `data/raw/surrogates/fuel_surrogates.csv` as `Mueller2016_<name>`, frac_type `mole`.
* Measured properties: `data/raw/surrogates/mueller2016_properties.csv`. Flash point (ASTM D93, Table 5, page 1457) for CF and the four surrogates; density at 20 C (ASTM D4052; values read from the bar labels of Figure 11, page 1455; CF has two measurements from two laboratories); cloud point (ASTM D5773 / NREL procedure; bar labels of Figure 14, page 1456). Only measured values; the paper's predicted densities are not transcribed.
* Palette compounds: 13 in total. Twelve have an InChIKey in `data/unified/compounds.csv`. **Missing: perhydrophenanthrene (PHP), used only in V2 (6.0 mol %).** Nothing is guessed for it; the API's own substitution rules (representatives, unmapped mass) decide what happens, and V2 is reported with its coverage.
* Pre-existing issue, not changed here: three older rows of `fuel_surrogates.csv` (KimVioli2014_UM1, UM2, Chen2016_JI) give heptamethylnonane a SMILES with 15 carbons; the correct compound has 16. The Mueller rows use the 16-carbon SMILES.

## Declared readings (committed before any prediction is made)

The cell is descriptive and `n_too_small` (n = 4 surrogates, one target fuel). No skill claim is made, and no claim of passing or failing a threshold.

Each surrogate is run through `predict_fuel(..., fuel_type='diesel', composition_kind='known')` at the measurement temperature (flash: no temperature; density: 20 C; cloud point: no temperature).

* **Flash point.** Per-surrogate error (predicted minus measured, C), designated route (overlay, as the API ships it). Also reported: the one-sided over-prediction direction, and the lab-margin warning that the API attaches.
* **Density.** Per-surrogate relative error at 20 C against the measured value (the CF target is compared with the mean of its two measurements).
* **Cloud point.** The API's designated freezing model (ideal SvL in mole fraction) against the measured cloud point, labelled "SvL vs cloud point, different property, indicative only". Cloud point is not freezing point, and no statement about the freezing model's skill follows.
* **Coverage dict.** For every surrogate and property: status, substituted mass, unmapped mass, and whether each component had own data. V2 is expected to show the effect of the missing PHP; this is reported, not corrected.
* **Target comparison.** A second table compares each surrogate's predicted values with the measured values of the target diesel CF, for the same three properties (this says how close the surrogate-as-recipe is to the target, not how accurate the model is).
* Nothing is selected, tuned or excluded after seeing the results. If a surrogate returns `None` for a property, it is listed as such with its status.
* The API is used as it is in the commit this task starts from; if `predict_route.py` changes before the run, the commit hash used is recorded in the Result section.
