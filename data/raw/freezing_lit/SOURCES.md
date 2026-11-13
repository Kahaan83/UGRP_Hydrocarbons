# freezing_lit: literature freezing points and fuel compositions (v5 paper ingest, 2026-11-13)

PDFs are staged in `extras/staging/lit/` (gitignored, not committed). Loader: `src/fuelblend/data/loaders_freezing_lit.py`;
test: `tests/test_loaders_freezing_lit.py`. Common columns: see `.superpowers/sdd/2026-11-12-v5-api-route-and-design/ingest-common.md`.

## vozka2018_gcxgc.csv, vozka2018_properties.csv
* Vozka, Simacek, Kilaz, Energy & Fuels 2018, supporting information for the HEFA/Jet A blend paper, doi 10.1021/acs.energyfuels.8b02787
  (doi taken from the batch brief; the SI itself prints none). ACS SI, licence not stated; numeric facts transcribed with citation.
* Table S2 (SI pp2-6): hydrocarbon type composition, wt %, class x carbon number, for Jet A, CAME, TALL, MFAT. 272 bins, zero bins kept as printed.
  Transcription: the SI is a text PDF, so values were read from the embedded text layer by a script (class headers + `C<n>` rows), not retyped.
  Closure check against the printed class totals: all match to 0.01 except Jet A isoparaffins (27.44 vs 27.46), MFAT isoparaffins (89.09 vs 89.11),
  Jet A monocycloparaffins (23.36 vs 23.39), Jet A cycloaromatics (7.11 vs 7.13); these differences are the paper's own rounding (13 bins each).
  Fuel totals as printed bins: Jet A 99.95, CAME 100.03, TALL 100.01, MFAT 99.98 (paper's total row says 100.00).
* Table S4 (SI pp9-10): measured freezing points (C) of Jet A, three neat HEFAs and six blends each. Minus signs: the text layer of this copy was
  NOT garbled (the survey expected '7' for '-'), but every sign and digit of the measured rows was checked against the rendered page images
  (p9 and p10). Cookson-equation rows are model values and are NOT ingested.
* Jet A -51.0 appears at the head of all three blocks (identical): entered once. Blends are `kind=real_fuel`, `components` = [Jet A, neat HEFA],
  `fractions` = [1 - f, f] with f = 0.1 ... 0.6 read from the labels C-10..C-60, T-10..., M-10....
* JUDGEMENT: blend basis is NOT confirmed. The SI never says vol % or mass %; the only vol % in the SI is the D86 definition in the Table S4
  footnote. The main paper (not staged) would say. Rows are tagged `fraction_basis=volume` as a working reading and `qualifier=blend_basis_unconfirmed`.
* Densities: Table S3 note (SI p9) quotes measured 0.7598, 0.7573, 0.7612 g/cm3 for CAME, TALL, MFAT, temperature not stated (`qualifier=T_unspecified`,
  stored x1000 as kg/m3). Table S3/S5 per-bin contributions are calculated, not ingested.
* Class mapping (`CLASS_MAP` in the loader): n-paraffins -> n_paraffin; isoparaffins -> iso_paraffin; monocycloparaffins -> monocycloparaffin / cycloalkane /
  cyclo_mono; di- and tricycloparaffins -> dicycloparaffin / cycloalkane / cyclo_di (the paper merges di and tri; tricyclics are 0 or negligible, AJFTD has
  them separate); alkylbenzenes -> alkylbenzene / aromatic; cycloaromatics -> cycloaromatic / aromatic / indan_tetralin; alkylnaphthalenes -> diaromatic /
  aromatic / aromatic_di (the AJFTD parser files naphthalenes under diaromatic). Same vocabulary as `data/real_fuels/v5/ajftd_gcxgc_composition.csv`.
  `vozka_as_ajftd_composition()` emits that layout; `blend_composition(sample)` mixes the neat tables linearly (exact for a mass blend, a small
  approximation for a volume blend).

## ascent065A_freezing.csv
* F2: Bell, Boehm, Heyne, "Freezing Point of Hydrocarbon Fuels from Single Species Concentrations", Energy Fuels 39, 4221 (2025), doi 10.1021/acs.energyfuels.4c06091,
  reproduced as Appendix A of the FAA ASCENT Project 065(A) annual report 2025 (WSU). No licence stated.
  Table A-1 (PDF p6): freezing points of 10 neat compounds (checked on the rendered page). Table A-2 + prose (PDF p10): mole % recipe of a 100 % n-alkane fuel
  (nC9-nC16, sums to 99.99, renormalised by the loader) and its measured freezing point -40.7 C (prose; model value -44.1 C not ingested).
* F3: ASCENT 065(A) slides, WSU, 2024 (no doi). Slide 21 (rendered, checked): six mixtures all containing 0.17 % nC16 with freezing points; slide 16: nC16 16.8 C,
  nC13 -4.2 C. JUDGEMENT: slide percentages do not state a basis; read as mole % because F2 states mole concentrations and uses mole-fraction axes.
  The slide nC16 (16.8) and nC13 (-4.2) pure values differ from F2 Table A-1 (18.1, -5.3); both kept, F3 rows `qualifier=replicate`. nC8, iso-octane, bicyclohexyl are
  identical in F3 and F2 and entered once. Slide 12 binary -18.4 C has no stated composition: not ingested. F2 Figs A-1, A-2, A-4, A-5 are plots: skipped.
* SMILES are written by hand (n-alkanes, 2,2,4-trimethylpentane, bicyclohexyl); `kind` pure / binary / multi.

## kittel2023_fuels.csv
* Kittel, Horsky, Simacek, Processes 11, 935 (2023), doi 10.3390/pr11030935, CC BY 4.0. Table 1 (PDF p4), seven neat samples, checked on the rendered page.
* Methods from the text (PDF p3): flash ASTM D56(A), aromatics ASTM D6379(A) (HPLC), smoke ASTM D1322, freezing ASTM D5972 (FPA-70X, range +20 to -80 C),
  net specific energy ASTM D4809 (measured value only), density at 15 C, distillation ASTM D86.
* Censoring: freezing '<-80' for Jet HC and FCC HN (`qualifier='<'`); PyrTIR is printed '-80' without '<' (kept as -80, note says it is at the instrument limit);
  smoke '>50' for HEFA Cam and HEFA 215 (`qualifier='>'`, sample did not smoke); FCC HN smoke point is printed '-' (not reported, no row);
  flash for PyrTIR and PyrPO printed as 'not measured' (no row); PyrPO residue and loss not measured (no row); PyrPO D86 T10 and end point are calculated from
  SimDist (`qualifier=from_simdist`).
* `distillation_T100_C` = 'end of distillation' of D86 (final boiling point). Aromatics (total, mono, di) in vol %. FCC HN's '9.3' is a naphthalene entry under a
  separate row (not ingested).
* Class composition: only aromatics vol % (mono/di). n-alkane contents appear in figures only (skipped); no GCxGC tables, so these seven are bulk-property
  validation rows with `is_hc` set True by assumption (petroleum or HEFA / pyrolysis streams; no molecule list).
