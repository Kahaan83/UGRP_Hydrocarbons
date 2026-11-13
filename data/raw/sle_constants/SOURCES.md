# Primary-source Tm and dHfus (data only)

File: `tm_dhfus_verified.csv`, 39 compounds. `list` = A (8 design-pool compounds) or B (the other 31: components of the ThermoML liquidus systems that
were unscored in `results/v5/BASE/delta_newdata_freeze_smoke/rows_freezing.csv` or whose system is jet-range, all components C8-C16, in
`data/processed/v5/thermoml_sle_hc.csv`; plus iso-octane from the ASCENT unscored list). Nothing is estimated: a missing value is blank and the note says so.

## Method
1. NIST WebBook phase-change pages (`cbook.cgi?InChI=<InChIKey>&Units=SI&Mask=4`, fetched 2026-11-13, one page per compound) were parsed for
   every `T triple`, `T fus`, `Delta fus H` and `Delta H trs` row with its original reference. Values in the CSV are copied from those rows; the
   citation text and DOI (where WebBook has one) come from the page's reference list. The original paper was not re-read.
2. Choice rule. dHfus = the adiabatic-calorimetry row (WebBook comment `DH`, crystal I -> liquid) of the best-documented study; Tm = the temperature
   printed with that row, which agrees with the TRC triple point to the digits shown. Where only a compilation or DSC row exists (1,4-/1,5-/1,6-DMN:
   Cheon and Kim 2007 DSC) that is used and flagged DSC. WebBook's averaged `T fus` rows are not used as Tm (they are not a single measurement).
3. Cross-checks (column `crosscheck`, not used for the values): pure-component rows in the two ThermoML archives used by
   `src/fuelblend/data/thermoml_sle.py` ("Normal melting temperature", "Triple point temperature", "Molar enthalpy of transition or fusion";
   mostly DSC from mixture papers, so noisier), and ASCENT 065A Table A-1 pure freezing points
   (`data/raw/freezing_lit/ascent065A_freezing.csv`, Tm only; n-C11..C16 agree with the chosen Tm to 0.1 K, iso-octane to 0.6 K).
4. `coutinho_*` columns come from calling `fuelblend.thermo.coutinho_wax.get_pure_thermal_properties` on a frame built exactly as
   `fuelblend.eval.composition_route.build_tmap` does (smiles, `hc_class`, carbon count). No code was edited. Blank = compound absent from that table.
5. `agrees` = |Tm - Tm_coutinho| <= 1.0 K and |dH - dH_coutinho| <= 5 % where dH is `svl_dhfus_kJ_mol` (see below). Blank when Coutinho has no entry.

## Fusion-only vs total (solid-solid) enthalpy
Odd n-alkanes (C11, C13, C15) and the long alkanes (C26, C36) have a rotator (crystal II -> I) transition below the melt. Columns:
`dhfus_kJ_mol` is always fusion only (I -> liquid); `dhfus_ss_transition_kJ_mol`, `ss_transition_T_K` the transition; `dhfus_total_incl_ss_kJ_mol` the sum.
`svl_dhfus_choice` / `svl_dhfus_kJ_mol` is what I recommend the Schroeder-van Laar term use.
Evidence for "total": the code's table formula dH = 4*nC - 10.6 kJ/mol gives 93.4 (C26) and 133.4 (C36), matching the measured totals 94.12 and 129.29,
not the fusion-only 60.70 and 88.83. For even alkanes (no transition near the melt) fusion = total and the formula matches (C12 37.4 vs 36.84, C14 45.4 vs 45.07,
C16 53.4 vs 53.36). The formula does not separate odd and even chains, so for odd chains it is 12-15 % above the measured total
(C11 33.4 vs 29.04, C13 41.4 vs 36.16, C15 49.4 vs 43.76) and 50-100 % above fusion only.
Caveat: I inferred the convention from the code's numbers; I did not read Coutinho's papers in this session. A user decision is needed.
Physically, between the transition temperature and Tm only the rotator phase is the equilibrium solid and dHfus alone applies; below the transition the
total applies. For the non-alkanes with a transition close to the melt (1,1'-bicyclohexyl 3.3 K, 1-methylnaphthalene 2.0 K, 2-methylnaphthalene 19 K) I kept
fusion only and flagged it; the totals are in the CSV.

## Judgement calls
* **Decalin.** WebBook cross-labels the McCullough et al. 1957 enthalpy rows: the stereo-unspecified page shows 14.414 kJ/mol at 242.78 K, the trans page shows
  9.489 kJ/mol at 230.18 K. The trans triple point on the trans page is 242.77 K and the cis average is 230.15 K (ThermoML
  paper 10.1016/j.fluid.2007.01.029 gives trans 242.71, cis 228.07), so I assigned 14.414 to trans and 9.489 to cis by temperature match. Parks and Hatton 1949 (WebBook prints 3.244 and 2.209 as
  kJ/mol, but the values are kcal/mol: 13.57 trans, 9.24 cis) support it within 6 % (trans) and 3 % (cis). The stereo-unspecified decahydronaphthalene has no entry (an isomer mixture has no single Tm).
* **1,8-dimethylnaphthalene.** WebBook's Finke 1977 enthalpy row shows 0.04687, which equals the entropy of fusion in kJ/mol/K (x 336.33 K = 15.76). I used the Acree 1991
  compilation value 15.77 kJ/mol at 336.3 K (it compiles that calorimetry). Cheon and Kim DSC gives 18.53 (17 % higher).
* **Hexadecane.** dHfus values 47.5-53.4 kJ/mol in the literature; Finke 1954 (53.359) kept, consistent with C10-C14 from the same study.
* **Hexatriacontane.** WebBook "Company, 1973" lists 130.666 for fusion, which is inconsistent with the Schaerer 1955 split (88.83 + 9.92 + 30.54); Schaerer used.
* **Anthracene, 2,3-DMN, 1,8-DMN.** Where DSC values spread from the adiabatic value, the calorimetric value is kept and the spread is in the note.
* **1-hexyne.** Tm only; no enthalpy of fusion in WebBook or ThermoML (`no primary source found`).
* **Not covered.** CRC/DIPPR tables were not used (not primary, not accessible). DTIC/AFRL was not consulted. ASCENT 065A is Tm-only, used as a check.

## Defect found in the existing Coutinho table (not fixed; code is out of scope)
`get_pure_thermal_properties` identifies aromatics by SMILES substring: `"c1ccccc1" in smi` matches fluorene and phenanthrene, so both get benzene's
Tm 278.7 K / dHfus 9.87 kJ/mol (true 387.94 K / 19.58 and 372.38 K / 16.46), and `"c1ccc2ccccc2c1"` matches 2-methylnaphthalene, which gets naphthalene's 353.4 K / 19.0
(true 307.73 K / 12.13). Cycloalkanes, 1,2,4-trimethylbenzene, methyl- and dimethylnaphthalenes (except via the substring above), pyrene, anthracene, acenaphthene,
biphenyl have no entry.
