# SLE primary-source candidates for Tm and dHfus (2026-11-14)

Deliverable: `candidate_rows.csv` (9 compounds). Nothing here changes `tm_dhfus_verified.csv`, any tier, `verification:` flag or margin; adopting values is the user's decision.

## What this is, honestly
All values come from the **NIST WebBook phase-change pages** (TRC / Domalski & Hearing 1996 compilations), one saved page per compound
(`nist_webbook_<name>/phase_change.html` + `SOURCES.md` with url, date, licence, sha256). Each row names the **original paper** the WebBook
cites, but I did **not** read those papers: Finke 1954, Scott 1962, Douslin & Huffman 1946, McCullough 1957, Pitzer & Scott 1943,
Taylor & Kilpatrick 1955 and Mears 1963 are paywalled or not open-access. The one open JPCRD route I tried (Acree and Chickos compilation) returned
403 from the publisher, so nothing from it is included. Same method and trust level as the existing `data/raw/sle_constants` file
(secondary compilation, adiabatic-calorimetry "DH" rows preferred). If primary-paper reading is required, that needs a library copy.

## Agreement with the registry (`pure_props_v5.csv`, source inverse_design_v2)
Registry values converted from degC to K.
| compound | Tm (K) WebBook vs registry | dHfus (kJ/mol) WebBook vs registry | verdict |
|---|---|---|---|
| toluene | 178.15 vs 178.15 | 6.636 vs 6.64 | agrees |
| n-hexane | 177.84 vs 177.85 | 13.079 vs 13.08 | agrees |
| n-octane | 216.38 vs 216.35 | 20.740 vs 20.73 | agrees |
| tetralin | 237.35 vs 237.35 | 12.447 vs 12.4 (registry "disputed") | agrees; supports 12.4 |
| o-xylene | 247.82 vs 247.95 | 13.598 vs 13.6 (registry "reference_only") | agrees (Tm 0.13 K, inside 0.25 K TRC unc) |
| mesitylene | 228.42 vs 228.45 | 9.5148 vs 9.51 | agrees |
| 2-methylpentane | 119.55 vs 119.45 | 6.268 vs 6.27 | agrees (Tm 0.10 K) |
| **hexylbenzene** | **209.7 vs 212.15** | no WebBook value (registry 17.5, "disputed") | **Tm differs by 2.45 K**; registry -61.0 degC matches no WebBook row (209.7, 209.7, 206.4 K) |
| n-nonane | 219.66 | 15.468 | **not in the registry** (no Tm/dHfus row in pure_props_v5.csv) |

## Caveats for whoever adopts these
- **n-nonane** has a solid-solid transition (crystal II -> I, 6.280 kJ/mol at 217.2 K, Finke 1954). 15.468 is fusion only from crystal I; Huffman 1931 (22.121) lumps both. Pick deliberately, as already done for the odd alkanes in `tm_dhfus_verified.csv` (`svl_dhfus_choice` columns).
- **Hexylbenzene dHfus** was not found in WebBook; it stays unsourced. Tm confidence is medium (TRC unc 0.3 K, but only 1963/1961/1929 sources).
- o-xylene has a tiny 0.031 kJ/mol II -> I transition at 208 K; negligible.
- Spread of dHfus across WebBook rows for toluene and n-hexane is ~0.1-0.7 kJ/mol; I took the adiabatic row with the best-documented Tm.
- Net: for 7 of 9 compounds the registry values are reproduced by a public compilation, so these would move from "old registry" to "WebBook-backed" without numerical change; hexylbenzene Tm is the one real discrepancy.

## Addendum: Acree & Chickos JPCRD 39, 043101 (2010) compendium (second source)
File: `candidate_rows_acree_chickos2010.csv`; source note `acree_chickos_jpcrd2010/SOURCES.md` (sha256 of the user's local PDF; PDF not copied into the repo).
**This is a compilation, not primary papers read.** Values were read off page images of Tables 6-9 (printed pages 043101-247, -248, -282, -337, -369, -398, -414, -437, -525). For fusion rows the compendium
overwhelmingly cites Domalski & Hearing 1996 (JPCRD 25, 1; key 1996DOM/HEA), the same compilation the WebBook uses, so it is **not independent** of the WebBook values; it only adds DSC rows (Mondieig et al. 2004 for octane and nonane) and confirms nothing beyond DOM/HEA.

| compound | compendium Tfus / fusH | vs WebBook | vs registry |
|---|---|---|---|
| toluene | 178 K / 6.61 (DOM/HEA, SOU/AND) | agrees | agrees |
| n-hexane | 177.8 / 13.08 | agrees | agrees |
| n-octane | 216.4 / 20.74 (DOM/HEA); 216.6 / 21.8 DSC | agrees | agrees (20.73) |
| tetralin | 237.4 / 12.45 | agrees | supports 12.4 |
| o-xylene | 247.8 / 13.6 | agrees | agrees (Tm +0.15 K vs registry) |
| mesitylene | 228.4 / 9.51 | agrees | agrees |
| 2-methylpentane | 119.6 / 6.27 | agrees | agrees |
| n-nonane | 219.7 / 15.48 fusion; separate row II->I 6.28 @217.2 K; DSC 219.5 / 15.0, trs 6.2 @218.2 | agrees | not in registry |
| hexylbenzene | **no fusion entry** (only vaporization 60.2 / 61.6 / 60.0 kJ/mol) | nothing to compare | registry Tm 212.15 K and dHfus 17.5 remain unsupported |

- **Hexylbenzene:** the compendium has no Tfus or dHfus, so the one remaining datum is the WebBook TRC Tm 209.7 K (Mears, Stanley et al. 1963; Jessup & Stanley 1961; Simon 1929 206.4 K). The registry's 212.15 K (-61.0 degC) and dHfus 17.5 are unsupported by both sources. A likely primary for Mears 1963 is J. Res. NBS 67A, 475 (doi 10.6028/jres.067a.049, "seven twelve-carbon hydrocarbons"); title match only, not confirmed to contain hexylbenzene.
- **n-nonane:** the compendium (like the WebBook) lists fusion and the solid-solid transition separately: fusion-only 15.48 kJ/mol at 219.7 K, crystal II->I 6.28 at 217.2 K. Total from below 217 K is 21.76 kJ/mol (DSC: 15.0 + 6.2 = 21.2). Which to use in SLE (fusion only vs total) stays your decision.

### Citation fixes (Crossref-checked 2026-11-14) in `candidate_rows.csv`
- tetralin: McCullough, Finke, Messerly, Todd, Kincheloe et al., J. Phys. Chem. 61, 1105 (1957), doi 10.1021/j150554a016 (not JACS 79, 4289).
- mesitylene: Taylor & Kilpatrick, J. Chem. Phys. 23, 1232 (1955), doi 10.1063/1.1742247 (not J. Res. NBS 55, 15).
- also confirmed and DOI added: Scott et al. 1962 (J. Phys. Chem. 66, 911), Douslin & Huffman 1946 (JACS 68, 1704), Pitzer & Scott 1943 (JACS 65, 803).
- The compendium's reference list does not carry the McCullough 1957 or Taylor & Kilpatrick 1955 papers for these rows (it cites DOM/HEA), so these two were confirmed by Crossref only, not by the compendium.
