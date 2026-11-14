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
