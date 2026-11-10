# r3 CRC 647 reconstruction report

Inputs per fuel: CRC 647 Table G.1 (ASTM D6379 aromatics and saturates in vol%, ASTM D2425 hydrocarbon-type shares), the D86 curve (mean of tests 1 and 2: IBP, T10, T20, T50, T90, FBP) and olefins vol% (D1319, else D1159). No target property is read.

CRC fuels written: 57 (of 57); skipped: none

## n / iso paraffin split (prior)

n-paraffin share of the paraffins = 0.4000 (n / (n + iso)). This is a prior from the same 4 GC jet fuels r2 used, not a fit to CRC; Table G.1 does not separate normal from branched paraffins.

## Mass share per sub-kind (over the fuels written)

| sub-kind | mean | min | max |
|---|---|---|---|
| n_paraffin | 0.2281 | 0.1458 | 0.3491 |
| iso_paraffin | 0.3421 | 0.2188 | 0.5237 |
| cyclo_mono | 0.1036 | 0.0000 | 0.2402 |
| cyclo_di | 0.0851 | 0.0000 | 0.2972 |
| alkylbenzene | 0.1447 | 0.0329 | 0.2029 |
| indan_tetralin | 0.0524 | 0.0119 | 0.1103 |
| aromatic_di | 0.0000 | 0.0000 | 0.0000 |
| alkene_alkyne | 0.0090 | 0.0042 | 0.0278 |
| fallback (total) | 0.0349 | 0.0000 | 0.0817 |

## Fallback

Fallback mass: mean 0.0349, max 0.0817 (worst fuel CRC647_40).

Mean fallback mass by bucket (fuels without that bucket count as 0):

| bucket | mean mass |
|---|---|
| no_rep_alkene_alkyne | 0.0071 |
| no_rep_aromatic_di | 0.0175 |
| unexplained_minor_aromatics | 0.0014 |
| unexplained_tricycloparaffins | 0.0089 |

Per-fuel detail: results/v4/R3/fallback_detail.csv.

## Representative per (sub-kind, carbon number) bin

Candidates require a usable registry entry, valid pure eta/rho, and an Antoine set plus LFL in the flash registry. Lowest InChIKey within a bin; `n cand` > 1 marks bins that hold more than one candidate.

| sub-kind | C# | n cand | inchikey | smiles |
|---|---|---|---|---|
| n_paraffin | 3 | 1 | ATUOYWHBWRKTHZ-UHFFFAOYSA-N | CCC |
| n_paraffin | 4 | 1 | IJDNQMDRQITEOD-UHFFFAOYSA-N | CCCC |
| n_paraffin | 5 | 1 | OFBQJSOFQDEBGM-UHFFFAOYSA-N | CCCCC |
| n_paraffin | 6 | 1 | VLKZOEOYAKHREP-UHFFFAOYSA-N | CCCCCC |
| n_paraffin | 7 | 1 | IMNFDUFMRHMDMM-UHFFFAOYSA-N | CCCCCCC |
| n_paraffin | 8 | 1 | TVMXDCGIABBOFY-UHFFFAOYSA-N | CCCCCCCC |
| n_paraffin | 9 | 1 | BKIMMITUMNQMOS-UHFFFAOYSA-N | CCCCCCCCC |
| n_paraffin | 10 | 1 | DIOQZVSQGTUSAI-UHFFFAOYSA-N | CCCCCCCCCC |
| n_paraffin | 11 | 1 | RSJKGSCJYJTIGS-UHFFFAOYSA-N | CCCCCCCCCCC |
| n_paraffin | 12 | 1 | SNRUBQQJIBEYMU-UHFFFAOYSA-N | CCCCCCCCCCCC |
| n_paraffin | 13 | 1 | IIYFAKIEWZDVMP-UHFFFAOYSA-N | CCCCCCCCCCCCC |
| n_paraffin | 14 | 1 | BGHCVCJVXZWKCC-UHFFFAOYSA-N | CCCCCCCCCCCCCC |
| n_paraffin | 15 | 1 | YCOZIPAWZNQLMR-UHFFFAOYSA-N | CCCCCCCCCCCCCCC |
| n_paraffin | 16 | 1 | DCAYPVUWAIABOU-UHFFFAOYSA-N | CCCCCCCCCCCCCCCC |
| n_paraffin | 17 | 1 | NDJKXXJCMXVBJW-UHFFFAOYSA-N | CCCCCCCCCCCCCCCCC |
| n_paraffin | 18 | 1 | RZJRJXONCZWCBN-UHFFFAOYSA-N | CCCCCCCCCCCCCCCCCC |
| iso_paraffin | 4 | 1 | NNPPMTNAJDCUHE-UHFFFAOYSA-N | CC(C)C |
| iso_paraffin | 5 | 1 | QWTDNUCVQCZILF-UHFFFAOYSA-N | CCC(C)C |
| iso_paraffin | 6 | 3 (>1) | AFABGHUZZDYHJO-UHFFFAOYSA-N | CCCC(C)C |
| iso_paraffin | 7 | 2 (>1) | BZHMBWZPUJHVEE-UHFFFAOYSA-N | CC(C)CC(C)C |
| iso_paraffin | 8 | 4 (>1) | HDGQICNBXPAKLR-UHFFFAOYSA-N | CCC(C)CC(C)C |
| iso_paraffin | 9 | 3 (>1) | DOGIHOCMZJUJNR-UHFFFAOYSA-N | CCCCC(C)CCC |
| iso_paraffin | 10 | 3 (>1) | IALRSQMWHFKJJA-UHFFFAOYSA-N | CCCCCC(C)CCC |
| iso_paraffin | 11 | 1 | CNPVJWYWYZMPDS-UHFFFAOYSA-N | CCCCCCCCC(C)C |
| iso_paraffin | 12 | 1 | HTZWVZNRDDOFEI-UHFFFAOYSA-N | CCCCCCCCC(C)CC |
| iso_paraffin | 16 | 1 | BANXPJUEBPWEOT-UHFFFAOYSA-N | CCCCCCCCCCCCCC(C)C |
| iso_paraffin | 17 | 1 | HZRSSYSIMDTFLT-UHFFFAOYSA-N | CCCCCCCCCC(C)CCCCCC |
| cyclo_mono | 5 | 1 | RGSFGYAAUTVSQA-UHFFFAOYSA-N | C1CCCC1 |
| cyclo_mono | 6 | 2 (>1) | GDOPTJXRTPNYNR-UHFFFAOYSA-N | CC1CCCC1 |
| cyclo_mono | 7 | 1 | UAEPNZWRGJTJPN-UHFFFAOYSA-N | CC1CCCCC1 |
| cyclo_mono | 8 | 1 | IIEWJVIFRVWJOD-UHFFFAOYSA-N | CCC1CCCCC1 |
| cyclo_mono | 9 | 2 (>1) | DEDZSLCZHWTGOR-UHFFFAOYSA-N | CCCC1CCCCC1 |
| cyclo_mono | 10 | 1 | GGBJHURWWWLEQH-UHFFFAOYSA-N | CCCCC1CCCCC1 |
| cyclo_mono | 11 | 1 | HLTMUYBTNSVOFY-UHFFFAOYSA-N | CCCCCC1CCCCC1 |
| cyclo_mono | 12 | 1 | QHWAQXOSHHKCFK-UHFFFAOYSA-N | CCCCCCC1CCCCC1 |
| cyclo_mono | 13 | 1 | MSTLSCNJAHAQNU-UHFFFAOYSA-N | CCCCCCCC1CCCCC1 |
| cyclo_mono | 14 | 1 | FBXWCEKQCVOOLT-UHFFFAOYSA-N | CCCCCCCCC1CCCCC1 |
| cyclo_mono | 15 | 1 | CLMFECCMAVQYQA-UHFFFAOYSA-N | CCCCCCCCCC1CCCCC1 |
| cyclo_mono | 16 | 1 | STWFZICHPLEOIC-UHFFFAOYSA-N | CCCCCCCCCCC1CCCCC1 |
| cyclo_mono | 18 | 1 | BLRBGKYYWDBAQQ-UHFFFAOYSA-N | CCCCCCCCCCCCC1CCCCC1 |
| cyclo_di | 10 | 2 (>1) | NNBZCPXTIHJBJL-MGCOHNPYSA-N | C1CC[C@H]2CCCC[C@@H]2C1 |
| alkylbenzene | 6 | 1 | UHOVQNZJYSORNB-UHFFFAOYSA-N | c1ccccc1 |
| alkylbenzene | 7 | 1 | YXFVVABEGXRONW-UHFFFAOYSA-N | Cc1ccccc1 |
| alkylbenzene | 8 | 4 (>1) | CTQNGGLPUBDAKN-UHFFFAOYSA-N | Cc1ccccc1C |
| alkylbenzene | 9 | 4 (>1) | AUHZEENZYGFFBQ-UHFFFAOYSA-N | Cc1cc(C)cc(C)c1 |
| alkylbenzene | 10 | 2 (>1) | HFPZCAJZSCWRBC-UHFFFAOYSA-N | Cc1ccc(C(C)C)cc1 |
| alkylbenzene | 12 | 1 | LTEQMZWBSYACLV-UHFFFAOYSA-N | CCCCCCc1ccccc1 |
| alkylbenzene | 14 | 1 | CDKDZKXSXLNROY-UHFFFAOYSA-N | CCCCCCCCc1ccccc1 |
| alkylbenzene | 16 | 1 | UZILCZKGXMQEQR-UHFFFAOYSA-N | CCCCCCCCCCc1ccccc1 |
| alkylbenzene | 17 | 1 | XBEADGFTLHRJRB-UHFFFAOYSA-N | CCCCCCCCCCCc1ccccc1 |
| alkylbenzene | 18 | 1 | KWKXNDCHNDYVRT-UHFFFAOYSA-N | CCCCCCCCCCCCc1ccccc1 |
| indan_tetralin | 9 | 1 | PQNFLJBBNBOBRQ-UHFFFAOYSA-N | c1ccc2c(c1)CCC2 |
| indan_tetralin | 10 | 1 | CXWXQJXEFPUFDZ-UHFFFAOYSA-N | c1ccc2c(c1)CCCC2 |
| alkene_alkyne | 6 | 1 | LIKMAJRDDDTEIG-UHFFFAOYSA-N | C=CCCCC |
| alkene_alkyne | 8 | 1 | FXNDIJDIPNCZQJ-UHFFFAOYSA-N | C=C(C)CC(C)(C)C |

## Two-ring classes: carbon numbers present in the registry

- cyclo_di: C10
- indan_tetralin: C9, C10
- aromatic_di: none

## Ring share of the saturates: measured (G.1) vs r2 fixed split

Measured ring share = (cycloparaffins + dicycloparaffins + tricycloparaffins) / sum of the four D2425 saturate types.

| | mean | min | max |
|---|---|---|---|
| G.1 measured (r3 input, 57 fuels) | 0.2575 | 0.0273 | 0.5588 |
| r2 fixed split (cycloalkane) | 0.5260 | 0.5260 | 0.5260 |

## Mean mass share per kind: r2 vs r3

Both computed from the composition CSVs, kinds derived from SMILES with `subkind` (r2 fuels found: 57; r3 fuels found: 57).

| kind | r2 mean | r3 mean |
|---|---|---|
| n_paraffin | 0.1085 | 0.2281 |
| iso_paraffin | 0.2680 | 0.3421 |
| cyclo_mono | 0.4178 | 0.1036 |
| cyclo_di | 0.0000 | 0.0851 |
| alkylbenzene | 0.0825 | 0.1447 |
| indan_tetralin | 0.0941 | 0.0524 |
| aromatic_di | 0.0000 | 0.0000 |
| alkene_alkyne | 0.0090 | 0.0090 |
| other (not a sub-kind) | 0.0000 | 0.0000 |
| fallback lump | 0.0202 | 0.0349 |

G.1 total aromatics minus the aromatics r2 used (D1319, else D6379), over the fuels: mean +2.49 vol%, min +0.10, max +4.39.

## Assumptions and known gaps

- vol% to mass uses the class densities 0.87 (aromatic) / 0.78 (olefin) / 0.79 (saturate) g/cm3, the same as r2.
- D2425 type shares are treated as mass shares inside each HPLC fraction (saturates, aromatics); the table prints no unit for them.
- Cells below the detection limit (printed with '<') are set to 0.
- The n / iso paraffin split is a prior (above); G.1 gives total paraffins only.
- One D86-based carbon-number distribution is applied to every class.
- cyclo_di, indan_tetralin and aromatic_di use the nearest carbon number available in the registry (max_gap=99), because the registry holds very few such molecules.
- Tricycloparaffins and the minor aromatic types (indenes, acenaphthenes, acenaphthylenes, tricyclic aromatics) have no realistic representative class and go to named fallback buckets.
- The G.1 (D6379) aromatics differ from the page-92 D1319/D6379 aromatics that r2 used (see the difference above), so the aromatic total also changes between r2 and r3.
- r2's report line "Bicyclic cycloalkane candidates available: 0" is wrong: it counted rings on a temporary RDKit Mol, whose RingInfo reports 0 rings. r3 keeps the Mol in a variable.
