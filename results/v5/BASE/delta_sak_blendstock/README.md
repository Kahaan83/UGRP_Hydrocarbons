# delta_sak_blendstock: an open, named SAK composition and its opt-in blendstock

**Result: a defensible composition exists** (one sample, one lab) and is added as an *opt-in* blendstock. Nothing else changed: default pools, `BLENDSTOCKS`,
`verification:` flags and margins are untouched. Everything here is model input, not a measurement of any blend; no lab data exist.

## 1. Did AJFTD already hold a labelled SAK / CTL fuel? No.
* AJFTD GCxGC fuels 47 (label 10326) and 139 (label 12918) are POSF numbers only. The repo carries no name, producer or pathway for them
  (`ajftd_gcxgc_fuels.csv`, `ajftd_fuels_wide.csv`; the raw `ajftd_fuel_details.csv` named in the build log is not in any checkout; a web search for POSF 10326
  returned nothing). The only labelled synthetic fuels in AJFTD are `NN% CtL / ATJ-SPK / HVO + fossil` *blends* (bulk properties only, no GCxGC).
* Their bins resemble a C9/C10 alkylbenzene stream but are not Virent SAK (alkylbenzene C8/C9/C10/C11 mass %):

| fuel | C8 | C9 | C10 | C11 | cycloaromatic | iso+cyclo+n paraffin |
|---|---|---|---|---|---|---|
| Virent SAK (Feldhausen 2022) | 3.7 | 62.0 | 27.5 | 0 | 2.0 | about 4.8 |
| AJFTD 47 (POSF 10326) | 0.67 | 65.3 | 22.8 | 3.6 | 4.5 | about 1.5 |
| AJFTD 139 (POSF 12918) | 0.67 | 65.3 | 22.6 | 3.5 | 4.7 | about 1.6 |

  They stay "unlabelled aromatic stream"; they are not used as SAK.

## 2. Source (open, CC BY)
Feldhausen, Bell, Yang, Faulhaber, Boehm, Heyne, *Synthetic aromatic kerosene property prediction improvements with isomer specific characterization via GCxGC and
vacuum ultraviolet spectroscopy*, Fuel 326 (2022) 125002, doi 10.1016/j.fuel.2022.125002. Open copy: US DOT ROSAP. Sample "provided courtesy of Virent, Inc".
URL, date, licence, sha256 in `data/raw/sak/SOURCES.md`. The Virent data-sheet PDF returned an empty body; nothing kept.

## 3. Composition table (transcribed from the paper's Table 1; mass %; 26 isomers, 93.581 % of the sample)
Machine-readable: `data/raw/sak/dot_64804/feldhausen2022_table1.csv` (with SMILES).

| group | species | % m | | group | species | % m |
|---|---|---|---|---|---|---|
| C9 cycloaromatic | indane | 0.799 | | C9 AB | 1,2,3-trimethylbenzene | 0.163 |
| C8 AB | p-xylene | 0.384 | | C10 AB | sec-butylbenzene | 0.156 |
| C8 AB | m-xylene | 1.274 | | C10 AB | 1,3-diethylbenzene | 4.445 |
| C8 AB | ethylbenzene | 0.181 | | C10 AB | 1-methyl-3-n-propylbenzene | 3.213 |
| C8 AB | o-xylene | 1.844 | | C10 AB | 1-methyl-4-n-propylbenzene | 1.535 |
| C9 AB | 1-methyl-3-ethylbenzene | 20.604 | | C10 AB | 1,3-dimethyl-5-ethylbenzene | 1.122 |
| C9 AB | 1-methyl-4-ethylbenzene | 10.879 | | C10 AB | 1,4-diethylbenzene | 0.811 |
| C9 AB | 1,3,5-trimethylbenzene | 0.148 | | C10 AB | 1-methyl-3-isopropylbenzene | 0.739 |
| C9 AB | 1,2,4-trimethylbenzene | 26.686 | | C10 AB | 1-methyl-4-isopropylbenzene | 0.332 |
| C9 AB | isopropylbenzene | 0.275 | | C10 AB | 1-methyl-2-n-propylbenzene | 0.145 |
| C9 AB | n-propylbenzene | 1.920 | | C10 AB | 1,4-dimethyl-2-ethylbenzene | 0.678 |
| C9 AB | 1-methyl-2-ethylbenzene | 1.356 | | C10 AB | 1,3-dimethyl-4-ethylbenzene | 2.443 |
| | | | | C10 AB | 1,2-dimethyl-4-ethylbenzene | 10.089 |
| | | | | C10 AB | 1,2,4,5-tetramethylbenzene | 1.360 |

Arithmetic checks (tests): C8 sum 3.683 (text 3.7), C9 AB 62.031 (text 62.0), C10 27.068 (text group 27.5), coeluting peak f* 31.631 (text 31.6), total 93.581 (text 93.6).
Group-type values from the main text (sec 3.1, GCxGC Striebich template): alkylbenzene C8 3.7 / C9 62.0 / C10 27.5; cycloaromatic C9 0.8 / C10 1.2; the remaining
about 4.8 % is "various cycloalkane and iso-alkane groups, each < 1.5 %, C8 to C10" and its split lives in the paper's Supplementary Table 4 (not obtained).

## 4. Mapping onto the existing representatives
`data/real_fuels/v5/sak_feldhausen2022_bins.csv` holds the group-type bins (the same class x carbon granularity as HEFA / ATJ / Jet A), fed through the same
`map_bins` + `SubReps` route (`blendstock_v5.load_blendstock("sak_virent_feldhausen2022")`). Result (mass fractions after renormalising the 4.8 % unmapped, which
has no class split and is dropped; this raises aromatics from 95.2 % to 100 %):

| representative | w | stands for |
|---|---|---|
| mesitylene | 0.651 | all C9 alkylbenzenes (really 1,2,4-TMB 27 % + ethyltoluenes 31 % + others; true 1,3,5-TMB is 0.15 %) |
| p-cymene | 0.289 | all C10 alkylbenzenes (diethylbenzenes, dimethylethylbenzenes ...) |
| o-xylene | 0.039 | C8 alkylbenzenes |
| tetralin / indane | 0.013 / 0.008 | C10 / C9 cycloaromatics |

`represented=True` (so recipes carry the existing `represented_blendstock_margin_not_validated` reason), `synthetic=True` (counts toward the 50 % synthetic cap; whether
D7566 caps SAK separately is not checked here: **needs user decision**). No `verification:` flag or margin was edited.

## 5. Confidence and limits
* **Composition: medium-high for what it is.** Peer-reviewed, isomer-resolved GCxGC-FID/VUV, CC BY, internally consistent to the third decimal. It is **one Virent sample**
  (a 2022 development batch); commercial SAK varies by batch and by producer, and Sasol-type FT aromatics were not found in open compositional form. Not labelled as the D7566 annex product.
* **Representative error, checked only where data allow.** For the 14 isomers present in the registry (73 % of mass) vs the 4-rep mapping: Refutas nu at 15 C 0.913 vs 0.912 cSt,
  density 872 vs 870 kg/m3; at 40 C 0.757 vs 0.681 cSt (10 % low), at -20 C 1.20 vs 1.44 cSt (reps 20 % high: over-prediction of cold viscosity is the conservative direction, but it is a mesitylene extrapolation artefact, not a safe margin), density -20 C +1.7 %. The partial isomer set is itself an approximation (the 12 missing isomers hold 27 % of mass).
  Mesitylene as the only C9 aromatic is the known weak spot (audit item 8).
* No measured SAK property in the repo was available to validate predictions (the paper reports density and viscosity only as figures; not digitised).
* 4.8 % of mass is dropped, not modelled. The renormalisation overstates aromatics by about 5 %.

## 6. Code and tests
`OPTIONAL_BLENDSTOCKS` in `src/fuelblend/design/blendstock_v5.py` (not iterated by `BLENDSTOCKS`, so the WS7 report script and default pools are untouched). `tests/test_v5_sak_blendstock.py`
(5 tests: published bins, Table 1 arithmetic, mapping, default-pool identity incl. `A` matrices bit-equal before/after loading SAK, existing members' columns unchanged when SAK joins a pool);
with `tests/test_v5_ws7_design_realism.py`: 18 passed. Note: in a fresh worktree the LightGBM CRLF model files abort pytest (known from WS7); I converted them locally and restored them, uncommitted.
`delta_ws7_design_realism.py cmd_blendstocks` was not re-run, so `blendstocks.csv` still lists SAK as blocked. Claims register not edited (follow-up: a CLM entry once the user accepts the source).
