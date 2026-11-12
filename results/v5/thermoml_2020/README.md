# ThermoML: nine J. Chem. Eng. Data 2020 files not in the bulk archive (v5)

Downloaded 2026-11-12 with the user's go-ahead from `https://trc.nist.gov/ThermoML/10.1021/acs.jced.<id>.xml`.
These files are served by NIST TRC but are absent from `ThermoML.v2020-09-30.tgz` and from NIST's browse index
(how they were found: `docs/data_audit_2026-11-11.md` section 5b). Raw files: `data/raw/thermoml_2020/` (gitignored).

Run (extractor unchanged):
`PYTHONPATH=src python -m fuelblend.data.thermoml_extract data/raw/thermoml_2020/thermoml_jced_2020_unlisted.tgz --out results/v5/thermoml_2020/thermoml_hc_2020.csv`

Output `thermoml_hc_2020.csv`: same schema as `results/v5/thermoml_reextract/thermoml_hc_v5.csv`, native ThermoML units
(Pa*s, m2/s, kg/m3, kPa). 9 files read, 0 errors, 3,768 hydrocarbon rows kept; `log.txt` has the not-kept counts.

| file (acs.jced.) | content kept | T, K |
|---|---|---|
| 0c00020 | iso-dodecane + n-propylbenzene / 1,3,5- / 1,2,4-trimethylbenzene: 198 binary viscosity, 264 binary density | 253.15-333.15 |
| 0c00280 | n-dodecane + four C10 aromatics: 220 binary viscosity, 308 binary density | 288.15-333.15 |
| 0c00581 | iso-dodecane + four C10 aromatics: 185 binary kinematic viscosity, 259 binary density | 288.15-333.15 |
| 0c00309 | isopropylcyclohexane + n-tridecane: 99 binary viscosity, 99 binary density | 293.15-333.15 |
| 9b01171 | exo-tetrahydrodicyclopentadiene + n-decane: 99 binary viscosity, 99 binary density | 293.15-333.15 |
| 9b01167 | alkyladamantanes + n-nonane / n-undecane (6 systems): 726 binary viscosity, 726 binary density | 293.15-343.15 |
| 9b01026 | benzene + n-heptane + n-octane: 108 ternary density, 16 binary density | 298.15-318.15 |
| 9b01208 | pure vapour pressure, 3 compounds: 40 points | 400-620 |
| 0c00111 | pure vapour pressure, 3 n-alkanes: 3 points (binary boiling temperatures are not a property the extractor keeps) | 377-411 |

Mixture viscosity: 1,527 rows in 19 systems; 17 of the 19 component sets do not occur among the 171 mixture-viscosity
systems of the 2019 archive. 33 mixture viscosity rows are at or below 263.15 K (all from 0c00020).
**Correction 2026-11-12:** that comparison is against the ThermoML archive only. 15 of the 19 systems are already in
`data/unified/measurements_dedup.csv` through NIST TDE and Leenhouts. Of the 1,278 mixture viscosity rows left after end points are
removed and repeats merged, 756 are at a (system, temperature) already in the unified data, 342 are new temperatures of known
systems and 180 belong to 4 new systems. Of 1,599 mixture density rows, 1,440 are in systems the unified data lacks.
Pure-compound rows in each file (the authors' own pure measurements) are kept as n_comp = 1.

## Status: extracted only

2026-11-12: converted to the lab-table schema (`data/real_fuels/v5/thermoml2020_lab.csv`) and scored once as a held-out test, see `results/v5/S1b/`.

NOT merged into `data/unified/`, `data/processed/v2/thermoml_clean_hc.csv`, any split, or any evaluation. No model was fitted or
scored on these rows. (Superseded 2026-11-12: the earlier sentence here asked for a held-out-versus-training decision on "17 unseen systems"; only 4 systems are new to the unified data, see the correction above and `results/v5/S1b/`.) The files are not in NIST's index, so they may have had less checking than the
archive; the NIST liability statement applies. ThermoML carries no flash point, so the flash points reported in the three
Luning Prak papers are not here.

## sha256 of the downloaded files
```
dd2e3e6b91d27516f2066a70cea49b6ada01aa7312d36cf5264213a79103dd36 acs.jced.acs.jced.0c00020.xml
9a86958b0f4cc96b347cf69ed3bdf5cf4985f93ee2be20b3f7c22533809083eb acs.jced.acs.jced.0c00111.xml
d53098793b7f6659a83fb07cf890898f390eb035318894f75cdbdcef7943c5df acs.jced.acs.jced.0c00280.xml
48432261f4300c487602457837f7b42c6fdbf3e07ce23e92a16f11f60f4cfa64 acs.jced.acs.jced.0c00309.xml
4f905e1db628e570e13cee2765855e6cbb7d2598ccf0db9235631c2a95760e9f acs.jced.acs.jced.0c00581.xml
5d207609ce2a73fbbac7cf774ce1ba09053211c7cfe0b0744aa50ac594b84cb8 acs.jced.acs.jced.9b01026.xml
c8829a5ed8f212559fd8d7b299060c37b312d95d1162b4a87fbbfa5d59e230a7 acs.jced.acs.jced.9b01167.xml
7c8f99d1c8791a7b84e16acafb85e567160ea093214527a67fa3e858e4c535ef acs.jced.acs.jced.9b01171.xml
ac6a4307798bf637fd81fd0efec401645252daec6628b2e87595411ed7594829 acs.jced.acs.jced.9b01208.xml
```
