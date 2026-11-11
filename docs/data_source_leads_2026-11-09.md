# Data source leads — web search, 2026-11-09

Scope: new data that fills the gaps in AGENTS.md §8.2/§12. The main gaps are real ≥3-component blends, LFL, Antoine for heavy odd n-alkanes, pour/freezing, and real-fuel compositions.
Items already on the v4 plan's Task 8 list are marked **(Task 8 overlap)**.

How far each lead was checked:
- **CONFIRMED** = the page or PDF was opened, or the local data was checked.
- **SNIPPET** = seen only in search results; still needs checking.

Access cost, from cheapest to most expensive: on disk < free download < login < manual/paywalled.

## Tier A — do first

| # | Source | Fills | Access | Status | Notes |
|---|---|---|---|---|---|
| 1 | **Růžička & Majer 1994**, *J. Phys. Chem. Ref. Data* 23, 1. Doi 10.1063/1.555942; free NIST reprint at https://srd.nist.gov/jpcrdreprint/1.555942.pdf | Antoine/Psat for **all n-alkanes C5–C20**, triple point to normal boiling point. That includes C11 at low T, C13, C15, C17, C18 and C19 (§8.5, §12 step 9) | free PDF | CONFIRMED: open PDF; Table 10 (p. 23) gives the recommended Cox-equation parameters, and Table 12 gives T at selected Psat | **Correction (same day): already ingested.** `extras/fuelblend_data/refit_antoine.py` reads `data/raw/ruzicka_majer/t01_transcription.md` and adds 13 Cox points over 313–373 K. AGENTS.md §8.5 is out of date. The Table 10 transcription was verified against the paper's Tables 12 and 15 (within 0.02 % / 0.005 K). The transcription's Table 14 was wrong and has been fixed (Pa, not kPa). **Open issue:** the stored C12/C16/C14 Antoine sets deviate from Cox by 17 %/17.5 %/6 % at 40–100 °C. The cause is one ThermoML reference, doi 10.1016/j.tca.2015.03.020, which reads 17–29 % low; DIPPR, the Lemmon–Huber EOS and VDI all agree with Cox. The sets for C13/15/17/18/19 come purely from Cox but are labelled `refit_measured`. |
| 1b | JPCRD 29, 1 (2000), "Critical properties and vapor pressure equation for alkanes CnH2n+2: normal alkanes with n≤36…" Doi 10.1063/1.556054; reprint at https://srd.nist.gov/jpcrdreprint/1.556054.pdf | Psat equation for n-alkanes up to C36 (C21+ wax-range components) | free PDF (likely) | SNIPPET (listed on the cover of #1) | Use it for C21–C24, which are now `check` status. |
| 2 | **AJFTD full database** (https://altjetfuels.illinois.edu/about_fuels) | Real-fuel properties plus **>160 fuel records with GCxGC composition** | login (you already have one) | CONFIRMED: the About page says "over 25,000 records" and "over 160 fuel records" with GCxGC; the properties listed include density, viscosity, flash point and heat of combustion | We have 81 fuels scraped. "Records" may mean test rows rather than fuels. **Your action:** with filters cleared, read the "Showing x of N" count (§12 step 5), and read `/about/database-access` for the terms before any re-scrape. The 160 GCxGC records would supply the missing real-fuel compositions. |
| 3 | `chemicals` Python library (`chemicals.safety`: `T_flash`, `LFL`). Docs at https://chemicals.readthedocs.io/_modules/chemicals/safety.html | LFL (currently 25 measured out of 695) and pure flash | pip, MIT, already in the v4 tech stack | SNIPPET | **Before you count these as measured, check where each value comes from.** If `LFL()` falls back to an estimation method, a value is not data. Keep only the tabulated lookups and record their source table. |

## Ingested 2026-11-09 (from PDFs the user downloaded)

Raw transcriptions and READMEs are in `extras/fuelblend_data/data/raw/{flash_mixtures,pires2018,vozka2019}/`. The loaders are `build_unified.load_ha_lee_2016()` and `load_paper_fuel_tables()`, tested in `tests/test_paper_ingest_2026_11.py`.
- **#7 Ha & Lee 2016:** 12 flash rows, 9 of them ternary. Multi-component HC flash records go from 3 to 12. Their pure n-dodecane (83 °C) reads about 9 °C high and is flagged as a conflict.
- **#5 Pires 2018:** 130 real-fuel property rows across 12 fuels, with the WSU and AFRL labs kept apart, plus class compositions (low quality, 48–180 % GC closure). There is no tabulated viscosity.
- **#4 Vozka 2019:** SI only. 27 rows: freezing points of 25 Jet A blends with FT/HEFA/SIP/ATJ, and the net heat of combustion of the neat blendstocks; class compositions for 11 blends. The main article (density, viscosity and flash of the blends) is still needed.

## Tier B — papers or supplements, mostly manual

| # | Source | Fills | Access | Status | Notes |
|---|---|---|---|---|---|
| 4 | Vozka, Vrtiška, Šimáček, Kilaz, *Energy & Fuels* 33(4) 3275 (2019), "Impact of alternative fuel blending components on fuel composition and properties in blends with Jet A". ACS figshare: https://acs.figshare.com/articles/journal_contribution/Impact_of_Alternative_Fuel_Blending_Components_on_Fuel_Composition_and_Properties_in_Blends_with_Jet_A/7844741 | Real blends: GCxGC composition, density, viscosity, flash, freezing, D86, heat of combustion | figshare may host the SI for free | SNIPPET **(Task 8 overlap)** | The figshare record is new: try it before using the paywalled article. The HEFA-feedstock companion paper (E&F 2018) has the same property set. |
| 5 | Pires, Han, Kramlich, Garcia-Pérez, *BioResources* 13(2) 2632 (2018), "Chemical composition and fuel properties of alternative jet fuels". https://rosap.ntl.bts.gov/view/dot/50415, data in `dot_50415_DS1.pdf` | 12 fuels (9 alternative and 3 commercial jet fuels): n/iso-paraffin, olefin, naphthene and aromatic fractions plus viscosity, density, flash point and surface tension | free (US DOT ROSAP) | CONFIRMED (abstract); DS1 not yet opened | The composition is by class only, so it fits the class × carbon-number pseudo-component route. |
| 6 | JCED 2017 doi 10.1021/acs.jced.7b00466, a ternary of n-butylcyclohexane + n-C12/n-C16 + toluene/butylbenzene/hexylbenzene; and JCED 2016 doi 10.1021/acs.jced.6b00516, the binary version | Real ternary hydrocarbon **flash points** | paywalled article | CONFIRMED locally: the density and viscosity are **already in our ThermoML tarball**. The flash points are **not**. Of 11,923 archive XML files, none has a flash property; it appears in the abstract only. | The flash values would have to come from the paper's tables, through the university library. The 12 ternary density/viscosity property sets are already usable. Check that the extractor kept them. |
| 7 | *Journal of the Korean Gas Society* (2016), n-nonane + n-decane + n-dodecane, Seta closed cup. https://koreascience.or.kr/article/JAKO201610364969655.page | Measured ternary n-alkane flash points | probably a free PDF | SNIPPET | Small, but these are exactly the ternary hydrocarbon flash data that are missing (only 3 rows now). |
| 8 | Liaw et al. ternary flash papers (*Combust. Flame* 138 (2004) 308, etc.) | Ternary flash | paywalled | SNIPPET **(Task 8 overlap)** | These are mostly model papers. Use only the measured tables they contain. |
| 9 | SLE / wax: Coutinho et al., *Fluid Phase Equilib.* 1995 (about 60 binary alkane SLE systems, >1000 points); Coutinho & Stenby, *IECR* 1996 (multicomponent n-alkanes); ACS Omega 2022, quaternary n-C11/C14/C16/C18 wax **disappearance** temperatures (PMC9134231, open) | Freezing/cloud (v4 Tasks 12–13) | 1995/1996 paywalled; 2022 open | SNIPPET | WDT is not the same as WAT/cloud point, so keep the measurement type in the `ref`/method column. The ThermoML SLE data already on disk comes first (plan item 6). |
| 10 | Edwards, *Jet Fuel Properties*, AFRL-RQ-WP-TR-2020-0017 (DTIC AD1093317), https://apps.dtic.mil/sti/pdfs/AD1093317.pdf | Reference fuels A-1/A-2/A-3 and alternative fuels: viscosity, flash, freezing | browser download (DTIC blocks scripts) | SNIPPET **(Task 8 overlap)** | Stanford HyChem's A2 page (https://web.stanford.edu/group/haiwanglab/HyChem/fuels/A2_spec.html) has a "key properties" table for POSF10325. |

## Checked and closed

- **No ThermoML release newer than 2020-09-30 exists.** The data.nist.gov record (mds2-2422) still points to `ThermoML.v2020-09-30.tgz` (data through 2019). This closes that Task 8 item.
- **ThermoML has no flash-point property at all:** 0 of 11,923 XML files. Flash data has to come from paper tables, Leenhouts, or Li.

## Not recommended

- **Shell.ai Hackathon 2025** (5-component blends, 10 targets, 2000 training rows): the components and properties are anonymized, the license is the competition's, and the provenance is unknown. It conflicts with the v4 no-fabrication/traceability rules.
- **Zenodo 10.5281/zenodo.10670349** (biodiesel cold flow vs FAME profile): the compounds are FAMEs, so they fail `is_hc`.
- **Heyne et al., *Fuel* 312 (2022)** (rosap dot/64805): one surrogate in a methods paper, not a dataset.
- **JETSCREEN/ASCENT JSON fuel database** (Blakey et al., *Front. Energy Res.* 2022; about 30k conventional and 400 alternative fuels): no public download found; parts are commercially restricted.
- **NIOSH Pocket Guide LEL as CSV**: no bulk export found, and n-dodecane has no entry.

## Sources
- https://srd.nist.gov/jpcrdreprint/1.555942.pdf
- https://altjetfuels.illinois.edu/about_fuels
- https://chemicals.readthedocs.io/_modules/chemicals/safety.html
- https://acs.figshare.com/articles/journal_contribution/Impact_of_Alternative_Fuel_Blending_Components_on_Fuel_Composition_and_Properties_in_Blends_with_Jet_A/7844741
- https://rosap.ntl.bts.gov/view/dot/50415
- https://rosap.ntl.bts.gov/view/dot/64805
- https://trc.nist.gov/ThermoML/10.1021/acs.jced.7b00466.html
- https://koreascience.or.kr/article/JAKO201610364969655.page
- https://rcastoragev2.blob.core.windows.net/fac5385c3077144472cce8760da0bd29/ao1c07072.PMC9134231.pdf
- https://apps.dtic.mil/sti/pdfs/AD1093317.pdf
- https://web.stanford.edu/group/haiwanglab/HyChem/fuels/A2_spec.html
- https://data.nist.gov/od/id/ark:/88434/mds2-2422
- https://alchemist-shellai-hackathon-2025.readthedocs.io/en/latest/
- https://www.frontiersin.org/articles/10.3389/fenrg.2022.771325/full
