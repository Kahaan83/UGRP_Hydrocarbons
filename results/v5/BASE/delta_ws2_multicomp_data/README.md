# delta_ws2_multicomp_data: can open data reach R-W2-2 (>= 8 independent multicomponent systems per class mix)?

Sourcing and ingest only (2026-11-14). The WS2 gate was NOT re-run, no margin was computed or changed, no `verification:` flag touched.
Reproduce: `PYTHONPATH=src python scripts/v5/base/ws2_multicomp_data_ingest.py` (reuses the loader/dedup/constant tiers of
`scripts/v5/base/delta_ws2_freeze_multicomp.py`).

## Answer
**No. The 8-per-class-mix count is not reachable from public open data I could find.** Even with every missing Tm/dHfus constant curated:

| class mix | systems in set | of which C18-C36 wax-appearance (Dauphin) | scorable now (verified or ThermoML-measured) | scorable, verified only | ceiling if all constants curated | ceiling without the 4 wax systems |
|---|---|---|---|---|---|---|
| alkane_only | 8 | 4 | 2 | 1 | 8 | 4 |
| alkane_naphthene | 1 | 0 | 1 | 1 | 1 | 1 |
| with_aromatic | 2 | 0 | 0 | 0 | 2 | 2 |

alkane_only reaches 8 only by counting the four Dauphin C18-C36 wax mixtures (one wax *appearance* point each, far above the jet range) and only
after curating constants for odd C19-C35 waxes. The other two classes have a ceiling of 1 and 2: the problem is a lack of systems, not of constants.
(Files: `counts_by_class_mix.csv`, `systems.csv` with per-system `missing_any` / `missing_verified`.)

## Ingest (canonical schema)
`measurements_multicomp_liquidus.csv`: 46 rows, 11 systems, `property = freezing_point_C` (liquidus / wax-disappearance / wax-appearance /
freezing point, label in `method`), `value` degC and `T_K`, `kind = multi`, `x_mole` and `w_mass` aligned to sorted `comp_inchikeys`,
`system_id` = md5[:12] of the sorted InChIKey set, `ref = source:doi:sample`. Dedup: the Oakley 2017 ternary equals the ThermoML point
(10.1016/j.jct.2016.12.030) and was dropped; exact (system, T to 0.5 K, x to 3 dp) duplicates are guarded. Nothing here is new data: it is the
union of the sources already in the repo (ThermoML SLE, ASCENT 065A, Shariatrad 2022, Dauphin/Coutinho SoftX) in one schema.

## Leads checked (this session)
| lead | outcome |
|---|---|
| ThermoML archive v2020-09-30 (liquidus, n>=3) | already exhausted by `thermoml_sle_hc.csv`: 2 ternaries (+ non-hydrocarbon ones such as hexacosane + dibenzofuran + biphenyl, JCED 2014, excluded by the hydrocarbon rule). Papers after 2019 are not in the archive; the live ThermoML web index needs JavaScript, so no crawl was possible. |
| Aftab, Javanmardi, Nasrifar 2015 (WDT of C11+C16+C18 and C14+C16+C18, 0.9 bar) | would add 2 alkane_only systems; Elsevier paywalled, not in ThermoML. **Blocked.** (Same group as the Shariatrad paper; an author data request is the only route.) |
| Bell/Boehm/Heyne 2025 (Energy Fuels) | DOT public copy read: mixture freezing points beyond ASCENT 065A are in figures only (not digitised). No new systems. `data/raw/ws2_bell2025_efuels_si/` |
| Coutinho SoftX database | already ingested (4 multicomponent Dauphin systems + binaries). |
| Coutinho/Dirand/Provost OGST 1998 (C26 + C28 + heptane ternary) | OGST site returned HTTP 403 to scripted download; C28 and heptane constants also unverified. **Blocked.** |
| Zhou et al., Processes 2026 (DSC, C38-C50 in C7/C9/C12) | binary only. `data/raw/ws2_zhou2026_processes/` |
| Johnsen, Fluid Phase Equilib. 2022 review (arXiv 1807.08514) | 43 binary data sets; cites but does not tabulate the Floter ternary. `data/raw/ws2_johnsen2022_review/` |
| Sequeira et al. 2024 (Int. J. Thermophys.) | binary DSC solidus/liquidus; no ternary. Not downloaded. |
| Vozka 2018 Jet A / HEFA blends, AJFTD GCxGC fuels | real fuels with lumped composition: not defined systems, no component constants. Excluded by the WS2 rule. |
| Affens 1981 NRL model jet fuels (ternary in iso-octane) | no digital full text. |

## Reading
* The only routes that could materially raise counts are (1) getting the two Aftab ternaries and Dirand's C26/C28/heptane ternary (alkane_only +3 at most),
  (2) new lab measurements on jet-range alkane + cycloalkane + aromatic blends, which no open source provides. Option (2) is outside this brief (no lab testing).
* Constants worth curating anyway (cheap, WebBook calorimetry like `tm_dhfus_verified.csv`): n-hexane, n-octane, n-nonane, p-xylene, cis/trans decalin
  (the unspecified isomer mixture cannot have a single Tm), which would make 2-3 more ASCENT/ThermoML systems scorable; the odd C19-C35 waxes only help the Dauphin set.
* Needs user decision: whether to keep the 16.98 C binary margin as a transferred, unvalidated value (current state), whether to request the Aftab data from the authors.
