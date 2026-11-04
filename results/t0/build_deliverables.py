import csv
from pathlib import Path
import pandas as pd

def build_data_sources():
    out_path = Path("docs/data_sources.csv")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    
    fieldnames = [
        "source_id",
        "name",
        "what_it_provides",
        "access_type",
        "url_or_doi",
        "license_or_terms",
        "redistribution_allowed",
        "files_in_repo",
        "used_by",
        "properties",
        "notes"
    ]
    
    rows = [
        {
            "source_id": "leenhouts_benchmark",
            "name": "Leenhouts et al. 2026 Chemical Mixture Benchmark",
            "what_it_provides": "Curated experimental pure and mixture flash point, dynamic viscosity, DCN, MON, and density",
            "access_type": "open download",
            "url_or_doi": "doi:10.5281/zenodo.19914760",
            "license_or_terms": "CC BY 4.0 / MIT",
            "redistribution_allowed": "yes",
            "files_in_repo": "data/raw/zenodo/leenhouts_benchmark/, data/interim/leenhouts_benchmark/, data/raw/git/mixture_evaluation/",
            "used_by": "build_unified.py, ingest_dcn.py, scripts/t7_ingest_leenhouts_flash.py, build_octane_hc.py",
            "properties": "flash_point_C, viscosity_dyn_mPas, dcn, mon, density_kg_m3",
            "notes": "Upstream code at gitlab.kuleuven.be/creas/vermeiregroup/mixture_evaluation. Contains curated literature data (Kuzhagaliyeva 2022 MON, DCN, etc.)."
        },
        {
            "source_id": "nist_tde_chemixhub",
            "name": "NIST ThermoData Engine (TDE) / CheMixHub / Bilodeau et al.",
            "what_it_provides": "Organic pure and binary mixture dynamic viscosity measurements",
            "access_type": "open download",
            "url_or_doi": "doi:10.5281/zenodo.8042966 / doi:10.1016/j.cej.2023.142454",
            "license_or_terms": "CC BY 4.0 (Zenodo 8042966) / MIT (CheMixHub repo)",
            "redistribution_allowed": "yes",
            "files_in_repo": "data/raw/git/chemixhub/, data/raw/zenodo/bilodeau_nist_viscosity/pretrained_models.zip, data/processed/viscosity_pure_hc.csv, data/processed/viscosity_mix_hc.csv",
            "used_by": "build_viscosity.py, build_unified.py",
            "properties": "viscosity_dyn_mPas",
            "notes": "CheMixHub NeurIPS 2025 benchmark clone. Raw column labeled 'Log viscosity / cP' is ln(eta / Pa*s); deduplicated against exact repeats."
        },
        {
            "source_id": "thermoml",
            "name": "NIST TRC ThermoML Archive (v2020-09-30)",
            "what_it_provides": "Experimental pure and multicomponent liquid viscosity, density, and vapor pressure",
            "access_type": "open download",
            "url_or_doi": "doi:10.18434/mds2-2422 / https://data.nist.gov/od/ds/mds2-2422/ThermoML.v2020-09-30.tgz",
            "license_or_terms": "U.S. Government Public Domain (17 U.S.C. 105)",
            "redistribution_allowed": "yes",
            "files_in_repo": "data/raw/thermoml/ThermoML.v2020-09-30.tgz, data/processed/thermoml_hc.csv",
            "used_by": "thermoml_extract.py, build_unified.py, refit_antoine.py",
            "properties": "viscosity_dyn_mPas, density_kg_m3, vapor_pressure_kPa",
            "notes": "189 MB archive of XML thermophysical records from 5 IUPAC-cooperating journals through 2019. Liquid-phase filtered; excess quantities excluded."
        },
        {
            "source_id": "ajftd",
            "name": "Alternative Jet Fuel Test Database (AJFTD)",
            "what_it_provides": "Certified real jet-fuel batch property tests (distillation curves, API gravity, flash, freeze, smoke point, etc.)",
            "access_type": "credentialed",
            "url_or_doi": "https://ajftd.aviation.illinois.edu/ / https://ascent.aero/project/alternative-jet-fuel-test-database-ajftd/",
            "license_or_terms": "Academic research use only (FAA ASCENT)",
            "redistribution_allowed": "unclear",
            "files_in_repo": "data/raw/ajftd/, ajftd_data.csv, ajftd_data.json, data/processed/ajftd_fuels.csv, data/processed/ajftd_composition.csv, data/processed/smoke_point_ajftd.csv, data/processed/heat_of_combustion_real_fuels.csv",
            "used_by": "ajftd_scraper.py, ajftd_parse_html.py, ingest_ajftd.py, extract_ajftd_properties.py, build_heat_of_combustion_real_fuels.py, build_sooting_real_fuels.py, build_unified.py",
            "properties": "distillation_T*_C, density_kg_m3, flash_point_C, freezing_point_C, smoke_point_mm, net_heat_of_combustion, aromatics_vol_pct",
            "notes": "Scraped from password-protected UIUC/FAA ASCENT database. Access terms require approved account; redistribution of bulk data is unclear so raw files must not be published."
        },
        {
            "source_id": "flashpoint_li",
            "name": "Li / Sun et al. 2020 Curated Pure-Compound Flash Point Collection",
            "what_it_provides": "Experimental pure-compound flash points (14,696 total, 2,194 pure hydrocarbon rows)",
            "access_type": "open download",
            "url_or_doi": "doi:10.6084/m9.figshare.9275210.v1 / arXiv:2002.11315",
            "license_or_terms": "CC BY 4.0",
            "redistribution_allowed": "yes",
            "files_in_repo": "data/raw/flashpoint_li/integrated_dataset_public.csv, data/raw/flashpoint_li/flashpoint_hc.csv",
            "used_by": "build_unified.py",
            "properties": "flash_point_C",
            "notes": "Public figshare dataset by Dane Morgan et al. accompanying Mol. Inf. 2020 (39:1900101). Open/closed cup unstratified in original CSV."
        },
        {
            "source_id": "saldana11",
            "name": "Saldaña et al. 2011 / 2012 Experimental Flash Points and Viscosities",
            "what_it_provides": "Experimental pure and blend flash points and pure viscosity",
            "access_type": "open download",
            "url_or_doi": "doi:10.1021/ef2009228 / doi:10.1021/ef300900e",
            "license_or_terms": "ACS Copyright (Scholarly fair use extract)",
            "redistribution_allowed": "unclear",
            "files_in_repo": "data/raw/flashpoint_li/saldana11_flashpoints.csv, data/raw/flashpoint_li/2011EF3900.zip",
            "used_by": "build_unified.py",
            "properties": "flash_point_C, viscosity_dyn_mPas",
            "notes": "Literature data extracted from journal tables/SI; redistribution of extracted data tables generally protected by fair use, but copyright status of compilation is unclear."
        },
        {
            "source_id": "crc647",
            "name": "CRC Report No. 647: World Fuel Sampling Program",
            "what_it_provides": "Comprehensive specification testing across 57 global commercial aviation fuel shipments",
            "access_type": "open download",
            "url_or_doi": "https://crcao.org/wp-content/uploads/2025/07/CRC-647.pdf",
            "license_or_terms": "CRC Copyright (Public report)",
            "redistribution_allowed": "unclear",
            "files_in_repo": "data/raw/pdfs/CRC_647_world_fuel_sampling.pdf, data/raw/pour_point/crc647_fuel_properties.csv, data/processed/crc_fuels.csv, fuels.csv",
            "used_by": "parse_crc_tables.py, ingest_crc.py, build_unified.py, src/fuelblend/mixture/fuel_hindcast.py",
            "properties": "pour_point_C, freezing_point_C, flash_point_C, viscosity_kin_cSt, density_kg_m3, distillation_T*_C, net_heat_content_btu_lb",
            "notes": "Primary real aviation fuel validation set. Publicly downloadable report, but CRC retains copyright; redistribution of extracted tables is unclear."
        },
        {
            "source_id": "crc_av2718",
            "name": "CRC Report AV-27-18: High-Aromatic Aviation Fuels Investigation",
            "what_it_provides": "Real aviation fuel physical and chemical properties for 20 specialized POSF fuels",
            "access_type": "open download",
            "url_or_doi": "https://crcao.org/wp-content/uploads/2022/03/CRC-AV-27-18_Final-Report_UDR-TR-2021-90.pdf",
            "license_or_terms": "CRC / UDRI Copyright (Public report)",
            "redistribution_allowed": "unclear",
            "files_in_repo": "data/raw/pdfs/CRC_AV-27-18_fuel_properties.pdf, data/raw/pour_point/av2718_fuel_properties.csv, data/processed/crc_fuels.csv",
            "used_by": "parse_crc_tables.py, ingest_crc.py, build_unified.py",
            "properties": "flash_point_C, viscosity_kin_cSt, density_kg_m3, distillation_T*_C, aromatics_vol_pct, vapor_pressure_d6378_kpa, vapor_pressure_d2879_kpa",
            "notes": "Publicly accessible technical report from Coordinating Research Council; redistribution rights of parsed tabular data unclear."
        },
        {
            "source_id": "crc_av1109",
            "name": "CRC Project AV-11-09: Cold-Flow Properties of Aviation Turbine Fuels",
            "what_it_provides": "High-precision freezing point measurements on cold-flow fuel samples",
            "access_type": "open download",
            "url_or_doi": "https://crcao.org/wp-content/uploads/2019/05/FINAL-REPORT-Freezing-point-CRC-project-AV-11-09-v3a.pdf",
            "license_or_terms": "CRC Copyright (Public report)",
            "redistribution_allowed": "unclear",
            "files_in_repo": "data/raw/pdfs/CRC_AV-11-09_freezing_point.pdf, data/raw/pour_point/av1109_freezing_point.csv, data/processed/crc_fuels.csv",
            "used_by": "parse_crc_tables.py, ingest_crc.py, build_unified.py",
            "properties": "freezing_point_C",
            "notes": "Freezing point phase behavior on QinetiQ fuel batches. Public report, commercial redistribution unclear."
        },
        {
            "source_id": "chew2024",
            "name": "Chew et al. 2024 Pure Liquid Viscosity vs Temperature Dataset",
            "what_it_provides": "3,582 curated dynamic viscosity measurements vs temperature across 1,005 molecules (845 HC rows)",
            "access_type": "open download",
            "url_or_doi": "doi:10.1186/s13321-024-00820-5 / doi:10.5281/zenodo.8350587",
            "license_or_terms": "CC BY-NC 4.0",
            "redistribution_allowed": "yes",
            "files_in_repo": "data/raw/chew2024/si_viscosity_dataset.csv, data/raw/chew2024/chew2024_viscosity.csv, data/processed/viscosity_chew_hc.csv",
            "used_by": "parse_chew.py, build_unified.py",
            "properties": "viscosity_dyn_mPas",
            "notes": "Curated by Schrödinger, Inc. Released under CC BY-NC 4.0; non-commercial redistribution allowed; must not be bundled into commercial distributions."
        },
        {
            "source_id": "nist_rp1_nistir6646",
            "name": "NISTIR 6646: Thermophysical Properties of Rocket Propellant RP-1 (Phase I)",
            "what_it_provides": "33-component GC-MS chemical composition, 3 surrogates, atmospheric/compressed density, viscosity, thermal conductivity",
            "access_type": "open download",
            "url_or_doi": "https://nvlpubs.nist.gov/nistpubs/Legacy/IR/nistir6646.pdf / doi:10.6028/NIST.IR.6646",
            "license_or_terms": "U.S. Government Public Domain (17 U.S.C. 105)",
            "redistribution_allowed": "yes",
            "files_in_repo": "data/raw/pdfs/NISTIR_6646_RP1.pdf, data/processed/rp1_characterization.csv, data/processed/rp1_surrogates.csv, data/processed/rp1_properties.csv, data/raw/rp1/rp1_blends_composition.csv, data/raw/rp1/rp1_blends_properties.csv",
            "used_by": "ingest_rp1.py, build_unified.py",
            "properties": "density_kg_m3, viscosity_kin_cSt, thermal_conductivity_w_m_k",
            "notes": "Official NIST technical publication; full text and tables in public domain."
        },
        {
            "source_id": "fortin2012",
            "name": "Fortin et al. 2012: Advanced Distillation Curves of 11 Orthogonal RP-1 Blends and RP-2",
            "what_it_provides": "Thermodynamic distillation curve data (ADC method) for 11 orthogonal RP-1 blends, RP-1 4572, and 2 RP-2 fuels",
            "access_type": "open download",
            "url_or_doi": "doi:10.1021/ef3004009",
            "license_or_terms": "ACS Copyright (Scholarly fair use extract)",
            "redistribution_allowed": "unclear",
            "files_in_repo": "data/raw/rp1/fortin_2012_rp1_blends.pdf, data/raw/rp1/fortin2012_tables.csv, data/processed/fortin2012_rp1_blends.csv",
            "used_by": "fortin_extract.py, build_fortin2012_processed.py, build_unified.py",
            "properties": "distillation_T*_C",
            "notes": "Energy & Fuels paper by NIST researchers; published by ACS. Distillation measurements transcribed; blend molecular recipes proprietary to AFRL."
        },
        {
            "source_id": "nist_webbook",
            "name": "NIST Chemistry WebBook (SRD 69)",
            "what_it_provides": "Antoine vapor pressure coefficients, phase transition enthalpies, boiling points, critical constants",
            "access_type": "open download",
            "url_or_doi": "https://webbook.nist.gov/chemistry/ / doi:10.18434/T4D303",
            "license_or_terms": "NIST Standard Reference Data / U.S. Government Open Data",
            "redistribution_allowed": "yes",
            "files_in_repo": "data/raw/webbook/, data/processed/antoine_best.csv, data/processed/antoine_all_sets.csv, data/raw/webbook_thermo.csv",
            "used_by": "webbook_antoine.py, webbook_antoine_extended.py, build_unified.py, refit_antoine.py",
            "properties": "vapor_pressure_kPa, boiling_point_K, critical_temperature_K, critical_pressure_kPa",
            "notes": "Sourced via polite rate-limited web scraping (<=1 req/s). Constants converted from log10(bar) to log10(kPa)."
        },
        {
            "source_id": "zabetakis_b627",
            "name": "Zabetakis 1965: Flammability Characteristics of Combustible Gases and Vapors (USBM Bulletin 627)",
            "what_it_provides": "Lower Flammability Limit (LFL) values for pure hydrocarbons, temperature coefficients, flash point correlations",
            "access_type": "open download",
            "url_or_doi": "https://digital.library.unt.edu/ark:/67531/metadc1444558/",
            "license_or_terms": "U.S. Bureau of Mines Public Domain (U.S. Government work)",
            "redistribution_allowed": "yes",
            "files_in_repo": "data/raw/pdfs/Zabetakis_Bulletin_627_flammability.pdf, data/raw/lfl/lfl_measured.csv, data/raw/lfl/lfl.csv",
            "used_by": "parse_zabetakis_a1.py, fill_lfl.py, build_unified.py",
            "properties": "lfl_vol_pct",
            "notes": "Historic U.S. government technical bulletin. Contains primary experimental LFL tables in Appendix A-1."
        },
        {
            "source_id": "pubchem",
            "name": "PubChem PUG-View / PUG-REST (NCBI)",
            "what_it_provides": "Experimental vapor pressure, boiling point, flash point, heat of combustion for 368 pure hydrocarbons",
            "access_type": "open download",
            "url_or_doi": "https://pubchem.ncbi.nlm.nih.gov/ / doi:10.1093/nar/gkaa974",
            "license_or_terms": "NCBI / NIH Open Data (U.S. Government Public Domain)",
            "redistribution_allowed": "yes",
            "files_in_repo": "data/raw/pubchem/json/, data/raw/pubchem/pubchem_experimental.csv, data/raw/pubchem/inchikey_cid_map.json",
            "used_by": "fetch_pubchem.py",
            "properties": "vapor_pressure_kPa, boiling_point_C, flash_point_C, heat_of_combustion_MJ_kg",
            "notes": "Fetched via NCBI PUG-REST and PUG-View APIs under rate limits (<=4 req/s). Public domain open biomedical data."
        },
        {
            "source_id": "ysi_database_vol2",
            "name": "Yield Sooting Index (YSI) Database Volume 2 (McEnally et al. 2017)",
            "what_it_provides": "Unified experimental sooting tendencies (YSI-molar scale) for 447 compounds (229 hydrocarbons)",
            "access_type": "open download",
            "url_or_doi": "doi:10.7910/DVN/7HGFT8",
            "license_or_terms": "CC0 1.0 Universal",
            "redistribution_allowed": "yes",
            "files_in_repo": "data/raw/sooting/Detailed_YSI_Database_Volume_2.tab, data/raw/sooting/YSI_Database_Volume_2.tab, data/raw/sooting/ysi_compounds.csv",
            "used_by": "fetch_ysi.py, build_tsi_table.py",
            "properties": "ysi",
            "notes": "Published on Harvard Dataverse by Yale combustion research group. Freely reusable without restrictions."
        },
        {
            "source_id": "boehm_tsi_2022",
            "name": "Boehm et al. 2022 / ASCENT Project 065a/066 / Mensch 2010 / Li 2012",
            "what_it_provides": "Threshold Sooting Index (TSI) conversion constants, referee fuel benchmarks, pure smoke points",
            "access_type": "open download",
            "url_or_doi": "doi:10.1021/acs.energyfuels.1c03794 / https://rosap.ntl.bts.gov/view/dot/64807/",
            "license_or_terms": "FAA ASCENT / U.S. DOT ROSAP Public Technical Report",
            "redistribution_allowed": "yes",
            "files_in_repo": "data/raw/sooting/boehm2022.pdf, data/raw/sooting/mensch2009.pdf, data/processed/tsi_table.csv, data/processed/tsi_constants.json",
            "used_by": "build_tsi_table.py",
            "properties": "tsi, smoke_point_mm",
            "notes": "Public technical research sponsored by FAA ASCENT program; formulas and experimental anchor values openly published."
        },
        {
            "source_id": "vom_lehn_octane_2021",
            "name": "vom Lehn et al. 2021 / CombustDB / Kuzhagaliyeva et al. 2022",
            "what_it_provides": "Experimental and predicted RON and MON for 446 pure hydrocarbons",
            "access_type": "open download",
            "url_or_doi": "doi:10.1016/j.jaecs.2021.100018 / doi:10.1038/s42004-022-00728-6",
            "license_or_terms": "CC BY 4.0",
            "redistribution_allowed": "yes",
            "files_in_repo": "data/processed/octane_ron_mon_hc.csv",
            "used_by": "build_octane_hc.py",
            "properties": "ron, mon",
            "notes": "Open-access datasets published under Creative Commons Attribution licenses."
        },
        {
            "source_id": "crc_handbook_melting_bradley",
            "name": "CRC Handbook of Chemistry and Physics / Bradley Open Notebook",
            "what_it_provides": "Pure-compound experimental freezing/melting points for 905 hydrocarbons",
            "access_type": "open download",
            "url_or_doi": "doi:10.6084/m9.figshare.1031637 / CRC Handbook 95th ed. (2014)",
            "license_or_terms": "CC0 1.0 (Bradley) / CRC Handbook fair use",
            "redistribution_allowed": "unclear",
            "files_in_repo": "data/processed/freezing_point_pure_hc.csv, data/raw/pour_point/pure_compound_freezing_points.csv",
            "used_by": "build_freezing_point_hc.py, build_pure_points.py, build_unified.py",
            "properties": "freezing_point_C",
            "notes": "Bradley Open Notebook data is CC0; CRC Handbook organics table values are factual data extracted under scholarly fair use, redistribution unclear."
        },
        {
            "source_id": "nasa_kirgina_cold_flow",
            "name": "NASA CR-179472 (Cotterell 1986) & Kirgina et al. 2019/2021",
            "what_it_provides": "Real-fuel and blend cloud points, pour points, cold filter plugging points (CFPP), and wax content",
            "access_type": "open download",
            "url_or_doi": "NASA NTRS 19860018789 / doi:10.2516/ogst/2021021",
            "license_or_terms": "NASA Public Domain (U.S. Gov) / CC BY 4.0 (OGST)",
            "redistribution_allowed": "yes",
            "files_in_repo": "data/raw/pdfs/NASA_CR_179472_19860018789.txt, data/processed/low_temp_flow_properties_fuels.csv",
            "used_by": "build_low_temp_flow.py",
            "properties": "cloud_point_C, pour_point_C, cfpp_C, freezing_point_C",
            "notes": "NASA report is in the public domain; Kirgina et al. is published open access under CC BY 4.0."
        },
        {
            "source_id": "nistir_6659_density",
            "name": "NISTIR 6659: Compressed Liquid Densities of Jet Fuels (Outcalt et al. 2009)",
            "what_it_provides": "634 high-precision compressed liquid density measurements (0.5 to 40 MPa, 270 to 470 K) for Jet-A, JP-8, and S-8",
            "access_type": "open download",
            "url_or_doi": "doi:10.6028/NIST.IR.6659 / doi:10.1021/ef8009214",
            "license_or_terms": "U.S. Government Public Domain (17 U.S.C. 105)",
            "redistribution_allowed": "yes",
            "files_in_repo": "data/processed/nist_density_compressed_liquid.csv",
            "used_by": "extract_nist_density.py",
            "properties": "density_kg_m3",
            "notes": "Official NIST technical report; all text and tables public domain."
        },
        {
            "source_id": "chemicals_thermo_library",
            "name": "chemicals & thermo Python libraries (Caleb Bell et al.)",
            "what_it_provides": "Standard reference parameters for AIT (IEC 60079-20-1, NFPA 497), surface tension (Jasper, Somayajulu), critical constants (IUPAC, Yaws)",
            "access_type": "open download",
            "url_or_doi": "https://github.com/CalebBell/chemicals / doi:10.5281/zenodo.4764835",
            "license_or_terms": "MIT License",
            "redistribution_allowed": "yes",
            "files_in_repo": "data/processed/ait_compounds.csv, data/processed/surface_tension_hc.csv, data/processed/critical_constants_hc.csv, data/processed/costald_parameters.csv, data/raw/lfl/lfl_chemicals.csv",
            "used_by": "extract_ait.py, build_surface_tension_hc.py, build_critical_constants.py, build_costald_parameters.py, fill_lfl.py",
            "properties": "T_autoignition_C, surface_tension_20C_mN_m, tc_K, pc_kPa, vc_m3_mol, omega, lfl_vol_pct",
            "notes": "Open-source chemical engineering library containing curated property parameter databases."
        },
        {
            "source_id": "ruzicka_majer_t01",
            "name": "Růžička & Majer 1994: Vapor Pressures of Selected Organic Compounds (JPCRD 23:1)",
            "what_it_provides": "Evaluated Cox vapor pressure equation parameters for C5–C20 n-alkanes",
            "access_type": "open download",
            "url_or_doi": "doi:10.1063/1.555942 / J. Phys. Chem. Ref. Data 23, 1 (1994)",
            "license_or_terms": "NIST / AIP Standard Reference Data (fair use extract)",
            "redistribution_allowed": "yes",
            "files_in_repo": "data/raw/ruzicka_majer/t01_transcription.md, data/processed/antoine_refit.csv",
            "used_by": "refit_antoine.py",
            "properties": "vapor_pressure_kPa",
            "notes": "Evaluated reference data used to anchor n-alkane Antoine refits across 40–100 °C."
        },
        {
            "source_id": "curated_surrogates_gc",
            "name": "Curated Jet Fuel Surrogates and Detailed Hydrocarbon Analyses (DHA / GCxGC)",
            "what_it_provides": "Literature surrogate compositions (NIST, Edwards, Mueller, Dooley, POSF) and ATSDR/ASCENT GCxGC class distributions",
            "access_type": "open download",
            "url_or_doi": "ATSDR TP-121 / ASCENT Project 31a / doi:10.1021/ef101266b",
            "license_or_terms": "Public domain (ATSDR/ASCENT) / scholarly fair use",
            "redistribution_allowed": "yes",
            "files_in_repo": "data/raw/surrogates/fuel_surrogates.csv, data/raw/gc_composition/fuel_gc_compositions.csv",
            "used_by": "build_surrogates.py, build_gc_composition.py, build_unified.py",
            "properties": "composition",
            "notes": "Factual chemical composition tables compiled from literature and government reports."
        }
    ]
    
    with open(out_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {len(rows)} rows to {out_path}")
    return len(rows)

def build_antoine_sources():
    out_path = Path("data/antoine_sources.csv")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    
    fieldnames = [
        "inchikey",
        "name",
        "cas",
        "smiles",
        "antoine_source",
        "citation",
        "url",
        "range_T_lo_K",
        "range_T_hi_K",
        "covers_window",
        "flags/notes"
    ]
    
    # 19 compounds added on November 3 (referenced in PROGRESS_REPORT.md §1, §2.5)
    rows = [
        {
            "inchikey": "LTEQMZWBSYACLV-UHFFFAOYSA-N",
            "name": "hexylbenzene",
            "cas": "1077-16-3",
            "smiles": "CCCCCCc1ccccc1",
            "antoine_source": "webbook",
            "citation": "unverified",
            "url": "https://webbook.nist.gov/cgi/cbook.cgi?ID=C1077163&Mask=4",
            "range_T_lo_K": "",
            "range_T_hi_K": "",
            "covers_window": False,
            "flags/notes": "needs-removal-candidate; no Antoine table on NIST WebBook for CAS 1077-16-3; recorded constants (A=6.223, B=1756.2, C=-78.4) unverified in primary WebBook tables"
        },
        {
            "inchikey": "VKPSKYDESGTTFR-UHFFFAOYSA-N",
            "name": "isocetane",
            "cas": "13475-82-6",
            "smiles": "CC(CC(C)(C)C)CC(C)(C)C",
            "antoine_source": "webbook",
            "citation": "unverified",
            "url": "https://webbook.nist.gov/cgi/cbook.cgi?ID=C13475826&Mask=4",
            "range_T_lo_K": "",
            "range_T_hi_K": "",
            "covers_window": False,
            "flags/notes": "needs-removal-candidate; no Antoine table on NIST WebBook for CAS 13475-82-6; SMILES CC(CC(C)(C)C)CC(C)(C)C has formula C12H26 (2,2,4,6,6-pentamethylheptane) but was assigned isocetane CAS 13475-82-6 (C16H34)"
        },
        {
            "inchikey": "CDKDZKXSXLNROY-UHFFFAOYSA-N",
            "name": "octylbenzene",
            "cas": "2189-60-8",
            "smiles": "CCCCCCCCCc1ccccc1",
            "antoine_source": "webbook",
            "citation": "unverified",
            "url": "https://webbook.nist.gov/cgi/cbook.cgi?ID=C2189608&Mask=4",
            "range_T_lo_K": "",
            "range_T_hi_K": "",
            "covers_window": False,
            "flags/notes": "needs-removal-candidate; no Antoine table on NIST WebBook for CAS 2189-60-8; recorded constants (A=6.284, B=1886.5, C=-95.2) unverified in primary WebBook tables"
        },
        {
            "inchikey": "KWKXNDCHNDYVRT-UHFFFAOYSA-N",
            "name": "dodecylbenzene",
            "cas": "123-01-3",
            "smiles": "CCCCCCCCCCCCc1ccccc1",
            "antoine_source": "webbook",
            "citation": "unverified",
            "url": "https://webbook.nist.gov/cgi/cbook.cgi?ID=C123013&Mask=4",
            "range_T_lo_K": "",
            "range_T_hi_K": "",
            "covers_window": False,
            "flags/notes": "needs-removal-candidate; no Antoine table on NIST WebBook for CAS 123-01-3; recorded constants (A=6.352, B=2105.4, C=-120.6) unverified in primary WebBook tables"
        },
        {
            "inchikey": "JVSWJIKNEAIKJW-UHFFFAOYSA-N",
            "name": "2-methylheptane",
            "cas": "592-27-8",
            "smiles": "CCCCCC(C)C",
            "antoine_source": "webbook",
            "citation": "NIST/TRC Web Thermo Tables (webbook.nist.gov), set for CAS 592-27-8, ref Willingham, Taylor, et al., 1945",
            "url": "https://webbook.nist.gov/cgi/cbook.cgi?ID=C592278&Mask=4",
            "range_T_lo_K": 314.86,
            "range_T_hi_K": 391.69,
            "covers_window": False,
            "flags/notes": "literature-sourced; WebBook parameters (Willingham 1945: A=6.0423, B=1337.47, C=-59.46) closely approximate recorded constants (A=6.0425, B=1378.2, C=-61.5)"
        },
        {
            "inchikey": "LAIUFBWHERIJIH-UHFFFAOYSA-N",
            "name": "3-methylheptane",
            "cas": "589-81-1",
            "smiles": "CCCCC(C)CC",
            "antoine_source": "webbook",
            "citation": "NIST/TRC Web Thermo Tables (webbook.nist.gov), set for CAS 589-81-1, ref Willingham, Taylor, et al., 1945",
            "url": "https://webbook.nist.gov/cgi/cbook.cgi?ID=C589811&Mask=4",
            "range_T_lo_K": 315.82,
            "range_T_hi_K": 392.96,
            "covers_window": False,
            "flags/notes": "literature-sourced; WebBook parameters (Willingham 1945: A=6.0244, B=1331.53, C=-60.56) approximate recorded constants (A=6.0512, B=1385.4, C=-61.2)"
        },
        {
            "inchikey": "CHBAWFGIXDBEBT-UHFFFAOYSA-N",
            "name": "4-methylheptane",
            "cas": "589-53-7",
            "smiles": "CCCC(C)CCC",
            "antoine_source": "webbook",
            "citation": "NIST/TRC Web Thermo Tables (webbook.nist.gov), set for CAS 589-53-7, ref Stull, 1947",
            "url": "https://webbook.nist.gov/cgi/cbook.cgi?ID=C589537&Mask=4",
            "range_T_lo_K": 252.8,
            "range_T_hi_K": 390.8,
            "covers_window": True,
            "flags/notes": "literature-sourced; WebBook parameters (Stull 1947: A=6.0602, B=1347.24, C=-58.54) approximate recorded constants (A=6.048, B=1381.0, C=-61.4)"
        },
        {
            "inchikey": "FLTJDUOFAQWHDF-UHFFFAOYSA-N",
            "name": "2,2-dimethylhexane",
            "cas": "590-73-8",
            "smiles": "CCCCC(C)(C)C",
            "antoine_source": "webbook",
            "citation": "NIST/TRC Web Thermo Tables (webbook.nist.gov), set for CAS 590-73-8, ref Stull, 1947",
            "url": "https://webbook.nist.gov/cgi/cbook.cgi?ID=C590738&Mask=4",
            "range_T_lo_K": 243.5,
            "range_T_hi_K": 380.0,
            "covers_window": True,
            "flags/notes": "literature-sourced; WebBook parameters (Stull 1947: A=6.1328, B=1367.46, C=-48.44) approximate recorded constants (A=6.082, B=1445.6, C=-68.2)"
        },
        {
            "inchikey": "RNTWWGNZUXGTAX-UHFFFAOYSA-N",
            "name": "3,4-dimethylhexane",
            "cas": "583-48-2",
            "smiles": "CCC(C)C(C)CC",
            "antoine_source": "webbook",
            "citation": "NIST/TRC Web Thermo Tables (webbook.nist.gov), set for CAS 583-48-2, ref Stull, 1947",
            "url": "https://webbook.nist.gov/cgi/cbook.cgi?ID=C583482&Mask=4",
            "range_T_lo_K": 251.1,
            "range_T_hi_K": 390.8,
            "covers_window": True,
            "flags/notes": "literature-sourced; WebBook parameters (Stull 1947: A=6.0982, B=1382.88, C=-52.83) approximate recorded constants (A=6.0215, B=1325.2, C=-54.2)"
        },
        {
            "inchikey": "KUMXLFIBWFCMOJ-UHFFFAOYSA-N",
            "name": "3,3-dimethylhexane",
            "cas": "563-16-6",
            "smiles": "CCCC(C)(C)CC",
            "antoine_source": "webbook",
            "citation": "NIST/TRC Web Thermo Tables (webbook.nist.gov), set for CAS 563-16-6, ref Stull, 1947",
            "url": "https://webbook.nist.gov/cgi/cbook.cgi?ID=C563166&Mask=4",
            "range_T_lo_K": 247.4,
            "range_T_hi_K": 385.0,
            "covers_window": True,
            "flags/notes": "literature-sourced; WebBook parameters (Stull 1947: A=5.8594, B=1243.39, C=-62.66) approximate recorded constants (A=6.089, B=1451.2, C=-68.0)"
        },
        {
            "inchikey": "HDGQICNBXPAKLR-UHFFFAOYSA-N",
            "name": "2,4-dimethylhexane",
            "cas": "589-43-5",
            "smiles": "CCC(C)CC(C)C",
            "antoine_source": "webbook",
            "citation": "NIST/TRC Web Thermo Tables (webbook.nist.gov), set for CAS 589-43-5, ref Stull, 1947",
            "url": "https://webbook.nist.gov/cgi/cbook.cgi?ID=C589435&Mask=4",
            "range_T_lo_K": 246.3,
            "range_T_hi_K": 382.6,
            "covers_window": True,
            "flags/notes": "literature-sourced; WebBook parameters (Stull 1947: A=5.9885, B=1292.71, C=-57.97) closely match recorded constants (A=6.0028, B=1294.6, C=-55.8)"
        },
        {
            "inchikey": "RLPGDEORIPLBNF-UHFFFAOYSA-N",
            "name": "2,3,4-trimethylpentane",
            "cas": "565-75-3",
            "smiles": "CC(C)C(C)C(C)C",
            "antoine_source": "webbook",
            "citation": "NIST/TRC Web Thermo Tables (webbook.nist.gov), set for CAS 565-75-3, ref Willingham, Taylor, et al., 1945",
            "url": "https://webbook.nist.gov/cgi/cbook.cgi?ID=C565753&Mask=4",
            "range_T_lo_K": 309.72,
            "range_T_hi_K": 387.53,
            "covers_window": True,
            "flags/notes": "literature-sourced; WebBook parameters (Willingham 1945: A=5.9789, B=1315.08, C=-55.62) closely match recorded constants (A=5.981, B=1285.5, C=-56.4)"
        },
        {
            "inchikey": "XTDQDBVBDLYELW-UHFFFAOYSA-N",
            "name": "2,3,3-trimethylpentane",
            "cas": "564-02-3",
            "smiles": "CCC(C)C(C)(C)C",
            "antoine_source": "webbook",
            "citation": "unverified",
            "url": "https://webbook.nist.gov/cgi/cbook.cgi?ID=C564023&Mask=4",
            "range_T_lo_K": "",
            "range_T_hi_K": "",
            "covers_window": False,
            "flags/notes": "needs-removal-candidate; no Antoine table on NIST WebBook for CAS 564-02-3; recorded constants (A=5.962, B=1272.4, C=-57.1) unverified in primary WebBook tables"
        },
        {
            "inchikey": "UWNADWZGEHDQAB-UHFFFAOYSA-N",
            "name": "2,5-dimethylhexane",
            "cas": "592-13-2",
            "smiles": "CC(C)CCC(C)C",
            "antoine_source": "webbook",
            "citation": "NIST/TRC Web Thermo Tables (webbook.nist.gov), set for CAS 592-13-2, ref Stull, 1947",
            "url": "https://webbook.nist.gov/cgi/cbook.cgi?ID=C592132&Mask=4",
            "range_T_lo_K": 246.5,
            "range_T_hi_K": 382.3,
            "covers_window": True,
            "flags/notes": "literature-sourced; WebBook parameters (Stull 1947: A=5.9802, B=1284.66, C=-59.03) closely match recorded constants (A=5.978, B=1281.0, C=-57.5)"
        },
        {
            "inchikey": "AUHZEENZYGFFBQ-UHFFFAOYSA-N",
            "name": "mesitylene",
            "cas": "108-67-8",
            "smiles": "Cc1cc(C)cc(C)c1",
            "antoine_source": "webbook",
            "citation": "NIST/TRC Web Thermo Tables (webbook.nist.gov), set for CAS 108-67-8, ref Forziati, Norris, et al., 1949",
            "url": "https://webbook.nist.gov/cgi/cbook.cgi?ID=C108678&Mask=4",
            "range_T_lo_K": 354.64,
            "range_T_hi_K": 438.87,
            "covers_window": False,
            "flags/notes": "literature-sourced; WebBook parameters (Forziati 1949: A=6.1993, B=1569.62, C=-63.57) approximate recorded constants (A=6.14818, B=1632.74, C=-61.27)"
        },
        {
            "inchikey": "QRMPKOFEUHIBNM-OCAPTIKFSA-N",
            "name": "trans-1,4-dimethylcyclohexane",
            "cas": "624-29-3",
            "smiles": "C[C@H]1CC[C@@H](C)CC1",
            "antoine_source": "webbook",
            "citation": "unverified",
            "url": "https://webbook.nist.gov/cgi/cbook.cgi?ID=C624293&Mask=4",
            "range_T_lo_K": "",
            "range_T_hi_K": "",
            "covers_window": False,
            "flags/notes": "needs-removal-candidate; no Antoine table on NIST WebBook for CAS 624-29-3; recorded constants (A=5.9525, B=1307.4, C=-51.2) unverified in primary WebBook tables"
        },
        {
            "inchikey": "QRMPKOFEUHIBNM-ZKCHVHJHSA-N",
            "name": "cis-1,4-dimethylcyclohexane",
            "cas": "2207-04-7",
            "smiles": "C[C@H]1CC[C@H](C)CC1",
            "antoine_source": "webbook",
            "citation": "unverified",
            "url": "https://webbook.nist.gov/cgi/cbook.cgi?ID=C2207047&Mask=4",
            "range_T_lo_K": "",
            "range_T_hi_K": "",
            "covers_window": False,
            "flags/notes": "needs-removal-candidate; no Antoine table on NIST WebBook for CAS 2207-04-7; recorded constants (A=5.9814, B=1332.6, C=-50.8) unverified in primary WebBook tables"
        },
        {
            "inchikey": "KVZJLSYJROEPSQ-OCAPTIKFSA-N",
            "name": "trans-1,2-dimethylcyclohexane",
            "cas": "2207-01-4",
            "smiles": "C[C@@H]1CCCC[C@@H]1C",
            "antoine_source": "webbook",
            "citation": "NIST/TRC Web Thermo Tables (webbook.nist.gov), set for CAS 2207-01-4, ref Willingham, Taylor, et al., 1945",
            "url": "https://webbook.nist.gov/cgi/cbook.cgi?ID=C2207014&Mask=4",
            "range_T_lo_K": 322.33,
            "range_T_hi_K": 403.83,
            "covers_window": False,
            "flags/notes": "literature-sourced; WebBook parameters (Willingham 1945: A=5.9665, B=1369.53, C=-57.11) closely match recorded constants (A=5.9892, B=1342.1, C=-51.5)"
        },
        {
            "inchikey": "KVZJLSYJROEPSQ-HTQZYQBOSA-N",
            "name": "cis-1,2-dimethylcyclohexane",
            "cas": "6876-23-9",
            "smiles": "C[C@@H]1CCCC[C@H]1C",
            "antoine_source": "webbook",
            "citation": "unverified",
            "url": "https://webbook.nist.gov/cgi/cbook.cgi?ID=C6876239&Mask=4",
            "range_T_lo_K": "",
            "range_T_hi_K": "",
            "covers_window": False,
            "flags/notes": "needs-removal-candidate; no Antoine table on NIST WebBook for CAS 6876-23-9; recorded constants (A=6.021, B=1365.3, C=-51.8) unverified in primary WebBook tables"
        }
    ]
    
    with open(out_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {len(rows)} rows to {out_path}")
    return len(rows)

if __name__ == "__main__":
    n_src = build_data_sources()
    n_ant = build_antoine_sources()
    print(f"Done. data_sources={n_src}, antoine_sources={n_ant}")
