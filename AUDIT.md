# Audit

FAIL: 0 | MISSING: 0 | WARN: 7 | OK: 123

## Needs attention

- **WARN** [env] project is inside OneDrive: sync of ~800 small files + 400 MB slows everything; move to a non-synced folder
- **WARN** [interim] inspect read errors: EmptyDataError:2, ValueError:1
- **WARN** [processed] Antoine compounds found: 22/24; missing: ['methylcyclohexane', 'n-eicosane']
- **WARN** [processed] Antoine range covers 40-100 C: 20 compounds do NOT: ['n-hexane', 'n-heptane', 'n-octane', 'n-nonane', 'n-decane', 'n-undecane', 'n-dodecane', 'n-tridecane', 'n-tetradecane', 'n-pentadecane', 'n-hexadecane', 'n-octadecane', 'cyclohexane', 'benzene', 'ethylbenzene', 'n-propylbenzene', 'n-butylbenzene', '1,2,4-trimethylbenzene', 'tetralin', 'naphthalene'] -> refit from vapour-pressure data
- **WARN** [processed] WebBook page name vs your name (heuristic): check by eye: ['n-hexane->NIST Chemistry\n  WebBook, SRD 69', 'n-heptane->NIST Chemistry\n  WebBook, SRD 69', 'n-octane->NIST Chemistry\n  WebBook, SRD 69', 'n-nonane->NIST Chemistry\n  WebBook, SRD 69', 'n-decane->NIST Chemistry\n  WebBook, SRD 69', 'n-undecane->NIST Chemistry\n  WebBook, SRD 69', 'n-dodecane->NIST Chemistry\n  WebBook, SRD 69', 'n-tridecane->NIST Chemistry\n  WebBook, SRD 69', 'n-tetradecane->NIST Chemistry\n  WebBook, SRD 69', 'n-pentadecane->NIST Chemistry\n  WebBook, SRD 69', 'n-hexadecane->NIST Chemistry\n  WebBook, SRD 69', 'n-octadecane->NIST Chemistry\n  WebBook, SRD 69']
- **WARN** [raw] pdf:AFRL_TR-2020-0017_jet_fuel_properties.pdf: optional - DTIC site; get later in a browser. save into data/raw/pdfs/ with this exact name
- **WARN** [raw] pdf:CRC_handbook_aviation_fuel_properties.pdf: optional - DTIC site; get later in a browser. save into data/raw/pdfs/ with this exact name

## Everything

| status | area | item | detail |
|---|---|---|---|
| OK | env | python 3.12.5 on Windows |  |
| OK | env | pandas | 3.0.6 |
| OK | env | numpy | 2.5.3 |
| OK | env | rdkit | 2026.03.6 |
| OK | env | requests | 2.34.2 |
| OK | env | lxml | 6.1.3 |
| OK | env | openpyxl | 3.1.5 |
| WARN | env | project is inside OneDrive | sync of ~800 small files + 400 MB slows everything; move to a non-synced folder |
| OK | raw | git:chemixhub |  |
| OK | raw | git:mixture_evaluation |  |
| OK | raw | CheMixHub nist-logV table | data\raw\git\chemixhub\datasets\nist-logV\processed_data\processed_NISTlogV.csv |
| OK | raw | zenodo:leenhouts_benchmark/datasets.zip | 123 MB |
| OK | raw | zenodo:leenhouts_benchmark/splits.zip | 201 MB |
| OK | raw | zenodo:bilodeau_nist_viscosity/pretrained_models.zip | 95 MB |
| OK | raw | pdf:NISTIR_6646_RP1.pdf | 5.1 MB |
| OK | raw | pdf:CRC_647_world_fuel_sampling.pdf | 1.3 MB |
| OK | raw | pdf:CRC_AV-27-18_fuel_properties.pdf | 1.7 MB |
| OK | raw | pdf:CRC_AV-11-09_freezing_point.pdf | 0.7 MB |
| WARN | raw | pdf:AFRL_TR-2020-0017_jet_fuel_properties.pdf | optional - DTIC site; get later in a browser. save into data/raw/pdfs/ with this exact name |
| WARN | raw | pdf:CRC_handbook_aviation_fuel_properties.pdf | optional - DTIC site; get later in a browser. save into data/raw/pdfs/ with this exact name |
| OK | raw | pdf:Zabetakis_Bulletin_627_flammability.pdf | 6.9 MB |
| OK | raw | pdf:Li_2020_flash_point_GNN_arXiv2002.11315.pdf | 0.7 MB |
| OK | raw | pdf:Leenhouts_2026_mixture_eval_arXiv2605.29698.pdf | 4.7 MB |
| OK | raw | pdf:Chew_2024_viscosity_SI.pdf | 2.0 MB |
| OK | raw | ThermoML ThermoML.v2020-09-30.tgz | 189 MB |
| OK | interim | flashpoint_pure_exp.csv | 585 rows |
| OK | interim | flashpoint_mix_exp.csv | 1006 rows |
| OK | interim | viscosity_pure_exp.csv | 3027 rows |
| OK | interim | viscosity_mix_exp.csv | 35938 rows |
| OK | interim | dcn_mix_exp.csv | 484 rows |
| OK | interim | mon_mix_exp.csv | 292 rows |
| OK | interim | NIST_Visc_Data.csv | 239200 rows |
| WARN | interim | inspect read errors | EmptyDataError:2, ValueError:1 |
| OK | processed | NIST viscosity (pure, hydrocarbon) | 4832 rows |
| OK | processed | NIST viscosity (mixtures, hydrocarbon) | 10947 rows |
| OK | processed | hydrocarbon compound table | 109 rows |
| OK | processed | thermoml_hc.csv | 46747 rows; props: 4; Kinematic viscosity, m2/s/1c=78; Kinematic viscosity, m2/s/2c=28; Mass density, kg/m3/1c=6651; Mass density, kg/m3/2c=13164; Mass density, kg/m3/3c=911; Vapor or sublimation pressur/1c=8017; Vapor or sublimation pressur/2c=5570; Vapor or sublimation pressur/3c=174; Viscosity, Pa*s/1c=3352; Viscosity, Pa*s/2c=8177; Viscosity, Pa*s/3c=625 |
| OK | processed | thermoml: rows with temperature | 100% |
| WARN | processed | Antoine compounds found | 22/24; missing: ['methylcyclohexane', 'n-eicosane'] |
| WARN | processed | Antoine range covers 40-100 C | 20 compounds do NOT: ['n-hexane', 'n-heptane', 'n-octane', 'n-nonane', 'n-decane', 'n-undecane', 'n-dodecane', 'n-tridecane', 'n-tetradecane', 'n-pentadecane', 'n-hexadecane', 'n-octadecane', 'cyclohexane', 'benzene', 'ethylbenzene', 'n-propylbenzene', 'n-butylbenzene', '1,2,4-trimethylbenzene', 'tetralin', 'naphthalene'] -> refit from vapour-pressure data |
| WARN | processed | WebBook page name vs your name (heuristic) | check by eye: ['n-hexane->NIST Chemistry\n  WebBook, SRD 69', 'n-heptane->NIST Chemistry\n  WebBook, SRD 69', 'n-octane->NIST Chemistry\n  WebBook, SRD 69', 'n-nonane->NIST Chemistry\n  WebBook, SRD 69', 'n-decane->NIST Chemistry\n  WebBook, SRD 69', 'n-undecane->NIST Chemistry\n  WebBook, SRD 69', 'n-dodecane->NIST Chemistry\n  WebBook, SRD 69', 'n-tridecane->NIST Chemistry\n  WebBook, SRD 69', 'n-tetradecane->NIST Chemistry\n  WebBook, SRD 69', 'n-pentadecane->NIST Chemistry\n  WebBook, SRD 69', 'n-hexadecane->NIST Chemistry\n  WebBook, SRD 69', 'n-octadecane->NIST Chemistry\n  WebBook, SRD 69'] |
| OK | coverage | flash point, pure | 244/585 rows all-hydrocarbon; n_comp among HC rows: {1: 244}; value -136 to 204 Celcius |
| OK | coverage | flash point, mixtures | 126/1006 rows all-hydrocarbon; n_comp among HC rows: {2: 123, 3: 3}; value -4.47 to 156 Celcius |
| OK | coverage | viscosity, pure | 692/3027 rows all-hydrocarbon; n_comp among HC rows: {1: 692}; T 293-323 K; value -4.89 to 3.44 ln(mPa*s) |
| OK | coverage | viscosity, mixtures | 2752/35938 rows all-hydrocarbon; n_comp among HC rows: {2: 2752}; T 293-323 K; value -1.54 to 3.21 ln(mPa*s) |
| OK | coverage | DCN, mixtures | 346/484 rows all-hydrocarbon; n_comp among HC rows: {2: 213, 3: 74, 4: 59}; value 7.25 to 102 - |
| OK | coverage | MON, mixtures | 106/292 rows all-hydrocarbon; n_comp among HC rows: {2: 30, 3: 54, 4: 15, 8: 1, 23: 1, 27: 1, 30: 1, 31: 1, 32: 1, 34: 1}; value 0 to 101 - |
| OK | unified | HC records: acid_number_mg_koh_g | 53 ({'real_fuel': 53}); T nan-nan K; systems=53 |
| OK | unified | HC records: aromatics_vol_pct | 77 ({'real_fuel': 77}); T nan-nan K; systems=77 |
| OK | unified | HC records: benzothiophenes_ppm_wt | 57 ({'real_fuel': 57}); T nan-nan K; systems=57 |
| OK | unified | HC records: cloud_point_C | 3 ({'real_fuel': 3}); T nan-nan K; systems=3 |
| OK | unified | HC records: conductivity_ps_m | 171 ({'real_fuel': 171}); T 279-311 K; systems=57 |
| OK | unified | HC records: d86_residue_pct | 57 ({'real_fuel': 57}); T nan-nan K; systems=57 |
| OK | unified | HC records: dcn | 426 ({'binary': 213, 'multi': 133, 'pure': 80}); T nan-nan K; systems=118 |
| OK | unified | HC records: density_kg_m3 | 18452 ({'binary': 12747, 'real_fuel': 3056, 'pure': 1752, 'multi': 897}); T -40-468 K; systems=2663 |
| OK | unified | HC records: dibenzothiophenes_ppm_wt | 57 ({'real_fuel': 57}); T nan-nan K; systems=57 |
| OK | unified | HC records: dielectric_constant | 340 ({'real_fuel': 340}); T 234-344 K; systems=57 |
| OK | unified | HC records: distillation_T0.025_C | 12 ({'real_fuel': 12}); T nan-nan K; systems=12 |
| OK | unified | HC records: distillation_T0_C | 1760 ({'real_fuel': 1760}); T nan-nan K; systems=1760 |
| OK | unified | HC records: distillation_T100_C | 2239 ({'real_fuel': 2239}); T nan-nan K; systems=2239 |
| OK | unified | HC records: distillation_T10_C | 2252 ({'real_fuel': 2252}); T nan-nan K; systems=2252 |
| OK | unified | HC records: distillation_T15_C | 33 ({'real_fuel': 33}); T nan-nan K; systems=33 |
| OK | unified | HC records: distillation_T20_C | 707 ({'real_fuel': 707}); T nan-nan K; systems=707 |
| OK | unified | HC records: distillation_T25_C | 33 ({'real_fuel': 33}); T nan-nan K; systems=33 |
| OK | unified | HC records: distillation_T30_C | 94 ({'real_fuel': 94}); T nan-nan K; systems=94 |
| OK | unified | HC records: distillation_T35_C | 33 ({'real_fuel': 33}); T nan-nan K; systems=33 |
| OK | unified | HC records: distillation_T40_C | 94 ({'real_fuel': 94}); T nan-nan K; systems=94 |
| OK | unified | HC records: distillation_T45_C | 33 ({'real_fuel': 33}); T nan-nan K; systems=33 |
| OK | unified | HC records: distillation_T50_C | 2254 ({'real_fuel': 2254}); T nan-nan K; systems=2254 |
| OK | unified | HC records: distillation_T55_C | 33 ({'real_fuel': 33}); T nan-nan K; systems=33 |
| OK | unified | HC records: distillation_T5_C | 264 ({'real_fuel': 264}); T nan-nan K; systems=264 |
| OK | unified | HC records: distillation_T60_C | 94 ({'real_fuel': 94}); T nan-nan K; systems=94 |
| OK | unified | HC records: distillation_T65_C | 33 ({'real_fuel': 33}); T nan-nan K; systems=33 |
| OK | unified | HC records: distillation_T70_C | 94 ({'real_fuel': 94}); T nan-nan K; systems=94 |
| OK | unified | HC records: distillation_T75_C | 33 ({'real_fuel': 33}); T nan-nan K; systems=33 |
| OK | unified | HC records: distillation_T80_C | 94 ({'real_fuel': 94}); T nan-nan K; systems=94 |
| OK | unified | HC records: distillation_T85_C | 34 ({'real_fuel': 34}); T nan-nan K; systems=34 |
| OK | unified | HC records: distillation_T90_C | 2239 ({'real_fuel': 2239}); T nan-nan K; systems=2239 |
| OK | unified | HC records: distillation_T95_C | 250 ({'real_fuel': 250}); T nan-nan K; systems=250 |
| OK | unified | HC records: distillation_loss_pct | 2504 ({'real_fuel': 2504}); T nan-nan K; systems=2504 |
| OK | unified | HC records: distillation_recovery_pct | 804 ({'real_fuel': 804}); T nan-nan K; systems=804 |
| OK | unified | HC records: distillation_residue_pct | 2512 ({'real_fuel': 2512}); T nan-nan K; systems=2512 |
| OK | unified | HC records: flash_point_C | 3264 ({'real_fuel': 2349, 'pure': 789, 'binary': 123, 'multi': 3}); T nan-nan K; systems=3169 |
| OK | unified | HC records: freezing_point_C | 2303 ({'real_fuel': 2278, 'pure': 25}); T nan-nan K; systems=2303 |
| OK | unified | HC records: gas_solubility_air_ppm_wt | 18 ({'real_fuel': 18}); T 293-293 K; systems=18 |
| OK | unified | HC records: gas_solubility_ar_ppm_wt | 18 ({'real_fuel': 18}); T 293-293 K; systems=18 |
| OK | unified | HC records: gas_solubility_n2_ppm_wt | 18 ({'real_fuel': 18}); T 293-293 K; systems=18 |
| OK | unified | HC records: gas_solubility_o2_ppm_wt | 18 ({'real_fuel': 18}); T 293-293 K; systems=18 |
| OK | unified | HC records: hydrogen_wt_pct | 57 ({'real_fuel': 57}); T nan-nan K; systems=57 |
| OK | unified | HC records: lubricity_wear_scar_mm | 57 ({'real_fuel': 57}); T nan-nan K; systems=57 |
| OK | unified | HC records: mon | 373 ({'pure': 267, 'multi': 76, 'binary': 30}); T nan-nan K; systems=289 |
| OK | unified | HC records: naphthalenes_vol_pct | 57 ({'real_fuel': 57}); T nan-nan K; systems=57 |
| OK | unified | HC records: net_heat_content_btu_lb | 57 ({'real_fuel': 57}); T nan-nan K; systems=57 |
| OK | unified | HC records: net_heat_content_calculated_method_btu_lb | 57 ({'real_fuel': 57}); T nan-nan K; systems=57 |
| OK | unified | HC records: nitrogen_solubility_ppm_wt | 19 ({'real_fuel': 19}); T 293-293 K; systems=19 |
| OK | unified | HC records: olefins_vol_pct | 57 ({'real_fuel': 57}); T nan-nan K; systems=57 |
| OK | unified | HC records: polar_species_final_0.6min_elution_mg_l | 19 ({'real_fuel': 19}); T nan-nan K; systems=19 |
| OK | unified | HC records: polar_species_total_mg_l | 20 ({'real_fuel': 20}); T nan-nan K; systems=20 |
| OK | unified | HC records: pour_point_C | 60 ({'real_fuel': 60}); T nan-nan K; systems=60 |
| OK | unified | HC records: refractive_index | 57 ({'real_fuel': 57}); T 293-294 K; systems=57 |
| OK | unified | HC records: saturates_vol_pct | 57 ({'real_fuel': 57}); T nan-nan K; systems=57 |
| OK | unified | HC records: smoke_point_mm | 2197 ({'real_fuel': 2197}); T nan-nan K; systems=2197 |
| OK | unified | HC records: specific_heat_capacity_j_g_k | 50 ({'real_fuel': 50}); T 273-273 K; systems=50 |
| OK | unified | HC records: specific_heat_capacity_slope | 50 ({'real_fuel': 50}); T nan-nan K; systems=50 |
| OK | unified | HC records: speed_of_sound_m_s | 170 ({'real_fuel': 170}); T 234-278 K; systems=57 |
| OK | unified | HC records: sulfur_mercaptan_wt_ppm | 57 ({'real_fuel': 57}); T nan-nan K; systems=57 |
| OK | unified | HC records: sulfur_total_by_hplc_class_ppm_wt | 57 ({'real_fuel': 57}); T nan-nan K; systems=57 |
| OK | unified | HC records: sulfur_wt_ppm | 57 ({'real_fuel': 57}); T nan-nan K; systems=57 |
| OK | unified | HC records: surface_tension_mn_m | 135 ({'real_fuel': 135}); T 263-313 K; systems=54 |
| OK | unified | HC records: thermal_stability_break_point_C | 53 ({'real_fuel': 53}); T nan-nan K; systems=53 |
| OK | unified | HC records: thermal_stability_spent_fuel_mL | 176 ({'real_fuel': 176}); T nan-nan K; systems=176 |
| OK | unified | HC records: thiols_sulfides_disulfides_ppm_wt | 57 ({'real_fuel': 57}); T nan-nan K; systems=57 |
| OK | unified | HC records: thiophenes_ppm_wt | 57 ({'real_fuel': 57}); T nan-nan K; systems=57 |
| OK | unified | HC records: vapor_pressure_d2879_kpa | 76 ({'real_fuel': 76}); T 298-339 K; systems=19 |
| OK | unified | HC records: vapor_pressure_d6378_kpa | 76 ({'real_fuel': 76}); T 298-339 K; systems=19 |
| OK | unified | HC records: vapor_pressure_kPa | 6500 ({'pure': 4415, 'binary': 2021, 'multi': 64}); T 230-520 K; systems=268 |
| OK | unified | HC records: viscosity_dyn_mPas | 18932 ({'binary': 16037, 'pure': 2270, 'multi': 625}); T 230-453 K; systems=812 |
| OK | unified | HC records: viscosity_kin_cSt | 2664 ({'real_fuel': 2569, 'pure': 67, 'binary': 28}); T -40-448 K; systems=2165 |
| OK | unified | multi-component (>=3) hydrocarbon records | 1798 |
| OK | unified | conflicting duplicate records | 1199/76175 |
| OK | unified | compounds with Antoine valid at 40-100 C | 81/1038 |
| OK | unified | compounds with LFL | 695/1038 (measured 25, estimated 670)  |
| OK | unified | rows removed by fuel-domain filter | 5327 (see measurements_out_of_domain.csv) |
| OK | unified | flash-point physics baseline report |  |
| OK | unified | real-fuel records (AJFTD) | 33513 |
| OK | manual | AJFTD jet-fuel export (needs registration) | data/raw/ajftd/ |
| OK | manual | NIST RP-1 11-blend data | data/raw/rp1/ |
| OK | manual | pour/cloud point tables (hand-extracted) | data/raw/pour_point/ |
| OK | manual | LFL values (from Zabetakis) | data/raw/lfl/ |
| OK | manual | Li et al. flash-point SI | data/raw/flashpoint_li/ |