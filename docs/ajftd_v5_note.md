# AJFTD v5 bulk-descriptor load (task D1)

Loader: `src/fuelblend/data/loaders_ajftd.py` (`python -m fuelblend.data.loaders_ajftd`, PYTHONPATH=src).
Outputs: `data/real_fuels/v5/ajftd_fuels_{long,wide}.csv`, `ajftd_build_log.json`. Provenance class `real_fuel_bulk_descriptors`.

Viscosity temperature: `temperature_value` (text deg C; -20 dominant, also -40, 40, 20, 90...). About 590 numeric rows
(methods `310`, `D445`, spec limit `8 Max`) have blank/"No Test" temperature; they are the -20 C spec-table entry and are loaded at
253.15 K with `T_source=assumed_spec_-20C`. Method `310` is a spec-table id, not an ASTM number.

Raw-data quirks handled: `density_api_by_meter` rows with Unit kg/m3 are really net heat of combustion in kJ/kg; some heat
values are 0.1004 (dropped as implausible); density sometimes in g/mL; distillation rows with blank method kept; "time_stamp" strings dropped.

Composition Metrics (29,807 rows, 15,315 numeric): NO hydrocarbon-type composition suitable for blend modelling. Bulk class numbers only:
aromatics_total (2,226 fuels), diaromatics (1,574), olefins (159), saturates (23), paraffins_total (54), cycloalkanes_total (51);
paraffins+cycloalkanes+aromatics together on 51 fuels, saturates+aromatics on 23. No iso/n-paraffin split, no carbon-number data;
`gas_chromatography` entries are free text. Labels JETSCREEN (2 fuels), CRC (3), HR (71) carry the same bulk descriptors only.
