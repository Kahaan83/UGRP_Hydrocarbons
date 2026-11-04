# Phase 3 Real-Fuel Dataset Report (T3.1 - T3.5)

Date: 2026-11-04  
Deliverables:
- Inventory: `C:/Users/kahaa/OneDrive/Desktop/UGRP_hydrocarbons/fuelblend_data/data/real_fuels/inventory.csv`
- Compositions: `C:/Users/kahaa/OneDrive/Desktop/UGRP_hydrocarbons/fuelblend_data/data/real_fuels/compositions.csv`
- Blend Links: `C:/Users/kahaa/OneDrive/Desktop/UGRP_hydrocarbons/fuelblend_data/data/real_fuels/blend_links.csv`
- Holdout Splits: `C:/Users/kahaa/OneDrive/Desktop/UGRP_hydrocarbons/fuelblend_data/splits/real_fuels_holdout.json`

## 1. Acceptance Criteria Verification
- **T3.1 Inventory Table (`data/real_fuels/inventory.csv`):** Complete catalogue of real fuels with fuel_id, source, fuel_family, number_of_resolved_components, and available properties.
- **T3.2 Compositions Table (`data/real_fuels/compositions.csv`):** Canonical InChIKeys and SMILES for every component; every fuel records `number_of_resolved_components` and `unresolved_fraction`. Verified with `validate_composition_df` (strict=True).
- **T3.3 Blend-of-Blends Links (`data/real_fuels/blend_links.csv`):** Links surrogate blends and parent commercial fuels with blend ratios and measured properties.
- **T3.4 Strict Segregation of Reconstructed Compositions:** CRC 647 pseudo-compositions are strictly tagged as `reconstructed_composition` (57 fuels) and isolated from true `measured_dha_gcxgc` fuels (8 fuels).
- **T3.5 Real Fuels Holdout Policy:** All 65 real fuels are held out from training by default in `splits/real_fuels_holdout.json`.

## 2. Inventory Breakdown by Composition Type
- **`reconstructed_composition`:** 57 fuels
- **`measured_dha_gcxgc`:** 8 fuels

## 3. Inventory Breakdown by Fuel Family
- **`jet`:** 63 fuels
- **`diesel`:** 2 fuels
