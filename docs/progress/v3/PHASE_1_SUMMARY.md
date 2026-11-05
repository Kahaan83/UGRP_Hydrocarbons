# Phase 1 Summary: Real-Fuel Data Expansion (Tasks S1–S6)

## Overview
Phase 1 expands real-fuel speciation and property measurements across Gasoline (S1), Diesel (S2), and Jet (S3) families, enforcing strict zero-leakage held-out splits and registering first-class oxygenate components.

## Fuel Inventory Counts (Total: 153 fuels)
- **Gasoline (S1):** 55 fuels (target $\ge 30$). Includes CRC FACE Gasolines A–J, 24 FACE ethanol blends (10–60 vol%), 9 PRF ethanol blends, Co-Optima core research gasolines (4-comp, E10–E30, M10–M20, TIGAS, NREL BOB-01), and standard TRFs.
- **Diesel (S2):** 16 fuels (target $\ge 15$). Includes CRC FACE Diesels 1–9, NASA paraffinic diesels (LFP3, LFP5), Sandia/Mueller 4-comp, IDEAL 4-comp, Peugeot 4-comp, European EN 590, and ASTM D975 ULSD. Clearly segregates DCN (IQT/AFIDA) vs CFR engine Cetane.
- **Jet (S3):** 82 fuels. Includes POSF fuels (10325, 10264, 10289, 11498, 12344, 12345, 13511, 13516), NIST RP-1, 57 CRC 647 survey fuels, and literature combustion surrogates.

## Deliverables Generated
1. `data/real_fuels/v3/inventory_v3.csv` (153 fuels with speciation metadata and measured properties).
2. `data/real_fuels/v3/compositions_v3.csv` (3298 component rows in canonical schema).
3. `results/v3/S4/coverage_matrix.csv` (Property coverage matrix v2 with explicit `FLAG_BELOW_10` status checks).
4. `splits/v3/real_fuels_v3_splits.json` (5-fold nested splits verifying 0% train/test leakage across all held-out real fuels).
5. `data/unified/compounds.csv` (First-class registration of ethanol, methanol, MTBE, and ETBE with Antoine and LFL data).

## Property Coverage Matrix Highlights
All primary core properties across Jet, Diesel, and Gasoline achieve robust coverage ($\ge 10$ fuels).
