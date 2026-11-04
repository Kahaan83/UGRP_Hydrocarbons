# Progress Report: NIST ThermoML Outlier Filtering & Tidy Dataset

**Date**: 2026-11-05  
**Source Archive**: `data/raw/thermoml/ThermoML.v2020-09-30.tgz` (189 MB)  
**Deliverable**: `data/processed/v2/thermoml_clean_hc.csv`  

---

## 1. Executive Summary

The raw NIST ThermoML archive compiles all published IUPAC data tables regardless of internal TRC experimental flags or measurement conditions. A strict physical filter was applied to remove:
- 60 unphysical non-positive viscosities ($\le 0$ Pa·s) resulting from negative excess viscosity placeholders.
- 130 gas-phase / near-critical densities ($< 550$ kg/m³).
- 1,426 high-pressure ($P > 500$ kPa) non-ambient measurements.
- 196 cryogenic ($T < 170$ K) and 57 thermal-cracking ($T > 550$ K) points.

Total rows evaluated: **46,747**  
Outliers dropped: **2172 (4.65%)**  
Clean liquid-phase rows retained: **44575 (95.35%)**

---

## 2. Clean Observations Breakdown

| Property Name | Clean Observations | Physical Liquid Bounds Applied |
|---|---|---|
| Mass density, kg/m³ | 19,404 | $550 \le ho \le 1100$ kg/m³, $P \le 500$ kPa |
| Vapor pressure, kPa | 13,310 | $P > 0$ kPa, liquid phase |
| Viscosity, Pa·s | 11,755 | $10^{-5} \le \eta \le 2.0$ Pa·s |
| Kinematic viscosity, m²/s | 106 | $
u > 0$ m²/s |
