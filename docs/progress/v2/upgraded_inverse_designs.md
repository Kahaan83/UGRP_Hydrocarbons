# Progress Report: Upgraded Inverse Blend Designs (Task 4)

**Date**: 2026-11-05  
**Deliverable**: `results/v2/improvement/upgraded_designs.csv`  

---

## 1. Executive Summary

Using the upgraded thermodynamic evaluators (anchored ambient vapor pressure, Coutinho multi-component wax depression, ASTM D341 Walther viscosity curvature, and recalibrated cetane/octane values), we evaluated and confirmed **100% constraint compliance** for all four priority designs under strict 90% conformal safety margins:

1. **Jet A Minimal 3-Component Blend**:
   - Formulation: 65% isododecane + 20% decalin + 15% tetralin
   - Freezing Point: -51.3 °C (passes <= -40.0 °C limit with 11.3 °C safety margin)
   - Kinematic Viscosity (-20 °C): 4.88 cSt (passes <= 8.0 cSt limit with 3.12 cSt margin)
   - Flash Point: 44.8 °C (passes >= 38.0 °C limit with 6.8 °C margin)
   - Density at 15 °C: 794.7 kg/m³ (passes [775, 840] envelope)
   - Net Heat of Combustion: 43.43 MJ/kg (passes >= 42.8 MJ/kg)
   - Aromatics: 15.0 wt% / 13.9 vol% (passes <= 25.0 wt% limit)
   - **All Constraints Passed with Calibrated Margins**: **TRUE**

2. **Jet A Realistic 4-Component Blend**:
   - Formulation: 58% isododecane + 12% n-dodecane + 15% decalin + 15% tetralin
   - Freezing Point: -47.1 °C (passes <= -40.0 °C with 7.1 °C margin)
   - Viscosity (-20 °C): 4.86 cSt (passes <= 8.0 cSt)
   - Flash Point: 45.1 °C (passes >= 38.0 °C)
   - Density: 794.6 kg/m³
   - **All Constraints Passed**: **TRUE**

3. **Diesel No. 2 Realistic 5-Component Blend**:
   - Formulation: 21% n-hexadecane + 7% n-octadecane + 27% HMN + 25% bicyclohexyl + 20% tetralin
   - DCN: 48.5 (passes >= 40.0 limit with 8.5 margin)
   - Cloud Point: -3.4 °C (passes <= 0.0 °C)
   - Viscosity (40 °C): 2.89 cSt (passes [1.9, 4.1] cSt envelope)
   - Density: 834.7 kg/m³ (passes [820, 845] kg/m³)
   - Flash Point: 83.2 °C (passes >= 52.0 °C)
   - **All Constraints Passed**: **TRUE**

4. **Gasoline Realistic 5-Component Blend**:
   - Formulation: 28% isooctane + 36% toluene + 20% isopentane + 5.7% n-butane + 10.3% diisobutylene
   - RON / MON / Sensitivity: 104.6 / 96.7 / 7.9 (passes >= 91 / >= 82 / [7, 14])
   - DVPE Volatility: 57.1 kPa (passes [45, 70] kPa summer envelope)
   - Density: 725.8 kg/m³ (passes [720, 775] kg/m³)
   - **All Constraints Passed**: **TRUE**
