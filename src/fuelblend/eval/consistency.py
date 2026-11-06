"""
src/fuelblend/eval/consistency.py

Multi-Property Physical Consistency Checker for Hydrocarbon Fuel Blends.

Validates that candidate blend property predictions satisfy fundamental thermodynamic
and fluid-mechanical coupling rules:
1. Flash Point vs Distillation (Mackison / API correlation; T_flash < T10).
2. Kinematic Viscosity vs Freezing Point (Physical state consistency: no fluid flow below WAT/freeze point).
3. Kinematic Viscosity vs Density (Walther-ASTM hydrocarbon envelope).
4. Derived Cetane Number (DCN) vs Aromatic Content (Combustion kinetic suppression).
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional


@dataclass
class ConsistencyResult:
    is_consistent: bool
    violations: List[str] = field(default_factory=list)
    warnings: List[str] = field(default_factory=list)
    details: Dict[str, Any] = field(default_factory=dict)


class PhysicalConsistencyChecker:
    """Enforces multi-property thermodynamic and physical consistency rules."""

    @staticmethod
    def check_flash_vs_distillation(
        flash_point_C: Optional[float],
        T10_C: Optional[float],
        T0_IBP_C: Optional[float] = None,
    ) -> List[str]:
        """Validate flash point against distillation front-end cuts."""
        violations = []
        if flash_point_C is None:
            return violations

        # Strict physical limit: Fuel cannot flash higher than its 10% recovery temperature
        if T10_C is not None:
            if flash_point_C >= T10_C:
                violations.append(
                    f"Thermodynamic violation: flash point ({flash_point_C:.1f} C) >= T10 distillation ({T10_C:.1f} C)."
                )
            elif flash_point_C >= (T10_C - 5.0):
                violations.append(
                    f"Physical warning: flash point ({flash_point_C:.1f} C) abnormally close to T10 ({T10_C:.1f} C)."
                )

            # Mackison / API correlation check: T_flash_expected ~ 0.683 * T10 - 66.8
            expected_flash = 0.683 * T10_C - 66.8
            residual = abs(flash_point_C - expected_flash)
            if residual > 35.0:
                violations.append(
                    f"Mackison correlation anomaly: flash point ({flash_point_C:.1f} C) deviates by {residual:.1f} C "
                    f"from API expectation ({expected_flash:.1f} C based on T10={T10_C:.1f} C)."
                )

        if T0_IBP_C is not None:
            if flash_point_C >= (T0_IBP_C + 15.0):
                violations.append(
                    f"Thermodynamic violation: flash point ({flash_point_C:.1f} C) exceeds IBP ({T0_IBP_C:.1f} C) + 15 C."
                )

        return violations

    @staticmethod
    def check_freezing_vs_viscosity(
        freezing_point_C: Optional[float],
        viscosity_minus20C_cSt: Optional[float],
    ) -> List[str]:
        """Validate that fluid flow can physically occur at -20 C."""
        violations = []
        if freezing_point_C is None or viscosity_minus20C_cSt is None:
            return violations

        # If fuel freezes above -20 C, it is solid/wax at -20 C
        if freezing_point_C > -20.0:
            violations.append(
                f"Physical state contradiction: fuel freezes at {freezing_point_C:.1f} C (> -20 C), "
                f"so liquid viscosity at -20 C ({viscosity_minus20C_cSt:.2f} cSt) cannot exist (fuel is solid/wax)."
            )
        elif freezing_point_C > -25.0 and viscosity_minus20C_cSt < 4.0:
            violations.append(
                f"Rheological anomaly: fuel freezes at {freezing_point_C:.1f} C (near -20 C) "
                f"yet claims ultra-low viscosity ({viscosity_minus20C_cSt:.2f} cSt) near cloud point."
            )

        return violations

    @staticmethod
    def check_viscosity_vs_density(
        viscosity_minus20C_cSt: Optional[float],
        density_15C_kg_m3: Optional[float],
    ) -> List[str]:
        """Validate hydrocarbon density-viscosity envelope."""
        violations = []
        if viscosity_minus20C_cSt is None or density_15C_kg_m3 is None:
            return violations

        # Hydrocarbon jet envelope: density 720 - 860 kg/m3
        if density_15C_kg_m3 < 730.0 and viscosity_minus20C_cSt > 10.0:
            violations.append(
                f"Hydrocarbon envelope anomaly: ultra-low density ({density_15C_kg_m3:.1f} kg/m3) "
                f"coupled with high cold viscosity ({viscosity_minus20C_cSt:.2f} cSt)."
            )
        elif density_15C_kg_m3 > 860.0 and viscosity_minus20C_cSt < 1.5:
            violations.append(
                f"Hydrocarbon envelope anomaly: high aromatic density ({density_15C_kg_m3:.1f} kg/m3) "
                f"coupled with unnaturally low viscosity ({viscosity_minus20C_cSt:.2f} cSt)."
            )

        return violations

    @staticmethod
    def check_cetane_vs_aromatics(
        dcn: Optional[float],
        aromatics_vol_pct: Optional[float],
    ) -> List[str]:
        """Validate combustion ignition delay against aromatic fraction."""
        violations = []
        if dcn is None or aromatics_vol_pct is None:
            return violations

        if aromatics_vol_pct > 30.0 and dcn > 50.0:
            violations.append(
                f"Combustion kinetics violation: aromatics ({aromatics_vol_pct:.1f} vol%) "
                f"cannot yield high cetane number (DCN={dcn:.1f}; aromatics suppress ignition)."
            )
        elif aromatics_vol_pct < 5.0 and dcn < 25.0:
            violations.append(
                f"Combustion kinetics anomaly: paraffinic blend (< 5% aromatics) "
                f"exhibits abnormally low cetane number (DCN={dcn:.1f})."
            )

        return violations

    def check_candidate_blend(self, properties: Dict[str, Any]) -> ConsistencyResult:
        """Run all physical consistency checks on a candidate blend dictionary."""
        all_violations = []
        all_warnings = []

        flash_C = properties.get("flash_point_C")
        t10_C = properties.get("distillation_T10_C") or properties.get("T10_C")
        ibp_C = properties.get("distillation_T0_C") or properties.get("T0_C")
        freeze_C = properties.get("freezing_point_C")
        visc20 = properties.get("viscosity_minus20C_cSt") or properties.get("viscosity_kin_minus20C")
        dens15 = properties.get("density_15C_kg_m3") or properties.get("density_kg_m3")
        dcn = properties.get("dcn") or properties.get("derived_cetane_number")
        arom = properties.get("aromatics_vol_pct") or properties.get("aromatics_pct")

        v_flash = self.check_flash_vs_distillation(flash_C, t10_C, ibp_C)
        v_freeze = self.check_freezing_vs_viscosity(freeze_C, visc20)
        v_visc = self.check_viscosity_vs_density(visc20, dens15)
        v_cetane = self.check_cetane_vs_aromatics(dcn, arom)

        for v in v_flash + v_freeze + v_visc + v_cetane:
            if "warning" in v.lower():
                all_warnings.append(v)
            else:
                all_violations.append(v)

        return ConsistencyResult(
            is_consistent=(len(all_violations) == 0),
            violations=all_violations,
            warnings=all_warnings,
            details={
                "flash_checked": flash_C is not None and t10_C is not None,
                "freezing_checked": freeze_C is not None and visc20 is not None,
                "viscosity_density_checked": visc20 is not None and dens15 is not None,
                "cetane_aromatics_checked": dcn is not None and arom is not None,
            },
        )
