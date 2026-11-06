"""Schema validators and definitions for measurements and compositions (T1.1, T1.2)."""
from __future__ import annotations

import argparse
import sys
from typing import Any, Dict, List
import numpy as np
import pandas as pd

# Canonical properties supported across all fuel families and consistency checks
CANONICAL_PROPERTIES = {
    # Core rheology & thermophysical
    "viscosity_kin_cSt",
    "viscosity_dyn_mPas",
    "density_kg_m3",
    "vapor_pressure_kPa",
    "surface_tension_mN_m",
    # Thermal & safety
    "flash_point_C",
    "autoignition_temp_C",
    # Low-temperature / cold flow
    "freezing_point_C",
    "cloud_point_C",
    "pour_point_C",
    "wax_appearance_temp_C",
    # Combustion & performance
    "heat_of_combustion_net_MJ_kg",
    "heat_of_combustion_gross_MJ_kg",
    "dcn",
    "ron",
    "mon",
    "octane_sensitivity",
    "smoke_point_mm",
    "yield_sooting_index",
    # Distillation checks (ASTM D86 / D2887)
    "distillation_TIBP_C",
    "distillation_T10_C",
    "distillation_T20_C",
    "distillation_T50_C",
    "distillation_T90_C",
    "distillation_TFBP_C",
}

CANONICAL_UNITS = {
    "cSt",
    "mPa*s",
    "mPa.s",
    "kg/m3",
    "kPa",
    "mN/m",
    "C",
    "K",
    "MJ/kg",
    "dimensionless",
    "mm",
}

PROPERTY_DEFAULT_UNITS = {
    "viscosity_kin_cSt": "cSt",
    "viscosity_dyn_mPas": "mPa*s",
    "density_kg_m3": "kg/m3",
    "vapor_pressure_kPa": "kPa",
    "surface_tension_mN_m": "mN/m",
    "flash_point_C": "C",
    "autoignition_temp_C": "C",
    "freezing_point_C": "C",
    "cloud_point_C": "C",
    "pour_point_C": "C",
    "wax_appearance_temp_C": "C",
    "heat_of_combustion_net_MJ_kg": "MJ/kg",
    "heat_of_combustion_gross_MJ_kg": "MJ/kg",
    "dcn": "dimensionless",
    "ron": "dimensionless",
    "mon": "dimensionless",
    "octane_sensitivity": "dimensionless",
    "smoke_point_mm": "mm",
    "yield_sooting_index": "dimensionless",
    "distillation_TIBP_C": "C",
    "distillation_T10_C": "C",
    "distillation_T20_C": "C",
    "distillation_T50_C": "C",
    "distillation_T90_C": "C",
    "distillation_TFBP_C": "C",
}

REQUIRED_MEASUREMENT_PROVENANCE = [
    "measurement_id",
    "source_id",
    "source_citation",
    "system_id",
    "property",
    "value",
    "unit",
    "method",
    "composition_basis",
    "quality_flag",
]

ALLOWED_COMPOSITION_BASES = {"mass", "volume", "mole", "none"}
ALLOWED_QUALITY_FLAGS = {"ok", "suspect", "conflict", "outlier"}


class ValidationError(Exception):
    """Raised when data fails schema validation rules."""
    pass


def validate_measurement_df(df: pd.DataFrame, strict: bool = True) -> Dict[str, Any]:
    """Validate long-format measurement DataFrame against schema rules."""
    errors = []

    # 1. Check required columns exist
    missing_cols = [c for c in REQUIRED_MEASUREMENT_PROVENANCE if c not in df.columns]
    if missing_cols:
        msg = f"Missing required columns in measurement table: {missing_cols}"
        if strict:
            raise ValidationError(msg)
        return {"valid": False, "errors": [msg]}

    # 2. Check for empty/null required provenance fields
    for col in ["source_id", "source_citation", "system_id", "property", "method", "composition_basis", "quality_flag"]:
        null_mask = df[col].isna() | (df[col].astype(str).str.strip() == "")
        if null_mask.any():
            bad_indices = df[null_mask].index.tolist()[:10]
            errors.append(f"Column '{col}' has empty values at rows: {bad_indices}")

    # 3. Check properties and units
    bad_props = set(df["property"].dropna()) - CANONICAL_PROPERTIES
    if bad_props:
        errors.append(f"Unknown properties detected: {sorted(list(bad_props))}")

    bad_units = set(df["unit"].dropna()) - CANONICAL_UNITS
    if bad_units:
        errors.append(f"Unknown units detected: {sorted(list(bad_units))}")

    # 4. Check values are numeric
    non_numeric = pd.to_numeric(df["value"], errors="coerce").isna() & df["value"].notna()
    if non_numeric.any():
        errors.append(f"Non-numeric values in 'value' column at rows: {df[non_numeric].index.tolist()[:10]}")

    # 5. Check quality flags and composition basis
    bad_flags = set(df["quality_flag"].dropna()) - ALLOWED_QUALITY_FLAGS
    if bad_flags:
        errors.append(f"Unknown quality flags: {sorted(list(bad_flags))}")

    bad_bases = set(df["composition_basis"].dropna()) - ALLOWED_COMPOSITION_BASES
    if bad_bases:
        errors.append(f"Unknown composition bases: {sorted(list(bad_bases))}")

    if errors and strict:
        raise ValidationError("Measurement validation failed:\n" + "\n".join(errors))

    return {"valid": len(errors) == 0, "errors": errors}


def validate_composition_df(df: pd.DataFrame, strict: bool = True, tolerance: float = 1e-3) -> Dict[str, Any]:
    """Validate mixture composition DataFrame."""
    errors = []
    required_cols = ["system_id", "inchikey", "fraction", "basis"]
    missing_cols = [c for c in required_cols if c not in df.columns]
    if missing_cols:
        msg = f"Missing required columns in composition table: {missing_cols}"
        if strict:
            raise ValidationError(msg)
        return {"valid": False, "errors": [msg]}

    if "unexplained_fraction" not in df.columns:
        df = df.copy()
        df["unexplained_fraction"] = 0.0

    # Group by system_id and basis to verify fraction summation
    for (sys_id, basis), group in df.groupby(["system_id", "basis"]):
        frac_sum = group["fraction"].sum()
        unexplained = group["unexplained_fraction"].iloc[0] if "unexplained_fraction" in group.columns else 0.0
        if pd.isna(unexplained):
            unexplained = 0.0

        total = frac_sum + unexplained
        if abs(total - 1.0) > tolerance and abs(frac_sum - 1.0) > tolerance:
            errors.append(
                f"System {sys_id} (basis={basis}) fraction sum is {frac_sum:.4f} "
                f"(unexplained={unexplained:.4f}, total={total:.4f}), expected 1.0 +- {tolerance}"
            )

    if errors and strict:
        raise ValidationError("Composition validation failed:\n" + "\n".join(errors[:20]))

    return {"valid": len(errors) == 0, "errors": errors}


def main():
    parser = argparse.ArgumentParser(description="Validate measurement and composition schemas.")
    parser.add_argument("--measurements", type=str, help="Path to measurement CSV")
    parser.add_argument("--compositions", type=str, help="Path to composition CSV")
    args = parser.parse_args()

    success = True
    if args.measurements:
        print(f"Validating measurements: {args.measurements}")
        df_m = pd.read_csv(args.measurements)
        res_m = validate_measurement_df(df_m, strict=False)
        if res_m["valid"]:
            print("OK: Measurements valid.")
        else:
            print("FAIL: Measurement errors:")
            for e in res_m["errors"]:
                print(f"  - {e}")
            success = False

    if args.compositions:
        print(f"Validating compositions: {args.compositions}")
        df_c = pd.read_csv(args.compositions)
        res_c = validate_composition_df(df_c, strict=False)
        if res_c["valid"]:
            print("OK: Compositions valid.")
        else:
            print("FAIL: Composition errors:")
            for e in res_c["errors"]:
                print(f"  - {e}")
            success = False

    sys.exit(0 if success else 1)


if __name__ == "__main__":
    main()
