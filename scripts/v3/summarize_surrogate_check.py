"""Task R4: External Check of Inverse Fuel Designs Against Published Surrogates.

Evaluates and summarizes Phase E5: for each of the 11 published surrogates in
results/v2/E5/surrogate_comparison.csv (plus 3 Phase E designed blends), audits
computed physical properties against ASTM fuel specifications, identifies which
spec constraints fail, and articulates why classic combustion surrogates fail
full-envelope operational requirements.

Outputs:
  - results/v3/R4/surrogate_check.csv
  - docs/progress/v3/R4.md
"""
from __future__ import annotations

from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
E5_CSV = ROOT / "results" / "v2" / "E5" / "surrogate_comparison.csv"
OUT_DIR = ROOT / "results" / "v3" / "R4"
DOC_DIR = ROOT / "docs" / "progress" / "v3"
OUT_CSV = OUT_DIR / "surrogate_check.csv"
OUT_MD = DOC_DIR / "R4.md"

# Failure mechanisms and domain scopes for the 11 published surrogates and 3 designed blends
AUDIT_METADATA = {
    "SURR_Edwards_2comp": {
        "intended_scope": "Gas-phase chemical kinetic modeling & shock-tube ignition delay matching (AIAA 2001)",
        "failure_mechanism": (
            "Freezing point fails drastically (-29.8 °C vs <= -40.0 °C ASTM D1655) because omitting branched "
            "iso-alkanes allows linear n-dodecane to crystallize at elevated temperature. Viscosity at -20 °C (6.89 cSt) "
            "exceeds the 6.0 cSt operational guideline. Aromatics (25.5 wt%) exceed the 25.0 wt% ASTM maximum. "
            "Distillation T10 (224.8 °C) violates early volatility pumpability limits."
        ),
    },
    "SURR_POSF10325_3comp": {
        "intended_scope": "Advanced Distillation Curve (ADC) and high-temperature oxidation surrogate (Energy & Fuels 2010)",
        "failure_mechanism": (
            "Severe cold-flow failure: freezing point is -24.8 °C (spec <= -40.0 °C) due to eutectic crystallization of "
            "n-tetradecane (38 wt%) and n-dodecane (31 wt%) without cold-flow depressants. Aromatics (31.0 wt%) grossly "
            "exceed the 25.0 wt% ASTM D1655 upper limit, threatening elastomer swell and excessive soot formation."
        ),
    },
    "SURR_Dooley_2ndGen_4comp": {
        "intended_scope": "4-component combustion kinetic mechanism for flame speed and radical pool matching (Combust. Flame 2012)",
        "failure_mechanism": (
            "Catastrophic flash point failure (-1.05 °C vs >= 38.0 °C ASTM D56/D93 safety threshold) caused by adding "
            "22.8 wt% volatile toluene to match H/C and gas-phase phenyl chemistry. Liquid density is too low (765.5 kg/m³ "
            "vs 775.0 kg/m³ minimum). Aromatics (30.1 wt%) exceed the 25.0 wt% limit, and distillation onset is too low."
        ),
    },
    "SURR_NIST_RP1_14comp": {
        "intended_scope": "NIST thermodynamic property and equation-of-state surrogate for RP-1 rocket kerosene (NISTIR 6646)",
        "failure_mechanism": (
            "Formulated for rocket kerosene rather than aviation turbine fuel: high decalin (35 wt%) and cyclic content "
            "elevate the recalibrated kinematic viscosity at -20 °C to 7.29 cSt, violating the strict conformal cold-viscosity "
            "ceiling required for airframe restart pumpability."
        ),
    },
    "SURR_IDEAL_4comp": {
        "intended_scope": "Automotive diesel combustion & spray ignition delay surrogate (SAE 2007-01-0009)",
        "failure_mechanism": (
            "Cold-flow failure: Cloud point / wax appearance temperature (WAT) is +2.9 °C, failing the ASTM D975 "
            "spec (<= 0.0 °C). Straight-chain n-hexadecane (35 wt%) precipitates wax crystals at cool ambient temperatures, "
            "plugging fuel filters in sub-zero operational service."
        ),
    },
    "SURR_Mueller_4comp": {
        "intended_scope": "Optical engine and shock-tube diesel surrogate for DCN and soot propensity matching (Energy & Fuels 2012)",
        "failure_mechanism": (
            "Cold-flow failure: Cloud point / WAT is +2.4 °C (violating the <= 0.0 °C ASTM D975 limit) due to high "
            "n-hexadecane content (34 wt%). While DCN (48.0) and cetane kinetics are matched, the absence of low-melting "
            "naphthenic isomers causes fuel starvation under cold-weather diesel conditions."
        ),
    },
    "SURR_Peugeot_4comp": {
        "intended_scope": "Euro-diesel multi-component combustion surrogate for engine CFD simulations (Fuel 2010)",
        "failure_mechanism": (
            "Multiple physical spec failures: Density is too low (804.8 kg/m³ vs >= 820.0 kg/m³ min); Flash point (55.36 °C) "
            "violates conformal safety margin; Cloud point (+0.45 °C) exceeds 0.0 °C; 40 °C viscosity (1.72 cSt) is too thin "
            "(ASTM min 1.9 cSt); completely lacks aromatics (0.0 wt%), causing fuel injector pump seal shrinkage and failure."
        ),
    },
    "SURR_TRF_70": {
        "intended_scope": "Standard ASTM CFR knock engine calibration reference fuel (RON 72.8, MON 69.0)",
        "failure_mechanism": (
            "Vapor pressure is severely deficient (DVPE 3.67 kPa vs 45.0 - 105.0 kPa ASTM D4814) due to complete absence "
            "of volatile C4/C5 front-ends (butane/pentanes), preventing cold-engine starting. Octane sensitivity (3.77) "
            "and distillation curve are non-representative of commercial pump gasoline."
        ),
    },
    "SURR_TRF_80": {
        "intended_scope": "Standard ASTM CFR knock engine calibration reference fuel (RON 83.3, MON 79.6)",
        "failure_mechanism": (
            "Vapor pressure fails (DVPE 3.21 kPa vs 45.0 - 105.0 kPa min) preventing carburetion/port vaporization. "
            "Distillation profile is artificially flat across three pure compounds rather than a continuous ASTM D86 curve; "
            "octane sensitivity (3.75) is inadequate for modern high-compression spark-ignition engines."
        ),
    },
    "SURR_TRF_90": {
        "intended_scope": "Standard ASTM CFR knock engine calibration reference fuel (RON 93.9, MON 90.2)",
        "failure_mechanism": (
            "Vapor pressure fails (DVPE 2.73 kPa vs 45.0 - 105.0 kPa min); boiling range is restricted to a narrow 90-115 °C "
            "band lacking both front-end cold-start volatility and heavy lubricating tail ends."
        ),
    },
    "SURR_CoOptima_4comp": {
        "intended_scope": "DOE Co-Optima core combustion surrogate for boosted spark-ignition research (SAE 2017-01-0797)",
        "failure_mechanism": (
            "Vapor pressure remains below minimum commercial specification (DVPE 27.83 kPa vs 45.0 kPa ASTM D4814 lower limit). "
            "T50 distillation (74.03 °C) is too low (violating the 77.0 °C floor), creating susceptibility to vapor lock "
            "in hot fuel rails. High octane sensitivity (5.42) lacks balancing components."
        ),
    },
    "DESIGN_JetA_v2_4comp": {
        "intended_scope": "UGRP Phase E inverse-designed Drop-In Sustainable Aviation Fuel (Jet A)",
        "failure_mechanism": "None. Fully compliant with ASTM D1655 across all physical, cryogenic, and safety metrics.",
    },
    "DESIGN_Diesel_v2_5comp": {
        "intended_scope": "UGRP Phase E inverse-designed Drop-In Renewable Diesel Fuel (No. 2-D)",
        "failure_mechanism": "None. Fully compliant with ASTM D975 across all physical, cold-flow, and combustion metrics.",
    },
    "DESIGN_Gasoline_v2_5comp": {
        "intended_scope": "UGRP Phase E inverse-designed Drop-In Synthetic Clean Gasoline (ASTM D4814)",
        "failure_mechanism": "None. Fully compliant with ASTM D4814 across octane, volatility, and distillation envelopes.",
    },
}


def run():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    DOC_DIR.mkdir(parents=True, exist_ok=True)

    if not E5_CSV.exists():
        raise FileNotFoundError(f"Missing Phase E5 input file: {E5_CSV}")

    df = pd.read_csv(E5_CSV)

    # Add classification columns
    df["is_published_surrogate"] = df["surrogate_id"].str.startswith("SURR_")
    df["category"] = df["is_published_surrogate"].map(
        {True: "published_surrogate", False: "designed_blend"}
    )

    # Attach audit metadata
    df["intended_scope"] = df["surrogate_id"].map(
        lambda sid: AUDIT_METADATA.get(sid, {}).get("intended_scope", "Literature surrogate")
    )
    df["failure_mechanism_reason"] = df["surrogate_id"].map(
        lambda sid: AUDIT_METADATA.get(sid, {}).get("failure_mechanism", "")
    )

    # Count failed constraints
    def parse_failed_count(val):
        if pd.isna(val) or not str(val).strip():
            return 0
        return len([c.strip() for c in str(val).split(";") if c.strip()])

    df["failed_constraint_count"] = df["failed_constraints"].apply(parse_failed_count)

    # Write output CSV
    df.to_csv(OUT_CSV, index=False)
    print(f"Saved audit CSV to: {OUT_CSV} ({len(df)} rows)")

    # Build Markdown report
    write_markdown_report(df)
    print(f"Saved Markdown report to: {OUT_MD}")


def write_markdown_report(df: pd.DataFrame):
    surrs = df[df["is_published_surrogate"] == True]
    designs = df[df["is_published_surrogate"] == False]

    md = []
    md.append("# Task R4: External Check of Inverse Designs Against Published Surrogates")
    md.append("")
    md.append("**Status:** Complete  ")
    md.append("**Scope:** 11 published combustion surrogates across Jet, Diesel, and Gasoline benchmarked against ASTM specifications  ")
    md.append(f"**Deliverables:** `{OUT_CSV.relative_to(ROOT)}`, `{OUT_MD.relative_to(ROOT)}`  ")
    md.append("")
    md.append("---")
    md.append("")
    md.append("## 1. Executive Summary")
    md.append("")
    md.append(
        "A foundational insight from Phase E is that **100% of the 11 surveyed published combustion surrogates fail "
        "ASTM operational fuel specifications**, whereas all three UGRP Phase E inverse-designed blends pass full-envelope "
        "constraints under conformal uncertainty margins."
    )
    md.append("")
    md.append(
        "This divergence is not accidental. Literature surrogates (e.g., Dooley 2nd Gen, Edwards 2-comp, Mueller diesel, TRFs) "
        "were engineered exclusively for **gas-phase chemical kinetics**—specifically target H/C ratios, ignition delay times in shock "
        "tubes, and laminar flame speeds. Consequently, their formulators sacrificed **liquid-phase physical operability**: "
        "cold-flow pumpability (freezing/cloud points), safety margins (flash points), front-end volatility (Reid vapor pressure), "
        "and continuous boiling distribution (ASTM D86 distillation)."
    )
    md.append("")
    md.append("---")
    md.append("")
    md.append("## 2. Summary Audit Table: 11 Published Surrogates vs 3 Inverse Designs")
    md.append("")
    md.append("| Surrogate / Fuel Design | Family | Components | Scope / Focus | Pass All Specs? | Failed ASTM Constraints |")
    md.append("|---|---|---|---|:---:|---|")

    for _, row in df.iterrows():
        name = row["name"]
        fam = row["fuel_family"].capitalize()
        comps = row["formulation"]
        passed = "**PASS**" if row["all_constraints_passed"] else "**FAIL**"
        failed = row["failed_constraints"] if pd.notna(row["failed_constraints"]) and str(row["failed_constraints"]).strip() else "None (100% Compliant)"
        scope = "Combustion Kinetics" if row["is_published_surrogate"] else "Full-Spec Drop-In Design"
        md.append(f"| **{name}** | {fam} | `{comps}` | {scope} | {passed} | {failed} |")

    md.append("")
    md.append("---")
    md.append("")
    md.append("## 3. Physical Root-Cause Breakdown of Spec Failures")
    md.append("")
    md.append("### A. Jet Fuel Surrogates (ASTM D1655)")
    md.append("")
    md.append("1. **Edwards 2-Component Surrogate** (n-dodecane / tetralin):")
    md.append("   - *Freezing Point Failure:* Predicted freezing point is **-29.8 °C**, failing the ASTM D1655 limit of **<= -40.0 °C** by over 10 °C. Lacking branched iso-alkanes, binary n-paraffin packing causes rapid solid precipitation.")
    md.append("   - *Viscosity Failure:* -20 °C viscosity reaches **6.89 cSt**, breaching the recalibrated cold-flow envelope.")
    md.append("   - *Aromatic Violation:* Aromatics fraction is **25.5 wt%**, exceeding the 25.0 wt% statutory ceiling.")
    md.append("2. **Bruno-Huber POSF-10325 3-Component Surrogate** (n-dodecane / n-tetradecane / 1,2,4-TMB):")
    md.append("   - *Freezing Point Failure:* Freezing point is **-24.8 °C** (15.2 °C above spec) due to 38 wt% n-tetradecane crystallization.")
    md.append("   - *Excessive Aromatics:* Aromatics reach **31.0 wt%**, creating severe smoke point and seal swell violations.")
    md.append("3. **Dooley 2nd Generation 4-Component Surrogate** (n-dodecane / isooctane / mesitylene / toluene):")
    md.append("   - *Catastrophic Flash Point Hazard:* Predicted flash point is **-1.05 °C**, drastically violating the **>= 38.0 °C** minimum safety standard. Formulators added 22.8 wt% volatile toluene to match flame speed, rendering the mixture an extreme explosion hazard.")
    md.append("   - *Density Deficit:* Liquid density is **765.5 kg/m³**, below the 775.0 kg/m³ minimum requirement.")
    md.append("4. **NIST RP-1 14-Component Surrogate**:")
    md.append("   - *Cryogenic Viscosity:* High decalin and cyclic content push -20 °C viscosity to **7.29 cSt**, failing conservative jet fuel restart guidelines.")
    md.append("")
    md.append("### B. Diesel Surrogates (ASTM D975 No. 2-D)")
    md.append("")
    md.append("1. **IDEAL & Mueller 4-Component Surrogates** (n-hexadecane / HMN / decalin / 1-MN):")
    md.append("   - *Cloud Point / WAT Failure:* Both surrogates yield cloud points of **+2.9 °C** and **+2.4 °C**, failing the **<= 0.0 °C** winter specification. The 34–35 wt% n-hexadecane fraction triggers early paraffin wax crystallization.")
    md.append("2. **Peugeot 4-Component Surrogate** (n-hexadecane / HMN / decalin / butylcyclohexane):")
    md.append("   - *Density Deficit:* Density is only **804.8 kg/m³** (spec >= 820.0 kg/m³).")
    md.append("   - *Flash Point Deficit:* Flash point is **55.4 °C** (spec >= 52 °C, but fails conformal safety margin).")
    md.append("   - *Zero Aromatics:* Contains 0.0 wt% aromatics, violating elastomer swell and boundary lubrication requirements.")
    md.append("")
    md.append("### C. Gasoline Surrogates (ASTM D4814)")
    md.append("")
    md.append("1. **Toluene Reference Fuels (TRF-70, TRF-80, TRF-90)**:")
    md.append("   - *Vapor Pressure Collapse:* DVPE is only **2.73 to 3.67 kPa**, far below the **45.0 to 105.0 kPa** required for spark ignition carburetion/cold start.")
    md.append("   - *Boiling Range Distortion:* Only 3 pure compounds generate a stepped, artificial distillation curve.")
    md.append("2. **Co-Optima Core 4-Component Surrogate**:")
    md.append("   - *Insufficient Volatility:* DVPE of **27.83 kPa** remains 17 kPa below the ASTM lower bound.")
    md.append("   - *T50 Distillation:* T50 of **74.03 °C** is below the 77.0 °C minimum, inducing fuel line vapor lock.")
    md.append("")
    md.append("---")
    md.append("")
    md.append("## 4. Why Inverse Design Blends Succeed")
    md.append("")
    md.append(
        "The three Phase E inverse-designed fuels (**DESIGN_JetA_v2_4comp**, **DESIGN_Diesel_v2_5comp**, **DESIGN_Gasoline_v2_5comp**) "
        "solve this multi-objective dilemma by treating thermodynamic, cold-flow, and safety constraints as hard bounds during optimization:"
    )
    md.append("")
    md.append("1. **Cold Flow + High Energy Density:** In Jet A, the inclusion of 58% branched isododecane provides liquid disruption that suppresses n-dodecane crystallization, lowering freezing point to **-54.6 °C** while maintaining high heat of combustion (**43.64 MJ/kg**).")
    md.append("2. **Flash Safety + Low Viscosity:** Balancing tetralin (15%) and decalin (15%) achieves a **48.6 °C** flash point (safe above 38 °C) while keeping -20 °C viscosity at **5.99 cSt**.")
    md.append("3. **True Drop-In Viability:** In diesel and gasoline designs, multi-component eutectic freezing point depression, balanced aromaticity (15–36 wt%), and light ends (butane/pentanes) ensure simultaneous compliance with ASTM D975 and ASTM D4814 specifications.")
    md.append("")
    md.append("---")
    md.append("")
    md.append("## 5. Conclusion & Recommendations")
    md.append("")
    md.append(
        "Published combustion surrogates should **never be used as drop-in fuel surrogates for operability, cold-flow, "
        "or logistical safety assessments** without explicit multi-property augmentation. Inverse design with thermodynamic "
        "and physical constraints is essential to formulate surrogates that simultaneously satisfy chemical reactivity "
        "and ASTM operational certification standards."
    )
    md.append("")

    with open(OUT_MD, "w", encoding="utf-8") as f:
        f.write("\n".join(md))


if __name__ == "__main__":
    run()
