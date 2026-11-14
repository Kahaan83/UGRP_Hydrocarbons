"""Generate configs/uncertainty_v5.yaml: results/v5/S1/uncertainty_v5.yaml with ONLY the flash lab ideal applied_margin replaced
by 2.206932808388615 (user decision Q2, 2026-11-13; row designated,old_plus_new_minus_mesitylene of delta_lab_flash_li2014/lab_margin.csv),
plus (2026-11-13, user decision "Yes, add it") the new entry smoke_point_mm.lab_mixtures.bulk_regression = 7.8009965628100595 mm
(row smoke_mixtures,bulk_designated of delta_newdata_freeze_smoke/summary.csv).  This script supersedes its first version (flash only)."""
from pathlib import Path
import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parents[3]
SRC = ROOT / "results/v5/S1/uncertainty_v5.yaml"
OUT = ROOT / "configs/uncertainty_v5.yaml"
CSV = ROOT / "results/v5/BASE/delta_lab_flash_li2014/lab_margin.csv"
KEY = "flash_point_C.lab_mixtures.ideal"
SMOKE_CSV = ROOT / "results/v5/BASE/delta_newdata_freeze_smoke/summary.csv"
SMOKE_KEY = "smoke_point_mm.lab_mixtures.bulk_regression"
SMOKE_MARGIN = 7.8009965628100595
FREEZE_DIR = ROOT / "results/v5/BASE/delta_freeze_margin"
FREEZE_KEY = "freezing_point_C.lab_liquidus.svl_verified"
FLASH_MARGIN = 2.206932808388615
SC_DIR = ROOT / "results/v5/BASE/delta_smoke_cert"
SC_KEY = "smoke_point_mm.ajftd_v5_in_domain.bulk_regression"


def add_smoke_cert(doc):
    """v5 delta_smoke_cert (user decision 2026-11-13, relayed by the parent session): real-fuel in-domain smoke margin, written only if
    the declared R-SC-2 gate passed.  It REPLACES smoke_point_mm.lab_mixtures.bulk_regression for in-domain jet recipes inside the class-mix range."""
    import json
    rd = json.loads((SC_DIR / "readings.json").read_text(encoding="utf-8"))
    if not rd["gate_passed"]:
        raise SystemExit("delta_smoke_cert R-SC-2 gate not passed: margin not written")
    m = float(rd["certifying_margin_final"])
    assert rd["direction"] == "over" and m > 0 and m < 10
    doc["margins"][SC_KEY] = {
        "value": m, "unit": "mm", "direction": "over", "alpha": 0.1, "n_groups": len(rd["eligible_groups"]), "group_unit": "source_group",
        "n_rows": int(rd["n_in_domain_fuels"]), "residual_file": "results/v5/BASE/delta_smoke_cert/margin_fit.csv", "applied_margin": m,
        "model": "smoke_point_from_bulk (aromatics only), margin = max(fuel-level, largest eligible per-group margin, production-predictor margin R-SC-2c)",
        "scope": "real AJFTD jet fuels, aromatics 8-23 vol %, uncensored smoke; applies only to recipes whose class mix lies inside the AJFTD GCxGC "
                 "range (src/fuelblend/design/smoke_margin_v5.py); otherwise smoke stays uncertified.",
        "residual_convention": "predicted - measured",
        "caveats": ["Smoke spec limit is secondary_sources_only: this margin does not make any design certified by itself.",
                    "Refit model is single-feature OLS; stored predictor differs by up to 1.87 mm (R-SC-2c); the larger margin is used."],
        "user_decision": "2026-11-13 'Yes, replace (Recommended)' to 'Should 4.11 mm REPLACE your earlier 7.8 mm lab-mixture margin for in-domain jet recipes "
                         "inside the real-fuel class-mix range?'; 'Yes, add it (Recommended)' to 'May WS1 add the key "
                         "smoke_point_mm.ajftd_v5_in_domain.bulk_regression = 4.11 to configs/uncertainty_v5.yaml?'",
        "source": "results/v5/BASE/delta_smoke_cert/readings.json certifying_margin_final"}


def add_freezing(doc):
    """v5 delta_freeze_margin: written ONLY if the declared rule passed (readings.json rule_write_margin).  Raises otherwise."""
    import json
    rd = json.loads((FREEZE_DIR / "readings.json").read_text(encoding="utf-8"))
    if not rd["rule_write_margin"]:
        raise SystemExit("delta_freeze_margin rule not passed: freezing margin not written")
    sm = pd.read_csv(FREEZE_DIR / "summary.csv").query("scope == 'S_F' and model == 'svl_verified'").iloc[0]
    m = float(sm["margin_alpha10"])
    assert m == rd["R-FM-3"]["margin_C"] and sm["margin_direction"] == "under" and int(sm["n_systems"]) >= 9
    assert FREEZE_KEY not in doc["margins"]
    doc["margins"][FREEZE_KEY] = {
        "value": m, "unit": "C", "direction": "under", "alpha": 0.1, "n_groups": int(sm["n_systems"]), "group_unit": "system_id",
        "n_rows": int(sm["n_scored"]), "residual_file": "results/v5/BASE/delta_freeze_margin/rows.csv", "applied_margin": m,
        "model": "ideal Schroeder-van Laar in mole fraction on verified Tm / dHfus (data/raw/sle_constants/tm_dhfus_verified.csv), designated freezing route",
        "scope": "ThermoML liquidus rows, all components C8-C16 hydrocarbons with verified constants, T <= 273.15 K (scope S_F): binaries of n-alkanes, "
                 "naphthalenes, decalins and bicyclohexyl, 9 systems. Valid for known compositions whose every component has a verified Tm and dHfus.",
        "residual_convention": "predicted - measured liquidus temperature (C); the margin is ADDED to the prediction (safe_bound = hi)",
        "caveats": ["9 systems is the minimum for alpha = 0.10; the margin is the largest group score (worst system).",
                    "Binaries only (2 components); 5 of the 9 systems were unscorable with the old Coutinho constants and were added after seeing that scoring.",
                    "A liquidus is not an ASTM D5972 freezing point; transfer to multi-component jet fuels is not demonstrated.",
                    "Real jet fuels (Vozka, reconstructed composition) are over-predicted by the same constants (safe side); the margin is not claimed for reconstructed real fuels."],
        "source": "results/v5/BASE/delta_freeze_margin/summary.csv row S_F,svl_verified margin_alpha10"}


def main():
    row = pd.read_csv(CSV).query("variant == 'designated' and scope == 'old_plus_new_minus_mesitylene'").iloc[0]
    new = float(row["margin"])
    assert new == 2.206932808388615 and int(row["n_groups"]) == 35
    doc = yaml.safe_load(SRC.read_text(encoding="utf-8"))
    e = doc["margins"][KEY]
    e["applied_margin_v4"] = e["applied_margin"]
    e["applied_margin"] = new
    e["applied_margin_source"] = "results/v5/BASE/delta_lab_flash_li2014/lab_margin.csv row designated,old_plus_new_minus_mesitylene"
    e["n_groups_applied"] = 35
    e["user_decision"] = "2026-11-13 'use 2.2 for Q2': mesitylene excluded (suspect Antoine data); margin set by the new ternary n-alkane system"
    sm = pd.read_csv(SMOKE_CSV).query("scope == 'smoke_mixtures' and model == 'bulk_designated'").iloc[0]
    assert float(sm["margin_alpha10"]) == SMOKE_MARGIN and sm["margin_direction"] == "over" and int(sm["n_systems"]) == 12 and int(sm["n_scored"]) == 56
    assert SMOKE_KEY not in doc["margins"]
    doc["margins"][SMOKE_KEY] = {
        "value": SMOKE_MARGIN, "unit": "mm", "direction": "over", "alpha": 0.1, "n_groups": 12, "group_unit": "system_id", "n_rows": 56,
        "residual_file": "results/v5/BASE/delta_newdata_freeze_smoke/rows_smoke.csv", "applied_margin": SMOKE_MARGIN,
        "model": "AJFTD bulk regression of smoke point on aromatics vol % (smoke_point_from_bulk), designated route",
        "scope": "smoke_mixtures: literature binaries and surrogates (12 systems, 56 scored points); not real jet fuels. Valid only inside the "
                 "regression domain (aromatics 8-23 vol %, result 19-25 mm); outside it the API and the designs stay uncertified.",
        "residual_convention": "predicted - measured",
        "caveats": ["Bound from binaries and surrogates; 37 of 56 scored points are outside the regression's 8-23 vol % aromatics domain.",
                    "The margin (7.8 mm) is wider than the 19-25 mm trusted range of the regression.",
                    "User chose it after seeing the numbers; probably certifies nothing (designs predict about 27.6 mm vs limit 25)."],
        "user_decision": "2026-11-13 'Yes, add it' to 'Add the 7.8 mm smoke margin to uncertainty_v5.yaml?'",
        "source": "results/v5/BASE/delta_newdata_freeze_smoke/summary.csv row smoke_mixtures,bulk_designated margin_alpha10"}
    assert doc["margins"][KEY]["applied_margin"] == FLASH_MARGIN
    assert doc["margins"][SMOKE_KEY]["applied_margin"] == SMOKE_MARGIN
    add_freezing(doc)
    add_smoke_cert(doc)
    doc["version"] = "5.0-S1+Q2+SM+FM+SC"
    doc["generated_by"] = "scripts/v5/base/build_uncertainty_v5_q2.py"
    head = ("# GENERATED - do not edit. Produced by scripts/v5/base/build_uncertainty_v5_q2.py from results/v5/S1/uncertainty_v5.yaml.\n"
            f"# Changes vs S1: {KEY}.applied_margin (3.698768343688585 -> {new}); new entries {SMOKE_KEY} (7.8009965628100595) and {FREEZE_KEY}. Read applied_margin only.\n")
    OUT.write_text(head + yaml.safe_dump(doc, sort_keys=False, width=140), encoding="utf-8")
    print("wrote", OUT)


if __name__ == "__main__":
    main()
