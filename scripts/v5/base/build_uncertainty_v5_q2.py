"""Generate configs/uncertainty_v5.yaml: results/v5/S1/uncertainty_v5.yaml with ONLY the flash lab ideal applied_margin replaced
by 2.206932808388615 (user decision Q2, 2026-11-13; row designated,old_plus_new_minus_mesitylene of delta_lab_flash_li2014/lab_margin.csv)."""
from pathlib import Path
import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parents[3]
SRC = ROOT / "results/v5/S1/uncertainty_v5.yaml"
OUT = ROOT / "configs/uncertainty_v5.yaml"
CSV = ROOT / "results/v5/BASE/delta_lab_flash_li2014/lab_margin.csv"
KEY = "flash_point_C.lab_mixtures.ideal"


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
    doc["version"] = "5.0-S1+Q2"
    doc["generated_by"] = "scripts/v5/base/build_uncertainty_v5_q2.py"
    head = ("# GENERATED - do not edit. Produced by scripts/v5/base/build_uncertainty_v5_q2.py from results/v5/S1/uncertainty_v5.yaml.\n"
            f"# Only {KEY}.applied_margin differs (3.698768343688585 -> {new}). Read applied_margin only.\n")
    OUT.write_text(head + yaml.safe_dump(doc, sort_keys=False, width=140), encoding="utf-8")
    print("wrote", OUT)


if __name__ == "__main__":
    main()
