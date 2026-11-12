"""v5 delta_diesel_surrogates: Mueller 2016 diesel surrogates (known composition, mole %) through predict_fuel.

Declared readings: results/v5/BASE/delta_diesel_surrogates/README.md. Descriptive only, n = 4, no skill claim.
Outputs: rows.csv (per surrogate x property) and target_vs_surrogate.csv in results/v5/BASE/delta_diesel_surrogates/.
"""
from __future__ import annotations

import subprocess
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "src"))
from fuelblend.api.predict_route import predict_fuel  # noqa: E402

OUT = ROOT / "results" / "v5" / "BASE" / "delta_diesel_surrogates"
T_K = 293.15            # density at 20 C; flash and cloud (freezing) are T-independent in the API
PROPS = [("flash_point_C", "flash_point_C"), ("density_kg_m3", "density_kg_m3"), ("cloud_point_C", "freezing_point_C")]
CLOUD_LABEL = "SvL vs cloud point, different property, indicative only"


def measured_table() -> dict:
    p = pd.read_csv(ROOT / "data/raw/surrogates/mueller2016_properties.csv")
    return {(c, pr): g.value.tolist() for (c, pr), g in p.groupby(["fuel_code", "property"])}


def main() -> None:
    warnings.simplefilter("ignore")
    OUT.mkdir(parents=True, exist_ok=True)
    meas = measured_table()
    d = pd.read_csv(ROOT / "data/raw/surrogates/fuel_surrogates.csv")
    d = d[d.fuel_code.str.startswith("Mueller2016_")]
    rows = []
    for code, g in d.groupby("fuel_code", sort=False):
        res = predict_fuel(list(zip(g.component_smiles, g.frac_value)), T_K=T_K, fuel_type="diesel",
                           mole_fractions=list(g.frac_value), composition_kind="known")
        for prop, key in PROPS:
            r = res[key]
            cov = r["coverage"]
            m = meas.get((code, prop), [])
            mv = float(np.mean(m)) if m else np.nan
            val = r["value"]
            err = (val - mv) if (val is not None and m) else np.nan
            rows.append(dict(
                fuel_code=code, property=prop, value=val, status=r["status"], method=r["method"],
                coverage_own_mass=cov["own_mass"], coverage_substituted_mass=cov["substituted_mass"],
                coverage_unmapped_mass=cov["unmapped_mass"], lo=r["lo"], hi=r["hi"], in_domain=r["in_domain"],
                measured=mv, n_measured=len(m), error=err,
                rel_error_pct=(100 * err / mv if prop == "density_kg_m3" and np.isfinite(err) else np.nan),
                label=CLOUD_LABEL if prop == "cloud_point_C" else "",
                warnings="; ".join(r["warnings"])))
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "rows.csv", index=False)

    # surrogate-predicted vs target diesel CF measured (density: mean of the two laboratory values, both listed)
    cf = {pr: meas[("Mueller2016_CF", pr)] for pr in ("flash_point_C", "density_kg_m3", "cloud_point_C")}
    t = []
    for r in df.itertuples():
        tgt = float(np.mean(cf[r.property]))
        t.append(dict(fuel_code=r.fuel_code, property=r.property, surrogate_predicted=r.value, status=r.status,
                      target_CF_measured=tgt, target_CF_values=";".join(str(v) for v in cf[r.property]),
                      target_basis="mean of two lab values" if len(cf[r.property]) > 1 else "single measurement",
                      diff_pred_minus_target=(r.value - tgt) if r.value is not None else np.nan,
                      surrogate_measured=r.measured, surrogate_measured_minus_target=r.measured - tgt, label=r.label))
    pd.DataFrame(t).to_csv(OUT / "target_vs_surrogate.csv", index=False)
    try:
        h = subprocess.check_output(["git", "rev-parse", "--short", "HEAD"], cwd=ROOT, text=True).strip()
    except Exception:
        h = "unknown"
    (OUT / "run_commit.txt").write_text(h + "\n")
    pd.set_option("display.width", 250); pd.set_option("display.max_columns", 30)
    print(df[["fuel_code", "property", "value", "status", "coverage_own_mass", "coverage_substituted_mass",
              "coverage_unmapped_mass", "measured", "error", "rel_error_pct"]].round(3).to_string())
    print("commit", h)


if __name__ == "__main__":
    main()
