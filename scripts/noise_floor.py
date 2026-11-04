"""Noise-floor estimate per property (WP2). Upper-bound proxy: observed min-max
spread of independent measurements for records with >=3 observations."""
import io, os, sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
import pandas as pd

DEDUP = "data/unified/measurements_dedup.csv"
PROPS = ["viscosity_dyn_mPas", "viscosity_kin_cSt", "density_kg_m3", "flash_point_C",
         "vapor_pressure_kPa", "pour_point_C", "freezing_point_C", "mon", "dcn",
         "smoke_point_mm", "net_heat_content_btu_lb"]


def main():
    df = pd.read_csv(DEDUP, low_memory=False)
    df["spread_v"] = df["vmax"] - df["vmin"]
    lines = [
        "# NOISE_FLOOR.md (WP2, auto-generated)",
        "",
        "Estimate: for records with >=3 observations, the observed min-max spread of independent",
        "measurements, median across records. Upper-bound proxy for irreducible noise (remaining",
        "T/P/composition differences within bins included). Model errors should be reported next to",
        "this floor. ASTM precision statements are added manually per method.",
        "",
        "| property | n(>=3 obs) | median spread | p90 spread | conflict rate | outlier rate |",
        "|---|---|---|---|---|---|",
    ]
    for prop in PROPS:
        sub = df[df["property"] == prop]
        g = sub[sub["n_obs"] >= 3]
        meds = g["spread_v"].median() if len(g) else float("nan")
        p90 = g["spread_v"].quantile(0.9) if len(g) else float("nan")
        cr = sub["conflict"].mean() * 100
        orr = (sub["n_outlier"] > 0).mean() * 100
        lines.append(f"| {prop} | {len(g)} | {meds:.4g} | {p90:.4g} | {cr:.1f}% | {orr:.1f}% |")
    out = os.path.join(os.path.dirname(__file__), "..", "reports", "NOISE_FLOOR.md")
    with io.open(out, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
