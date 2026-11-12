"""v5 BASE Task 10: flash input package (LFL overlay, mesitylene Antoine set) against the frozen scoreboard.

Lab mixtures are re-predicted per row through RouteCtx overrides (in memory; registry files untouched) for four variants and the lab
flash margin is recomputed with one_sided_margin, with and without the mesitylene system.  CRC 647 r3 (overlay) comes from the committed
L1/F_r3 per-fuel file, T2 rows from the frozen scoreboard.  Output: results/v5/BASE/delta_flash_inputs/{delta.csv, lab_margin.csv, readings.json}
"""
from __future__ import annotations

import importlib.util
import json
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "src"))
from fuelblend.baselines.flash_unifac import flash_point_unifac_detail  # noqa: E402
from fuelblend.eval.cellkit import frame, make_row  # noqa: E402
from fuelblend.eval.composition_route import RouteCtx  # noqa: E402
from fuelblend.uq.conformal import coverage, one_sided_margin  # noqa: E402

OUT = ROOT / "results" / "v5" / "BASE" / "delta_flash_inputs"
V4_MARGIN = 3.698768343688585
VARIANTS = ["designated", "lfl_overlay", "mesitylene_fixed", "lfl_overlay+mesitylene_fixed"]
over5 = lambda e: e > 5.0                                                                   # noqa: E731

_spec = importlib.util.spec_from_file_location("cells_t2_flash_mod", ROOT / "scripts/v5/base/cells_t2_flash.py")
t2f = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(t2f)


def lab_predictions() -> pd.DataFrame:
    lab = pd.read_csv(ROOT / "results/v4/F/flash_lab_loso.csv")
    lfl, mes = t2f.lfl_overlay_dict(), t2f.mesitylene_dict()
    ov = {"designated": {}, "lfl_overlay": dict(lfl_override=lfl), "mesitylene_fixed": dict(antoine_override=mes),
          "lfl_overlay+mesitylene_fixed": dict(lfl_override=lfl, antoine_override=mes)}
    smi = [json.loads(s) for s in lab.comp_smiles]
    xs = [json.loads(s) for s in lab.x_mole]
    out = lab[["system_id", "measured_C"]].copy()
    for v, kw in ov.items():
        ctx = RouteCtx(None, None, **kw)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            out[v] = [flash_point_unifac_detail(s, x, registry=ctx.flash_reg, use_unifac=False)["flash_C"] for s, x in zip(smi, xs)]
    mes_sys = set(pd.read_csv(ROOT / "results/v5/S2/mesitylene_flash_rerun.csv").query("has_mesitylene").system_id)
    out["has_mesitylene"] = out.system_id.isin(mes_sys)
    return out


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    sb = pd.read_csv(ROOT / "results/v5/BASE/scoreboard.csv")
    sf = sb[(sb.property == "flash_point_C") & (sb.model == "ideal_le_chatelier")]
    lab = lab_predictions()
    rows, margins = [], []
    for v in VARIANTS:
        err = lab[v] - lab.measured_C
        for scope, m in (("all_31_systems", np.ones(len(lab), bool)), ("leave_mesitylene_system_out", ~lab.has_mesitylene.to_numpy())):
            marg = one_sided_margin(err[m], lab.system_id[m], 0.10, "over")
            margins.append(dict(variant=v, scope=scope, margin=marg, n_groups=int(lab.system_id[m].nunique()), v4_margin=V4_MARGIN,
                                diff_vs_v4=marg - V4_MARGIN))
        marg_all = [r["margin"] for r in margins if r["variant"] == v and r["scope"] == "all_31_systems"][0]
        k, n = coverage(err, lab.system_id, marg_all, direction="over")
        rows.append(make_row("flash_point_C", "T1_lab", "lab_mixtures_31_systems", "ideal_le_chatelier", v, "recomputed_v5", "dev",
                             "results/v4/F/flash_lab_loso.csv (re-predicted via RouteCtx)", err=err, groups=lab.system_id, unsafe=over5,
                             margin90=marg_all, coverage90=k / n))
    # CRC 647 r3: overlay per-fuel (committed L1 run) against the frozen designated aggregates
    f = pd.read_csv(ROOT / "results/v5/L1/F_r3/flash_real_fuels.csv")
    f = f[f.statistics_allowed.astype(bool)]
    for lab_, mask in (("dev", f.split_role != "locked"), ("fold_b", f.split_role == "locked")):
        g = f[mask]
        rows.append(make_row("flash_point_C", "CRC647_r3", "crc647_r3", "ideal_le_chatelier", "lfl_overlay", "aggregated_v5", lab_,
                             "results/v5/L1/F_r3/flash_real_fuels.csv", err=g.ideal_err_C, groups=g.fuel_id, unsafe=over5))
    # T2: frozen scoreboard variants, copied
    t2 = sf[(sf.tier == "T2_measured") & sf.variant.isin(VARIANTS)]
    base_rows = pd.concat([pd.DataFrame(rows), t2], ignore_index=True)
    d = frame(base_rows.to_dict("records"))
    keys = ["tier", "set", "label"]
    base = sf[(sf.variant == "designated")][["tier", "set", "label", "mae", "bias", "unsafe_rate"]].rename(
        columns={"mae": "baseline_mae", "bias": "baseline_bias", "unsafe_rate": "baseline_unsafe_rate"})
    # lab designated baseline = the re-predicted designated row itself (frozen scoreboard row is an aggregate of the same file)
    d = d.merge(base, on=keys, how="left")
    dm = d.variant == "designated"
    d.loc[dm & d.baseline_mae.isna(), ["baseline_mae", "baseline_bias", "baseline_unsafe_rate"]] = d.loc[dm & d.baseline_mae.isna(), ["mae", "bias", "unsafe_rate"]].to_numpy()
    d.to_csv(OUT / "delta.csv", index=False)
    mt = pd.DataFrame(margins)
    mt.to_csv(OUT / "lab_margin.csv", index=False)

    o = d[(d.variant == "lfl_overlay") & d.baseline_unsafe_rate.notna()]
    mm = mt.set_index(["variant", "scope"]).margin
    r1 = bool((o.unsafe_rate <= o.baseline_unsafe_rate + 1e-9).all())
    r2 = {v: float(abs(mm[(v, "leave_mesitylene_system_out")] - mm[(v, "all_31_systems")])) for v in ("designated", "lfl_overlay")}
    r3 = float(mm[("lfl_overlay", "all_31_systems")] - V4_MARGIN)
    res = dict(reading1_unsafe_rate_never_rises=r1, reading2_loo_mesitylene_margin_diff=r2, reading2_pass=all(x <= 0.5 for x in r2.values()),
               reading3_overlay_margin_minus_v4=r3, reading3_pass=abs(r3) < 0.5,
               lab_margin_overlay=float(mm[("lfl_overlay", "all_31_systems")]), v4_margin=V4_MARGIN)
    (OUT / "readings.json").write_text(json.dumps(res, indent=2))
    print(mt.to_string())
    print(d[["tier", "set", "variant", "label", "n_groups", "bias", "mae", "unsafe_rate", "baseline_unsafe_rate", "margin90"]].to_string())
    print(json.dumps(res, indent=2))


if __name__ == "__main__":
    main()
