"""v5 BASE Task 15: fresh holdout (score once).  Readings and frozen model version: results/v5/BASE/delta_fresh_holdout/README.md.

Outputs: data/real_fuels/v5/nistir6659_density.csv; results/v5/BASE/delta_fresh_holdout/{label_search.csv, rows_scored.csv, delta.csv, readings.json, scored_once.json}
"""
from __future__ import annotations

import importlib.util
import json
import re
import sys
import warnings
from datetime import date
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "src"))
from fuelblend.data.loaders_nistir6659 import load_nistir6659  # noqa: E402
from fuelblend.eval.cellkit import frame, make_row  # noqa: E402

BASE = ROOT / "results" / "v5" / "BASE"
OUT = BASE / "delta_fresh_holdout"
V5 = ROOT / "data" / "real_fuels" / "v5"
SRC = ROOT / "extras/fuelblend_data/data/processed/nist_density_compressed_liquid.csv"
FROZEN = "5a62e20"
SEARCH_DIRS = [ROOT / "data" / "real_fuels", ROOT / "data" / "processed", ROOT / "data" / "unified"]
SKIP_NAME = re.compile(r"ajftd", re.I)               # the AJFTD tables are where the 'no property row' came from


def guard(out: Path) -> None:
    if (out / "scored_once.json").exists():
        raise SystemExit(f"refusing: {out / 'scored_once.json'} exists (this holdout is scored once)")


def label_search(labels: list[str]) -> pd.DataFrame:
    pats = {l: re.compile(rf"POSF[\s_\-]?0*{l}\b", re.I) for l in labels}
    rows = []
    for d in SEARCH_DIRS:
        for p in d.rglob("*.csv"):
            if SKIP_NAME.search(str(p.relative_to(ROOT))) or p.stat().st_size > 50e6:
                continue
            txt = p.read_text(encoding="utf-8", errors="ignore")
            for l, pat in pats.items():
                n = len(pat.findall(txt))
                if n:
                    rows.append(dict(label=l, file=str(p.relative_to(ROOT)).replace("\\", "/"), n_matches=n))
    return pd.DataFrame(rows, columns=["label", "file", "n_matches"])


def main() -> None:
    guard(OUT)
    if not SRC.exists():
        print("SKIPPED: extras/ source file absent")
        return
    warnings.simplefilter("ignore")
    OUT.mkdir(parents=True, exist_ok=True)
    nd = load_nistir6659(SRC)
    nd.to_csv(V5 / "nistir6659_density.csv", index=False)
    fu = pd.read_csv(BASE / "t2_fuels.csv")
    blank = fu[(fu.source == "ajftd") & (fu.n_target_rows == 0) & ~fu.low_coverage.astype(bool)]
    ls = label_search(sorted(set(blank.label.astype(str))))
    ls.to_csv(OUT / "label_search.csv", index=False)

    # POSF 4658 = AJFTD:23 (already a T2 fuel; new temperatures, not a new fuel)
    x = nd[(nd.posf == "4658") & (nd.P_kPa <= 1000.0)].copy()
    ent = fu[(fu.source == "ajftd") & (fu.label.astype(str) == "4658")].iloc[0]
    grp = pd.read_csv(ROOT / "results/v5/S3/corrected_2026-11-12_groups.csv")
    members = grp[grp.members.astype(str).str.contains(ent.entity + r"\b", regex=True)]
    from fuelblend.baselines.flash_unifac import _lookup
    from fuelblend.data.reconstruction_r3 import SubReps
    from fuelblend.eval.composition_route import RouteCtx, kind_bins, map_bins, predict_properties
    from fuelblend.mixture.common import CompoundRegistry, PurePropertyProvider
    reg = CompoundRegistry(ROOT / "data/unified/compounds.csv")
    prov = PurePropertyProvider(reg, ROOT / "data/cache/pure_grid.npz", verbose=False)
    reps = SubReps(reg, prov, vp_lookup=lambda s: _lookup(s, None))
    c = pd.read_csv(V5 / "ajftd_gcxgc_composition.csv")
    mf = map_bins(kind_bins("ajftd", c[c.fuel_id == int(ent.fuel_id)]), reps)
    p = predict_properties(RouteCtx(reg, prov), mf, pd.DataFrame(dict(property="density_kg_m3", T_K=x.T_K.to_numpy())))
    x["pred"], x["status"] = p.pred.to_numpy(), p.status.to_numpy()
    x["err"] = x.pred - x.value
    x["entity"], x["group"] = ent.entity, ent.group
    x.to_csv(OUT / "rows_scored.csv", index=False, float_format="%.4f")
    ok = x[(x.status == "ok") & np.isfinite(x.err)]

    bk = pd.read_csv(BASE / "rows_t2_bulk.csv")
    _s = importlib.util.spec_from_file_location("cells_t2_bulk_fh", ROOT / "scripts/v5/base/cells_t2_bulk.py")
    bulk = importlib.util.module_from_spec(_s)
    _s.loader.exec_module(bulk)
    ref = bk[(bk.property == "density_kg_m3") & (bk.variant == "designated") & (np.round(bk.T_K, 2) == 288.15) & (bk.status == "ok") & ~bk.low_coverage.astype(bool)]
    ref = bulk.one_per_group(ref)
    e_ref = ref.pred - ref.measured
    lo, hi = float(np.percentile(e_ref, 2.5)), float(np.percentile(e_ref, 97.5))
    r = make_row("density_kg_m3", "T2_measured", "nistir6659_posf4658_P_le_1000kPa", "ideal_volume", "designated", "recomputed_v5", "exploratory",
                 "data/real_fuels/v5/nistir6659_density.csv", err=ok.err, groups=ok.group,
                 extra=dict(frozen_model=FROZEN, n_T=len(ok), T_range_K=f"{ok.T_K.min():.0f}-{ok.T_K.max():.0f}",
                            t2_descriptive_interval_lo=lo, t2_descriptive_interval_hi=hi,
                            baseline_mae=float(np.abs(e_ref).mean()), baseline_bias=float(e_ref.mean()), baseline_unsafe_rate=np.nan))
    frame([r]).to_csv(OUT / "delta.csv", index=False)
    res = dict(n_fresh_fuels=0, n_fuels_scored=1, note="POSF 4658 is AJFTD:23, already read at 288.15 K; new temperatures only",
               n_rows=len(ok), bias=r["bias"], mae=r["mae"], interval=[lo, hi], reading_pass=bool(lo <= r["bias"] <= hi),
               bias_by_T={f"{t:.0f}": float(e) for t, e in zip(ok.T_K, ok.err)},
               double_count_groups_with_this_fuel=members.group.tolist(), label_search_hits=int(len(ls)),
               labels_searched=int(blank.label.nunique()))
    (OUT / "readings.json").write_text(json.dumps(res, indent=2))
    (OUT / "scored_once.json").write_text(json.dumps(dict(date=str(date.today()), frozen_model=FROZEN), indent=2))
    print(json.dumps(res, indent=2))
    print(ls.to_string())


if __name__ == "__main__":
    main()
