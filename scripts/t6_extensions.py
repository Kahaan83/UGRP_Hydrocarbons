#!/usr/bin/env python
"""T6: extension-property decision + honest baselines for implementable ones.

    python scripts/t6_extensions.py

Decides implement/exclude for GCV, MON, DCN, freezing, pour, cloud with documented
reasons. For MON/DCN (the implementable pair): (a) null S1-median baseline and
(b) pure-anchored mass-weighted ideal mixing, both scored on S2 with a
system-grouped bootstrap, so later models have numbers to beat.
Writes reports/T6_EXTENSIONS.md.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

BASE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(BASE / "src"))
from fuelblend.eval import metrics as evm  # noqa: E402

MT = BASE / "data/processed/model_table.csv"
SPLIT_FILE = BASE / "splits/s1s2s3.json"
REPORT = BASE / "reports/T6_EXTENSIONS.md"
LOG: list[str] = []


def log(msg) -> None:
    print(msg, flush=True)
    LOG.append(str(msg))


def main() -> None:
    log("T6 extension decision - implement or exclude per property")
    mt = pd.read_csv(MT, low_memory=False)
    raw = json.loads(SPLIT_FILE.read_text())
    S = {k: set(int(i) for i in v) for k, v in raw.items() if isinstance(v, list) and v and isinstance(v[0], int)}
    mt["split"] = ""
    mt.loc[sorted(S["S1"]), "split"] = "S1"
    mt.loc[sorted(S["S2"]), "split"] = "S2"
    mt.loc[sorted(S["S3"]), "split"] = "S3"
    log(f"model_table rows={len(mt)} splits S1={len(S['S1'])} S2={len(S['S2'])} S3={len(S['S3'])}")
    rd = pd.read_csv(BASE / "data/processed/excluded_by_R1.csv", low_memory=False) \
        if (BASE / "data/processed/excluded_by_R1.csv").exists() else pd.DataFrame()

    log("")
    log("## data inventory per candidate")
    for p in ["net_heat_content", "gross_heat_of_combustion", "gcv", "mon", "dcn",
              "freezing_point_C", "pour_point_C", "cloud_point_C"]:
        rows = mt[mt["property"] == p]
        excl = len(rd[rd["property"] == p]) if len(rd) else 0
        log(f"{p:34s} pure={int((rows['kind'] == 'pure').sum()):5d} "
            f"mix={int(rows['kind'].isin(['binary', 'multi']).sum()):5d} excluded_by_R1={excl}")

    log("")
    log("## decisions")
    log("- GCV: EXCLUDE - zero pure-compound heats of combustion in open HC data (T0); the physics chain cannot start.")
    log("- freezing_point_C: EXCLUDE for now - 25 pure rows < 200-row training floor (same rule raised it in WP4).")
    log("- pour_point_C / cloud_point_C: EXCLUDE - exists only as excluded L0 real-fuel rows, no molecular anchor.")
    log("- mon / dcn: IMPLEMENT - pure and mixture rows both exist (Leenhouts), so provide null + pure-anchored baselines.")

    comp = pd.read_csv(BASE / "data/unified/compounds.csv", low_memory=False)
    for p in ["mon", "dcn"]:
        g = mt[(mt["property"] == p) & (~mt["conflict"].astype(bool)) & mt["value"].notna()]
        s1 = g[g["split"] == "S1"]
        s2 = g[g["split"] == "S2"]
        log("")
        log(f"== {p}: S1 n={len(s1)} S2 n={len(s2)} (S2 systems={s2['system_id'].nunique()}), "
            f"range [{g['value'].min():.1f}, {g['value'].max():.1f}]")
        if not len(s1) or not len(s2):
            continue
        cl = s2["system_id"].values
        med = float(np.median(s1["value"]))
        st = evm.report(s2["value"].values, np.full(len(s2), med), cl, mode="generic")
        log(f"   null baseline  (S1 median = {med:.1f}): MAE={st['mae']:.2f} RMSE={st['rmse']:.2f}")

        # pure-anchored ideal mixing: per-inchikey median pure value, mass-weighted blend
        pure = g[(g["split"] == "S1") & (g["kind"] == "pure")]
        val_ik: dict[str, float] = {}
        for r in pure.itertuples():
            try:
                iks = json.loads(r.comp_inchikeys)
                if len(iks) == 1:
                    val_ik.setdefault(str(iks[0]), []).append(float(r.value))
            except Exception:
                pass
        val_ik = {k: float(np.median(v)) for k, v in val_ik.items() if v}
        log(f"   pure anchors: {len(val_ik)} molecules (S1 pure median per inchikey)")
        mix = g[(g["split"] == "S2") & (g["kind"] != "pure")]
        preds, meas, cls = [], [], []
        miss = 0
        for r in mix.itertuples():
            try:
                iks = [str(k) for k in json.loads(r.comp_inchikeys)]
                w = json.loads(r.w_mass)
            except Exception:
                miss += 1
                continue
            if not iks or len(iks) != len(w) or any(k not in val_ik for k in iks):
                miss += 1
                continue
            w = np.asarray(w, float)
            preds.append(float(np.sum(w * np.asarray([val_ik[k] for k in iks]) / w.sum())))
            meas.append(float(r.value))
            cls.append(str(r.system_id))
        if preds:
            xx = pd.Series(meas, dtype=float)
            pp = pd.Series(preds, dtype=float)
            st2 = evm.report(xx.values, pp.values, np.asarray(cls), mode="generic")
            log(f"   ideal-mix (pure-anchored mass-weighted): n={len(preds)} skipped={miss} "
                f"MAE={st2['mae']:.2f} RMSE={st2['rmse']:.2f} "
                f"(S2 mixture systems={len(set(cls))})")
        else:
            log(f"   ideal-mix: no scoreable S2 mixture rows (skipped={miss})")

    REPORT.write_text("\n".join(LOG) + "\n", encoding="utf-8")
    print(f"\nwrote {REPORT}")


if __name__ == "__main__":
    main()
