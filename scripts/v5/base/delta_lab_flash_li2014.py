"""v5 BASE Task 10 (U2): measured n-alkane mixture flash points (Li et al. 2014, Ha and Lee 2016) scored with the lab-flash route.

Per-row prediction is the one in delta_flash_inputs.lab_predictions (ideal Le Chatelier, use_unifac=False), designated and LFL-overlay
variants.  New points are labelled dev.  Output: results/v5/BASE/delta_lab_flash_li2014/{delta.csv, lab_margin.csv, readings.json}
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
from fuelblend.data.loaders_li2014 import load_li2014  # noqa: E402
from fuelblend.eval.composition_route import RouteCtx  # noqa: E402
from fuelblend.uq.conformal import one_sided_margin  # noqa: E402

OUT = ROOT / "results" / "v5" / "BASE" / "delta_lab_flash_li2014"
VARIANTS = ["designated", "lfl_overlay"]


def _load(name: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / "scripts/v5/base" / f"{name}.py")
    m = importlib.util.module_from_spec(spec)
    sys.modules[name] = m
    spec.loader.exec_module(m)
    return m


def predict(smi_json: pd.Series, x_json: pd.Series) -> pd.DataFrame:
    t2f = _load("cells_t2_flash")
    ov = {"designated": {}, "lfl_overlay": dict(lfl_override=t2f.lfl_overlay_dict())}
    smi = [json.loads(s) for s in smi_json]
    xs = [json.loads(s) for s in x_json]
    out = {}
    for v, kw in ov.items():
        ctx = RouteCtx(None, None, **kw)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            out[v] = [flash_point_unifac_detail(s, x, registry=ctx.flash_reg, use_unifac=False)["flash_C"] for s, x in zip(smi, xs)]
    return pd.DataFrame(out)


def is_duplicate(row, old: pd.DataFrame) -> bool:
    """Same system_id, component-aligned |dx|<=0.005 and |dflash|<=0.5 C."""
    a = dict(zip(json.loads(row.comp_smiles), json.loads(row.x_mole)))
    for _, q in old[old.system_id == row.system_id].iterrows():
        b = dict(zip(json.loads(q.comp_smiles), json.loads(q.x_mole)))
        if set(a) == set(b) and max(abs(a[k] - b[k]) for k in a) <= 0.005 and abs(row.measured_C - q.measured_C) <= 0.5:
            return True
    return False


def summarise(err: np.ndarray) -> dict:
    return dict(n=len(err), bias=float(np.mean(err)), mae=float(np.mean(np.abs(err))), max_over=float(np.max(err)),
                unsafe_rate_gt5=float(np.mean(err > 5.0)))


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    old = pd.read_csv(ROOT / "results/v4/F/flash_lab_loso.csv")
    new_all = load_li2014()
    dup = np.array([is_duplicate(r, old) for r in new_all.itertuples()])
    new = new_all[~dup].reset_index(drop=True)
    print(f"transcribed {len(new_all)}  duplicates {int(dup.sum())}  new {len(new)}")
    if len(new) == 0:
        (OUT / "readings.json").write_text(json.dumps(dict(R_LI_1="not evaluable: no new points", R_LI_2="not evaluable: no new points"), indent=2))
        return
    po, pn = predict(old.comp_smiles, old.x_mole), predict(new.comp_smiles, new.x_mole)
    mes_sys = set(pd.read_csv(ROOT / "results/v5/S2/mesitylene_flash_rerun.csv").query("has_mesitylene").system_id)
    old = old.assign(has_mesitylene=old.system_id.isin(mes_sys), label="dev_old")
    new = new.assign(has_mesitylene=False, label="dev")
    allr = pd.concat([old[["system_id", "measured_C", "has_mesitylene", "label"]], new[["system_id", "measured_C", "has_mesitylene", "label"]]],
                     ignore_index=True)
    for v in VARIANTS:
        allr[v] = np.concatenate([po[v], pn[v]])
        allr["err_" + v] = allr[v] - allr.measured_C
    allr["source"] = list(old.source) + list(new.source)
    allr["system_label"] = [""] * len(old) + list(new.system_label)
    allr.to_csv(OUT / "points.csv", index=False)

    rows, margins = [], []
    scopes = {
        "old_31_systems": (allr.label == "dev_old"),
        "old_minus_mesitylene": (allr.label == "dev_old") & ~allr.has_mesitylene,
        "old_plus_new": pd.Series(True, index=allr.index),
        "old_plus_new_minus_mesitylene": ~allr.has_mesitylene,
        "new_only": (allr.label == "dev"),
    }
    for v in VARIANTS:
        e = allr["err_" + v]
        for scope, m in scopes.items():
            g = allr.system_id[m]
            marg = one_sided_margin(e[m], g, 0.10, "over")
            margins.append(dict(variant=v, scope=scope, margin=marg, n_groups=int(g.nunique()), n_points=int(m.sum())))
            rows.append(dict(variant=v, scope=scope, label="dev", **summarise(e[m].to_numpy())))
        for sysid, g in allr[allr.label == "dev"].groupby("system_id"):
            rows.append(dict(variant=v, scope=f"new_system:{sysid}", label="dev", **summarise(g["err_" + v].to_numpy())))
    pd.DataFrame(rows).to_csv(OUT / "delta.csv", index=False)
    mt = pd.DataFrame(margins)
    mt.to_csv(OUT / "lab_margin.csv", index=False)
    mm = mt.set_index(["variant", "scope"]).margin
    en = allr.loc[allr.label == "dev", "err_lfl_overlay"].to_numpy()
    r1 = dict(bias=float(en.mean()), abs_bias_le_2=bool(abs(en.mean()) <= 2.0), over5_rate=float((en > 5).mean()), over5_le_10pct=bool((en > 5).mean() <= 0.10))
    r1["pass"] = r1["abs_bias_le_2"] and r1["over5_le_10pct"]
    r2 = {v: float(abs(mm[(v, "old_plus_new_minus_mesitylene")] - mm[(v, "old_plus_new")])) for v in VARIANTS}
    res = dict(n_transcribed=len(new_all), n_duplicates=int(dup.sum()), n_new=len(new), n_new_systems=int(new.system_id.nunique()),
               R_LI_1=r1, R_LI_2_diff_with_vs_without_mesitylene=r2, R_LI_2_pass_overlay=bool(r2["lfl_overlay"] < 0.5),
               margins={f"{v}|{s}": float(mm[(v, s)]) for v, s in mm.index})
    (OUT / "readings.json").write_text(json.dumps(res, indent=2))
    print(mt.to_string())
    print(pd.DataFrame(rows).to_string())
    print(json.dumps(res, indent=2))


if __name__ == "__main__":
    main()
