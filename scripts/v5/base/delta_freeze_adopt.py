"""delta_freeze_adopt (results/v5/BASE/delta_freeze_adopt/README.md): ideal and UNIFAC SvL margins re-derived on the new constants / data.

Sub-readings: (A) old SLE file + old constants (d9309d3), must reproduce delta_freeze_unifac R-FU-1; (B) same 67 rows, new constants;
(C) full new scope incl. eutectic points (headline); (C-ps) B + the Pitzer & Scott eutectics only; R-FA-2 descriptive.

    .venv/Scripts/python.exe scripts/v5/base/delta_freeze_adopt.py
Writes only into results/v5/BASE/delta_freeze_adopt/ (never configs/).
"""
from __future__ import annotations

import importlib.util
import json
import subprocess
import sys
import tempfile
import warnings
from pathlib import Path
from unittest import mock

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "src"))
warnings.simplefilter("ignore")
OUT = ROOT / "results" / "v5" / "BASE" / "delta_freeze_adopt"

from fuelblend.baselines.flash_unifac import make_unifac_gamma_fn  # noqa: E402
from fuelblend.eval import freeze_multicomp as FM  # noqa: E402
from fuelblend.eval.composition_route import build_tmap_verified  # noqa: E402


def _load(name, rel):
    s = importlib.util.spec_from_file_location(name, ROOT / rel)
    m = importlib.util.module_from_spec(s)
    s.loader.exec_module(m)
    return m


FMG = _load("fm_margin", "scripts/v5/base/delta_freeze_margin.py")
OLDM = FMG.OLDM
SLE_REL = "data/processed/v5/thermoml_sle_hc.csv"
TM_REL = "data/raw/sle_constants/tm_dhfus_verified.csv"
OLD_SLE_REV = "1f9b89d~1"      # SLE file before the Pitzer & Scott eutectics
OLD_TM_REV = "d9309d3"         # constants before Tasks 0-3 (bicyclohexyl 6.86, no n-nonane / 2,3-DMB / m- / p-xylene)
PS_FILE = "calorimetry_primary:PS43_TableX"
REF_A = dict(margin_ideal=16.9755, margin_unifac=13.8427, mae_ideal=3.3271, mae_unifac=2.6238)   # delta_freeze_unifac/summary.csv S_F


def _git_show(rev, rel, dst):
    dst.write_bytes(subprocess.run(["git", "show", f"{rev}:{rel}"], cwd=ROOT, check=True, capture_output=True).stdout)
    return dst


def _xkey(r):
    x = tuple(round(float(v), 6) for v in json.loads(r.x_liquid)) if isinstance(r.x_liquid, str) else None
    return (r.system_id, round(float(r.T_K), 4), x)


def prepare_sle(src, dst, kinds):
    """Copy of the SLE csv with the kept rows: liquidus, plus (if 'eutectic' in kinds) eutectic rows relabelled 'liquidus' (orig_kind kept),
    except a eutectic identical (system_id, T_K, x_liquid) to a liquidus row (one row per point).  Returns the kept frame."""
    d = pd.read_csv(src)
    d["orig_kind"] = d.point_kind
    liq = d[d.point_kind == "liquidus"]
    keep = [liq]
    if "eutectic" in kinds:
        lk = {_xkey(r) for r in liq.itertuples()}
        eu = d[(d.point_kind == "eutectic") & d.x_liquid.notna()]       # 2 ThermoML eutectics (T > 321 K, out of S_F) carry no composition
        eu = eu[[_xkey(r) not in lk for r in eu.itertuples()]].copy()
        if kinds.get("eutectic") == "ps_only":
            eu = eu[eu.file == PS_FILE]
        eu["point_kind"] = "liquidus"
        keep.append(eu)
    out = pd.concat(keep).sort_index()
    out.to_csv(dst, index=False)
    return out.reset_index(drop=True)


def score_rows(sle_csv, tmap):
    """delta_freeze_unifac.main()'s row loop on (sle_csv, tmap): scope_rows (ideal) + UNIFAC prediction and fallback flag."""
    with mock.patch.object(FMG, "SLE_CSV", Path(sle_csv)):
        rows = FMG.scope_rows(tmap)
    sle = pd.read_csv(sle_csv)
    sle = sle[sle.point_kind == "liquidus"].reset_index(drop=True)
    assert len(rows) == len(sle)
    pu, fbs = [], []
    for r, s in zip(rows.itertuples(), sle.itertuples()):
        if not r.all_verified:
            pu.append(np.nan)
            fbs.append(False)
            continue
        smi, iks = json.loads(s.smiles), json.loads(s.inchikeys)
        x = [float(v) for v in json.loads(s.x_liquid)]
        x = list(np.asarray(x) / np.sum(x))
        cc = [tmap[k] for k in iks]
        cc = [(float(c["tm_K"]), float(c["dh_fus_J_mol"])) for c in cc]
        gfn, fb = make_unifac_gamma_fn(smi)
        pu.append(FM.svl_activity_C(x, cc, gfn))
        fbs.append(bool(fb))
    rows["pred_unifac"] = pu
    rows["unifac_fallback"] = fbs
    rows["orig_kind"] = sle.orig_kind.to_numpy() if "orig_kind" in sle else "liquidus"
    rows["file"] = sle.file.to_numpy()
    rows["x_liquid"] = sle.x_liquid.to_numpy()
    return rows


def pair(d, scope):
    return [OLDM.summarize(d, scope, m, c, "under") for m, c in (("svl_verified", "pred_verified"), ("svl_unifac", "pred_unifac"))]


def reading(S, scope, d):
    a = next(r for r in S if r["scope"] == scope and r["model"] == "svl_verified")
    b = next(r for r in S if r["scope"] == scope and r["model"] == "svl_unifac")
    nfb = int(d.unifac_fallback.sum())
    pref = bool(b["margin_alpha10"] < a["margin_alpha10"] and b["mae"] <= a["mae"] and nfb == 0)
    return dict(margin_ideal_C=float(a["margin_alpha10"]), margin_unifac_C=float(b["margin_alpha10"]), mae_ideal_C=float(a["mae"]),
                mae_unifac_C=float(b["mae"]), bias_ideal_C=float(a["bias"]), bias_unifac_C=float(b["bias"]), n_systems=int(b["n_systems"]),
                n_rows=int(b["n_scored"]), unifac_fallbacks=nfb, rule_prefers_unifac=pref,
                systems=sorted(d.group.unique().tolist()))


def per_system(sf):
    sf = sf.copy()
    sf["r_v"] = sf.pred_verified - sf.measured
    sf["r_u"] = sf.pred_unifac - sf.measured
    return sf.groupby(["group", "sample"]).agg(n=("r_v", "size"), n_eutectic=("orig_kind", lambda v: int((v == "eutectic").sum())),
                                               worst_ideal=("r_v", "min"), worst_unifac=("r_u", "min"),
                                               mae_ideal=("r_v", lambda v: v.abs().mean()), mae_unifac=("r_u", lambda v: v.abs().mean()),
                                               fallback=("unifac_fallback", "any")).sort_values("worst_ideal")


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    tmp = Path(tempfile.mkdtemp(prefix="freeze_adopt_"))
    old_sle = _git_show(OLD_SLE_REV, SLE_REL, tmp / "sle_old.csv")
    tmap_old = build_tmap_verified(_git_show(OLD_TM_REV, TM_REL, tmp / "tm_old.csv"))
    tmap_new = build_tmap_verified()
    assert abs(tmap_old["WVIIMZNLDWSIRH-UHFFFAOYSA-N"]["dh_fus_J_mol"] - 6860.0) < 1e-6
    assert abs(tmap_new["WVIIMZNLDWSIRH-UHFFFAOYSA-N"]["dh_fus_J_mol"] - 13940.0) < 1e-6

    # (A) old file, liquidus only, old constants
    prepare_sle(old_sle, tmp / "A.csv", {})
    rA = score_rows(tmp / "A.csv", tmap_old)
    sfA = rA[rA.in_SF].copy()
    assert len(sfA) == 67 and sfA.group.nunique() == 9, (len(sfA), sfA.group.nunique())
    # (B) the same 67 row identities, new constants (old file rows are a prefix of the new file)
    rB_all = score_rows(tmp / "A.csv", tmap_new)
    sfB = rB_all[rA.in_SF.to_numpy()].copy()
    assert len(sfB) == 67 and bool(sfB.in_SF.all())
    # (C) new file, liquidus + eutectic, new constants
    prepare_sle(ROOT / SLE_REL, tmp / "C.csv", {"eutectic": "all"})
    rC = score_rows(tmp / "C.csv", tmap_new)
    sfC = rC[rC.in_SF].copy()
    # (C-ps) B rows + Pitzer & Scott eutectics only
    prepare_sle(ROOT / SLE_REL, tmp / "Cps.csv", {"eutectic": "ps_only"})
    rP = score_rows(tmp / "Cps.csv", tmap_new)
    sfP = rP[rP.in_SF].copy()
    allv = rC[rC.all_verified & rC.all_hc].copy()

    S = pair(sfC, "S_F") + pair(sfA, "S_F_A_old") + pair(sfB, "S_F_B_newconst") + pair(sfP, "S_F_Cps") + \
        pair(allv, "thermoml_all_verified_any_T_any_nC")
    summ = pd.DataFrame(S)
    summ.to_csv(OUT / "summary.csv", index=False)   # full precision: the yaml margin is read from here
    rC.to_csv(OUT / "rows.csv", index=False, float_format="%.5f")
    sfA.to_csv(OUT / "rows_A.csv", index=False, float_format="%.5f")
    sfB.to_csv(OUT / "rows_B.csv", index=False, float_format="%.5f")
    per_system(sfC).to_csv(OUT / "per_system_S_F.csv", float_format="%.4f")

    rd_A, rd_B, rd_C, rd_P = (reading(S, "S_F_A_old", sfA), reading(S, "S_F_B_newconst", sfB), reading(S, "S_F", sfC),
                              reading(S, "S_F_Cps", sfP))
    rd_R2 = reading(S, "thermoml_all_verified_any_T_any_nC", allv)
    rd_R2.pop("systems")
    a_ok = (abs(rd_A["margin_ideal_C"] - REF_A["margin_ideal"]) < 1e-3 and abs(rd_A["margin_unifac_C"] - REF_A["margin_unifac"]) < 1e-3
            and abs(rd_A["mae_ideal_C"] - REF_A["mae_ideal"]) < 1e-3 and abs(rd_A["mae_unifac_C"] - REF_A["mae_unifac"]) < 1e-3)
    stop = None
    if not a_ok:
        stop = "sub-reading A does not reproduce delta_freeze_unifac R-FU-1; attribution invalid"
    elif rd_C["n_systems"] < 9:
        stop = f"n_systems {rd_C['n_systems']} < 9"
    elif not rd_C["rule_prefers_unifac"]:
        stop = "preference rule does not prefer UNIFAC on (C): user decision needed (decision was made on 13.84 vs 16.98)"
    new_sys = sorted(set(sfC.group) - set(sfA.group))
    readings = {k: rd_C[k] for k in ("margin_ideal_C", "margin_unifac_C", "mae_ideal_C", "mae_unifac_C", "n_systems", "n_rows",
                                     "unifac_fallbacks", "rule_prefers_unifac")}
    readings.update(stop_reason=stop, scope="S_F (C): liquidus + eutectic points, new constants", new_systems=new_sys,
                    n_rows_old_67_set=len(sfA), n_eutectic_rows_C=int((sfC.orig_kind == "eutectic").sum()),
                    A_reproduces_R_FU_1=bool(a_ok), A_old=rd_A, B_newconst=rd_B, C=rd_C, C_ps=rd_P, R_FA_2=rd_R2)
    (OUT / "readings.json").write_text(json.dumps(readings, indent=1), encoding="utf-8")
    pd.set_option("display.width", 250)
    print(summ[["scope", "model", "n_scored", "n_systems", "mae", "bias", "margin_alpha10", "frac_unsafe", "max_group_unsafe"]].to_string())
    print(per_system(sfC).to_string())
    print(json.dumps({k: v for k, v in readings.items() if k not in ("A_old", "B_newconst", "C", "C_ps")}, indent=1))


if __name__ == "__main__":
    main()
