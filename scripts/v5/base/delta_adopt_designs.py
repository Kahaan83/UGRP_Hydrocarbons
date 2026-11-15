"""delta_adopt_designs runner (declared in results/v5/BASE/delta_adopt_designs/README.md before any run).

Design re-runs under the adopted defaults of docs/decisions_2026-11-15.md: pools P0, P0x-a, P0x-b, P0y, P0y-b, P1, P2, P3 x one spec
per process.  Protocol (declared): n_starts=10, seed=42, Engler D86, lambda 0.25, alt_standin on, strict_constraints on, synthetic cap
0.50 for P2/P3 only.  ``--n-starts`` / ``--seed`` overrides are for smoke tests and are refused for the default output folder.

    python scripts/v5/base/delta_adopt_designs.py --spec astm:jet_a                         # all eight pools, resumable
    python scripts/v5/base/delta_adopt_designs.py --spec bis:jet_a1 --pool P0yb P1          # a subset (parallel processes are safe)
    python scripts/v5/base/delta_adopt_designs.py --rebuild                                 # designs.csv + compare.csv from per-design files
    python scripts/v5/base/delta_adopt_designs.py --timing --out <dir>                      # per-evaluation timing, svl_unifac vs ideal_svl

Unit of record: designs/<pool>/<design_id>__{designs,constraints,readings}.csv (written atomically).  A pool whose model build or strict
guard raises ``ModelUnavailable`` writes a ``model_unavailable`` row (a reading) and the run continues.  ``designs.csv`` and
``compare.csv`` are rebuilt from the per-design files after every pool, so an interrupted run keeps finished rows; a pool already on disk
for the spec is skipped unless ``--force``.
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import math
import os
import subprocess
import sys
import time
import traceback
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts" / "v5" / "base"))
warnings.simplefilter("ignore")

OVERLAY_ENV = "FUELBLEND_FLASH_REGISTRY_ADDITIONS"
if os.environ.get(OVERLAY_ENV) == "0":
    raise SystemExit(f"{OVERLAY_ENV}=0: the adopted default (decision 3) is the overlay ON; P0y / P0y-b need bicyclohexyl's entry")

from fuelblend.design import blendstock_v5 as B  # noqa: E402
from fuelblend.design import engine_v4 as E  # noqa: E402
from fuelblend.eval import composition_route as CR  # noqa: E402

_s = importlib.util.spec_from_file_location("ws7run", ROOT / "scripts" / "v5" / "base" / "delta_ws7_design_realism.py")
W = importlib.util.module_from_spec(_s)
_s.loader.exec_module(W)
import delta_smoke_branch_opt as SB  # noqa: E402
import delta_pool_bicyclohexyl as BC  # noqa: E402

DEFAULT_OUT = ROOT / "results" / "v5" / "BASE" / "delta_adopt_designs"
POOLS = ["P0", "P0xa", "P0xb", "P0y", "P0yb", "P1", "P2", "P3"]
ALIASES = {"P0x-a": "P0xa", "P0x-b": "P0xb", "P0y-b": "P0yb"}
SPECS = ["astm:jet_a", "bis:jet_a1"]
N_STARTS, SEED = 10, 42                                  # declared; never reduced for the results folder
RAYLEIGH_MARGINS = ROOT / "results" / "v5" / "BASE" / "delta_d86_rayleigh" / "d86_rayleigh_margins.json"
IDEAL_KEY = CR.FREEZING_MODELS["ideal_svl"]
COMMITTED = {   # pool -> (folder, design-id tag) of the committed reading compared against (README table)
    "P0": (ROOT / "results/v5/BASE/delta_ws7_design_realism/designs/P0", "P0"),
    "P1": (ROOT / "results/v5/BASE/delta_ws7_design_realism/designs/P1", "P1"),
    "P0xa": (ROOT / "results/v5/BASE/delta_smoke_branch_opt/designs/P0xa", "P0xa_standin"),
    "P0xb": (ROOT / "results/v5/BASE/delta_smoke_branch_opt/designs/P0xb", "P0xb_standin"),
    "P2": (ROOT / "results/v5/BASE/delta_smoke_branch_opt/designs/P2", "P2_standin"),
    "P3": (ROOT / "results/v5/BASE/delta_smoke_branch_opt/designs/P3", "P3_standin"),
    "P0y": (ROOT / "results/v5/BASE/delta_pool_bicyclohexyl/designs/P0y", "P0y_standin"),
    "P0yb": (ROOT / "results/v5/BASE/delta_pool_bicyclohexyl/designs/P0yb", "P0yb_standin"),
}
SPEC_LIMIT_FREEZE = {"astm:jet_a": -40.0, "bis:jet_a1": -47.0}


# ------------------------------------------------------------------------------------------------ pools / provenance
def make_pool(pool: str) -> B.BlendPool:
    if pool in ("P0", "P1", "P2", "P3"):
        return W.make_pool(pool)
    if pool in ("P0xa", "P0xb"):
        return SB.make_pool(pool)
    if pool in ("P0y", "P0yb"):
        return B.BlendPool([B.pure_component(n) for n in BC.pool_members(pool)])
    raise ValueError(pool)


def synthetic_cap(pool: str):
    return B.SYNTHETIC_CAP if pool in ("P2", "P3") else None


def design_id(pool: str, spec: str) -> str:
    return f"{spec.replace(':', '_')}__{pool}_standin"


def provenance() -> dict:
    def git(*a):
        try:
            return subprocess.run(["git", *a], cwd=ROOT, capture_output=True, text=True, timeout=30).stdout.strip()
        except Exception:
            return ""
    dirty = git("status", "--porcelain", "--", "src", "configs", "data", "scripts/v5/base")
    return dict(git_head=git("rev-parse", "HEAD"), code_dirty=bool(dirty), freezing_model=CR.designated_freezing_model(),
                registry_overlay=os.environ.get(OVERLAY_ENV) != "0")


def _atomic_csv(df: pd.DataFrame, path: Path, **kw):
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    df.to_csv(tmp, index=False, **kw)
    os.replace(tmp, path)


# ------------------------------------------------------------------------------------------------ readings 4-7
def _dicyclo(smiles: str) -> bool:
    from rdkit import Chem
    m = Chem.MolFromSmiles(smiles)
    if m is None or any(a.GetIsAromatic() for a in m.GetAtoms()):
        return False
    if any(b.GetBondType() != Chem.BondType.SINGLE for b in m.GetBonds()):
        return False
    return m.GetRingInfo().NumRings() == 2


def readings(bp, models, u, res, spec, margins) -> dict:
    L = models["library"]
    pool = list(models["pool"])
    w = bp.expand(u)
    smi = [L.smiles[n] for n in pool]
    c = {x.key: x for x in res.constraints}
    r: dict = {}
    for key, tag, sign in (("flash_point_C", "flash", -1.0), ("freezing_point_C", "freezing", +1.0)):
        x = c.get(key)
        if x is None:
            continue
        lim = x.lo if sign < 0 else x.hi
        m = x.applied_margin
        r.update({f"{tag}_status": x.status, f"{tag}_model": x.model, f"{tag}_pred_C": x.prediction, f"{tag}_margin_C": m,
                  f"{tag}_bound_C": (x.prediction + sign * m) if (m is not None and np.isfinite(m)) else np.nan,
                  f"{tag}_limit_C": lim, f"{tag}_slack": x.slack})
    # ideal_svl variant: same constants (designated freezing tmap), same RDKit-MW mole fractions as the judge
    im, _, _, _ = E.read_applied_margin(margins.get(IDEAL_KEY))
    try:
        tmap = L._pr._tmap_freezing()
        iks = [L.ik(n) for n in pool]
        mw = np.array([L._ctx.mw(s) for s in smi])
        x = (w / w.sum()) / mw
        x = x / x.sum()
        v, _ = CR.freezing_C(x, iks, smi, tmap, model="ideal_svl")
        r.update(freezing_ideal_svl_C=v, freezing_ideal_svl_margin_C=im, freezing_ideal_svl_bound_C=v + im if im is not None else np.nan)
    except Exception as e:
        r.update(freezing_ideal_svl_C=np.nan, freezing_ideal_svl_note=f"unavailable: {type(e).__name__}: {e}")
    r["freezing_spec_limit_C"] = SPEC_LIMIT_FREEZE.get(spec, np.nan)
    # the judge's freezing PropertyModel keeps the name 'ideal_sle_mole_fraction' (library_route dataclasses.replace); its value is the
    # designated route, recorded here explicitly
    r["freezing_route"] = CR.designated_freezing_model()
    # Engler D86 (designated; the judge's numbers)
    for k, x in c.items():
        if k.startswith("distillation_"):
            r[f"engler_{k[len('distillation_'):]}_pred"] = x.prediction
            r[f"engler_{k[len('distillation_'):]}_slack"] = x.slack
            r[f"engler_{k[len('distillation_'):]}_status"] = x.status
    # Rayleigh D86 (check)
    rm = json.load(open(RAYLEIGH_MARGINS, encoding="utf-8"))["margins"]
    try:
        rho = np.array([B._rho15(L, n) for n in pool])
        if not np.all(np.isfinite(rho[w > 0])):
            raise E.ModelUnavailable("no grid rho at 288.15 K for a recipe member")
        q = B.rayleigh_quantities_smiles(smi, w, np.where(np.isfinite(rho), rho, 1.0))
        r.update(rayleigh_T10=q["T10"], rayleigh_T50=q["T50"], rayleigh_T90=q["T90"], rayleigh_FBP=q["FBP"],
                 rayleigh_T50_T10=q["T50_T10"], rayleigh_T90_T10=q["T90_T10"],
                 rayleigh_slack_T50_T10_C=q["T50_T10"] - 15.0 - float(rm["T50_T10"]),
                 rayleigh_slack_T90_T10_C=q["T90_T10"] - 40.0 - float(rm["T90_T10"]), rayleigh_status="ok")
    except Exception as e:
        r.update(rayleigh_status=f"unavailable: {type(e).__name__}: {e}")
    # dicycloparaffin mass fraction (reported only, no gate)
    r["dicycloparaffin_wt_frac"] = float(sum(wi for s, wi in zip(smi, w) if wi > 0 and _dicyclo(s)) / w.sum())
    # declared labels
    labels = []
    if res.passed:
        if r.get("rayleigh_status") == "ok" and (r["rayleigh_slack_T50_T10_C"] < 0 or r["rayleigh_slack_T90_T10_C"] < 0):
            labels.append("d86_model_dependent")
        if np.isfinite(r.get("freezing_ideal_svl_bound_C", np.nan)) and r["freezing_ideal_svl_bound_C"] > r["freezing_spec_limit_C"]:
            labels.append("freezing_route_dependent")
        if any(wi > 0 and n in ("bicyclohexyl", "decalin_cis") for n, wi in zip(pool, w)):
            labels.append("overlay_dependent")
    r["labels"] = ",".join(labels)
    return r


# ------------------------------------------------------------------------------------------------ one design
def run_one(pool: str, spec: str, out: Path, n_starts: int, seed: int, log) -> dict:
    did = design_id(pool, spec)
    prov = provenance()
    log(f"{did}: start {time.strftime('%Y-%m-%d %H:%M:%S')} n_starts={n_starts} seed={seed} {json.dumps(prov)}")
    t0 = time.time()
    base = dict(design_id=did, spec=spec, pool=pool, n_starts_declared=n_starts, seed_declared=seed, alt_standin=True, strict_constraints=True,
                synthetic_cap=synthetic_cap(pool), d86=True, **prov)
    stage = "build"
    try:
        bp = make_pool(pool)
        models = B.add_d86_models(B.build_models(bp))
        t1 = time.time()
        log(f"{did}: built {len(bp.comp_ids)} compounds in {t1 - t0:.1f} s")
        stage = "design"
        res, top, info = B.design_pool(spec, bp, models, W.MARGINS, n_starts=n_starts, seed=seed, d86=True, lam=B.LAMBDA_UNC,
                                       synthetic_cap=synthetic_cap(pool), alt_standin=True, strict_constraints=True)
    except E.ModelUnavailable as e:
        d = dict(base, status="model_unavailable", passed=False, stage=stage, message=str(e), build_s=round(time.time() - t0, 1))
        _atomic_csv(pd.DataFrame([d]), out / "designs" / pool / f"{did}__designs.csv")
        log(f"{did}: MODEL UNAVAILABLE at {stage}: {e}")
        return d
    d, cons = W._rows(res, top, f"{pool}_standin", spec)
    sl = np.array([x["slack"] for x in cons], float)
    mn = float(np.nanmin(sl)) if len(sl) else float("nan")
    u = np.array([top[n] for n in bp.names], float)
    margins = E.load_margins(W.MARGINS)
    rd = readings(bp, models, u, res, spec, margins)
    d = dict(base, **d, status="designed", stage="done", message="", top_weights_full=json.dumps(top),
             min_slack=mn, binding=",".join(x["constraint"] for x in cons if x["slack"] <= mn + 0.003),
             n_pass=sum(x["status"] == "pass" for x in cons), n_fail=sum(x["status"] == "fail" for x in cons),
             n_uncertified=sum(x["status"] == "uncertified" for x in cons),
             failing=",".join(x["constraint"] for x in cons if x["status"] == "fail"),
             uncertified=",".join(x["constraint"] for x in cons if x["status"] == "uncertified"),
             n_compounds=len(bp.comp_ids), build_s=round(t1 - t0, 1), design_s=round(time.time() - t1, 1), labels=rd["labels"])
    folder = out / "designs" / pool
    _atomic_csv(pd.DataFrame(cons), folder / f"{did}__constraints.csv")
    _atomic_csv(pd.DataFrame([dict(design_id=did, **rd)]), folder / f"{did}__readings.csv")
    _atomic_csv(pd.DataFrame([d]), folder / f"{did}__designs.csv")      # last: its presence marks the design as finished
    log(json.dumps({k: v for k, v in d.items() if k not in ("weights",)}, indent=1, default=str))
    for x in cons:
        log(f"  {x['constraint']:30s} {x['status']:16s} pred {x['prediction']:.5g} slack {x['slack']:+.5f} margin {x['applied_margin']}  "
            f"{str(x['reason'])[:150]}")
    log("readings: " + json.dumps(rd, default=str))
    return d


# ------------------------------------------------------------------------------------------------ tables
def compare_one(did_new: Path, pool: str, spec: str) -> list:
    folder, tag = COMMITTED[pool]
    old_id = f"{spec.replace(':', '_')}__{tag}"
    fo = folder / f"{old_id}__designs.csv"
    dn = pd.read_csv(did_new).iloc[0]
    rows = []
    if dn.get("status") != "designed" or not fo.exists():
        return [dict(design_id=dn.design_id, committed_id=old_id, field="status", new=dn.get("status"),
                     committed="present" if fo.exists() else "absent", delta=np.nan)]
    do = pd.read_csv(fo).iloc[0]
    cn = pd.read_csv(str(did_new).replace("__designs.csv", "__constraints.csv")).set_index("constraint")
    co = pd.read_csv(fo.with_name(f"{old_id}__constraints.csv")).set_index("constraint")
    mo = float(co.slack.astype(float).min())
    rows.append(dict(design_id=dn.design_id, committed_id=old_id, field="passed", new=bool(dn.passed), committed=bool(do.passed), delta=np.nan))
    rows.append(dict(design_id=dn.design_id, committed_id=old_id, field="min_slack", new=float(dn.min_slack), committed=mo,
                     delta=float(dn.min_slack) - mo))
    for k in sorted(set(cn.index) | set(co.index)):
        for f in ("status", "prediction", "slack", "applied_margin", "model"):
            a = cn.loc[k, f] if k in cn.index else np.nan
            b = co.loc[k, f] if k in co.index else np.nan
            num = f in ("prediction", "slack", "applied_margin")
            rows.append(dict(design_id=dn.design_id, committed_id=old_id, field=f"{k}.{f}", new=a, committed=b,
                             delta=(float(a) - float(b)) if num and pd.notna(a) and pd.notna(b) else np.nan))
    wn, wo = json.loads(dn.top_weights_full), json.loads(do.top_weights)
    for k in sorted(set(wn) | set(wo)):
        a, b = float(wn.get(k, 0.0)), float(wo.get(k, 0.0))
        rows.append(dict(design_id=dn.design_id, committed_id=old_id, field=f"weight.{k}", new=a, committed=b, delta=a - b))
    rows.append(dict(design_id=dn.design_id, committed_id=old_id, field="weight.L1", new=np.nan, committed=np.nan,
                     delta=float(sum(abs(float(wn.get(k, 0.0)) - float(wo.get(k, 0.0))) for k in set(wn) | set(wo)))))
    return rows


def rebuild(out: Path):
    files = sorted((out / "designs").glob("*/*__designs.csv"))
    if not files:
        return
    d = pd.concat([pd.read_csv(f) for f in files], ignore_index=True)
    rfiles = sorted((out / "designs").glob("*/*__readings.csv"))
    if rfiles:
        r = pd.concat([pd.read_csv(f) for f in rfiles], ignore_index=True)
        d = d.merge(r.drop(columns=[c for c in ("labels",) if c in r.columns]), on="design_id", how="left")
    order = {p: i for i, p in enumerate(POOLS)}
    d = d.sort_values(["spec", "pool"], key=lambda s: s.map(order) if s.name == "pool" else s)
    _atomic_csv(d, out / "designs.csv")
    comp = []
    for f in files:
        x = pd.read_csv(f).iloc[0]
        comp += compare_one(f, x.pool, x.spec)
    _atomic_csv(pd.DataFrame(comp), out / "compare.csv")


# ------------------------------------------------------------------------------------------------ timing (per optimiser evaluation)
def timing(out: Path, n_evals: int = 40):
    """One full property evaluation of the committed P0y-b Jet A recipe as the optimiser calls it (``_Sides.vector`` on the
    recipe-margin models plus the lambda penalty ``U``), under DESIGNATED_FREEZING_MODEL svl_unifac vs ideal_svl; models rebuilt for
    each.  ``n_evals`` distinct small perturbations (the D86 / viscosity / U caches key on weight bytes) after one warm-up."""
    from fuelblend.design import library_route as LR
    spec = "astm:jet_a"
    d = pd.read_csv(COMMITTED["P0yb"][0] / "astm_jet_a__P0yb_standin__designs.csv").iloc[0]
    tw = json.loads(d.top_weights)
    rows = []
    for model in ("svl_unifac", "ideal_svl"):
        CR.DESIGNATED_FREEZING_MODEL = model
        t0 = time.time()
        bp = make_pool("P0yb")
        models = B.add_d86_models(B.build_models(bp))
        build = time.time() - t0
        u0 = np.array([tw.get(n, 0.0) for n in bp.names], float) + 1e-4
        u0 /= u0.sum()
        sp = E.load_spec(spec, synthetic_blend=True, alt_limits=True, d86=True)
        margins = E.load_margins(W.MARGINS)
        use = LR.with_recipe_margins(models, margins)
        sides = E._Sides(sp, list(models["pool"]), use, margins, graded=True, alt_point_standin=True)
        U = B.alternate_spread_fn(models, sp)
        pool = list(models["pool"])
        rng = np.random.default_rng(0)

        def ev(u):
            v = sides.vector(bp.expand(u))
            return float(np.min(v)) - B.LAMBDA_UNC * float(U(pool, bp.expand(u)))
        ev(np.full(len(u0), 1.0 / len(u0)))                        # warm-up (gamma cache, lazy loads)
        us = [np.clip(u0 * (1 + 0.01 * rng.standard_normal(len(u0))), 1e-9, None) for _ in range(n_evals)]
        us = [x / x.sum() for x in us]
        t = time.perf_counter()
        vals = [ev(x) for x in us]
        per = (time.perf_counter() - t) / n_evals
        fz = float(models["props"]["freezing"].fn(pool, bp.expand(u0), None))
        rows.append(dict(freezing_model=model, build_s=round(build, 2), n_evals=n_evals, per_eval_ms=per * 1e3,
                         objective_at_recipe=vals[0], freezing_at_recipe_C=fz))
        print(rows[-1], flush=True)
    CR.DESIGNATED_FREEZING_MODEL = "svl_unifac"
    t = pd.DataFrame(rows)
    t["ratio_to_ideal"] = t.per_eval_ms / float(t.loc[t.freezing_model == "ideal_svl", "per_eval_ms"].iloc[0])
    _atomic_csv(t, out / "timing_per_eval.csv", float_format="%.4f")
    print(t.to_string())


# ------------------------------------------------------------------------------------------------ main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--spec", choices=SPECS)
    ap.add_argument("--pool", nargs="*", default=None, help="default: all eight pools")
    ap.add_argument("--n-starts", type=int, default=N_STARTS, help="smoke tests only; refused for the results folder")
    ap.add_argument("--seed", type=int, default=SEED, help="smoke tests only; refused for the results folder")
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT)
    ap.add_argument("--force", action="store_true", help="re-run pools already on disk for this spec")
    ap.add_argument("--rebuild", action="store_true", help="only rebuild designs.csv and compare.csv")
    ap.add_argument("--timing", action="store_true", help="per-evaluation timing svl_unifac vs ideal_svl (P0y-b Jet A recipe)")
    a = ap.parse_args()
    out = a.out.resolve()
    if a.timing:
        timing(out)
        return
    if a.rebuild:
        rebuild(out)
        return
    if a.spec is None:
        ap.error("--spec is required")
    if out == DEFAULT_OUT.resolve() and (a.n_starts != N_STARTS or a.seed != SEED):
        raise SystemExit(f"declared protocol is n_starts={N_STARTS}, seed={SEED}; overrides only with --out outside the results folder")
    pools = [ALIASES.get(p, p) for p in (a.pool or POOLS)]
    bad = [p for p in pools if p not in POOLS]
    if bad:
        ap.error(f"unknown pool(s) {bad}; one of {POOLS} (aliases {ALIASES})")
    out.mkdir(parents=True, exist_ok=True)
    for pool in pools:
        did = design_id(pool, a.spec)
        if (out / "designs" / pool / f"{did}__designs.csv").exists() and not a.force:
            print(f"{did}: on disk, skipped (--force to re-run)", flush=True)
            continue
        with open(out / f"log_{pool}_{a.spec.replace(':', '_')}.txt", "a", encoding="utf-8") as fh:
            def log(msg, _fh=fh):
                print(msg, flush=True)
                _fh.write(msg + "\n")
                _fh.flush()
            try:
                run_one(pool, a.spec, out, a.n_starts, a.seed, log)
            except Exception:                          # anything else: recorded in the log, the run continues
                log(f"{did}: ERROR\n{traceback.format_exc()}")
        rebuild(out)


if __name__ == "__main__":
    main()
