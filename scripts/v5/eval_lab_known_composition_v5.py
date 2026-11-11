"""v5 S1: cold-bin lab viscosity margins with the Boehm 2022 table added (declaration: results/v5/S1/README.md).

Reuses the v4 logic (scripts/v4/eval_lab_known_composition.py, build_viscosity_leaderboard_v4.refutas_nu,
build_uncertainty_v4.add); writes ONLY to results/v5/S1/ (v4 inputs and outputs are read, never written).

Methods per (system, T) record, residual = ln(pred / measured) of DYNAMIC viscosity:
  refutas    v4 raw Refutas (mass-fraction VBN blending, ideal-volume density)
  arrhenius  ln(eta_mix) = sum_i x_i ln(eta_i)  (mole fractions, pure eta_i(T) from the PurePropertyProvider)
  walther    Refutas nu at warm reference T (288.15..318.15 K) -> ASTM D341 line (walther_viscosity.py) -> extrapolated
             to T (only for T < 283.15 K), times ideal-volume density
Each also as <method>_biascorr: a per-bin constant offset fit on DEV systems only (dev rows: leave-one-system-out).

Split: every v4 system keeps its v4 role (splits/v4/lab_systems_splits.json); Boehm systems not in that file are assigned by
md5 hash BEFORE any residual is computed (assign_new_roles).  The split file refuses to be overwritten by a differing split.

Run:  .venv/Scripts/python.exe scripts/v5/eval_lab_known_composition_v5.py
"""
from __future__ import annotations

import hashlib
import importlib.util
import json
import math
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

OUT = ROOT / "results" / "v5" / "S1"
V4_SPLIT = ROOT / "splits" / "v4" / "lab_systems_splits.json"
BOEHM = ROOT / "data" / "real_fuels" / "v5" / "boehm2022_viscosity.csv"
MEAS = ROOT / "data" / "unified" / "measurements_dedup.csv"
V4_YAML = ROOT / "configs" / "uncertainty_v4.yaml"
SALT = "v5-S1-boehm2022"
LOCK_FRAC = 0.20
BINS = ("T_lt_263K", "T_263_283K", "T_ge_283K")
COLD_BINS = BINS[:2]
WALTHER_TREF_K = (288.15, 298.15, 308.15, 318.15)
WALTHER_MAX_T_K = 283.15                      # Walther is extrapolation; evaluated for the two cold bins only
METHODS = ("refutas", "arrhenius", "walther")
MIN_DEV_SYSTEMS = 9
LAB_VISC_KEY = "viscosity_kin_cSt.lab_mixtures"          # v4 key stem; refutas keeps '.raw_refutas.<bin>'
KEY_NAME = {"refutas": "raw_refutas", "arrhenius": "arrhenius", "walther": "walther_extrap"}


def _load(name: str, rel: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / rel)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


def t_bin(T_K: float) -> str:
    """Same bins as v4 (263.15 and 283.15 K belong to the warmer bin)."""
    return "T_lt_263K" if T_K < 263.15 else ("T_263_283K" if T_K < 283.15 else "T_ge_283K")


# ------------------------------------------------------------------------------------------------- split
def hash_role(system_id: str, salt: str = SALT, frac: float = LOCK_FRAC) -> str:
    """Deterministic role of a new system: locked iff md5(salt|system_id) / 2**128 < frac."""
    h = int(hashlib.md5(f"{salt}|{system_id}".encode()).hexdigest(), 16) / 2 ** 128
    return "locked" if h < frac else "dev"


def assign_new_roles(boehm_ids, v4_split: dict, salt: str = SALT, frac: float = LOCK_FRAC) -> dict:
    """-> {'dev': [...], 'locked': [...], 'boehm_new': {'dev': [...], 'locked': [...]}, 'boehm_in_v4': {...}}.
    Systems already in the v4 split keep their v4 role; the others get hash_role."""
    v4_dev, v4_lock = set(v4_split["dev"]), set(v4_split["locked"])
    assert not (v4_dev & v4_lock)
    new = {"dev": [], "locked": []}
    old = {"dev": [], "locked": []}
    for s in sorted(set(boehm_ids)):
        if s in v4_dev:
            old["dev"].append(s)
        elif s in v4_lock:
            old["locked"].append(s)
        else:
            new[hash_role(s, salt, frac)].append(s)
    return dict(dev=sorted(v4_dev | set(new["dev"])), locked=sorted(v4_lock | set(new["locked"])),
                boehm_new=new, boehm_in_v4=old)


def write_split(path: Path, sp: dict) -> dict:
    if path.exists():
        old = json.loads(path.read_text())
        if old != json.loads(json.dumps(sp)):
            raise SystemExit(f"{path} exists and differs from the recomputed split; refusing to overwrite")
        return old
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(sp, indent=1) + "\n")
    return sp


# --------------------------------------------------------------------------------------------- statistics
def loo_system_offsets(resid, systems):
    """Per-system leave-one-system-out offset: mean over the OTHER systems of each system's mean residual.
    -> (offset per row, all-system offset).  Needs >= 2 systems, else NaN."""
    resid, systems = np.asarray(resid, float), np.asarray(systems)
    if resid.size == 0:
        return np.array([]), float("nan")
    u, inv = np.unique(systems, return_inverse=True)
    b = np.array([resid[inv == k].mean() for k in range(len(u))])
    full = float(b.mean())
    if len(u) < 2:
        return np.full(resid.shape, np.nan), full
    off = (b.sum() - b) / (len(u) - 1)
    return off[inv], full


def walther_extrapolate(T_ref_K, nu_ref_cSt, T_K) -> float:
    from fuelblend.thermo.walther_viscosity import fit_walther_parameters, walther_viscosity
    A, B = fit_walther_parameters(list(T_ref_K), list(nu_ref_cSt))
    return walther_viscosity(float(T_K), A, B)


def arrhenius_ln_eta(x, ln_eta_i) -> float:
    """ln(eta_mix) = sum x_i ln(eta_i); x re-normalised."""
    x = np.asarray(x, float)
    return float((x / x.sum() * np.asarray(ln_eta_i, float)).sum())


# ----------------------------------------------------------------------------------------------- scoring
def lab_selection(meas: pd.DataFrame) -> pd.DataFrame:
    """Same rows as v4 build_lab without the temperature filter."""
    return meas[meas.kind.isin(["binary", "multi"]) & meas.is_hc.astype(bool) & (meas.property == "viscosity_dyn_mPas")
                & ~meas.conflict.astype(bool)]


def score_rows(meas: pd.DataFrame, reg, provider, cls, nc, lb, classes, source_col=None) -> pd.DataFrame:
    """-> one row per (record, method) with predicted/measured dynamic viscosity (mPa s) and ln residual."""
    w = lab_selection(meas)
    rows = []
    for r in w.itertuples():
        ks = json.loads(r.comp_inchikeys)
        ws = np.array(json.loads(r.w_mass), float)
        xs = np.array(json.loads(r.x_mole), float)
        ii = [reg.ik2i.get(k) for k in ks]
        if any(i is None or not provider.valid[i] or not reg.usable[i] or cls[i] not in classes or not np.isfinite(nc[i])
               for i in ii) or not np.all(np.isfinite(ws)) or ws.sum() <= 0 or not np.all(np.isfinite(xs)):
            continue
        ii = np.array(ii)
        ws = ws / ws.sum()
        T = float(r.T_K)
        nu, rho_mix, ln_eta, _ = lb.refutas_nu(provider, ii, ws, T)
        pred = {"refutas": nu * rho_mix / 1000.0, "arrhenius": math.exp(arrhenius_ln_eta(xs, ln_eta))}
        if T < WALTHER_MAX_T_K:
            nu_ref = [lb.refutas_nu(provider, ii, ws, Tr)[0] for Tr in WALTHER_TREF_K]
            pred["walther"] = walther_extrapolate(WALTHER_TREF_K, nu_ref, T) * rho_mix / 1000.0
        for m, p in pred.items():
            if np.isfinite(p) and p > 0:
                rows.append(dict(system_id=r.system_id, source=getattr(r, "sources", ""), T_K=T, t_bin=t_bin(T), n_comp=len(ks),
                                 method=m, predicted=p, measured=float(r.value), resid=math.log(p) - math.log(float(r.value))))
    return pd.DataFrame(rows)


def add_bias_corrected(res: pd.DataFrame) -> pd.DataFrame:
    """Append <method>_biascorr rows.  Offset per (method, bin) from DEV systems only; dev rows use leave-one-system-out
    offsets, locked rows the all-dev offset.  Bins with < 2 dev systems get no corrected rows."""
    out = []
    for (m, b), g in res.groupby(["method", "t_bin"]):
        dev = g[g.split_role == "dev"]
        if dev.system_id.nunique() < 2:
            continue
        off_dev, full = loo_system_offsets(dev.resid.to_numpy(), dev.system_id.to_numpy())
        c_dev = dev.assign(method=f"{m}_biascorr", resid_raw=dev.resid, resid=dev.resid - off_dev, offset=off_dev)
        lk = g[g.split_role == "locked"]
        c_lk = lk.assign(method=f"{m}_biascorr", resid_raw=lk.resid, resid=lk.resid - full, offset=full)
        out += [c_dev, c_lk]
    base = res.assign(resid_raw=res.resid, offset=0.0)
    return pd.concat([base] + out, ignore_index=True) if out else base


def summarize(res: pd.DataFrame) -> pd.DataFrame:
    recs = []
    for (m, b, role), g in res.groupby(["method", "t_bin", "split_role"]):
        recs.append(dict(method=m, t_bin=b, split_role=role, n_systems=g.system_id.nunique(), n_records=len(g),
                         bias_ln=g.resid.mean(), rmse_ln=float(np.sqrt((g.resid ** 2).mean())),
                         n_boehm_systems=g[g.source == "boehm2022"].system_id.nunique()))
    return pd.DataFrame(recs)


def method_key(method: str, b: str) -> str:
    corr = method.endswith("_biascorr")
    base = method[: -len("_biascorr")] if corr else method
    return f"{LAB_VISC_KEY}.{KEY_NAME[base]}{'_biascorr' if corr else ''}.{b}"


# ---------------------------------------------------------------------------------------------------- main
def build_yaml(res: pd.DataFrame, bu, fv_rel: str) -> dict:
    """Run the existing v4 `add` for each (method, bin); then copy the other v4 entries unchanged."""
    import yaml
    bu.ENTRIES.clear()
    bu.LOCKED_ROWS.clear()
    descr = {"refutas": "raw Refutas (kinematic, mass-fraction VBN blending) with ideal-volume density, pure data from the PurePropertyProvider",
             "arrhenius": "Arrhenius log-mixing rule ln eta = sum x_i ln eta_i (mole fractions), pure eta_i(T) from the PurePropertyProvider",
             "walther": "Refutas at 288-318 K extrapolated along the ASTM D341 (Walther) line to T, ideal-volume density"}
    for method in sorted(res.method.unique(), key=lambda s: (s.endswith("_biascorr"), METHODS.index(s.split("_")[0]))):
        base = method.replace("_biascorr", "")
        for b in BINS:
            g = res[(res.method == method) & (res.t_bin == b)]
            vd, vl = g[g.split_role == "dev"], g[g.split_role == "locked"]
            if len(vd) == 0:
                continue
            key = method_key(method, b)
            if method == "refutas" and b not in COLD_BINS:
                continue                              # warm bin: the v4 entry is kept verbatim (S1 scope = cold bins)
            bu.add(key, vd.resid, vd.system_id, kind="two", direction="over", unit="ln-ratio",
                   model=descr[base] + (" + constant ln offset fit on dev systems (dev rows leave-one-system-out)" if method.endswith("_biascorr") else ""),
                   residual_file=fv_rel,
                   scope=(f"lab hydrocarbon binary/multi dynamic-viscosity mixtures with measured composition (v4 lab table + Boehm 2022), "
                          f"temperature bin {b} (edges belong to the warmer bin); dev systems only"),
                   group_unit="lab system (distinct InChIKey set)", bias=vd.resid.mean(), locked=(vl.resid, vl.system_id), log_space=True,
                   caveats=["Known compositions; Boehm rows: volume fractions converted with 15 C densities, nested pre-blends assumed ideally additive.",
                            "Pure eta_i(T), rho_i(T) at cold T come from the provider (fits of warm data, extrapolated).",
                            "Binaries/ternaries of ~18 hydrocarbons; transfer to 100-component jet fuels is not demonstrated.",
                            "Two-sided ln-ratio margin; inf when fewer than 9 independent dev systems exist (no substitute number)."])
            ENTRY = bu.ENTRIES[key]
            ENTRY["recompute"] = dict(column="resid", group="system_id",
                                      where=[["method", "eq", method], ["split_role", "ne", "locked"], ["t_bin", "eq", b]])
            ENTRY["residual_derivation"] = f"residual = ln(pred) - ln(measured) of {method}, column resid, dev systems in t_bin {b}; group = system_id"
    v4 = yaml.safe_load(V4_YAML.read_text(encoding="utf-8"))
    old = {f"{LAB_VISC_KEY}.raw_refutas.{b}" for b in COLD_BINS}
    margins = {k: v for k, v in v4["margins"].items() if k not in old}
    margins.update(bu.ENTRIES)
    doc = {k: v for k, v in v4.items() if k not in ("margins", "version", "generated_by")}
    doc = dict(version="5.0-S1", generated_by="scripts/v5/eval_lab_known_composition_v5.py", **doc)
    doc["margins"] = margins
    doc["v5_note"] = ("Derived from configs/uncertainty_v4.yaml: the two COLD lab raw_refutas viscosity entries (T_lt_263K, T_263_283K) are REPLACED by entries recomputed on the v4 lab "
                      "table + Boehm 2022; the warm T_ge_283K raw_refutas entry is kept verbatim from v4 (same conformal rule); extra keys .arrhenius/.walther_extrap/.*_biascorr are comparators. All other entries are copied unchanged.")
    return doc


def main():
    import yaml
    OUT.mkdir(parents=True, exist_ok=True)
    meas4 = pd.read_csv(MEAS, low_memory=False)
    boehm = pd.read_csv(BOEHM, low_memory=False)
    assert list(boehm.columns) == list(meas4.columns), "Boehm table is not in the measurements_dedup schema"
    boehm["sources"] = boehm["sources"].fillna("boehm2022")
    meas = pd.concat([meas4, boehm], ignore_index=True)

    # 1. split FIRST (ids only; no residual exists yet)
    v4_split = json.loads(V4_SPLIT.read_text())
    b_ids = lab_selection(boehm).system_id.unique()
    sp = assign_new_roles(b_ids, v4_split)
    sp_doc = dict(salt=SALT, locked_fraction=LOCK_FRAC, rule="locked iff int(md5(salt|system_id),16)/2**128 < locked_fraction; v4 systems keep v4 roles",
                  v4_split_file="splits/v4/lab_systems_splits.json", n_boehm_systems=len(b_ids), **sp)
    sp_doc = write_split(OUT / "lab_systems_splits_v5.json", sp_doc)
    role = {**{s: "dev" for s in sp_doc["dev"]}, **{s: "locked" for s in sp_doc["locked"]}}
    print(f"split: {len(sp_doc['dev'])} dev / {len(sp_doc['locked'])} locked; Boehm systems {len(b_ids)} "
          f"(new dev {len(sp_doc['boehm_new']['dev'])}, new locked {len(sp_doc['boehm_new']['locked'])}, already in v4 split "
          f"{sum(len(v) for v in sp_doc['boehm_in_v4'].values())})", flush=True)

    # 2. residuals
    lb = _load("build_viscosity_leaderboard_v4", "scripts/v4/build_viscosity_leaderboard_v4.py")
    from fuelblend.mixture.cold_flow import compound_meta
    from fuelblend.mixture.common import CLASSES, CompoundRegistry, PurePropertyProvider
    reg = CompoundRegistry(ROOT / "data" / "unified" / "compounds.csv")
    _, cls, nc = compound_meta(ROOT / "data" / "unified" / "compounds.csv")
    provider = PurePropertyProvider(reg, ROOT / "data" / "cache" / "pure_grid.npz", verbose=False)
    res = score_rows(meas, reg, provider, cls, nc, lb, CLASSES)
    missing = set(res.system_id) - set(role)
    assert not missing, f"{len(missing)} scored systems without a split role"
    res["split_role"] = res.system_id.map(role)

    # regression: v4 refutas rows reproduce results/v4/L/viscosity_lab_residuals.csv
    old = pd.read_csv(ROOT / "results" / "v4" / "L" / "viscosity_lab_residuals.csv")
    v4only = res[(res.method == "refutas") & (res.source != "boehm2022")]
    # same selection order as build_lab => compare positionally (keys system/T/value are not unique across compositions)
    v4only = v4only.reset_index(drop=True)
    n_match = min(len(old), len(v4only))
    assert len(old) == len(v4only), f"row count differs: v4 {len(old)} vs v5 refutas non-Boehm {len(v4only)}"
    assert (old.system_id.to_numpy() == v4only.system_id.to_numpy()).all()
    reg_diff = float(np.abs(old.resid.to_numpy() - v4only.resid.to_numpy()).max())
    common, chk = range(n_match), range(len(old))
    print(f"v4 regression: {len(common)}/{len(chk)} v4 refutas residuals matched, max |diff| = {reg_diff:.2e}", flush=True)
    assert reg_diff < 1e-6, "v5 refutas does not reproduce the v4 lab residuals"

    res = add_bias_corrected(res)
    fv = OUT / "viscosity_lab_residuals_v5.csv"
    res.to_csv(fv, index=False, float_format="%.8g")
    summ = summarize(res)
    summ.to_csv(OUT / "summary_by_bin_method.csv", index=False, float_format="%.5g")

    # 3. generated yaml via the existing v4 `add`
    bu = _load("build_uncertainty_v4", "scripts/v4/build_uncertainty_v4.py")
    doc = build_yaml(res, bu, str(fv.relative_to(ROOT)).replace("\\", "/"))
    text = ("# GENERATED - do not edit. Produced by scripts/v5/eval_lab_known_composition_v5.py (v5 S1) via scripts/v4/build_uncertainty_v4.py::add.\n"
            "# Same residual convention and conformal rule as configs/uncertainty_v4.yaml, which is NOT modified. 'inf' = no valid margin; read applied_margin only.\n"
            + yaml.safe_dump(doc, sort_keys=False, default_flow_style=False, width=140))
    (OUT / "uncertainty_v5.yaml").write_text(text, encoding="utf-8")

    # 4. results report
    L = ["# v5 S1 results (generated by scripts/v5/eval_lab_known_composition_v5.py)", "",
         f"Reading declared in README.md before scoring. Split: {len(sp_doc['dev'])} dev / {len(sp_doc['locked'])} locked systems "
         f"(Boehm: {len(b_ids)} systems; new dev {len(sp_doc['boehm_new']['dev'])}, new locked {len(sp_doc['boehm_new']['locked'])}, already in v4 split "
         f"{sum(len(v) for v in sp_doc['boehm_in_v4'].values())}). v4 regression on refutas: {len(common)}/{len(chk)} rows matched, max |diff| {reg_diff:.1e}.", "",
         "## Bias and RMSE (ln(pred/meas), dynamic viscosity) per bin x method x role", "",
         "| bin | method | role | systems | Boehm systems | records | bias | RMSE |", "|---|---|---|---|---|---|---|---|"]
    for b in BINS:
        for m in sorted(summ.method.unique(), key=lambda s: (s.endswith("_biascorr"), s)):
            for role_ in ("dev", "locked"):
                s = summ[(summ.t_bin == b) & (summ.method == m) & (summ.split_role == role_)]
                if len(s):
                    s = s.iloc[0]
                    L.append(f"| {b} | {m} | {role_} | {s.n_systems} | {s.n_boehm_systems} | {s.n_records} | {s.bias_ln:+.4f} | {s.rmse_ln:.4f} |")
    L += ["", "## Conformal margins (dev systems, alpha 0.1, two-sided |ln|)", "",
          "| key | dev systems | k | margin (value) | applied_margin | locked covered |", "|---|---|---|---|---|---|"]
    for k, e in doc["margins"].items():
        if k.startswith(LAB_VISC_KEY):
            lc = e.get("locked_check")
            lcs = f"{lc['groups_covered']}/{lc['n_locked_groups']}" if lc else "n/a"
            L.append(f"| {k} | {e['n_groups']} | {e['conformal_rank_k']} | {e['value'] if isinstance(e['value'], str) else round(e['value'], 4)} | "
                     f"{e.get('applied_margin', 'none (inf)')} | {lcs} |")
    (OUT / "RESULTS.md").write_text("\n".join(L) + "\n", encoding="utf-8")
    print("\n".join(L), flush=True)


if __name__ == "__main__":
    main()
