"""v4 Task 13: solid-liquid equilibrium (liquidus) temperatures from ThermoML, ideal Schroeder-van Laar vs library Coutinho.

Data      : data/processed/v4/thermoml_sle_hc.csv (extras/fuelblend_data/thermoml_extract.py --sle; all-hydrocarbon binaries/ternaries,
            LIQUID composition, ePresentation 'Direct value, X' only, P <= 150 kPa).
Constants : data/processed/v4/pure_props_v4.csv, kind == 'pure' AND status == 'verified' only (R16: no disputed / unverified /
            reference_only / not_comparable values, no automatic substitution).  Tm in C, dHfus in kJ/mol.
Models
  ideal_svl      : ln x_i = -(dHfus_i/R)(1/T - 1/Tm_i) per component (mole fractions, no dCp, no solid solution),
                   liquidus T = max_i T_i.  Needs verified Tm AND dHfus for EVERY component of the system, else the row is not predicted.
  coutinho_lib   : src/fuelblend/thermo/coutinho_wax.predict_freezing_point_coutinho, library code unmodified, called as the
                   repo calls it (mass fractions, library's own Tm/dHfus table, eutectic depression subtracted).
  coutinho_lib_x : same function fed MOLE fractions (diagnostic: the library puts the fractions it is given inside the SvL log).
  *_offset       : model + a constant fitted on the TRAINING folds only (mean of obs - pred).
CV        : GroupKFold(5) by system_id (md5[:12] of sorted InChIKeys); offsets are refit inside each training fold.
Coutinho predictions are kept only where the library has Tm/dHfus for the crystallising component named in the ThermoML record
(it silently skips components without constants, which would otherwise return another solute's liquidus).
Scored rows exclude pure-component end points (x = 0 or 1), which are melting points, not mixture liquidus points.
Outputs   : results/v4/SLE/sle_eval.csv (one row per SLE point), sle_summary.csv, sle_margin.csv
"""
from __future__ import annotations

import hashlib
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem.inchi import InchiToInchiKey
from sklearn.model_selection import GroupKFold

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from fuelblend.thermo.coutinho_wax import R_GAS, get_pure_thermal_properties, predict_freezing_point_coutinho  # noqa: E402

SLE = ROOT / "data" / "processed" / "v4" / "thermoml_sle_hc.csv"
PURE = ROOT / "data" / "processed" / "v4" / "pure_props_v4.csv"
OUT = ROOT / "results" / "v4" / "SLE"
P_MAX_KPA = 150.0
N_SPLITS = 5


def _hc_class():
    spec = importlib.util.spec_from_file_location("fb_common", ROOT / "extras" / "fuelblend_data" / "common.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m.hc_class, m.carbon_number


HC_CLASS, CARBON_NUMBER = _hc_class()


def verified_constants() -> dict[str, dict[str, float]]:
    """inchikey -> {tm_K, dh_J_mol} from verified, pure-kind rows of the v4 audited table."""
    p = pd.read_csv(PURE)
    p = p[(p.kind == "pure") & (p.status == "verified") & p.prop.isin(["Tm", "dHfus"])]
    out: dict[str, dict[str, float]] = {}
    for _, r in p.iterrows():
        m = Chem.MolFromSmiles(r.smiles)
        if m is None:
            continue
        ik = Chem.inchi.MolToInchiKey(m)[:14]            # connectivity block: ignores stereo (decalin cis/trans stays distinct by smiles only if flat)
        d = out.setdefault(ik, {})
        if r.prop == "Tm":
            d["tm_K"] = float(r.value) + 273.15
        else:
            d["dh_J_mol"] = float(r.value) * 1000.0
    return {k: v for k, v in out.items() if "tm_K" in v and "dh_J_mol" in v}


def ideal_svl(x, iks, const) -> float:
    ts = []
    for xi, ik in zip(x, iks):
        c = const.get(ik[:14])
        if c is None:
            return np.nan                                  # a component without verified constants: no prediction, no substitution
        if xi <= 0:
            continue
        inv = 1.0 / c["tm_K"] - (R_GAS / c["dh_J_mol"]) * np.log(xi)
        if inv <= 0:
            return np.nan
        ts.append(1.0 / inv)
    return max(ts) if ts else np.nan


def build_rows() -> pd.DataFrame:
    d = pd.read_csv(SLE)
    d = d[d.presentation == "Direct value, X"]
    d = d[d.P_kPa.isna() | (d.P_kPa <= P_MAX_KPA)].reset_index(drop=True)
    const = verified_constants()
    recs, comp_rows = [], {}
    for _, r in d.iterrows():
        inchis = json.loads(r.inchis)
        names = json.loads(r.names)
        vals = np.array(json.loads(r.comp_values), float)
        mols = [Chem.MolFromInchi(i) for i in inchis]
        mw = np.array([Descriptors.MolWt(m) for m in mols])
        iks = [InchiToInchiKey(i) for i in inchis]
        if r.comp_type == "mass fraction":
            w = vals / vals.sum()
            x = (w / mw) / np.sum(w / mw)
        else:
            x = vals / vals.sum()
            w = x * mw / np.sum(x * mw)
        cls = [HC_CLASS(m) for m in mols]
        for ik, m, c in zip(iks, mols, cls):
            comp_rows.setdefault(ik, dict(inchikey=ik, smiles=Chem.MolToSmiles(m), n_carbon=float(CARBON_NUMBER(m)), cls=c))
        endpoint = bool(np.any(x < 1e-6) or np.any(x > 1 - 1e-6))
        cry = [InchiToInchiKey(i) for i in json.loads(r.crystal_inchis) if i]
        recs.append(dict(
            file=r.file, names="+".join(names), n_comp=len(names), comp_type=r.comp_type,
            system_id=hashlib.md5("|".join(sorted(iks)).encode()).hexdigest()[:12],
            system_class="+".join(sorted(cls)), inchikeys=json.dumps(iks), x_mole=json.dumps([float(v) for v in np.round(x, 6)]),
            w_mass=json.dumps([float(v) for v in np.round(w, 6)]), crystal=";".join(str(c) for c in json.loads(r.crystal_names)),
            crystal_ik=json.dumps(cry), P_kPa=r.P_kPa, T_obs_K=float(r.value), endpoint=endpoint,
            ideal_svl_K=ideal_svl(x, iks, const)))
    df = pd.DataFrame(recs)
    # library Coutinho (unmodified); needs a comp_df keyed by inchikey exactly as get_pure_thermal_properties expects
    comp_df = pd.DataFrame(comp_rows.values()).set_index("inchikey")
    tmap = get_pure_thermal_properties(comp_df)
    # The library silently skips components it has no Tm/dHfus for.  If the component that actually crystallises (ThermoML
    # PhaseID) is one of them, the library answer is the liquidus of a DIFFERENT solute, so such rows are not eligible.
    def lib_covers(ik):
        i = tmap.get(ik)
        return bool(i) and np.isfinite(i["tm_K"]) and np.isfinite(i["dh_fus_J_mol"])
    df["coutinho_crystallizer_covered"] = [bool(json.loads(c)) and all(lib_covers(k) for k in json.loads(c)) for c in df.crystal_ik]
    for col, key in (("coutinho_lib_K", "w_mass"), ("coutinho_lib_x_K", "x_mole")):
        out = []
        for _, r in df.iterrows():
            t, _, meta = predict_freezing_point_coutinho(json.loads(r[key]), json.loads(r.inchikeys), tmap)
            out.append(np.nan if (t == -999.0 or not np.isfinite(t) or not r.coutinho_crystallizer_covered) else t + 273.15)
        df[col] = out
    df["scored"] = ~df.endpoint
    return df


def add_cv(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    models = ["ideal_svl", "coutinho_lib", "coutinho_lib_x"]
    for m in models:
        df[f"{m}_offset_K"] = np.nan
    df["fold"] = -1
    sc = df.index[df.scored]
    groups = df.loc[sc, "system_id"]
    gkf = GroupKFold(n_splits=min(N_SPLITS, groups.nunique()))
    for k, (tr, te) in enumerate(gkf.split(sc, groups=groups)):
        tri, tei = sc[tr], sc[te]
        df.loc[tei, "fold"] = k
        for m in models:
            ok = df.loc[tri, f"{m}_K"].notna()
            if ok.sum() == 0:
                continue
            off = float((df.loc[tri, "T_obs_K"] - df.loc[tri, f"{m}_K"])[ok].mean())    # training folds only
            df.loc[tei, f"{m}_offset_K"] = df.loc[tei, f"{m}_K"] + off
    return df


def metrics(err: pd.Series) -> dict:
    e = err.dropna()
    if not len(e):
        return dict(n=0, mae=np.nan, bias=np.nan, rmse=np.nan, p90_abs=np.nan)
    return dict(n=len(e), mae=float(e.abs().mean()), bias=float(e.mean()), rmse=float(np.sqrt((e ** 2).mean())),
                p90_abs=float(np.quantile(e.abs(), 0.9)))


def summarise(df: pd.DataFrame) -> pd.DataFrame:
    s = df[df.scored]
    both = s.ideal_svl_K.notna() & s.coutinho_lib_K.notna()
    rows = []
    for subset, mask in (("all", slice(None)), ("common (ideal & coutinho both defined)", both)):
        sub = s if subset == "all" else s[mask]
        for m in ("ideal_svl", "ideal_svl_offset", "coutinho_lib", "coutinho_lib_offset", "coutinho_lib_x", "coutinho_lib_x_offset"):
            col = f"{m}_K" if not m.endswith("offset") else f"{m}_K"
            err = sub[col] - sub.T_obs_K                     # pred - obs, K (positive = predicts warmer liquidus than observed)
            for cls, g in [("ALL", sub)] + list(sub.groupby("system_class")):
                e = err.loc[g.index]
                rows.append(dict(subset=subset, model=m, system_class=cls, n_systems=g.loc[e.notna(), "system_id"].nunique(),
                                 **metrics(e)))
    out = pd.DataFrame(rows)
    return out[out.n > 0].reset_index(drop=True)


def margin_table(df: pd.DataFrame) -> pd.DataFrame:
    """Unsafe-side residual (observed liquidus - predicted, K; positive = the real mixture freezes WARMER than predicted)
    for the offset-free models on the rows where both models are defined.  A margin proposal is read off the point-level and
    the system-level (mean per system) quantiles; no parameter is fitted here, so the residuals are out-of-fold by construction."""
    s = df[df.scored & df.ideal_svl_K.notna() & df.coutinho_lib_K.notna()]
    rows = []
    for m in ("ideal_svl", "coutinho_lib"):
        r = s.T_obs_K - s[f"{m}_K"]
        per_sys = r.groupby(s.system_id).mean()
        rows.append(dict(model=m, n_points=len(r), n_systems=per_sys.size,
                         point_q50=r.quantile(0.5), point_q90=r.quantile(0.9), point_q95=r.quantile(0.95), point_max=r.max(),
                         system_mean_max=per_sys.max(), system_mean_q90=per_sys.quantile(0.9)))
    return pd.DataFrame(rows)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    df = add_cv(build_rows())
    df.to_csv(OUT / "sle_eval.csv", index=False)
    summ = summarise(df)
    summ.to_csv(OUT / "sle_summary.csv", index=False)
    mg = margin_table(df)
    mg.to_csv(OUT / "sle_margin.csv", index=False)
    s = df[df.scored]
    print(f"SLE rows (direct, P<={P_MAX_KPA:g}): {len(df)}  scored (no pure end points): {len(s)}  systems: {s.system_id.nunique()}")
    print(f"ideal_svl defined: {s.ideal_svl_K.notna().sum()} pts / {s[s.ideal_svl_K.notna()].system_id.nunique()} systems;"
          f" coutinho defined: {s.coutinho_lib_K.notna().sum()} pts / {s[s.coutinho_lib_K.notna()].system_id.nunique()} systems")
    print(summ[summ.system_class == "ALL"].to_string(index=False))
    print(mg.round(2).to_string(index=False))


if __name__ == "__main__":
    main()
