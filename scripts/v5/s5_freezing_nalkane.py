"""v5 S5: freezing point from the measured n-alkane distribution (first-species-to-freeze) versus a constant and ideal SvL.

Models, readings and data are declared in results/v5/S5/README.md BEFORE this script was first run.  Nothing here edits the
README; the Results section is filled by hand from the CSVs written below.

Models (all leave-one-group-out, LOGO; the scored group's freezing value never enters its own prediction):
  C   mean of training freezing points
  I   ideal Schroeder-van Laar on the full mapped composition (G1 code: ideal_sle_C and the coutinho Tm/dHfus table)
  F0  no fit: T_i from ln x_i = -(dHfus_i/R)(1/T - 1/Tm_i) per n-alkane with audited constants; prediction = max_i T_i
  F1  a + b * T_first(F0), a and b fitted on training groups
  P   ridge on [n-alkane mass share C9-C12, C13-C16, C17+, total n-alkane share, aromatic share] (alpha by inner leave-one-out)
Sensitivity (not part of the declared readings): F0x / F1x use the repo's complete Coutinho n-alkane table (C6-C36) instead of
the audited table, which has no C8, C9, C10, C17, C19.

Outputs: results/v5/S5/{s5_fuels.csv, s5_predictions.csv, s5_metrics.csv, s5_bootstrap.csv, s5_nalkane_constants.csv,
                        s5_excluded.csv, s5_readings.csv}
"""
from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem.inchi import InchiToInchiKey
from scipy.stats import spearmanr
from sklearn.linear_model import RidgeCV

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
V5 = ROOT / "data" / "real_fuels" / "v5"
S2 = ROOT / "results" / "v5" / "S2"
S3 = ROOT / "results" / "v5" / "S3"
G1D = ROOT / "results" / "v5" / "G1"
OUT = ROOT / "results" / "v5" / "S5"
R_GAS = 8.314462618
MAX_UNMAPPED = 0.10
MIN_GROUPS = 20
N_BOOT = 5000
SEED = 0
PRIORITY = {"AV-33-22": 0, "FuelLib": 1, "AJFTD": 2}          # same tie-break as S3 select_one_per_group
ALPHAS = np.logspace(-3, 3, 13)
AROM = {"alkylbenzene", "cycloaromatic", "diaromatic", "aromatic"}


# ---------------------------------------------------------------------------------------------------------------- physics
def nalkane_mw(n: int) -> float:
    return 12.011 * n + 1.008 * (2 * n + 2)


def nalkane_ik(n: int) -> str:
    return InchiToInchiKey(Chem.MolToInchi(Chem.MolFromSmiles("C" * n)))


def audited_nalkane_constants(path: Path | None = None) -> dict[int, dict]:
    """carbon number -> {tm_K, dh_J_mol, tier} for n-alkanes in the audited SLE constants table (results/v5/SLE2)."""
    c = pd.read_csv(path or ROOT / "results" / "v5" / "SLE2" / "constants_extended.csv")
    c = c[c.tier.notna() & (c.tier != "") & c.tm_K.notna() & c.dh_J_mol.notna()]
    by_ik = {nalkane_ik(n): n for n in range(1, 61)}
    return {by_ik[r.inchikey]: dict(tm_K=float(r.tm_K), dh_J_mol=float(r.dh_J_mol), tier=r.tier)
            for r in c.itertuples() if r.inchikey in by_ik}


def coutinho_nalkane_constants(g1) -> dict[int, dict]:
    ik = {n: nalkane_ik(n) for n in range(6, 37)}
    info = pd.DataFrame({"smiles": ["C" * n for n in ik], "cls": "n_paraffin", "n_carbon": [float(n) for n in ik]}, index=list(ik.values()))
    t = g1.get_pure_thermal_properties(info)
    return {n: dict(tm_K=t[k]["tm_K"], dh_J_mol=t[k]["dh_fus_J_mol"], tier="coutinho_table") for n, k in ik.items()}


def t_first(moles_by_cn: dict[int, float], total_moles: float, consts: dict[int, dict]) -> tuple[float, int | None]:
    """F0: first species to freeze.  moles_by_cn: n-alkane carbon number -> moles in the WHOLE liquid; total_moles: all moles
    (every other component is an inert solvent).  Returns (max_i T_i in K, carbon number of that species); (nan, None) if no
    n-alkane with constants is present.  T_i solves ln x_i = -(dH/R)(1/T - 1/Tm)  ->  1/T = 1/Tm - (R/dH) ln x_i."""
    best, arg = np.nan, None
    for n, mol in moles_by_cn.items():
        c = consts.get(n)
        if c is None or mol <= 0 or total_moles <= 0:
            continue
        x = mol / total_moles
        inv = 1.0 / c["tm_K"] - (R_GAS / c["dh_J_mol"]) * np.log(x)
        if inv <= 0:
            continue
        T = 1.0 / inv
        if not np.isfinite(best) or T > best:
            best, arg = T, n
    return best, arg


# ---------------------------------------------------------------------------------------------------------------- targets
def freezing_ajftd(long_df: pd.DataFrame) -> pd.DataFrame:
    """Per fuel: median of UNCENSORED freezing values.  Rows with any qualifier (e.g. '<') are censored and dropped; a fuel with
    no uncensored value is returned with value NaN and censored=True."""
    d = long_df[long_df.property == "freezing_point_C"].copy()
    d["fuel_id"] = d.fuel_id.astype(str)
    d["cens"] = d.qualifier.notna()
    rows = []
    for fid, g in d.groupby("fuel_id"):
        u = g[~g.cens & g.value.notna()]
        rows.append(dict(fuel_id=fid, value=float(u.value.median()) if len(u) else np.nan, n_tests=len(u), n_censored=int(g.cens.sum()),
                         censored_only=len(u) == 0, method="+".join(sorted({str(m).replace(" ", "") for m in u.method})) if len(u) else ""))
    return pd.DataFrame(rows)


def freezing_av3322(bulk: pd.DataFrame) -> pd.DataFrame:
    q = bulk.coa_freezing_point_C_qualifier.notna()
    return pd.DataFrame(dict(fuel_id=bulk.fuel_id.astype(str), value=np.where(q, np.nan, bulk.coa_freezing_point_C), n_tests=1,
                             n_censored=q.astype(int), censored_only=q, method="supplier CoA"))


def freezing_fuellib(props: pd.DataFrame) -> pd.DataFrame:
    d = props[(props.property == "freezing_point_C") & props.value.notna()]
    return pd.DataFrame(dict(fuel_id=d.fuel_id.astype(str), value=d.value.astype(float), n_tests=1, n_censored=0, censored_only=False,
                             method=d.method.astype(str)))


# ---------------------------------------------------------------------------------------------------------------- LOGO models
def feature_vector(bins: pd.DataFrame) -> np.ndarray:
    """bins: cls ('n_paraffin' | 'aromatic' | 'other'), n_carbon, mass.  Shares of the total GC mass."""
    tot = bins.mass.sum()
    na = bins[bins.cls == "n_paraffin"]
    f = lambda lo, hi: float(na[(na.n_carbon >= lo) & (na.n_carbon <= hi)].mass.sum() / tot)  # noqa: E731
    return np.array([f(9, 12), f(13, 16), f(17, 99), float(na.mass.sum() / tot), float(bins[bins.cls == "aromatic"].mass.sum() / tot)])


def logo_predict(y: np.ndarray, groups: np.ndarray, tfirst: np.ndarray, X: np.ndarray, train_mask: np.ndarray | None = None,
                 pred_idx: np.ndarray | None = None) -> pd.DataFrame:
    """Leave-one-group-out predictions of C, F1, P.  Training rows = train_mask & groups != scored group (so a scored group's own
    y is never used).  pred_idx: rows to predict (default all); rows outside train_mask (e.g. locked fuels) are scored with models
    trained on all of train_mask."""
    n = len(y)
    tm = np.ones(n, bool) if train_mask is None else train_mask
    idx = np.arange(n) if pred_idx is None else pred_idx
    out = {k: np.full(n, np.nan) for k in ("C", "F1", "P", "F1_a", "F1_b", "P_alpha")}
    for i in idx:
        tr = tm & (groups != groups[i])
        ytr = y[tr]
        out["C"][i] = ytr.mean()
        ok = tr & np.isfinite(tfirst)
        if ok.sum() >= 3 and np.isfinite(tfirst[i]) and np.ptp(tfirst[ok]) > 0:
            b, a = np.polyfit(tfirst[ok], y[ok], 1)
            out["F1"][i], out["F1_a"][i], out["F1_b"][i] = a + b * tfirst[i], a, b
        mu, sd = X[tr].mean(0), X[tr].std(0)
        sd[sd == 0] = 1.0
        rg = RidgeCV(alphas=ALPHAS).fit((X[tr] - mu) / sd, ytr)
        out["P"][i], out["P_alpha"][i] = float(rg.predict(((X[i] - mu) / sd)[None])[0]), rg.alpha_
    return pd.DataFrame(out)


# ---------------------------------------------------------------------------------------------------------------- statistics
def metrics_row(pred: np.ndarray, obs: np.ndarray) -> dict:
    m = np.isfinite(pred) & np.isfinite(obs)
    e = pred[m] - obs[m]
    if m.sum() == 0:
        return dict(n=0, bias=np.nan, mae=np.nan, rmse=np.nan, spearman=np.nan)
    rho = spearmanr(pred[m], obs[m]).statistic if m.sum() >= 3 and np.ptp(pred[m]) > 0 else np.nan
    return dict(n=int(m.sum()), bias=float(e.mean()), mae=float(np.abs(e).mean()), rmse=float(np.sqrt((e ** 2).mean())), spearman=float(rho))


def paired_bootstrap(abs_ref: np.ndarray, abs_mod: np.ndarray, n_boot: int = N_BOOT, seed: int = SEED) -> tuple[float, float, float]:
    """MAE(ref) - MAE(model), groups resampled with replacement (one row per group).  -> (diff, lo95, hi95)."""
    d = abs_ref - abs_mod
    rng = np.random.default_rng(seed)
    ix = rng.integers(0, len(d), size=(n_boot, len(d)))
    bs = d[ix].mean(1)
    return float(d.mean()), float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5))


def verdict(n: int, lo: float, hi: float) -> str:
    if n < MIN_GROUPS:
        return f"not evaluable (n={n} < {MIN_GROUPS})"
    if lo > 0:
        return "model better"
    if hi < 0:
        return "model worse"
    return "no detectable difference"


# ---------------------------------------------------------------------------------------------------------------- data
def load_module(name: str, rel: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / rel)
    m = importlib.util.module_from_spec(spec)
    sys.modules[name] = m
    spec.loader.exec_module(m)
    return m


def cls3(c: str) -> str:
    return "n_paraffin" if c in ("n_paraffin", "n-alkane") else ("aromatic" if c in AROM else "other")


def build_entities() -> dict[str, dict]:
    """entity id -> dict(source, fid, label, bins DataFrame[cls, n_carbon, mass], mapped DataFrame[ik, smiles, mass], unmapped)."""
    ents: dict[str, dict] = {}
    # AJFTD (G1 mapping, committed)
    comp = pd.read_csv(V5 / "ajftd_gcxgc_composition.csv")
    mapped = pd.read_csv(V5 / "ajftd_gcxgc_mapped.csv")
    msum = pd.read_csv(G1D / "mapping_summary.csv").set_index("fuel_id")
    for fid, g in comp.groupby("fuel_id"):
        if fid not in msum.index:
            continue
        b = pd.DataFrame(dict(cls=[cls3(c) for c in g.hc_class], n_carbon=g.n_carbon.astype(int), mass=g.mass_pct.astype(float)))
        mp = mapped[mapped.fuel_id == fid]
        ents[f"AJFTD:{fid}"] = dict(source="AJFTD", fid=str(fid), label=str(msum.loc[fid, "fuel_label"]), bins=b,
                                    mapped=pd.DataFrame(dict(ik=mp.inchikey, smiles=mp.smiles, frac=mp.fraction)),
                                    unmapped=float(msum.loc[fid, "unmapped_fraction"]))
    # CRC AV-33-22 (S2 representative map)
    comp = pd.read_csv(V5 / "crc_av3322_composition.csv")
    amap = pd.read_csv(S2 / "av3322_representative_map.csv")
    for fid, g in comp.groupby("fuel_id"):
        m = g.merge(amap[["group", "carbon_number", "rep_smiles", "rep_inchikey", "mapped"]], on=["group", "carbon_number"], how="left")
        tot = m.mass_pct.sum()
        mm = m[m.mapped.fillna(False).astype(bool)]
        agg = mm.groupby(["rep_inchikey", "rep_smiles"]).mass_pct.sum().reset_index()
        b = pd.DataFrame(dict(cls=[cls3("n_paraffin" if x == "n-alkane" else ("aromatic" if x in AROM else "other")) for x in g.group],
                              n_carbon=g.carbon_number.astype(int), mass=g.mass_pct.astype(float)))
        ents[f"AV-33-22:{fid}"] = dict(source="AV-33-22", fid=str(fid), label=str(fid), bins=b, unmapped=float(1 - agg.mass_pct.sum() / tot),
                                       mapped=pd.DataFrame(dict(ik=agg.rep_inchikey, smiles=agg.rep_smiles, frac=agg.mass_pct / agg.mass_pct.sum())))
    # FuelLib
    fl = pd.read_csv(V5 / "fuellib_gcxgc_compositions.csv")
    fmap = pd.read_csv(S2 / "fuellib_representative_map.csv")
    for fid, g in fl[fl.fuel_id.str.contains("posf")].groupby("fuel_id"):
        m = g.merge(fmap[["gc_bin", "n_carbon", "rep_smiles", "rep_inchikey", "mapped"]], on=["gc_bin", "n_carbon"], how="left")
        tot = m.fraction.sum()
        mm = m[m.mapped.fillna(False).astype(bool)]
        agg = mm.groupby(["rep_inchikey", "rep_smiles"]).fraction.sum().reset_index()
        b = pd.DataFrame(dict(cls=[cls3(c) for c in g.hc_class], n_carbon=g.n_carbon.astype(int), mass=g.fraction.astype(float)))
        ents[f"FuelLib:{fid}"] = dict(source="FuelLib", fid=fid, label=fid.replace("FUELLIB_posf", ""), bins=b,
                                      unmapped=float(1 - agg.fraction.sum() / tot),
                                      mapped=pd.DataFrame(dict(ik=agg.rep_inchikey, smiles=agg.rep_smiles, frac=agg.fraction / agg.fraction.sum())))
    return ents


def species_moles(ent: dict) -> tuple[dict[int, float], float, float]:
    """-> (n-alkane moles by true carbon number from the bins, total moles of the whole mapped liquid, non-n-alkane mapped moles).
    n-alkane bins are taken at their own carbon number (not the representative, which can sit up to 3 carbons away).  Other
    components come from the mapped representatives; unmapped non-n-alkane mass is not in the liquid (same as G1)."""
    b = ent["bins"]
    tot_mass = b.mass.sum()
    na = b[b.cls == "n_paraffin"].groupby("n_carbon").mass.sum()
    mol_na = {int(n): float(m / nalkane_mw(int(n))) for n, m in na.items() if m > 0}
    mapped_mass = tot_mass * (1 - ent["unmapped"])
    mp = ent["mapped"]
    other = 0.0
    for r in mp.itertuples():
        mol = Chem.MolFromSmiles(r.smiles)
        if set(r.smiles) <= {"C"}:
            continue                                          # n-alkane representative: already counted from the bins
        other += r.frac * mapped_mass / Descriptors.MolWt(mol)
    return mol_na, sum(mol_na.values()) + other, other


def ideal_svl_I(ent: dict, g1) -> float:
    """G1 model I: ideal SvL on the mapped composition (representatives), returns deg C (nan if no component qualifies)."""
    mp = ent["mapped"]
    mols = {r.ik: Chem.MolFromSmiles(r.smiles) for r in mp.itertuples()}
    info = pd.DataFrame(dict(smiles=[r.smiles for r in mp.itertuples()],
                             cls=[g1.hc_class(mols[r.ik]) for r in mp.itertuples()],
                             n_carbon=[float(sum(a.GetSymbol() == "C" for a in mols[r.ik].GetAtoms())) for r in mp.itertuples()]),
                        index=list(mp.ik)).loc[lambda d: ~d.index.duplicated()]
    tmap = g1.get_pure_thermal_properties(info)
    w = mp.frac.to_numpy(float)
    x = w / np.array([Descriptors.MolWt(mols[k]) for k in mp.ik])
    x = x / x.sum()
    return float(g1.ideal_sle_C(x, list(mp.ik), tmap))


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    g1 = load_module("g1_eval_s5", "scripts/v5/eval_ajftd_gcxgc_v5.py")
    ents = build_entities()
    cons = audited_nalkane_constants()
    consx = coutinho_nalkane_constants(g1)
    pd.DataFrame([dict(carbon_number=n, in_audited_table=n in cons, **(cons.get(n) or {}),
                       coutinho_tm_K=consx.get(n, {}).get("tm_K"), coutinho_dh_J_mol=consx.get(n, {}).get("dh_J_mol"))
                  for n in range(6, 27)]).to_csv(OUT / "s5_nalkane_constants.csv", index=False)

    # freezing targets
    long_aj = pd.read_csv(V5 / "ajftd_fuels_long.csv", low_memory=False)
    fz = {"AJFTD": freezing_ajftd(long_aj), "AV-33-22": freezing_av3322(pd.read_csv(V5 / "crc_av3322_bulk.csv")),
          "FuelLib": freezing_fuellib(pd.read_csv(V5 / "fuellib_properties.csv"))}
    fz = {k: v.set_index("fuel_id") for k, v in fz.items()}

    # groups and locks from S3
    s3f = pd.read_csv(S3 / "s3_fuels.csv")
    grp = dict(zip(s3f.entity, s3f.group))
    glock = dict(zip(s3f.entity, s3f.group_locked))

    rows, excl = [], []
    for eid, e in ents.items():
        src, fid = e["source"], e["fid"]
        if eid not in grp:
            continue
        f = fz[src].loc[fid] if fid in fz[src].index else None
        if f is None:
            continue                                           # no freezing record at all: not a candidate, not listed
        rec = dict(entity=eid, source=src, fuel_id=fid, label=e["label"], group=grp[eid], locked=bool(glock[eid]),
                   unmapped_fraction=e["unmapped"], freezing_C=f.value, n_tests=int(f.n_tests), n_censored=int(f.n_censored),
                   test_method=f.method)
        if bool(f.censored_only) or not np.isfinite(f.value):
            excl.append(dict(rec, reason="freezing value censored (qualifier) - no uncensored measurement"))
            continue
        if e["unmapped"] > MAX_UNMAPPED:
            excl.append(dict(rec, reason=f"unmapped composition {e['unmapped']:.3f} > {MAX_UNMAPPED}"))
            continue
        mol_na, tot_mol, other_mol = species_moles(e)
        for tag, c in (("", cons), ("x", consx)):
            T, arg = t_first(mol_na, tot_mol, c)
            rec[f"T_first{tag}_C"] = T - 273.15 if np.isfinite(T) else np.nan
            rec[f"first_species{tag}"] = arg
        rec["nalkane_mass_pct"] = float(e["bins"][e["bins"].cls == "n_paraffin"].mass.sum() / e["bins"].mass.sum() * 100)
        miss = [n for n in mol_na if n not in cons]
        rec["nalkane_mass_without_constants_pct"] = float(sum(e["bins"][(e["bins"].cls == "n_paraffin") & e["bins"].n_carbon.isin(miss)].mass) / e["bins"].mass.sum() * 100)
        rec["I_C"] = ideal_svl_I(e, g1)
        rec["_X"] = feature_vector(e["bins"])
        rows.append(rec)

    d = pd.DataFrame(rows)
    # one entity per group: scoreable members only, S3 source priority, then entity id; groups with a locked fuel are 'locked'
    d["prio"] = d.source.map(PRIORITY)
    d = d.sort_values(["group", "prio", "entity"])
    d["chosen"] = ~d.duplicated("group", keep="first")
    dup = d[~d.chosen][["entity", "group"]].rename(columns={"entity": "dropped_as_double_count"})
    for r in d[~d.chosen].itertuples():
        excl.append({**{k: getattr(r, k) for k in ("entity", "source", "fuel_id", "label", "group", "locked", "unmapped_fraction", "freezing_C",
                                                   "n_tests", "n_censored", "test_method")},
                     "reason": f"double-count group {r.group}: another member scored (S3 priority)"})
    d = d[d.chosen].drop(columns=["chosen", "prio"]).reset_index(drop=True)
    d["role"] = np.where(d.locked, "locked", "headline")

    y = d.freezing_C.to_numpy(float)
    g = d.group.to_numpy()
    X = np.vstack(d["_X"].to_numpy())
    head = (d.role == "headline").to_numpy()
    locked_idx = np.where(~head)[0]
    # pooled LOGO (train on all other headline groups); locked fuels scored from the models trained on all headline groups
    for tag, tcol in (("", "T_first_C"), ("x", "T_firstx_C")):
        res = logo_predict(y, g, d[tcol].to_numpy(float), X, train_mask=head)
        d[f"F1{tag}"] = res.F1
        d[f"F1{tag}_a"], d[f"F1{tag}_b"] = res.F1_a, res.F1_b
        if tag == "":
            d["C"], d["P"], d["P_alpha"] = res.C, res.P, res.P_alpha
    d["F0"] = d.T_first_C
    d["F0x"] = d.T_firstx_C
    d["I"] = d.I_C
    # within-source LOGO (sensitivity)
    for tag in ("C", "F1", "P"):
        d[f"{tag}_ws"] = np.nan
    for src in ("AJFTD", "AV-33-22"):
        m = (d.source == src).to_numpy() & head
        if m.sum() >= 5:
            sub = d[m]
            r = logo_predict(sub.freezing_C.to_numpy(float), sub.group.to_numpy(), sub.T_first_C.to_numpy(float), np.vstack(sub["_X"]))
            d.loc[m, "C_ws"], d.loc[m, "F1_ws"], d.loc[m, "P_ws"] = r.C.to_numpy(), r.F1.to_numpy(), r.P.to_numpy()
    # sensitivity: refit and rescore without the one very waxy fuel (AJFTD:131, 51 % n-alkane, freezing -8 C), a high-leverage point
    lev = (d.entity == "AJFTD:131").to_numpy()
    h2 = head & ~lev
    r2 = logo_predict(y, g, d.F0.to_numpy(float), X, train_mask=h2, pred_idx=np.where(h2)[0])
    d["C_x131"], d["F1_x131"], d["P_x131"] = r2.C, r2.F1, r2.P
    X_cols = ["x_nalk_C9_12", "x_nalk_C13_16", "x_nalk_C17p", "x_nalk_total", "x_aromatic"]
    for i, c in enumerate(X_cols):
        d[c] = X[:, i]
    d = d.drop(columns=["_X", "I_C", "T_first_C", "T_firstx_C"])
    d.to_csv(OUT / "s5_predictions.csv", index=False)
    pd.DataFrame(excl).to_csv(OUT / "s5_excluded.csv", index=False)
    d[["entity", "source", "fuel_id", "label", "group", "role", "freezing_C", "n_tests", "test_method", "unmapped_fraction", "nalkane_mass_pct",
       "nalkane_mass_without_constants_pct", "first_species", "first_speciesx"]
      ].to_csv(OUT / "s5_fuels.csv", index=False)

    # metrics
    sets = {"headline": d.role == "headline", "headline AJFTD": (d.role == "headline") & (d.source == "AJFTD"),
            "headline AV-33-22": (d.role == "headline") & (d.source == "AV-33-22"),
            "headline FuelLib": (d.role == "headline") & (d.source == "FuelLib"), "locked AJFTD": d.role == "locked",
            "headline excl. AJFTD:131": (d.role == "headline") & (d.entity != "AJFTD:131")}
    models = ["C", "I", "F0", "F1", "P", "F0x", "F1x", "C_ws", "F1_ws", "P_ws", "C_x131", "F1_x131", "P_x131"]
    mrows = []
    for sn, m in sets.items():
        s = d[m]
        for mod in models:
            ex = sn == "headline excl. AJFTD:131"
            if mod.endswith("_ws") and sn not in ("headline AJFTD", "headline AV-33-22"):
                continue
            if mod.endswith("_x131") != False and not ex:
                continue
            if ex and mod in ("C", "F1", "P", "C_ws", "F1_ws", "P_ws", "F1x"):
                continue                                            # this set uses the refitted *_x131 columns
            row = metrics_row(s[mod].to_numpy(float), s.freezing_C.to_numpy(float))
            if mod in ("C", "C_ws", "C_x131"):
                row["spearman"] = np.nan                            # a constant carries no ranking; LOGO mean is anti-correlated by construction
            mrows.append(dict(set=sn, model=mod, **row))
    pd.DataFrame(mrows).to_csv(OUT / "s5_metrics.csv", index=False)

    # paired bootstrap
    comps = [("R-S5-1", "C", m) for m in ("F0", "F1", "P")] + [("R-S5-2", "I", m) for m in ("F0", "F1")] + \
            [("context", "C", "I"), ("sensitivity", "C", "F0x"), ("sensitivity", "C", "F1x"), ("sensitivity", "I", "F0x"), ("sensitivity", "I", "F1x"),
             ("sensitivity within-source", "C_ws", "F1_ws"), ("sensitivity within-source", "C_ws", "P_ws")]
    brows = []
    for sn in ("headline", "headline AJFTD", "headline AV-33-22", "headline excl. AJFTD:131"):
        s = d[sets[sn]]
        for rd, ref, mod in comps:
            if (ref.endswith("_ws") or mod.endswith("_ws")) and sn in ("headline", "headline excl. AJFTD:131"):
                continue
            if sn == "headline excl. AJFTD:131":
                if rd not in ("R-S5-1", "R-S5-2"):
                    continue
                ref = ref + "_x131" if ref == "C" else ref
                mod = mod + "_x131" if mod in ("F1", "P") else mod
            ok = s[ref].notna() & s[mod].notna()
            ss = s[ok]
            n = len(ss)
            if n < 2:
                continue
            ar = (ss[ref] - ss.freezing_C).abs().to_numpy(float)
            am = (ss[mod] - ss.freezing_C).abs().to_numpy(float)
            diff, lo, hi = paired_bootstrap(ar, am)
            brows.append(dict(reading=rd, set=sn, reference=ref, model=mod, n_groups=n, mae_ref=ar.mean(), mae_model=am.mean(),
                              diff_ref_minus_model=diff, ci_lo=lo, ci_hi=hi, verdict=verdict(n, lo, hi)))
    b = pd.DataFrame(brows)
    b.to_csv(OUT / "s5_bootstrap.csv", index=False)

    # readings
    def get(sn, ref, mod):
        r = b[(b.set == sn) & (b.reference == ref) & (b.model == mod)]
        return r.iloc[0] if len(r) else None
    rr = []
    for sn in ("headline", "headline AJFTD", "headline AV-33-22"):
        for mod in ("F0", "F1", "P"):
            r = get(sn, "C", mod)
            if r is not None:
                rr.append(dict(reading="R-S5-1", set=sn, comparison=f"C vs {mod}", n_groups=r.n_groups, diff=r.diff_ref_minus_model, lo=r.ci_lo, hi=r.ci_hi,
                               met="MET" if r.verdict == "model better" else ("NOT EVALUABLE" if r.verdict.startswith("not evaluable") else "NOT MET"),
                               detail=r.verdict))
        for mod in ("F0", "F1"):
            r = get(sn, "I", mod)
            if r is not None:
                rr.append(dict(reading="R-S5-2", set=sn, comparison=f"I vs {mod}", n_groups=r.n_groups, diff=r.diff_ref_minus_model, lo=r.ci_lo, hi=r.ci_hi,
                               met="MET" if r.verdict == "model better" else ("NOT EVALUABLE" if r.verdict.startswith("not evaluable") else "NOT MET"),
                               detail=r.verdict))
    for mod in ("F1", "P"):
        a, v = get("headline AJFTD", "C", mod), get("headline AV-33-22", "C", mod)
        if a is None or v is None:
            continue
        ev = a.n_groups >= MIN_GROUPS and v.n_groups >= MIN_GROUPS
        met = a.verdict == "model better" and v.verdict == "model better"
        rr.append(dict(reading="R-S5-3 component", set="AJFTD and AV-33-22 separately", comparison=f"C vs {mod}", n_groups=f"{a.n_groups}/{v.n_groups}",
                       diff=f"{a.diff_ref_minus_model:+.2f}/{v.diff_ref_minus_model:+.2f}", lo=f"{a.ci_lo:+.2f}/{v.ci_lo:+.2f}", hi=f"{a.ci_hi:+.2f}/{v.ci_hi:+.2f}",
                       met="MET" if met else ("NOT EVALUABLE" if not ev else "NOT MET"), detail=f"AJFTD: {a.verdict}; AV-33-22: {v.verdict}"))
    pd.DataFrame(rr).to_csv(OUT / "s5_readings.csv", index=False)

    pd.set_option("display.width", 250, "display.max_columns", 40, "display.max_rows", 200)
    print(f"groups scored: headline {int(head.sum())} (AJFTD {int(((d.source == 'AJFTD') & head).sum())}, AV-33-22 {int(((d.source == 'AV-33-22') & head).sum())}, "
          f"FuelLib {int(((d.source == 'FuelLib') & head).sum())}); locked {len(locked_idx)}; excluded rows {len(excl)}")
    print(pd.DataFrame(mrows).round(3).to_string())
    print(b.round(3).to_string())
    print(pd.DataFrame(rr).to_string())


if __name__ == "__main__":
    main()
