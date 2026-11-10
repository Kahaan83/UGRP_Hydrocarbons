"""v4 Task 12: density, net heat of combustion and freezing point re-run on the v4 real-fuel corpus.

Same models and the same metric set as v3, but on v4 inputs only:
  * targets   : jet_properties_v4.csv rows with provenance == "raw_file" (CRC647, ASTM D4052/D4809-type tables)
  * mixtures  : compositions_v4.csv (CRC647 compositions are RECONSTRUCTED from measured aromatics + D86, R10, so
                every error below mixes composition-reconstruction error and physics error)
  * nothing in THIS script is fitted on the locked fuels, but see caveats: the legacy -14.58 C offset (v2 eval_cold_flow_v2) was
    fit on all 57 CRC fuels incl. locked and the library Coutinho parameters were calibrated on "Dauphin & CRC" data, so
    neither is an out-of-sample test.  The dev-only offset below is the clean one.  Pure-component inputs come from pure-compound data only (ThermoML pure-liquid
    densities, the v2 pure-GCV table, the library Tm/dHfus table).  No threshold or margin is chosen here, and the
    locked fuels (splits/v4) are scored as a separate never-tuned slice.

Models
  density  : ideal volume mixing  rho = 1/sum(w_i/rho_i(T)); rho_i(T) = linear-in-T fit of measured ThermoML pure-liquid
             densities (tier 1) or a class/carbon-number regression on those fits (tier 2).  No constant fallback.
  hoc      : mass-weighted group-additivity GCV (v2 table) -> NCV with the v2 conversion; primary variant uses the
             v2 constant 13.8 wt% H (identical to v2 C1), plus composition-derived and measured-H diagnostic variants.
  freezing : src/fuelblend/thermo/coutinho_wax.py as shipped: ideal WAT (mass fractions inside the log, as v3 R5 used it),
             WAT-14.58 offset, Coutinho eutectic; plus an ideal mole-fraction Schroeder-van Laar variant (this script).

Outputs: results/v4/B/{density,hoc,freezing}_v4.csv (per fuel [x T]) and {density,hoc,freezing}_summary_v4.csv.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from fuelblend.baselines.heat_of_combustion_v2 import (  # noqa: E402
    GroupAdditivityGCVModel, mixture_heat_of_combustion_v2)
from fuelblend.data.loaders_jet import included_jet_fuel_ids  # noqa: E402
from fuelblend.thermo.coutinho_wax import (  # noqa: E402
    R_GAS, get_pure_thermal_properties, predict_freezing_point_coutinho)

V4 = ROOT / "data" / "real_fuels" / "v4"
COMP = V4 / "compositions_v4.csv"
PROPS = V4 / "jet_properties_v4.csv"
INV = V4 / "inventory_v4.csv"
JETCLS = V4 / "jet_class_compositions.csv"
SPLITS = ROOT / "splits" / "v4" / "real_fuels_v4_splits.json"
THERMOML = ROOT / "data" / "processed" / "v2" / "thermoml_clean_hc.csv"
GCV_V2 = ROOT / "data" / "processed" / "v2" / "gcv_v2_MJ_kg.csv"
OUT = ROOT / "results" / "v4" / "B"

BOOT_SEED = 20261110
BOOT_B = 2000
T_REF_K = 288.71          # 15.56 C (60 F), the ASTM D4052 reference of the CRC647 densities
OFFSET_C = 14.58          # v3 R5 empirical WAT offset (kept for the like-for-like comparison)
DROPPED = ROOT / "data" / "quarantine" / "dropped_fuels.csv"
ALLOWED_COMPOSITION_PROVENANCE = {"raw_file", "reconstructed", "cited_typed", "reconstructed_r2"}   # anything else (e.g. quarantined) fails


# ----------------------------------------------------------------------------- metrics (same set as v2 C1 / v3 R5)
def bias_stats(pred, meas, B=BOOT_B, seed=BOOT_SEED) -> dict:
    """MAE, RMSE, MAPE %, bias (pred - meas) with a fuel-cluster percentile bootstrap (one cluster per fuel)."""
    p, m = np.asarray(pred, float), np.asarray(meas, float)
    ok = np.isfinite(p) & np.isfinite(m)
    p, m = p[ok], m[ok]
    n = int(p.size)
    if n == 0:
        return dict(N=0)
    e = p - m
    out = dict(N=n, MAE=float(np.mean(np.abs(e))), RMSE=float(np.sqrt(np.mean(e ** 2))),
               MAPE_pct=float(np.mean(np.abs(e / m)) * 100.0) if np.all(m != 0) else np.nan, bias=float(np.mean(e)))
    if n >= 5:
        rng = np.random.default_rng(seed)
        idx = rng.integers(0, n, size=(B, n))
        mae_b = np.abs(e[idx]).mean(axis=1)
        bias_b = e[idx].mean(axis=1)
        out.update(MAE_lo95=float(np.percentile(mae_b, 2.5)), MAE_hi95=float(np.percentile(mae_b, 97.5)),
                   bias_lo95=float(np.percentile(bias_b, 2.5)), bias_hi95=float(np.percentile(bias_b, 97.5)))
    return out


# ----------------------------------------------------------------------------- v4 inputs
def norm_mol(smiles: str):
    """Parse a SMILES and strip isotope labels.  The reconstructed CRC647 compositions use perdeuterated alkanes as
    GC-standard stand-ins for the n-alkanes; for bulk properties the protium molecule is the physical component."""
    m = Chem.MolFromSmiles(smiles)
    deut = any(a.GetIsotope() != 0 for a in m.GetAtoms())
    for a in m.GetAtoms():
        a.SetIsotope(0)
    return Chem.MolFromSmiles(Chem.MolToSmiles(m)), deut


def hc_class(m) -> str:  # same rules as extras/fuelblend_data/common.hc_class (not importable from src)
    if any(a.GetIsAromatic() for a in m.GetAtoms()):
        return "aromatic"
    if any(b.GetBondType() != Chem.BondType.SINGLE for b in m.GetBonds()):
        return "alkene_alkyne"
    if m.GetRingInfo().NumRings() > 0:
        return "cycloalkane"
    branched = any(a.GetSymbol() == "C" and a.GetDegree() > 2 for a in m.GetAtoms())
    return "iso_paraffin" if branched else "n_paraffin"


def select_fuels() -> tuple[list[str], dict]:
    """Fuels that may be evaluated: v4 inventory CRC647 fuels plus the included GC jet fuels, in the v4 split universe.
    The GC fuels (class x carbon-number compositions, no molecule list, n=1 headline) carry no bulk raw_file target
    that this module can score, so in practice the set is the CRC647 inventory fuels."""
    inv = pd.read_csv(INV)
    allowed = set(inv.fuel_id) | set(included_jet_fuel_ids(JETCLS))
    sp = json.load(open(SPLITS))
    role = sp["fuel_role"]
    return sorted(allowed), role


def load_inputs():
    allowed, role = select_fuels()
    comp = pd.read_csv(COMP)
    comp = comp[~comp.is_lumped.astype(bool) & (comp.inchikey != "LUMP_NO_REALISTIC_REP")].copy()   # lump = unexplained mass, renormalised away
    comp["fuel_id"] = comp.system_id.str.replace("fuel:", "", regex=False)
    props = pd.read_csv(PROPS)
    props = props[props.provenance == "raw_file"]
    fuels = sorted(set(comp.fuel_id) & set(props.fuel_id) & set(allowed))
    fuels = [f for f in fuels if f.startswith("CRC647")]          # GC fuels have no molecule list (see select_fuels)
    comp = comp[comp.fuel_id.isin(fuels)].copy()
    rows = []
    for r in comp.itertuples():
        m, deut = norm_mol(r.smiles)
        rows.append(dict(smiles_norm=Chem.MolToSmiles(m), inchikey_norm=Chem.MolToInchiKey(m), deuterated=deut,
                         cls=hc_class(m), n_carbon=sum(a.GetSymbol() == "C" for a in m.GetAtoms()),
                         mw=Descriptors.MolWt(m)))
    comp = pd.concat([comp.reset_index(drop=True), pd.DataFrame(rows)], axis=1)
    unex = comp.groupby("fuel_id").unexplained_fraction.first()
    return fuels, role, comp, props[props.fuel_id.isin(fuels)], unex


def assert_inputs_clean(df: pd.DataFrame, headline=("POSF-10325",)) -> None:
    """Every evaluated row: target provenance raw_file, composition provenance in the allowed v4 set, fuel not dropped
    by the Task 3 quarantine, fuel not in the n=1 GC headline set, and fuel in the v4 split universe."""
    bad_t = df[df.target_provenance != "raw_file"]
    assert bad_t.empty, f"non-raw_file targets evaluated: {bad_t.fuel_id.unique()}"
    for cell in df.composition_provenance.unique():
        assert set(cell.split(",")) <= ALLOWED_COMPOSITION_PROVENANCE, f"disallowed composition provenance {cell}"
    dropped = set(pd.read_csv(DROPPED).fuel_id)
    assert not (set(df.fuel_id) & dropped), "quarantine-dropped fuel evaluated"
    assert not (set(df.fuel_id) & set(headline)), "n=1 headline fuel must not enter statistics"
    assert (df.split_role != "not_in_split").all(), "fuel outside the v4 split universe"


def fuel_meta(f, role):
    return dict(fuel_id=f, split_role=role.get(f, "not_in_split"),
                slice="locked" if role.get(f) == "locked" else "dev")


def mass_frac(g):
    w = g.fraction.to_numpy(float)
    return w / w.sum()


# ----------------------------------------------------------------------------- pure density model
class PureDensity:
    """rho_i(T) = a_i + b_i (T - T_REF) from measured ThermoML pure liquids (tier 1), else class regression (tier 2)."""

    def __init__(self, comps: pd.DataFrame):
        t = pd.read_csv(THERMOML, low_memory=False)
        t = t[(t.n_comp == 1) & t.prop.str.contains("Mass density") & (t.prop_phase == "Liquid")
              & (t.P_kPa.isna() | (t.P_kPa <= 110)) & t["T_K"].between(250, 380)].copy()
        t["ik"] = t.inchis.map(lambda s: Chem.inchi.InchiToInchiKey(json.loads(s)[0]))
        self.t1 = {}
        for ik, g in t.groupby("ik"):
            T, y = g["T_K"].to_numpy(float), g["value"].to_numpy(float)
            if np.unique(np.round(T)).size < 3 or np.ptp(T) < 10:
                continue
            keep = np.ones(T.size, bool)
            for _ in range(2):          # robust: drop >3 MAD outliers (a few ThermoML records are off by 100s of kg/m3)
                b, a = np.polyfit(T[keep] - T_REF_K, y[keep], 1)
                res = y - (a + b * (T - T_REF_K))
                mad = np.median(np.abs(res[keep])) + 1e-9
                keep = np.abs(res) < 3 * 1.4826 * mad + 1.0
            b, a = np.polyfit(T[keep] - T_REF_K, y[keep], 1)
            self.t1[ik] = (a, b, float(T[keep].min()), float(T[keep].max()), int(keep.sum()))
        # tier 2: per-class regression of a_i on [1, 1/nC] over the tier-1 compounds, mean slope per class
        info = comps.drop_duplicates("inchikey_norm").set_index("inchikey_norm")
        rows = [(info.loc[k, "cls"], info.loc[k, "n_carbon"], v[0], v[1]) for k, v in self.t1.items() if k in info.index]
        self.fit_universe = pd.DataFrame(rows, columns=["cls", "n_carbon", "a", "b"])
        self.cls_fit = {}
        pooled = self.fit_universe
        for c in pooled.cls.unique():
            d = pooled[pooled.cls == c]
            self.cls_fit[c] = self._ridge(d if len(d) >= 4 else pooled)
        self.cls_fit["_all"] = self._ridge(pooled)

    @staticmethod
    def _ridge(d):
        X = np.c_[np.ones(len(d)), 1.0 / d.n_carbon.to_numpy(float)]
        beta = np.linalg.solve(X.T @ X + 1e-6 * np.eye(2), X.T @ d.a.to_numpy(float))
        return beta, float(d.b.mean())

    def rho(self, ik, cls, nc, T_K):
        """-> (rho kg/m3, tier).  tier 1 = measured compound fit (extrapolation <= 30 K beyond data, else tier 2)."""
        v = self.t1.get(ik)
        if v is not None and v[2] - 30 <= T_K <= v[3] + 30:
            return v[0] + v[1] * (T_K - T_REF_K), 1
        beta, b = self.cls_fit.get(cls, self.cls_fit["_all"])
        return float(beta[0] + beta[1] / nc + b * (T_K - T_REF_K)), 2


def eval_density(comp, props, role, fuels, pd_model):
    dens = props[props.property == "density_kg_m3"]
    rows = []
    for r in dens.itertuples():
        g = comp[comp.fuel_id == r.fuel_id]
        if g.empty or not np.isfinite(r.T_C):
            continue
        w = mass_frac(g)
        T_K = r.T_C + 273.15
        rho, tier = zip(*[pd_model.rho(ik, c, n, T_K) for ik, c, n in zip(g.inchikey_norm, g.cls, g.n_carbon)])
        rho, tier = np.array(rho), np.array(tier)
        pred = 1.0 / np.sum(w / rho)
        rows.append(dict(**fuel_meta(r.fuel_id, role), T_C=r.T_C, method=r.method, measured_kg_m3=r.value,
                         ideal_volume_kg_m3=pred, err_kg_m3=pred - r.value, rel_err_pct=(pred - r.value) / r.value * 100,
                         source_file=r.source_file, target_provenance=r.provenance,
                         composition_provenance=",".join(sorted(set(g.provenance))),
                         n_components=len(g), unexplained_fraction_dropped=float(g.unexplained_fraction.iloc[0]),
                         mass_frac_measured_rho=float(w[tier == 1].sum()), mass_frac_class_rho=float(w[tier == 2].sum()),
                         mass_frac_deuterated_standin=float(w[g.deuterated.to_numpy()].sum())))
    return pd.DataFrame(rows)


# ----------------------------------------------------------------------------- heat of combustion
def pure_gcv_table(comp):
    g = pd.read_csv(GCV_V2).set_index("inchikey")
    info = comp.drop_duplicates("inchikey_norm")
    # models for compounds that are not in the v2 table: refit the v2 group-additivity model on its Tier A rows
    tier_a = g[g.tier == "A"].reset_index()
    mw = {k: Descriptors.MolWt(Chem.MolFromSmiles(s)) for k, s in zip(tier_a.inchikey, tier_a.smiles)}
    model = GroupAdditivityGCVModel()
    model.fit_from_tier_a(tier_a.assign(mw=tier_a.inchikey.map(mw)))
    out = {}
    for r in info.itertuples():
        if r.inchikey_norm in g.index:
            out[r.inchikey_norm] = (float(g.loc[r.inchikey_norm, "value"]), "v2_table_" + g.loc[r.inchikey_norm, "tier"])
        else:
            out[r.inchikey_norm] = (model.predict_gcv(r.smiles_norm, r.mw), "refit_group_additivity")
    return out


def eval_hoc(comp, props, role, fuels):
    gcv = pure_gcv_table(comp)
    tgt = props[props.property == "net_heat_of_combustion_MJ_kg"].set_index("fuel_id")
    hyd = props[props.property == "hydrogen_wt_pct"].set_index("fuel_id")
    rows = []
    for f in fuels:
        if f not in tgt.index:
            continue
        g = comp[comp.fuel_id == f]
        w = mass_frac(g)
        gv = np.array([gcv[k][0] for k in g.inchikey_norm])
        tiers = np.array([gcv[k][1] for k in g.inchikey_norm])
        mw = g.mw.to_numpy(float)
        nc = g.n_carbon.to_numpy(float)
        std = mixture_heat_of_combustion_v2(w, gv, mw, comp_n_carbon=nc, use_standard_jet_h=True)
        cmp_ = mixture_heat_of_combustion_v2(w, gv, mw, comp_n_carbon=nc, use_standard_jet_h=False)
        h_meas = float(hyd.loc[f, "value"]) if f in hyd.index else np.nan
        ncv_meas_h = std["gcv_MJ_kg"] - 0.2122 * h_meas if np.isfinite(h_meas) else np.nan
        meas = float(tgt.loc[f, "value"])
        rows.append(dict(**fuel_meta(f, role), measured_ncv_MJ_kg=meas, ncv_v2_const_H_MJ_kg=std["ncv_MJ_kg"],
                         ncv_comp_H_MJ_kg=cmp_["ncv_MJ_kg"], ncv_measured_H_diag_MJ_kg=ncv_meas_h,
                         gcv_mix_MJ_kg=std["gcv_MJ_kg"], H_wt_pct_const=std["w_H_wt_pct"],
                         H_wt_pct_from_composition=cmp_["w_H_wt_pct"], H_wt_pct_measured=h_meas,
                         err_v2_const_H=std["ncv_MJ_kg"] - meas, err_comp_H=cmp_["ncv_MJ_kg"] - meas,
                         err_measured_H_diag=ncv_meas_h - meas if np.isfinite(ncv_meas_h) else np.nan,
                         source_file=tgt.loc[f, "source_file"], target_provenance=tgt.loc[f, "provenance"],
                         composition_provenance=",".join(sorted(set(g.provenance))), n_components=len(g),
                         unexplained_fraction_dropped=float(g.unexplained_fraction.iloc[0]),
                         mass_frac_tierA_gcv=float(w[tiers == "v2_table_A"].sum()),
                         mass_frac_deuterated_standin=float(w[g.deuterated.to_numpy()].sum())))
    return pd.DataFrame(rows)


# ----------------------------------------------------------------------------- freezing point
def thermal_map_for(comp):
    info = comp.drop_duplicates("inchikey_norm").set_index("inchikey_norm")[["smiles_norm", "cls", "n_carbon"]]
    info = info.rename(columns={"smiles_norm": "smiles"})
    return get_pure_thermal_properties(info)


def ideal_sle_mole_C(w, iks, mws, tmap):
    """Ideal Schroeder-van Laar liquidus with true MOLE fractions: max over crystallising solutes (no eutectic term)."""
    x = (np.asarray(w, float) / np.asarray(mws, float))
    x = x / x.sum()
    best = np.nan
    for xi, ik in zip(x, iks):
        t = tmap.get(ik)
        if t and np.isfinite(t["tm_K"]) and np.isfinite(t["dh_fus_J_mol"]) and xi > 1e-4:
            inv = 1.0 / t["tm_K"] - R_GAS / t["dh_fus_J_mol"] * np.log(xi)
            if inv > 0:
                T = 1.0 / inv - 273.15
                best = T if not np.isfinite(best) else max(best, T)
    return best


def eval_freezing(comp, props, role, fuels):
    tmap = thermal_map_for(comp)
    tgt = props[props.property == "freezing_point_C"].set_index("fuel_id")
    rows = []
    for f in fuels:
        if f not in tgt.index:
            continue
        g = comp[comp.fuel_id == f]
        w = mass_frac(g)
        iks = g.inchikey_norm.tolist()
        t_cout, de, meta = predict_freezing_point_coutinho(w, iks, tmap)
        no_cryst = t_cout <= -900
        wat = meta.get("ideal_wat_C", np.nan)
        meas = float(tgt.loc[f, "value"])
        sle_mole = ideal_sle_mole_C(w, iks, g.mw.to_numpy(float), tmap)
        n_par = g.cls.eq("n_paraffin").to_numpy()
        rows.append(dict(**fuel_meta(f, role), measured_C=meas, wat_ideal_library_C=wat,
                         wat_minus_legacy_offset_C=wat - OFFSET_C if np.isfinite(wat) else np.nan,
                         coutinho_C=np.nan if no_cryst else t_cout, ideal_sle_mole_C=sle_mole,
                         delta_eutectic_C=de, no_crystallizing_component=bool(no_cryst),
                         paraffin_mass_fraction=float(w[n_par].sum()), paraffin_std_nc=meta.get("paraffin_std_nc", np.nan),
                         err_wat_ideal=wat - meas, err_wat_minus_legacy_offset=wat - OFFSET_C - meas,
                         err_coutinho=np.nan if no_cryst else t_cout - meas, err_ideal_sle_mole=sle_mole - meas,
                         source_file=tgt.loc[f, "source_file"], target_provenance=tgt.loc[f, "provenance"],
                         composition_provenance=",".join(sorted(set(g.provenance))), n_components=len(g),
                         unexplained_fraction_dropped=float(g.unexplained_fraction.iloc[0]),
                         mass_frac_deuterated_standin=float(w[g.deuterated.to_numpy()].sum())))
    df = pd.DataFrame(rows)
    # Dev-only offset: mean (ideal WAT - measured) over the DEV fuels; locked fuels are scored with this constant and
    # never contribute to it.  (Dev rows are then in-sample for this one constant; locked rows are out-of-sample.)
    dv = df[df.slice == "dev"]
    off = float((dv.wat_ideal_library_C - dv.measured_C).mean())
    df["dev_fit_offset_C"] = off
    df["wat_minus_devfit_offset_C"] = df.wat_ideal_library_C - off
    df["err_wat_minus_devfit_offset"] = df.wat_minus_devfit_offset_C - df.measured_C
    return df


# ----------------------------------------------------------------------------- summaries
def summarize(df, cols: dict, extra=None):
    """cols: model -> error column; one row per (slice, model); slices: dev, locked, dev+locked (reported separately)."""
    out = []
    for sl, d in [("dev", df[df.slice == "dev"]), ("locked", df[df.slice == "locked"]), ("all_dev+locked", df)]:
        for model, (ecol, mcol) in cols.items():
            e = d[ecol].to_numpy(float)
            meas = d[mcol].to_numpy(float)
            s = bias_stats(e + meas, meas)
            s.update(slice=sl, model=model)
            if extra:
                s.update(extra(d, ecol))
            out.append(s)
    cols_first = ["slice", "model", "N"]
    s = pd.DataFrame(out)
    return s[cols_first + [c for c in s.columns if c not in cols_first]]


def freeze_extra(d, ecol):
    e = d[ecol].to_numpy(float)
    e = e[np.isfinite(e)]
    return dict(optimistic_gt3C_rate=float(np.mean(e < -3.0)) if e.size else np.nan,
                n_optimistic_gt3C=int(np.sum(e < -3.0)), worst_optimistic_C=float(e.min()) if e.size else np.nan)


def main(argv=None):
    global COMP, OUT
    ap = argparse.ArgumentParser()
    ap.add_argument("--compositions", type=Path, default=COMP, help="composition table (default: the v4 file)")
    ap.add_argument("--out", type=Path, default=OUT, help="output directory (default: results/v4/B)")
    a = ap.parse_args(argv)
    COMP, OUT = a.compositions, a.out
    OUT.mkdir(parents=True, exist_ok=True)
    fuels, role, comp, props, unex = load_inputs()
    print(f"{len(fuels)} v4 CRC647 fuels with composition + raw_file targets; "
          f"locked={sum(role.get(f) == 'locked' for f in fuels)} dev={sum(role.get(f) != 'locked' for f in fuels)}")

    dens = eval_density(comp, props, role, fuels, PureDensity(comp))
    assert_inputs_clean(dens)
    dens.to_csv(OUT / "density_v4.csv", index=False)
    prim = dens[np.isclose(dens.T_C, 15.56)]
    ds = pd.concat([
        summarize(prim, {"ideal_volume": ("err_kg_m3", "measured_kg_m3")}).assign(scope="15.56C D4052 (v3-comparable)"),
        summarize(dens, {"ideal_volume": ("err_kg_m3", "measured_kg_m3")}).assign(scope="all T (-0.2..69.5C; rows not fuels)"),
        summarize(prim[prim.mass_frac_measured_rho >= 0.8], {"ideal_volume": ("err_kg_m3", "measured_kg_m3")}).assign(
            scope="15.56C, >=80% of mass has a measured pure density"),
    ])
    ds.to_csv(OUT / "density_summary_v4.csv", index=False)

    hoc = eval_hoc(comp, props, role, fuels)
    assert_inputs_clean(hoc)
    hoc.to_csv(OUT / "hoc_v4.csv", index=False)
    hs = summarize(hoc, {"group_additivity_v2_constH13.8 (v2-comparable)": ("err_v2_const_H", "measured_ncv_MJ_kg"),
                         "group_additivity_composition_H": ("err_comp_H", "measured_ncv_MJ_kg"),
                         "group_additivity_measuredH_DIAGNOSTIC": ("err_measured_H_diag", "measured_ncv_MJ_kg")})
    hs.to_csv(OUT / "hoc_summary_v4.csv", index=False)

    fr = eval_freezing(comp, props, role, fuels)
    assert_inputs_clean(fr)
    fr.to_csv(OUT / "freezing_v4.csv", index=False)
    fs = summarize(fr, {"wat_ideal_library (mass-fraction ln)": ("err_wat_ideal", "measured_C"),
                        "wat_minus_devfit_offset (fit on 45 dev only; locked out-of-sample)": ("err_wat_minus_devfit_offset", "measured_C"),
                        "wat_minus_14.58 (LEGACY constant, fit on all 57 incl. locked; NOT out-of-sample)": ("err_wat_minus_legacy_offset", "measured_C"),
                        "coutinho_library (params calibrated on Dauphin+CRC; locked NOT a clean holdout)": ("err_coutinho", "measured_C"),
                        "ideal_sle_mole_fraction": ("err_ideal_sle_mole", "measured_C")}, extra=freeze_extra)
    fs.to_csv(OUT / "freezing_summary_v4.csv", index=False)

    pd.set_option("display.width", 250)
    pd.set_option("display.max_columns", 30)
    for name, s in [("density", ds), ("hoc", hs), ("freezing", fs)]:
        print(f"\n== {name}")
        print(s.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
