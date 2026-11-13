"""v4 Task 17: property-model wrappers for ``engine_v4`` built from the Task 4-15 library code.

Nothing here substitutes a number.  Pure-component constants come from ``data/processed/v4/pure_props_v4.csv`` and
ONLY rows with ``status == 'verified'`` are used (ruling R16); any other status raises ``ModelUnavailable`` naming the
compound, property and status, which the engine reports as ``uncertified`` for that constraint.

Primary model per property and the yaml margin it is judged by (all keys are in ``configs/uncertainty_v4.yaml``):

  density     ideal volume mixing of verified rho15              density_kg_m3.lab_mixtures.ideal_volume.over / .under (side-specific; hydrocarbons only)
  flash       Le Chatelier + Antoine + LFL, gamma = 1 (Task 9)   flash_point_C.lab_mixtures.ideal  MAX  D93.ideal
  viscosity   raw Refutas (known compositions; review fix 4)    viscosity_kin_cSt.lab_mixtures.raw_refutas.<T bin> (cold bins inf -> uncertified;
              recalibrated Refutas and ridge are alternates)    hydrocarbons only
  freezing    ideal Schroeder-van Laar (library, mass fractions) freezing_point_C.sle_alkane_blends.ideal_svl  (inf)
  nhoc        mass-weighted group-additivity GCV (Task 12)       heat_of_combustion_net_MJ_kg...constH13.8
  smoke       AJFTD bulk regression on aromatics (Task 14)       smoke_point_mm.ajftd100_bulk.selected_in_fold  (inf)
  ron/mon/dcn linear-by-volume of verified pure values           NO margin entry (octane models are base + one
                                                                  blendstock only; Task 15) -> never certified
  dvpe        ideal-Raoult sum x_i Psat(37.8 C)                   NO margin (DVPE data units-unconfirmed, excluded)

Encoder (Task 16) lost to the structured ridge and is NOT used.
"""
from __future__ import annotations

import math
from functools import lru_cache
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, inchi as rdinchi

from fuelblend.baselines import flash_unifac as fu
from fuelblend.baselines import sooting
from fuelblend.baselines.heat_of_combustion_v2 import GroupAdditivityGCVModel, mixture_heat_of_combustion_v2
from fuelblend.design.engine_v4 import ModelUnavailable, PropertyModel, viscosity_t_bin
from fuelblend.mixture.cold_flow import compound_meta
from fuelblend.mixture.common import CLASSES, CompoundRegistry, PurePropertyProvider
from fuelblend.mixture.structured_residual import StructuredResidual
from fuelblend.thermo.coutinho_wax import ideal_wat_schroeder, predict_freezing_point_coutinho

ROOT = Path(__file__).resolve().parents[3]
PURE_PROPS = ROOT / "data" / "processed" / "v4" / "pure_props_v4.csv"
COMPOUNDS = ROOT / "data" / "unified" / "compounds.csv"
POOL_GRID = ROOT / "data" / "cache" / "pure_grid.npz"
VISC_PRED = ROOT / "results" / "v4" / "V" / "viscosity_predictions_v4.csv"
VISC_HP = ROOT / "results" / "v4" / "V" / "selected_hyperparameters_v4.csv"
GCV_V2 = ROOT / "data" / "processed" / "v2" / "gcv_v2_MJ_kg.csv"
R_GAS = 8.314462618
FEATS = [f"frac_{c}" for c in CLASSES] + ["nC_mean", "nC_std"]
_NAPHTHALENE_CORE = Chem.MolFromSmarts("c1ccc2ccccc2c1")
CRC_CLASSES = set(CLASSES)          # classes present in the CRC 647 reconstructed compositions (validation data)


class Library:
    """All static data the model wrappers need; build once per process."""

    def __init__(self, pure_props=None):
        # pure_props: alternative constants table with the pure_props_v4.csv columns (v5 overlay); default is v4
        pp = pd.read_csv(PURE_PROPS if pure_props is None else pure_props)
        self.pp = pp
        self.smiles = pp.drop_duplicates("compound").set_index("compound").smiles.to_dict()
        self._tab = {(r.compound, r.prop): (r.value, r.status) for r in pp.itertuples()}
        self.reg = CompoundRegistry(COMPOUNDS)
        _, self.cls_all, self.nc_all = compound_meta(COMPOUNDS)
        self.provider = PurePropertyProvider(self.reg, POOL_GRID, verbose=False)
        self._ik = {}
        self._fit_viscosity()
        self._fit_gcv()

    # ---------------------------------------------------------------- compound metadata
    def ik(self, name):
        if name not in self._ik:
            self._ik[name] = rdinchi.MolToInchiKey(Chem.MolFromSmiles(self.smiles[name]))
        return self._ik[name]

    def idx(self, name):
        return self.reg.ik2i.get(self.ik(name))

    def cls(self, name):
        j = self.idx(name)
        return self.cls_all[j] if j is not None else "unknown"

    def nc(self, name):
        j = self.idx(name)
        return float(self.nc_all[j]) if j is not None else float("nan")

    def v(self, name, prop):
        """Verified pure constant, else ModelUnavailable (ruling R16: no silent substitution)."""
        val, st = self._tab.get((name, prop), (float("nan"), "absent"))
        if st != "verified" or not math.isfinite(float(val)):
            raise ModelUnavailable(f"{name}.{prop} status={st}")
        return float(val)

    def classes(self, pool):
        return {n: self.cls(n) for n in pool}

    def is_hydrocarbon(self, name) -> bool:
        """True iff the compound contains only C and H (the lab margins were measured on hydrocarbon mixtures only)."""
        mol = Chem.MolFromSmiles(self.smiles[name])
        return mol is not None and all(a.GetSymbol() in ("C", "H") for a in Chem.AddHs(mol).GetAtoms())

    def require_hydrocarbons(self, pool, w, prop):
        """Lab margins are hydrocarbon-only: any component with w > 0 that is not a hydrocarbon -> ModelUnavailable."""
        bad = [n for n, wi in zip(pool, np.asarray(w, float)) if wi > 0 and not self.is_hydrocarbon(n)]
        if bad:
            raise ModelUnavailable(f"{prop}: lab-mixture margins are hydrocarbon-only; non-hydrocarbon component(s) {bad} "
                                   "have no validated margin")

    # ---------------------------------------------------------------- fits from stored Task 11 outputs
    def _fit_viscosity(self):
        d = pd.read_csv(VISC_PRED)
        d = d[~d.locked.astype(bool)]                       # dev fuels only: locked fuels never fitted
        resid = np.log(d.y.to_numpy()) - np.log(d.nu_refutas.to_numpy())
        hp = pd.read_csv(VISC_HP)
        alpha = float(hp[(hp.model == "Structured ridge") & (hp.fold == "FINAL(dev->locked)")].hp.iloc[0])
        self.ridge = StructuredResidual(FEATS, "ridge", alpha=alpha).fit(d[FEATS].to_numpy(float), d.T_K.to_numpy(), resid)
        B = np.column_stack([np.ones(len(d)), 1000.0 / d.T_K.to_numpy()])
        self.recal_beta = np.linalg.lstsq(B, resid, rcond=None)[0]
        self.ridge_env = {k: (float(d[k].min()), float(d[k].max())) for k in FEATS}
        self.visc_fit_info = dict(n_rows=int(len(d)), n_fuels=int(d.fuel_id.nunique()), ridge_alpha=alpha)

    def _fit_gcv(self):
        g = pd.read_csv(GCV_V2)
        g = g[g.property == "gcv_MJ_kg"].drop_duplicates("inchikey").set_index("inchikey")
        self.gcv_tab = g
        tier_a = g[g.tier == "A"].reset_index()
        mw = {k: Descriptors.MolWt(Chem.MolFromSmiles(s)) for k, s in zip(tier_a.inchikey, tier_a.smiles)}
        self.gcv_model = GroupAdditivityGCVModel()
        self.gcv_model.fit_from_tier_a(tier_a.assign(mw=tier_a.inchikey.map(mw)))

    def gcv(self, name):
        ik = self.ik(name)
        if ik in self.gcv_tab.index:
            return float(self.gcv_tab.loc[ik, "value"])
        return self.gcv_model.predict_gcv(self.smiles[name], self.v(name, "MW"))

    def thermal_map(self, pool):
        out = {}
        for n in pool:
            tm = self.v(n, "Tm")
            dh = self.v(n, "dHfus")
            out[self.ik(n)] = dict(tm_K=tm + 273.15, dh_fus_J_mol=dh * 1000.0, n_carbon=self.nc(n), cls=self.cls(n))
        return out


@lru_cache(maxsize=1)
def library(pure_props=None) -> Library:
    return Library(pure_props)


# ------------------------------------------------------------------------------------------------ helpers
def _norm(w):
    w = np.asarray(w, float)
    return w / w.sum()


def _vol_frac(L, pool, w):
    rho = np.array([L.v(n, "rho15") for n in pool])
    v = _norm(w) / rho
    return v / v.sum()


def _mole_frac(L, pool, w):
    mw = np.array([L.v(n, "MW") for n in pool])
    x = _norm(w) / mw
    return x / x.sum()


# ------------------------------------------------------------------------------------------------ model fns
def make_fns(L: Library) -> dict:
    def density(pool, w, T_C):
        L.require_hydrocarbons(pool, w, "density")
        rho = np.array([L.v(n, "rho15") for n in pool])
        return float(1.0 / np.sum(_norm(w) / rho))

    def _flash(pool, w, T_C, unifac):
        sm = [L.smiles[n] for n in pool]
        bad = [n for n, s in zip(pool, sm) if fu._lookup(s, None) is None]
        if bad:
            raise ModelUnavailable(f"no Antoine/LFL entry in compounds_anchored_vp.csv for {bad}")
        x = _mole_frac(L, pool, w)
        with np.errstate(all="ignore"):
            T = fu.flash_point_unifac(sm, x) if unifac else fu.flash_point_ideal(sm, x)
        return float(T)

    def flash_ideal(pool, w, T_C):
        return _flash(pool, w, T_C, False)

    def flash_unifac(pool, w, T_C):
        return _flash(pool, w, T_C, True)

    def _visc_inputs(pool, w, T_C):
        T_K = T_C + 273.15
        idx = []
        for n in pool:
            j = L.idx(n)
            if j is None or not L.provider.valid[j] or L.cls(n) not in CLASSES or not np.isfinite(L.nc(n)):
                raise ModelUnavailable(f"{n}: no pure eta(T)/rho(T) grid or class outside the viscosity model")
            idx.append(j)
        idx = np.array(idx)
        ln_eta, rho = L.provider.lookup(idx[None, :], T_K)
        ln_eta, rho = ln_eta[0], rho[0]
        nu_i = 1000.0 * np.exp(ln_eta) / rho
        vbn = 14.534 * np.log(np.log(nu_i + 0.8)) + 10.975
        wn = _norm(w)
        vbn_m = float((wn * vbn).sum())
        nu = float(np.exp(np.exp((vbn_m - 10.975) / 14.534)) - 0.8)
        cls = [L.cls(n) for n in pool]
        feats = {f"frac_{c}": float(sum(wi for wi, ci in zip(wn, cls) if ci == c)) for c in CLASSES}
        nc = np.array([L.nc(n) for n in pool])
        m = float((wn * nc).sum())
        feats["nC_mean"], feats["nC_std"] = m, float(np.sqrt((wn * (nc - m) ** 2).sum()))
        return nu, feats, T_K

    def visc_raw(pool, w, T_C):
        L.require_hydrocarbons(pool, w, "viscosity")
        return _visc_inputs(pool, w, T_C)[0]

    def visc_recal(pool, w, T_C):
        nu, _, T_K = _visc_inputs(pool, w, T_C)
        return float(nu * math.exp(L.recal_beta[0] + L.recal_beta[1] * 1000.0 / T_K))

    def visc_ridge_unguarded(pool, w, T_C):
        """DIAGNOSTIC ONLY (never a design model): the Task 11 ridge with no envelope check."""
        nu, f, T_K = _visc_inputs(pool, w, T_C)
        x = np.array([[f[k] for k in FEATS]])
        return float(nu * math.exp(float(L.ridge.predict(x, np.array([T_K]))[0])))

    def visc_ridge(pool, w, T_C):
        nu, f, T_K = _visc_inputs(pool, w, T_C)
        out = {k: (round(f[k], 3), L.ridge_env[k]) for k in FEATS
               if not (L.ridge_env[k][0] - 1e-9 <= f[k] <= L.ridge_env[k][1] + 1e-9)}
        if out:       # a linear-in-features residual extrapolates without bound; no value is returned outside the envelope
            raise ModelUnavailable(f"structured ridge applied outside its CRC dev-fuel feature envelope: {out}")
        return visc_ridge_unguarded(pool, w, T_C)

    def freezing_ideal(pool, w, T_C):
        tmap = L.thermal_map(pool)
        T = ideal_wat_schroeder(_norm(w), [L.ik(n) for n in pool], tmap)
        if not math.isfinite(T):
            raise ModelUnavailable("no crystallising component with verified Tm and dHfus")
        return float(T)

    def freezing_coutinho(pool, w, T_C):
        tmap = L.thermal_map(pool)
        T, _, meta = predict_freezing_point_coutinho(_norm(w), [L.ik(n) for n in pool], tmap)
        if not math.isfinite(T) or T <= -900:
            raise ModelUnavailable("no crystallising component with verified Tm and dHfus")
        return float(T)

    def freezing_ideal_mole(pool, w, T_C):
        tmap = L.thermal_map(pool)
        x = _mole_frac(L, pool, w)
        best = None
        for xi, n in zip(x, pool):
            t = tmap[L.ik(n)]
            if xi > 1e-4:
                inv = 1.0 / t["tm_K"] - R_GAS / t["dh_fus_J_mol"] * math.log(xi)
                if inv > 0:
                    T = 1.0 / inv - 273.15
                    best = T if best is None else max(best, T)
        if best is None:
            raise ModelUnavailable("no crystallising component")
        return float(best)

    def _nhoc(pool, w, const_h):
        gcv = [L.gcv(n) for n in pool]
        mw = [L.v(n, "MW") for n in pool]
        nc = [L.nc(n) for n in pool]
        return mixture_heat_of_combustion_v2(_norm(w), gcv, mw, comp_n_carbon=nc, use_standard_jet_h=const_h)["ncv_MJ_kg"]

    def nhoc_const_h(pool, w, T_C):
        return float(_nhoc(pool, w, True))

    def nhoc_comp_h(pool, w, T_C):
        return float(_nhoc(pool, w, False))

    def arom_wt(pool, w, T_C):
        wn = _norm(w)
        return float(100.0 * sum(wi for wi, n in zip(wn, pool) if L.cls(n) == "aromatic"))

    def arom_vol(pool, w, T_C):
        v = _vol_frac(L, pool, w)
        return float(100.0 * sum(vi for vi, n in zip(v, pool) if L.cls(n) == "aromatic"))

    def naph_vol(pool, w, T_C):
        """Vol % of naphthalene-type components (naphthalene core, i.e. two fused benzene rings; D1840 'naphthalenes').

        Tetralin, indane, biphenyl and decalin are NOT naphthalene-type.  A component whose class is unknown or whose SMILES
        does not parse makes the quantity unknown (ModelUnavailable): it must never read as 0."""
        for n in pool:
            if L.cls(n) == "unknown" or Chem.MolFromSmiles(L.smiles[n]) is None:
                raise ModelUnavailable(f"naphthalenes: no class/ring data for {n}")
        v = _vol_frac(L, pool, w)
        return float(100.0 * sum(vi for vi, n in zip(v, pool) if Chem.MolFromSmiles(L.smiles[n]).HasSubstructMatch(_NAPHTHALENE_CORE)))

    def ethanol_vol(pool, w, T_C):
        v = _vol_frac(L, pool, w)
        return float(100.0 * sum(vi for vi, n in zip(v, pool) if n == "ethanol"))

    def smoke_bulk(pool, w, T_C):
        return float(sooting.smoke_point_from_bulk(arom_vol(pool, w, T_C)))

    def smoke_ysi(pool, w, T_C):
        sm = [L.smiles[n] for n in pool]
        x = _mole_frac(L, pool, w)
        try:
            y = sooting.ysi_mix(sm, x)
        except KeyError as e:
            raise ModelUnavailable(str(e))
        return float(sooting.smoke_point_from_ysi(y, sooting.mixture_mw(sm, x)))

    def _lin_vol(prop):
        def fn(pool, w, T_C):
            v = _vol_frac(L, pool, w)
            return float(sum(vi * L.v(n, prop) for vi, n in zip(v, pool)))
        return fn

    def dvpe(pool, w, T_C):
        x = _mole_frac(L, pool, w)
        return float(sum(xi * L.v(n, "Psat37_8") for xi, n in zip(x, pool)))

    return dict(density=density, flash_ideal=flash_ideal, flash_unifac=flash_unifac, visc_raw=visc_raw,
                visc_recal=visc_recal, visc_ridge=visc_ridge, visc_ridge_unguarded=visc_ridge_unguarded, freezing_ideal=freezing_ideal,
                freezing_coutinho=freezing_coutinho, freezing_ideal_mole=freezing_ideal_mole,
                nhoc_const_h=nhoc_const_h, nhoc_comp_h=nhoc_comp_h, arom_wt=arom_wt, arom_vol=arom_vol, naph_vol=naph_vol,
                ethanol_vol=ethanol_vol, smoke_bulk=smoke_bulk, smoke_ysi=smoke_ysi, ron=_lin_vol("RON"),
                mon=_lin_vol("MON"), dcn=_lin_vol("DCN"), dvpe=dvpe)


def build_models(L: Library | None = None) -> dict:
    """The ``models`` dict consumed by ``engine_v4.design`` / ``robustness_table``."""
    L = L or library()
    f = make_fns(L)
    FL = ["flash_point_C.lab_mixtures.ideal", "flash_point_C.real_fuels_D93_reconstructed.ideal"]
    FLU = ["flash_point_C.lab_mixtures.unifac", "flash_point_C.real_fuels_D93_reconstructed.unifac"]
    oct_why = ("octane margins exist only for base-gasoline + ONE blendstock models (Task 15, 18 blendstock groups); "
               "a multi-component pure-compound recipe has no validated octane model")
    JET = {"jet"}                   # every yaml margin below was validated on jet fuels (CRC 647 / AJFTD / jet-range lab blends)

    def P(*a, **kw):
        m = PropertyModel(*a, **kw)
        if not m.by_construction and "certified_kinds" not in kw:
            m.certified_kinds = JET
        return m

    props = {
        # review fix 4: the pools are blends of KNOWN pure compounds, so the lab-mixture margins apply (not the CRC
        # reconstruction margins); they are not jet specific (certified_kinds=None) and are hydrocarbon-only.
        "density": P("ideal_volume_rho15", f["density"],
                     margin_keys_by_side={"over": ["density_kg_m3.lab_mixtures.ideal_volume.over"],
                                          "under": ["density_kg_m3.lab_mixtures.ideal_volume.under"]},
                     validated_classes=CRC_CLASSES, certified_kinds=None,
                     notes="lab-mixture side-specific margins; non-hydrocarbon component with w>0 -> unavailable"),
        "flash": P("le_chatelier_ideal", f["flash_ideal"], FL, validated_classes=CRC_CLASSES,
                   notes="margin = max(lab ideal, D93 ideal) applied_margin"),
        "viscosity": P("refutas_raw", f["visc_raw"],
                       margin_keys_fn=lambda T_C: [f"viscosity_kin_cSt.lab_mixtures.raw_refutas.{viscosity_t_bin(T_C)}"],
                       T_valid_C=None, validated_classes=CRC_CLASSES, certified_kinds=None,
                       notes="raw Refutas of known pure compounds + lab-mixture margin of the temperature bin "
                             "(cold bins have no valid margin -> uncertified); recalibration/ridge are alternates"),
        "freezing": P("ideal_svl_library", f["freezing_ideal"], ["freezing_point_C.sle_alkane_blends.ideal_svl"],
                      validated_classes=CRC_CLASSES),
        "nhoc": P("gcv_group_additivity_constH", f["nhoc_const_h"],
                  ["heat_of_combustion_net_MJ_kg.crc647.group_additivity_v2_constH13.8"], validated_classes=CRC_CLASSES),
        "aromatics_wt": P("class_mass_fraction", f["arom_wt"], by_construction=True),
        "aromatics_vol": P("class_volume_fraction", f["arom_vol"], by_construction=True),
        "naphthalenes_vol_pct": P("naphthalene_core_volume_fraction", f["naph_vol"], by_construction=True),
        "ethanol_vol": P("ethanol_volume_fraction", f["ethanol_vol"], by_construction=True),
        "smoke_point": P("ajftd_bulk_regression", f["smoke_bulk"], ["smoke_point_mm.ajftd100_bulk.selected_in_fold"],
                         validated_classes=CRC_CLASSES),
        "ron": P("linear_volume_pure_RON", f["ron"], margin_missing_reason=oct_why),
        "mon": P("linear_volume_pure_MON", f["mon"], margin_missing_reason=oct_why),
        "dcn": P("linear_volume_pure_DCN", f["dcn"],
                 margin_missing_reason="no diesel residuals exist (0 diesel fuels qualify; yaml 'excluded')"),
        "dvpe": P("ideal_raoult_psat37.8", f["dvpe"],
                  margin_missing_reason="DVPE targets are units-unconfirmed and excluded from uncertainty_v4.yaml"),
    }
    alternates = {
        "flash": [P("le_chatelier_unifac", f["flash_unifac"], FLU, validated_classes=CRC_CLASSES)],
        "viscosity": [P("refutas_recalibrated", f["visc_recal"], ["viscosity_kin_cSt.crc647.recalibrated_refutas.all_T"],
                        T_valid_C=(-40.0, 20.0), validated_classes=CRC_CLASSES,
                        notes="CRC 647 2-parameter recalibration (reconstructed-composition margin); alternate only"),
                      P("refutas+structured_ridge", f["visc_ridge"], ["viscosity_kin_cSt.crc647.structured_ridge.all_T"],
                        T_valid_C=(-40.0, 20.0), validated_classes=CRC_CLASSES)],
        "freezing": [P("coutinho_library", f["freezing_coutinho"], ["freezing_point_C.crc647.coutinho_library"],
                       validated_classes=CRC_CLASSES),
                     P("ideal_sle_mole_fraction", f["freezing_ideal_mole"],
                       margin_missing_reason="mole-fraction ideal SLE has no yaml margin entry")],
        "nhoc": [P("gcv_group_additivity_compH", f["nhoc_comp_h"],
                   ["heat_of_combustion_net_MJ_kg.crc647.group_additivity_composition_H"], validated_classes=CRC_CLASSES)],
        "smoke_point": [P("ysi_route", f["smoke_ysi"], margin_missing_reason="YSI route not evaluated on real fuels")],
    }
    return dict(props=props, alternates=alternates, classes={n: L.cls(n) for n in L.smiles},
                validated_classes=None)
