"""v5 delta_ws2_freeze_multicomp: multicomponent (n_comp >= 3) liquidus / wax-disappearance data, class-mix gate, model comparison.

Declared readings: results/v5/BASE/delta_ws2_freeze_multicomp/README.md (committed before this script existed, incl. the addendum).
Run: PYTHONPATH=src .venv/Scripts/python.exe scripts/v5/base/delta_ws2_freeze_multicomp.py
Writes rows.csv, summary.csv, systems.csv, readings.json in the delta folder.  Writes nothing to configs/.
"""
from __future__ import annotations

import json
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from rdkit import Chem
from rdkit.Chem import Descriptors

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "src"))
from fuelblend.baselines.flash_unifac import make_unifac_gamma_fn  # noqa: E402
from fuelblend.data.thermoml_sle import system_id  # noqa: E402
from fuelblend.eval import freeze_multicomp as FM  # noqa: E402
from fuelblend.thermo.coutinho_wax import predict_freezing_point_coutinho  # noqa: E402
from fuelblend.uq.conformal import one_sided_margin  # noqa: E402

OUT = ROOT / "results" / "v5" / "BASE" / "delta_ws2_freeze_multicomp"
GATE_SYSTEMS = 8
MODELS = ["svl_ideal", "svl_unifac", "ideal_solid_solution_nalkanes", "coutinho_empirical_depression"]


def can(s):
    return Chem.MolToSmiles(Chem.MolFromSmiles(s))


def ik(s):
    return Chem.MolToInchiKey(Chem.MolFromSmiles(s))


def mw(s):
    return Descriptors.MolWt(Chem.MolFromSmiles(s))


def load_rows() -> pd.DataFrame:
    recs = []
    t = pd.read_csv(ROOT / "data/processed/v5/thermoml_sle_hc.csv")
    t = t[(t.point_kind == "liquidus") & (t.n_comp >= 3)]
    for r in t.itertuples():
        recs.append(dict(source="thermoml", doi=r.doi, sample=" + ".join(json.loads(r.names)), smiles=[can(s) for s in json.loads(r.smiles)],
                         x=[float(v) for v in json.loads(r.x_liquid)], T_K=float(r.T_K), method="liquidus"))
    for f, src, meth in [("data/raw/freezing_lit/oakley2017_liquidus.csv", "oakley2017", "liquidus"),
                         ("data/raw/freezing_lit/ascent065A_freezing.csv", "ascent065A", "freezing_point"),
                         ("data/raw/wdt_quaternary_shariatrad2022/shariatrad2022_wdt.csv", "shariatrad2022", "wax_disappearance")]:
        d = pd.read_csv(ROOT / f)
        d = d[d.kind == "multi"]
        for r in d.itertuples():
            smi = [can(s) for s in json.loads(r.smiles)]
            fr = np.array(json.loads(r.fractions), float)
            if r.fraction_basis == "mass":
                fr = fr / np.array([mw(s) for s in smi])
            elif r.fraction_basis != "mole":
                raise ValueError(r.fraction_basis)
            fr = fr / fr.sum()
            recs.append(dict(source=src, doi=getattr(r, "doi", None), sample=r.sample_id, smiles=smi, x=list(fr),
                             T_K=float(r.value) + 273.15, method=meth))
    w = pd.read_csv(ROOT / "data/processed/v2/wax_dauphin_coutinho.csv")
    w = w[w.n_comp >= 3]
    for r in w.itertuples():
        smi = ["C" * int(n[1:]) for n in json.loads(r.components)]
        assert [ik(s) for s in smi] == json.loads(r.inchikeys), r.system_id
        x = np.array(json.loads(r.mole_fractions), float)
        recs.append(dict(source="dauphin_softx2018", doi="10.1016/j.softx.2018.06.003", sample=r.system_name, smiles=smi, x=list(x / x.sum()),
                         T_K=float(r.wat_measured_K), method="wax_appearance"))
    d = pd.DataFrame(recs)
    d["inchikeys"] = d.smiles.map(lambda ss: [ik(s) for s in ss])
    d["system_id"] = d.inchikeys.map(system_id)
    d["n_comp"] = d.smiles.map(len)
    d["class_mix"] = d.smiles.map(FM.class_mix)
    return d


def dedup(d: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Drop a non-ThermoML row that matches a ThermoML row of the same system within |dT| < 0.5 K and max |dx| < 0.005."""
    th = d[d.source == "thermoml"]
    drop = []
    for i, r in d[d.source != "thermoml"].iterrows():
        for _, q in th[th.system_id == r.system_id].iterrows():
            xr = dict(zip(r.inchikeys, r.x))
            if abs(r.T_K - q.T_K) < 0.5 and max(abs(xr[k] - v) for k, v in zip(q.inchikeys, q.x)) < 0.005:
                drop.append(i)
                break
    return d.drop(index=drop).reset_index(drop=True), d.loc[drop]


def constants():
    v = pd.read_csv(ROOT / "data/raw/sle_constants/tm_dhfus_verified.csv")
    t = pd.read_csv(ROOT / "data/processed/v5/thermoml_fusion_pure_hc_by_compound.csv")
    tier, ver = {}, {}
    for r in t.itertuples():
        if np.isfinite(r.tm_K) and np.isfinite(r.dhfus_kJ_mol) and r.dhfus_kJ_mol > 0:
            tier[r.inchikey] = (float(r.tm_K), float(r.dhfus_kJ_mol) * 1000.0, "thermoml_measured")
    for r in v.itertuples():
        if np.isfinite(r.tm_K) and np.isfinite(r.svl_dhfus_kJ_mol) and r.svl_dhfus_kJ_mol > 0:
            tier[r.inchikey] = ver[r.inchikey] = (float(r.tm_K), float(r.svl_dhfus_kJ_mol) * 1000.0, "verified")
    return tier, ver


def predict_row(r, cmap):
    consts = [cmap.get(k) for k in r.inchikeys]
    if any(c is None for c in consts):
        return {m: np.nan for m in MODELS}, "no_constant", [s for s, c in zip(r.smiles, consts) if c is None], False
    cc = [(c[0], c[1]) for c in consts]
    out = {"svl_ideal": FM.svl_ideal_C(r.x, cc)}
    gfn, fb = make_unifac_gamma_fn(r.smiles)
    out["svl_unifac"] = FM.svl_activity_C(r.x, cc, gfn)
    out["ideal_solid_solution_nalkanes"] = FM.solid_solution_nalkanes_C(r.x, r.smiles, cc)
    tmap = {}
    for k, s, c in zip(r.inchikeys, r.smiles, cc):
        m = Chem.MolFromSmiles(s)
        cls = "n_paraffin" if FM.is_n_alkane(s) else "other"
        tmap[k] = dict(tm_K=c[0], dh_fus_J_mol=c[1], n_carbon=float(sum(a.GetSymbol() == "C" for a in m.GetAtoms())), cls=cls)
    v, _, _ = predict_freezing_point_coutinho(r.x, r.inchikeys, tmap)
    out["coutinho_empirical_depression"] = float(v) if v is not None and v > -900 else np.nan
    return out, "ok", [], bool(fb)


def score(d, cmap, label):
    rows = []
    for r in d.itertuples():
        p, status, missing, fb = predict_row(r, cmap)
        for m in MODELS:
            rows.append(dict(constants=label, source=r.source, doi=r.doi, sample=r.sample, system_id=r.system_id, class_mix=r.class_mix,
                             n_comp=r.n_comp, method=r.method, measured_C=r.T_K - 273.15, model=m, pred_C=p[m], status=status,
                             missing=json.dumps(missing), unifac_fallback=fb))
    return pd.DataFrame(rows)


def summarize(rows, margin_binary):
    out = []
    for (const, model), g0 in rows.groupby(["constants", "model"]):
        for scope, g1 in [("all", g0), ("excluding_wax_appearance", g0[g0.method != "wax_appearance"])]:
            for cm, g in list(g1.groupby("class_mix")) + [("pooled", g1)]:
                ok = g[(g.status == "ok") & np.isfinite(g.pred_C)]
                e = ok.pred_C - ok.measured_C
                rec = dict(constants=const, model=model, scope=scope, class_mix=cm, n_rows=len(g), n_scored=len(ok),
                           n_systems_in_M=g.system_id.nunique(), n_systems_scored=ok.system_id.nunique(),
                           mae=e.abs().mean(), bias=e.mean(), frac_under=float((e < 0).mean()) if len(e) else np.nan,
                           max_under=float(-e.min()) if len(e) else np.nan,
                           covered_by_binary_margin=int(((ok.pred_C + margin_binary) >= ok.measured_C).sum()))
                rec["gate_passed"] = bool(rec["n_systems_scored"] >= GATE_SYSTEMS) and cm != "pooled"
                rec["margin_alpha10_under"] = (float(one_sided_margin(e.to_numpy(), ok.system_id.to_numpy(), alpha=0.10, direction="under"))
                                               if rec["gate_passed"] and model == "svl_ideal" else np.nan)
                out.append(rec)
    return pd.DataFrame(out)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    d, dropped = dedup(load_rows())
    tier, ver = constants()
    rows = pd.concat([score(d, tier, "verified_then_thermoml"), score(d, ver, "verified_only")], ignore_index=True)
    ym = yaml.safe_load((ROOT / "configs/uncertainty_v5.yaml").read_text())
    mb = float(ym["margins"]["freezing_point_C.lab_liquidus.svl_verified"]["value"])   # read only
    summ = summarize(rows, mb)
    systems = (d.groupby(["system_id", "class_mix", "source", "method"]).agg(n_comp=("n_comp", "first"), n_rows=("T_K", "size"),
               sample=("sample", "first"), T_min_C=("T_K", lambda s: s.min() - 273.15), T_max_C=("T_K", lambda s: s.max() - 273.15))
               .reset_index())
    sc = rows[(rows.constants == "verified_then_thermoml") & (rows.model == "svl_ideal")]
    systems["scorable"] = systems.system_id.map(lambda s: bool((sc[sc.system_id == s].status == "ok").all()))
    systems["missing_constants"] = systems.system_id.map(lambda s: sorted(set(sum((json.loads(m) for m in sc[sc.system_id == s].missing), []))))
    rows.to_csv(OUT / "rows.csv", index=False)
    summ.to_csv(OUT / "summary.csv", index=False)
    systems.to_csv(OUT / "systems.csv", index=False)
    main_s = summ[(summ.constants == "verified_then_thermoml") & (summ.model == "svl_ideal") & (summ.scope == "all")]
    readings = {
        "dedup_dropped": dropped[["source", "sample", "T_K"]].to_dict("records"),
        "binary_margin_C": mb,
        "R-W2-1": {r.class_mix: dict(n_systems_in_M=int(r.n_systems_in_M), n_systems_scored=int(r.n_systems_scored)) for r in main_s.itertuples()},
        "R-W2-2": {r.class_mix: dict(gate_systems=GATE_SYSTEMS, passed=bool(r.gate_passed), margin=r.margin_alpha10_under)
                   for r in main_s.itertuples() if r.class_mix != "pooled"},
        "R-W2-3": summ[(summ.constants == "verified_then_thermoml") & (summ.scope == "all")][
            ["model", "class_mix", "n_scored", "n_systems_scored", "mae", "bias", "frac_under", "max_under"]].to_dict("records"),
        "R-W2-4": {r.class_mix: dict(n_scored=int(r.n_scored), covered=int(r.covered_by_binary_margin)) for r in main_s.itertuples()},
    }
    (OUT / "readings.json").write_text(json.dumps(readings, indent=2, default=float), encoding="utf-8")
    pd.set_option("display.width", 250, "display.max_columns", 30)
    print(systems.to_string())
    print(summ[summ.scope == "all"].round(2).to_string())
    print(json.dumps({k: readings[k] for k in ("dedup_dropped", "R-W2-1", "R-W2-2", "R-W2-4")}, indent=2, default=float))


if __name__ == "__main__":
    warnings.simplefilter("ignore")
    main()
