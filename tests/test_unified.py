"""Integration test for build_unified.py. Run: python tests/test_unified.py"""
import json, shutil, subprocess, sys, tempfile
from pathlib import Path
import numpy as np, pandas as pd
from rdkit import Chem

HERE = Path(__file__).resolve().parent.parent
W = Path(tempfile.mkdtemp())
for f in HERE.glob("*.py"): shutil.copy(f, W / f.name)
shutil.copy(HERE / "compounds.csv", W / "compounds.csv")
inchi = lambda smi: Chem.MolToInchi(Chem.MolFromSmiles(smi))
key = lambda smi: Chem.MolToInchiKey(Chem.MolFromSmiles(smi))
HEX, HEP, DEC, DOD, TOL, ETH = (inchi(s) for s in ("CCCCCC", "CCCCCCC", "CCCCCCCCCC", "CCCCCCCCCCCC", "Cc1ccccc1", "CCO"))

# ---- Leenhouts-like tables (real column schema from the published inventory)
lp = W / "data/interim/leenhouts_benchmark/datasets/datasets"
(lp / "flashpoint/processed_data").mkdir(parents=True); (lp / "viscosity/processed_data").mkdir(parents=True)
pd.DataFrame({"value": [-22.0, 71.0, 13.0], "unit": ["Celcius"] * 3, "component_inchi_0": [HEX, DEC, ETH]}) \
  .to_csv(lp / "flashpoint/processed_data/flashpoint_pure_exp.csv", index=False)
# ternary with UNSORTED order: dec 0.251, dod 0.210, hex 0.539 ; plus an ethanol mixture that must be dropped
pd.DataFrame({"value": [71.85, 20.0], "unit": ["Celcius"] * 2, "fraction_type": ["molar"] * 2,
              "component_inchi_0": [DEC, HEX], "component_inchi_1": [DOD, ETH], "component_inchi_2": [HEX, None],
              "component_inchi_3": [None, None], "fraction_0": [0.251, 0.5], "fraction_1": [0.210, 0.5],
              "fraction_2": [0.539, np.nan], "fraction_3": [np.nan, np.nan]}) \
  .to_csv(lp / "flashpoint/processed_data/flashpoint_mix_exp.csv", index=False)

nist_pure = pd.read_csv(HERE / "data/processed/viscosity_pure_hc.csv")
nist_mix = pd.read_csv(HERE / "data/processed/viscosity_mix_hc.csv")
(W / "data/processed").mkdir(parents=True)
shutil.copy(HERE / "data/processed/viscosity_pure_hc.csv", W / "data/processed/")
shutil.copy(HERE / "data/processed/viscosity_mix_hc.csv", W / "data/processed/")
def single_rows(df, keycols):
    d = df.copy(); d["Tk"] = (d.T_K * 2).round() / 2
    d = d.drop_duplicates(keycols + ["Tk", "Ref_ID", "eta_mPas"])
    n = d.groupby(keycols + ["Tk"]).Ref_ID.transform("size")
    return d[n == 1]
sp = single_rows(nist_pure, ["smiles"]); sm_ = single_rows(nist_mix, ["smiles_list", "x_mole"])
p0, p1 = sp.iloc[0], sp.iloc[1]          # REAL NIST rows seen exactly once: p0 re-reported identically, p1 re-reported 2x too high
m0 = sm_.iloc[0]
sm = json.loads(m0.smiles_list); xm = json.loads(m0.x_mole)
vis_rows = [dict(value=np.log(p0.eta_mPas), unit="ln(mPa*s)", temperature=p0.T_K, component_inchi_0=inchi(p0.smiles)),
            dict(value=np.log(p1.eta_mPas * 2.0), unit="ln(mPa*s)", temperature=p1.T_K, component_inchi_0=inchi(p1.smiles))]
pd.DataFrame(vis_rows).to_csv(lp / "viscosity/processed_data/viscosity_pure_exp.csv", index=False)
pd.DataFrame([dict(value=np.log(m0.eta_mPas), unit="ln(mPa*s)", fraction_type="molar", temperature=m0.T_K,
                   component_inchi_0=inchi(sm[0]), component_inchi_1=inchi(sm[1]), fraction_0=xm[0], fraction_1=xm[1])]) \
  .to_csv(lp / "viscosity/processed_data/viscosity_mix_exp.csv", index=False)
# make a CONFLICT: same NIST pure point at +10 K was given with a value 2x too large in Leenhouts but NIST has 1x there?
# (only conflicts if NIST has that T; instead force one explicitly below)

# ---- ThermoML-like (mass-fraction mixture, kinematic + dynamic + vapour pressure)
pd.DataFrame([
    dict(file="a.xml", prop="Viscosity, Pa*s", value=0.00031, T_K=298.15, P_kPa=101.3, n_comp=1, inchis=json.dumps([HEX]), names="[]", comp_type="pure", comp_values=json.dumps([1.0]), prop_phase="Liquid"),
    dict(file="a.xml", prop="Mass density, kg/m3", value=654.8, T_K=298.15, P_kPa=101.3, n_comp=1, inchis=json.dumps([HEX]), names="[]", comp_type="pure", comp_values=json.dumps([1.0]), prop_phase="Liquid"),
    dict(file="b.xml", prop="Vapor or sublimation pressure, kPa", value=0.05, T_K=313.15, P_kPa=np.nan, n_comp=1, inchis=json.dumps([DOD]), names="[]", comp_type="pure", comp_values=json.dumps([1.0]), prop_phase="Liquid"),
    dict(file="c.xml", prop="Mass density, kg/m3", value=780.0, T_K=298.15, P_kPa=101.3, n_comp=2, inchis=json.dumps([HEX, TOL]), names="[]", comp_type="mass fraction", comp_values=json.dumps([0.5, 0.5]), prop_phase="Liquid"),
    dict(file="d.xml", prop="Mass density, kg/m3", value=700.0, T_K=298.15, P_kPa=101.3, n_comp=2, inchis=json.dumps([HEX, TOL]), names="[]", comp_type=None, comp_values=json.dumps([None, None]), prop_phase="Liquid"),
]).to_csv(W / "data/processed/thermoml_hc.csv", index=False)

# same hexane density at the same T but 50 MPa: must stay a SEPARATE record, not merge with the atmospheric one
extra = pd.DataFrame([dict(file="a.xml", prop="Mass density, kg/m3", value=700.0, T_K=298.15, P_kPa=50000.0, n_comp=1, inchis=json.dumps([HEX]),
                           names="[]", comp_type="pure", comp_values=json.dumps([1.0]), prop_phase="Liquid"),
                      dict(file="a.xml", prop="Excess viscosity, Pa*s", value=-0.0001, T_K=298.15, P_kPa=101.0, n_comp=1, inchis=json.dumps([HEX]),
                           names="[]", comp_type="pure", comp_values=json.dumps([1.0]), prop_phase="Liquid")])
extra = pd.concat([extra, pd.DataFrame([
    dict(file="z.xml", prop="Mass density, kg/m3", value=25.3, T_K=777.0, P_kPa=500.0, n_comp=1, inchis=json.dumps([HEX]), names="[]", comp_type="pure", comp_values=json.dumps([1.0]), prop_phase="Liquid"),
    dict(file="z.xml", prop="Viscosity, Pa*s", value=-0.0001, T_K=298.15, P_kPa=101.0, n_comp=1, inchis=json.dumps([HEX]), names="[]", comp_type="pure", comp_values=json.dumps([1.0]), prop_phase="Liquid"),
    dict(file="z.xml", prop="Vapor or sublimation pressure, kPa", value=1.4e6, T_K=338.0, P_kPa=np.nan, n_comp=1, inchis=json.dumps([HEX]), names="[]", comp_type="pure", comp_values=json.dumps([1.0]), prop_phase="Liquid")])])
pd.concat([pd.read_csv(W / "data/processed/thermoml_hc.csv"), extra]).to_csv(W / "data/processed/thermoml_hc.csv", index=False)

# ---- AJFTD-like real fuel rows, Antoine, LFL
pd.DataFrame([dict(fuel_id="POSF10325", property="flash_point_C", value=43.0, unit_std="C", T_K=np.nan, method="D93"),
              dict(fuel_id="POSF10325", property="viscosity_kin_cSt", value=7.9, unit_std="cSt", T_K=253.15, method="D445")]) \
  .to_csv(W / "data/processed/ajftd_fuels.csv", index=False)
pd.DataFrame([dict(name="n-hexane", cas="110-54-3", smiles="CCCCCC", T_lo=286.0, T_hi=380.0, A=3.45, B=1044.0, C=-53.9, coverage=1.0, covers_window=True),
              dict(name="n-decane", cas="124-18-5", smiles="CCCCCCCCCC", T_lo=367.0, T_hi=448.0, A=4.07, B=1501.0, C=-80.0, coverage=0.09, covers_window=False)]) \
  .to_csv(W / "data/processed/antoine_best.csv", index=False)

def run():
    r = subprocess.run([sys.executable, "build_unified.py"], cwd=W, capture_output=True, text=True)
    assert r.returncode == 0, r.stdout + r.stderr
    return r.stdout

out = run()
U = W / "data/unified"
m = pd.read_csv(U / "measurements.csv"); d = pd.read_csv(U / "measurements_dedup.csv"); c = pd.read_csv(U / "compounds.csv")

# non-hydrocarbon rows dropped
assert not m.comp_inchikeys.str.contains(key("CCO")).any()
assert key("CCO") not in set(c.inchikey)
# ternary: order sorted by InChIKey, fractions stay attached to the right molecule
t = m[(m.source == "leenhouts") & (m.n_comp == 3)].iloc[0]
ks, xs = json.loads(t.comp_inchikeys), json.loads(t.x_mole)
assert ks == sorted(ks) and abs(sum(xs) - 1) < 1e-6
assert abs(xs[ks.index(key("CCCCCCCCCC"))] - 0.251) < 1e-6 and abs(xs[ks.index(key("CCCCCC"))] - 0.539) < 1e-6
w = json.loads(t.w_mass); Mdec, Mdod, Mhex = 142.286, 170.34, 86.178
exp_w = np.array([0.251 * Mdec, 0.210 * Mdod, 0.539 * Mhex]); exp_w /= exp_w.sum()
assert abs(w[ks.index(key("CCCCCCCCCC"))] - exp_w[0]) < 2e-3
# units: ln(mPa*s)->mPa*s ; Pa*s->mPa*s ; flash stays in C
lv = m[(m.source == "leenhouts") & (m.property == "viscosity_dyn_mPas") & (m.n_comp == 1)].sort_values("value").iloc[0]
if abs(lv.value - p0.eta_mPas) > 1e-6: lv = m[(m.source == "leenhouts") & (m.n_comp == 1) & (m.property == "viscosity_dyn_mPas")].iloc[0]
assert any(abs(v - p0.eta_mPas) / p0.eta_mPas < 1e-6 for v in m[(m.source == "leenhouts") & (m.n_comp == 1) & (m.property == "viscosity_dyn_mPas")].value)
tv = m[(m.source == "thermoml") & (m.property == "viscosity_dyn_mPas")].iloc[0]; assert abs(tv.value - 0.31) < 1e-9
# ThermoML mass-fraction mixture is stored with mass fractions preserved, mole fractions derived
tm = m[(m.source == "thermoml") & (m.n_comp == 2)]; assert len(tm) == 1            # the row with no composition was skipped
wm = json.loads(tm.iloc[0].w_mass); assert np.allclose(wm, [0.5, 0.5], atol=1e-4)
assert not np.allclose(json.loads(tm.iloc[0].x_mole), [0.5, 0.5], atol=1e-2)       # mole != mass for hexane/toluene
# dedupe: the real NIST pure point reported by both sources merges to one record with 2 sources and no conflict
g = d[(d.property == "viscosity_dyn_mPas") & (d.kind == "pure") & d.sources.str.contains("leenhouts") & d.sources.str.contains("nist_tde")]
assert len(g) == 2, g[["T_K", "value", "n_obs", "conflict"]]
assert (g.n_sources == 2).all()
assert g.conflict.sum() == 1                      # identical re-report merges cleanly; the 2x one is flagged
good = g[~g.conflict].iloc[0]; assert abs(good.value - p0.eta_mPas) / p0.eta_mPas < 1e-6
gm = d[(d.kind == "binary") & d.sources.str.contains("leenhouts") & d.sources.str.contains("nist_tde")]
assert len(gm) >= 1 and not gm.conflict.any()
assert len(d) < len(m)
# pressure kept apart; excess viscosity never becomes a viscosity
hx = d[(d.property == "density_kg_m3") & (d.kind == "pure") & (d.system_id == m[m.comp_inchikeys.str.contains(key("CCCCCC")) & (m.n_comp == 1)].system_id.iloc[0])]
assert sorted(hx.P_kPa.round(0).unique())[-1] == 50000 and len(hx) == 2, hx[["T_K", "P_kPa", "value"]]
assert (m[m.property == "viscosity_dyn_mPas"].value > 0).all()
# domain filter: gas-like 'liquid' density, negative viscosity and a GPa 'vapour pressure' are removed and reported
assert m.value[m.property == "density_kg_m3"].min() > 400 and m.value[m.property == "vapor_pressure_kPa"].max() <= 1000
ood = pd.read_csv(U / "measurements_out_of_domain.csv"); oo = ood[ood.source == "thermoml"]          # (the real NIST rows it also drops are genuine outliers: 8.8e5 mPa.s and ~4e-4 mPa.s)
assert len(oo) == 3 and set(oo.property) == {"density_kg_m3", "viscosity_dyn_mPas", "vapor_pressure_kPa"}
assert "domain filter: dropped 1 density_kg_m3" in out
# real fuels carried through without molecules
rf = d[d.kind == "real_fuel"]; assert len(rf) == 2 and (rf.n_comp == 0).all()
# physics join: hexane Antoine covers window, decane does not; no LFL yet -> template written
ch = c.set_index("inchikey")
assert bool(ch.loc[key("CCCCCC"), "antoine_covers_window"]) and not bool(ch.loc[key("CCCCCCCCCC"), "antoine_covers_window"])
assert (W / "data/raw/lfl/lfl_template.csv").exists()
assert abs(ch.loc[key("CCCCCC"), "mw"] - 86.178) < 0.01
assert abs(ch.loc[key("CCCCCC"), "antoine_A"] - 5.45) < 1e-9          # 3.45 (bar) + 2 -> kPa
# a filled LFL file is picked up and the report counts runnable flash-point records
tpl = pd.read_csv(W / "data/raw/lfl/lfl_template.csv"); assert key("CCCCCC") in set(tpl.inchikey)
pd.DataFrame([dict(smiles="CCCCCC", lfl_vol_pct=1.2, T_C=25, source="Zabetakis")]).to_csv(W / "data/raw/lfl/lfl.csv", index=False)
out2 = run(); c2 = pd.read_csv(U / "compounds.csv").set_index("inchikey")
assert abs(c2.loc[key("CCCCCC"), "lfl_vol_pct"] - 1.2) < 1e-9 and c2.loc[key("CCCCCC"), "lfl_source"] == "Zabetakis"
assert "can run" in out2 and "UNIFIED" in out2.upper()
print("ALL UNIFIED TESTS PASSED")

# ---- an OLD-format thermoml_hc.csv (no prop_phase) must be refused loudly
old = pd.read_csv(W / "data/processed/thermoml_hc.csv").drop(columns=["prop_phase"]); old.to_csv(W / "data/processed/thermoml_hc.csv", index=False)
out3 = run(); assert "OLD extractor" in out3 and "thermoml: no rows" in out3
print("OLD-FORMAT REFUSAL OK")
