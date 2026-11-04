"""Flash-point physics baseline on a synthetic unified dataset built from REAL refit Antoine constants. Run: python tests/test_baselines_flash.py"""
import json, shutil, subprocess, sys, tempfile
from pathlib import Path
import numpy as np, pandas as pd
from rdkit import Chem
HERE = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(HERE))
import baselines_flash as bf
W = Path(tempfile.mkdtemp()); [shutil.copy(f, W / f.name) for f in HERE.glob("*.py")]
(W / "data/unified").mkdir(parents=True)
key = lambda s: Chem.MolToInchiKey(Chem.MolFromSmiles(s))
# log10(P/kPa) = A - B/(T+C): constants refitted from measured ThermoML vapour pressures (user's data)
ANT = {"C" * 10: (6.199, 1580.177, -70.46, 276.8, 475.4, 0.75, 10), "C" * 14: (6.4, 1900.0, -88.0, 283.0, 498.0, 0.50, 14)}
# derive consistent decane/tetradecane constants: decane from the user's refit; tetradecane approximated for the test only
rows = []
for s, (A, B, C, lo, hi, lfl, n) in ANT.items():
    rows.append(dict(inchikey=key(s), smiles=s, is_hc=True, cls="n_paraffin", n_carbon=n, antoine_A=A, antoine_B=B, antoine_C=C,
                     antoine_Tlo=lo, antoine_Thi=hi, lfl_vol_pct=lfl, lfl_source="measured: Zabetakis627_Table2"))
rows.append(dict(inchikey=key("Cc1ccccc1"), smiles="Cc1ccccc1", is_hc=True, cls="aromatic", n_carbon=7, antoine_A=6.08, antoine_B=1346.0, antoine_C=-53.5,
                 antoine_Tlo=280, antoine_Thi=410, lfl_vol_pct=1.0, lfl_source="estimated: 0.55*Cst (Jones rule)"))
rows.append(dict(inchikey=key("CCCCCCCCCCCCCCCC"), smiles="C" * 16, is_hc=True, cls="n_paraffin", n_carbon=16, antoine_A=np.nan, antoine_B=np.nan,
                 antoine_C=np.nan, antoine_Tlo=np.nan, antoine_Thi=np.nan, lfl_vol_pct=0.43, lfl_source="measured"))
pd.DataFrame(rows).to_csv(W / "data/unified/compounds.csv", index=False)
k10, k14, kt, k16 = key("C" * 10), key("C" * 14), key("Cc1ccccc1"), key("C" * 16)
def rec(sys_, keys, x, val):
    o = np.argsort(keys); keys = [keys[i] for i in o]; x = [x[i] for i in o]
    return dict(property="flash_point_C", is_hc=True, kind="pure" if len(keys) == 1 else "binary", n_comp=len(keys), system_id=sys_,
                comp_inchikeys=json.dumps(keys), x_mole=json.dumps(x), value=val, n_obs=1, conflict=False, T_K=np.nan)
pd.DataFrame([rec("a", [k10], [1.0], 46.0), rec("b", [k14], [1.0], 100.0), rec("c", [k10, k14], [0.5, 0.5], 62.0),
              rec("d", [kt], [1.0], 4.0), rec("e", [k16], [1.0], 135.0)]).to_csv(W / "data/unified/measurements_dedup.csv", index=False)

# unit-level: the solver returns a root of H(T) = 1, and the decane value matches the real-data sanity check (47.4 C)
A, B, C = ANT["C" * 10][:3]
p, why = bf.predict_flash_C([1.0], [A], [B], [C], [0.75]); assert why == "" and abs(p - 47.4) < 0.3, p
assert abs(bf.hazard(p + 273.15, [1.0], [A], [B], [C], np.array([0.0075]), False)) < 1e-6
p_lt, _ = bf.predict_flash_C([1.0], [A], [B], [C], [0.75], lfl_temp=True); assert p_lt < p            # LFL falls with T -> flash point lower

r = subprocess.run([sys.executable, "baselines_flash.py"], cwd=W, capture_output=True, text=True); assert r.returncode == 0, r.stdout + r.stderr
o = pd.read_csv(W / "data/unified/baseline_flash.csv").set_index("system_id")
assert o.loc["e", "status"].startswith("missing antoine")                                         # hexadecane without Antoine: reported, not guessed
d, t14, mix = o.loc["a", "predicted_C"], o.loc["b", "predicted_C"], o.loc["c", "predicted_C"]
assert abs(d - 47.4) < 0.3 and d < mix < t14                                                       # ideal mixture lies between its pure components
assert abs(o.loc["a", "err_C"] - (d - 46.0)) < 1e-9 and bool(o.loc["a", "antoine_in_range"]) and bool(o.loc["a", "lfl_all_measured"])
assert o.loc["d", "lfl_all_measured"] == False and bool(o.loc["d", "has_aromatic"])               # estimated LFL is flagged
rep = (W / "data/unified/BASELINE_FLASH_REPORT.md").read_text()
assert "Overall" in rep and "missing antoine_A" in rep and "ideal-mixing" in rep.lower()
print("ALL FLASH BASELINE TESTS PASSED")
