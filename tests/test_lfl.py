"""fill_lfl.py: priority (user > Zabetakis > estimate), estimate maths, pickup by build_unified. Run: python tests/test_lfl.py"""
import shutil, subprocess, sys, tempfile
from pathlib import Path
import pandas as pd
from rdkit import Chem
HERE = Path(__file__).resolve().parent.parent
W = Path(tempfile.mkdtemp()); [shutil.copy(f, W / f.name) for f in HERE.glob("*.py")]
key = lambda s: Chem.MolToInchiKey(Chem.MolFromSmiles(s))
smis = {"n-hexane": "CCCCCC", "n-decane": "CCCCCCCCCC", "n-hexadecane": "C" * 16, "isooctane": "CC(C)CC(C)(C)C", "toluene": "Cc1ccccc1", "ethanol": "CCO"}
(W / "data/unified").mkdir(parents=True)
pd.DataFrame([dict(inchikey=key(s), smiles=s, is_hc=(n != "ethanol"), cls="x") for n, s in smis.items()]).to_csv(W / "data/unified/compounds.csv", index=False)
(W / "data/raw/lfl").mkdir(parents=True)
pd.DataFrame([dict(smiles="Cc1ccccc1", lfl_vol_pct=1.1, T_C=25, source="my handbook")]).to_csv(W / "data/raw/lfl/lfl_measured.csv", index=False)
r = subprocess.run([sys.executable, "fill_lfl.py"], cwd=W, capture_output=True, text=True); assert r.returncode == 0, r.stdout + r.stderr
d = pd.read_csv(W / "data/raw/lfl/lfl.csv").set_index("inchikey")
assert len(d) == 5 and key("CCO") not in d.index                                  # non-hydrocarbon skipped
assert d.loc[key("CCCCCC"), "lfl_vol_pct"] == 1.2 and d.loc[key("CCCCCC"), "source"].startswith("measured: Zabetakis")
assert d.loc[key("C" * 10), "lfl_vol_pct"] == 0.75 and d.loc[key("C" * 16), "lfl_vol_pct"] == 0.43
assert d.loc[key("Cc1ccccc1"), "lfl_vol_pct"] == 1.1 and d.loc[key("Cc1ccccc1"), "source"] == "measured: my handbook"   # user override wins
est = d.loc[key("CC(C)CC(C)(C)C")]                                                 # isooctane: C8H18 -> Cst 1.648 -> 0.55*Cst = 0.907
assert est.source.startswith("estimated") and abs(est.lfl_vol_pct - 0.907) < 0.005
# the estimate rule reproduces the measured n-alkanes (hexane..hexadecane) within 10 %
sys.path.insert(0, str(HERE)); import fill_lfl as fl
for n in range(6, 17):
    e = 0.55 * fl.cst_vol_pct(Chem.MolFromSmiles("C" * n)); m = fl.ZABETAKIS_T2[n]
    assert abs(e - m) / m < 0.10, (n, e, m)
print("ALL LFL TESTS PASSED")
