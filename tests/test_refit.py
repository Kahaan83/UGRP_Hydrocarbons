"""Antoine refit on synthetic vapour pressures with known truth. Run: python tests/test_refit.py"""
import json, shutil, subprocess, sys, tempfile
from pathlib import Path
import numpy as np, pandas as pd
from rdkit import Chem
HERE = Path(__file__).resolve().parent.parent
W = Path(tempfile.mkdtemp()); [shutil.copy(f, W / f.name) for f in HERE.glob("*.py")]
(W / "data/unified").mkdir(parents=True); (W / "data/processed").mkdir(parents=True)
rng = np.random.default_rng(1)
key = lambda s: Chem.MolToInchiKey(Chem.MolFromSmiles(s))
true = {"CCCCCCCCCC": (4.07857 + 2, 1501.268, -78.67),      # decane-like (kPa form)
        "CCCCCCCCCCCC": (4.10 + 2, 1625.0, -92.0)}          # dodecane-like
rows = []
for smi, (A, B, C) in true.items():
    for T in np.linspace(300, 440, 25):
        P = 10 ** (A - B / (T + C)) * 10 ** rng.normal(0, 0.004)          # ~1 % noise
        rows.append(dict(source="thermoml", property="vapor_pressure_kPa", value=P, T_K=T, n_comp=1, comp_inchikeys=json.dumps([key(smi)]), ref="r%d" % int(T)))
rows.append(dict(source="thermoml", property="vapor_pressure_kPa", value=500.0, T_K=300.0, n_comp=1, comp_inchikeys=json.dumps([key(smi)]), ref="bad"))  # gross outlier (also P>300 -> filtered)
for T in (310, 320, 330):                                                    # hexadecane: too little data
    rows.append(dict(source="thermoml", property="vapor_pressure_kPa", value=0.01, T_K=T, n_comp=1, comp_inchikeys=json.dumps([key("CCCCCCCCCCCCCCCC")]), ref="x"))
for T, P in zip(np.linspace(300, 360, 8), np.linspace(50, 10, 8)):       # tetradecane with BAD data: pressure falls as T rises (used to crash the run)
    rows.append(dict(source="thermoml", property="vapor_pressure_kPa", value=P, T_K=T, n_comp=1, comp_inchikeys=json.dumps([key("CCCCCCCCCCCCCC")]), ref="bad%d" % int(T)))
pd.DataFrame(rows).to_csv(W / "data/unified/measurements.csv", index=False)
pd.DataFrame([dict(inchikey=key(s), smiles=s) for s in list(true) + ["CCCCCCCCCCCCCCCC", "CCCCCCCCCCCCCC"]]).to_csv(W / "data/unified/compounds.csv", index=False)
# WebBook decane set (bar) valid only 367-448 K: refit must NOT disagree with it in-range
pd.DataFrame([dict(name="n-decane", smiles="CCCCCCCCCC", T_lo=367.63, T_hi=448.27, A=4.07857, B=1501.268, C=-78.67)]).to_csv(W / "data/processed/antoine_best.csv", index=False)

r = subprocess.run([sys.executable, "refit_antoine.py"], cwd=W, capture_output=True, text=True); assert r.returncode == 0, r.stdout + r.stderr
f = pd.read_csv(W / "data/processed/antoine_refit.csv").set_index("inchikey")
for smi, (A, B, C) in true.items():
    q = f.loc[key(smi)]; assert q.status == "ok", q
    for T in (313.15, 373.15):                                                # the window the WebBook sets could not cover
        err = abs(np.log10(10 ** (q.A - q.B / (T + q.C))) - (A - B / (T + C)))
        assert err < 0.01, (smi, T, err)                                      # < 2.5 % in pressure
    assert bool(q.covers_window) is False or q["T_max"] >= 373.15            # measured span 300-440 K -> covers
assert f.loc[key("CCCCCCCCCC"), "anchored"] and f.loc[key("CCCCCCCCCC"), "webbook_maxdiff"] < 0.02
assert "skipped" in f.loc[key("CCCCCCCCCCCCCCCC"), "status"]
assert f.loc[key("CCCCCCCCCCCCCC"), "status"].startswith("failed")          # bad compound is reported, the run continues
assert f.loc[key("CCCCCCCCCC"), "status"] == "ok"                           # and the good ones are unaffected
print("ALL REFIT TESTS PASSED")
