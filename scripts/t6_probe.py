import pandas as pd, json

mt = pd.read_csv("data/processed/model_table.csv", low_memory=False)
g = mt[(mt["property"].isin(["mon", "dcn"])) & (mt["kind"] != "pure")]
r = g.iloc[0]
print("comp_inchikeys:", type(r.comp_inchikeys).__name__, repr(r.comp_inchikeys)[:120])
print("w_mass:", type(r.w_mass).__name__, repr(r.w_mass)[:120])
try:
    ik = json.loads(r.comp_inchikeys) if isinstance(r.comp_inchikeys, str) else r.comp_inchikeys
    print("parsed ik:", ik)
except Exception as e:
    print("parse err:", e)
try:
    w = json.loads(r.w_mass) if isinstance(r.w_mass, str) else r.w_mass
    print("parsed w:", w)
except Exception as e:
    print("parse err:", e)