import json
import pandas as pd

df = pd.read_csv("data/processed/model_table.csv")
comp = pd.read_csv("data/unified/compounds.csv").set_index("inchikey")

mix = df[df["kind"].isin(["binary", "multi"])].drop_duplicates("system_id")
multi = mix[mix["n_comp"] >= 3]
print("Total multi systems:", len(multi))
for idx, r in multi.iterrows():
    sys_id = r["system_id"]
    keys = json.loads(r["comp_inchikeys"])
    c_nums = [comp.loc[k, "n_carbon"] for k in keys if k in comp.index]
    min_c = min(c_nums) if c_nums else 0
    max_c = max(c_nums) if c_nums else 0
    mean_c = sum(c_nums) / len(c_nums) if c_nums else 0
    props = df[df["system_id"] == sys_id]["property"].unique().tolist()
    print(f"{sys_id} n={r['n_comp']:2d} C:[{min_c:2d}-{max_c:2d}, avg={mean_c:4.1f}] props={props}")
