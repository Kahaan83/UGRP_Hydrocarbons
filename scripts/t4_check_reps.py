"""T4 prep: find surrogate representatives with Antoine cover + LFL per class."""
import pandas as pd

comp = pd.read_csv("data/unified/compounds.csv", low_memory=False)
with_antoine = comp[comp["antoine_A"].notna() & comp["lfl_vol_pct"].notna() & (comp["is_hc"] == True)]
print("classes with antoine+lfl, ordered by n_carbon:")
for cls in ["cycloalkane", "aromatic", "n_paraffin", "iso_paraffin"]:
    sub = with_antoine[with_antoine["cls"] == cls].sort_values("n_carbon")
    cols = ["smiles", "n_carbon", "mw", "antoine_Tlo", "antoine_Thi", "antoine_source", "lfl_vol_pct"]
    show = sub[cols].head(12).to_string(index=False) if len(sub) else "(none)"
    print(f"\n--- {cls} (n={len(sub)})\n{show}")
