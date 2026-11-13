"""R-SM-3: every non-smoke record of predict_fuel is unchanged by the smoke margin.
'Before' = the same code with the smoke key removed from the v5 margins cache (smoke then stays uncertified, the previous behaviour);
'after' = the shipped yaml.  Blends: the 60-component demo (seed 0) and a 5-compound known blend, each at 233.15, 288.15, 313.15 K."""
import sys
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT / "src"))
from fuelblend.api import predict_route as pr  # noqa: E402

ctx = pr.get_ctx()
P = pr.predicates(ctx)
cand = [i for i in range(len(ctx.reg.ik)) if (P["HAS_GRID"](i) or P["HAS_FLASH"](i)) and pr.subkind(ctx.reg.smiles[i])
        and 8 <= pr._mol_info(ctx.reg.smiles[i])[2] <= 16 and pr._mol_info(ctx.reg.smiles[i])[3]]
pick = np.random.default_rng(0).choice(cand, 60, replace=False)
B60 = {ctx.reg.ik[int(i)]: 1.0 for i in pick}
B5 = {"CCCCCCCCCC": 0.4, "CCCCCCCCCCCC": 0.25, "CCCCCCCCCCCCCCCC": 0.1, "Cc1ccccc1": 0.15, "Cc1cc(C)cc(C)c1": 0.1}
TS = (233.15, 288.15, 313.15)


def snap():
    out = {}
    for bn, b in (("b60", B60), ("b5", B5)):
        for T in TS:
            r = pr.predict_fuel(b, T_K=T, fuel_type="jet")
            for k, v in r.items():
                out[(bn, T, k)] = repr(v)
    return out


full = pr._v5_margins()
pr._flash_margins_cache = {k: v for k, v in full.items() if k != pr.SMOKE_MARGIN_KEY}
before = snap()
pr._flash_margins_cache = full
after = snap()
diff = sorted(k for k in before if before[k] != after[k])
print("records", len(before), "differing", len(diff))
for k in diff:
    print(k)
nonsmoke = [k for k in diff if k[2] != "smoke_point_mm"]
print("non-smoke differing:", len(nonsmoke), nonsmoke)
s = pr.predict_fuel(B5, T_K=288.15, fuel_type="jet")["smoke_point_mm"]
print("b5 smoke value", s["value"], "in_domain", s["in_domain"], "margin", s["margin"])
s = pr.predict_fuel(B60, T_K=288.15, fuel_type="jet")["smoke_point_mm"]
print("b60 smoke value", s["value"], "in_domain", s["in_domain"], "margin", s["margin"], "cov", s["coverage"]["substituted_mass"], s["coverage"]["unmapped_mass"], s["warnings"])
