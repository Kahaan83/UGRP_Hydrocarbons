import sys, os, json
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
import pandas as pd

def _fake_df():
    rows = [
        # system made of molecules A only (warm)
        dict(system_id="A", property="viscosity_kin_cSt", kind="binary", T_key=300.0,
             T_K=300.0, P_kPa=101.325, value=1.0, x_mole="[0.5,0.5]",
             comp_inchikeys=json.dumps(["A","B"]), w_mass="[0.5,0.5]", conflict=False,
             n_obs=1, n_refs=1, n_outlier=0, vmin=0.9, vmax=1.1),
        dict(system_id="H2", property="viscosity_kin_cSt", kind="binary", T_key=300.0,
             T_K=300.0, P_kPa=101.325, value=2.0, x_mole="[0.5,0.5]",
             comp_inchikeys=json.dumps(["C","D"]), w_mass="[0.5,0.5]", conflict=False,
             n_obs=1, n_refs=1, n_outlier=0, vmin=1.9, vmax=2.1),
        dict(system_id="H3", property="viscosity_kin_cSt", kind="pure", T_key=300.0,
             T_K=300.0, P_kPa=101.325, value=3.0, x_mole="[1.0]",
             comp_inchikeys=json.dumps(["E"]), w_mass="[1.0]", conflict=False,
             n_obs=1, n_refs=1, n_outlier=0, vmin=2.9, vmax=3.1),
    ]
    return pd.DataFrame(rows)

def test_split_no_leak():
    from fuelblend.splits.split_generator import split_s1_s2_s3, leakage_check
    from fuelblend.eval import metrics
    import numpy as np
    df = _fake_df()
    res = split_s1_s2_s3(df, seed=20261030, hold_frac=0.4)
    leakage_check(df, res)
    mols_s2 = set(res["holdout_molecules_s2"]); mols_s3 = set(res["holdout_molecules_s3"])
    assert not (mols_s2 & mols_s3)
    assigned = set(res["S1"] + res["S2"] + res["S3"])
    for i in df.index:
        ms = set(json.loads(df.loc[i, "comp_inchikeys"]))
        if ms and (ms <= mols_s2 or ms <= mols_s3):
            assert i in assigned, f"held-out-molecule row {i} lost"
    y = np.array([1.0, 2.0, 3.0]); p = np.array([1.1, 1.8, 3.2]); c = np.array(["a","b","c"])
    rep = metrics.report(y, p, c, mode="viscosity")
    assert abs(rep["mae"] - np.mean([.1, .2, .2])) < 1e-9
    rep = metrics.report(y, p, c, mode="flash")
    assert rep["overpred_rate_5"] == 0.0  # no row overpredicted by >5
    print("ALL HARNESS TESTS PASSED")

def test_locked_no_overlap():
    import json as _json
    from fuelblend.splits.locked_sets import locked_sets
    from fuelblend.splits.split_generator import split_s1_s2_s3
    df = pd.read_csv("data/unified/measurements_dedup.csv", low_memory=False)
    locked_sets()
    lk = json.load(open("splits/locked.json"))
    res = split_s1_s2_s3(df)  # excludes locked rows by construction
    train = set(res["S1"]) | set(res["S2"])
    locked_fuels = set(lk["real_fuel_ids_locked"])
    locked_multi = set(lk["multi_system_ids_locked"])
    overlap = 0
    for i in train:
        r = df.loc[i]
        if r["kind"] == "real_fuel" and str(r["fuel_id"]) in locked_fuels:
            overlap += 1
        if r["system_id"] in locked_multi:
            overlap += 1
    assert overlap == 0, f"locked rows present in training split: {overlap}"
    print("LOCKED-SET LEAKAGE TEST PASSED")

if __name__ == "__main__":
    test_split_no_leak()
    test_locked_no_overlap()
