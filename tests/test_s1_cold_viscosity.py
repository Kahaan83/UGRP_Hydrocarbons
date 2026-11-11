"""v5 S1: cold-bin lab viscosity margins (scripts/v5/eval_lab_known_composition_v5.py)."""
import importlib.util
import json
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import yaml

ROOT = Path(__file__).resolve().parents[1]
_spec = importlib.util.spec_from_file_location("s1v5", ROOT / "scripts" / "v5" / "eval_lab_known_composition_v5.py")
s1 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(s1)
S1 = ROOT / "results" / "v5" / "S1"


def test_t_bin_edges_match_v4():
    assert s1.t_bin(263.14) == "T_lt_263K" and s1.t_bin(263.15) == "T_263_283K"
    assert s1.t_bin(283.14) == "T_263_283K" and s1.t_bin(283.15) == "T_ge_283K"


def test_hash_role_deterministic_and_about_20_percent():
    ids = [f"sys{i:05d}" for i in range(4000)]
    roles = [s1.hash_role(i) for i in ids]
    assert roles == [s1.hash_role(i) for i in ids]
    frac = roles.count("locked") / len(roles)
    assert 0.17 < frac < 0.23
    assert s1.hash_role("abc", frac=0.0) == "dev" and s1.hash_role("abc", frac=1.0) == "locked"


def test_assign_new_roles_keeps_v4_roles_and_hashes_the_rest():
    v4 = {"dev": ["a", "b"], "locked": ["c"]}
    sp = s1.assign_new_roles(["a", "c", "n1", "n2", "n3", "n4", "n5", "n6"], v4)
    assert sp["boehm_in_v4"] == {"dev": ["a"], "locked": ["c"]}
    new = sp["boehm_new"]
    for s in new["dev"]:
        assert s1.hash_role(s) == "dev"
    for s in new["locked"]:
        assert s1.hash_role(s) == "locked"
    assert set(new["dev"]) | set(new["locked"]) == {"n1", "n2", "n3", "n4", "n5", "n6"}
    assert "a" in sp["dev"] and "b" in sp["dev"] and "c" in sp["locked"]
    assert not set(sp["dev"]) & set(sp["locked"])


def test_write_split_refuses_to_overwrite_a_differing_split(tmp_path):
    p = tmp_path / "split.json"
    sp = {"dev": ["a"], "locked": ["b"], "x": 1}
    s1.write_split(p, sp)
    assert s1.write_split(p, sp) == sp
    with pytest.raises(SystemExit):
        s1.write_split(p, {"dev": ["a", "z"], "locked": ["b"], "x": 1})


def test_loo_system_offsets():
    resid = np.array([1.0, 1.0, 3.0, 5.0, 5.0, 5.0])
    sysid = np.array(["a", "a", "b", "c", "c", "c"])
    off, full = s1.loo_system_offsets(resid, sysid)
    assert full == pytest.approx((1 + 3 + 5) / 3)
    assert off[0] == off[1] == pytest.approx((3 + 5) / 2)
    assert off[2] == pytest.approx((1 + 5) / 2)
    assert off[3] == pytest.approx((1 + 3) / 2)
    off1, _ = s1.loo_system_offsets(np.array([1.0]), np.array(["a"]))
    assert np.isnan(off1).all()


def test_bias_correction_uses_dev_systems_only():
    def frame(locked_resid):
        rows = [dict(system_id=s, split_role="dev", method="refutas", t_bin="T_lt_263K", resid=r, source="x")
                for s, r in (("d1", -0.2), ("d2", -0.1), ("d3", -0.3))]
        rows.append(dict(system_id="l1", split_role="locked", method="refutas", t_bin="T_lt_263K", resid=locked_resid, source="x"))
        return pd.DataFrame(rows)

    a = s1.add_bias_corrected(frame(0.5))
    b = s1.add_bias_corrected(frame(-9.0))                   # wildly different locked residual
    ca, cb = (x[(x.method == "refutas_biascorr") & (x.split_role == "dev")].sort_values("system_id") for x in (a, b))
    assert np.allclose(ca.resid.to_numpy(), cb.resid.to_numpy())          # dev correction independent of locked rows
    lk = a[(a.method == "refutas_biascorr") & (a.split_role == "locked")]
    assert lk.offset.iloc[0] == pytest.approx(-0.2)                       # all-dev mean offset
    assert lk.resid.iloc[0] == pytest.approx(0.5 + 0.2)
    d1 = ca[ca.system_id == "d1"].iloc[0]
    assert d1.offset == pytest.approx((-0.1 - 0.3) / 2) and d1.resid == pytest.approx(-0.2 + 0.2)


def test_walther_extrapolation_recovers_an_exact_d341_line():
    A, B = 9.0, 3.4
    T = np.array(s1.WALTHER_TREF_K)
    nu = 10.0 ** (10.0 ** (A - B * np.log10(T))) - 0.7
    got = s1.walther_extrapolate(T, nu, 253.15)
    want = 10.0 ** (10.0 ** (A - B * np.log10(253.15))) - 0.7
    assert got == pytest.approx(want, rel=1e-6)


def test_arrhenius_log_mixing_rule():
    assert s1.arrhenius_ln_eta([0.5, 0.5], [np.log(1.0), np.log(4.0)]) == pytest.approx(np.log(2.0))
    assert s1.arrhenius_ln_eta([2, 2], [np.log(1.0), np.log(4.0)]) == pytest.approx(np.log(2.0))   # re-normalised


def test_margin_is_finite_iff_nine_dev_systems():
    from fuelblend.uq.conformal import two_sided_margin
    rng = np.random.default_rng(0)
    for n, finite in ((8, False), (9, True), (20, True)):
        g = np.repeat([f"s{i}" for i in range(n)], 3)
        r = rng.normal(0, 0.1, size=g.size)
        assert bool(np.isfinite(two_sided_margin(r, g, 0.1))) == finite


def test_method_key_names():
    assert s1.method_key("refutas", "T_lt_263K") == "viscosity_kin_cSt.lab_mixtures.raw_refutas.T_lt_263K"
    assert s1.method_key("walther_biascorr", "T_263_283K") == "viscosity_kin_cSt.lab_mixtures.walther_extrap_biascorr.T_263_283K"


def test_run_designs_margins_flag_is_optional():
    src = (ROOT / "scripts" / "v4" / "run_designs_v4.py").read_text(encoding="utf-8")
    assert "def main(out=None, margins=None, pure_props=None)" in src and "--margins" in src
    assert "mpath = None if margins is None" in src          # default behaviour (v4 yaml) unchanged


def test_boehm_table_matches_measurements_schema():
    head = pd.read_csv(ROOT / "data" / "unified" / "measurements_dedup.csv", nrows=1)
    b = pd.read_csv(ROOT / "data" / "real_fuels" / "v5" / "boehm2022_viscosity.csv")
    assert list(b.columns) == list(head.columns) and len(b) == 1470


def test_v4_inputs_and_outputs_untouched_by_s1():
    out = subprocess.run(["git", "diff", "--quiet", "HEAD", "--", "configs/uncertainty_v4.yaml", "results/v4", "splits/v4"], cwd=ROOT)
    assert out.returncode == 0


@pytest.mark.skipif(not (S1 / "uncertainty_v5.yaml").exists(), reason="run scripts/v5/eval_lab_known_composition_v5.py first")
def test_generated_v5_yaml_and_split():
    y = yaml.safe_load((S1 / "uncertainty_v5.yaml").read_text(encoding="utf-8"))
    m = y["margins"]
    for k, e in m.items():
        if k.startswith(s1.LAB_VISC_KEY):                    # inf iff fewer than 9 groups (generated, not hand-set)
            assert (e["value"] == "inf") == (e["n_groups"] < s1.MIN_DEV_SYSTEMS)
            assert ("applied_margin" in e) == (e["value"] != "inf")
    for b in s1.COLD_BINS:
        e = m[f"viscosity_kin_cSt.lab_mixtures.raw_refutas.{b}"]
        assert e["n_groups"] >= s1.MIN_DEV_SYSTEMS and np.isfinite(e["applied_margin"])
    v4 = yaml.safe_load((ROOT / "configs" / "uncertainty_v4.yaml").read_text(encoding="utf-8"))["margins"]
    for k, e in v4.items():
        if k not in {f"viscosity_kin_cSt.lab_mixtures.raw_refutas.{b}" for b in s1.COLD_BINS}:
            assert m[k] == e                                  # everything else copied unchanged
    sp = json.loads((S1 / "lab_systems_splits_v5.json").read_text())
    v4sp = json.loads((ROOT / "splits" / "v4" / "lab_systems_splits.json").read_text())
    assert set(v4sp["dev"]) <= set(sp["dev"]) and set(v4sp["locked"]) <= set(sp["locked"])
    for s in sp["boehm_new"]["dev"]:
        assert s1.hash_role(s) == "dev"
    for s in sp["boehm_new"]["locked"]:
        assert s1.hash_role(s) == "locked"
    from fuelblend.uq.conformal import two_sided_margin
    res = pd.read_csv(S1 / "viscosity_lab_residuals_v5.csv")
    g = res[(res.method == "refutas") & (res.split_role == "dev") & (res.t_bin == "T_lt_263K")]
    assert two_sided_margin(g.resid, g.system_id, 0.1) == pytest.approx(
        m["viscosity_kin_cSt.lab_mixtures.raw_refutas.T_lt_263K"]["applied_margin"])
