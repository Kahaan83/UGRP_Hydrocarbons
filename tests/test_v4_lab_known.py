"""v4 review fix 3: known-composition lab residuals, lab split and conformal margins."""
import hashlib
import importlib.util
import json
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import yaml
from rdkit import Chem

ROOT = Path(__file__).resolve().parents[1]
_spec = importlib.util.spec_from_file_location("eval_lab_known_composition",
                                               ROOT / "scripts" / "v4" / "eval_lab_known_composition.py")
lk = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(lk)
t_bin = lk.t_bin

BASE_COMMIT = "65212cc"           # last commit before this plan added keys (HEAD would be vacuous after the commit)
NEW_KEYS = ["density_kg_m3.lab_mixtures.ideal_volume.over", "density_kg_m3.lab_mixtures.ideal_volume.under",
            "viscosity_kin_cSt.lab_mixtures.raw_refutas.T_lt_263K", "viscosity_kin_cSt.lab_mixtures.raw_refutas.T_263_283K",
            "viscosity_kin_cSt.lab_mixtures.raw_refutas.T_ge_283K"]


def _yaml():
    return yaml.safe_load((ROOT / "configs" / "uncertainty_v4.yaml").read_text(encoding="utf-8"))


def test_t_bin_edges():
    assert t_bin(263.14) == "T_lt_263K"
    assert t_bin(263.15) == "T_263_283K"
    assert t_bin(283.14) == "T_263_283K"
    assert t_bin(283.15) == "T_ge_283K"


def test_cold_viscosity_bin_has_no_substitute():
    m = _yaml()["margins"]["viscosity_kin_cSt.lab_mixtures.raw_refutas.T_lt_263K"]
    if m["n_groups"] < 9:
        assert m["value"] == "inf" and "applied_margin" not in m


def test_new_keys_present_with_expected_shape():
    mg = _yaml()["margins"]
    for k in NEW_KEYS:
        assert k in mg, k
    assert mg[NEW_KEYS[0]]["direction"] == "over" and mg[NEW_KEYS[1]]["direction"] == "under"
    for k in NEW_KEYS[2:]:
        assert mg[k]["unit"] == "ln-ratio" and mg[k]["direction"] == "two_sided"


def test_existing_margin_keys_unchanged():
    old = subprocess.run(["git", "show", f"{BASE_COMMIT}:configs/uncertainty_v4.yaml"], cwd=ROOT, capture_output=True,
                         text=True, check=True).stdout
    old = yaml.safe_load(old)["margins"]
    new = _yaml()["margins"]
    assert len(old) > 20
    for k, v in old.items():
        assert k in new and new[k] == v, k


class _StubPD:
    """Tier-1 only for the ik's in t1; rho() mimics PureDensity.rho."""

    def __init__(self, iks):
        self.t1 = {ik: (700.0, -0.8, 280.0, 320.0, 5) for ik in iks}

    def rho(self, ik, cls, nc, T_K):
        v = self.t1.get(ik)
        if v is not None and v[2] - 30 <= T_K <= v[3] + 30:
            return v[0] + v[1] * (T_K - 288.71), 1
        return 750.0, 2


def _row(inchis, vals, value, T_K=298.15, ctype="mass fraction"):
    return dict(file="f.xml", prop="Mass density, kg/m3", prop_phase="Liquid", value=value, T_K=T_K, P_kPa=101.0,
                n_comp=len(inchis), inchis=json.dumps(inchis), names="[]", comp_type=ctype, comp_values=json.dumps(vals))


def test_system_without_pure_density_is_excluded():
    hexane = Chem.MolToInchi(Chem.MolFromSmiles("CCCCCC"))
    heptane = Chem.MolToInchi(Chem.MolFromSmiles("CCCCCCC"))
    octane = Chem.MolToInchi(Chem.MolFromSmiles("CCCCCCCC"))
    ik = lambda i: Chem.inchi.InchiToInchiKey(i)
    pdm = _StubPD([ik(hexane), ik(heptane)])               # octane has no tier-1 fit
    df = pd.DataFrame([_row([hexane, heptane], [0.5, 0.5], 650.0),
                       _row([hexane, octane], [0.5, 0.5], 660.0)])
    rows, stats = lk.build_density_rows(df, pdm)
    assert rows.system_id.nunique() == 1
    assert stats["n_excluded_no_pure_density"] == 1
    r = rows.iloc[0]
    assert r.predicted == pytest.approx(1.0 / (0.5 / (700 - 0.8 * (298.15 - 288.71)) * 2), rel=1e-9)
    assert r.resid == pytest.approx(r.predicted - r.measured)
    assert r.system_id == hashlib.md5("|".join(sorted([ik(hexane), ik(heptane)])).encode()).hexdigest()[:12]


def test_mole_fraction_converted_to_mass():
    hexane = Chem.MolToInchi(Chem.MolFromSmiles("CCCCCC"))
    decane = Chem.MolToInchi(Chem.MolFromSmiles("CCCCCCCCCC"))
    ik = Chem.inchi.InchiToInchiKey
    pdm = _StubPD([ik(hexane), ik(decane)])
    pdm.t1[ik(decane)] = (730.0, -0.7, 280.0, 320.0, 5)
    df = pd.DataFrame([_row([hexane, decane], [0.5, 0.5], 700.0, ctype="mole fraction")])
    rows, _ = lk.build_density_rows(df, pdm)
    mw6, mw10 = 86.178, 142.286
    w6 = mw6 / (mw6 + mw10)
    rho6 = 700 - 0.8 * (298.15 - 288.71); rho10 = 730 - 0.7 * (298.15 - 288.71)
    assert rows.predicted.iloc[0] == pytest.approx(1 / (w6 / rho6 + (1 - w6) / rho10), rel=1e-3)


def test_split_written_before_residuals_and_locked_disjoint():
    sp = json.loads((ROOT / "splits" / "v4" / "lab_systems_splits.json").read_text())
    assert not set(sp["locked"]) & set(sp["dev"])
    n = len(sp["locked"]) + len(sp["dev"])
    assert abs(len(sp["locked"]) - round(0.2 * n)) <= 1
    assert sp["seed"] == 20261111


def test_lab_locked_systems_not_in_margin_rows():
    cfg = _yaml()["margins"]
    sp = json.loads((ROOT / "splits" / "v4" / "lab_systems_splits.json").read_text())
    locked = set(sp["locked"])
    ops = {"eq": lambda s, v: s == v, "ne": lambda s, v: s != v}
    for k in NEW_KEYS:
        rc = cfg[k]["recompute"]
        df = pd.read_csv(ROOT / cfg[k]["residual_file"])
        sel = df
        for col, op, val in rc["where"]:
            sel = sel[ops[op](sel[col], val)]
        assert len(sel) == cfg[k]["n_rows"], k
        assert not (sel.split_role == "locked").any(), k
        assert not set(sel.system_id) & locked, k
        assert sel[rc["group"]].nunique() == cfg[k]["n_groups"], k
