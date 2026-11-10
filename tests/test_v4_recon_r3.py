"""r3 CRC 647 reconstruction (Table G.1 hydrocarbon types): sub-kinds, class masses, spreading (stub reps)."""
import importlib.util
import inspect
import pathlib

import numpy as np
import pandas as pd
import pytest

from fuelblend.data import reconstruction_r3 as r3
from fuelblend.data.reconstruction_r3 import SUBKINDS, UNEXPLAINED, class_masses, reconstruct_r3, subkind
from fuelblend.mixture.fuel_hindcast import RHO_CLASS


@pytest.mark.parametrize("smiles, expected", [
    ("C1CCC2CCCCC2C1", "cyclo_di"),                       # decalin
    ("CCCCC1CCCCC1", "cyclo_mono"),                       # butylcyclohexane
    ("c1ccc2c(c1)CCCC2", "indan_tetralin"),               # tetralin
    ("c1ccc2c(c1)CCC2", "indan_tetralin"),                # indane
    ("CCCCCc1ccccc1", "alkylbenzene"),
    ("C=Cc1ccccc1", None),                                # styrene: non-aromatic double bond
    ("c1ccc(C2CCCCC2)cc1", None),                         # phenylcyclohexane: rings not fused
    ("Cc1cccc2ccccc12", "aromatic_di"),                   # 1-methylnaphthalene
    ("C1C2CC3C4CC5CC(C14)C(C2)C3C5", None),               # diamantane cage
    ("CCCCCCCC", "n_paraffin"),
    ("CCCCCCCC(C)C", "iso_paraffin"),                     # 2-methylnonane
    ("CCCCCC=C", "alkene_alkyne"),
    ("not a smiles", None),
])
def test_subkind(smiles, expected):
    assert subkind(smiles) == expected


def g1_row(**over):
    row = dict(aromatics_vol_pct=20.0, saturates_vol_pct=80.0,
               hc_type_paraffins_pct=50.0, hc_type_cycloparaffins_pct=30.0, hc_type_dicycloparaffins_pct=10.0,
               hc_type_tricycloparaffins_pct=10.0, hc_type_alkylbenzenes_pct=60.0, hc_type_indans_tetralins_pct=20.0,
               hc_type_indenes_pct=5.0, hc_type_naphthalene_pct=5.0, hc_type_naphthalenes_pct=5.0,
               hc_type_acenaphthenes_pct=0.0, hc_type_acenaphthylenes_pct=0.0, hc_type_tricyclic_aromatics_pct=5.0)
    row.update(over)
    return row


def test_class_masses_hand_computed():
    cm = class_masses(g1_row(), 0.0, 0.25)
    sat, aro = 80.0 * 0.79, 20.0 * 0.87                    # 63.2 and 17.4 (units: vol% x class factor)
    tot = sat + aro                                         # 80.6
    assert sum(cm.values()) == pytest.approx(1.0, abs=1e-12)
    assert cm["n_paraffin"] == pytest.approx(sat / tot * 0.50 * 0.25)
    assert cm["iso_paraffin"] == pytest.approx(sat / tot * 0.50 * 0.75)
    assert cm["cyclo_mono"] == pytest.approx(sat / tot * 0.30)
    assert cm["cyclo_di"] == pytest.approx(sat / tot * 0.10)
    assert cm["unexplained_tricycloparaffins"] == pytest.approx(sat / tot * 0.10)
    assert cm["alkylbenzene"] == pytest.approx(aro / tot * 0.60)
    assert cm["indan_tetralin"] == pytest.approx(aro / tot * 0.20)
    assert cm["aromatic_di"] == pytest.approx(aro / tot * 0.10)                 # naphthalene + naphthalenes
    assert cm["unexplained_minor_aromatics"] == pytest.approx(aro / tot * 0.10)  # indenes 5 + tricyclic 5
    assert cm["alkene_alkyne"] == 0.0
    assert set(cm) == set(SUBKINDS) | set(UNEXPLAINED)


def test_class_masses_olefins_and_type_shares_absorb_closure():
    cm = class_masses(g1_row(), 5.0, 0.3)
    s, a, o = 75.0 * RHO_CLASS["sat"], 20.0 * RHO_CLASS["aromatic"], 5.0 * RHO_CLASS["olefin"]
    assert sum(cm.values()) == pytest.approx(1.0, abs=1e-12)
    assert cm["alkene_alkyne"] == pytest.approx(o / (s + a + o))
    # type cells summing to 90 instead of 100 give identical masses: shares inside each fraction are used
    scaled = g1_row(**{k: v * 0.9 for k, v in g1_row().items() if k.startswith("hc_type_")})
    cm2 = class_masses(scaled, 5.0, 0.3)
    assert all(cm2[k] == pytest.approx(cm[k]) for k in cm)


def test_class_masses_degenerate_returns_none():
    assert class_masses(g1_row(aromatics_vol_pct=np.nan), 0.0, 0.3) is None
    sat_zero = g1_row(hc_type_paraffins_pct=0.0, hc_type_cycloparaffins_pct=0.0, hc_type_dicycloparaffins_pct=0.0,
                      hc_type_tricycloparaffins_pct=0.0)
    assert class_masses(sat_zero, 0.0, 0.3) is None
    aro_zero = {k: (0.0 if k in ("hc_type_alkylbenzenes_pct", "hc_type_indans_tetralins_pct", "hc_type_indenes_pct",
                                 "hc_type_naphthalene_pct", "hc_type_naphthalenes_pct", "hc_type_acenaphthenes_pct",
                                 "hc_type_acenaphthylenes_pct", "hc_type_tricyclic_aromatics_pct") else v)
                for k, v in g1_row().items()}
    assert class_masses(aro_zero, 0.0, 0.3) is None
    assert class_masses(aro_zero | {"aromatics_vol_pct": 0.0}, 0.0, 0.3) is not None   # no aromatics: nothing to split
    assert class_masses({k: v for k, v in g1_row().items() if k != "hc_type_paraffins_pct"}, 0.0, 0.3) is None


class StubReps:
    """get() returns a fake registry index (None for kinds in ``missing``) and records every call."""

    def __init__(self, missing=()):
        self.missing = set(missing)
        self.ik = [f"IK{i}" for i in range(500)]
        self.smiles = [f"S{i}" for i in range(500)]
        self.calls = []

    def get(self, kind, n, max_gap=3):
        self.calls.append((kind, n, max_gap))
        if kind in self.missing:
            return None
        return SUBKINDS.index(kind) * 50 + int(n)


ROW = dict(olefins_vol_pct=1.0, d86_ibp=157.0, d86_t10=173.0, d86_t20=180.0, d86_t50=199.0, d86_t90=241.0,
           d86_fbp=257.0)


def test_fractions_and_fallback_sum_to_one():
    reps = StubReps()
    fr, fb, det = reconstruct_r3(ROW, g1_row(), 0.29, reps)
    assert list(fr.columns) == ["inchikey", "smiles", "fraction", "kind"]
    assert fr.fraction.sum() + fb == pytest.approx(1.0, abs=1e-9)
    assert sum(det.values()) == pytest.approx(fb, abs=1e-12)
    assert set(det) == set(UNEXPLAINED)                        # every sub-kind has a representative here
    assert set(fr.kind) <= set(SUBKINDS)


def test_missing_kind_sends_exactly_its_mass_to_no_rep():
    g1 = g1_row(hc_type_tricycloparaffins_pct=0.0, hc_type_indenes_pct=0.0, hc_type_tricyclic_aromatics_pct=0.0)
    cm = class_masses(g1, ROW["olefins_vol_pct"], 0.29)
    fr, fb, det = reconstruct_r3(ROW, g1, 0.29, StubReps(missing=["cyclo_di"]))
    assert set(det) == {"no_rep_cyclo_di"}
    assert det["no_rep_cyclo_di"] == pytest.approx(cm["cyclo_di"], abs=1e-12)
    assert fb == pytest.approx(cm["cyclo_di"], abs=1e-12)
    assert not (fr.kind == "cyclo_di").any()
    assert fr.fraction.sum() + fb == pytest.approx(1.0, abs=1e-9)


def test_two_ring_kinds_use_nearest_carbon_number():
    reps = StubReps()
    reconstruct_r3(ROW, g1_row(), 0.29, reps)
    gaps = {}
    for kind, _, gap in reps.calls:
        gaps.setdefault(kind, set()).add(gap)
    for k in ("cyclo_di", "indan_tetralin", "aromatic_di"):
        assert gaps[k] == {99}, k
    for k in ("n_paraffin", "iso_paraffin", "cyclo_mono", "alkylbenzene", "alkene_alkyne"):
        assert gaps[k] == {3}, k


def test_short_d86_or_bad_g1_is_skipped():
    row = dict(ROW)
    for k in ("d86_t20", "d86_t50", "d86_t90"):
        row.pop(k)
    fr, fb, det = reconstruct_r3(row, g1_row(), 0.29, StubReps())
    assert fr.empty and fb == 1.0 and det == {"skipped": 1.0}
    fr, fb, det = reconstruct_r3(ROW, g1_row(aromatics_vol_pct=np.nan), 0.29, StubReps())
    assert fr.empty and fb == 1.0 and det == {"skipped": 1.0}


def test_reconstruction_reads_no_target():
    src = inspect.getsource(r3)
    body = "\n".join(l for l in src.splitlines() if "RHO_CLASS" not in l or "import" not in l).lower()
    for w in ("density", "flash", "freez", "visc", "heat_of_combustion"):
        assert w not in body, w


def test_sub_reps_binning_and_lowest_inchikey():
    smi = ["CCCCCCCC", "C1CCC2CCCCC2C1", "c1ccc2c(c1)CCCC2", "CCCCC1CCCCC1", "C=Cc1ccccc1", "CCCCCCCC(C)C"]
    ik = [f"K{c}" for c in "FEDCBA"]                    # deliberately unsorted

    class Reg:
        pass
    reg = Reg()
    reg.ik, reg.smiles, reg.usable = ik, smi, np.ones(len(smi), bool)

    class Prov:
        valid = np.ones(len(smi), bool)
        def ensure(self, idx):
            pass
    reps = r3.SubReps(reg, Prov(), vp_lookup=lambda s: 1)
    assert reps.get("cyclo_di", 14, max_gap=99) == 1
    assert reps.get("cyclo_di", 14) is None
    assert reps.get("indan_tetralin", 10) == 2
    assert reps.get("cyclo_mono", 10) == 3
    assert reps.get("alkylbenzene", 8) is None          # styrene is not an alkylbenzene
    assert reps.get("iso_paraffin", 9) == 5


def _load_script():
    spec = importlib.util.spec_from_file_location(
        "rebuild_r3", pathlib.Path(__file__).parents[1] / "scripts/v4/rebuild_crc_reconstruction_r3.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_skip_is_logged_and_excluded_end_to_end(caplog, tmp_path):
    mod = _load_script()
    recs = []
    for fid, props in (("CRC647_1", {"distillation_T0_C": 157.0, "distillation_T10_C": 173.0, "distillation_T20_C": 180.0,
                                     "distillation_T50_C": 199.0, "distillation_T90_C": 241.0, "distillation_T100_C": 257.0}),
                       ("CRC647_2", {"distillation_T0_C": 157.0, "distillation_T10_C": 173.0, "distillation_T50_C": 199.0})):
        recs.append(dict(fuel_id=fid, property="aromatics_vol_pct", value=18.0, method="ASTM D1319", role="descriptor"))
        recs.append(dict(fuel_id=fid, property="olefins_vol_pct", value=1.0, method="ASTM D1319", role="descriptor"))
        for k, v in props.items():
            recs.append(dict(fuel_id=fid, property=k, value=v, method="ASTM D86 (test 1)", role="descriptor"))
    rows = mod.descriptor_rows(pd.DataFrame(recs))
    g1 = {"CRC647_1": g1_row(), "CRC647_2": g1_row()}
    with caplog.at_level("WARNING"):
        out, stats, kind_mass, skipped, details = mod.build_compositions(rows, g1, 0.29, StubReps(missing=["aromatic_di"]))
    assert skipped == ["CRC647_2"]
    assert any("fewer than 4 distillation points" in r.message for r in caplog.records)
    df = pd.DataFrame(out, columns=mod.COLS)
    df.to_csv(tmp_path / "c.csv", index=False)
    assert set(df.system_id) == {"fuel:CRC647_1"}
    assert [s[0] for s in stats] == ["CRC647_1"] and list(details) == ["CRC647_1"]
    assert all(len(v) == 1 for v in kind_mass.values())
    lump = df[df.inchikey == "LUMP_NO_REALISTIC_REP"]
    assert len(lump) == 1 and bool(lump.is_lumped.iloc[0]) and lump.smiles.iloc[0] == ""
    assert lump.fraction.iloc[0] == pytest.approx(stats[0][1]) and lump.fraction.iloc[0] < 1.0
    assert (df.provenance == "reconstructed_r3").all() and (df.basis == "mass").all() and not df.verified.any()
    assert df.fraction.sum() == pytest.approx(1.0, abs=1e-9)
    assert sum(details["CRC647_1"].values()) == pytest.approx(stats[0][1], abs=1e-12)


def test_load_g1_maps_ids_and_zeroes_censored_cells(tmp_path):
    mod = _load_script()
    p = tmp_path / "g1.csv"
    pd.DataFrame([
        dict(fuel_code="008", property="aromatics_vol_pct", value_raw="18.5", value_std=18.5),
        dict(fuel_code="008", property="hc_type_indenes_pct", value_raw="<0.2", value_std=0.2),
        dict(fuel_code="101", property="aromatics_vol_pct", value_raw="20.0", value_std=20.0),
    ]).to_csv(p, index=False)
    got = mod.load_g1(p, ["CRC647_8", "CRC647_101"])
    assert got["CRC647_8"] == {"aromatics_vol_pct": 18.5, "hc_type_indenes_pct": 0.0}
    assert got["CRC647_101"]["aromatics_vol_pct"] == 20.0
    with pytest.raises(SystemExit):
        mod.load_g1(p, ["CRC647_8", "CRC647_999"])
