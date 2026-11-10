"""v4 review fix 1: realistic representatives + r2 CRC reconstruction (stub reps; no registry needed)."""
import inspect

import numpy as np
import pandas as pd
import pytest

from fuelblend.data import reconstruction_r2 as r2
from fuelblend.data.reconstruction_r2 import is_realistic, reconstruct


def test_cage_and_polycyclic_rejected():
    assert not is_realistic("C1C[C@@H]2[C@H]3CC[C@H](C3)[C@@H]2C1", "cycloalkane")   # JP-10
    assert not is_realistic("C1C2CC3C4CC5CC(C14)C(C2)C3C5", "cycloalkane")          # diamantane
    assert not is_realistic("c1ccc2cc3ccccc3cc2c1", "aromatic_mono")                # anthracene
    assert not is_realistic("c1ccc(-c2ccccc2)cc1", "aromatic_mono")                 # biphenyl
    assert not is_realistic("CC(CC(C)(C)C)CC(C)(C)C", "iso_paraffin")               # quaternary carbons


def test_realistic_accepted():
    assert is_realistic("CCCCC1CCCCC1", "cycloalkane")
    assert is_realistic("C1CCC2CCCCC2C1", "cycloalkane")                            # decalin
    assert is_realistic("CCCCCc1ccccc1", "aromatic_mono")
    assert is_realistic("c1ccc2c(c1)CCCC2", "aromatic_mono")                        # tetralin
    assert is_realistic("Cc1cccc2ccccc12", "aromatic_di")
    assert is_realistic("CCCCCCCC(C)C", "iso_paraffin")
    assert is_realistic("CCCCCCCC", "n_paraffin")
    assert is_realistic("CCCCCC=C", "alkene_alkyne")
    assert not is_realistic("CCCCCC#C", "alkene_alkyne")


def test_deuterated_rejected_everywhere():
    d = "[2H]C([2H])([2H])C([2H])([2H])C([2H])([2H])[2H]"
    assert not is_realistic(d, "n_paraffin") and not is_realistic(d, "iso_paraffin")


class StubReps:
    """get() returns a fake registry index; None for the kinds in ``missing``."""

    def __init__(self, missing=()):
        self.missing = set(missing)
        self.ik = [f"IK{i}" for i in range(400)]
        self.smiles = [f"S{i}" for i in range(400)]
        self.kinds = ["n_paraffin", "iso_paraffin", "cycloalkane", "aromatic_mono", "aromatic_di", "alkene_alkyne"]

    def get(self, kind, n, max_gap=3):
        if kind in self.missing:
            return None
        return self.kinds.index(kind) * 50 + int(n)


ROW = dict(aromatics_vol_pct=18.0, olefins_vol_pct=1.0, diaromatics_vol_pct=2.0,
           d86_ibp=157.0, d86_t10=173.0, d86_t20=180.0, d86_t50=199.0, d86_t90=241.0, d86_fbp=257.0)
SPLIT = (0.3, 0.4, 0.3)


def test_missing_bin_goes_to_fallback():
    fr, fb = reconstruct(ROW, SPLIT, StubReps(missing=["cycloalkane"]))
    full, fb0 = reconstruct(ROW, SPLIT, StubReps())
    assert fb0 == pytest.approx(0.0, abs=1e-12)
    cyc = full[full.kind == "cycloalkane"].fraction.sum()
    assert cyc > 0.1
    assert fb == pytest.approx(cyc, abs=1e-9)
    assert not (fr.kind == "cycloalkane").any()


def test_short_d86_is_skipped():
    row = dict(ROW)
    for k in ("d86_t20", "d86_t50", "d86_t90"):
        row.pop(k)
    fr, fb = reconstruct(row, SPLIT, StubReps())
    assert fr.empty and fb == 1.0
    fr, fb = reconstruct({**ROW, "aromatics_vol_pct": np.nan}, SPLIT, StubReps())
    assert fr.empty and fb == 1.0


def test_fractions_sum():
    fr, fb = reconstruct(ROW, SPLIT, StubReps())
    assert list(fr.columns) == ["inchikey", "smiles", "fraction", "kind"]
    assert fr.fraction.sum() + fb == pytest.approx(1.0, abs=1e-9)
    fr, fb = reconstruct(ROW, SPLIT, StubReps(missing=["alkene_alkyne", "aromatic_di"]))
    assert fb > 0 and fr.fraction.sum() + fb == pytest.approx(1.0, abs=1e-9)


def test_reconstruction_reads_no_target():
    src = inspect.getsource(r2)
    body = "\n".join(l for l in src.splitlines() if "RHO_CLASS" not in l or "import" not in l).lower()
    for w in ("density", "flash", "freez", "visc", "heat_of_combustion"):
        assert w not in body, w


def test_nominal_split_sums_to_one(tmp_path):
    p = tmp_path / "jet.csv"
    pd.DataFrame([
        dict(fuel_id="A", corpus_status="included",
             class_fractions='{"n_paraffin": 0.2, "iso_paraffin": 0.3, "cycloalkane": 0.3, "aromatic": 0.2}'),
        dict(fuel_id="B", corpus_status="included",
             class_fractions='{"n_paraffin": 0.1, "iso_paraffin": 0.2, "cycloalkane": 0.5, "aromatic": 0.2}'),
        dict(fuel_id="C", corpus_status="included", class_fractions='{"paraffin_unsplit": 0.8, "aromatic": 0.2}'),
        dict(fuel_id="D", corpus_status="excluded_partial_composition",
             class_fractions='{"n_paraffin": 0.9, "aromatic": 0.0}'),
    ]).to_csv(p, index=False)
    s = r2.nominal_sat_split(p)
    assert sum(s) == pytest.approx(1.0, abs=1e-12)
    assert s[0] == pytest.approx((0.25 + 0.125) / 2)


def test_realistic_reps_requires_pure_props_and_flash_entry():
    smi = ["CCCCCCCC", "CCCCCCCCC", "CCCCCCCCCC", "C1CC2CC1C2", "[2H]C([2H])([2H])CCCCCC"]
    ik = [f"K{c}" for c in "EDCBA"]                      # deliberately not in sorted order

    class Reg:
        pass
    reg = Reg()
    reg.ik, reg.smiles, reg.usable = ik, smi, np.ones(5, bool)

    class Prov:
        valid = np.array([True, True, False, True, True])
        def ensure(self, idx):
            pass
    flash_ok = {smi[0], smi[1], smi[2], smi[3]}
    reps = r2.RealisticReps(reg, Prov(), vp_lookup=lambda s: 1 if s in flash_ok else None)
    assert reps.get("n_paraffin", 8) == 0
    assert reps.get("n_paraffin", 9) == 1
    assert reps.get("n_paraffin", 10) == 1             # decane has invalid pure props -> nearest within gap is C9
    assert reps.get("cycloalkane", 7) is None          # bicyclo[2.2.1]-type cage is rejected
    assert reps.get("n_paraffin", 8) != 4              # deuterated octane rejected
    assert reps.get("n_paraffin", 20) is None          # beyond max_gap


def test_skip_is_logged_and_excluded_end_to_end(caplog):
    import importlib.util, pathlib
    spec = importlib.util.spec_from_file_location("rebuild", pathlib.Path(__file__).parents[1] / "scripts/v4/rebuild_crc_reconstruction.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    p = pd.DataFrame([
        dict(fuel_id="CRC647_1", property="aromatics_vol_pct", value=18.0, method="ASTM D1319", role="descriptor"),
        dict(fuel_id="CRC647_1", property="distillation_T0_C", value=157.0, method="ASTM D86 (test 1)", role="descriptor"),
        dict(fuel_id="CRC647_1", property="distillation_T10_C", value=173.0, method="ASTM D86 (test 1)", role="descriptor"),
        dict(fuel_id="CRC647_1", property="distillation_T50_C", value=199.0, method="ASTM D86 (test 1)", role="descriptor"),
    ])
    rows = mod.descriptor_rows(p)
    with caplog.at_level("WARNING"):
        fr, fb = reconstruct(rows["CRC647_1"], SPLIT, StubReps())
    assert fr.empty and fb == 1.0
    assert any("fewer than 4 distillation points" in r.message for r in caplog.records)
    assert r2.diaromatic_mass(ROW) == pytest.approx(r2.diaromatic_mass(ROW))
    assert 0 < r2.diaromatic_mass(ROW) < 0.2
