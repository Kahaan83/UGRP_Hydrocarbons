"""v4 Task 15: octane refit on real bRON blends (R17 per-property targets, R19 blendstock groups, R6 in-fold selection)."""
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from fuelblend.models import octane_blend_v4 as O
from fuelblend.models import octane_ghosh_lumps as G

ROOT = Path(__file__).resolve().parents[1]
CORPUS = ROOT / "data" / "real_fuels" / "v4" / "gasoline_octane_blends.csv"
SPLITS = ROOT / "splits" / "v4" / "real_fuels_v4_splits.json"
real = pytest.mark.skipif(not (CORPUS.exists() and SPLITS.exists()), reason="v4 corpus absent")


def test_prf_linear_by_volume_exact():
    heptane, isooctane = "CCCCCCC", "CC(C)CC(C)(C)C"
    for pct in (0, 10, 25, 50, 87.5, 100):
        for prop in ("ron", "mon"):
            v = O.linear_by_volume_pure([(isooctane, pct / 100), (heptane, 1 - pct / 100)], prop)
            assert abs(v - pct) <= 1.0                         # PRF scale: ON == vol % iso-octane
            assert v == pytest.approx(pct)


def _toy_corpus():
    comps = lambda cas, f: json.dumps([{"name": "x", "inchikey": f"IK{cas}", "frac": f, "basis": "volume"},
                                       {"name": "L", "inchikey": "LUMP:BRON_BASE_RON90.3_MON84.7", "frac": 1 - f, "basis": "volume"}])
    base = dict(base_ron=90.3, base_mon=84.7, mon=80.0, target_ok=True, target_ok_mon=True)
    return pd.DataFrame([
        dict(blend_id="BRON:64-17-5:0.1", components=comps("a", 0.1), ron=92.0, target_ok_ron=True, **base),
        dict(blend_id="BRON:64-17-5:0.2", components=comps("a", 0.2), ron=93.0, target_ok_ron=True, **base),
        # OR column True but RON flag False and RON missing -> must never be a RON target
        dict(blend_id="BRON:534-22-5:0.3", components=comps("b", 0.3), ron=np.nan, target_ok_ron=False, **base),
        # FACE-like row with a RON value but flag False -> never a target
        dict(blend_id="PRF+10", components=comps("c", 0.1), ron=55.6, target_ok_ron=False,
             **{**base, "target_ok": False, "target_ok_mon": False}),
        # unresolved component: dropped, not imputed
        dict(blend_id="BRON:67-56-1:0.1", components=json.dumps([{"name": "x", "inchikey": "IKd", "frac": None, "basis": "unresolved"},
             {"name": "L", "inchikey": "LUMP:BRON_BASE_RON90.3_MON84.7", "frac": None, "basis": "unresolved"}]),
             ron=91.0, target_ok_ron=True, **base),
    ])


def test_rows_with_target_ok_ron_false_never_ron_targets():
    t, rep = O.prepare_blends(_toy_corpus(), "ron")
    assert list(t.blend_id) == ["BRON:64-17-5:0.1", "BRON:64-17-5:0.2"]
    assert rep["n_dropped_frac_none"] == 1 and rep["n_flag_false"] == 2
    tm, _ = O.prepare_blends(_toy_corpus(), "mon")           # MON flag is separate: the NaN-RON row IS a MON target
    assert "BRON:534-22-5:0.3" in set(tm.blend_id) and "PRF+10" not in set(tm.blend_id)


@real
def test_real_corpus_targets_respect_per_property_flags():
    df = pd.read_csv(CORPUS)
    for prop in ("ron", "mon"):
        t, rep = O.prepare_blends(df, prop)
        ok = set(df.loc[df[f"target_ok_{prop}"].astype(bool), "blend_id"])
        assert set(t.blend_id) <= ok and len(t) == rep["n_used"] == len(ok)
        assert not any(b.startswith(("PRF", "TPRF", "FACE", "TMB", "CP", "Hex", "MC")) for b in t.blend_id)
        assert t.y.notna().all() and t.base_on.notna().all()
    # the row whose RON is missing but target_ok (OR) is True must not be a RON target
    r = df[df.blend_id.str.startswith("BRON:534-22-5:0.3109")]
    if len(r):
        assert "BRON:534-22-5:0.3109" not in set(O.prepare_blends(df, "ron")[0].blend_id)


def _table(n_groups=8, per=3, seed=0):
    rng = np.random.default_rng(seed)
    rows = []
    for g in range(n_groups):
        bov = rng.normal(100, 8)
        for f in (0.1, 0.2, 0.3)[:per]:
            lump = ["alcohol", "ester", "ketone", "ethanol"][g % 4]
            rows.append(dict(group=f"G{g}", f=f, base_on=90.3, base_ron=90.3, base_mon=84.7, base_sens=5.6,
                             y=(1 - f) * 90.3 + f * bov, lump=lump, mw=60 + 5 * g, rho=800.0, n_c=3 + g,
                             o_mass_frac=0.2, unsat=0.0, ring=0.0, arom=0.0))
    return pd.DataFrame(rows)


@pytest.mark.parametrize("make", [lambda: O.GBMResidual(), lambda: O.GPResidual(), lambda: O.GhoshLumps("ron")])
def test_selection_inside_folds_ignores_test_labels_and_groups(make, monkeypatch):
    seen = []
    orig = G.grouped_inner_scores

    def spy(factory, params, train, k=5):
        seen.append(train.copy())
        return orig(factory, params, train, k)

    monkeypatch.setattr(G, "grouped_inner_scores", spy)
    t = _table()
    tr, te = t[t.group.isin([f"G{i}" for i in range(6)])], t[t.group.isin(["G6", "G7"])]
    p1 = O.fit_predict(make(), tr, te)
    te_bad = te.assign(y=te.y + 1e3)                       # garbage test labels must not change predictions
    p2 = O.fit_predict(make(), tr, te_bad)
    assert np.allclose(p1, p2)
    m = make().fit(tr)
    if not isinstance(m, O.GPResidual):          # GP selects by marginal likelihood on train inside fit()
        assert seen, "selection routine was never called"
    for s_ in seen:                              # every selection call saw training rows/groups only
        assert set(s_.group) <= set(tr.group) and not (set(s_.group) & set(te.group))
        assert len(s_) <= len(tr)
    assert m.selection_groups_ == set(tr.group) and not (m.selection_groups_ & set(te.group))


def test_reused_model_reselects_per_fit():
    t = _table()
    m = O.GBMResidual()
    m.fit(t[t.group.isin(["G0", "G1", "G2", "G3", "G4"])])
    first = dict(m.params_)
    m.fit(t[t.group.isin(["G3", "G4", "G5", "G6", "G7"])])
    assert m.params is None and m.selection_groups_ == {"G3", "G4", "G5", "G6", "G7"}
    assert first in O.GBMResidual.GRID and m.params_ in O.GBMResidual.GRID


def test_ghosh_fit_path_recovers_blending_values_and_keeps_old_api():
    t = _table(n_groups=8)
    mdl = G.fit(t, "ron")
    assert np.abs(mdl.predict(t) - t.y).mean() < 2.0
    # backward compatible helpers untouched
    assert G.compute_blending_octane_value(90.0, 92.0, 0.2) == pytest.approx(100.0)
    assert G.classify_component_to_lump("CCO", None) == "ethanol"
    assert G.blendstock_lump("CCCCO") == "alcohol" and G.blendstock_lump("CCCC(=O)OCC") == "ester"
    assert G.blendstock_lump("CC(C)C(C)(C)C") == "isooctane"


def test_cluster_bootstrap_paired_detects_clear_difference():
    g = np.repeat(np.arange(10), 3)
    a, b = np.full(30, 1.0), np.full(30, 3.0)
    d, lo, hi, ng = O.cluster_boot(g, a, b)
    assert d == pytest.approx(-2.0) and hi < 0 and ng == 10


# ---- DVPE diagnostic (Task 15 fix) -------------------------------------------------------
import sys  # noqa: E402
sys.path.insert(0, str(ROOT / "scripts" / "v4"))
import eval_dvpe_v4 as DV  # noqa: E402


def test_dvpe_cas_typos_and_row_filters():
    assert DV.fix_cas("137.32-6") == "137-32-6" and DV.fix_cas("141-78*-6") == "141-78-6" and DV.fix_cas("71.36-3") == "71-36-3"
    df = pd.DataFrame([
        dict(Blend_Name="10% ethanol in 2016-NREL BOB-01", Blend_DVPE="6.44", Blend_DVPE_temp=37.8, Blend_DVPE_source="s",
             **{"Blend_Primary Component": "10.1 wt% ethanol", "Blend component primary CAS": "64-17-5"}),
        dict(Blend_Name="10% methylfuran in 2016-NREL BOB-01", Blend_DVPE="5.71", Blend_DVPE_temp=25.0, Blend_DVPE_source="s",
             **{"Blend_Primary Component": "10.7 wt% methylfuran", "Blend component primary CAS": "534-22-5"}),
        dict(Blend_Name="Mixed ketones", Blend_DVPE="<6.89", Blend_DVPE_temp=37.8, Blend_DVPE_source="s",
             **{"Blend_Primary Component": "2-pentanone", "Blend component primary CAS": "107-87-9"}),
        dict(Blend_Name="10% ketone mix in 2016-NREL BOB-01", Blend_DVPE="4.97", Blend_DVPE_temp=37.8, Blend_DVPE_source="s",
             **{"Blend_Primary Component": "2-pentanone", "Blend component primary CAS": "107-87-9"}),
        dict(Blend_Name="Vertifuel", Blend_DVPE="35.4", Blend_DVPE_temp=38.0, Blend_DVPE_source="s",
             **{"Blend_Primary Component": np.nan, "Blend component primary CAS": np.nan}),
    ])
    t, ex = DV.parse_rows(df)
    assert list(t.cas) == ["64-17-5"] and t.f_vol.iloc[0] == 0.10 and t.y.iloc[0] == 6.44
    assert ex["excluded_censored"] == 1 and ex["excluded_temp_not_37p8"] == 1
    assert ex["excluded_named_mixture"] == 1 and ex["excluded_not_nrel_bob_oxygenate"] == 1


@pytest.mark.skipif(not DV.RAW.exists(), reason="raw Co-Optima file absent")
def test_dvpe_real_rows_map_to_bron_groups_and_have_no_neat_bob():
    t, ex = DV.parse_rows(pd.read_csv(DV.RAW, low_memory=False, encoding="latin-1"))
    sp = json.loads(SPLITS.read_text())
    assert set(t.group) <= set(sp["fuel_group"].values())
    assert (t.f_vol > 0).all()                      # no 0 % oxygenate (neat BOB) row exists
    assert ex["n_resolvable"] == len(t) > 40
