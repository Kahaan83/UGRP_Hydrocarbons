"""Outputs of scripts/v5/base/delta_newdata_freeze_smoke.py exist, readings keys are present and row counts match the loaders."""
import json
from pathlib import Path

import pandas as pd

from fuelblend.data import loaders_freezing_lit as LF
from fuelblend.data import loaders_smoke_lit as LS

D = Path(__file__).resolve().parents[1] / "results" / "v5" / "BASE" / "delta_newdata_freeze_smoke"
SLE = Path(__file__).resolve().parents[1] / "data" / "processed" / "v5" / "thermoml_sle_hc.csv"


def test_outputs_exist():
    for f in ("README.md", "rows_freezing.csv", "rows_smoke.csv", "summary.csv", "readings.json", "rows_kittel_freezing_not_scorable.csv"):
        assert (D / f).exists(), f


def test_readings_keys():
    r = json.loads((D / "readings.json").read_text())
    for k in ("R-NF-1", "R-NF-2", "R-NF-3", "R-NF-4", "counts"):
        assert k in r
    assert {"thermoml_jet_cold", "thermoml_all", "ascent_mixtures", "vozka"} <= set(r["R-NF-2"])
    assert set(r["R-NF-4"]) >= {"bulk", "ysi", "direction"} and r["R-NF-4"]["direction"] == "over"


def test_freezing_row_counts():
    f = pd.read_csv(D / "rows_freezing.csv")
    sle = pd.read_csv(SLE)
    assert (f.dataset == "thermoml").sum() == (sle.point_kind == "liquidus").sum() == 838          # 11 eutectic rows excluded
    p = LF.load_properties()
    assert (f.dataset == "vozka").sum() == ((p.source_id == "vozka2018_si") & (p.property == "freezing_point_C")).sum() == 22
    asc = p[p.source_id.str.startswith("ascent065A") & (p.property == "freezing_point_C")]
    assert (f.dataset == "ascent").sum() == 17 and len(asc) == 19            # 19 raw values, 2 replicate pairs averaged
    k = pd.read_csv(D / "rows_kittel_freezing_not_scorable.csv")
    assert len(k) == 7 and int(k.censored.sum()) == 2


def test_smoke_row_counts():
    s = pd.read_csv(D / "rows_smoke.csv")
    d = LS.load_smoke_lit()
    d = d[(d.property == "smoke_point_mm") & ~d.qualifier.isin(["reference", "reported_mean"])]
    n_k = (s.source == "kittel2023").sum()
    assert n_k == 7
    assert (s.source != "kittel2023").sum() <= len(d)                        # replicates / methods averaged into fewer points
    assert (s.scope == "excluded_no_recipe").sum() == 3                      # Mensch JP-8 drum, Rainville Jet-A, Botero gasoline
    assert set(s.scope) <= {"smoke_mixtures", "smoke_pure", "smoke_rainville_calib", "smoke_kittel", "excluded_no_recipe"}
    summ = pd.read_csv(D / "summary.csv")
    assert {"mae", "mae_const", "margin_alpha10", "n_systems"} <= set(summ.columns)
