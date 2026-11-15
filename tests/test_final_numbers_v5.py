"""The final numbers pack is built from committed files only, and its headline values equal the sources.

Runs the builder into a temp dir (never results/v5/BASE) and recomputes each checked value from the source files
read with ``git show`` at the same pinned commit; no hand-typed numbers.
"""
import importlib.util
import io
import json
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import yaml

REPO = Path(__file__).resolve().parents[1]
SCRIPT = REPO / "scripts/v5/final/build_final_numbers_v5.py"


def _load():
    spec = importlib.util.spec_from_file_location("build_final_numbers_v5", SCRIPT)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


MOD = _load()
COMMIT = MOD.DEFAULT_COMMIT


def git(path):
    return subprocess.run(["git", "-C", str(REPO), "show", f"{COMMIT}:{path}"], capture_output=True, check=True).stdout


def gcsv(path):
    return pd.read_csv(io.BytesIO(git(path)), low_memory=False)


@pytest.fixture(scope="module")
def pack(tmp_path_factory):
    out = tmp_path_factory.mktemp("final")
    MOD.main(["--out", str(out), "--doc", str(out / "doc.md"), "--commit", COMMIT])
    return out


def test_outputs_exist(pack):
    for n in ["table_A_headline", "table_B_drift_by_cause", "table_C_margins", "table_D_ml_gates", "table_E_designs",
              "table_F_negative_results", "table_G_limitations"]:
        assert (pack / f"{n}.csv").exists()
    assert len(list((pack / "figures").glob("*.png"))) == 7
    assert (pack / "README.md").exists() and (pack / "doc.md").exists()


def test_refuses_base(tmp_path):
    with pytest.raises(SystemExit):
        MOD.main(["--out", str(REPO / "results/v5/BASE/final_x"), "--doc", str(tmp_path / "d.md")])


def test_headline_equals_scoreboards(pack):
    A = pd.read_csv(pack / "table_A_headline.csv")
    cur, frz = gcsv(MOD.SRC["sb_rerun"]), gcsv(MOD.SRC["sb_frozen"])
    for _, r in A.iterrows():
        def pick(sb):
            q = sb[(sb.property == r.property) & (sb.tier == r.tier) & (sb.set == r.set) & (sb.model == r.model)
                   & (sb.variant == r.variant) & (sb.label == r.label) & sb.T_K.isna() & (sb.status == r.status)
                   & (sb.cells_file == r.cells_file) & (sb.evidence_path.astype(str) == str(r.evidence_path))]
            assert len(q) == 1, (r.property, r.set, r.model)
            return q.iloc[0]
        c, f = pick(cur), pick(frz)
        for col, src in [("bias_current", c.bias), ("mae_current", c.mae), ("bias_frozen", f.bias),
                         ("mae_frozen", f.mae)]:
            assert np.isclose(r[col], src, equal_nan=True, rtol=1e-12), (col, r.property, r.set)
        assert str(r.verdict_current) == str(c.verdict) and str(r.verdict_frozen) == str(f.verdict)
    # the AJFTD flash example quoted in section B
    A2 = pd.read_csv(pack / "table_A2_t2_per_source.csv")
    aj = A2[(A2.property == "flash_point_C") & (A2.set == "ajftd_gcxgc") & (A2.model == "ideal_le_chatelier")].iloc[0]
    q = cur[(cur.property == "flash_point_C") & (cur.set == "ajftd_gcxgc") & (cur.model == "ideal_le_chatelier")
            & (cur.variant == "designated")].iloc[0]
    assert aj.bias_current == pytest.approx(q.bias)


def test_drift_and_margins(pack):
    B = pd.read_csv(pack / "table_B_drift_by_cause.csv")
    dr = gcsv(MOD.SRC["drift"])
    s = dr[(dr.path == MOD.SRC["sb_frozen"]) & (dr.kind == "cell")]
    assert B.set_index("cause").n_cells.to_dict() == s.cause.value_counts().to_dict()
    C = pd.read_csv(pack / "table_C_margins.csv").set_index("margin_key")
    y = yaml.safe_load(git(MOD.SRC["yaml"]))["margins"]
    for k in ["flash_point_C.lab_mixtures.ideal", "freezing_point_C.lab_liquidus.svl_unifac",
              "freezing_point_C.lab_liquidus.svl_verified"]:
        assert C.loc[k, "applied_margin"] == pytest.approx(y[k]["applied_margin"])
    v = json.loads(git(MOD.SRC["mlv"]))["R-MLV-6"]["ml"]
    for b in ["T_lt_263K", "T_263_283K", "T_ge_283K"]:
        assert C.loc[f"ml_pair_residual.{b}", "applied_margin"] == pytest.approx(v[b]["margin"])


def test_ml_and_designs(pack):
    D = pd.read_csv(pack / "table_D_ml_gates.csv")
    v = json.loads(git(MOD.SRC["mlv"]))
    d = json.loads(git(MOD.SRC["mld"]))
    for g, src in [("R-MLV-1", v["R-MLV-1"]["ml"]), ("R-MLV-2", v["R-MLV-2"]["ml"]), ("R-MLV-4", v["R-MLV-4"]),
                   ("R-MLD-1", d["R-MLD-1"]["ml"]), ("R-MLD-2", d["R-MLD-2"]["ml"])]:
        r = D[D.gate == g].iloc[0]
        assert r.mae_baseline == pytest.approx(src["mae_a"]) and r.mae_ml == pytest.approx(src["mae_b"])
    r = D[(D.gate == "R-MLD-3")].iloc[0]
    assert str(r["pass"]) == str(d["R-MLD-3"]["pass"])
    E = pd.read_csv(pack / "table_E_designs.csv")
    rt = gcsv(MOD.SRC["designs_table"])
    assert len(E) == len(rt)
    assert (E.certified == "Y").sum() == (rt.certified == "Y").sum()
    assert np.allclose(E.min_slack, rt.min_slack)
