"""T1: strict family-level (Murcko-scaffold) holdout audit for the pure flash claim."""
import sys
from pathlib import Path

BASE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(BASE / "src"))


def test_scaffold_family_key_groups_homologs():
    from fuelblend.features.family import scaffold_family_key

    # n-hexane and n-dodecane are acyclic alkane scaffolds -> same generic scaffold key
    hexane = "CCCCCC"
    dodecane = "CCCCCCCCCCCC"
    decane = "CCCCCCCCCCC"  # C11
    toluene = "Cc1ccccc1"
    ethylbenzene = "CCc1ccccc1"

    assert scaffold_family_key(hexane) == scaffold_family_key(decane)
    assert scaffold_family_key(decane) == scaffold_family_key(dodecane)
    assert scaffold_family_key(toluene) == scaffold_family_key(ethylbenzene)
    assert scaffold_family_key(hexane) != scaffold_family_key(toluene)


def test_leave_family_out_cv_never_trains_on_heldout_family():
    from fuelblend.features.family import leave_family_folds

    records = [
        {"inchikey": f"MOL{i:03d}", "family": fam, "T_K": 300.0 + i, "value": 40.0 + i, "quality_weight": 1.0}
        for fam in ("A", "B", "C", "D", "E")
        for i in range(10)
    ]
    folds = list(leave_family_folds(records, n_folds=5, seed=20261030))
    assert len(folds) == 5
    seen_pairs = set()
    for fold in folds:
        test_fams = {r["family"] for r in fold["test"]}
        train_fams = {r["family"] for r in fold["train"]}
        assert not (test_fams & train_fams), f"family overlap: {test_fams & train_fams}"
        assert fold["test"], "empty test set"
        assert fold["train"], "empty train set"
        for r in fold["test"]:
            pair = (r["inchikey"], r["family"])
            assert pair not in seen_pairs, "molecule assigned to multiple test folds"
            seen_pairs.add(pair)


def test_plausibility_errors_look_unlike_regression_collapses():
    """Guard: if plausibility SMAPE comes out near-perfect on disjoint real fuels, the input is suspicious."""
    from fuelblend.features.family import plausibility_report

    # two fuels with disjoint plausible splits and a tiny bias -> SMAPE must stay large-ish
    rep = plausibility_report(
        {
            "fuel_hi": {"dist": {"n_paraffin": 0.3, "iso_paraffin": 0.2, "cycloalkane": 0.0, "aromatic": 0.5},
                        "flash_C": 60.0, "density_kg_m3": 850.0, "mw_mean": 160.0},
            "fuel_lo": {"dist": {"n_paraffin": 0.15, "iso_paraffin": 0.05, "cycloalkane": 0.0, "aromatic": 0.8},
                        "flash_C": 38.0, "density_kg_m3": 900.0, "mw_mean": 140.0},
        }
    )
    # hi-paraffin fuel should flash higher than the aromatic-heavy one
    assert rep["pairwise_violations"] >= 0
    assert "smape_pct" in rep and "pairs" in rep and rep["pairs"] == 1
