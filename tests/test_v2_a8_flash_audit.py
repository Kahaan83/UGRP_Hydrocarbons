from pathlib import Path
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
CSV_PATH = ROOT / "results" / "v2" / "A8" / "flash_usage_audit.csv"
DOC_PATH = ROOT / "docs" / "progress" / "v2" / "A8.md"


def test_flash_usage_audit_output_exists():
    assert CSV_PATH.exists(), f"Missing {CSV_PATH}"
    assert DOC_PATH.exists(), f"Missing {DOC_PATH}"


def test_flash_usage_audit_content():
    df = pd.read_csv(CSV_PATH)
    required_cols = [
        "category", "source_id", "sample_kind", "test_method",
        "n_observations", "pct_of_total", "pipeline_usage"
    ]
    for col in required_cols:
        assert col in df.columns

    # Verify counts by sample kind
    kinds = df[df["category"] == "sample_kind"].set_index("sample_kind")
    assert kinds.loc["pure", "n_observations"] == 2726
    assert kinds.loc["real_fuel", "n_observations"] == 2487
    assert kinds.loc["binary", "n_observations"] == 123
    assert kinds.loc["multi", "n_observations"] == 3

    # Total defined lab mixtures is exactly 126
    total_mix = kinds.loc["binary", "n_observations"] + kinds.loc["multi", "n_observations"]
    assert total_mix == 126

    # Verify total sum
    total_obs = kinds["n_observations"].sum()
    assert total_obs == 5339
