"""Unit tests for execution descriptions."""

import re
from pathlib import Path

from crocodeel.execution_description import ExecutionDescription


# ---------------------------------------------------------------------------
# ExecutionDescription.__str__()
# ---------------------------------------------------------------------------


def test_execution_description() -> None:
    """Test execution description with a single abundance table."""
    description = ExecutionDescription(
        species_ab_table_fp=Path("abundance.tsv"),
        species_ab_table_2_fp=None,
        rf_model_fp=Path("model.joblib"),
        filtering_ab_thr_factor=2.0,
        probability_cutoff=0.5,
        rate_cutoff=0.001,
    )

    result = str(description)

    assert result.startswith("# crocodeel version:")
    assert "rf_model: model.joblib" in result
    assert "species_ab_table: abundance.tsv" in result
    assert "filtering_ab_thr_factor: 2.0" in result
    assert "probability_cutoff: 0.5" in result
    assert "rate_cutoff: 0.001" in result

    # The second abundance table is only included when provided.
    assert "species_ab_table_2:" not in result


def test_execution_description_with_second_abundance_table() -> None:
    """Test execution description with two abundance tables."""
    description = ExecutionDescription(
        species_ab_table_fp=Path("abundance.tsv"),
        species_ab_table_2_fp=Path("abundance_2.tsv"),
        rf_model_fp=Path("model.joblib"),
        filtering_ab_thr_factor=None,
        probability_cutoff=0.5,
        rate_cutoff=0.001,
    )

    result = str(description)

    assert "rf_model: model.joblib" in result
    assert "species_ab_table: abundance.tsv" in result
    assert "species_ab_table_2: abundance_2.tsv" in result
    assert "filtering_ab_thr_factor: None" in result
    assert "probability_cutoff: 0.5" in result
    assert "rate_cutoff: 0.001" in result


# ---------------------------------------------------------------------------
# ExecutionDescription.datetime
# ---------------------------------------------------------------------------


def test_execution_description_datetime() -> None:
    """Test that the execution datetime has the expected format."""
    description = ExecutionDescription(
        species_ab_table_fp=Path("abundance.tsv"),
        species_ab_table_2_fp=None,
        rf_model_fp=Path("model.joblib"),
        filtering_ab_thr_factor=None,
        probability_cutoff=0.5,
        rate_cutoff=0.001,
    )

    assert re.fullmatch(
        r"\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}",
        description.datetime,
    )
