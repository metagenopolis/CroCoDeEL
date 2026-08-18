"""Unit tests for contamination search."""

from unittest.mock import MagicMock

import numpy as np
import pandas as pd
import pytest

from crocodeel.conta_event import ContaminationEvent
from crocodeel.search_conta import (
    ContaminationSearcherWorker,
    _passes_cutoffs,
    _prepare_search,
)


@pytest.fixture
def species_ab_table() -> pd.DataFrame:
    """Return a small species abundance table."""
    return pd.DataFrame(
        {
            "sample1": [0.1, 0.2],
            "sample2": [0.3, 0.4],
        },
        index=["species1", "species2"],
    )


@pytest.fixture
def conta_event() -> ContaminationEvent:
    """Return a contamination event for testing."""
    return ContaminationEvent(
        source="source",
        target="target",
        rate=0.1,
        probability=0.9,
        conta_line_species=["species1"],
    )


def test_prepare_search_without_second_table(
    species_ab_table: pd.DataFrame,
) -> None:
    """Test preparation of a search using a single abundance table."""
    table, sample_pairs, num_sample_pairs = _prepare_search(
        species_ab_table,
        None,
    )

    pd.testing.assert_frame_equal(table, species_ab_table)

    assert num_sample_pairs == 4
    assert list(sample_pairs) == [
        ("sample1", "sample1"),
        ("sample1", "sample2"),
        ("sample2", "sample1"),
        ("sample2", "sample2"),
    ]


def test_prepare_search_with_second_table(
    species_ab_table: pd.DataFrame,
) -> None:
    """Test preparation of a search using two abundance tables."""
    species_ab_table_2 = pd.DataFrame(
        {
            "sample3": [0.5, 0.6],
            "sample4": [0.7, 0.8],
        },
        index=["species1", "species2"],
    )

    table, sample_pairs, num_sample_pairs = _prepare_search(
        species_ab_table,
        species_ab_table_2,
    )

    expected_table = pd.DataFrame(
        {
            "sample1": [0.1, 0.2],
            "sample2": [0.3, 0.4],
            "sample3": [0.5, 0.6],
            "sample4": [0.7, 0.8],
        },
        index=["species1", "species2"],
    )

    pd.testing.assert_frame_equal(table, expected_table)

    assert num_sample_pairs == 4
    assert list(sample_pairs) == [
        ("sample1", "sample3"),
        ("sample1", "sample4"),
        ("sample2", "sample3"),
        ("sample2", "sample4"),
    ]


def test_prepare_search_with_missing_species(
    species_ab_table: pd.DataFrame,
) -> None:
    """Test that missing species are filled with negative infinity."""
    species_ab_table_2 = pd.DataFrame(
        {
            "sample3": [0.5],
        },
        index=["species1"],
    )

    table, sample_pairs, num_sample_pairs = _prepare_search(
        species_ab_table,
        species_ab_table_2,
    )

    assert table.loc["species1", "sample3"] == 0.5
    assert table.loc["species2", "sample3"] == -np.inf

    assert num_sample_pairs == 2
    assert list(sample_pairs) == [
        ("sample1", "sample3"),
        ("sample2", "sample3"),
    ]


@pytest.mark.parametrize(
    "probability, rate, probability_cutoff, rate_cutoff, expected",
    [
        (0.9, 0.1, 0.5, 0.0, True),
        (0.5, 0.1, 0.5, 0.0, True),
        (0.9, 0.1, 0.9, 0.1, True),
        (0.49, 0.1, 0.5, 0.0, False),
        (0.9, 0.1, 0.95, 0.0, False),
        (0.9, 0.1, 0.5, 0.2, False),
        (0.5, 0.0, 0.5, 0.0, True),
    ],
)
def test_passes_cutoffs(
    probability: float,
    rate: float,
    probability_cutoff: float,
    rate_cutoff: float,
    expected: bool,
) -> None:
    """Test contamination event probability and rate cutoffs."""
    event = ContaminationEvent(
        source="source",
        target="target",
        rate=rate,
        probability=probability,
        conta_line_species=[],
    )

    assert (
        _passes_cutoffs(
            event,
            probability_cutoff,
            rate_cutoff,
        )
        is expected
    )


def test_worker_ignores_same_sample(
    species_ab_table: pd.DataFrame,
) -> None:
    """Test that a sample is never reported as contaminating itself."""
    classifier = MagicMock()

    worker = ContaminationSearcherWorker(
        species_ab_table,
        classifier,
    )

    result = worker.classify_sample_pair(
        ("sample1", "sample1"),
    )

    assert result is None
    classifier.predict_proba.assert_not_called()


def test_worker_returns_none_when_features_cannot_be_extracted(
    species_ab_table: pd.DataFrame,
) -> None:
    """Test that failed feature extraction produces no event."""
    classifier = MagicMock()

    worker = ContaminationSearcherWorker(
        species_ab_table,
        classifier,
    )

    worker.feature_extractor = MagicMock()
    worker.feature_extractor.extract.return_value = None

    result = worker.classify_sample_pair(
        ("sample1", "sample2"),
    )

    assert result is None
    classifier.predict_proba.assert_not_called()


def test_worker_classifies_sample_pair(
    species_ab_table: pd.DataFrame,
) -> None:
    """Test classification of a sample pair."""
    classifier = MagicMock()
    classifier.predict_proba.return_value = np.array([[0.1, 0.9]])

    worker = ContaminationSearcherWorker(
        species_ab_table,
        classifier,
    )

    features = MagicMock()
    features.values = np.array([1.0, 2.0, 3.0])
    features.conta_line_offset = 2.0
    features.conta_line_species = ["species1", "species2"]

    worker.feature_extractor = MagicMock()
    worker.feature_extractor.extract.return_value = features

    result = worker.classify_sample_pair(
        ("sample1", "sample2"),
    )

    assert result is not None
    assert result.source == "sample1"
    assert result.target == "sample2"
    assert result.probability == pytest.approx(0.9)
    assert result.rate == pytest.approx(0.01)
    assert result.conta_line_species == ["species1", "species2"]

    classifier.predict_proba.assert_called_once()

    prediction_input = classifier.predict_proba.call_args.args[0]

    assert prediction_input.shape == (1, 3)

    np.testing.assert_array_equal(
        prediction_input,
        np.array([[1.0, 2.0, 3.0]]),
    )
