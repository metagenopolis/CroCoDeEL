"""Unit tests for contamination search."""

import importlib
import logging
from unittest.mock import MagicMock, patch

import numpy as np
import pandas as pd
import pytest
from sklearn.ensemble import RandomForestClassifier

import crocodeel.search_conta as search_conta
from crocodeel.conta_event import ContaminationEvent
from crocodeel.conta_features import ContaminationFeatureExtractor
from crocodeel.search_conta import (
    ContaminationSearcherDriver,
    ContaminationSearcherWorker,
    Defaults,
    _load_rf_model,
    _log_search_results,
    _passes_cutoffs,
    _prepare_search,
    run_search_conta,
)

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


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


# ---------------------------------------------------------------------------
# _load_rf_model()
# ---------------------------------------------------------------------------


def test_load_rf_model() -> None:
    """Test loading the packaged Random Forest model."""
    model_fp = importlib.resources.files("crocodeel").joinpath(
        "models",
        "crocodeel_rf_Feb2026_2.joblib",
    )

    with model_fp.open("rb") as model_fh:
        model = _load_rf_model(model_fh)

    assert isinstance(model, RandomForestClassifier)


# ---------------------------------------------------------------------------
# _prepare_search()
# ---------------------------------------------------------------------------


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


# ---------------------------------------------------------------------------
# Defaults and cutoff handling
# ---------------------------------------------------------------------------


def test_search_conta_defaults() -> None:
    """Test default contamination search cutoffs."""
    assert Defaults.PROBABILITY_CUTOFF == 0.5
    assert Defaults.RATE_CUTOFF == 0.0


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


# ---------------------------------------------------------------------------
# ContaminationSearcherWorker
# ---------------------------------------------------------------------------


def test_worker_init(species_ab_table: pd.DataFrame) -> None:
    """Test that the worker initializes its feature extractor and classifier."""
    classifier = MagicMock()

    worker = ContaminationSearcherWorker(
        species_ab_table,
        classifier,
    )

    assert worker.rf_classifier is classifier
    assert isinstance(
        worker.feature_extractor,
        ContaminationFeatureExtractor,
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


# ---------------------------------------------------------------------------
# _log_search_results()
# ---------------------------------------------------------------------------


def test_log_search_results(caplog: pytest.LogCaptureFixture) -> None:
    """Test logging of contamination search results."""
    events = [
        ContaminationEvent(
            source="source1",
            target="target1",
            rate=0.1,
            probability=0.9,
            conta_line_species=[],
        ),
        ContaminationEvent(
            source="source2",
            target="target2",
            rate=0.2,
            probability=0.9,
            conta_line_species=[],
        ),
    ]

    with caplog.at_level(logging.INFO):
        _log_search_results(events)

    assert "2 contamination events detected" in caplog.text
    assert "2 samples contaminated" in caplog.text


def test_log_search_results_one_event(caplog: pytest.LogCaptureFixture) -> None:
    """Test logging of a single contamination event."""
    events = [
        ContaminationEvent(
            source="source1",
            target="target1",
            rate=0.1,
            probability=0.9,
            conta_line_species=[],
        ),
    ]

    with caplog.at_level(logging.INFO):
        _log_search_results(events)

    assert "1 contamination event detected" in caplog.text
    assert "1 sample contaminated" in caplog.text


def test_log_search_results_empty(caplog: pytest.LogCaptureFixture) -> None:
    """Test logging when no contamination events are detected."""
    with caplog.at_level(logging.INFO):
        _log_search_results([])

    assert "0 contamination events detected" in caplog.text
    assert "0 samples contaminated" in caplog.text


# ---------------------------------------------------------------------------
# run_search_conta()
# ---------------------------------------------------------------------------


def test_run_search_conta(
    species_ab_table: pd.DataFrame,
) -> None:
    """Test the contamination search orchestration."""
    rf_model = MagicMock()

    prepared_table = species_ab_table.copy()
    sample_pairs = iter(
        [
            ("sample1", "sample2"),
        ]
    )
    num_sample_pairs = 1

    low_rate_event = ContaminationEvent(
        source="source1",
        target="target1",
        rate=0.1,
        probability=0.8,
        conta_line_species=[],
    )

    high_rate_event = ContaminationEvent(
        source="source2",
        target="target2",
        rate=0.2,
        probability=0.8,
        conta_line_species=[],
    )

    # Return events unsorted to verify that run_search_conta()
    # sorts them by decreasing contamination rate.
    driver = MagicMock()
    driver.search_contamination.return_value = [
        low_rate_event,
        high_rate_event,
    ]

    rf_model_fh = MagicMock()

    with (
        patch(
            "crocodeel.search_conta._load_rf_model",
            return_value=rf_model,
        ) as mock_load_rf_model,
        patch(
            "crocodeel.search_conta._prepare_search",
            return_value=(
                prepared_table,
                sample_pairs,
                num_sample_pairs,
            ),
        ) as mock_prepare_search,
        patch(
            "crocodeel.search_conta.ContaminationSearcherDriver",
            return_value=driver,
        ) as mock_driver,
        patch(
            "crocodeel.search_conta._log_search_results",
        ) as mock_log_results,
        patch(
            "crocodeel.search_conta.perf_counter",
            side_effect=[0.0, 1.5],
        ),
    ):
        result = run_search_conta(
            species_ab_table=species_ab_table,
            species_ab_table_2=None,
            rf_model_fh=rf_model_fh,
            probability_cutoff=0.5,
            rate_cutoff=0.01,
            nproc=2,
        )

    mock_load_rf_model.assert_called_once_with(
        rf_model_fh,
    )

    mock_prepare_search.assert_called_once_with(
        species_ab_table,
        None,
    )

    mock_driver.assert_called_once_with(
        species_ab_table=prepared_table,
        all_sample_pairs=sample_pairs,
        num_sample_pairs=num_sample_pairs,
        rf_model=rf_model,
        probability_cutoff=0.5,
        rate_cutoff=0.01,
        nproc=2,
    )

    driver.search_contamination.assert_called_once_with()

    # Logging happens after sorting, so logged events follow the return order.
    mock_log_results.assert_called_once_with(
        [high_rate_event, low_rate_event],
    )

    # The returned events are sorted by decreasing contamination rate.
    assert result == [
        high_rate_event,
        low_rate_event,
    ]


def test_run_search_conta_no_events(
    species_ab_table: pd.DataFrame,
) -> None:
    """Test that a search with no detected events returns an empty list."""
    rf_model = MagicMock()
    rf_model_fh = MagicMock()

    prepared_table = species_ab_table.copy()
    sample_pairs = iter(
        [
            ("sample1", "sample2"),
        ]
    )
    num_sample_pairs = 1

    driver = MagicMock()
    driver.search_contamination.return_value = []

    with (
        patch(
            "crocodeel.search_conta._load_rf_model",
            return_value=rf_model,
        ),
        patch(
            "crocodeel.search_conta._prepare_search",
            return_value=(
                prepared_table,
                sample_pairs,
                num_sample_pairs,
            ),
        ),
        patch(
            "crocodeel.search_conta.ContaminationSearcherDriver",
            return_value=driver,
        ),
        patch(
            "crocodeel.search_conta._log_search_results",
        ) as mock_log_results,
        patch(
            "crocodeel.search_conta.perf_counter",
            side_effect=[0.0, 1.5],
        ),
    ):
        result = run_search_conta(
            species_ab_table=species_ab_table,
            species_ab_table_2=None,
            rf_model_fh=rf_model_fh,
            probability_cutoff=0.5,
            rate_cutoff=0.01,
            nproc=2,
        )

    assert result == []

    driver.search_contamination.assert_called_once_with()
    mock_log_results.assert_called_once_with([])


@pytest.mark.parametrize(
    "events",
    [
        (
            ContaminationEvent(
                source="source_b",
                target="target_a",
                rate=0.2,
                probability=0.8,
                conta_line_species=[],
            ),
            ContaminationEvent(
                source="source_a",
                target="target_b",
                rate=0.2,
                probability=0.8,
                conta_line_species=[],
            ),
            ContaminationEvent(
                source="source_a",
                target="target_a",
                rate=0.2,
                probability=0.8,
                conta_line_species=[],
            ),
        ),
        (
            ContaminationEvent(
                source="source_a",
                target="target_a",
                rate=0.2,
                probability=0.8,
                conta_line_species=[],
            ),
            ContaminationEvent(
                source="source_a",
                target="target_b",
                rate=0.2,
                probability=0.8,
                conta_line_species=[],
            ),
            ContaminationEvent(
                source="source_b",
                target="target_a",
                rate=0.2,
                probability=0.8,
                conta_line_species=[],
            ),
        ),
    ],
    ids=["reverse_order", "sorted_order"],
)
def test_run_search_conta_sorts_events_deterministically(
    species_ab_table: pd.DataFrame,
    events: tuple[ContaminationEvent, ...],
) -> None:
    """Test deterministic sorting of events with identical rates."""
    rf_model = MagicMock()
    rf_model_fh = MagicMock()

    prepared_table = species_ab_table.copy()
    sample_pairs = iter([("sample1", "sample2")])

    driver = MagicMock()
    driver.search_contamination.return_value = list(events)

    with (
        patch(
            "crocodeel.search_conta._load_rf_model",
            return_value=rf_model,
        ),
        patch(
            "crocodeel.search_conta._prepare_search",
            return_value=(
                prepared_table,
                sample_pairs,
                1,
            ),
        ),
        patch(
            "crocodeel.search_conta.ContaminationSearcherDriver",
            return_value=driver,
        ),
        patch(
            "crocodeel.search_conta._log_search_results",
        ) as mock_log_results,
        patch(
            "crocodeel.search_conta.perf_counter",
            side_effect=[0.0, 1.5],
        ),
    ):
        result = run_search_conta(
            species_ab_table=species_ab_table,
            species_ab_table_2=None,
            rf_model_fh=rf_model_fh,
            probability_cutoff=0.5,
            rate_cutoff=0.01,
            nproc=2,
        )

    expected = [
        ContaminationEvent(
            source="source_a",
            target="target_a",
            rate=0.2,
            probability=0.8,
            conta_line_species=[],
        ),
        ContaminationEvent(
            source="source_a",
            target="target_b",
            rate=0.2,
            probability=0.8,
            conta_line_species=[],
        ),
        ContaminationEvent(
            source="source_b",
            target="target_a",
            rate=0.2,
            probability=0.8,
            conta_line_species=[],
        ),
    ]

    # Events with identical rates are sorted by source, then target.
    assert result == expected

    # Logging happens after sorting and therefore uses the same order.
    mock_log_results.assert_called_once_with(expected)


# ---------------------------------------------------------------------------
# ContaminationSearcherDriver
# ---------------------------------------------------------------------------


def test_contamination_searcher_driver_init(
    species_ab_table: pd.DataFrame,
) -> None:
    """Test that ContaminationSearcherDriver stores its parameters."""
    all_sample_pairs = iter(
        [
            ("sample1", "sample2"),
        ]
    )
    rf_model = MagicMock()

    driver = ContaminationSearcherDriver(
        species_ab_table=species_ab_table,
        all_sample_pairs=all_sample_pairs,
        num_sample_pairs=1,
        rf_model=rf_model,
        probability_cutoff=0.5,
        rate_cutoff=0.01,
        nproc=2,
    )

    assert driver.species_ab_table is species_ab_table
    assert driver.all_sample_pairs is all_sample_pairs
    assert driver.num_sample_pairs == 1
    assert driver.rf_model is rf_model
    assert driver.probability_cutoff == 0.5
    assert driver.rate_cutoff == 0.01
    assert driver.nproc == 2


def test_contamination_searcher_driver_search_contamination(
    species_ab_table: pd.DataFrame,
) -> None:
    """Test that the driver filters and collects contamination events."""
    all_sample_pairs = iter(
        [
            ("sample1", "sample2"),
            ("sample2", "sample1"),
            ("sample1", "sample1"),
        ]
    )
    rf_model = MagicMock()

    valid_event = ContaminationEvent(
        source="source",
        target="target",
        rate=0.1,
        probability=0.9,
        conta_line_species=["species1"],
    )
    low_probability_event = ContaminationEvent(
        source="source2",
        target="target2",
        rate=0.2,
        probability=0.4,
        conta_line_species=[],
    )
    low_rate_event = ContaminationEvent(
        source="source3",
        target="target3",
        rate=0.005,
        probability=0.9,
        conta_line_species=[],
    )

    driver = ContaminationSearcherDriver(
        species_ab_table=species_ab_table,
        all_sample_pairs=all_sample_pairs,
        num_sample_pairs=3,
        rf_model=rf_model,
        probability_cutoff=0.5,
        rate_cutoff=0.01,
        nproc=2,
    )

    pool = MagicMock()

    pbar = MagicMock()
    pbar.__iter__.return_value = iter(
        [
            valid_event,
            None,
            low_probability_event,
            low_rate_event,
        ]
    )

    with (
        patch(
            "crocodeel.search_conta.Pool",
            return_value=pool,
        ) as mock_pool,
        patch(
            "crocodeel.search_conta.tqdm",
            return_value=pbar,
        ) as mock_tqdm,
    ):
        result = driver.search_contamination()

    assert result == [valid_event]

    mock_pool.assert_called_once_with(
        processes=2,
        initializer=search_conta._init_worker,
        initargs=(species_ab_table, rf_model),
    )

    mock_tqdm.assert_called_once()
    assert mock_tqdm.call_args.kwargs["total"] == 3
    assert mock_tqdm.call_args.kwargs["leave"] is False

    pbar.set_postfix_str.assert_called_once_with(
        "1 conta events found",
    )


# ---------------------------------------------------------------------------
# _init_worker() and _classify_sample_pair()
# ---------------------------------------------------------------------------


def test_init_worker_builds_module_level_worker(
    species_ab_table: pd.DataFrame,
) -> None:
    """Test that _init_worker builds the module-level worker."""
    classifier = MagicMock()

    search_conta._init_worker(species_ab_table, classifier)

    assert isinstance(
        search_conta._worker,
        ContaminationSearcherWorker,
    )
    assert search_conta._worker.rf_classifier is classifier


def test_classify_sample_pair_delegates_to_module_level_worker(
    species_ab_table: pd.DataFrame,
) -> None:
    """Test that _classify_sample_pair delegates to the module-level worker."""
    classifier = MagicMock()

    search_conta._init_worker(species_ab_table, classifier)

    # Same sample on both sides always returns None, regardless of the
    # classifier - this is enough to confirm delegation to the worker
    # built by _init_worker() without duplicating the more thorough
    # ContaminationSearcherWorker tests above.
    result = search_conta._classify_sample_pair(("sample1", "sample1"))

    assert result is None
