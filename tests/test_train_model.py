"""Unit tests for Random Forest model training."""

import numpy as np
import pandas as pd
import pytest
from sklearn.ensemble import RandomForestClassifier

from crocodeel.conta_features import ContaminationFeatures
from crocodeel.train_model import (
    FeaturesComputerWorker,
    _build_performance_report,
    _filter_invalid_features,
    _reconstruct_sample_pairs,
    _train_model,
)
from crocodeel.exceptions import InputDataError


@pytest.fixture
def training_species_ab_table() -> pd.DataFrame:
    """Return a small abundance table with valid training sample names."""
    return pd.DataFrame(
        {
            "conta_source_case_1": [0.1, 0.2],
            "conta_target_case_1": [0.3, 0.4],
            "non_conta_source_case_2": [0.5, 0.6],
            "non_conta_target_case_2": [0.7, 0.8],
        },
        index=["species1", "species2"],
    )


@pytest.fixture
def sample_pairs() -> list[tuple[str, str]]:
    """Return sample pairs for testing."""
    return [
        ("conta_source_case_1", "conta_target_case_1"),
        ("non_conta_source_case_2", "non_conta_target_case_2"),
    ]


def test_reconstruct_sample_pairs(
    training_species_ab_table: pd.DataFrame,
) -> None:
    """Test reconstruction of source-target sample pairs."""
    sample_pairs = _reconstruct_sample_pairs(
        training_species_ab_table,
    )

    assert sample_pairs == [
        ("conta_source_case_1", "conta_target_case_1"),
        ("non_conta_source_case_2", "non_conta_target_case_2"),
    ]


def test_reconstruct_sample_pairs_rejects_invalid_sample_names(
    training_species_ab_table: pd.DataFrame,
) -> None:
    """Test that invalid sample names are rejected."""
    training_species_ab_table = training_species_ab_table.rename(
        columns={"conta_source_case_1": "invalid_sample"}
    )

    with pytest.raises(
        InputDataError,
        match="do not match the expected pattern",
    ):
        _reconstruct_sample_pairs(training_species_ab_table)


def test_reconstruct_sample_pairs_rejects_source_without_target() -> None:
    """Test that a source without a corresponding target is rejected."""
    species_ab_table = pd.DataFrame(
        {
            "conta_source_case_1": [0.1, 0.2],
        }
    )

    with pytest.raises(
        InputDataError,
        match="source samples without corresponding targets",
    ):
        _reconstruct_sample_pairs(species_ab_table)


def test_reconstruct_sample_pairs_rejects_target_without_source() -> None:
    """Test that a target without a corresponding source is rejected."""
    species_ab_table = pd.DataFrame(
        {
            "conta_target_case_1": [0.1, 0.2],
        }
    )

    with pytest.raises(
        InputDataError,
        match="target samples without corresponding sources",
    ):
        _reconstruct_sample_pairs(species_ab_table)


def test_filter_invalid_features() -> None:
    """Test removal of sample pairs with missing features."""
    all_features = np.array(
        [
            [1.0, 2.0, 3.0],
            [np.nan, np.nan, np.nan],
            [4.0, 5.0, 6.0],
            [7.0, np.nan, 9.0],
        ]
    )

    is_contaminated = np.array(
        [True, False, False, True],
    )

    filtered_features, filtered_labels = _filter_invalid_features(
        all_features,
        is_contaminated,
    )

    np.testing.assert_array_equal(
        filtered_features,
        np.array(
            [
                [1.0, 2.0, 3.0],
                [4.0, 5.0, 6.0],
            ]
        ),
    )

    np.testing.assert_array_equal(
        filtered_labels,
        np.array([True, False]),
    )


def test_filter_invalid_features_keeps_all_valid_features() -> None:
    """Test that valid features are preserved."""
    all_features = np.array(
        [
            [1.0, 2.0],
            [3.0, 4.0],
        ]
    )

    is_contaminated = np.array([True, False])

    filtered_features, filtered_labels = _filter_invalid_features(
        all_features,
        is_contaminated,
    )

    np.testing.assert_array_equal(
        filtered_features,
        all_features,
    )
    np.testing.assert_array_equal(
        filtered_labels,
        is_contaminated,
    )


def test_train_model() -> None:
    """Test Random Forest model training."""
    features = np.array(
        [
            [0.0, 0.0],
            [0.1, 0.2],
            [0.9, 0.8],
            [1.0, 1.0],
        ]
    )

    is_contaminated = np.array(
        [False, False, True, True],
    )

    model = _train_model(
        features,
        is_contaminated,
        ntrees=10,
        rng_seed=0,
        nproc=1,
    )

    assert isinstance(model, RandomForestClassifier)
    assert model.n_estimators == 10
    assert model.random_state == 0
    assert model.n_jobs == 1

    predictions = model.predict(features)

    assert predictions.shape == (4,)
    assert set(predictions).issubset({False, True})


def test_build_performance_report() -> None:
    """Test generation of the model performance report."""
    features_train = np.array(
        [
            [0.0, 0.0],
            [0.1, 0.1],
            [0.9, 0.9],
            [1.0, 1.0],
        ]
    )

    features_test = np.array(
        [
            [0.0, 0.1],
            [0.9, 1.0],
        ]
    )

    is_contaminated_train = np.array(
        [False, False, True, True],
    )

    is_contaminated_test = np.array(
        [False, True],
    )

    model = RandomForestClassifier(
        n_estimators=10,
        random_state=0,
    )
    model.fit(
        features_train,
        is_contaminated_train,
    )

    report = _build_performance_report(
        model,
        features_train,
        features_test,
        is_contaminated_train,
        is_contaminated_test,
    )

    assert set(report) == {"train", "test"}

    assert "classification_report" in report["train"]
    assert "classification_report" in report["test"]

    train_report = report["train"]["classification_report"]
    test_report = report["test"]["classification_report"]

    assert isinstance(train_report, dict)
    assert isinstance(test_report, dict)

    assert "accuracy" in train_report
    assert "accuracy" in test_report


def test_features_computer_worker_returns_none_when_features_cannot_be_extracted(
    training_species_ab_table: pd.DataFrame,
) -> None:
    """Test that failed feature extraction returns None."""
    worker = FeaturesComputerWorker(
        training_species_ab_table,
    )

    worker.feature_extractor.extract = lambda source, target: None

    result = worker.compute_features_sample_pair(
        (
            "conta_source_case_1",
            "conta_target_case_1",
        )
    )

    assert result is None


def test_features_computer_worker_returns_features(
    training_species_ab_table: pd.DataFrame,
) -> None:
    """Test that extracted features are returned as a NumPy array."""
    worker = FeaturesComputerWorker(
        training_species_ab_table,
    )

    expected_features = np.arange(
        ContaminationFeatures.NUM_FEATURES,
        dtype=float,
    )

    class FakeFeatures:
        values = expected_features

    worker.feature_extractor.extract = lambda source, target: FakeFeatures()

    result = worker.compute_features_sample_pair(
        (
            "conta_source_case_1",
            "conta_target_case_1",
        )
    )

    assert result is not None
    np.testing.assert_array_equal(
        result,
        expected_features,
    )
