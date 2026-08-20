"""Unit tests for Random Forest model training."""

from unittest.mock import MagicMock, patch

import numpy as np
import pandas as pd
import pytest
from sklearn.ensemble import RandomForestClassifier

from crocodeel.conta_features import ContaminationFeatures
from crocodeel.exceptions import InputDataError
from crocodeel.train_model import (FeaturesComputerDriver,
                                   FeaturesComputerWorker,
                                   _build_performance_report,
                                   _compute_features, _filter_invalid_features,
                                   _reconstruct_sample_pairs, _train_model,
                                   run_train_model)

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


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


# ---------------------------------------------------------------------------
# _reconstruct_sample_pairs()
# ---------------------------------------------------------------------------


def test_reconstruct_sample_pairs(
    training_species_ab_table: pd.DataFrame,
) -> None:
    """Test reconstruction of source-target sample pairs."""
    result = _reconstruct_sample_pairs(
        training_species_ab_table,
    )

    assert result == [
        ("conta_source_case_1", "conta_target_case_1"),
        ("non_conta_source_case_2", "non_conta_target_case_2"),
    ]


def test_reconstruct_sample_pairs_rejects_invalid_sample_names(
    training_species_ab_table: pd.DataFrame,
) -> None:
    """Test that invalid sample names are rejected."""
    training_species_ab_table = training_species_ab_table.rename(
        columns={
            "conta_source_case_1": "invalid_sample",
        }
    )

    with pytest.raises(
        InputDataError,
        match="do not match the expected pattern",
    ):
        _reconstruct_sample_pairs(
            training_species_ab_table,
        )


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
        _reconstruct_sample_pairs(
            species_ab_table,
        )


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
        _reconstruct_sample_pairs(
            species_ab_table,
        )


def test_reconstruct_sample_pairs_rejects_source_and_target_without_pair() -> None:
    """Test that unpaired source and target samples are rejected."""
    species_ab_table = pd.DataFrame(
        {
            "conta_source_case_1": [0.1, 0.2],
            "conta_target_case_2": [0.3, 0.4],
        }
    )

    with pytest.raises(
        InputDataError,
        match=(
            "source samples without corresponding targets.*"
            "target samples without corresponding sources"
        ),
    ):
        _reconstruct_sample_pairs(
            species_ab_table,
        )


# ---------------------------------------------------------------------------
# _compute_features()
# ---------------------------------------------------------------------------


def test_compute_features(
    training_species_ab_table: pd.DataFrame,
    sample_pairs: list[tuple[str, str]],
) -> None:
    """Test that _compute_features delegates to FeaturesComputerDriver."""
    expected_features = np.array(
        [
            [1.0, 2.0],
            [3.0, 4.0],
        ]
    )

    driver = MagicMock()
    driver.compute_all_features.return_value = expected_features

    with patch(
        "crocodeel.train_model.FeaturesComputerDriver",
        return_value=driver,
    ) as mock_driver:
        result = _compute_features(
            training_species_ab_table,
            sample_pairs,
            nproc=2,
        )

    mock_driver.assert_called_once_with(
        training_species_ab_table,
        sample_pairs,
        2,
    )

    driver.compute_all_features.assert_called_once_with()

    np.testing.assert_array_equal(
        result,
        expected_features,
    )


# ---------------------------------------------------------------------------
# _filter_invalid_features()
# ---------------------------------------------------------------------------


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

    is_contaminated = np.array(
        [True, False],
    )

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


# ---------------------------------------------------------------------------
# _train_model()
# ---------------------------------------------------------------------------


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

    assert isinstance(
        model,
        RandomForestClassifier,
    )
    assert model.n_estimators == 10
    assert model.random_state == 0
    assert model.n_jobs == 1

    predictions = model.predict(features)

    assert predictions.shape == (4,)


# ---------------------------------------------------------------------------
# _build_performance_report()
# ---------------------------------------------------------------------------


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


# ---------------------------------------------------------------------------
# FeaturesComputerWorker
# ---------------------------------------------------------------------------


def test_features_computer_worker_returns_none_when_features_cannot_be_extracted(
    training_species_ab_table: pd.DataFrame,
) -> None:
    """Test that failed feature extraction returns None."""
    worker = FeaturesComputerWorker(
        training_species_ab_table,
    )

    worker.feature_extractor.extract = (
        lambda source, target: None
    )

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

    worker.feature_extractor.extract = (
        lambda source, target: FakeFeatures()
    )

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


# ---------------------------------------------------------------------------
# FeaturesComputerDriver
# ---------------------------------------------------------------------------


def test_features_computer_driver() -> None:
    """Test that feature results are collected in sample-pair order."""
    species_ab_table = pd.DataFrame(
        {
            "sample1": [0.1, 0.2],
            "sample2": [0.3, 0.4],
        },
        index=["species1", "species2"],
    )

    sample_pairs = [
        ("sample1", "sample2"),
        ("sample2", "sample1"),
        ("sample1", "sample2"),
    ]

    features_1 = np.arange(
        ContaminationFeatures.NUM_FEATURES,
        dtype=float,
    )

    features_2 = (
        np.arange(
            ContaminationFeatures.NUM_FEATURES,
            dtype=float,
        )
        + 10
    )

    worker = MagicMock()
    pool = MagicMock()
    pool.__enter__.return_value = pool
    pool.imap.return_value = [
        features_1,
        None,
        features_2,
    ]

    with (
        patch(
            "crocodeel.train_model.FeaturesComputerWorker",
            return_value=worker,
        ) as mock_worker,
        patch(
            "crocodeel.train_model.Pool",
            return_value=pool,
        ) as mock_pool,
        patch(
            "crocodeel.train_model.tqdm",
            side_effect=lambda tasks, **kwargs: tasks,
        ),
    ):
        driver = FeaturesComputerDriver(
            species_ab_table=species_ab_table,
            sample_pairs=sample_pairs,
            nproc=2,
        )

        result = driver.compute_all_features()

    mock_worker.assert_called_once_with(
        species_ab_table,
    )

    mock_pool.assert_called_once_with(
        processes=2,
    )

    pool.imap.assert_called_once_with(
        worker.compute_features_sample_pair,
        sample_pairs,
        chunksize=FeaturesComputerDriver.DEFAULT_CHUNKSIZE,
    )

    np.testing.assert_array_equal(
        result[0],
        features_1,
    )

    # Failed feature extraction leaves the corresponding row as NaN.
    assert np.all(
        np.isnan(result[1])
    )

    np.testing.assert_array_equal(
        result[2],
        features_2,
    )


# ---------------------------------------------------------------------------
# run_train_model()
# ---------------------------------------------------------------------------


def test_run_train_model(
    training_species_ab_table: pd.DataFrame,
) -> None:
    """Test the orchestration of Random Forest model training."""
    sample_pairs = [
        ("conta_source_case_1", "conta_target_case_1"),
        ("non_conta_source_case_2", "non_conta_target_case_2"),
    ]

    all_features = np.array(
        [
            [1.0, 2.0],
            [3.0, 4.0],
        ]
    )

    features_train = np.array(
        [[1.0, 2.0]]
    )
    features_test = np.array(
        [[3.0, 4.0]]
    )

    is_contaminated_train = np.array(
        [True]
    )
    is_contaminated_test = np.array(
        [False]
    )

    rf_model = MagicMock()

    performance_report = {
        "train": {
            "classification_report": {},
        },
        "test": {
            "classification_report": {},
        },
    }

    model_fh = MagicMock()
    model_fh.name = "model.pkl"

    json_report_fh = MagicMock()

    with (
        patch(
            "crocodeel.train_model._reconstruct_sample_pairs",
            return_value=sample_pairs,
        ) as mock_reconstruct,
        patch(
            "crocodeel.train_model._compute_features",
            return_value=all_features,
        ) as mock_compute_features,
        patch(
            "crocodeel.train_model._filter_invalid_features",
            return_value=(
                all_features,
                np.array([True, False]),
            ),
        ) as mock_filter,
        patch(
            "crocodeel.train_model.train_test_split",
            return_value=(
                features_train,
                features_test,
                is_contaminated_train,
                is_contaminated_test,
            ),
        ) as mock_split,
        patch(
            "crocodeel.train_model._train_model",
            return_value=rf_model,
        ) as mock_train,
        patch(
            "crocodeel.train_model.joblib.dump",
        ) as mock_dump,
        patch(
            "crocodeel.train_model._build_performance_report",
            return_value=performance_report,
        ) as mock_build_report,
        patch(
            "crocodeel.train_model.json.dump",
        ) as mock_json_dump,
    ):
        run_train_model(
            species_ab_table=training_species_ab_table,
            model_fh=model_fh,
            json_report_fh=json_report_fh,
            test_size=0.3,
            ntrees=100,
            rng_seed=42,
            nproc=2,
        )

    mock_reconstruct.assert_called_once_with(
        training_species_ab_table,
    )

    mock_compute_features.assert_called_once_with(
        training_species_ab_table,
        sample_pairs,
        2,
    )

    # The source sample name determines the contamination label.
    mock_filter.assert_called_once()

    filter_args = mock_filter.call_args.args

    np.testing.assert_array_equal(
        filter_args[0],
        all_features,
    )

    np.testing.assert_array_equal(
        filter_args[1],
        np.array([True, False]),
    )

    mock_split.assert_called_once()

    split_args = mock_split.call_args.args
    split_kwargs = mock_split.call_args.kwargs

    np.testing.assert_array_equal(
        split_args[0],
        all_features,
    )

    np.testing.assert_array_equal(
        split_args[1],
        np.array([True, False]),
    )

    assert split_kwargs == {
        "test_size": 0.3,
        "random_state": 42,
    }

    mock_train.assert_called_once_with(
        features_train,
        is_contaminated_train,
        100,
        42,
        2,
    )

    rf_model.set_params.assert_called_once_with(
        n_jobs=1,
    )

    mock_dump.assert_called_once_with(
        rf_model,
        model_fh,
        compress=3,
    )

    mock_build_report.assert_called_once_with(
        rf_model,
        features_train,
        features_test,
        is_contaminated_train,
        is_contaminated_test,
    )

    mock_json_dump.assert_called_once_with(
        performance_report,
        json_report_fh,
        indent=4,
    )
