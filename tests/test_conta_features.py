import numpy as np
import pandas as pd
import pytest

from crocodeel.conta_features import (ContaminationFeatureExtractor,
                                      _UnitSlopeRegression)

# ---------------------------------------------------------------------------
# _UnitSlopeRegression
# ---------------------------------------------------------------------------


def test_unit_slope_regression_fit():
    """Test fitting a unit-slope regression on log-transformed abundances."""
    model = _UnitSlopeRegression()

    X = np.array([[-6.0], [-5.0], [-4.0]])
    y = np.array([[-4.0], [-3.0], [-2.0]])

    result = model.fit(X, y)

    # The data follow y = x + 2.
    assert result is model
    assert model.coef_ == 1
    assert model.intercept_ == pytest.approx(2.0)


def test_unit_slope_regression_predict():
    """Test predictions on log-transformed abundances."""
    model = _UnitSlopeRegression()

    X = np.array([[-6.0], [-5.0], [-4.0]])
    y = np.array([[-4.0], [-3.0], [-2.0]])

    model.fit(X, y)

    result = model.predict(np.array([[-3.0], [-2.0]]))

    # The fitted contamination line is y = x + 2.
    expected = np.array([[-1.0], [0.0]])

    np.testing.assert_array_equal(result, expected)


def test_unit_slope_regression_score():
    """Test the negative mean squared error score."""
    model = _UnitSlopeRegression()

    X_train = np.array([[-6.0], [-5.0], [-4.0]])
    y_train = np.array([[-4.0], [-3.0], [-2.0]])

    model.fit(X_train, y_train)

    X = np.array([[-3.0], [-2.0]])
    y = np.array([[-1.0], [0.5]])

    # Predictions are [-1.0, 0.0], giving an MSE of 0.125.
    # score() returns the negative MSE.
    result = model.score(X, y)

    assert result == pytest.approx(-0.125)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def feature_extractor():
    """Create a ContaminationFeatureExtractor for testing."""
    species_ab_table = pd.DataFrame(
        {
            "sample1": [0.1, 0.01, 0.001, 0.0001, 0.00001],
            "sample2": [0.2, 0.02, 0.002, 0.0002, 0.00002],
        },
        index=[
            "species_1",
            "species_2",
            "species_3",
            "species_4",
            "species_5",
        ],
    )

    return ContaminationFeatureExtractor(species_ab_table)


# ---------------------------------------------------------------------------
# extract()
# ---------------------------------------------------------------------------


def test_extract_returns_none_with_too_few_candidates(
    feature_extractor,
    monkeypatch,
) -> None:
    """Test that extraction fails when too few candidate species are found."""
    candidate_species = np.zeros(
        (
            ContaminationFeatureExtractor.CONTA_LINE_MIN_NUM_SPECIES - 1,
            2,
        )
    )

    monkeypatch.setattr(
        feature_extractor,
        "_get_conta_line_candidate_species",
        lambda sample_pair_species_ab: candidate_species,
    )

    result = feature_extractor.extract("sample1", "sample2")

    assert result is None


def test_extract_returns_none_when_ransac_fails(
    feature_extractor,
    monkeypatch,
) -> None:
    """Test that extraction fails when RANSAC raises a ValueError."""
    candidate_species = np.zeros(
        (
            ContaminationFeatureExtractor.CONTA_LINE_MIN_NUM_SPECIES,
            2,
        )
    )

    def raise_value_error(conta_line_candidate_species_ab):
        raise ValueError("RANSAC failed")

    monkeypatch.setattr(
        feature_extractor,
        "_get_conta_line_candidate_species",
        lambda sample_pair_species_ab: candidate_species,
    )
    monkeypatch.setattr(
        feature_extractor,
        "_estimate_conta_line_offset",
        raise_value_error,
    )

    result = feature_extractor.extract("sample1", "sample2")

    assert result is None


def test_extract_returns_none_with_too_few_ransac_inliers(
    feature_extractor,
    monkeypatch,
) -> None:
    """Test that extraction fails when too few species belong to the line."""
    candidate_species = np.zeros(
        (
            ContaminationFeatureExtractor.CONTA_LINE_MIN_NUM_SPECIES,
            2,
        )
    )
    conta_line_species = np.zeros(
        (
            ContaminationFeatureExtractor.CONTA_LINE_MIN_NUM_SPECIES - 1,
            2,
        )
    )

    monkeypatch.setattr(
        feature_extractor,
        "_get_conta_line_candidate_species",
        lambda sample_pair_species_ab: candidate_species,
    )
    monkeypatch.setattr(
        feature_extractor,
        "_estimate_conta_line_offset",
        lambda conta_line_candidate_species_ab: (
            conta_line_species,
            2.0,
        ),
    )

    result = feature_extractor.extract("sample1", "sample2")

    assert result is None


def test_extract_returns_contamination_features(
    feature_extractor,
    monkeypatch,
) -> None:
    """Test successful extraction of contamination features."""
    candidate_species = np.zeros(
        (
            ContaminationFeatureExtractor.CONTA_LINE_MIN_NUM_SPECIES,
            2,
        )
    )
    conta_line_species = np.zeros(
        (
            ContaminationFeatureExtractor.CONTA_LINE_MIN_NUM_SPECIES,
            2,
        )
    )
    conta_line_offset = 2.0

    expected_features = np.arange(6, dtype=float)
    expected_species = [
        "species_1",
        "species_2",
    ]

    monkeypatch.setattr(
        feature_extractor,
        "_get_conta_line_candidate_species",
        lambda sample_pair_species_ab: candidate_species,
    )
    monkeypatch.setattr(
        feature_extractor,
        "_estimate_conta_line_offset",
        lambda conta_line_candidate_species_ab: (
            conta_line_species,
            conta_line_offset,
        ),
    )
    monkeypatch.setattr(
        feature_extractor,
        "_compute_features",
        lambda sample_pair_species_ab, conta_line_species_ab, offset:
            expected_features,
    )
    monkeypatch.setattr(
        feature_extractor,
        "_refine_conta_line",
        lambda sample_pair_species_ab, conta_line_species_ab:
            expected_species,
    )

    result = feature_extractor.extract("sample1", "sample2")

    assert result is not None
    np.testing.assert_array_equal(result.values, expected_features)
    assert result.conta_line_offset == conta_line_offset
    assert result.conta_line_species == expected_species


# ---------------------------------------------------------------------------
# _compute_features()
# ---------------------------------------------------------------------------


def test_compute_features(
    feature_extractor,
    monkeypatch,
) -> None:
    """Test computation of contamination features."""
    sample_pair_species_ab = np.array(
        [
            [-np.inf, 1.0],
            [-np.inf, 2.0],
            [3.0, 4.0],
            [5.0, 6.0],
        ]
    )

    conta_line_species_ab = np.array(
        [
            [1.0, 2.0],
            [2.0, 3.0],
            [3.0, 4.0],
        ]
    )

    conta_line_offset = 1.0

    monkeypatch.setattr(
        feature_extractor,
        "_get_mean_distance_to_nearest_neighbors",
        lambda data: 0.25,
    )
    monkeypatch.setattr(
        feature_extractor,
        "_get_mean_ab_top_source_specific_species",
        lambda data: 10.0,
    )
    monkeypatch.setattr(
        feature_extractor,
        "_get_diff_mean_ab_top10_source_species_vs_ab_cutoff1",
        lambda mean_ab, data, offset: 0.5,
    )
    monkeypatch.setattr(
        feature_extractor,
        "_get_diff_mean_ab_top10_source_species_vs_ab_cutoff2",
        lambda mean_ab, data: 0.75,
    )

    feature_extractor.ransac.n_trials_ = 100

    features = feature_extractor._compute_features(
        sample_pair_species_ab,
        conta_line_species_ab,
        conta_line_offset,
    )

    expected_distances = np.abs(
        conta_line_species_ab[:, 1]
        - conta_line_species_ab[:, 0]
        - conta_line_offset
    ) / np.sqrt(2)

    expected_mean_distance = expected_distances.mean()

    np.testing.assert_array_equal(
        features,
        np.array(
            [
                100,
                3,
                expected_mean_distance,
                0.25,
                0.5,
                0.75,
            ]
        ),
    )


def test_compute_features_without_source_specific_species(
    feature_extractor,
    monkeypatch,
) -> None:
    """Test feature computation when no source-specific species exist."""
    sample_pair_species_ab = np.array(
        [
            [1.0, 2.0],
            [2.0, 3.0],
            [3.0, 4.0],
        ]
    )

    conta_line_species_ab = np.array(
        [
            [1.0, 2.0],
            [2.0, 3.0],
        ]
    )

    monkeypatch.setattr(
        feature_extractor,
        "_get_mean_distance_to_nearest_neighbors",
        lambda data: 0.25,
    )

    feature_extractor.ransac.n_trials_ = 100

    features = feature_extractor._compute_features(
        sample_pair_species_ab,
        conta_line_species_ab,
        conta_line_offset=1.0,
    )

    expected_distances = np.abs(
        conta_line_species_ab[:, 1]
        - conta_line_species_ab[:, 0]
        - 1.0
    ) / np.sqrt(2)

    expected_mean_distance = expected_distances.mean()

    np.testing.assert_array_equal(
        features,
        np.array(
            [
                100,
                2,
                expected_mean_distance,
                0.25,
                0,
                0,
            ]
        ),
    )


# ---------------------------------------------------------------------------
# _get_conta_line_candidate_species()
# ---------------------------------------------------------------------------


def test_get_conta_line_candidate_species(feature_extractor):
    """Test selection of candidate species."""
    sample_pair_species_ab = np.array(
        [
            [-6.0, -4.0],
            [-5.0, -4.5],
            [-4.0, -5.0],
            [-3.0, -6.0],
        ]
    )

    result = feature_extractor._get_conta_line_candidate_species(
        sample_pair_species_ab
    )

    # Species [2] and [3] are excluded because source abundance
    # is lower than target abundance.
    expected = np.array(
        [
            [-6.0, -4.0],
            [-5.0, -4.5],
        ]
    )

    np.testing.assert_array_equal(result, expected)


def test_get_conta_line_candidate_species_excludes_target_specific_species(
    feature_extractor,
):
    """Test that species absent from the target are excluded."""
    sample_pair_species_ab = np.array(
        [
            [-np.inf, -4.0],
            [-6.0, -4.5],
            [-5.0, -3.0],
        ]
    )

    result = feature_extractor._get_conta_line_candidate_species(
        sample_pair_species_ab
    )

    # The first species is absent from the target (-inf), so it is excluded.
    expected = np.array(
        [
            [-6.0, -4.5],
            [-5.0, -3.0],
        ]
    )

    np.testing.assert_array_equal(result, expected)


def test_get_conta_line_candidate_species_keeps_equal_abundances(
    feature_extractor,
):
    """Test that species with equal source and target abundances are retained."""
    sample_pair_species_ab = np.array(
        [
            [-5.0, -5.0],
            [-4.0, -4.0],
        ]
    )

    result = feature_extractor._get_conta_line_candidate_species(
        sample_pair_species_ab
    )

    # Both species satisfy source >= target.
    expected = np.array(
        [
            [-5.0, -5.0],
            [-4.0, -4.0],
        ]
    )

    np.testing.assert_array_equal(result, expected)


def test_get_conta_line_candidate_species_keeps_species_at_quadrant_threshold(
    feature_extractor,
):
    """Test that species at the upper-left quadrant threshold are retained."""
    sample_pair_species_ab = np.array(
        [
            [-6.0, -3.0],
            [-5.0, -3.5],
            [-4.0, -3.9],
        ]
    )

    result = feature_extractor._get_conta_line_candidate_species(
        sample_pair_species_ab
    )

    # The last species has exactly two species in its upper-left quadrant,
    # which is the maximum allowed. All three species are therefore retained.
    expected = np.array(
        [
            [-6.0, -3.0],
            [-5.0, -3.5],
            [-4.0, -3.9],
        ]
    )

    np.testing.assert_array_equal(result, expected)


def test_get_conta_line_candidate_species_excludes_species_above_quadrant_threshold(
    feature_extractor,
):
    """Test that species above the upper-left quadrant threshold are excluded."""
    sample_pair_species_ab = np.array(
        [
            [-6.0, -2.0],
            [-5.0, -2.5],
            [-4.0, -3.0],
            [-3.0, -3.2],
        ]
    )

    result = feature_extractor._get_conta_line_candidate_species(
        sample_pair_species_ab
    )

    # The last species has three species in its upper-left quadrant,
    # exceeding the maximum of two, so it is excluded.
    expected = np.array(
        [
            [-6.0, -2.0],
            [-5.0, -2.5],
            [-4.0, -3.0],
        ]
    )

    np.testing.assert_array_equal(result, expected)


# ---------------------------------------------------------------------------
# _estimate_conta_line_offset()
# ---------------------------------------------------------------------------


def test_estimate_conta_line_offset(feature_extractor):
    """Test estimation of a contamination line offset."""
    conta_line_candidate_species_ab = np.array(
        [
            [-5.0, -3.0],
            [-4.0, -2.0],
            [-3.0, -1.0],
            [-2.0, 0.0],
            [-1.5, 0.5],
        ]
    )

    result_species_ab, result_offset = (
        feature_extractor._estimate_conta_line_offset(
            conta_line_candidate_species_ab
        )
    )

    # All species follow source = target + 2.
    # An offset of 2 corresponds to a contamination rate of 1%.
    np.testing.assert_allclose(result_offset, 2.0)

    # All species belong to the contamination line and should therefore
    # be identified as RANSAC inliers.
    np.testing.assert_array_equal(
        result_species_ab,
        conta_line_candidate_species_ab,
    )


def test_estimate_conta_line_offset_rejects_outlier(feature_extractor):
    """Test that RANSAC excludes an outlier from the contamination line."""
    conta_line_candidate_species_ab = np.array(
        [
            [-5.0, -3.0],
            [-4.0, -2.0],
            [-3.0, -1.0],
            [-2.0, 0.0],
            [-1.0, 3.0],  # Outlier
        ]
    )

    result_species_ab, result_offset = (
        feature_extractor._estimate_conta_line_offset(
            conta_line_candidate_species_ab
        )
    )

    # The first four species follow source = target + 2.
    # The last species is an outlier and should not affect the estimate.
    np.testing.assert_allclose(result_offset, 2.0)

    expected_species_ab = np.array(
        [
            [-5.0, -3.0],
            [-4.0, -2.0],
            [-3.0, -1.0],
            [-2.0, 0.0],
        ]
    )

    np.testing.assert_array_equal(
        result_species_ab,
        expected_species_ab,
    )


# ---------------------------------------------------------------------------
# _refine_conta_line()
# ---------------------------------------------------------------------------


def test_refine_conta_line_keeps_species_within_iqr(feature_extractor):
    """Test that species within the IQR bounds are retained."""
    sample_pair_species_ab = np.array(
        [
            [-6.0, -4.0],
            [-5.0, -3.1],
            [-4.0, -2.0],
            [-3.0, -1.1],
            [-2.0, 0.0],
        ]
    )

    conta_line_species_ab = sample_pair_species_ab.copy()

    result = feature_extractor._refine_conta_line(
        sample_pair_species_ab,
        conta_line_species_ab,
    )

    # All source-target offsets are between 1.9 and 2.0,
    # so all species are retained.
    expected = [
        "species_1",
        "species_2",
        "species_3",
        "species_4",
        "species_5",
    ]

    assert result == expected


def test_refine_conta_line_removes_low_offset_outlier(
    feature_extractor,
):
    """Test that a species below the lower IQR bound is removed."""
    sample_pair_species_ab = np.array(
        [
            [-6.0, -4.0],  # offset = 2.0
            [-5.0, -3.1],  # offset = 1.9
            [-4.0, -2.0],  # offset = 2.0
            [-3.0, -1.1],  # offset = 1.9
            [-2.0, -3.0],  # offset = -1.0
        ]
    )

    conta_line_species_ab = sample_pair_species_ab[:4]

    result = feature_extractor._refine_conta_line(
        sample_pair_species_ab,
        conta_line_species_ab,
    )

    # The last species has an offset of -1.0 and is below
    # the lower IQR bound.
    expected = [
        "species_1",
        "species_2",
        "species_3",
        "species_4",
    ]

    assert result == expected


def test_refine_conta_line_removes_high_offset_outlier(
    feature_extractor,
):
    """Test that a species above the upper IQR bound is removed."""
    sample_pair_species_ab = np.array(
        [
            [-6.0, -4.0],  # offset = 2.0
            [-5.0, -3.1],  # offset = 1.9
            [-4.0, -2.0],  # offset = 2.0
            [-3.0, -1.1],  # offset = 1.9
            [-2.0, 1.0],  # offset = 3.0
        ]
    )

    conta_line_species_ab = sample_pair_species_ab[:4]

    result = feature_extractor._refine_conta_line(
        sample_pair_species_ab,
        conta_line_species_ab,
    )

    # The last species has an offset of 3.0 and is above
    # the upper IQR bound.
    expected = [
        "species_1",
        "species_2",
        "species_3",
        "species_4",
    ]

    assert result == expected


# ---------------------------------------------------------------------------
# _get_mean_ab_top_source_specific_species()
# ---------------------------------------------------------------------------


def test_get_mean_ab_top_source_specific_species_with_fewer_than_10_species(
    feature_extractor,
):
    """Test mean abundance when fewer than 10 species are available."""
    source_specific_species_ab = np.array(
        [
            [-6.0, -3.0],
            [-5.0, -4.0],
            [-4.0, -2.0],
            [-3.0, -5.0],
        ]
    )

    result = feature_extractor._get_mean_ab_top_source_specific_species(
        source_specific_species_ab
    )

    # Fewer than 10 species are available, so all source abundances
    # are included in the mean.
    expected = np.mean([-3.0, -4.0, -2.0, -5.0])

    assert result == pytest.approx(expected)


def test_get_mean_ab_top_source_specific_species_with_exactly_10_species(
    feature_extractor,
):
    """Test mean abundance when exactly 10 species are available."""
    source_abundances = np.array(
        [
            -5.0,
            -4.5,
            -4.0,
            -3.5,
            -3.0,
            -2.5,
            -2.0,
            -1.5,
            -1.0,
            -0.5,
        ]
    )

    source_specific_species_ab = np.column_stack(
        [
            np.full(len(source_abundances), -np.inf),
            source_abundances,
        ]
    )

    result = feature_extractor._get_mean_ab_top_source_specific_species(
        source_specific_species_ab
    )

    # Exactly 10 species are available, so all 10 are included.
    expected = source_abundances.mean()

    assert result == pytest.approx(expected)


def test_get_mean_ab_top_source_specific_species_with_more_than_10_species(
    feature_extractor,
):
    """Test mean abundance of the 10 most abundant species."""
    source_abundances = np.array(
        [
            -6.0,
            -5.5,
            -5.0,
            -4.5,
            -4.0,
            -3.5,
            -3.0,
            -2.5,
            -2.0,
            -1.5,
            -1.0,
            -0.5,
        ]
    )

    source_specific_species_ab = np.column_stack(
        [
            np.full(len(source_abundances), -np.inf),
            source_abundances,
        ]
    )

    result = feature_extractor._get_mean_ab_top_source_specific_species(
        source_specific_species_ab
    )

    # In log10 space, the highest abundances are the least negative values:
    # -0.5, -1.0, ..., -5.0.
    expected = np.mean(
        [
            -0.5,
            -1.0,
            -1.5,
            -2.0,
            -2.5,
            -3.0,
            -3.5,
            -4.0,
            -4.5,
            -5.0,
        ]
    )

    assert result == pytest.approx(expected)


def test_get_mean_ab_top_source_specific_species_custom_num_species(
    feature_extractor,
):
    """Test selecting a custom number of most abundant species."""
    source_specific_species_ab = np.array(
        [
            [-np.inf, -5.0],
            [-np.inf, -3.0],
            [-np.inf, -1.0],
            [-np.inf, -2.0],
            [-np.inf, -4.0],
        ]
    )

    result = feature_extractor._get_mean_ab_top_source_specific_species(
        source_specific_species_ab,
        num_species=2,
    )

    # The two most abundant species are -1.0 and -2.0.
    expected = np.mean([-1.0, -2.0])

    assert result == pytest.approx(expected)
