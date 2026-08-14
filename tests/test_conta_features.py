import numpy as np
import pandas as pd
import pytest

from crocodeel.conta_features import (
    ContaminationFeatureExtractor,
)


@pytest.fixture
def feature_extractor():
    """Create a ContaminationFeatureExtractor for testing."""
    species_ab_table = pd.DataFrame(
        {"sample1": [1.0]},
        index=["species_1"],
    )

    return ContaminationFeatureExtractor(species_ab_table)


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

    result = feature_extractor._get_conta_line_candidate_species(sample_pair_species_ab)

    # [2] and [3] are excluded because the source abundance is lower
    # than the target abundance.
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

    result = feature_extractor._get_conta_line_candidate_species(sample_pair_species_ab)

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

    result = feature_extractor._get_conta_line_candidate_species(sample_pair_species_ab)

    # Both species satisfy source >= target, so both are retained.
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

    result = feature_extractor._get_conta_line_candidate_species(sample_pair_species_ab)

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

    result = feature_extractor._get_conta_line_candidate_species(sample_pair_species_ab)

    # The last species has three species in its upper-left quadrant,
    # exceeding the maximum of two, so it is excluded.
    # The other species have at most two.
    expected = np.array(
        [
            [-6.0, -2.0],
            [-5.0, -2.5],
            [-4.0, -3.0],
        ]
    )

    np.testing.assert_array_equal(result, expected)


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

    result_species_ab, result_offset = feature_extractor._estimate_conta_line_offset(
        conta_line_candidate_species_ab
    )

    # All species follow source = target + 2.
    # An offset of 2 corresponds to a contamination rate of 1%.
    np.testing.assert_allclose(result_offset, 2.0)

    # All species belong to the contamination line and should
    # therefore be identified as RANSAC inliers.
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

    result_species_ab, result_offset = feature_extractor._estimate_conta_line_offset(
        conta_line_candidate_species_ab
    )

    # The first four species follow source = target + 2.
    # The last species is an outlier and should not affect the
    # estimated contamination line.
    np.testing.assert_allclose(result_offset, 2.0)

    expected_species_ab = np.array(
        [
            [-5.0, -3.0],
            [-4.0, -2.0],
            [-3.0, -1.0],
            [-2.0, 0.0],
        ]
    )

    np.testing.assert_array_equal(result_species_ab, expected_species_ab)
