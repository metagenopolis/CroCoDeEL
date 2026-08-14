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
            [0.0, 3.0],
            [1.0, 2.0],
            [2.0, 1.0],
            [3.0, 0.0],
        ]
    )

    result = feature_extractor._get_conta_line_candidate_species(sample_pair_species_ab)

    # [2.0, 1.0] and [3.0, 0.0] are excluded because the source
    # abundance is lower than the target abundance.
    expected = np.array(
        [
            [0.0, 3.0],
            [1.0, 2.0],
        ]
    )

    np.testing.assert_array_equal(result, expected)


def test_get_conta_line_candidate_species_excludes_target_specific_species(
    feature_extractor,
):
    """Test that species absent from the target are excluded."""
    sample_pair_species_ab = np.array(
        [
            [-np.inf, 1.0],
            [0.0, 2.0],
            [1.0, 3.0],
        ]
    )

    result = feature_extractor._get_conta_line_candidate_species(sample_pair_species_ab)

    # The first species is absent from the target (-inf), so it is excluded.
    expected = np.array(
        [
            [0.0, 2.0],
            [1.0, 3.0],
        ]
    )

    np.testing.assert_array_equal(result, expected)


def test_get_conta_line_candidate_species_keeps_equal_abundances(
    feature_extractor,
):
    """Test that species with equal source and target abundances are retained."""
    sample_pair_species_ab = np.array(
        [
            [1.0, 1.0],
            [2.0, 2.0],
        ]
    )

    result = feature_extractor._get_conta_line_candidate_species(sample_pair_species_ab)

    # Both species satisfy source >= target, so both are retained.
    expected = np.array(
        [
            [1.0, 1.0],
            [2.0, 2.0],
        ]
    )

    np.testing.assert_array_equal(result, expected)


def test_get_conta_line_candidate_species_keeps_species_at_quadrant_threshold(
    feature_extractor,
):
    """Test that species at the upper-left quadrant threshold are retained."""
    sample_pair_species_ab = np.array(
        [
            [0.0, 3.0],
            [1.0, 2.5],
            [2.0, 2.1],
        ]
    )

    result = feature_extractor._get_conta_line_candidate_species(sample_pair_species_ab)

    # The first species has exactly two other species in its
    # upper-left quadrant, which is the maximum allowed.
    # All three species therefore remain candidates.
    expected = np.array(
        [
            [0.0, 3.0],
            [1.0, 2.5],
            [2.0, 2.1],
        ]
    )

    np.testing.assert_array_equal(result, expected)


def test_get_conta_line_candidate_species_excludes_species_above_quadrant_threshold(
    feature_extractor,
):
    """Test that species above the upper-left quadrant threshold are excluded."""
    sample_pair_species_ab = np.array(
        [
            [0.0, 5.0],
            [1.0, 4.0],
            [2.0, 3.0],
            [3.0, 2.5],
        ]
    )

    result = feature_extractor._get_conta_line_candidate_species(
        sample_pair_species_ab
    )

    # The last species has three species in its upper-left quadrant,
    # exceeding the maximum of two, so it is excluded.
    # The other species have at most two.
    expected = np.array(
        [
            [0.0, 5.0],
            [1.0, 4.0],
            [2.0, 3.0],
        ]
    )

    np.testing.assert_array_equal(result, expected)
