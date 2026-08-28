"""Unit tests for the specialised unit-slope RANSAC."""

from collections.abc import Iterator

import numpy as np
import pytest
from sklearn.base import BaseEstimator, RegressorMixin
from sklearn.linear_model import RANSACRegressor

from crocodeel.conta_features import ContaminationFeatureExtractor
from crocodeel.ransac import _dynamic_max_trials, fit_unit_slope_ransac

MAX_TRIALS = ContaminationFeatureExtractor.RANSAC_MAX_TRIALS
RESIDUAL_THRESHOLD = ContaminationFeatureExtractor.RANSAC_RESIDUAL_THRESHOLD
RANDOM_STATE = ContaminationFeatureExtractor.RANSAC_RANDOM_STATE


def fit(
    abscissa: np.ndarray,
    ordinate: np.ndarray,
) -> tuple[np.ndarray, float, int]:
    """Fit using the RANSAC parameters used by CroCoDeEL."""
    return fit_unit_slope_ransac(
        abscissa,
        ordinate,
        max_trials=MAX_TRIALS,
        residual_threshold=RESIDUAL_THRESHOLD,
        random_state=RANDOM_STATE,
    )


# ---------------------------------------------------------------------------
# Behaviour
# ---------------------------------------------------------------------------


def test_fits_the_offset_of_a_perfect_line() -> None:
    """Test that a noiseless unit-slope line is recovered exactly."""
    abscissa = np.array([-6.0, -5.0, -4.0, -3.0, -2.0, -1.0, 0.0, 1.0])
    ordinate = abscissa + 2.0

    inlier_mask, intercept, num_trials = fit(abscissa, ordinate)

    assert inlier_mask.all()
    assert intercept == pytest.approx(2.0)

    # A consensus set containing every point is found on the first trial,
    # so the dynamic stopping criterion ends the search immediately.
    assert num_trials == 1


def test_excludes_outliers_from_the_consensus_set() -> None:
    """Test that points outside the consensus set are excluded."""
    abscissa = np.arange(-8.0, 0.0)
    ordinate = abscissa + 1.0
    ordinate[2] += 5.0
    ordinate[5] -= 4.0

    inlier_mask, intercept, _ = fit(abscissa, ordinate)

    assert not inlier_mask[2]
    assert not inlier_mask[5]
    assert inlier_mask.sum() == len(abscissa) - 2
    assert intercept == pytest.approx(1.0)


def test_keeps_points_at_the_residual_threshold() -> None:
    """Test that the residual threshold is inclusive."""
    abscissa = np.zeros(10)
    ordinate = np.array(
        [0.0] * 8
        + [RESIDUAL_THRESHOLD, RESIDUAL_THRESHOLD + 1.0]
    )

    inlier_mask, _, _ = fit(abscissa, ordinate)

    # The boundary is inclusive, matching sklearn's RANSAC implementation.
    assert inlier_mask[8]
    assert not inlier_mask[9]


def test_raises_when_no_consensus_set_is_found() -> None:
    """Test that scattered points produce no valid consensus set."""
    # The residuals are deliberately irregular so that no sampled pair
    # produces even a single-point consensus set.
    abscissa = np.zeros(12)
    ordinate = np.array(
        [
            6.328,
            10.517,
            16.488,
            20.241,
            26.125,
            32.343,
            39.987,
            45.349,
            48.467,
            54.634,
            58.865,
            66.313,
        ]
    )

    with pytest.raises(ValueError, match="consensus set"):
        fit(abscissa, ordinate)

    with pytest.raises(ValueError):
        reference_fit(abscissa, ordinate)


def test_accepts_a_single_point_consensus_set() -> None:
    """Test that a single inlier is accepted, matching sklearn."""
    # Evenly spaced residuals make some sampled-pair midpoints coincide with
    # another observation. sklearn starts with one as the best inlier count,
    # so such a candidate can become the selected consensus set.
    abscissa = np.zeros(12)
    ordinate = np.arange(12.0) * 5.0

    inlier_mask, intercept, num_trials = fit(abscissa, ordinate)
    expected = reference_fit(abscissa, ordinate)

    assert inlier_mask.sum() == 1
    np.testing.assert_array_equal(inlier_mask, expected[0])
    assert intercept == expected[1]
    assert num_trials == expected[2]


def test_intercept_is_refitted_on_the_whole_consensus_set() -> None:
    """Test that the final intercept is fitted using all consensus points."""
    abscissa = np.zeros(9)
    offsets = np.array(
        [0.0, 0.05, -0.05, 0.1, -0.1, 0.02, -0.02, 0.08, -0.08]
    )
    ordinate = abscissa + 1.0 + offsets

    inlier_mask, intercept, _ = fit(abscissa, ordinate)

    assert inlier_mask.all()
    assert intercept == pytest.approx((ordinate - abscissa).mean())


def test_stop_probability_controls_when_the_search_ends() -> None:
    """Test that the confidence requirement decides how long the search runs.

    Half of these points lie on the identity line and half on y = x + 1.5.
    Demanding no confidence stops as soon as any candidate is accepted, on
    whichever line the first sampled pair happened to fall; the default keeps
    looking and finds the contamination line.
    """
    abscissa = np.array([-9.0, -8.0, -7.0, -2.0, -2.0, -1.0, 0.0, 0.0, 1.0, 1.0])
    ordinate = np.array([-9.0, -8.0, -5.5, -0.5, -2.0, 0.5, 1.5, 1.5, 1.0, 1.0])

    _, thorough_intercept, thorough_trials = fit_unit_slope_ransac(
        abscissa,
        ordinate,
        max_trials=MAX_TRIALS,
        residual_threshold=RESIDUAL_THRESHOLD,
        random_state=RANDOM_STATE,
    )
    _, immediate_intercept, immediate_trials = fit_unit_slope_ransac(
        abscissa,
        ordinate,
        max_trials=MAX_TRIALS,
        residual_threshold=RESIDUAL_THRESHOLD,
        random_state=RANDOM_STATE,
        stop_probability=0.0,
    )

    assert immediate_trials == 1
    assert immediate_intercept == pytest.approx(0.0)

    assert thorough_trials > immediate_trials
    assert thorough_intercept == pytest.approx(1.5)


def test_dynamic_max_trials_never_divides_by_zero() -> None:
    """Test the guard for a consensus set with no inliers.

    fit_unit_slope_ransac() cannot reach this: the stopping criterion is only
    consulted once a candidate has been accepted, which requires at least one
    inlier. The guard is kept because sklearn has it, and because without it
    the ratio would divide by log(1) and yield -inf rather than a trial count.
    """
    assert _dynamic_max_trials(0, 8, 2, 0.99) == float("inf")


# ---------------------------------------------------------------------------
# Parity with scikit-learn
# ---------------------------------------------------------------------------


class _UnitSlopeRegression(RegressorMixin, BaseEstimator):
    """Unit-slope estimator used by the sklearn RANSAC reference."""

    def __init__(self) -> None:
        self.coef_ = None
        self.intercept_ = None

    def fit(self, X, y):  # pylint: disable=invalid-name
        """Fit the regression by estimating the intercept."""
        self.coef_ = 1
        self.intercept_ = np.mean(y - X)
        return self

    def predict(self, X):  # pylint: disable=invalid-name
        """Predict values using the fitted unit-slope regression."""
        return X + self.intercept_

    def score(self, X, y, sample_weight=None):  # pylint: disable=invalid-name
        """Return the negative mean squared error of predictions."""
        squared_errors = (y - self.predict(X)) ** 2

        if sample_weight is None:
            return -float(np.mean(squared_errors))

        return -float(np.average(squared_errors, weights=sample_weight))


def reference_fit(
    abscissa: np.ndarray,
    ordinate: np.ndarray,
) -> tuple[np.ndarray, float, int]:
    """Fit the same problem using sklearn's RANSACRegressor."""
    regressor = RANSACRegressor(
        estimator=_UnitSlopeRegression(),
        min_samples=2,
        max_trials=MAX_TRIALS,
        random_state=RANDOM_STATE,
        residual_threshold=RESIDUAL_THRESHOLD,
    )
    regressor.fit(
        abscissa.reshape(-1, 1),
        ordinate.reshape(-1, 1),
    )

    return (
        regressor.inlier_mask_,
        float(regressor.estimator_.intercept_),
        regressor.n_trials_,
    )


def random_problems(
    num_problems: int,
) -> Iterator[tuple[np.ndarray, np.ndarray]]:
    """Generate contamination-line-shaped problems of varying difficulty."""
    rng = np.random.default_rng(0)

    for _ in range(num_problems):
        num_points = int(rng.integers(8, 120))
        abscissa = np.sort(rng.normal(-6.0, 2.0, num_points))
        ordinate = (
            abscissa
            + rng.choice([0.0, 1.5]) * (rng.random(num_points) < 0.7)
            + rng.normal(0.0, rng.choice([0.02, 0.5]), num_points)
        )
        yield abscissa, ordinate


def test_matches_sklearn_on_random_problems() -> None:
    """Test exact parity with sklearn across random problem instances."""
    for abscissa, ordinate in random_problems(200):
        expected = reference_fit(abscissa, ordinate)
        inlier_mask, intercept, num_trials = fit(abscissa, ordinate)

        # The specialized implementation must preserve all outputs that can
        # affect downstream CroCoDeEL features and contamination results.
        np.testing.assert_array_equal(inlier_mask, expected[0])
        assert intercept == expected[1]
        assert num_trials == expected[2]


def test_matches_sklearn_when_inlier_counts_tie() -> None:
    """Test sklearn's score-based tie-breaking between consensus sets."""
    abscissa = np.array(
        [0.0, 0.0, 1.0, 1.0, 2.0, 2.0, 3.0, 3.0, 4.0, 4.0]
    )
    ordinate = abscissa + np.tile([0.0, 1.0], 5)

    expected = reference_fit(abscissa, ordinate)
    inlier_mask, intercept, num_trials = fit(abscissa, ordinate)

    np.testing.assert_array_equal(inlier_mask, expected[0])
    assert intercept == expected[1]
    assert num_trials == expected[2]
