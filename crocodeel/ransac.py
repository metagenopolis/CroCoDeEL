"""RANSAC fitting specialised for contamination lines.

CroCoDeEL's contamination line is constrained to a slope of one in log10
abundance space, so only the intercept needs to be estimated. This module
implements that specialised case directly instead of using sklearn's
general-purpose ``RANSACRegressor``.
"""

import numpy as np
from sklearn.utils import check_random_state
from sklearn.utils.random import sample_without_replacement

# Smallest representable increment above 1.0. This is used by sklearn's
# RANSAC stopping calculation to avoid taking a logarithm of zero.
_EPSILON = np.spacing(1)


def _dynamic_max_trials(
    num_inliers: int,
    num_samples: int,
    min_samples: int,
    probability: float,
) -> float:
    """Return the number of trials required by RANSAC's stopping criterion.

    This reproduces sklearn's calculation of the number of trials needed to
    obtain, with the requested probability, at least one all-inlier subset.
    """
    inlier_ratio = num_inliers / float(num_samples)
    numerator = max(_EPSILON, 1 - probability)
    denominator = max(_EPSILON, 1 - inlier_ratio**min_samples)

    if numerator == 1:
        return 0

    if denominator == 1:
        return float("inf")

    return abs(float(np.ceil(np.log(numerator) / np.log(denominator))))


def fit_unit_slope_ransac(
    abscissa: np.ndarray,
    ordinate: np.ndarray,
    max_trials: int,
    residual_threshold: float,
    random_state: int,
    min_samples: int = 2,
    stop_probability: float = 0.99,
) -> tuple[np.ndarray, float, int]:
    """Fit ``ordinate = abscissa + intercept`` with RANSAC.

    Because the slope is fixed at one, each RANSAC trial only needs to
    estimate an intercept. On the residuals ``ordinate - abscissa``, this
    amounts to selecting ``min_samples`` residuals and taking their mean.

    The candidate models are generated using the same random sampler as
    sklearn and are independent of the results of previous trials. They can
    therefore all be scored together with vectorised NumPy operations.
    sklearn's sequential model-selection, tie-breaking, and early-stopping
    rules are then replayed over those precomputed scores.

    This preserves the behaviour of the corresponding ``RANSACRegressor``
    while avoiding the overhead of fitting and scoring a general-purpose
    estimator for every trial.

    The number of trials is returned because it is used as one of the
    contamination features.

    Args:
        abscissa: Values of the target sample.
        ordinate: Values of the source sample.
        max_trials: Maximum number of RANSAC trials.
        residual_threshold: Maximum absolute residual for an inlier.
        random_state: Random seed used by sklearn's RANSAC implementation.
        min_samples: Number of observations sampled for each trial.
        stop_probability: Probability used by RANSAC's stopping criterion.

    Returns:
        A tuple containing the best model's inlier mask, fitted intercept,
        and number of trials performed.

    Raises:
        ValueError: If no trial produces a valid consensus set.
    """
    residuals = ordinate - abscissa
    num_samples = residuals.shape[0]

    # Generate the same random subsets that sklearn would generate.
    # Candidate generation does not depend on previous trial scores, so the
    # subsets can safely be generated before any candidates are evaluated.
    random_generator = check_random_state(random_state)
    subsets = np.empty((max_trials, min_samples), dtype=np.intp)

    for trial in range(max_trials):
        subsets[trial] = sample_without_replacement(
            num_samples,
            min_samples,
            random_state=random_generator,
        )

    # For a unit-slope model, the intercept is simply the mean residual of
    # the sampled observations. Score every candidate simultaneously by
    # comparing its residuals with the inlier threshold.
    candidate_offsets = residuals[subsets].mean(axis=1)
    candidate_inliers = (
        np.abs(
            residuals[np.newaxis, :] - candidate_offsets[:, np.newaxis]
        )
        <= residual_threshold
    )
    candidate_num_inliers = candidate_inliers.sum(axis=1)

    # The expensive per-trial calculations are now complete. Keep the
    # remaining loop sequential because sklearn's choice of the best model
    # and its dynamic stopping criterion depend on previous trials.
    best_trial = -1
    best_num_inliers = 1
    best_score = -np.inf

    # This is the current upper bound on the number of trials that need to be
    # replayed. It decreases whenever a better consensus set is found.
    dynamic_max_trials: float = max_trials
    num_trials = 0

    while num_trials < dynamic_max_trials:
        num_trials += 1
        trial = num_trials - 1

        if candidate_num_inliers[trial] < best_num_inliers:
            continue

        # sklearn's unit-slope estimator uses negative mean squared error as
        # its score. Reproduce that score on the current consensus set for
        # the same tie-breaking behaviour.
        inlier_residuals = residuals[candidate_inliers[trial]]
        score = -float(
            np.mean((inlier_residuals - candidate_offsets[trial]) ** 2)
        )

        # Keep the previous model when this candidate has the same number of
        # inliers but a worse score.
        if (
            candidate_num_inliers[trial] == best_num_inliers
            and score < best_score
        ):
            continue

        best_trial = trial
        best_num_inliers = candidate_num_inliers[trial]
        best_score = score

        # A larger consensus set increases the estimated probability that
        # the correct model has already been sampled, allowing RANSAC to stop
        # before reaching max_trials.
        dynamic_max_trials = min(
            dynamic_max_trials,
            _dynamic_max_trials(
                best_num_inliers,
                num_samples,
                min_samples,
                stop_probability,
            ),
        )

    if best_trial < 0:
        raise ValueError("RANSAC could not find a valid consensus set.")

    inlier_mask = candidate_inliers[best_trial]

    # As in sklearn's RANSACRegressor, refit the final model using every
    # observation belonging to the selected consensus set.
    intercept = float(residuals[inlier_mask].mean())

    return inlier_mask, intercept, num_trials
