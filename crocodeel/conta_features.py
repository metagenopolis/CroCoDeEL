from typing import Final, ClassVar, Optional
from dataclasses import dataclass
import numpy as np
import pandas as pd
from sklearn.linear_model import RANSACRegressor
from sklearn.base import RegressorMixin, BaseEstimator
from sklearn.metrics import mean_squared_error
from sklearn.neighbors import NearestNeighbors

@dataclass
class ContaminationFeatures:
    NUM_FEATURES: ClassVar[int] = 6

    values: np.ndarray
    conta_line_offset: float
    conta_line_species: list[str]

class _UnitSlopeRegression(RegressorMixin, BaseEstimator):
    def __init__(self):
        self.coef_ = None
        self.intercept_ = None

    def fit(self, X, y):
        self.coef_ = 1
        self.intercept_ = np.mean(y - X)
        return self

    def predict(self, X):
        return X + self.intercept_

    def score(self, X, y, sample_weight=None):
        return -mean_squared_error(y, self.predict(X), sample_weight=sample_weight)

class ContaminationFeatureExtractor:
    CONTA_LINE_MIN_NUM_SPECIES: Final[int] = 8
    UPPER_LEFT_QUADRANT_MAX_NUM_SPECIES: Final[int] = 2
    RANSAC_MAX_TRIALS: Final[int] = 30
    RANSAC_RANDOM_STATE: Final[int] = 42
    RANSAC_RESIDUAL_THRESHOLD: Final[float] = 0.2

    def __init__(self, species_ab_table: pd.DataFrame):
        self.species_ab_table = species_ab_table

        self.ransac = RANSACRegressor(
            estimator=_UnitSlopeRegression(),
            min_samples=2,
            max_trials=self.RANSAC_MAX_TRIALS,
            random_state=self.RANSAC_RANDOM_STATE,
            residual_threshold=self.RANSAC_RESIDUAL_THRESHOLD,
        )

    def extract(self, source: str, target: str) -> Optional[ContaminationFeatures]:
        # Select abundance of all species from the current sample pair
        sample_pair_species_ab = self.species_ab_table[
            [target, source]
        ].to_numpy()

        # Step 1: Selection of candidate species for a contamination line
        conta_line_candidate_species_ab = self._get_conta_line_candidate_species(
            sample_pair_species_ab
        )

        # Not enough candidates species for a contamination line
        if conta_line_candidate_species_ab.shape[0] < self.CONTA_LINE_MIN_NUM_SPECIES:
            return None

        # Step 2: Search for a potential contamination line
        # Use RANSAC regressor to estimate its offset
        try:
            conta_line_species_ab, conta_line_offset = self._estimate_conta_line_offset(
                conta_line_candidate_species_ab
            )
        except ValueError:
            # RANSAC could not find a valid consensus set
            # no contamination found, exit
            return None

        # Not enough species in the contamination line
        # no contamination found, exit
        if conta_line_species_ab.shape[0] < self.CONTA_LINE_MIN_NUM_SPECIES:
            return None

        # Step 3: Compute features describing the potential contamination line
        conta_line_features = self._compute_features(
            sample_pair_species_ab, conta_line_species_ab, conta_line_offset,
        )

        # Refine contamination line
        refined_conta_line_species_names = (
            self._refine_conta_line(sample_pair_species_ab, conta_line_species_ab)
        )

        return ContaminationFeatures(
            values=conta_line_features,
            conta_line_offset=conta_line_offset,
            conta_line_species = refined_conta_line_species_names)

    def _get_conta_line_candidate_species(
        self, sample_pair_species_ab: np.ndarray
    ) -> np.ndarray:
        # Select species shared by both samples but more abundant in the source
        mask_upper_triangle = (
            sample_pair_species_ab[:, 1] >= sample_pair_species_ab[:, 0]
        ) & (sample_pair_species_ab[:, 0] != -np.inf)

        upper_triangle_species_ab = sample_pair_species_ab[mask_upper_triangle, :]
        # Get candidate species for a contamination line
        x = upper_triangle_species_ab[:, 0]  # target
        y = upper_triangle_species_ab[:, 1]  # source

        # Broadcasting: Compare every element to every other element
        # Shapes: (N, 1) vs (1, N) -> Result (N, N) matrices
        is_left = x[:, np.newaxis] <= x[np.newaxis, :]   # j is to the left of i
        is_upper = y[:, np.newaxis] >= y[np.newaxis, :]  # j is above i

        # Sum columns to get count for each species i
        # Subtract 1 because a species is always in its own upper-left quadrant
        upper_left_quadrant_num_species = np.sum(is_left & is_upper, axis=0) - 1

        # Filter candidates
        mask_candidate_species = (
            upper_left_quadrant_num_species <= self.UPPER_LEFT_QUADRANT_MAX_NUM_SPECIES
        )
        candidates_species_ab = upper_triangle_species_ab[mask_candidate_species]

        return candidates_species_ab

    def _refine_conta_line(
        self, sample_pair_species_ab: np.ndarray, conta_line_species_ab: np.ndarray
    ) -> list[str]:
        # Select species shared by both samples but more abundant in the source
        mask_upper_triangle = (
            sample_pair_species_ab[:, 1] >= sample_pair_species_ab[:, 0]
        ) & (sample_pair_species_ab[:, 0] != -np.inf)

        conta_line_species_ab_ratio = conta_line_species_ab[:, 1]- conta_line_species_ab[:, 0]

        q1_conta_line_offset = np.quantile(conta_line_species_ab_ratio, 0.25)
        q3_conta_line_offset = np.quantile(conta_line_species_ab_ratio, 0.75)
        iqr_conta_line_offset = q3_conta_line_offset - q1_conta_line_offset

        lower_bound_conta_line_offset = max(
            q1_conta_line_offset - 1.5 * iqr_conta_line_offset,
            np.min(conta_line_species_ab_ratio),
        )

        upper_bound_conta_line_offset = min(
            q3_conta_line_offset + 1.5 * iqr_conta_line_offset,
            np.max(conta_line_species_ab_ratio),
        )

        upper_triangle_species_ab_ratio = (
            sample_pair_species_ab[mask_upper_triangle, 1]
            - sample_pair_species_ab[mask_upper_triangle, 0]
        )

        mask_refined_conta_line_species = (
            upper_triangle_species_ab_ratio >= lower_bound_conta_line_offset
        ) & (upper_triangle_species_ab_ratio <= upper_bound_conta_line_offset)

        refined_conta_line_species_names = self.species_ab_table.index.values[
            mask_upper_triangle
        ][mask_refined_conta_line_species].tolist()

        return refined_conta_line_species_names

    def _estimate_conta_line_offset(
        self, conta_line_candidate_species_ab: np.ndarray
    ) -> tuple[np.ndarray, float]:
        self.ransac.fit(
            conta_line_candidate_species_ab[:, [0]],  # target
            conta_line_candidate_species_ab[:, [1]],  # source
        )

        conta_line_species_ab = conta_line_candidate_species_ab[
            self.ransac.inlier_mask_
        ]
        conta_line_offset = self.ransac.estimator_.intercept_

        return conta_line_species_ab, conta_line_offset

    def _compute_features(
        self,
        sample_pair_species_ab: np.ndarray,
        conta_line_species_ab: np.ndarray,
        conta_line_offset: float,
    ) -> np.ndarray:
        # Species detected only in the source sample
        source_specific_species_ab = sample_pair_species_ab[
            (sample_pair_species_ab[:, 0] == -np.inf)
            & (sample_pair_species_ab[:, 1] != -np.inf),
            :,
        ]

        num_species_conta_line = conta_line_species_ab.shape[0]

        mean_distance_to_nearest_neighbors = self._get_mean_distance_to_nearest_neighbors(
            conta_line_species_ab
        )

        distances = np.abs(
            conta_line_species_ab[:, 1]
            - conta_line_species_ab[:, 0]
            - conta_line_offset
        ) / np.sqrt(2)
        mean_distance_to_the_contamination_line = distances.mean()

        # Case where features cannot be calculated
        if source_specific_species_ab.size == 0:
            diff_mean_ab_top10_source_species_vs_ab_cutoff1 = 0
            diff_mean_ab_top10_source_species_vs_ab_cutoff2 = 0
        else:
            # Mean abundance in the source of the 10 most abundant species in the source
            # Named 'm' in the paper
            mean_ab_top10_source_specific_species = (
                self._get_mean_ab_top_source_specific_species(
                    source_specific_species_ab
                )
            )

            diff_mean_ab_top10_source_species_vs_ab_cutoff1 = (
                self._get_diff_mean_ab_top10_source_species_vs_ab_cutoff1(
                    mean_ab_top10_source_specific_species,
                    sample_pair_species_ab,
                    conta_line_offset,
                )
            )

            diff_mean_ab_top10_source_species_vs_ab_cutoff2 = (
                self._get_diff_mean_ab_top10_source_species_vs_ab_cutoff2(
                    mean_ab_top10_source_specific_species,
                    conta_line_species_ab,
                )
            )

        return np.array(
            [
                self.ransac.n_trials_,
                num_species_conta_line,
                mean_distance_to_the_contamination_line,
                mean_distance_to_nearest_neighbors,
                diff_mean_ab_top10_source_species_vs_ab_cutoff1,
                diff_mean_ab_top10_source_species_vs_ab_cutoff2,
            ]
        )

    def _get_mean_ab_top_source_specific_species(
        self,
        source_specific_species_ab,
        num_species=10,
    ):
        """"""
        assert source_specific_species_ab.size != 0

        # Abundance of source specific species in the source
        values = source_specific_species_ab[:, 1]

        if len(values) <= num_species:
            return values.mean()

        top_n_idx = np.argpartition(-values, num_species)[:num_species]
        return values[top_n_idx].mean()

    def _get_diff_mean_ab_top10_source_species_vs_ab_cutoff1(
        self,
        mean_ab_top10_source_specific_species,
        sample_pair_species_ab,
        conta_line_offset
    ):
        # Define a pseudo zero
        min_non_zero = np.min(
            sample_pair_species_ab[sample_pair_species_ab[:, 0] != -np.inf, 0]
        )
        pseudo_zero = min_non_zero - 1

        # Named 'c1' in the paper
        ab_cutoff1 = pseudo_zero + conta_line_offset

        # '|m-c1|' in the paper
        return np.abs(mean_ab_top10_source_specific_species - ab_cutoff1)

    def _get_diff_mean_ab_top10_source_species_vs_ab_cutoff2(
        self,
        mean_ab_top10_source_specific_species,
        candidate_species_inliers,
    ):
        # Abundance of inliers in the source
        values = candidate_species_inliers[:, 1]

        # Select the 10% least abundant species
        num_species = max(int((0.1 * len(values))), 1)
        idx = np.argpartition(values, num_species)[:num_species]

        # Mean abundance of these species
        # Named 'c2' in the paper
        ab_cutoff2 = values[idx].mean()

        # '|m-c2|' in the paper
        return np.abs(mean_ab_top10_source_specific_species - ab_cutoff2)

    def _get_mean_distance_to_nearest_neighbors(self, data, num_neighbors=5):
        """"""
        if data.shape[0] < num_neighbors:
            num_neighbors = data.shape[0]
        neighbors_model = NearestNeighbors(n_neighbors=num_neighbors)
        neighbors_model.fit(data)
        nearest_neighbors_distances, _ = neighbors_model.kneighbors(data)
        return nearest_neighbors_distances.mean().mean()
