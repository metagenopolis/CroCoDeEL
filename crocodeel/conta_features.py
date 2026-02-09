from typing import Final, ClassVar, Optional
from dataclasses import dataclass
import numpy as np
import pandas as pd
from sklearn.linear_model import RANSACRegressor, LinearRegression
from sklearn.metrics import mean_squared_error
from sklearn.neighbors import NearestNeighbors
from scipy.stats import spearmanr

@dataclass
class ContaminationFeatures:
    NUM_FEATURES: ClassVar[int] = 10

    values: np.ndarray
    conta_line_offset: float
    conta_line_species: list[str]

class _UnitSlopeRegression(LinearRegression):
    def fit(self, X, y, sample_weight=None):
        self.coeffs = (1, np.mean(y) - np.mean(X))
        return super().fit(X, y)

    def predict(self, X):
        y_hat = X * self.coeffs[0] + self.coeffs[1]
        return y_hat

    def score(self, X, y, sample_weight=None):
        return mean_squared_error(y, self.predict(X))

class ContaminationFeatureExtractor:
    CONTA_LINE_MIN_NUM_SPECIES: Final[int] = 6

    def __init__(self, species_ab_table: pd.DataFrame):
        self.species_ab_table = species_ab_table

    def extract(self, source: str, target: str) -> Optional[ContaminationFeatures]:
        # Step 1: Selection of candidate species for a contamination line
        (
            cur_sample_pair_species_ab,
            candidate_species_conta_line,
            candidate_species_conta_line_idxs,
        ) = self._select_candidate_species_conta_line(
            source, target
        )

        # Not enough candidates species for a contamination line
        if candidate_species_conta_line.shape[0] < self.CONTA_LINE_MIN_NUM_SPECIES:
            return None

        # Step 2: Search for a potential contamination line
        # Use RANSAC regressor to estimate its offset
        candidate_species_inliers, conta_line_offset = self._search_potential_conta_line(
            candidate_species_conta_line
        )

        # Not enough inlier species in the potential contamination line
        # no contamination found, exit
        if np.sum(candidate_species_inliers) < self.CONTA_LINE_MIN_NUM_SPECIES:
            return None

        candidate_species_inliers_idxs = candidate_species_conta_line_idxs[
            candidate_species_inliers
        ]
        candidate_species_inliers = candidate_species_conta_line[
            candidate_species_inliers
        ]

        # Step 3: Compute features describing the potential contamination line
        conta_line_features = self._compute_features(
            conta_line_offset, candidate_species_inliers, cur_sample_pair_species_ab
        )

        return ContaminationFeatures(
            values=conta_line_features,
            conta_line_offset=conta_line_offset,
            conta_line_species = candidate_species_inliers_idxs.tolist())

    def _select_candidate_species_conta_line(self, source: str, target: str):
        # Select abundance of all species from the current sample pair
        cur_sample_pair_species_ab = self.species_ab_table[
            [target, source]
        ]
        # Select species shared by both samples but more abundant in the source
        shared_species_upper_triangle_ab = (
            cur_sample_pair_species_ab[source]
            >= cur_sample_pair_species_ab[target]
        ) & (cur_sample_pair_species_ab[target] != -np.inf)

        # convert to numpy array
        shared_species_upper_triangle_ab = cur_sample_pair_species_ab[
            shared_species_upper_triangle_ab
        ]
        shared_species_upper_triangle_idxs = shared_species_upper_triangle_ab.index.values
        cur_sample_pair_species_ab = cur_sample_pair_species_ab.to_numpy()
        shared_species_upper_triangle_ab = shared_species_upper_triangle_ab.to_numpy()

        # Test each species in the upper left triangle as a candidate
        candidates_species = []
        for species_id, species_ab in enumerate(shared_species_upper_triangle_ab):
            num_species_in_upper_left_quadrant = np.sum(
                (shared_species_upper_triangle_ab[:, 0] <= species_ab[0])
                & (shared_species_upper_triangle_ab[:, 1] >= species_ab[1])
            )
            num_species_in_upper_left_quadrant -= 1
            # The species is selected if no more than 2 other species are in upper left quadrant
            if num_species_in_upper_left_quadrant <= 2:
                candidates_species.append(species_id)
        candidates_species_idxs = shared_species_upper_triangle_idxs[candidates_species]
        candidates_species = shared_species_upper_triangle_ab[candidates_species]

        return (
            cur_sample_pair_species_ab,
            candidates_species,
            candidates_species_idxs,
        )

    def _search_potential_conta_line(self, candidate_species_conta_line):
        ransac = RANSACRegressor(
            estimator=_UnitSlopeRegression(), random_state=42, residual_threshold=0.2
        )
        ransac.fit(
            candidate_species_conta_line[:, [0]],
            candidate_species_conta_line[:, [1]],
        )

        candidate_species_inliers = ransac.inlier_mask_
        conta_line_offset = ransac.estimator_.coeffs[1]

        return candidate_species_inliers, conta_line_offset

    def _compute_features(
        self, conta_line_offset, candidate_species_inliers, cur_sample_pair_species_ab
    ):
        # Species detected only in the source sample
        source_specific_species_ab = cur_sample_pair_species_ab[
            (cur_sample_pair_species_ab[:, 0] == -np.inf)
            & (cur_sample_pair_species_ab[:, 1] != -np.inf),
            :,
        ]

        # Shared species between source and target samples
        shared_species_ab = cur_sample_pair_species_ab[
            (cur_sample_pair_species_ab[:, 0] != -np.inf)
            & (cur_sample_pair_species_ab[:, 1] != -np.inf)
        ]
        num_shared_species = shared_species_ab.shape[0]

        # Feature 1
        num_species_conta_line = candidate_species_inliers.shape[0]

        # Feature 2
        ratio_species_conta_line_to_shared_species = (
            num_species_conta_line / num_shared_species
        )

        # Feature 8
        num_species_above_conta_line = np.sum(
            shared_species_ab[:, 1] > shared_species_ab[:, 0] + conta_line_offset + 0.2
        )

        # Feature 7
        ratio_species_above_line_to_shared_species = (
            num_species_above_conta_line / num_shared_species
        )

        # Feature 3
        mean_distance_to_nearest_neighbors = self._get_mean_distance_to_nearest_neighbors(
            candidate_species_inliers
        )

        # Feature 4
        mean_distance_to_farthest_neighbors = self._get_mean_distance_to_farthest_neighbors(
            candidate_species_inliers
        )

        # Feature 5
        spearman_corr_all_species = spearmanr(
            cur_sample_pair_species_ab[:, 0], cur_sample_pair_species_ab[:, 1]
        )[0]

        # Feature 6
        distances = np.abs(
            candidate_species_inliers[:, 1]
            - candidate_species_inliers[:, 0]
            - conta_line_offset
        ) / np.sqrt(2)
        mean_distance_to_the_contamination_line = distances.mean()

        # Case where features 9 and 10 cannot be calculated
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

            # Feature 9
            diff_mean_ab_top10_source_species_vs_ab_cutoff1 = (
                self._get_diff_mean_ab_top10_source_species_vs_ab_cutoff1(
                    mean_ab_top10_source_specific_species,
                    cur_sample_pair_species_ab,
                    conta_line_offset,
                )
            )

            # Feature 10
            diff_mean_ab_top10_source_species_vs_ab_cutoff2 = (
                self._get_diff_mean_ab_top10_source_species_vs_ab_cutoff2(
                    mean_ab_top10_source_specific_species,
                    candidate_species_inliers,
                )
            )

        return np.array(
            [
                ratio_species_conta_line_to_shared_species,
                ratio_species_above_line_to_shared_species,
                num_species_conta_line,
                num_species_above_conta_line,
                spearman_corr_all_species,
                mean_distance_to_the_contamination_line,
                mean_distance_to_nearest_neighbors,
                mean_distance_to_farthest_neighbors,
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
        cur_sample_pair_species_ab,
        conta_line_offset
    ):
        # Define a pseudo zero
        min_non_zero = np.min(
            cur_sample_pair_species_ab[cur_sample_pair_species_ab[:, 0] != -np.inf, 0]
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

    def _get_mean_distance_to_farthest_neighbors(self, data, num_neighbors=5):
        """"""
        neighbors_model = NearestNeighbors(n_neighbors=data.shape[0])
        neighbors_model.fit(data)
        distances, _ = neighbors_model.kneighbors(data)
        farthest_neighbors_distances = distances[:, -num_neighbors:]
        return farthest_neighbors_distances.mean().mean()
