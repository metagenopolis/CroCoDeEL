import numpy as np
from sklearn.linear_model import RANSACRegressor, LinearRegression
from sklearn.metrics import mean_squared_error
from sklearn.neighbors import NearestNeighbors
from scipy.stats import spearmanr


def _get_mean_ab_top_source_specific_species(
    source_specific_species_ab,
    num_species=10,
):
    """"""
    assert source_specific_species_ab.size != 0
    source_specific_species_ab_sorted = source_specific_species_ab[
        source_specific_species_ab[:, 1].argsort()[::-1]
    ]

    return (source_specific_species_ab_sorted[:num_species, 1]).mean()


def _get_diff_mean_ab_top10_source_specices_vs_ab_cutoff1(
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
    return np.abs(
        mean_ab_top10_source_specific_species - ab_cutoff1
    )


def _get_diff_mean_ab_top10_source_specices_vs_ab_cutoff2(
    mean_ab_top10_source_specific_species,
    candidate_species_inliers,
):
    # Select the 10% least abundant inlier species in the source
    sorted_points = candidate_species_inliers[candidate_species_inliers[:, 1].argsort()]
    num_points_to_select = max(int((0.1 * len(candidate_species_inliers))),1)
    selected_points = sorted_points[:num_points_to_select]

    # Mean abundance of the 10% least abundant inlier species in the source
    # Named 'c2' in the paper
    ab_cutoff2 = np.mean(selected_points[:, 1])

    # '|m-c2|' in the paper
    return np.abs(
        mean_ab_top10_source_specific_species - ab_cutoff2
    )


def _get_mean_distance_to_nearest_neighbors(data, num_neighbors=5):
    """"""
    if data.shape[0] < num_neighbors:
        num_neighbors = data.shape[0]
    neighbors_model = NearestNeighbors(n_neighbors=num_neighbors)
    neighbors_model.fit(data)
    nearest_neighbors_distances, _ = neighbors_model.kneighbors(data)
    return nearest_neighbors_distances.mean().mean()

def _get_mean_distance_to_farthest_neighbors(data, num_neighbors=5):
    """"""
    neighbors_model = NearestNeighbors(n_neighbors=data.shape[0])
    neighbors_model.fit(data)
    distances, _ = neighbors_model.kneighbors(data)
    sorted_distances = np.sort(distances)
    farthest_neighbors_distances = sorted_distances[:, -num_neighbors:]
    return farthest_neighbors_distances.mean().mean()


def compute_conta_line_features(
    conta_line_offset, candidate_species_inliers, cur_sample_pair_species_ab
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
    mean_distance_to_nearest_neighbors = _get_mean_distance_to_nearest_neighbors(
        candidate_species_inliers
    )

    # Feature 4
    mean_distance_to_farthest_neighbors = _get_mean_distance_to_farthest_neighbors(
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
        diff_mean_ab_top10_source_specices_vs_ab_cutoff1 = 0
        diff_mean_ab_top10_source_specices_vs_ab_cutoff2 = 0
    else:
        # Mean abundance in the source of the 10 most abundant species in the source
        # Named 'm' in the paper
        mean_ab_top10_source_specific_species = (
            _get_mean_ab_top_source_specific_species(
                source_specific_species_ab
            )
        )

        # Feature 9
        diff_mean_ab_top10_source_specices_vs_ab_cutoff1 = (
            _get_diff_mean_ab_top10_source_specices_vs_ab_cutoff1(
                mean_ab_top10_source_specific_species,
                cur_sample_pair_species_ab,
                conta_line_offset,
            )
        )

        # Feature 10
        diff_mean_ab_top10_source_specices_vs_ab_cutoff2 = (
            _get_diff_mean_ab_top10_source_specices_vs_ab_cutoff2(
                mean_ab_top10_source_specific_species,
                candidate_species_inliers,
            )
        )

    return (
        ratio_species_conta_line_to_shared_species,
        ratio_species_above_line_to_shared_species,
        num_species_conta_line,
        num_species_above_conta_line,
        spearman_corr_all_species,
        mean_distance_to_the_contamination_line,
        mean_distance_to_nearest_neighbors,
        mean_distance_to_farthest_neighbors,
        diff_mean_ab_top10_source_specices_vs_ab_cutoff1,
        diff_mean_ab_top10_source_specices_vs_ab_cutoff2,
    )


class _UnitSlopeRegression(LinearRegression):
    def fit(self, X, y, sample_weight=None):
        self.coeffs = (1, np.mean(y) - np.mean(X))
        return super().fit(X, y)

    def predict(self, X):
        y_hat = X * self.coeffs[0] + self.coeffs[1]
        return y_hat

    def score(self, X, y, sample_weight=None):
        return mean_squared_error(y, self.predict(X))

def search_potential_conta_line(candidate_species_conta_line):
    """ """
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

def select_candidate_species_conta_line(species_ab_table, source_sample_name, target_sample_name):
    # Select abundance of all species from the current sample pair
    cur_sample_pair_species_ab = species_ab_table[
        [target_sample_name, source_sample_name]
    ]
    # Select species shared by both samples but more abundant in the source
    shared_species_upper_triangle_ab = (
        cur_sample_pair_species_ab[source_sample_name]
        >= cur_sample_pair_species_ab[target_sample_name]
    ) & (cur_sample_pair_species_ab[target_sample_name] != -np.inf)

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
