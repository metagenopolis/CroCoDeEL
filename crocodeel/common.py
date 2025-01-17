import numpy as np
from sklearn.linear_model import RANSACRegressor, LinearRegression
from sklearn.metrics import mean_squared_error
from sklearn.neighbors import NearestNeighbors
from scipy.stats import spearmanr


def _get_mean_abundance_of_most_abundant_species_specific_to_source_sample(
    source_specific_species_ab,
    intercept_specific_species_to_source_sample,
    number_of_species=10,
):
    """"""
    source_specific_species_ab_sorted = source_specific_species_ab[
        source_specific_species_ab[:, 1].argsort()[::-1]
    ]
    if source_specific_species_ab_sorted.shape[0] == 0:
        return intercept_specific_species_to_source_sample
    return (source_specific_species_ab_sorted[:number_of_species, 1]).mean()


def _get_distance_between_mean_abundance_of_specific_species_and_contamination_line(
    source_specific_species_ab,
    intercept_specific_species_to_source_sample,
    num_species=10,
):
    mean_abundance_of_most_abundant_species_specific_to_source_sample = (
        _get_mean_abundance_of_most_abundant_species_specific_to_source_sample(
            source_specific_species_ab, intercept_specific_species_to_source_sample, num_species
        )
    )
    return np.abs(
        mean_abundance_of_most_abundant_species_specific_to_source_sample
        - intercept_specific_species_to_source_sample
    )

def _get_distance_between_mean_abundance_of_specific_species_and_contamination_line2(
    source_specific_species_ab,
    intercept_specific_species_to_source_sample,
    candidate_species_inliers,
    num_species=10,
):
    mean_abundance_of_most_abundant_species_specific_to_source_sample = (
        _get_mean_abundance_of_most_abundant_species_specific_to_source_sample(
            source_specific_species_ab, intercept_specific_species_to_source_sample, num_species
        )
    )
    if (
        mean_abundance_of_most_abundant_species_specific_to_source_sample
        == intercept_specific_species_to_source_sample
    ):
        return 0
    sorted_points = candidate_species_inliers[
        candidate_species_inliers[:, 1].argsort()
    ]
    if int(0.1 * len(candidate_species_inliers)) > 1:
        num_points_to_select = int(0.1 * len(candidate_species_inliers))
    else:
        num_points_to_select = 1
    selected_points = sorted_points[:num_points_to_select]
    return np.abs(
        mean_abundance_of_most_abundant_species_specific_to_source_sample - np.mean(selected_points[:, 1])
    )

def _get_mean_distance_to_nearest_neighbors(data, number_of_neighbors=5):
    """"""
    if data.shape[0] < number_of_neighbors:
        number_of_neighbors = data.shape[0]
    neighbors_model = NearestNeighbors(n_neighbors=number_of_neighbors)
    neighbors_model.fit(data)
    nearest_neighbors_distances, _ = neighbors_model.kneighbors(data)
    return nearest_neighbors_distances.mean().mean()

def _get_mean_distance_to_farthest_neighbors(data, number_of_neighbors=5):
    """"""
    neighbors_model = NearestNeighbors(n_neighbors=data.shape[0])
    neighbors_model.fit(data)
    distances, _ = neighbors_model.kneighbors(data)
    sorted_distances = np.sort(distances)
    farthest_neighbors_distances = sorted_distances[:, -number_of_neighbors:]
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

    #
    num_species_conta_line = candidate_species_inliers.shape[0]
    ratio_species_conta_line_to_shared_species = (
        num_species_conta_line / num_shared_species
    )
    num_species_above_conta_line = np.sum(
        shared_species_ab[:, 1] > shared_species_ab[:, 0] + conta_line_offset + 0.2
    )
    ratio_species_above_line_to_shared_species = (
        num_species_above_conta_line / num_shared_species
    )

    #
    mean_distance_to_nearest_neighbors = _get_mean_distance_to_nearest_neighbors(
        candidate_species_inliers
    )
    mean_distance_to_farthest_neighbors = _get_mean_distance_to_farthest_neighbors(
        candidate_species_inliers
    )

    #
    pseudo_zero = np.min(cur_sample_pair_species_ab[cur_sample_pair_species_ab[:, 0] != -np.inf, 0]) - 1
    intercept_specific_species_to_source_sample = pseudo_zero + conta_line_offset
    distance_between_mean_abundance_of_specific_species_and_contamination_line = (
        _get_distance_between_mean_abundance_of_specific_species_and_contamination_line(
            source_specific_species_ab, intercept_specific_species_to_source_sample
        )
    )
    distance_between_mean_abundance_of_specific_species_and_contamination_line2 = (
        _get_distance_between_mean_abundance_of_specific_species_and_contamination_line2(
            source_specific_species_ab,
            intercept_specific_species_to_source_sample,
            candidate_species_inliers,
        )
    )

    #
    correlation_spearman_all_species = spearmanr(
        cur_sample_pair_species_ab[:, 0], cur_sample_pair_species_ab[:, 1]
    )[0]

    #
    distances = np.abs(
        candidate_species_inliers[:, 1]
        - candidate_species_inliers[:, 0]
        - conta_line_offset
    ) / np.sqrt(2)
    mean_distance_to_the_contamination_line = distances.mean()

    return (
        ratio_species_conta_line_to_shared_species,
        ratio_species_above_line_to_shared_species,
        num_species_conta_line,
        num_species_above_conta_line,
        correlation_spearman_all_species,
        mean_distance_to_the_contamination_line,
        mean_distance_to_nearest_neighbors,
        mean_distance_to_farthest_neighbors,
        distance_between_mean_abundance_of_specific_species_and_contamination_line,
        distance_between_mean_abundance_of_specific_species_and_contamination_line2,
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
