import numpy as np
from sklearn.linear_model import RANSACRegressor, LinearRegression
from sklearn.metrics import mean_squared_error
from sklearn.neighbors import NearestNeighbors
from scipy.stats import spearmanr


def _get_mean_abundance_of_most_abundant_species_specific_to_source_sample(
    specific_species_to_source_sample,
    intercept_specific_species_to_source_sample,
    number_of_species=10,
):
    """"""
    specific_species_to_source_sample_sorted = specific_species_to_source_sample[
        specific_species_to_source_sample[:, 1].argsort()[::-1]
    ]
    if specific_species_to_source_sample_sorted.shape[0] == 0:
        return intercept_specific_species_to_source_sample
    return (specific_species_to_source_sample_sorted[:number_of_species, 1]).mean()


def _get_distance_between_mean_abundance_of_specific_species_and_contamination_line(
    specific_species_to_source_sample,
    intercept_specific_species_to_source_sample,
    number_of_species=10,
):
    mean_abundance_of_most_abundant_species_specific_to_source_sample = (
        _get_mean_abundance_of_most_abundant_species_specific_to_source_sample(
            specific_species_to_source_sample, intercept_specific_species_to_source_sample, number_of_species
        )
    )
    return np.abs(
        mean_abundance_of_most_abundant_species_specific_to_source_sample
        - intercept_specific_species_to_source_sample
    )

def _get_distance_between_mean_abundance_of_specific_species_and_contamination_line2(
    specific_species_to_source_sample,
    intercept_specific_species_to_source_sample,
    species_potentially_in_contamination_line_inliers,
    number_of_species=10,
):
    mean_abundance_of_most_abundant_species_specific_to_source_sample = (
        _get_mean_abundance_of_most_abundant_species_specific_to_source_sample(
            specific_species_to_source_sample, intercept_specific_species_to_source_sample, number_of_species
        )
    )
    if (
        mean_abundance_of_most_abundant_species_specific_to_source_sample
        == intercept_specific_species_to_source_sample
    ):
        return 0
    sorted_points = species_potentially_in_contamination_line_inliers[
        species_potentially_in_contamination_line_inliers[:, 1].argsort()
    ]
    if int(0.1 * len(species_potentially_in_contamination_line_inliers)) > 1:
        num_points_to_select = int(0.1 * len(species_potentially_in_contamination_line_inliers))
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

def compute_features(intercept, species_potentially_in_contamination_line_inliers, not_filtered_data):
    # Specific species to the source sample
    specific_species_to_source_sample = not_filtered_data[
        (not_filtered_data[:, 0] == -np.inf) & (not_filtered_data[:, 1] != -np.inf), :
    ]

    # Shared species between source and target samples
    shared_species = not_filtered_data[(not_filtered_data[:, 0] != -np.inf) & (not_filtered_data[:, 1] != -np.inf)]
    number_of_shared_species = shared_species.shape[0]

    #
    number_of_species_in_contamination_line = species_potentially_in_contamination_line_inliers.shape[0]
    ratio_species_in_contamination_line_to_shared_species = (
        number_of_species_in_contamination_line / number_of_shared_species
    )
    number_of_species_above_line = np.sum(shared_species[:, 1] > shared_species[:, 0] + intercept + 0.2)
    ratio_species_above_line_to_shared_species = number_of_species_above_line / number_of_shared_species

    #
    mean_distance_to_nearest_neighbors = _get_mean_distance_to_nearest_neighbors(
        species_potentially_in_contamination_line_inliers
    )
    mean_distance_to_farthest_neighbors = _get_mean_distance_to_farthest_neighbors(
        species_potentially_in_contamination_line_inliers
    )

    #
    pseudo_zero = np.min(not_filtered_data[not_filtered_data[:, 0] != -np.inf, 0]) - 1
    intercept_specific_species_to_source_sample = pseudo_zero + intercept
    distance_between_mean_abundance_of_specific_species_and_contamination_line = (
        _get_distance_between_mean_abundance_of_specific_species_and_contamination_line(
            specific_species_to_source_sample, intercept_specific_species_to_source_sample
        )
    )
    distance_between_mean_abundance_of_specific_species_and_contamination_line2 = (
        _get_distance_between_mean_abundance_of_specific_species_and_contamination_line2(
            specific_species_to_source_sample,
            intercept_specific_species_to_source_sample,
            species_potentially_in_contamination_line_inliers,
        )
    )

    #
    correlation_spearman_all_species = spearmanr(not_filtered_data[:, 0], not_filtered_data[:, 1])[0]

    #
    distances = np.abs(
        species_potentially_in_contamination_line_inliers[:, 1]
        - species_potentially_in_contamination_line_inliers[:, 0]
        - intercept
    ) / np.sqrt(2)
    mean_distance_to_the_contamination_line = distances.mean()

    return (
        ratio_species_in_contamination_line_to_shared_species,
        ratio_species_above_line_to_shared_species,
        number_of_species_in_contamination_line,
        number_of_species_above_line,
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

def get_coefficients_of_potential_contamination_line(species_potentially_in_contamination_line):
    """ """
    ransac = RANSACRegressor(
        estimator=_UnitSlopeRegression(), random_state=42, residual_threshold=0.2
    )
    ransac.fit(
        species_potentially_in_contamination_line[:, [0]],
        species_potentially_in_contamination_line[:, [1]],
    )

    species_inliers = ransac.inlier_mask_
    intercept = ransac.estimator_.coeffs[1]

    return species_inliers, intercept

def select_species_potentially_in_contamination_line(species_ab_table, source_sample_name, target_sample_name):
    """Return"""
    # Select all species or only those in upper triangle
    not_filtered_data = species_ab_table[[target_sample_name, source_sample_name]]
    shared_species_upper_triangle = (
        not_filtered_data[source_sample_name] >= not_filtered_data[target_sample_name]
    ) & (not_filtered_data[target_sample_name] != -np.inf)

    # convert to numpy array
    shared_species_upper_triangle = not_filtered_data[shared_species_upper_triangle]
    shared_species_upper_triangle_indexes = shared_species_upper_triangle.index.values
    not_filtered_data = not_filtered_data.to_numpy()
    shared_species_upper_triangle = shared_species_upper_triangle.to_numpy()

    # search species potentially in the contamination line
    species_potentially_in_contamination_line = []
    for species_id, point in enumerate(shared_species_upper_triangle):
        number_of_points_in_upper_left_triangle = np.sum(
            (shared_species_upper_triangle[:, 0] <= point[0]) & (shared_species_upper_triangle[:, 1] >= point[1])
        )
        number_of_points_in_upper_left_triangle -= 1
        # In potential contamination line if no more than 2 species in upper left triangle
        if number_of_points_in_upper_left_triangle <= 2:
            species_potentially_in_contamination_line.append(species_id)
    species_potentially_in_contamination_line_indexes = shared_species_upper_triangle_indexes[
        species_potentially_in_contamination_line
    ]
    species_potentially_in_contamination_line = shared_species_upper_triangle[
        species_potentially_in_contamination_line
    ]

    return (
        not_filtered_data,
        species_potentially_in_contamination_line,
        species_potentially_in_contamination_line_indexes,
    )
