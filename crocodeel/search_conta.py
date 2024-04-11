import os
os.environ['OMP_NUM_THREADS'] = '1'
from dataclasses import dataclass, field
from multiprocessing import Pool
from functools import partial
from itertools import product
import numpy as np
import pandas as pd
import tqdm
from sklearn.linear_model import RANSACRegressor, LinearRegression
from sklearn.metrics import mean_squared_error
from sklearn.neighbors import NearestNeighbors
from scipy.stats import spearmanr
from conta_event import ContaminationEvent
from rf_model import RandomForestModel

class UnitSlopeRegression(LinearRegression):
    def fit(self, X, y, sample_weight=None):
        self.coeffs = (1, np.mean(y) - np.mean(X))
        return super().fit(X, y)

    def predict(self, X):
        y_hat = X * self.coeffs[0] + self.coeffs[1]
        return y_hat

    def score(self, X, y, sample_weight=None):
        return mean_squared_error(y, self.predict(X))

class ContaminationSearcherWorker:
    UPPER_LEFT_TRIANGLE_LIMIT = 2
    RESIDUAL_THRESHOLD = 0.2
    NUMBER_NEAREST_NEIGHBORS = 5
    NUMBER_FARTHEST_NEIGHBORS = 5
    NUMBER_SPECIFIC_SPECIES_TO_CONSIDER = 10
    PROBABILITY_CUTOFF = 0.5

    def __init__(self, mgs_profiles, rf_classifier):
        self.mgs_profiles = mgs_profiles.div(mgs_profiles.sum(axis=0), axis=1)
        with np.errstate(divide="ignore"):
            self.mgs_profiles = self.mgs_profiles.apply(np.log10)
        # Make sure that species names are strings
        self.mgs_profiles.index = self.mgs_profiles.index.astype(str)
        self.rf_classifier = rf_classifier

    def get_number_of_points_in_upper_left_triangle(self, point, other_points):
        """Return the number of points present in the upper left triangle of a point."""
        number_of_points_in_upper_left_triangle = np.sum((other_points[:, 0] <= point[0]) & (other_points[:, 1] >= point[1]))
        return number_of_points_in_upper_left_triangle-1

    def is_potentially_in_contamination_line(self, point, other_points):
        """Return if the point is potentially in the contamination line or not."""
        number_of_points_in_upper_left_triangle = self.get_number_of_points_in_upper_left_triangle(point, other_points)
        return number_of_points_in_upper_left_triangle <= self.UPPER_LEFT_TRIANGLE_LIMIT

    def select_species_potentially_in_contamination_line(self, source_sample_name, target_sample_name):
        """Return """
        
        # Select all species or only those in upper triangle
        not_filtered_data = self.mgs_profiles[[target_sample_name, source_sample_name]]
        common_species_upper_triangle = (not_filtered_data[source_sample_name] >= not_filtered_data[target_sample_name]) & (not_filtered_data[target_sample_name] != -np.inf)

        # deal with pseudo zero
        minimum_abundance_target_sample = not_filtered_data[not_filtered_data[target_sample_name] != -np.inf][target_sample_name].min()
        pseudo_zero = minimum_abundance_target_sample-1
        not_filtered_data = not_filtered_data.replace(-np.inf, pseudo_zero)

        # convert to numpy array
        common_species_upper_triangle_df= not_filtered_data[common_species_upper_triangle]
        not_filtered_data = not_filtered_data.to_numpy()
        common_species_upper_triangle = common_species_upper_triangle_df.to_numpy()

        # search species potentially in the contamination line
        species_potentially_in_contamination_line = []
        for species_id, point in enumerate(common_species_upper_triangle):
            if self.is_potentially_in_contamination_line(point, common_species_upper_triangle):
                species_potentially_in_contamination_line.append(species_id)
        species_potentially_in_contamination_line_indexes = common_species_upper_triangle_df.index[species_potentially_in_contamination_line]
        species_potentially_in_contamination_line = common_species_upper_triangle[species_potentially_in_contamination_line]

        return not_filtered_data, common_species_upper_triangle, species_potentially_in_contamination_line, species_potentially_in_contamination_line_indexes, pseudo_zero

    def get_coefficients_of_potential_contamination_line(self, species_potentially_in_contamination_line, species_potentially_in_contamination_line_indexes):
        """ """
        X = species_potentially_in_contamination_line[:,0].reshape(-1, 1)
        y = species_potentially_in_contamination_line[:,1].reshape(-1, 1)

        X_indexes = np.array(species_potentially_in_contamination_line_indexes)

        try:
            ransac = RANSACRegressor(estimator=UnitSlopeRegression(), random_state=42, residual_threshold=self.RESIDUAL_THRESHOLD)
            ransac.fit(X, y)

            inliers = ransac.inlier_mask_
            outliers = np.logical_not(inliers)
            _, intercept = ransac.estimator_.coeffs

            species_inliers = np.hstack([X[inliers], y[inliers]])
            species_inliers_indexes = X_indexes[inliers]
            species_outliers = np.hstack([X[outliers], y[outliers]])
            species_outliers_indexes = X_indexes[outliers]

        except ValueError:
            species_inliers, species_outliers, species_inliers_indexes, species_outliers_indexes, intercept = np.empty((0,0)), np.empty((0,0)), np.empty((0,0)), np.empty((0,0)), 0

        return species_inliers, species_outliers, species_inliers_indexes, species_outliers_indexes, intercept

    def select_specific_species_to_source_sample(self, not_filtered_data, pseudo_zero):
        """ """
        specific_species_to_source_sample = not_filtered_data[(not_filtered_data[:, 0] == pseudo_zero) & (not_filtered_data[:, 1] != pseudo_zero), :]
        return specific_species_to_source_sample

    def select_shared_species(self, not_filtered_data, pseudo_zero):
        """ """
        shared_species = not_filtered_data[(not_filtered_data[:,0] != pseudo_zero) & (not_filtered_data[:,1] != pseudo_zero)]
        return shared_species

    def get_number_of_species_in_contamination_line(self, data):
        return data.shape[0]

    def get_intercept_specific_species_to_source_sample(self, intercept, pseudo_zero):
        return pseudo_zero + intercept

    def get_number_of_specific_species_to_source_sample_above_line(self, intercept, specific_species_to_source_sample, pseudo_zero):
        """ """
        intercept_specific_species_to_source_sample = self.get_intercept_specific_species_to_source_sample(intercept, pseudo_zero)
        return np.count_nonzero(specific_species_to_source_sample > intercept_specific_species_to_source_sample)

    def get_spearman_correlation(self, data):
        if data.shape[0] < 2:
            return 0, 0
        return spearmanr(data[:,0], data[:,1])

    def get_number_of_species_above_line(self, shared_species, intercept):
        """"""
        number_of_species_above_line = 0
        for x, y in shared_species:
            if y > (x + intercept+0.2):
                number_of_species_above_line += 1
        return number_of_species_above_line

    def get_mean_abundance_of_most_abundant_species_specific_to_source_sample(self, specific_species_to_source_sample, intercept_specific_species_to_source_sample, number_of_species=NUMBER_SPECIFIC_SPECIES_TO_CONSIDER):
        """"""
        specific_species_to_source_sample_sorted = specific_species_to_source_sample[specific_species_to_source_sample[:,1].argsort()[::-1]]
        if specific_species_to_source_sample_sorted.shape[0] == 0:
            return intercept_specific_species_to_source_sample
        return (specific_species_to_source_sample_sorted[:number_of_species,1]).mean()

    def get_distance_between_mean_abundance_of_specific_species_and_contamination_line(self, specific_species_to_source_sample, intercept_specific_species_to_source_sample, number_of_species=NUMBER_SPECIFIC_SPECIES_TO_CONSIDER):
        mean_abundance_of_most_abundant_species_specific_to_source_sample = self.get_mean_abundance_of_most_abundant_species_specific_to_source_sample(specific_species_to_source_sample, intercept_specific_species_to_source_sample, number_of_species)
        return np.abs(mean_abundance_of_most_abundant_species_specific_to_source_sample - intercept_specific_species_to_source_sample)

    def get_distance_between_mean_abundance_of_specific_species_and_contamination_line2(self, specific_species_to_source_sample, intercept_specific_species_to_source_sample, species_potentially_in_contamination_line_inliers, number_of_species=10):
        mean_abundance_of_most_abundant_species_specific_to_source_sample = self.get_mean_abundance_of_most_abundant_species_specific_to_source_sample(specific_species_to_source_sample, intercept_specific_species_to_source_sample, number_of_species)
        if mean_abundance_of_most_abundant_species_specific_to_source_sample == intercept_specific_species_to_source_sample:
            return 0
        sorted_points = species_potentially_in_contamination_line_inliers[species_potentially_in_contamination_line_inliers[:, 1].argsort()]
        if int(0.1 * len(species_potentially_in_contamination_line_inliers)) > 1:
            num_points_to_select = int(0.1 * len(species_potentially_in_contamination_line_inliers))
        else :
            num_points_to_select = 1
        selected_points = sorted_points[:num_points_to_select]
        return np.abs(mean_abundance_of_most_abundant_species_specific_to_source_sample - np.mean(selected_points[:,1]))

    def get_mean_distance_to_nearest_neighbors(self, data, number_of_neighbors=NUMBER_NEAREST_NEIGHBORS):
        """"""
        if data.shape[0] < number_of_neighbors:
            number_of_neighbors = data.shape[0]
        neighbors_model = NearestNeighbors(n_neighbors=number_of_neighbors)
        neighbors_model.fit(data)
        nearest_neighbors_distances, _ = neighbors_model.kneighbors(data)
        return nearest_neighbors_distances.mean().mean()

    def get_mean_distance_to_farthest_neighbors(self, data, number_of_neighbors=NUMBER_FARTHEST_NEIGHBORS):
        """"""
        neighbors_model = NearestNeighbors(n_neighbors=data.shape[0])
        neighbors_model.fit(data)
        distances, _ = neighbors_model.kneighbors(data)
        sorted_distances = np.sort(distances)
        farthest_neighbors_distances = sorted_distances[:, -number_of_neighbors:]
        return farthest_neighbors_distances.mean().mean()

    def get_mean_distance_to_the_contamination_line(self, data, intercept):
        distances = np.abs(-data[:,0]+data[:,1]-intercept)/np.sqrt(2)
        return distances.mean()

    def get_metrics(self, intercept, species_potentially_in_contamination_line_inliers, not_filtered_data, pseudo_zero):
        # Specific species to the source sample
        specific_species_to_source_sample = self.select_specific_species_to_source_sample(not_filtered_data, pseudo_zero)

        # Shared species between source and target samples
        shared_species = self.select_shared_species(not_filtered_data, pseudo_zero)
        number_of_shared_species = shared_species.shape[0]

        #
        number_of_species_in_contamination_line = self.get_number_of_species_in_contamination_line(species_potentially_in_contamination_line_inliers)
        ratio_species_in_contamination_line_to_shared_species = number_of_species_in_contamination_line/number_of_shared_species
        number_of_species_above_line = self.get_number_of_species_above_line(shared_species, intercept)
        ratio_species_above_line_to_shared_species = number_of_species_above_line/number_of_shared_species

        #
        mean_distance_to_nearest_neighbors = self.get_mean_distance_to_nearest_neighbors(species_potentially_in_contamination_line_inliers)
        mean_distance_to_farthest_neighbors = self.get_mean_distance_to_farthest_neighbors(species_potentially_in_contamination_line_inliers)

        #
        intercept_specific_species_to_source_sample = self.get_intercept_specific_species_to_source_sample(intercept, pseudo_zero)
        distance_between_mean_abundance_of_specific_species_and_contamination_line = self.get_distance_between_mean_abundance_of_specific_species_and_contamination_line(specific_species_to_source_sample, intercept_specific_species_to_source_sample)
        distance_between_mean_abundance_of_specific_species_and_contamination_line2 = self.get_distance_between_mean_abundance_of_specific_species_and_contamination_line2(specific_species_to_source_sample, intercept_specific_species_to_source_sample, species_potentially_in_contamination_line_inliers)

        #
        correlation_spearman_all_species, _ = self.get_spearman_correlation(not_filtered_data)

        #
        mean_distance_to_the_contamination_line = self.get_mean_distance_to_the_contamination_line(species_potentially_in_contamination_line_inliers, intercept)

        return (ratio_species_in_contamination_line_to_shared_species,
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

    def crocodeel(self, species_potentially_in_contamination_line,
                species_potentially_in_contamination_line_indexes,
                not_filtered_data,
                pseudo_zero):
        (species_potentially_in_contamination_line_inliers, species_potentially_in_contamination_line_outliers, inliers_indexes, outliers_indexes,
        intercept) = self.get_coefficients_of_potential_contamination_line(
            species_potentially_in_contamination_line, species_potentially_in_contamination_line_indexes)

        # Not enough species in the contamination line
        if species_potentially_in_contamination_line_inliers.shape[0] < 5:
            contamination_rate = 0
            contamination_probability = 0
            inliers_indexes = np.empty((0, 0))

            return contamination_probability, contamination_rate, inliers_indexes

        X = np.array([self.get_metrics(
            intercept,
            species_potentially_in_contamination_line_inliers,
            not_filtered_data,
            pseudo_zero)])
        contamination_probability = self.rf_classifier.predict_proba(X)[0,1]

        if contamination_probability >= self.PROBABILITY_CUTOFF:
            contamination_rate = np.round(10**(-intercept), 4)
            return contamination_probability, contamination_rate, inliers_indexes
        else:
            return self.crocodeel(species_potentially_in_contamination_line_outliers,outliers_indexes,
                            not_filtered_data,
                            pseudo_zero)
    
    def classify_sample_pair(self, sample_pair):
        source, target = sample_pair
        contamination_probability = 0
        contamination_rate = 0
        inliers_indexes = np.empty((0, 0))

        if source == target:
            return ContaminationEvent(source, target)

        not_filtered_data, common_species_upper_triangle, species_potentially_in_contamination_line, species_potentially_in_contamination_line_indexes, pseudo_zero = self.select_species_potentially_in_contamination_line(source, target)

        if common_species_upper_triangle.shape[0] <= 5:
            return ContaminationEvent(source, target)

        contamination_probability, contamination_rate, inliers_indexes = self.crocodeel(species_potentially_in_contamination_line,
                                                                                    species_potentially_in_contamination_line_indexes,
                                                                                    not_filtered_data,
                                                                                    pseudo_zero)
        return ContaminationEvent(
            source,
            target,
            rate=contamination_rate,
            probability=contamination_probability,
            contamination_specific_species=inliers_indexes.tolist(),
        )

@dataclass
class ContaminationSearcherDriver:
    mgs_profiles: pd.DataFrame
    nproc: int = field(default=1)
    chunksize: int = field(default=50)

    def search_contamination(self):
        all_samples = self.mgs_profiles.columns
        all_sample_pairs = product(all_samples, repeat=2)
        num_sample_pairs = len(all_samples) ** 2

        rf_classifier = RandomForestModel.load()
        worker = ContaminationSearcherWorker(self.mgs_profiles,rf_classifier)

        all_conta_events = []

        with Pool(processes=self.nproc) as pool:
            all_tasks = pool.imap_unordered(
                worker.classify_sample_pair, all_sample_pairs, chunksize=self.chunksize
            )
            pbar = partial(
                tqdm.tqdm,
                total=num_sample_pairs,
                bar_format="{l_bar}{bar}| {n_fmt}/{total_fmt} sample pairs inspected",
            )

            for conta_event in pbar(all_tasks):
                if conta_event.probability >= ContaminationSearcherWorker.PROBABILITY_CUTOFF:
                    all_conta_events.append(conta_event)

        return all_conta_events
