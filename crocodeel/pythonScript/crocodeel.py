import numpy as np
import pandas as pd
import argparse
import sys
import os
import joblib
from sklearn.linear_model import RANSACRegressor, LinearRegression
from sklearn.metrics import mean_squared_error
from sklearn.neighbors import NearestNeighbors
from scipy.stats import spearmanr

UPPER_LEFT_TRIANGLE_LIMIT = 2
RESIDUAL_THRESHOLD = 0.2
NUMBER_NEAREST_NEIGHBORS = 5
NUMBER_FARTHEST_NEIGHBORS = 5
NUMBER_SPECIFIC_SPECIES_TO_CONSIDER = 10

def get_arguments(args=sys.argv[1:]):
    parser = argparse.ArgumentParser(description="Parses command.")
    parser.add_argument("-i", "--mgs_profiles_path", type=str, help="MGS profiles path.", required=True)
    parser.add_argument("-s", "--source", type=str, help="Name of the source sample.", required=True)
    arguments = parser.parse_args(args)
    return arguments

class CustomLinearRegression(LinearRegression):
    def fit(self, X, y):
        self.coeffs = (1, np.mean(y)-np.mean(X))
        return super().fit(X, y)
    
    def predict(self, X):
        y_hat = X * self.coeffs[0] + self.coeffs[1]
        return y_hat
    
    def score(self, X, y):
        return mean_squared_error(y, self.predict(X))

def get_number_of_points_in_upper_left_triangle(point, other_points):
    """Return the number of points present in the upper left triangle of a point."""
    number_of_points_in_upper_left_triangle = 0
    for other_point in other_points:
        if other_point[0] <= point[0] and other_point[1] >= point[1]:
            number_of_points_in_upper_left_triangle += 1
    return number_of_points_in_upper_left_triangle

def is_potentially_in_contamination_line(point, other_points, limit):
    """Return if the point is potentially in the contamination line or not."""
    number_of_points_in_upper_left_triangle = get_number_of_points_in_upper_left_triangle(point, other_points)
    if number_of_points_in_upper_left_triangle <= limit:
        return True
    
def select_species_potentially_in_contamination_line(df, source_sample_name, target_sample_name, limit):
    """Return """
    source_sample = df.loc[:, source_sample_name]
    target_sample = df.loc[:, target_sample_name]
    data = pd.concat([target_sample, source_sample], axis=1)
    minimum_abundance_target_sample = data[data[target_sample_name] != 0][target_sample_name].min()
    
    # Not filtered data, log
    not_filtered_data = data.replace(0, minimum_abundance_target_sample/10)
    not_filtered_data = np.log10(not_filtered_data)
    not_filtered_data = np.hstack([not_filtered_data[target_sample_name], not_filtered_data[source_sample_name]]).reshape(2,-1).T
    
    # Filtered data, only common species in the upper triangle, log
    common_species_upper_triangle_df = data[(data[source_sample_name] != 0) & (data[target_sample_name] != 0)]
    common_species_upper_triangle_df = common_species_upper_triangle_df[common_species_upper_triangle_df[source_sample_name] >= common_species_upper_triangle_df[target_sample_name]] 
    common_species_upper_triangle_df = common_species_upper_triangle_df.replace(0, minimum_abundance_target_sample/10)
    common_species_upper_triangle_df = np.log10(common_species_upper_triangle_df)
    common_species_upper_triangle = np.hstack([common_species_upper_triangle_df[target_sample_name], common_species_upper_triangle_df[source_sample_name]]).reshape(2,-1).T
    
    # Filtered data, only species potentially in the contamination line, log
    species_potentially_in_contamination_line = []
    species_potentially_in_contamination_line_indexes = []
    for i, point in enumerate(common_species_upper_triangle):
        other_points = np.delete(common_species_upper_triangle, i, axis=0) 
        if is_potentially_in_contamination_line(point, other_points, limit):
            species_potentially_in_contamination_line.append(point)  
            species_potentially_in_contamination_line_indexes.append(common_species_upper_triangle_df.index[i])
    species_potentially_in_contamination_line = np.array(species_potentially_in_contamination_line)
    
    return not_filtered_data, common_species_upper_triangle, species_potentially_in_contamination_line, species_potentially_in_contamination_line_indexes, minimum_abundance_target_sample

def get_coefficients_of_potential_contamination_line(species_potentially_in_contamination_line, species_potentially_in_contamination_line_indexes, estimator=CustomLinearRegression()):
    """ """
    X = species_potentially_in_contamination_line[:,0].reshape(-1, 1)
    y = species_potentially_in_contamination_line[:,1].reshape(-1, 1)
    
    X_indexes = np.array(species_potentially_in_contamination_line_indexes)

    try: 
        ransac = RANSACRegressor(estimator=estimator, random_state=42, residual_threshold=RESIDUAL_THRESHOLD)
        ransac.fit(X, y)

        inliers = ransac.inlier_mask_
        outliers = np.logical_not(inliers)
        slope, intercept = ransac.estimator_.coeffs
        
        species_inliers = np.vstack([X[inliers], y[inliers]]).reshape(2,-1).T
        species_inliers_indexes = X_indexes[inliers]
        species_outliers = np.vstack([X[outliers], y[outliers]]).reshape(2,-1).T
        species_outliers_indexes = X_indexes[outliers]
        
    except ValueError as e:
        species_inliers, species_outliers, species_inliers_indexes, species_outliers_indexes, slope, intercept = np.empty((0,0)), np.empty((0,0)), np.empty((0,0)), np.empty((0,0)), 0, 0
        
    return species_inliers, species_outliers, species_inliers_indexes, species_outliers_indexes, slope, intercept

def select_specific_species_to_source_sample(not_filtered_data, minimum_abundance_target_sample):
    """ """
    specific_species_to_source_sample = not_filtered_data[(not_filtered_data[:, 0] == np.log10(minimum_abundance_target_sample/10)) & (not_filtered_data[:, 1] != np.log10(minimum_abundance_target_sample/10)), :]
    return specific_species_to_source_sample

def select_shared_species(not_filtered_data, minimum_abundance_target_sample):
    """ """
    shared_species = not_filtered_data[(not_filtered_data[:,0] != np.log10(minimum_abundance_target_sample/10)) & (not_filtered_data[:,1] != np.log10(minimum_abundance_target_sample/10))]
    return shared_species

def get_number_of_species_in_contamination_line(data):
    return data.shape[0]

def get_intercept_specific_species_to_source_sample(slope, intercept, minimum_abundance_target_sample):
    return slope * np.log10(minimum_abundance_target_sample/10) + intercept

def get_number_of_specific_species_to_source_sample_above_line(slope, intercept, specific_species_to_source_sample, minimum_abundance_target_sample):
    """ """
    intercept_specific_species_to_source_sample = get_intercept_specific_species_to_source_sample(slope, intercept, minimum_abundance_target_sample)
    return np.count_nonzero(specific_species_to_source_sample > intercept_specific_species_to_source_sample)

def get_spearman_correlation(data):
    if data.shape[0] < 2:
        return 0, 0
    return spearmanr(data[:,0], data[:,1])

def get_number_of_species_above_line(shared_species, slope, intercept):
    """"""
    number_of_species_above_line = 0
    for x, y in shared_species:
        if y > (slope * x + intercept+0.2):
            number_of_species_above_line += 1
    return number_of_species_above_line

def get_mean_abundance_of_most_abundant_species_specific_to_source_sample(specific_species_to_source_sample, intercept_specific_species_to_source_sample, number_of_species=NUMBER_SPECIFIC_SPECIES_TO_CONSIDER):
    """"""
    specific_species_to_source_sample_sorted = specific_species_to_source_sample[specific_species_to_source_sample[:,1].argsort()[::-1]]
    if specific_species_to_source_sample_sorted.shape[0] == 0:
        return intercept_specific_species_to_source_sample
    else :
        return (specific_species_to_source_sample_sorted[:number_of_species,1]).mean()
        
def get_distance_between_mean_abundance_of_specific_species_and_contamination_line(specific_species_to_source_sample, intercept_specific_species_to_source_sample, number_of_species=NUMBER_SPECIFIC_SPECIES_TO_CONSIDER):
    mean_abundance_of_most_abundant_species_specific_to_source_sample = get_mean_abundance_of_most_abundant_species_specific_to_source_sample(specific_species_to_source_sample, intercept_specific_species_to_source_sample, number_of_species)
    return np.abs(mean_abundance_of_most_abundant_species_specific_to_source_sample - intercept_specific_species_to_source_sample)

def get_distance_between_mean_abundance_of_specific_species_and_contamination_line2(specific_species_to_source_sample, intercept_specific_species_to_source_sample, species_potentially_in_contamination_line_inliers, number_of_species=10):
    mean_abundance_of_most_abundant_species_specific_to_source_sample = get_mean_abundance_of_most_abundant_species_specific_to_source_sample(specific_species_to_source_sample, intercept_specific_species_to_source_sample, number_of_species)
    if mean_abundance_of_most_abundant_species_specific_to_source_sample == intercept_specific_species_to_source_sample:
        return 0
    sorted_points = species_potentially_in_contamination_line_inliers[species_potentially_in_contamination_line_inliers[:, 1].argsort()]
    if int(0.1 * len(species_potentially_in_contamination_line_inliers)) > 1:
        num_points_to_select = int(0.1 * len(species_potentially_in_contamination_line_inliers))
    else : 
        num_points_to_select = 1
    selected_points = sorted_points[:num_points_to_select]
    return np.abs(mean_abundance_of_most_abundant_species_specific_to_source_sample - np.mean(selected_points[:,1]))

def get_mean_distance_to_nearest_neighbors(data, number_of_neighbors=NUMBER_NEAREST_NEIGHBORS):
    """"""
    if data.shape[0] < number_of_neighbors:
        number_of_neighbors = data.shape[0]
    neighbors_model = NearestNeighbors(n_neighbors=number_of_neighbors)
    neighbors_model.fit(data)
    nearest_neighbors_distances, indexes = neighbors_model.kneighbors(data)
    return nearest_neighbors_distances.mean().mean()

def get_mean_distance_to_farthest_neighbors(data, number_of_neighbors=NUMBER_FARTHEST_NEIGHBORS):
    """"""
    neighbors_model = NearestNeighbors(n_neighbors=data.shape[0])
    neighbors_model.fit(data)
    distances, indexes = neighbors_model.kneighbors(data)
    sorted_distances = np.sort(distances)
    farthest_neighbors_distances = sorted_distances[:, -number_of_neighbors:]
    return farthest_neighbors_distances.mean().mean()

def get_mean_distance_to_the_contamination_line(data, intercept):
    distances = np.abs(-data[:,0]+data[:,1]-intercept)/np.sqrt(2)
    return distances.mean()
    
def get_metrics(slope, intercept, species_potentially_in_contamination_line_inliers, not_filtered_data, minimum_abundance_target_sample):
    # Specific species to the source sample
    specific_species_to_source_sample = select_specific_species_to_source_sample(not_filtered_data, minimum_abundance_target_sample)
    
    # Shared species between source and target samples
    shared_species = select_shared_species(not_filtered_data, minimum_abundance_target_sample)
    number_of_shared_species = shared_species.shape[0]
    
    # 
    number_of_species_in_contamination_line = get_number_of_species_in_contamination_line(species_potentially_in_contamination_line_inliers)
    ratio_species_in_contamination_line_to_shared_species = number_of_species_in_contamination_line/number_of_shared_species
    number_of_species_above_line = get_number_of_species_above_line(shared_species, slope, intercept)
    ratio_species_above_line_to_shared_species = number_of_species_above_line/number_of_shared_species
        
    # 
    mean_distance_to_nearest_neighbors = get_mean_distance_to_nearest_neighbors(species_potentially_in_contamination_line_inliers)
    mean_distance_to_farthest_neighbors = get_mean_distance_to_farthest_neighbors(species_potentially_in_contamination_line_inliers)
    
    #
    intercept_specific_species_to_source_sample = get_intercept_specific_species_to_source_sample(slope, intercept, minimum_abundance_target_sample)
    distance_between_mean_abundance_of_specific_species_and_contamination_line = get_distance_between_mean_abundance_of_specific_species_and_contamination_line(specific_species_to_source_sample, intercept_specific_species_to_source_sample)
    distance_between_mean_abundance_of_specific_species_and_contamination_line2 = get_distance_between_mean_abundance_of_specific_species_and_contamination_line2(specific_species_to_source_sample, intercept_specific_species_to_source_sample, species_potentially_in_contamination_line_inliers)

    #
    correlation_spearman_all_species, p_values_all_species = get_spearman_correlation(not_filtered_data)
    
    #
    mean_distance_to_the_contamination_line = get_mean_distance_to_the_contamination_line(species_potentially_in_contamination_line_inliers, intercept)

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

def crocodeel(species_potentially_in_contamination_line,
              species_potentially_in_contamination_line_indexes,
              not_filtered_data,
              minimum_abundance_target_sample, 
              source_sample_name, 
              target_sample_name,
              rf_classifier):
    (species_potentially_in_contamination_line_inliers, species_potentially_in_contamination_line_outliers, inliers_indexes, outliers_indexes,
     slope, intercept) = get_coefficients_of_potential_contamination_line(
         species_potentially_in_contamination_line, species_potentially_in_contamination_line_indexes)
    
    if species_potentially_in_contamination_line_inliers.shape[0] != 0:
        X = np.array([get_metrics(
            slope,
            intercept,
            species_potentially_in_contamination_line_inliers,
            not_filtered_data,
            minimum_abundance_target_sample)])
        y_prediction = rf_classifier.predict(X)[0]
        probabilities = rf_classifier.predict_proba(X)
        contamination_probability = probabilities[:, 1][0]
        
        if y_prediction == 1 :
            contamination_rate = np.round(1/(10**intercept), 4)
            if species_potentially_in_contamination_line_inliers.shape[0] < 5: # Predicted as contaminated but not enough species in the contamination line
                y_prediction = 0
                contamination_rate = 0
                contamination_probability = 0
                inliers_indexes = np.empty((0, 0))
            
            return y_prediction, contamination_probability, contamination_rate, inliers_indexes
        
        else :
            if species_potentially_in_contamination_line_outliers.shape[0] >= 5: 
                return crocodeel(species_potentially_in_contamination_line_outliers, 
                                 outliers_indexes,
                                 not_filtered_data,
                                 minimum_abundance_target_sample, 
                                 source_sample_name, 
                                 target_sample_name,
                                 rf_classifier)
            else:
                contamination_probability = 0
                contamination_rate = 0
                inliers_indexes = np.empty((0, 0))
                return y_prediction, contamination_probability, contamination_rate, inliers_indexes
    else:
        y_prediction = 0
        contamination_probability = 0
        contamination_rate = 0
        inliers_indexes = np.empty((0, 0))
        return y_prediction, contamination_probability, contamination_rate, inliers_indexes
                 

def main():
    arguments = get_arguments()
    dict = {}

    mgs_profiles = pd.read_csv(arguments.mgs_profiles_path, sep='\t', header=0, index_col=0)
    mgs_profiles = mgs_profiles.div(mgs_profiles.sum(axis=0), axis=1)

    rf_classifier = joblib.load(os.path.dirname(__file__) + '/../models/crocodeel_last_version.joblib')

    samples = mgs_profiles.columns.tolist()
    other_samples = samples.copy()
    other_samples.remove(arguments.source)
    
    for target in other_samples:
        not_filtered_data, common_species_upper_triangle, species_potentially_in_contamination_line, species_potentially_in_contamination_line_indexes, minimum_abundance_target_sample = select_species_potentially_in_contamination_line(
            mgs_profiles, arguments.source, target, limit=UPPER_LEFT_TRIANGLE_LIMIT)
        
        if common_species_upper_triangle.shape[0] > 5:
            y_prediction, contamination_probability, contamination_rate, inliers_indexes = crocodeel(species_potentially_in_contamination_line,
                                                                                    species_potentially_in_contamination_line_indexes,
                                                                                    not_filtered_data,
                                                                                    minimum_abundance_target_sample, 
                                                                                    arguments.source, 
                                                                                    target,
                                                                                    rf_classifier)
        else :
            y_prediction = 0
            contamination_probability = 0
            contamination_rate = 0
            inliers_indexes = np.empty((0, 0))
        
        inliers_indexes_list = ",".join(map(str, inliers_indexes))
        dict[target] = [y_prediction, contamination_probability, contamination_rate, inliers_indexes_list]
        
    with open(arguments.source + '.txt', 'a') as fichier:
        fichier.write("source\ttarget\tcontamination\tprobability\trate\tmsp_list\n")
        for key, value in dict.items():
            if value[0] == 1:
                fichier.write(arguments.source + '\t' + key + '\t' + "\t".join(map(str, value)) + "\n")
        
if __name__ == '__main__':
    main()
