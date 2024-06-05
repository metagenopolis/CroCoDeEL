from multiprocessing import Pool
from functools import partial
from itertools import product
import logging
from time import perf_counter
from typing import Optional, Final
import numpy as np
import pandas as pd
import tqdm
from crocodeel.conta_event import ContaminationEvent
from crocodeel.ressources import RandomForestModel
from crocodeel.common import (
    select_species_potentially_in_contamination_line,
    get_coefficients_of_potential_contamination_line,
    compute_features
)


def run_search_conta(
    species_ab_table: pd.DataFrame, species_ab_table_2: Optional[pd.DataFrame], nproc=int
) -> list[ContaminationEvent]:
    if species_ab_table_2 is not None:
        all_samples = species_ab_table.columns
        all_samples_2 = species_ab_table_2.columns
        all_sample_pairs = product(all_samples, all_samples_2)
        num_sample_pairs = len(all_samples) * len(all_samples_2)
        species_ab_table = species_ab_table.join(species_ab_table_2, how="outer").fillna(0.0)
    else:
        all_samples = species_ab_table.columns
        all_sample_pairs = product(all_samples, all_samples)
        num_sample_pairs = len(all_samples) ** 2

    contamination_searcher = ContaminationSearcherDriver(
        species_ab_table,
        all_sample_pairs,
        num_sample_pairs,
        nproc)

    start = perf_counter()
    logging.info("Search for contaminations started")
    conta_events = contamination_searcher.search_contamination()
    logging.info("Search completed in %.1f seconds", np.round(perf_counter() - start, 1))

    contaminated_samples = {conta_event.target for conta_event in conta_events}
    logging.info("%d contamination events detected", len(conta_events))
    logging.info("%d samples contaminated", len(contaminated_samples))

    conta_events.sort(key=lambda e: e.rate, reverse=True)

    return conta_events


class ContaminationSearcherWorker:
    PROBABILITY_CUTOFF: Final[float] = 0.5

    def __init__(self, species_ab_table, rf_classifier):
        self.species_ab_table = species_ab_table
        self.rf_classifier = rf_classifier

    def classify_sample_pair(self, sample_pair):
        source, target = sample_pair

        if source == target:
            return ContaminationEvent(source, target)

        # Search a potential contamination line
        (
            not_filtered_data,
            species_potentially_in_contamination_line,
            species_potentially_in_contamination_line_indexes,
        ) = select_species_potentially_in_contamination_line(self.species_ab_table, source, target)

        while True:
            # Not enough species in the potential contamination line
            # no contamination found, exit loop
            if species_potentially_in_contamination_line.shape[0] <= 5:
                return ContaminationEvent(source, target)

            # Run RANSAC to estimate the intercept of the potential contamination line
            # and get inlier and outlier species
            species_inliers, intercept = get_coefficients_of_potential_contamination_line(
                species_potentially_in_contamination_line
            )

            # Not enough inlier species in the potential contamination line
            # no contamination found, exit loop
            if np.sum(species_inliers) <= 5:
                return ContaminationEvent(source, target)

            species_outliers = np.logical_not(species_inliers)
            species_inliers_indexes = species_potentially_in_contamination_line_indexes[species_inliers]
            species_inliers = species_potentially_in_contamination_line[species_inliers]

            # Extract features describing the current sample pair and the potential contamination line
            features = compute_features(intercept, species_inliers, not_filtered_data)
            features = np.array([features])

            # Apply the random forest model with the extracted features
            contamination_probability = self.rf_classifier.predict_proba(features)[0, 1]

            # contamination found, exit loop
            if contamination_probability >= self.PROBABILITY_CUTOFF:
                contamination_rate = np.round(10 ** (-intercept), 4)
                return ContaminationEvent(
                    source,
                    target,
                    rate=contamination_rate,
                    probability=contamination_probability,
                    contamination_specific_species=species_inliers_indexes.tolist(),
                )

            # no contamination found with inliers
            # try with remaining outliers
            species_potentially_in_contamination_line = species_potentially_in_contamination_line[species_outliers]
            species_potentially_in_contamination_line_indexes = species_potentially_in_contamination_line_indexes[
                species_outliers
            ]


class ContaminationSearcherDriver:
    DEFAULT_CHUNKSIZE: Final[int] = 50

    def __init__(
        self, species_ab_table: pd.DataFrame, all_sample_pairs, num_sample_pairs: int, nproc: int = 1
    ):
        self.species_ab_table = species_ab_table
        self.all_sample_pairs = all_sample_pairs
        self.num_sample_pairs = num_sample_pairs
        self.nproc = nproc

    def search_contamination(self):
        rf_classifier = RandomForestModel.load()
        worker = ContaminationSearcherWorker(self.species_ab_table, rf_classifier)

        all_conta_events = []

        with Pool(processes=self.nproc) as pool:
            all_tasks = pool.imap_unordered(
                worker.classify_sample_pair, self.all_sample_pairs, chunksize=self.DEFAULT_CHUNKSIZE
            )
            pbar = partial(
                tqdm.tqdm,
                total=self.num_sample_pairs,
                bar_format="{l_bar}{bar}| {n_fmt}/{total_fmt} sample pairs inspected",
            )
            for conta_event in pbar(all_tasks):
                if conta_event.probability >= ContaminationSearcherWorker.PROBABILITY_CUTOFF:
                    all_conta_events.append(conta_event)

        return all_conta_events
