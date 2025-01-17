from multiprocessing import Pool
from itertools import product
import logging
from time import perf_counter
from typing import Optional, Final
import numpy as np
import pandas as pd
from tqdm import tqdm
from crocodeel.conta_event import ContaminationEvent
from crocodeel.ressources import RandomForestModel
from crocodeel.common import (
    select_candidate_species_conta_line,
    search_potential_conta_line,
    compute_conta_line_features
)


def run_search_conta(
    species_ab_table: pd.DataFrame, species_ab_table_2: Optional[pd.DataFrame], nproc=int
) -> list[ContaminationEvent]:
    if species_ab_table_2 is not None:
        all_samples = species_ab_table.columns
        all_samples_2 = species_ab_table_2.columns
        all_sample_pairs = product(all_samples, all_samples_2)
        num_sample_pairs = len(all_samples) * len(all_samples_2)
        species_ab_table = species_ab_table.join(species_ab_table_2, how="outer").fillna(-np.inf)
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

    def classify_sample_pair(self, sample_pair) -> ContaminationEvent:
        source_sample_name, target_sample_name = sample_pair

        if source_sample_name == target_sample_name:
            return ContaminationEvent(source_sample_name, target_sample_name)

        # Step 1: Selection of candidate species for a contamination line
        (
            cur_sample_pair_species_ab,
            candidate_species_conta_line,
            candidate_species_conta_line_idxs,
        ) = select_candidate_species_conta_line(
            self.species_ab_table, source_sample_name, target_sample_name
        )

        while True:
            # Not enough candidates species for a contamination line
            # no contamination found, exit loop
            if candidate_species_conta_line.shape[0] <= 5:
                return ContaminationEvent(source_sample_name, target_sample_name)

            # Step 2: Search for a potential contamination line
            # Use RANSAC regressor to estimate its offset
            candidate_species_inliers, conta_line_offset = search_potential_conta_line(
                candidate_species_conta_line
            )

            # Not enough inlier species in the potential contamination line
            # no contamination found, exit loop
            if np.sum(candidate_species_inliers) <= 5:
                return ContaminationEvent(source_sample_name, target_sample_name)

            candidate_species_outliers = np.logical_not(candidate_species_inliers)
            candidate_species_inliers_idxs = candidate_species_conta_line_idxs[
                candidate_species_inliers
            ]
            candidate_species_inliers = candidate_species_conta_line[
                candidate_species_inliers
            ]

            # Step 3: Compute features describing the potential contamination line
            conta_line_features = compute_conta_line_features(
                conta_line_offset, candidate_species_inliers, cur_sample_pair_species_ab
            )
            conta_line_features = np.array([conta_line_features])

            # Step 4: Apply the Random Forest model to confirm the contamination event
            conta_probability = self.rf_classifier.predict_proba(conta_line_features)
            conta_probability = conta_probability[0, 1]

            # contamination found, exit loop
            if conta_probability >= self.PROBABILITY_CUTOFF:
                contamination_rate = np.round(10 ** (-conta_line_offset), 4)
                return ContaminationEvent(
                    source_sample_name,
                    target_sample_name,
                    rate=contamination_rate,
                    probability=conta_probability,
                    contamination_specific_species=candidate_species_inliers_idxs.tolist(),
                )

            # no contamination found with inlier species
            # try with remaining outlier species
            # not needed with latest version of the algorithm
            # to be suppressed
            candidate_species_conta_line = candidate_species_conta_line[
                candidate_species_outliers
            ]
            candidate_species_conta_line_idxs = candidate_species_conta_line_idxs[
                candidate_species_outliers
            ]


class ContaminationSearcherDriver:
    DEFAULT_CHUNKSIZE: Final[int] = 50

    def __init__(
        self,
        species_ab_table: pd.DataFrame,
        all_sample_pairs,
        num_sample_pairs: int,
        nproc: int = 1,
    ):
        self.species_ab_table = species_ab_table
        self.all_sample_pairs = all_sample_pairs
        self.num_sample_pairs = num_sample_pairs
        self.nproc = nproc

    def search_contamination(self) -> list[ContaminationEvent]:
        rf_classifier = RandomForestModel.load()
        worker = ContaminationSearcherWorker(self.species_ab_table, rf_classifier)

        all_conta_events = []

        with Pool(processes=self.nproc) as pool:
            all_tasks = pool.imap_unordered(
                worker.classify_sample_pair,
                self.all_sample_pairs,
                chunksize=self.DEFAULT_CHUNKSIZE,
            )
            pbar = tqdm(
                all_tasks,
                total=self.num_sample_pairs,
                leave=False,
                bar_format="{l_bar}{bar}| {n_fmt}/{total_fmt} sample pairs inspected{postfix}",
            )

            for conta_event in pbar:
                if conta_event.probability >= ContaminationSearcherWorker.PROBABILITY_CUTOFF:
                    all_conta_events.append(conta_event)
                    pbar.set_postfix_str(f'{len(all_conta_events)} conta events found')

        return all_conta_events
