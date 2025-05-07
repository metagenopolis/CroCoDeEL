from multiprocessing import Pool
import logging
from pathlib import Path
from time import perf_counter
from typing import BinaryIO, TextIO, Final, Iterator
import json
import numpy as np
import pandas as pd
from tqdm import tqdm
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import classification_report
import joblib
from crocodeel.common import (
    select_candidate_species_conta_line,
    search_potential_conta_line,
    compute_conta_line_features
)


def run_train_model(
    species_ab_table: pd.DataFrame,
    model_fh: BinaryIO,
    json_report_fh: TextIO,
    test_size: float,
    ntrees: int,
    rng_seed: int,
    nproc: int,
):
    all_samples = species_ab_table.columns
    source_samples = [sample for sample in all_samples if "source" in sample]
    target_samples = [sample.replace("source", "target") for sample in source_samples]
    all_sample_pairs = zip(source_samples, target_samples)
    num_sample_pairs = len(source_samples)
    is_contaminated = np.char.startswith(source_samples, "conta")

    logging.info(
        "Abundance table contains %i sample pairs: %i contaminated and %i non-contaminated",
        len(is_contaminated),
        sum(is_contaminated),
        sum(~is_contaminated),
    )

    features_computer = FeaturesComputerDriver(
        species_ab_table, all_sample_pairs, num_sample_pairs, nproc
    )
    start = perf_counter()
    logging.info(
        "Computing features using %d process%s...",
        nproc,
        "" if nproc == 1 else "es",
    )
    all_features = features_computer.compute_all_features()
    logging.info(
        "Feature computation completed in %.1f seconds", perf_counter() - start
    )

    # Check sample pairs where feature computation failed
    # due to missing candidate contamination line
    sample_pairs_no_features = np.any(np.isnan(all_features), axis=1)
    num_conta_no_features = sum(is_contaminated[sample_pairs_no_features])
    num_no_conta_no_features = sum(sample_pairs_no_features) - num_conta_no_features
    if num_conta_no_features != 0:
        logging.info(
            "Feature computation failed for %i contaminated sample pairs",
            num_conta_no_features,
        )
    if num_no_conta_no_features != 0:
        logging.info(
            "Feature computation failed for %i non-contaminated sample pairs",
            num_no_conta_no_features,
        )

    # Remove these problematic sample pairs
    is_contaminated = is_contaminated[~sample_pairs_no_features]
    all_features = all_features[~sample_pairs_no_features]
    logging.info(
        "%i sample pairs remain: %i contaminated, %i non-contaminated",
        len(is_contaminated),
        sum(is_contaminated),
        sum(~is_contaminated),
    )

    # Create train and test splits
    features_train, features_test, is_contaminated_train, is_contaminated_test = (
        train_test_split(
            all_features, is_contaminated, test_size=test_size, random_state=rng_seed
        )
    )
    logging.info(
        "Dataset split: %.0f%% for training, %.0f%% for testing",
        100 * (1 - test_size),
        100 * test_size,
    )

    # Train the model
    rf_model = RandomForestClassifier(
        n_estimators=ntrees, n_jobs=nproc, random_state=rng_seed
    )
    start = perf_counter()
    logging.info(
        "Training Random Forest model with %d trees using %d process%s...",
        ntrees,
        nproc,
        "" if nproc == 1 else "es",
    )
    rf_model.fit(features_train, is_contaminated_train)
    logging.info("Training completed in %.1f seconds", perf_counter() - start)

    # Save the model
    joblib.dump(rf_model, model_fh)
    logging.info("Model saved to %s", Path(model_fh.name).resolve())
    model_fh.close()

    # Create and save the performance report
    pred_is_contaminated_train = rf_model.predict(features_train)
    pred_is_contaminated_test = rf_model.predict(features_test)
    performance_report = {
        "train": {
            "classification_report": classification_report(
                is_contaminated_train,
                pred_is_contaminated_train,
                output_dict=True,
            )
        },
        "test": {
            "classification_report": classification_report(
                is_contaminated_test,
                pred_is_contaminated_test,
                output_dict=True,
            )
        }
    }
    json.dump(performance_report, json_report_fh, indent=4)
    logging.info(
        "Model performance report saved to %s", Path(json_report_fh.name).resolve()
    )
    json_report_fh.close()

class Defaults:
    TEST_SIZE: Final[float] = 0.3
    NTREES: Final[int] = 1000
    RNG_SEED: Final[int] = 0

class FeaturesComputerWorker:

    def __init__(self, species_ab_table: pd.DataFrame):
        self.species_ab_table = species_ab_table

    def compute_features_sample_pair(self, sample_pair: tuple[int, tuple[str, str]]):
        sample_pair_id, (source, target) = sample_pair
        # Search a potential contamination line
        (cur_sample_pair_species_ab, candidate_species_conta_line, _) = (
            select_candidate_species_conta_line(self.species_ab_table, source, target)
        )

        if candidate_species_conta_line.shape[0] <= 5:
            return sample_pair_id, None

        candidate_species_inliers, conta_line_offset = search_potential_conta_line(
            candidate_species_conta_line
        )
        candidate_species_inliers = candidate_species_conta_line[
            candidate_species_inliers
        ]
        conta_line_features = compute_conta_line_features(
            conta_line_offset, candidate_species_inliers, cur_sample_pair_species_ab
        )

        return sample_pair_id, conta_line_features


class FeaturesComputerDriver:
    DEFAULT_CHUNKSIZE: Final[int] = 50

    def __init__(
        self,
        species_ab_table: pd.DataFrame,
        all_sample_pairs: Iterator[tuple[str, str]],
        num_sample_pairs: int,
        nproc: int = 1,
    ):
        self.species_ab_table = species_ab_table
        self.all_sample_pairs = list(enumerate(all_sample_pairs))
        self.num_sample_pairs = num_sample_pairs
        self.nproc = nproc

    def compute_all_features(self):
        worker = FeaturesComputerWorker(self.species_ab_table)

        # TODO: replace hardcoded number of features by a constant
        all_features = np.empty((self.num_sample_pairs, 10))

        with Pool(processes=self.nproc) as pool:
            all_tasks = pool.imap_unordered(
                worker.compute_features_sample_pair,
                self.all_sample_pairs,
                chunksize=self.DEFAULT_CHUNKSIZE,
            )
            pbar = tqdm(
                all_tasks,
                total=self.num_sample_pairs,
                leave=False,
                bar_format="{l_bar}{bar}| {n_fmt}/{total_fmt} sample pairs processed",
            )

            for cur_sample_pair_id, cur_sample_pair_features in pbar:
                all_features[cur_sample_pair_id] = cur_sample_pair_features

        return all_features
