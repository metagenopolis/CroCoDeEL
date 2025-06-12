from multiprocessing import Pool
import logging
from pathlib import Path
from time import perf_counter
from typing import BinaryIO, TextIO, Final, Optional
import json
import re
import sys
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
from sklearn.model_selection import cross_val_score, StratifiedKFold
from sklearn.metrics import make_scorer, precision_score, recall_score, f1_score, accuracy_score


def run_train_model(
    species_ab_table: pd.DataFrame,
    model_fh: BinaryIO,
    json_report_fh: TextIO,
    test_size: float,
    ntrees: int,
    rng_seed: int,
    nproc: int,
) -> None:
    sample_pairs = _reconstruct_sample_pairs(species_ab_table)
    is_contaminated = np.array(
        [source_sample.startswith("conta") for source_sample, _ in sample_pairs]
    )

    logging.info(
        "Abundance table contains %i sample pairs: %i contaminated and %i non-contaminated",
        len(sample_pairs),
        sum(is_contaminated),
        sum(~is_contaminated),
    )

    features_computer = FeaturesComputerDriver(
        species_ab_table, sample_pairs, nproc
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

    # Cross-validation setup
    rf_model = RandomForestClassifier(
        n_estimators=ntrees, n_jobs=nproc, random_state=rng_seed
    )
    cv = StratifiedKFold(n_splits=10, shuffle=True, random_state=rng_seed)

    logging.info("Performing 10-fold cross-validation...")

    # Evaluate with multiple metrics
    precision_scores = cross_val_score(rf_model, all_features, is_contaminated, cv=cv,
                                       scoring=make_scorer(precision_score))
    recall_scores = cross_val_score(rf_model, all_features, is_contaminated, cv=cv,
                                    scoring=make_scorer(recall_score))
    f1_scores = cross_val_score(rf_model, all_features, is_contaminated, cv=cv,
                                scoring=make_scorer(f1_score))
    accuracy_scores = cross_val_score(rf_model, all_features, is_contaminated, cv=cv,
                                scoring=make_scorer(accuracy_score))

    logging.info("10-fold CV results:")
    logging.info("  Precision: %.3f ± %.3f", np.mean(precision_scores), np.std(precision_scores))
    logging.info("  Recall:    %.3f ± %.3f", np.mean(recall_scores), np.std(recall_scores))
    logging.info("  F1-score:  %.3f ± %.3f", np.mean(f1_scores), np.std(f1_scores))
    logging.info("  Accuracy:  %.3f ± %.3f", np.mean(accuracy_scores), np.std(accuracy_scores))


def _reconstruct_sample_pairs(
    species_ab_table: pd.DataFrame,
) -> list[tuple[str, str]] :
    all_samples = species_ab_table.columns

    # Identify invalid sample names
    sample_name_pattern = re.compile(r"^(conta_|non_conta_)(source_|target_)case_\d+$")
    invalid_sample_names = [
        name for name in all_samples if not sample_name_pattern.match(name)
    ]
    if invalid_sample_names:
        logging.error(
            "The following sample names do not match the expected pattern '%s': %s",
            sample_name_pattern.pattern,
            ", ".join(invalid_sample_names),
        )
        sys.exit(1)

    # Identify source and target samples
    source_samples = [sample for sample in all_samples if "source" in sample]
    target_samples = [sample for sample in all_samples if "target" in sample]

    # Check if each source has a corresponding target
    sources_without_targets = [
        source_sample
        for source_sample in source_samples
        if source_sample.replace("source", "target") not in target_samples
    ]
    if sources_without_targets:
        logging.error(
            "The following source samples have no corresponding targets: %s",
            ", ".join(sources_without_targets),
        )
    # Check if each target has a corresponding source
    targets_without_sources = [
        target_sample
        for target_sample in target_samples
        if target_sample.replace("target", "source") not in source_samples
    ]
    if targets_without_sources:
        logging.error(
            "The following target samples have no corresponding sources: %s",
            ", ".join(targets_without_sources),
        )
    if sources_without_targets or targets_without_sources:
        sys.exit(1)

    # Everything is OK
    # We can reconstruct sample pairs
    sample_pairs = [
        (source_sample, source_sample.replace("source", "target"))
        for source_sample in source_samples
    ]

    return sample_pairs


class Defaults:
    TEST_SIZE: Final[float] = 0.3
    NTREES: Final[int] = 1000
    RNG_SEED: Final[int] = 0


class FeaturesComputerWorker:

    def __init__(self, species_ab_table: pd.DataFrame) -> None:
        self.species_ab_table = species_ab_table

    def compute_features_sample_pair(
        self, sample_pair: tuple[str, str]
    ) -> Optional[np.ndarray]:
        source, target = sample_pair
        # Search a potential contamination line
        (cur_sample_pair_species_ab, candidate_species_conta_line, _) = (
            select_candidate_species_conta_line(self.species_ab_table, source, target)
        )

        if candidate_species_conta_line.shape[0] <= 5:
            return None

        candidate_species_inliers, conta_line_offset = search_potential_conta_line(
            candidate_species_conta_line
        )
        candidate_species_inliers = candidate_species_conta_line[
            candidate_species_inliers
        ]
        conta_line_features = compute_conta_line_features(
            conta_line_offset, candidate_species_inliers, cur_sample_pair_species_ab
        )

        return conta_line_features


class FeaturesComputerDriver:
    DEFAULT_CHUNKSIZE: Final[int] = 50

    def __init__(
        self,
        species_ab_table: pd.DataFrame,
        sample_pairs: list[tuple[str, str]],
        nproc: int = 1,
    ) -> None:
        self.species_ab_table = species_ab_table
        self.sample_pairs = sample_pairs
        self.num_sample_pairs = len(sample_pairs)
        self.nproc = nproc

    def compute_all_features(self) -> np.ndarray:
        worker = FeaturesComputerWorker(self.species_ab_table)

        # TODO: replace hardcoded number of features by a constant
        all_features = np.empty((self.num_sample_pairs, 10))

        with Pool(processes=self.nproc) as pool:
            all_tasks = pool.imap(
                worker.compute_features_sample_pair,
                self.sample_pairs,
                chunksize=self.DEFAULT_CHUNKSIZE,
            )
            pbar = tqdm(
                all_tasks,
                total=self.num_sample_pairs,
                leave=False,
                bar_format="{l_bar}{bar}| {n_fmt}/{total_fmt} sample pairs processed",
            )

            for cur_sample_pair_id, cur_sample_pair_features in enumerate(pbar):
                all_features[cur_sample_pair_id] = cur_sample_pair_features

        return all_features
