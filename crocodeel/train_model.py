"""Train a Random Forest model for cross-sample contamination detection."""

import json
import logging
import re
from multiprocessing import Pool
from time import perf_counter
from typing import Any, BinaryIO, Final, TextIO

import joblib
import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import classification_report
from sklearn.model_selection import train_test_split
from tqdm import tqdm

from crocodeel.conta_features import (ContaminationFeatureExtractor,
                                      ContaminationFeatures)
from crocodeel.exceptions import InputDataError

SamplePair = tuple[str, str]


def run_train_model(
    species_ab_table: pd.DataFrame,
    model_fh: BinaryIO,
    json_report_fh: TextIO,
    test_size: float,
    ntrees: int,
    rng_seed: int,
    nproc: int,
) -> None:
    """Train a Random Forest model and save its performance report."""
    sample_pairs = _reconstruct_sample_pairs(species_ab_table)

    is_contaminated = np.array(
        [source_sample.startswith("conta") for source_sample, _ in sample_pairs],
        dtype=bool,
    )

    logging.info(
        "Abundance table contains %i sample pairs: "
        "%i contaminated and %i non-contaminated",
        len(sample_pairs),
        int(is_contaminated.sum()),
        int((~is_contaminated).sum()),
    )

    all_features = _compute_features(
        species_ab_table,
        sample_pairs,
        nproc,
    )

    all_features, is_contaminated = _filter_invalid_features(
        all_features,
        is_contaminated,
    )

    logging.info(
        "%i sample pairs remain: %i contaminated, %i non-contaminated",
        len(is_contaminated),
        int(is_contaminated.sum()),
        int((~is_contaminated).sum()),
    )

    (
        features_train,
        features_test,
        is_contaminated_train,
        is_contaminated_test,
    ) = train_test_split(
        all_features,
        is_contaminated,
        test_size=test_size,
        random_state=rng_seed,
    )

    logging.info(
        "Dataset split: %.0f%% for training, %.0f%% for testing",
        100 * (1 - test_size),
        100 * test_size,
    )

    rf_model = _train_model(
        features_train,
        is_contaminated_train,
        ntrees,
        rng_seed,
        nproc,
    )

    rf_model.set_params(n_jobs=1)

    joblib.dump(
        rf_model,
        model_fh,
        compress=3,
    )
    logging.info("Model saved to %s", model_fh.name)

    performance_report = _build_performance_report(
        rf_model,
        features_train,
        features_test,
        is_contaminated_train,
        is_contaminated_test,
    )

    json.dump(
        performance_report,
        json_report_fh,
        indent=4,
    )

    logging.info(
        "Model performance report saved to %s",
        json_report_fh.name,
    )


def _reconstruct_sample_pairs(
    species_ab_table: pd.DataFrame,
) -> list[SamplePair]:
    """Reconstruct source-target sample pairs from sample names."""
    sample_names = list(species_ab_table.columns)

    sample_name_pattern = re.compile(r"^(conta_|non_conta_)(source_|target_)case_\d+$")

    invalid_sample_names = [
        name for name in sample_names if sample_name_pattern.fullmatch(name) is None
    ]

    if invalid_sample_names:
        raise InputDataError(
            "The following sample names do not match the expected pattern "
            f"'{sample_name_pattern.pattern}': " + ", ".join(invalid_sample_names)
        )
    source_samples = [sample for sample in sample_names if "_source_" in sample]

    target_samples = [sample for sample in sample_names if "_target_" in sample]

    source_sample_names = set(source_samples)
    target_sample_names = set(target_samples)

    sources_without_targets = [
        source_sample
        for source_sample in source_samples
        if source_sample.replace("_source_", "_target_") not in target_sample_names
    ]

    targets_without_sources = [
        target_sample
        for target_sample in target_samples
        if target_sample.replace("_target_", "_source_") not in source_sample_names
    ]

    if sources_without_targets or targets_without_sources:
        error_messages = []

        if sources_without_targets:
            error_messages.append(
                "source samples without corresponding targets: "
                + ", ".join(sources_without_targets)
            )

        if targets_without_sources:
            error_messages.append(
                "target samples without corresponding sources: "
                + ", ".join(targets_without_sources)
            )

        raise InputDataError("; ".join(error_messages))

    return [
        (
            source_sample,
            source_sample.replace("_source_", "_target_"),
        )
        for source_sample in source_samples
    ]


def _compute_features(
    species_ab_table: pd.DataFrame,
    sample_pairs: list[SamplePair],
    nproc: int,
) -> np.ndarray:
    """Compute contamination features for all sample pairs."""
    start = perf_counter()

    logging.info(
        "Computing features using %d process%s...",
        nproc,
        "" if nproc == 1 else "es",
    )

    features_computer = FeaturesComputerDriver(
        species_ab_table,
        sample_pairs,
        nproc,
    )

    all_features = features_computer.compute_all_features()

    logging.info(
        "Feature computation completed in %.1f seconds",
        perf_counter() - start,
    )

    return all_features


def _filter_invalid_features(
    all_features: np.ndarray,
    is_contaminated: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Remove sample pairs for which feature extraction failed."""
    invalid_features = np.any(
        np.isnan(all_features),
        axis=1,
    )

    num_conta_invalid = int(np.sum(is_contaminated & invalid_features))
    num_no_conta_invalid = int(np.sum(~is_contaminated & invalid_features))

    if num_conta_invalid:
        logging.info(
            "Feature computation failed for %i contaminated sample pairs",
            num_conta_invalid,
        )

    if num_no_conta_invalid:
        logging.info(
            "Feature computation failed for %i non-contaminated sample pairs",
            num_no_conta_invalid,
        )

    return (
        all_features[~invalid_features],
        is_contaminated[~invalid_features],
    )


def _train_model(
    features_train: np.ndarray,
    is_contaminated_train: np.ndarray,
    ntrees: int,
    rng_seed: int,
    nproc: int,
) -> RandomForestClassifier:
    """Train a Random Forest classifier."""
    rf_model = RandomForestClassifier(
        n_estimators=ntrees,
        n_jobs=nproc,
        random_state=rng_seed,
    )

    start = perf_counter()

    logging.info(
        "Training Random Forest model with %d trees using %d process%s...",
        ntrees,
        nproc,
        "" if nproc == 1 else "es",
    )

    rf_model.fit(
        features_train,
        is_contaminated_train,
    )

    logging.info(
        "Training completed in %.1f seconds",
        perf_counter() - start,
    )

    return rf_model


def _build_performance_report(
    rf_model: RandomForestClassifier,
    features_train: np.ndarray,
    features_test: np.ndarray,
    is_contaminated_train: np.ndarray,
    is_contaminated_test: np.ndarray,
) -> dict[str, Any]:
    """Build the model performance report."""
    pred_is_contaminated_train = rf_model.predict(features_train)
    pred_is_contaminated_test = rf_model.predict(features_test)

    return {
        "train": {
            "classification_report": classification_report(
                is_contaminated_train,
                pred_is_contaminated_train,
                output_dict=True,
            ),
        },
        "test": {
            "classification_report": classification_report(
                is_contaminated_test,
                pred_is_contaminated_test,
                output_dict=True,
            ),
        },
    }


class Defaults:
    """Default parameters for Random Forest model training."""

    TEST_SIZE: Final[float] = 0.3
    NTREES: Final[int] = 100
    RNG_SEED: Final[int] = 0


class FeaturesComputerWorker:
    """Compute contamination features for individual sample pairs."""

    def __init__(
        self,
        species_ab_table: pd.DataFrame,
    ) -> None:
        self.feature_extractor = ContaminationFeatureExtractor(species_ab_table)

    def compute_features_sample_pair(
        self,
        sample_pair: SamplePair,
    ) -> np.ndarray | None:
        """Compute features for one sample pair."""
        source, target = sample_pair

        features = self.feature_extractor.extract(
            source,
            target,
        )

        if features is None:
            return None

        return features.values


# Initialize the worker once per process so the abundance table
# is sent to each worker only once, rather than being re-serialized per chunk.
_worker: FeaturesComputerWorker | None = None


def _init_worker(species_ab_table: pd.DataFrame) -> None:
    """Build the per-process feature-computation worker.

    Used as a ``multiprocessing.Pool`` initializer, so it runs exactly once
    per worker process.
    """
    global _worker
    _worker = FeaturesComputerWorker(species_ab_table)


def _compute_features_sample_pair(sample_pair: SamplePair) -> np.ndarray | None:
    """Compute features for one sample pair using the current process's worker."""
    assert _worker is not None
    return _worker.compute_features_sample_pair(sample_pair)


class FeaturesComputerDriver:
    """Compute contamination features using one or more processes."""

    DEFAULT_CHUNKSIZE: Final[int] = 50

    def __init__(
        self,
        species_ab_table: pd.DataFrame,
        sample_pairs: list[SamplePair],
        nproc: int = 1,
    ) -> None:
        self.species_ab_table = species_ab_table
        self.sample_pairs = sample_pairs
        self.num_sample_pairs = len(sample_pairs)
        self.nproc = nproc

    def compute_all_features(self) -> np.ndarray:
        """Compute features for all sample pairs."""
        # Use NaN to mark sample pairs for which feature extraction failed.
        # These rows are removed later by _filter_invalid_features().
        all_features = np.full(
            (
                self.num_sample_pairs,
                ContaminationFeatures.NUM_FEATURES,
            ),
            np.nan,
        )

        with Pool(
            processes=self.nproc,
            initializer=_init_worker,
            initargs=(self.species_ab_table,),
        ) as pool:
            all_tasks = pool.imap(
                _compute_features_sample_pair,
                self.sample_pairs,
                chunksize=self.DEFAULT_CHUNKSIZE,
            )

            pbar = tqdm(
                all_tasks,
                total=self.num_sample_pairs,
                leave=False,
                bar_format=(
                    "{l_bar}{bar}| " "{n_fmt}/{total_fmt} sample pairs processed"
                ),
            )

            for sample_pair_id, sample_pair_features in enumerate(pbar):
                if sample_pair_features is not None:
                    all_features[sample_pair_id] = sample_pair_features

        return all_features
