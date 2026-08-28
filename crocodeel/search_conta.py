"""Search for cross-sample contamination in metagenomic abundance data."""

import importlib.resources
import logging
import warnings
from itertools import batched, product
from multiprocessing import Pool
from pathlib import Path
from time import perf_counter
from typing import BinaryIO, Final, Iterator, Sequence

import joblib
import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.exceptions import InconsistentVersionWarning  # type: ignore[attr-defined]
from tqdm import tqdm

from crocodeel.conta_event import ContaminationEvent
from crocodeel.conta_features import (
    ContaminationFeatureExtractor,
    ContaminationFeatures,
)
from crocodeel.exceptions import InputDataError

SamplePair = tuple[str, str]


def run_search_conta(
    species_ab_table: pd.DataFrame,
    species_ab_table_2: pd.DataFrame | None,
    rf_model_fh: BinaryIO,
    probability_cutoff: float,
    rate_cutoff: float,
    nproc: int,
) -> list[ContaminationEvent]:
    """Search for cross-sample contamination."""
    rf_model = _load_rf_model(rf_model_fh)

    (
        species_ab_table,
        all_sample_pairs,
        num_sample_pairs,
    ) = _prepare_search(
        species_ab_table,
        species_ab_table_2,
    )

    start = perf_counter()

    logging.info(
        "Searching contamination using %d process%s...",
        nproc,
        "" if nproc == 1 else "es",
    )

    conta_events = ContaminationSearcherDriver(
        species_ab_table=species_ab_table,
        all_sample_pairs=all_sample_pairs,
        num_sample_pairs=num_sample_pairs,
        rf_model=rf_model,
        probability_cutoff=probability_cutoff,
        rate_cutoff=rate_cutoff,
        nproc=nproc,
    ).search_contamination()

    logging.info(
        "Search completed in %.1f seconds",
        perf_counter() - start,
    )

    # Sort by decreasing contamination rate, then by source and target
    # to ensure deterministic output when rates are identical.
    conta_events.sort(
        key=lambda event: (-event.rate, event.source, event.target),
    )

    _log_search_results(conta_events)

    return conta_events


def _load_rf_model(
    rf_model_fh: BinaryIO,
) -> RandomForestClassifier:
    """Load and validate the Random Forest model."""
    with warnings.catch_warnings():
        warnings.filterwarnings(
            action="ignore",
            category=InconsistentVersionWarning,
        )
        rf_model = joblib.load(rf_model_fh)

    if not isinstance(rf_model, RandomForestClassifier):
        raise InputDataError(
            f"The model file {rf_model_fh.name} is not an "
            "sklearn RandomForestClassifier."
        )

    if rf_model.n_features_in_ != ContaminationFeatures.NUM_FEATURES:
        raise InputDataError(
            f"The Random Forest model {rf_model_fh.name} expects "
            f"{rf_model.n_features_in_} features, but "
            f"{ContaminationFeatures.NUM_FEATURES} are required."
        )

    return rf_model


def _prepare_search(
    species_ab_table: pd.DataFrame,
    species_ab_table_2: pd.DataFrame | None,
) -> tuple[
    pd.DataFrame,
    Iterator[SamplePair],
    int,
]:
    """Prepare the abundance table and sample pairs for contamination search."""
    if species_ab_table_2 is None:
        sample_names = species_ab_table.columns
        sample_pairs = product(sample_names, sample_names)
        num_sample_pairs = len(sample_names) ** 2

        return (
            species_ab_table,
            sample_pairs,
            num_sample_pairs,
        )

    sample_names_1 = species_ab_table.columns
    sample_names_2 = species_ab_table_2.columns

    sample_pairs = product(
        sample_names_1,
        sample_names_2,
    )
    num_sample_pairs = len(sample_names_1) * len(sample_names_2)

    combined_table = species_ab_table.join(
        species_ab_table_2,
        how="outer",
    ).fillna(-np.inf)

    return (
        combined_table,
        sample_pairs,
        num_sample_pairs,
    )


def _log_search_results(
    conta_events: list[ContaminationEvent],
) -> None:
    """Log a summary of contamination search results."""
    contaminated_samples = {event.target for event in conta_events}

    num_events = len(conta_events)
    num_samples = len(contaminated_samples)

    logging.info(
        "%d contamination event%s detected",
        num_events,
        "" if num_events == 1 else "s",
    )

    logging.info(
        "%d sample%s contaminated",
        num_samples,
        "" if num_samples == 1 else "s",
    )


def _passes_cutoffs(
    conta_event: ContaminationEvent,
    probability_cutoff: float,
    rate_cutoff: float,
) -> bool:
    """Return whether an event passes the probability and rate cutoffs."""
    return (
        conta_event.probability >= probability_cutoff
        and conta_event.rate > rate_cutoff
    )


class Defaults:
    """Default parameters for contamination search."""

    MODEL_FILE: Final[Path] = Path(
        str(
            importlib.resources.files().joinpath(
                "models",
                "crocodeel_rf_Feb2026_2.joblib",
            )
        )
    )

    PROBABILITY_CUTOFF: Final[float] = 0.5
    RATE_CUTOFF: Final[float] = 0.0


class ContaminationSearcherWorker:
    """Classify batches of sample pairs for contamination."""

    def __init__(
        self,
        species_ab_table: pd.DataFrame,
        rf_classifier: RandomForestClassifier,
    ) -> None:
        self.feature_extractor = ContaminationFeatureExtractor(species_ab_table)
        self.rf_classifier = rf_classifier

    def classify_sample_pairs(
        self,
        sample_pairs: Sequence[SamplePair],
    ) -> list[ContaminationEvent]:
        """Classify a batch of sample pairs.

        Features are extracted individually because feature extraction depends
        on the source and target samples. The resulting feature matrix is then
        passed to the Random Forest in a single call to avoid paying the
        prediction overhead once per sample pair.
        """
        candidates: list[tuple[SamplePair, ContaminationFeatures]] = []

        for sample_pair in sample_pairs:
            source, target = sample_pair

            if source == target:
                continue

            features = self.feature_extractor.extract(
                source,
                target,
            )

            if features is not None:
                candidates.append((sample_pair, features))

        if not candidates:
            return []

        feature_matrix = np.vstack(
            [features.values for _, features in candidates],
        )
        conta_probabilities = self.rf_classifier.predict_proba(
            feature_matrix,
        )[:, 1]

        return [
            ContaminationEvent(
                source=source,
                target=target,
                rate=10 ** (-features.conta_line_offset),
                probability=float(conta_probability),
                conta_line_species=features.conta_line_species,
            )
            for ((source, target), features), conta_probability in zip(
                candidates,
                conta_probabilities,
            )
        ]


# Each worker creates its feature extractor and keeps its own model instance.
# This avoids repeatedly serializing these objects when batches are submitted.
_worker: ContaminationSearcherWorker | None = None  # pylint: disable=invalid-name


def _init_worker(
    species_ab_table: pd.DataFrame,
    rf_model: RandomForestClassifier,
) -> None:
    """Initialize the contamination-search worker for a process.

    The pool calls this function once when each worker process starts. The
    abundance table and Random Forest are therefore initialized once per
    process rather than being sent with every batch of sample pairs.
    """
    global _worker  # pylint: disable=global-statement
    _worker = ContaminationSearcherWorker(
        species_ab_table,
        rf_model,
    )


def _classify_sample_pairs(
    sample_pairs: Sequence[SamplePair],
) -> tuple[int, list[ContaminationEvent]]:
    """Classify a batch of sample pairs using the current process's worker.

    Return the batch size alongside the events so the progress bar can track
    the number of sample pairs inspected rather than the number of batches.
    """
    assert _worker is not None

    return (
        len(sample_pairs),
        _worker.classify_sample_pairs(sample_pairs),
    )


class ContaminationSearcherDriver:
    """Run contamination searches using one or more worker processes."""

    DEFAULT_CHUNKSIZE: Final[int] = 50

    def __init__(
        self,
        species_ab_table: pd.DataFrame,
        all_sample_pairs: Iterator[SamplePair],
        num_sample_pairs: int,
        rf_model: RandomForestClassifier,
        probability_cutoff: float,
        rate_cutoff: float,
        nproc: int,
    ) -> None:
        self.species_ab_table = species_ab_table
        self.all_sample_pairs = all_sample_pairs
        self.num_sample_pairs = num_sample_pairs
        self.rf_model = rf_model
        self.probability_cutoff = probability_cutoff
        self.rate_cutoff = rate_cutoff
        self.nproc = nproc

    def search_contamination(
        self,
    ) -> list[ContaminationEvent]:
        """Search all sample pairs and return detected contamination events."""
        conta_events: list[ContaminationEvent] = []

        # Group sample pairs before submitting them to the pool so that each
        # worker can extract features and run the Random Forest on one batch.
        sample_pair_batches = batched(
            self.all_sample_pairs,
            self.DEFAULT_CHUNKSIZE,
        )

        with Pool(
            processes=self.nproc,
            initializer=_init_worker,
            initargs=(self.species_ab_table, self.rf_model),
        ) as pool:
            # Each pool task is already one batch, so use chunksize=1 to avoid
            # combining several batches into a larger scheduling unit.
            results = pool.imap_unordered(
                _classify_sample_pairs,
                sample_pair_batches,
                chunksize=1,
            )

            pbar = tqdm(
                total=self.num_sample_pairs,
                leave=False,
                bar_format=(
                    "{l_bar}{bar}| "
                    "{n_fmt}/{total_fmt} sample pairs inspected"
                    "{postfix}"
                ),
            )

            # The pool returns batches as soon as workers finish them. Update
            # the progress bar by the number of pairs in each completed batch.
            with pbar:
                for num_sample_pairs, batch_conta_events in results:
                    pbar.update(num_sample_pairs)

                    passing_conta_events = [
                        conta_event
                        for conta_event in batch_conta_events
                        if _passes_cutoffs(
                            conta_event,
                            self.probability_cutoff,
                            self.rate_cutoff,
                        )
                    ]

                    if passing_conta_events:
                        conta_events.extend(passing_conta_events)
                        pbar.set_postfix_str(
                            f"{len(conta_events)} conta events found",
                        )

        return conta_events
