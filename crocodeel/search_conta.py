"""Search for cross-sample contamination in metagenomic abundance data."""

import importlib.resources
import logging
import warnings
from itertools import product
from multiprocessing import Pool
from pathlib import Path
from time import perf_counter
from typing import BinaryIO, Final, Iterator

import joblib
import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.exceptions import InconsistentVersionWarning  # type: ignore[attr-defined]
from tqdm import tqdm

from crocodeel.conta_event import ContaminationEvent
from crocodeel.conta_features import (ContaminationFeatureExtractor,
                                      ContaminationFeatures)
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
    """Classify individual sample pairs for contamination."""

    def __init__(
        self,
        species_ab_table: pd.DataFrame,
        rf_classifier: RandomForestClassifier,
    ) -> None:
        self.feature_extractor = ContaminationFeatureExtractor(species_ab_table)
        self.rf_classifier = rf_classifier

    def classify_sample_pair(
        self,
        sample_pair: SamplePair,
    ) -> ContaminationEvent | None:
        """Classify a sample pair and return a contamination event if detected."""
        source, target = sample_pair

        if source == target:
            return None

        features: ContaminationFeatures | None = self.feature_extractor.extract(
            source,
            target,
        )

        if features is None:
            return None

        probabilities = self.rf_classifier.predict_proba(features.values.reshape(1, -1))
        conta_probability = float(probabilities[0, 1])

        return ContaminationEvent(
            source=source,
            target=target,
            rate=10 ** (-features.conta_line_offset),
            probability=conta_probability,
            conta_line_species=features.conta_line_species,
        )


# Initialize the worker once per process so the abundance table and model
# are sent to each worker only once, rather than being re-serialized per chunk.
_worker: ContaminationSearcherWorker | None = None  # pylint: disable=invalid-name


def _init_worker(
    species_ab_table: pd.DataFrame,
    rf_model: RandomForestClassifier,
) -> None:
    """Build the per-process contamination-search worker.

    Used as a ``multiprocessing.Pool`` initializer, so it runs exactly once
    per worker process.
    """
    global _worker  # pylint: disable=global-statement
    _worker = ContaminationSearcherWorker(species_ab_table, rf_model)


def _classify_sample_pair(sample_pair: SamplePair) -> ContaminationEvent | None:
    """Classify one sample pair using the current process's worker."""
    assert _worker is not None
    return _worker.classify_sample_pair(sample_pair)


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

        with Pool(
            processes=self.nproc,
            initializer=_init_worker,
            initargs=(self.species_ab_table, self.rf_model),
        ) as pool:
            results = pool.imap_unordered(
                _classify_sample_pair,
                self.all_sample_pairs,
                chunksize=self.DEFAULT_CHUNKSIZE,
            )

            pbar = tqdm(
                results,
                total=self.num_sample_pairs,
                leave=False,
                bar_format=(
                    "{l_bar}{bar}| "
                    "{n_fmt}/{total_fmt} sample pairs inspected"
                    "{postfix}"
                ),
            )

            for conta_event in pbar:
                if conta_event is not None and _passes_cutoffs(
                    conta_event,
                    self.probability_cutoff,
                    self.rate_cutoff,
                ):
                    conta_events.append(conta_event)
                    pbar.set_postfix_str(f"{len(conta_events)} conta events found")

        return conta_events
