"""Utilities for reading and preprocessing species abundance tables."""

import logging
from typing import Optional, TextIO

import numpy as np
import pandas as pd

from crocodeel.exceptions import InputDataError
from crocodeel.utils import format_sample_names


def _validate_and_normalize_species_names(species_ab_table: pd.DataFrame) -> None:
    """Check and normalize species names to strings."""
    species_names_type = species_ab_table.index.inferred_type

    if species_names_type not in ("integer", "string"):
        raise InputDataError(
            "Species names in first column are of the "
            f"'{species_names_type}' type but should be 'string' or 'integer'"
        )

    species_ab_table.index = species_ab_table.index.astype(str)


def _check_numeric_abundances(species_ab_table: pd.DataFrame) -> None:
    """Check that species abundances are numeric."""
    bad_format_samples = [
        sample
        for sample in species_ab_table.columns
        if not pd.api.types.is_numeric_dtype(species_ab_table[sample].dtype)
    ]

    if bad_format_samples:
        raise InputDataError(
            "Species abundance in the following samples is not numeric: "
            + format_sample_names(bad_format_samples)
        )


def _check_non_negative_abundances(
    species_ab_table: pd.DataFrame,
) -> None:
    """Check that species abundances are non-negative."""
    negative_samples = species_ab_table.columns[(species_ab_table < 0).any(axis=0)]

    if not negative_samples.empty:
        raise InputDataError(
            "The following samples contain negative species abundances: "
            + format_sample_names(negative_samples.tolist())
        )


def _check_no_missing_abundances(
    species_ab_table: pd.DataFrame,
) -> None:
    """Check that species abundances contain no missing (NA) values."""
    samples_with_missing_values = species_ab_table.columns[
        species_ab_table.isna().any(axis=0)
    ]

    if not samples_with_missing_values.empty:
        raise InputDataError(
            "The following samples contain missing species abundances: "
            + format_sample_names(samples_with_missing_values.tolist())
        )


def _check_non_empty_samples(
    species_ab_table: pd.DataFrame,
) -> None:
    """Check that each sample contains at least one non-zero abundance."""
    empty_samples = species_ab_table.columns[species_ab_table.sum(axis=0) <= 0]

    if not empty_samples.empty:
        raise InputDataError(
            "The following samples have no non-zero species abundances: "
            + format_sample_names(empty_samples.tolist())
        )


def read(fh: TextIO) -> pd.DataFrame:
    """Read a species abundance table from a TSV file."""
    logging.info("Reading %s", fh.name)
    species_ab_table = pd.read_csv(fh, sep="\t", header=0, index_col=0, comment="#")
    num_species = species_ab_table.shape[0]
    num_samples = species_ab_table.shape[1]

    logging.info(
        "Abundance table quantifies %d species in %d sample%s",
        num_species,
        num_samples,
        "s" if num_samples > 1 else "",
    )

    return species_ab_table


def normalize(species_ab_table: pd.DataFrame) -> pd.DataFrame:
    """Normalize species abundances to relative abundances."""
    species_ab_table = species_ab_table.div(species_ab_table.sum(axis=0), axis=1)

    logging.info("Species abundance table normalized")
    return species_ab_table


def log_transform(species_ab_table: pd.DataFrame) -> pd.DataFrame:
    """Apply a log10 transformation to species abundances.

    Zeros are intentionally converted to -inf.
    """
    with np.errstate(divide="ignore"):
        log_species_ab_table = pd.DataFrame(
            np.log10(species_ab_table),
            index=species_ab_table.index,
            columns=species_ab_table.columns,
        )

    logging.info("Species abundance table log-transformed")
    return log_species_ab_table


def filter_low_ab(
    species_ab_table: pd.DataFrame, filtering_ab_thr_factor: float
) -> pd.DataFrame:
    """Remove species with low abundances from each sample."""
    values = species_ab_table.to_numpy()

    min_positive = np.min(
        values,
        axis=0,
        where=values > 0,
        initial=np.inf,
    )

    filtered = np.where(
        values <= filtering_ab_thr_factor * min_positive,
        0.0,
        values,
    )

    logging.info("Low-abundance species filtered out")

    return pd.DataFrame(
        filtered,
        index=species_ab_table.index,
        columns=species_ab_table.columns,
    )


def read_filter_normalize(
    fh: TextIO, filtering_ab_thr_factor: Optional[float] = None
) -> pd.DataFrame:
    """Read, validate, optionally filter, and normalize a species abundance table."""
    species_ab_table = read(fh)
    _validate_and_normalize_species_names(species_ab_table)
    _check_numeric_abundances(species_ab_table)

    # Use a floating-point representation for all subsequent preprocessing.
    species_ab_table = species_ab_table.astype(float)

    _check_no_missing_abundances(species_ab_table)
    _check_non_negative_abundances(species_ab_table)
    _check_non_empty_samples(species_ab_table)

    if filtering_ab_thr_factor is not None:
        species_ab_table = filter_low_ab(species_ab_table, filtering_ab_thr_factor)
        _check_non_empty_samples(species_ab_table)

    species_ab_table = normalize(species_ab_table)
    species_ab_table = log_transform(species_ab_table)

    return species_ab_table


def compare_species_names(species_ab_table: pd.DataFrame, species_ab_table_2: pd.DataFrame) -> None:
    """Compare species names between two abundance tables.
    
    Logs warnings when the tables do not contain the same set of species.
    """
    species_names = set(species_ab_table.index)
    species_names_2 = set(species_ab_table_2.index)

    if species_names != species_names_2:
        logging.warning(
            "Abundance tables have only %d species names in common",
            len(species_names.intersection(species_names_2)),
        )
        logging.warning(
            "Make sure the abundance tables were generated with the same tool and database"
        )
        logging.warning(
            "Missing abundance values will be filled with zeros for non-shared species"
        )
