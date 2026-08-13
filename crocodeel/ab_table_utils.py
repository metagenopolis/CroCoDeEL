"""Utilities for reading and preprocessing species abundance tables."""

import logging
import sys
from typing import TextIO, Optional
import numpy as np
import pandas as pd


def _validate_and_normalize_species_names(species_ab_table: pd.DataFrame) -> None:
    """Check and normalize species names to strings."""
    species_names_type = species_ab_table.index.inferred_type

    if species_names_type not in ("integer", "string"):
        logging.error(
            "Species names in first column are of the '%s' type "
            "but should be 'string' or 'integer'",
            species_names_type,
        )
        sys.exit(1)

    species_ab_table.index = species_ab_table.index.astype(str)


def _check_numeric_abundances(species_ab_table: pd.DataFrame) -> None:
    """Check that species abundances are numeric."""
    bad_format_samples = [
        sample
        for sample in species_ab_table.columns
        if not pd.api.types.is_numeric_dtype(species_ab_table[sample].dtype)
    ]

    if bad_format_samples:
        logging.error(
            "Species abundance in the following samples is not numeric: %s",
            " ".join(bad_format_samples),
        )
        sys.exit(1)


def _check_non_negative_abundances(
    species_ab_table: pd.DataFrame,
) -> None:
    """Check that species abundances are non-negative."""
    negative_samples = species_ab_table.columns[(species_ab_table < 0).any(axis=0)]

    if not negative_samples.empty:
        logging.error(
            "The following samples contain negative species abundances: %s",
            " ".join(negative_samples),
        )
        sys.exit(1)


def _check_non_empty_samples(
    species_ab_table: pd.DataFrame,
) -> None:
    """Check that each sample contains at least one non-zero abundance."""
    empty_samples = species_ab_table.columns[species_ab_table.sum(axis=0) <= 0]

    if not empty_samples.empty:
        logging.error(
            "The following samples have no non-zero species abundances: %s",
            " ".join(empty_samples),
        )
        sys.exit(1)


def read(fh: TextIO) -> pd.DataFrame:
    """Read a species abundance table from a TSV file."""
    # Read table
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
    # Normalize to relative abundance
    species_ab_table = species_ab_table.div(species_ab_table.sum(axis=0), axis=1)

    logging.info("Species abundance table normalized")
    return species_ab_table


def filter_low_ab(
    species_ab_table: pd.DataFrame, filtering_ab_thr_factor: float
) -> pd.DataFrame:
    """Remove species with low abundances from each sample."""
    species_ab_table = species_ab_table.apply(
        lambda ab: np.where(ab <= filtering_ab_thr_factor * ab[ab > 0].min(), 0, ab)
    )
    logging.info("Low-abundance species filtered out")

    return species_ab_table


def read_filter_normalize(
    fh: TextIO, filtering_ab_thr_factor: Optional[float] = None
) -> pd.DataFrame:
    """Read, validate, optionally filter, and normalize a species abundance table."""
    species_ab_table = read(fh)
    _validate_and_normalize_species_names(species_ab_table)
    _check_numeric_abundances(species_ab_table)
    _check_non_negative_abundances(species_ab_table)
    _check_non_empty_samples(species_ab_table)

    if filtering_ab_thr_factor is not None:
        species_ab_table = filter_low_ab(species_ab_table, filtering_ab_thr_factor)
        _check_non_empty_samples(species_ab_table)

    species_ab_table = normalize(species_ab_table)

    return species_ab_table
