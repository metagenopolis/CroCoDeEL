import logging
import sys
from typing import TextIO, Optional
from pathlib import Path
import numpy as np
import pandas as pd


def read(fh: TextIO) -> pd.DataFrame:
    # Read table
    logging.info("Reading %s", Path(fh.name).resolve())
    species_ab_table = pd.read_csv(fh, sep="\t", header=0, index_col=0, comment="#")
    logging.info(
        "Abundance table quantifies %d species in %d samples",
        species_ab_table.shape[0],
        species_ab_table.shape[1],
    )

    # Check species names
    species_names_type = species_ab_table.index.inferred_type
    if species_names_type not in ("integer", "string"):
        logging.error(
            "Species names in first column are of the '%s' type "
            "but should be 'string' or 'integer'",
            species_names_type,
        )
        sys.exit(1)

    # Convert species names to strings when they are integers
    species_ab_table.index = species_ab_table.index.astype(str)

    # Check species abundance type
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

    return species_ab_table


def normalize(species_ab_table: pd.DataFrame) -> pd.DataFrame :
    # Normalize to relative abundance
    species_ab_table = species_ab_table.div(species_ab_table.sum(axis=0), axis=1)

    # Perform log10 transformation
    with np.errstate(divide="ignore"):
        species_ab_table = species_ab_table.apply(np.log10)

    logging.info("Species abundance table normalized and log-transformed")
    return species_ab_table


def filter_low_ab(
    species_ab_table: pd.DataFrame, filtering_ab_thr_factor: float
) -> pd.DataFrame:
    species_ab_table = species_ab_table.apply(
        lambda ab: np.where(ab <= filtering_ab_thr_factor * ab[ab > 0].min(), 0, ab)
    )
    logging.info("Low-abundance species filtered out")

    return species_ab_table


def read_filter_normalize(
    fh: TextIO, filtering_ab_thr_factor: Optional[float]
) -> pd.DataFrame:
    species_ab_table = read(fh)
    if filtering_ab_thr_factor is not None:
        species_ab_table = filter_low_ab(species_ab_table, filtering_ab_thr_factor)
    species_ab_table = normalize(species_ab_table)

    return species_ab_table


def compare_species_names(species_ab_table: pd.DataFrame, species_ab_table_2: pd.DataFrame) -> None:
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
