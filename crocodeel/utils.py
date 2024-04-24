import numpy as np
import pandas as pd
import logging
import sys
from typing import TextIO


def load_species_ab_table(fh: TextIO) -> pd.DataFrame:
    # Read table
    species_ab_table = pd.read_csv(fh, sep="\t", header=0, index_col=0)
    logging.info(
        "Abundance table quantifies %d species in %d samples", species_ab_table.shape[0], species_ab_table.shape[1]
    )

    # Check species names
    species_names_type = species_ab_table.index.inferred_type
    if species_names_type not in ("integer", "string"):
        logging.error(
            "Species names in first column are of the '%s' type but should be 'string' or 'integer'",
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
        logging.error("Species abundance in the following samples is not numeric: %s", " ".join(bad_format_samples))
        sys.exit(1)

    # Check that there are no missing values
    samples_missing_values = species_ab_table.columns[species_ab_table.isna().any()].tolist()
    if samples_missing_values:
        logging.error("The following samples have missing values: %s", " ".join(samples_missing_values))
        sys.exit(1)

    # Check that species are detected in all samples
    samples_sum = species_ab_table.sum(axis=0)
    samples_no_species = samples_sum[samples_sum == 0].index.tolist()
    if samples_no_species:
        logging.error("No species detected in the following samples: %s", " ".join(samples_no_species))
        sys.exit(1)

    # Normalize to relative abundance
    species_ab_table = species_ab_table.div(samples_sum, axis=1)

    # Perform log10 transformation
    with np.errstate(divide="ignore"):
        species_ab_table = species_ab_table.apply(np.log10)

    logging.info("Species abundance table normalized and log-transformed")

    return species_ab_table
