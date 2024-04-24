import numpy as np
import pandas as pd
from typing import TextIO


def load_species_ab_table(fh: TextIO) -> pd.DataFrame:
    # Read table
    species_ab_table = pd.read_csv(fh, sep="\t", header=0, index_col=0)

    # Normalize to relative abundance
    species_ab_table = species_ab_table.div(species_ab_table.sum(axis=0), axis=1)

    #Perform log10 transformation
    with np.errstate(divide="ignore"):
        species_ab_table = species_ab_table.apply(np.log10)

    # Make sure that species names are strings
    species_ab_table.index = species_ab_table.index.astype(str)

    return species_ab_table
