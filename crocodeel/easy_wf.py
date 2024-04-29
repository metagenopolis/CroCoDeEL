from typing import Any
from species_ab_table import SpeciesAbTableUtils
from search_conta import run_search_conta
from plot_conta import run_plot_conta, Defaults as plot_conta_defaults

def run_easy_wf(args: dict[str,Any]):
    species_ab_table = SpeciesAbTableUtils.load(args["species_ab_table_fh"])
    args["species_ab_table_fh"].close()
    args["species_ab_table"] = SpeciesAbTableUtils.normalize(species_ab_table)

    run_search_conta(args)

    args["conta_events_fh"].seek(0)
    args["nrow"] = plot_conta_defaults.NROW
    args["ncol"] = plot_conta_defaults.NCOL
    args["no_conta_line"] = False
    args["color_conta_species"] = False
    run_plot_conta(args)



