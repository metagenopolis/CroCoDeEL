from typing import BinaryIO, Any, Final
from functools import partial
import logging
from pathlib import Path
from time import perf_counter
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import tqdm
from crocodeel.conta_event import ContaminationEvent, ContaminationEventIO
from crocodeel import ab_table_utils


def run_plot_conta(args: dict[str, Any]):
    if "species_ab_table" in args:
        species_ab_table = args["species_ab_table"]
    else:
        species_ab_table = ab_table_utils.read(args["species_ab_table_fh"])
        args["species_ab_table_fh"].close()
        species_ab_table = ab_table_utils.normalize(species_ab_table)

    conta_events = list(ContaminationEventIO.read_tsv(args["conta_events_fh"]))
    if args["conta_events_fh"].mode == 'r':
        logging.info(
            "%d contamination events loaded from %s",
            len(conta_events),
            Path(args["conta_events_fh"].name).resolve(),
        )
    args["conta_events_fh"].close()

    start = perf_counter()
    logging.info("Generation of the PDF report started")
    ContaminationPlotsReport(
        species_ab_table=species_ab_table,
        conta_events=conta_events,
        nrow=args["nrow"],
        ncol=args["ncol"],
        no_contamination_line=args["no_conta_line"],
        color_contamination_specific_species=args["color_conta_species"],
    ).save_to_pdf(args["pdf_report_fh"])
    logging.info("PDF report generated in %.1f seconds", np.round(perf_counter() - start, 1))
    logging.info("PDF report saved in %s", Path(args["pdf_report_fh"].name).resolve())
    args["pdf_report_fh"].close()


class Defaults:
    MIN_NROW: Final[int] = 1
    NROW: Final[int] = 4
    MAX_NROW: Final[int] = 11
    MIN_NCOL: Final[int] = 1
    NCOL: Final[int] = 4
    MAX_NCOL: Final[int] = 11


class ContaminationPlotsReport:

    def __init__(
        self,
        species_ab_table: pd.DataFrame,
        conta_events: list[ContaminationEvent],
        nrow: int,
        ncol: int,
        no_contamination_line: bool,
        color_contamination_specific_species: bool,
    ):
        self.species_ab_table = species_ab_table
        self.conta_events = conta_events
        self.nrow = nrow
        self.ncol = ncol
        self.no_contamination_line = no_contamination_line
        self.color_contamination_specific_species = color_contamination_specific_species

        # Add pseudo_zero
        min_non_zero = self.species_ab_table[self.species_ab_table > -np.inf].min().min()
        self.pseudo_zero = int(np.floor(min_non_zero))
        self.species_ab_table.replace(-np.inf, self.pseudo_zero, inplace=True)

    def _create_plot(self, conta_event: ContaminationEvent, ax) -> None:
        spearman_rho = self.species_ab_table[conta_event.target].corr(
            self.species_ab_table[conta_event.source], method="spearman"
        )

        # Do not show species absent in both samples
        # for faster rendering and reduce PDF file size
        non_zero_species = (
            self.species_ab_table[conta_event.target] > self.pseudo_zero
        ) | (self.species_ab_table[conta_event.source] > self.pseudo_zero)

        scatterplot = ax.scatter(
            x=self.species_ab_table.loc[non_zero_species, conta_event.target],
            y=self.species_ab_table.loc[non_zero_species, conta_event.source],
            s=10,
            facecolor="none",
        )

        if self.color_contamination_specific_species:
            edge_colors = [
                (
                    "orange"
                    if species in conta_event.contamination_specific_species
                    else "black"
                )
                for species in self.species_ab_table.index[non_zero_species]
            ]
            scatterplot.set_edgecolor(edge_colors)
        else:
            scatterplot.set_edgecolor("black")

        # Add identity line line
        ax.axline(
            [self.pseudo_zero, self.pseudo_zero],
            slope=1,
            color="grey",
            linestyle="-",
            linewidth=0.5,
        )

        # Add contamination line
        if not self.no_contamination_line:
            ax.axline(
                (
                    self.pseudo_zero,
                    self.pseudo_zero - np.log10(conta_event.rate),
                ),
                slope=1,
                color="red",
                linestyle="-",
                linewidth=0.1,
                alpha=0.5,
            )

        # Set labels and title
        ticks = np.arange(self.pseudo_zero, 1)
        ax.set_xticks(ticks)
        ax.set_yticks(ticks)
        ax.tick_params(direction="in")

        ticks_labels = [f"$10^{{{t}}}$" for t in ticks]
        ticks_labels[0] = "0"
        ticks_labels[-1] = "1"
        ax.set_xticklabels(ticks_labels)
        ax.set_yticklabels(ticks_labels)

        ax.set_xlabel(conta_event.target)
        ax.set_ylabel(conta_event.source)

        ax.set_title(
            f"prob = {conta_event.probability}, "
            f"rate = {round(conta_event.rate * 100, 2)}%, "
            f"rho = {round(spearman_rho, 2)}"
        )

    def save_to_pdf(self, pdf_fh: BinaryIO) -> None:
        num_plots_per_page = self.nrow * self.ncol
        num_pages = int(
            np.ceil(float(len(self.conta_events) / num_plots_per_page))
        )
        with PdfPages(pdf_fh) as pdf:
            pbar = partial(
                tqdm.tqdm,
                total=num_pages,
                bar_format="{l_bar}{bar}| {n_fmt}/{total_fmt} pages generated",
            )
            for page in pbar(range(num_pages)):
                fig = plt.figure()
                fig, axs = plt.subplots(
                    self.nrow, self.ncol, figsize=(4 * self.ncol, 4 * self.nrow)
                )
                axs = axs.flatten()
                for plot_id in range(num_plots_per_page):
                    global_plot_id = (page * num_plots_per_page) + plot_id
                    if global_plot_id < len(self.conta_events):
                        self._create_plot(
                            self.conta_events[global_plot_id], axs[plot_id]
                        )
                    else:
                        axs[plot_id].axis("off")
                plt.tight_layout()
                fig.savefig(pdf, format="pdf")
                plt.close('all')
