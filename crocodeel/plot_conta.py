from dataclasses import dataclass, field
from typing import BinaryIO
from functools import partial
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import tqdm
from conta_event import ContaminationEvent


@dataclass
class ContaminationPlotsReport:
    mgs_profiles: pd.DataFrame
    conta_events: list[ContaminationEvent]
    nrow: int = field(default=4)
    ncol: int = field(default=4)
    no_contamination_line: bool = field(default=False)
    color_contamination_specific_species: bool = field(default=False)
    pseudo_zero: float = field(init=False)

    def __post_init__(self):
        self.mgs_profiles = self.mgs_profiles.div(self.mgs_profiles.sum(axis=0), axis=1)

        # log10 transformation
        min_non_zero = self.mgs_profiles[self.mgs_profiles > 0].min().min()
        self.pseudo_zero = round(np.log10(min_non_zero))
        with np.errstate(divide="ignore"):
            self.mgs_profiles = self.mgs_profiles.apply(np.log10)
        self.mgs_profiles.replace(-np.inf, self.pseudo_zero, inplace=True)

        # Make sure that species names are strings
        self.mgs_profiles.index = self.mgs_profiles.index.astype(str)

    def _create_plot(self, conta_event: ContaminationEvent, ax) -> None:
        spearman_rho = self.mgs_profiles[conta_event.target].corr(
            self.mgs_profiles[conta_event.source], method="spearman"
        )

        # Do not show species absent in both samples
        # for faster rendering and reduce PDF file size
        non_zero_species = (
            self.mgs_profiles[conta_event.target] > self.pseudo_zero
        ) | (self.mgs_profiles[conta_event.source] > self.pseudo_zero)

        scatterplot = ax.scatter(
            x=self.mgs_profiles.loc[non_zero_species, conta_event.target],
            y=self.mgs_profiles.loc[non_zero_species, conta_event.source],
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
                for species in self.mgs_profiles.index[non_zero_species]
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
                plt.close(fig)
