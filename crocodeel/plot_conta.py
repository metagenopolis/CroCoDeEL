from dataclasses import dataclass, field
from typing import BinaryIO
from functools import partial
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import tqdm
from contamination_case import ContaminationCase


@dataclass
class ContaminationPlotsReport:
    mgs_profiles: pd.DataFrame
    contamination_cases: list[ContaminationCase]
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

    def _create_plot(self, contamination_case: ContaminationCase, ax) -> None:
        spearman_rho = self.mgs_profiles[contamination_case.target].corr(
            self.mgs_profiles[contamination_case.source], method="spearman"
        )

        # Do not show species absent in both samples
        # for faster rendering and reduce PDF file size
        non_zero_species = (
            self.mgs_profiles[contamination_case.target] > self.pseudo_zero
        ) | (self.mgs_profiles[contamination_case.source] > self.pseudo_zero)

        scatterplot = ax.scatter(
            x=self.mgs_profiles.loc[non_zero_species, contamination_case.target],
            y=self.mgs_profiles.loc[non_zero_species, contamination_case.source],
            s=10,
            facecolor="none",
        )

        if self.color_contamination_specific_species:
            edge_colors = [
                (
                    "orange"
                    if species in contamination_case.contamination_specific_species
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
                    self.pseudo_zero - np.log10(contamination_case.rate),
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

        ax.set_xlabel(contamination_case.target)
        ax.set_ylabel(contamination_case.source)

        ax.set_title(
            f"prob = {contamination_case.probability}, "
            f"rate = {round(contamination_case.rate * 100, 2)}%, "
            f"rho = {round(spearman_rho, 2)}"
        )

    def save_to_pdf(self, pdf_fh: BinaryIO) -> None:
        num_plots_per_page = self.nrow * self.ncol
        num_pages = int(
            np.ceil(float(len(self.contamination_cases) / num_plots_per_page))
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
                    if global_plot_id < len(self.contamination_cases):
                        self._create_plot(
                            self.contamination_cases[global_plot_id], axs[plot_id]
                        )
                    else:
                        axs[plot_id].axis("off")
                plt.tight_layout()
                fig.savefig(pdf, format="pdf")
                plt.close(fig)
