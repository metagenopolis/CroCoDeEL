"""Command-line interface for CroCoDeEL."""

import argparse
import logging
import math
import multiprocessing
import os
from importlib.metadata import version
from pathlib import Path

import pandas as pd

from crocodeel import ab_table_utils
from crocodeel.conta_event import ContaminationEvent, ContaminationEventIO
from crocodeel.execution_description import ExecutionDescription
from crocodeel.plot_conta import (
    Defaults as plot_conta_defaults,
    run_plot_conta,
)
from crocodeel.search_conta import (
    Defaults as search_conta_defaults,
    run_search_conta,
)
from crocodeel.self_test import SelfTest
from crocodeel.train_model import (
    Defaults as train_model_defaults,
    run_train_model,
)
from crocodeel.exceptions import InputDataError, SelfTestError


def set_logging() -> None:
    """Configure logging for the command-line application."""
    logging.basicConfig(
        format="%(asctime)s :: %(levelname)s :: %(message)s",
        level=logging.INFO,
    )


def readable_file(fp_str: str) -> Path:
    """Validate that a path points to a readable regular file."""
    fp = Path(fp_str).resolve()

    if not fp.exists():
        raise argparse.ArgumentTypeError(f"{fp} does not exist")

    if not fp.is_file():
        raise argparse.ArgumentTypeError(f"{fp} is not a regular file.")

    if not os.access(fp, os.R_OK):
        raise argparse.ArgumentTypeError(f"{fp} is not readable.")

    return fp


def writable_file(fp_str: str) -> Path:
    """Validate that a path can be used as an output file."""
    fp = Path(fp_str).resolve()

    if fp.exists():
        if fp.is_dir():
            raise argparse.ArgumentTypeError(
                f"{fp} is a directory, not a file."
            )

        if not os.access(fp, os.W_OK):
            raise argparse.ArgumentTypeError(f"{fp} is not writable.")

        return fp

    parent_dir = fp.parent

    if not parent_dir.exists():
        raise argparse.ArgumentTypeError(
            f"directory {parent_dir} does not exist."
        )

    if not parent_dir.is_dir():
        raise argparse.ArgumentTypeError(
            f"{parent_dir} is not a directory."
        )

    if not os.access(parent_dir, os.W_OK):
        raise argparse.ArgumentTypeError(
            f"directory {parent_dir} is not writable."
        )

    return fp


def nproc(value: str) -> int:
    """Validate the number of requested parallel processes."""
    max_nproc = multiprocessing.cpu_count()

    try:
        ivalue = int(value)
    except ValueError as value_err:
        raise argparse.ArgumentTypeError(
            f"{value} is not an integer"
        ) from value_err

    if ivalue <= 0:
        raise argparse.ArgumentTypeError("minimum value is 1")

    if ivalue > max_nproc:
        raise argparse.ArgumentTypeError(
            f"maximum value is {max_nproc}"
        )

    return ivalue


def positive_int(value: str) -> int:
    """Validate a strictly positive integer."""
    try:
        ivalue = int(value)
    except ValueError as value_err:
        raise argparse.ArgumentTypeError(
            f"{value} is not an integer"
        ) from value_err

    if ivalue <= 0:
        raise argparse.ArgumentTypeError(
            "value must be greater than 0"
        )

    return ivalue


def bounded_float_01(value: str) -> float:
    """Parse a string as a finite float between 0 and 1."""
    try:
        fvalue = float(value)
    except ValueError as value_err:
        raise argparse.ArgumentTypeError(
            f"{value} is not a valid float"
        ) from value_err

    if not math.isfinite(fvalue) or not 0.0 <= fvalue <= 1.0:
        raise argparse.ArgumentTypeError(
            "value must be a finite float between 0 and 1"
        )

    return fvalue


def positive_float(value: str) -> float:
    """Validate a strictly positive floating-point value."""
    try:
        fvalue = float(value)
    except ValueError as value_err:
        raise argparse.ArgumentTypeError(
            f"{value} is not a number"
        ) from value_err

    if not math.isfinite(fvalue) or fvalue <= 0:
        raise argparse.ArgumentTypeError(
            "value must be a finite number greater than 0"
        )

    return fvalue


def add_abundance_table_arguments(
    parser: argparse.ArgumentParser,
) -> None:
    """Add arguments used to specify species abundance tables."""
    parser.add_argument(
        "-s",
        dest="species_ab_table_fp",
        type=readable_file,
        required=True,
        metavar="SPECIES_ABUNDANCE_TABLE",
        help="Input TSV file corresponding to the species abundance table",
    )

    parser.add_argument(
        "-s2",
        dest="species_ab_table_2_fp",
        type=readable_file,
        required=False,
        metavar="SPECIES_ABUNDANCE_TABLE_2",
        help=(
            "Optional input TSV file corresponding to another species "
            "abundance table. If provided, samples from this table will "
            "be considered as contamination targets while those from the "
            "first table as contamination sources."
        ),
    )

    parser.add_argument(
        "--filter-low-ab",
        dest="filtering_ab_thr_factor",
        type=positive_float,
        default=None,
        metavar="AB_THRESHOLD_FACTOR",
        help=(
            "Filter out low-abundance species that may be inaccurately "
            "quantified. In each sample, set the abundance of species to "
            "zero if they are up to %(metavar)s times more abundant than "
            "the least abundant species. Recommended value for "
            "MetaPhlAn4: 20 (default: None)"
        ),
    )


def add_train_model_abundance_table_arguments(
    parser: argparse.ArgumentParser,
) -> None:
    """Add abundance table arguments used for model training."""
    parser.add_argument(
        "-s",
        dest="species_ab_table_fp",
        type=readable_file,
        required=True,
        metavar="SPECIES_ABUNDANCE_TABLE",
        help="Input TSV file corresponding to the species abundance table",
    )

    parser.add_argument(
        "--filter-low-ab",
        dest="filtering_ab_thr_factor",
        type=positive_float,
        default=None,
        metavar="AB_THRESHOLD_FACTOR",
        help=(
            "Filter out low-abundance species that may be inaccurately "
            "quantified. In each sample, set the abundance of species to "
            "zero if they are up to %(metavar)s times more abundant than "
            "the least abundant species. Recommended value for "
            "MetaPhlAn4: 20 (default: None)"
        ),
    )


def add_search_arguments(
    parser: argparse.ArgumentParser,
    include_conta_events: bool = True,
) -> None:
    """Add arguments used to search for contamination."""
    parser.add_argument(
        "-m",
        dest="rf_model_fp",
        type=readable_file,
        default=search_conta_defaults.MODEL_FILE,
        metavar="RF_MODEL_FILE",
        help=(
            "Joblib file containing the pre-trained Random Forest model "
            "used to detect contamination events "
            "(default: %(default)s)"
        ),
    )

    parser.add_argument(
        "--probability-cutoff",
        dest="probability_cutoff",
        type=bounded_float_01,
        default=search_conta_defaults.PROBABILITY_CUTOFF,
        metavar="PROBABILITY_CUTOFF",
        help=(
            "Only report contamination events with a probability "
            "greater than %(metavar)s "
            "(default: %(default).2f)"
        ),
    )

    parser.add_argument(
        "--rate-cutoff",
        dest="rate_cutoff",
        type=bounded_float_01,
        default=search_conta_defaults.RATE_CUTOFF,
        metavar="RATE_CUTOFF",
        help=(
            "Only report events with a contamination rate greater "
            "than %(metavar)s "
            "(default: %(default).0f)"
        ),
    )

    if include_conta_events:
        parser.add_argument(
            "-c",
            dest="conta_events_fp",
            type=writable_file,
            required=True,
            metavar="CONTAMINATION_EVENTS_FILE",
            help="Output TSV file listing all contamination events",
        )

    parser.add_argument(
        "--nproc",
        dest="nproc",
        type=nproc,
        default=multiprocessing.cpu_count(),
        help=(
            "Number of parallel processes to search contaminations "
            "(default: %(default)d)"
        ),
    )


def add_plot_arguments(
    parser: argparse.ArgumentParser,
    include_conta_events: bool = True,
) -> None:
    """Add arguments used to generate the contamination report."""
    if include_conta_events:
        parser.add_argument(
            "-c",
            dest="conta_events_fp",
            type=readable_file,
            required=True,
            metavar="CONTAMINATION_EVENTS_FILE",
            help="Input TSV file listing all contamination events.",
        )

    parser.add_argument(
        "-r",
        dest="pdf_report_fp",
        type=writable_file,
        required=True,
        metavar="PDF_REPORT_FILE",
        help="Output PDF file with scatterplots for all contamination events",
    )

    parser.add_argument(
        "--nrow",
        dest="nrow",
        type=positive_int,
        choices=range(
            plot_conta_defaults.MIN_NROW,
            plot_conta_defaults.MAX_NROW + 1,
        ),
        default=plot_conta_defaults.NROW,
        metavar="NROW",
        help=(
            "Number of scatterplots to draw vertically on each page "
            "(default: %(default)d)"
        ),
    )

    parser.add_argument(
        "--ncol",
        dest="ncol",
        type=positive_int,
        choices=range(
            plot_conta_defaults.MIN_NCOL,
            plot_conta_defaults.MAX_NCOL + 1,
        ),
        default=plot_conta_defaults.NCOL,
        metavar="NCOL",
        help=(
            "Number of scatterplots to draw horizontally on each page "
            "(default: %(default)d)"
        ),
    )

    parser.add_argument(
        "--no-conta-line",
        dest="no_conta_line",
        action="store_true",
        help="Do not show contamination line in scatterplots",
    )

    parser.add_argument(
        "--color-conta-species",
        dest="color_conta_species",
        action="store_true",
        help="Use a different color for species introduced by contamination",
    )


def add_train_model_arguments(
    parser: argparse.ArgumentParser,
) -> None:
    """Add arguments used to train the Random Forest model."""
    parser.add_argument(
        "-m",
        dest="model_fp",
        type=writable_file,
        required=True,
        metavar="MODEL_FILE",
        help="Output file storing the trained Random Forest model",
    )

    parser.add_argument(
        "-r",
        dest="json_report_fp",
        type=writable_file,
        required=True,
        metavar="JSON_REPORT_FILE",
        help=(
            "Output JSON file storing classification performance metrics "
            "for train and test splits"
        ),
    )

    parser.add_argument(
        "--test-size",
        dest="test_size",
        type=bounded_float_01,
        default=train_model_defaults.TEST_SIZE,
        metavar="TEST_SIZE",
        help=(
            "Proportion of dataset to include in test split "
            "(default: %(default).2f)"
        ),
    )

    parser.add_argument(
        "--ntrees",
        dest="ntrees",
        type=positive_int,
        default=train_model_defaults.NTREES,
        metavar="NTREES",
        help=(
            "Number of trees in the RandomForest model "
            "(default: %(default)d)"
        ),
    )

    parser.add_argument(
        "--rng-seed",
        dest="rng_seed",
        type=int,
        default=train_model_defaults.RNG_SEED,
        metavar="RNG_SEED",
        help=(
            "Seed of the random number generator for reproducibility "
            "(default: %(default)d)"
        ),
    )

    parser.add_argument(
        "--nproc",
        dest="nproc",
        type=nproc,
        default=multiprocessing.cpu_count(),
        help=(
            "Number of parallel processes to train the model "
            "(default: %(default)d)"
        ),
    )


def get_arguments() -> argparse.Namespace:
    """Parse and validate command-line arguments."""
    prog_name = "CroCoDeEL"
    prog_version = version(prog_name.lower())

    parser = argparse.ArgumentParser(
        description=(
            f"{prog_name} is a tool that detects cross-sample contamination "
            "in shotgun metagenomic data"
        ),
    )

    parser.add_argument(
        "-v",
        "--version",
        action="version",
        version=f"{prog_name} version {prog_version}",
    )

    subparsers = parser.add_subparsers(
        title=f"{prog_name} subcommands",
        dest="command",
        required=True,
    )

    easy_wf_parser = subparsers.add_parser(
        "easy_wf",
        help=(
            "Run search_conta and plot_conta in one command "
            "to detect cross-sample contamination and generate a PDF report"
        ),
    )

    search_conta_parser = subparsers.add_parser(
        "search_conta",
        help="Search for cross-sample contamination",
    )

    plot_conta_parser = subparsers.add_parser(
        "plot_conta",
        help=(
            "Create a PDF report with scatterplots representing "
            "species abundance profiles for each contamination event"
        ),
    )

    self_test_parser = subparsers.add_parser(
        "self_test",
        help=f"Run an end-to-end test of the {prog_name} workflow on a toy dataset",
    )

    train_model_parser = subparsers.add_parser(
        "train_model",
        help=(
            f"Train the Random Forest model used by {prog_name} "
            "to classify sample pairs (for advanced users)"
        ),
    )

    # search_conta
    add_abundance_table_arguments(search_conta_parser)
    add_search_arguments(search_conta_parser)

    # plot_conta
    add_abundance_table_arguments(plot_conta_parser)
    add_plot_arguments(plot_conta_parser)

    # easy_wf
    add_abundance_table_arguments(easy_wf_parser)
    add_search_arguments(easy_wf_parser)
    add_plot_arguments(
        easy_wf_parser,
        include_conta_events=False,
    )

    # train_model
    add_train_model_abundance_table_arguments(train_model_parser)
    add_train_model_arguments(train_model_parser)

    # self_test
    self_test_parser.add_argument(
        "--keep-results",
        dest="keep_results",
        action="store_true",
        help="Keep all temporary result files.",
    )

    return parser.parse_args()


def load_abundance_tables(
    args: argparse.Namespace,
) -> tuple[pd.DataFrame, pd.DataFrame | None]:
    """Load and validate the species abundance tables."""
    with args.species_ab_table_fp.open(
        "r",
        encoding="utf8",
    ) as species_ab_table_fh:
        species_ab_table = ab_table_utils.read_filter_normalize(
            species_ab_table_fh,
            args.filtering_ab_thr_factor,
        )

    species_ab_table_2 = None

    if (
        args.species_ab_table_2_fp is not None
        and args.species_ab_table_fp != args.species_ab_table_2_fp
    ):
        with args.species_ab_table_2_fp.open(
            "r",
            encoding="utf8",
        ) as species_ab_table_2_fh:
            species_ab_table_2 = ab_table_utils.read_filter_normalize(
                species_ab_table_2_fh,
                args.filtering_ab_thr_factor,
            )

        ab_table_utils.compare_species_names(
            species_ab_table,
            species_ab_table_2,
        )

    return species_ab_table, species_ab_table_2


def run_search(
    args: argparse.Namespace,
    species_ab_table: pd.DataFrame,
    species_ab_table_2: pd.DataFrame | None,
) -> list[ContaminationEvent]:
    """Search for contamination events and save them to a TSV file."""
    exec_desc = ExecutionDescription(
        args.species_ab_table_fp,
        args.species_ab_table_2_fp,
        args.rf_model_fp,
        args.filtering_ab_thr_factor,
        args.probability_cutoff,
        args.rate_cutoff,
    )

    with (
        args.rf_model_fp.open("rb") as rf_model_fh,
        args.conta_events_fp.open(
            "w",
            encoding="utf8",
        ) as conta_events_fh,
    ):
        print(exec_desc, file=conta_events_fh)

        conta_events = run_search_conta(
            species_ab_table,
            species_ab_table_2,
            rf_model_fh,
            args.probability_cutoff,
            args.rate_cutoff,
            args.nproc,
        )

        ContaminationEventIO.write_tsv(
            conta_events,
            conta_events_fh,
        )

    return conta_events


def load_contamination_events(
    args: argparse.Namespace,
) -> list[ContaminationEvent]:
    """Load contamination events from a TSV file."""
    with args.conta_events_fp.open(
        "r",
        encoding="utf8",
    ) as conta_events_fh:
        return ContaminationEventIO.read_tsv(
            conta_events_fh,
        )


def generate_pdf_report(
    args: argparse.Namespace,
    species_ab_table: pd.DataFrame,
    species_ab_table_2: pd.DataFrame | None,
    conta_events: list[ContaminationEvent],
) -> None:
    """Generate the PDF contamination report."""
    with args.pdf_report_fp.open("wb") as pdf_report_fh:
        run_plot_conta(
            species_ab_table,
            species_ab_table_2,
            conta_events,
            pdf_report_fh,
            args.nrow,
            args.ncol,
            args.no_conta_line,
            args.color_conta_species,
        )


def log_contamination_warnings(
    conta_events: list[ContaminationEvent],
) -> None:
    """Log warnings associated with detected contamination events."""
    if not conta_events:
        return

    logging.warning(
        "\033[31m\033[1mCroCoDeEL is a decision-support tool and should not be "
        "considered a definitive contamination classification system.\033[0m"
    )

    logging.warning(
        "Reported contamination events may include false positives, "
        "particularly for samples with similar species abundance profiles "
        "(e.g. longitudinal samples)."
    )

    logging.warning(
        "We strongly recommend manually reviewing the scatterplots "
        "associated with each predicted contamination event to identify "
        "potential false positives."
    )

    logging.warning(
        "For easier exploration and validation of contamination events, "
        "consider using the CroCoDeEL Interpretation Interface: "
        "https://metagenopolis.github.io/CroCoDeEL_interpreter/"
    )


def run_search_conta_command(
    args: argparse.Namespace,
) -> None:
    """Run the search_conta command."""
    species_ab_table, species_ab_table_2 = load_abundance_tables(args)

    conta_events = run_search(
        args,
        species_ab_table,
        species_ab_table_2,
    )

    log_contamination_warnings(conta_events)


def run_plot_conta_command(
    args: argparse.Namespace,
) -> None:
    """Run the plot_conta command."""
    species_ab_table, species_ab_table_2 = load_abundance_tables(args)

    conta_events = load_contamination_events(args)

    generate_pdf_report(
        args,
        species_ab_table,
        species_ab_table_2,
        conta_events,
    )


def run_easy_workflow(
    args: argparse.Namespace,
) -> None:
    """Run contamination detection and generate the PDF report."""
    species_ab_table, species_ab_table_2 = load_abundance_tables(args)

    conta_events = run_search(
        args,
        species_ab_table,
        species_ab_table_2,
    )

    generate_pdf_report(
        args,
        species_ab_table,
        species_ab_table_2,
        conta_events,
    )

    log_contamination_warnings(conta_events)


def run_train_model_command(
    args: argparse.Namespace,
) -> None:
    """Run the train_model command."""
    with args.species_ab_table_fp.open(
        "r",
        encoding="utf8",
    ) as species_ab_table_fh:
        species_ab_table = ab_table_utils.read_filter_normalize(
            species_ab_table_fh,
            args.filtering_ab_thr_factor,
        )

    with (
        args.model_fp.open("wb") as model_fh,
        args.json_report_fp.open(
            "w",
            encoding="utf8",
        ) as json_report_fh,
    ):
        run_train_model(
            species_ab_table,
            model_fh,
            json_report_fh,
            args.test_size,
            args.ntrees,
            args.rng_seed,
            args.nproc,
        )


def main() -> int:
    """Run the CroCoDeEL command-line application."""
    set_logging()
    args = get_arguments()

    try:
        if args.command == "self_test":
            SelfTest(args.keep_results).run()
        elif args.command == "search_conta":
            run_search_conta_command(args)
        elif args.command == "plot_conta":
            run_plot_conta_command(args)
        elif args.command == "easy_wf":
            run_easy_workflow(args)
        elif args.command == "train_model":
            run_train_model_command(args)

        return 0

    except (InputDataError, SelfTestError) as error:
        logging.error("%s", error)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
