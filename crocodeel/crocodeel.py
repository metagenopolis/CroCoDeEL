#!/bin/env python3
# -*- coding: utf-8 -*-

import sys
import argparse
import logging
import multiprocessing
from pathlib import Path
import os
from importlib.metadata import version
from crocodeel.execution_description import ExecutionDescription
from crocodeel import ab_table_utils
from crocodeel.conta_event import ContaminationEventIO
from crocodeel.search_conta import run_search_conta, Defaults as search_conta_defaults
from crocodeel.plot_conta import run_plot_conta, Defaults as plot_conta_defaults
from crocodeel.train_model import run_train_model, Defaults as train_model_defaults
from crocodeel.test_install import TestInstall


def set_logging() -> None:
    logging.basicConfig(
        format="%(asctime)s :: %(levelname)s :: %(message)s", level=logging.INFO
    )

def readable_file(fp_str: str) -> Path:
    fp = Path(fp_str).resolve()

    if not fp.exists():
        raise argparse.ArgumentTypeError(f"{fp} does not exist")
    if not fp.is_file():
        raise argparse.ArgumentTypeError(f"{fp} is not a regular file.")
    if not os.access(fp, os.R_OK):
        raise argparse.ArgumentTypeError(f"{fp} is not readable.")

    return fp


def writable_file(fp_str: str) -> Path:
    fp = Path(fp_str).resolve()

    if fp.exists():
        if fp.is_dir():
            raise argparse.ArgumentTypeError(f"{fp} is a directory, not a file.")
        if not os.access(fp, os.W_OK):
            raise argparse.ArgumentTypeError(f"{fp} is not writable.")
        return fp

    parent_dir = fp.parent or Path(".")
    if not parent_dir.exists():
        raise argparse.ArgumentTypeError(f"directory {parent_dir} does not exist.")
    if not parent_dir.is_dir():
        raise argparse.ArgumentTypeError(f"{parent_dir} is not a directory.")
    if not os.access(parent_dir, os.W_OK):
        raise argparse.ArgumentTypeError(f"directory {parent_dir} is not writable.")

    return fp

def nproc(value: str) -> int:
    max_nproc = multiprocessing.cpu_count()

    try:
        ivalue = int(value)
    except ValueError as value_err:
        raise argparse.ArgumentTypeError(f"{value} is not an integer") from value_err

    if ivalue <= 0:
        raise argparse.ArgumentTypeError("minimum value is 1")
    if ivalue > max_nproc:
        raise argparse.ArgumentTypeError(f"maximum value is {max_nproc}")

    return ivalue


def bounded_float_01(value: str) -> float:
    try:
        fvalue = float(value)
    except ValueError as value_err:
        raise argparse.ArgumentTypeError(f"{value} is not a valid float") from value_err

    if not 0.0 <= fvalue <= 1.0:
        raise argparse.ArgumentTypeError("value must be a float between 0 and 1")

    return fvalue


def get_arguments() -> argparse.Namespace:
    prog_name = "CroCoDeEL"
    prog_version = version(prog_name.lower())
    parser = argparse.ArgumentParser(
        description=f"{prog_name} is a tool that detects cross-sample contamination "
        "in shotgun metagenomic data",
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
        help="Run search_conta and plot_conta in one command "
        "to detect cross-sample contamination and generate a PDF report",
    )
    search_conta_parser = subparsers.add_parser(
        "search_conta", help="Search for cross-sample contamination"
    )
    plot_conta_parser = subparsers.add_parser(
        "plot_conta",
        help="Create a PDF report with scatterplots representing "
        "species abundance profiles for each contamination event",
    )
    test_install_parser = subparsers.add_parser(
        "test_install",
        help=f"Test if {prog_name} is correctly installed "
        "and generates the expected results",
    )
    train_model_parser = subparsers.add_parser(
        "train_model",
        help=f"Train the Random Forest model used by {prog_name} "
        "to classify sample pairs (for advanced users)",
    )

    for cur_parser in (
        search_conta_parser,
        easy_wf_parser,
        plot_conta_parser,
        test_install_parser,
        train_model_parser
    ):
        cur_parser.add_argument(
            "-s",
            dest="species_ab_table_fp",
            type=readable_file,
            required=cur_parser != test_install_parser,
            metavar="SPECIES_ABUNDANCE_TABLE",
            help=argparse.SUPPRESS
            if cur_parser == test_install_parser
            else "Input TSV file corresponding to the species abundance table",
        )
        cur_parser.add_argument(
            "-s2",
            dest="species_ab_table_2_fp",
            type=readable_file,
            required=False,
            metavar="SPECIES_ABUNDANCE_TABLE_2",
            help=argparse.SUPPRESS
            if cur_parser in (test_install_parser, train_model_parser)
            else "Optional input TSV file corresponding to another species abundance table. "
            "If provided, samples from this table will be considered as contamination targets "
            "while those from the first table as contamination sources.",
        )
        cur_parser.add_argument(
            "--filter-low-ab",
            dest="filtering_ab_thr_factor",
            type=float,
            required=False,
            default=None,
            metavar="AB_THRESHOLD_FACTOR",
            help=argparse.SUPPRESS
            if cur_parser in (test_install_parser, train_model_parser)
            else "Filter out low-abundance species that may be inaccurately quantified. "
            "In each sample, set the abundance of species to zero if they are up to "
            "%(metavar)s times more abundant than the least abundant species. "
            "Recommended value for MetaPhlAn4: 20 (default: None)",
        )

    for cur_parser in search_conta_parser, easy_wf_parser, test_install_parser:
        cur_parser.add_argument(
            "-m",
            dest="rf_model_fp",
            type=readable_file,
            required=False,
            default=search_conta_defaults.MODEL_FILE,
            metavar="RF_MODEL_FILE",
            help=argparse.SUPPRESS
            if cur_parser != search_conta_parser
            else "Joblib file containing the pre-trained Random Forest model "
            "used to detect contamination events (default: %(default)s)",
        )
        cur_parser.add_argument(
            "--probability-cutoff",
            dest="probability_cutoff",
            type=bounded_float_01,
            default=search_conta_defaults.PROBABILITY_CUTOFF,
            metavar="PROBABILITY_CUTOFF",
            help=argparse.SUPPRESS
            if cur_parser != search_conta_parser
            else "Only report contamination events with a probability greater than "
            "%(metavar)s (default: %(default).2f)",
        )
        cur_parser.add_argument(
            "--rate-cutoff",
            dest="rate_cutoff",
            type=bounded_float_01,
            default=search_conta_defaults.RATE_CUTOFF,
            metavar="RATE_CUTOFF",
            help=argparse.SUPPRESS
            if cur_parser != search_conta_parser
            else "Only report events with a contamination rate greater than "
            "%(metavar)s (default: %(default).0f)",
        )
        cur_parser.add_argument(
            "-c",
            dest="conta_events_fp",
            type=writable_file,
            required=cur_parser != test_install_parser,
            metavar="CONTAMINATION_EVENTS_FILE",
            help=argparse.SUPPRESS
            if cur_parser == test_install_parser
            else "Output TSV file listing all contamination events",
        )

    plot_conta_parser.add_argument(
        "-c",
        dest="conta_events_fp",
        type=readable_file,
        required=True,
        metavar="CONTAMINATION_EVENTS_FILE",
        help="Input TSV file listing all contaminations events.",
    )

    for cur_parser in plot_conta_parser, easy_wf_parser, test_install_parser:
        cur_parser.add_argument(
            "-r",
            dest="pdf_report_fp",
            type=writable_file,
            required=cur_parser != test_install_parser,
            metavar="PDF_REPORT_FILE",
            help=argparse.SUPPRESS
            if cur_parser == test_install_parser
            else "Output PDF file with scatterplots for all contamination events",
        )

    train_model_parser.add_argument(
        "-m",
        dest="model_fp",
        type=writable_file,
        required=True,
        metavar="MODEL_FILE",
        help="Output file storing the trained Random Forest model",
    )
    train_model_parser.add_argument(
        "-r",
        dest="json_report_fp",
        type=writable_file,
        required=True,
        metavar="JSON_REPORT_FILE",
        help="Output JSON file storing classification performance metrics "
        "for train and test splits",
    )
    train_model_parser.add_argument(
        "--test-size",
        dest="test_size",
        type=bounded_float_01,
        default=train_model_defaults.TEST_SIZE,
        metavar="TEST_SIZE",
        help="Proportion of dataset to include in test split (default: %(default).2f)",
    )
    train_model_parser.add_argument(
        "--ntrees",
        dest="ntrees",
        type=int,
        default=train_model_defaults.NTREES,
        metavar="NTREES",
        help="Number of trees in the RandomForest model (default: %(default)d)",
    )
    train_model_parser.add_argument(
        "--rng-seed",
        dest="rng_seed",
        type=int,
        default=train_model_defaults.RNG_SEED,
        metavar="RNG_SEED",
        help="Seed of the random number generator for reproducibility (default: %(default)d)",
    )

    for cur_parser in (
        search_conta_parser,
        easy_wf_parser,
        test_install_parser,
        train_model_parser,
    ):
        cur_parser.add_argument(
            "--nproc",
            dest="nproc",
            type=nproc,
            default=1
            if cur_parser == test_install_parser
            else multiprocessing.cpu_count(),
            help=argparse.SUPPRESS
            if cur_parser == test_install_parser
            else "Number of parallel processes to search contaminations (default: %(default)d)",
        )

    for cur_parser in plot_conta_parser, easy_wf_parser, test_install_parser:
        cur_parser.add_argument(
            "--nrow",
            dest="nrow",
            type=int,
            choices=range(plot_conta_defaults.MIN_NROW, plot_conta_defaults.MAX_NROW),
            default=plot_conta_defaults.NROW,
            metavar="NROW",
            help=argparse.SUPPRESS
            if cur_parser != plot_conta_parser
            else "Number of scatterplots to draw vertically on each page (default: %(default)d)",
        )
        cur_parser.add_argument(
            "--ncol",
            dest="ncol",
            type=int,
            choices=range(plot_conta_defaults.MIN_NCOL, plot_conta_defaults.MAX_NCOL),
            default=plot_conta_defaults.NCOL,
            metavar="NCOL",
            help=argparse.SUPPRESS
            if cur_parser != plot_conta_parser
            else "Number of scatterplots to draw horizontally on each page (default: %(default)d)",
        )
        cur_parser.add_argument(
            "--no-conta-line",
            dest="no_conta_line",
            action="store_true",
            help=argparse.SUPPRESS
            if cur_parser != plot_conta_parser
            else "Do not show contamination line in scatterplots",
        )
        cur_parser.add_argument(
            "--color-conta-species",
            dest="color_conta_species",
            action="store_true",
            help=argparse.SUPPRESS
            if cur_parser != plot_conta_parser
            else "Use a different color for species introduced by contamination",
        )

    test_install_parser.add_argument(
        "--keep-results",
        dest="keep_results",
        action="store_true",
        help="Keep all temporary results files.",
    )

    return parser.parse_args(args=sys.argv[1:] or ["--help"])


def main() -> None:
    set_logging()
    args = get_arguments()

    if args.command == "test_install":
        test_install = TestInstall(args.keep_results)
        args.species_ab_table_fp = test_install.species_ab_table_fp
        args.conta_events_fp = test_install.conta_events_fp
        args.pdf_report_fp = test_install.pdf_report_fp

    # Add comment line in output file describing execution context
    if args.command in ("easy_wf", "search_conta"):
        exec_desc = ExecutionDescription(
            args.species_ab_table_fp,
            args.species_ab_table_2_fp,
            args.rf_model_fp,
            args.filtering_ab_thr_factor,
            args.probability_cutoff,
            args.rate_cutoff,
        )

        with open(args.conta_events_fp, "w", encoding="utf8") as conta_events_fh:
            print(exec_desc, file=conta_events_fh)

    # Load first species abundance table
    with open(args.species_ab_table_fp, "r", encoding="utf8") as species_ab_table_fh:
        species_ab_table = ab_table_utils.read_filter_normalize(
            species_ab_table_fh,
            args.filtering_ab_thr_factor,
        )

    # Load second species abundance table if necessary
    species_ab_table_2 = None
    if (
        args.species_ab_table_2_fp is not None
        and args.species_ab_table_fp != args.species_ab_table_2_fp
    ):
        with open(args.species_ab_table_2_fp, "r", encoding="utf8") as species_ab_table_2_fh:
            species_ab_table_2 = ab_table_utils.read_filter_normalize(
                species_ab_table_2_fh,
                args.filtering_ab_thr_factor,
            )
        ab_table_utils.compare_species_names(
            species_ab_table,
            species_ab_table_2,
        )

    # Search and save contamination events
    if args.command in ("easy_wf", "search_conta", "test_install"):
        with open(args.rf_model_fp, "rb") as rf_model_fh:
            conta_events = run_search_conta(
                species_ab_table,
                species_ab_table_2,
                rf_model_fh,
                args.probability_cutoff,
                args.rate_cutoff,
                args.nproc,
            )
        with open(args.conta_events_fp, "a", encoding="utf8") as conta_events_fh:
            ContaminationEventIO.write_tsv(conta_events, conta_events_fh)

    # Or load them
    if args.command == "plot_conta":
        with open(args.conta_events_fp, "r", encoding="utf8") as conta_events_fh:
            conta_events = ContaminationEventIO.read_tsv(conta_events_fh)

    # Create PDF with scatterplots
    if args.command in ("plot_conta", "easy_wf", "test_install"):
        with open(args.pdf_report_fp, "wb") as pdf_report_fh:
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

    if args.command in ("easy_wf", "search_conta") and conta_events:
        logging.warning("Contamination events may be false positives, especially "
                            "when dealing with samples with similar species abundance profiles "
                            "(longitudinal data, animals raised together)")
        if args.command == "search_conta":
            logging.warning("Check each reported contamination event by inspecting "
                            "scatterplots generated by the plot_conta subcommand")
        else:
            logging.warning("Check each reported contamination event by inspecting "
                "scatterplots in the PDF report %s", args.pdf_report_fp)

    if args.command == "test_install":
        test_install.check_results()

    if args.command == "train_model":
        with open(args.model_fp, "wb") as model_fh, \
             open(args.json_report_fp, "w", encoding="utf8") as json_report_fh:
            run_train_model(
                species_ab_table,
                model_fh,
                json_report_fh,
                args.test_size,
                args.ntrees,
                args.rng_seed,
                args.nproc,
            )

if __name__ == "__main__":
    main()
