#!/bin/env python3
# -*- coding: utf-8 -*-

import sys
import argparse
import logging
import multiprocessing
from pathlib import Path
import filecmp
from tempfile import NamedTemporaryFile
from importlib.metadata import version
from crocodeel.execution_description import ExecutionDescription
from crocodeel import ab_table_utils
from crocodeel.conta_event import ContaminationEventIO
from crocodeel.search_conta import run_search_conta
from crocodeel.plot_conta import run_plot_conta, Defaults as plot_conta_defaults
from crocodeel.ressources import TestData


def set_logging() -> None:
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    logging.basicConfig(format="%(asctime)s :: %(levelname)s :: %(message)s")


def nproc(value) -> int:
    max_nproc = multiprocessing.cpu_count()

    try:
        value = int(value)
    except ValueError as value_err:
        raise argparse.ArgumentTypeError("NPROC is not an integer") from value_err

    if value <= 0:
        raise argparse.ArgumentTypeError("minimum NPROC is 1")
    if value > max_nproc:
        raise argparse.ArgumentTypeError(f"maximum NPROC is {max_nproc}")

    return value


def get_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(prog="crocodeel")
    parser.add_argument(
       "-v", "--version", action="version", version=f"%(prog)s version { version("crocodeel")}"
    )

    subparsers = parser.add_subparsers(
        title="positional arguments",
        help="Select command",
        dest="command",
        required=True,
    )
    easy_wf_parser = subparsers.add_parser(
        "easy_wf",
        help="Search cross-sample contamination and "
        "create a PDF report in one command",
    )
    search_conta_parser = subparsers.add_parser(
        "search_conta", help="Search cross-sample contamination"
    )
    plot_conta_parser = subparsers.add_parser(
        "plot_conta",
        help="Create a PDF report with scatterplots representing "
        "species abundance profiles for each contamination event",
    )
    test_install_parser = subparsers.add_parser(
        "test_install",
        help="Test if %(prog)s is correctly installed "
        "and generates expected results",
    )

    for cur_parser in (
        search_conta_parser,
        easy_wf_parser,
        plot_conta_parser,
        test_install_parser,
    ):
        cur_parser.add_argument(
            "-s",
            dest="species_ab_table_fh",
            type=argparse.FileType("r"),
            required=cur_parser != test_install_parser,
            default=str(TestData.SPECIES_ABUNDANCE_TABLE)
            if cur_parser == test_install_parser
            else None,
            metavar="SPECIES_ABUNDANCE_TABLE",
            help=argparse.SUPPRESS
            if cur_parser == test_install_parser
            else "Input TSV file corresponding to the species abundance table",
        )
        cur_parser.add_argument(
            "-s2",
            dest="species_ab_table_fh_2",
            type=argparse.FileType("r"),
            required=False,
            metavar="SPECIES_ABUNDANCE_TABLE_2",
            help=argparse.SUPPRESS
            if cur_parser == test_install_parser
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
            if cur_parser == test_install_parser
            else "Filter out low-abundance species that may be inaccurately quantified. "
            "In each sample, set the abundance of species to zero if they are up to "
            "%(metavar)s times more abundant than the least abundant species .  "
            "Recommended value for MetaPhlAn4: 20 (default: None)",
        )

    for cur_parser in search_conta_parser, easy_wf_parser, test_install_parser:
        cur_parser.add_argument(
            "-c",
            dest="conta_events_fh",
            type=argparse.FileType("w"),
            required=cur_parser != test_install_parser,
            default=NamedTemporaryFile(
                mode="w",
                prefix="contamination_events_",
                suffix=".tsv",
                delete=False,
                delete_on_close=False,
            )
            if cur_parser == test_install_parser
            else None,
            metavar="CONTAMINATION_EVENTS_FILE",
            help=argparse.SUPPRESS
            if cur_parser == test_install_parser
            else "Output TSV file listing all contamination events",
        )
    plot_conta_parser.add_argument(
        "-c",
        dest="conta_events_fh",
        type=argparse.FileType("r"),
        required=True,
        metavar="CONTAMINATION_EVENTS_FILE",
        help="Input TSV file listing all contaminations events.",
    )

    for cur_parser in plot_conta_parser, easy_wf_parser, test_install_parser:
        cur_parser.add_argument(
            "-r",
            dest="pdf_report_fh",
            type=argparse.FileType("wb"),
            required=cur_parser != test_install_parser,
            default=NamedTemporaryFile(
                mode="wb",
                prefix="contamination_events_",
                suffix=".pdf",
                delete=False,
                delete_on_close=False,
            )
            if cur_parser == test_install_parser
            else None,
            metavar="PDF_REPORT_FILE",
            help=argparse.SUPPRESS
            if cur_parser == test_install_parser
            else "Output PDF file with scatterplots for all contamination events",
        )

    for cur_parser in search_conta_parser, easy_wf_parser, test_install_parser:
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

    if args.command == 'test_install':
        logging.info("Running tests on the toy dataset")

    # Add comment line in output file describing execution context
    if args.command in ("easy_wf", "search_conta"):
        exec_desc = ExecutionDescription(
            args.species_ab_table_fh,
            args.species_ab_table_fh_2,
            args.filtering_ab_thr_factor,
        )
        print(exec_desc, file=args.conta_events_fh)

    # Load first species abundance table
    species_ab_table = ab_table_utils.read_filter_normalize(
        args.species_ab_table_fh,
        args.filtering_ab_thr_factor,
    )
    args.species_ab_table_fh.close()

    # Load second species abundance table if necessary
    species_ab_table_2 = None
    if args.species_ab_table_fh_2 is not None:
        if args.species_ab_table_fh.name != args.species_ab_table_fh_2.name:
            species_ab_table_2 = ab_table_utils.read_filter_normalize(
                args.species_ab_table_fh_2,
                args.filtering_ab_thr_factor,
            )
            ab_table_utils.compare_species_names(
                species_ab_table,
                species_ab_table_2,
            )
        args.species_ab_table_fh_2.close()

    # Search and save contamination events
    if args.command in ("easy_wf", "search_conta", "test_install"):
        conta_events = run_search_conta(
            species_ab_table, species_ab_table_2, args.nproc
        )
        ContaminationEventIO.write_tsv(conta_events, args.conta_events_fh)
        args.conta_events_fh.close()

    # Or load them
    if args.command == "plot_conta":
        conta_events = ContaminationEventIO.read_tsv(args.conta_events_fh)
        args.conta_events_fh.close()

    # Create PDF with scatterplots
    if args.command in ("plot_conta", "easy_wf", "test_install"):
        run_plot_conta(
            species_ab_table,
            species_ab_table_2,
            conta_events,
            args.pdf_report_fh,
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
                "scatterplots in the PDF report %s", Path(args.pdf_report_fh.name).resolve())

    if args.command == "test_install":
        if filecmp.cmp(TestData.EXPECTED_CONTA_EVENTS_FILE, args.conta_events_fh.name):
            logging.info(
                "All contamination events with expected rates and probabilities found"
            )
        else:
            logging.error("Contamination events found are not those expected")
            sys.exit(1)

        if not args.keep_results:
            Path(args.pdf_report_fh.name).unlink()
            Path(args.conta_events_fh.name).unlink()
            logging.info("Temporary result files deleted")

        logging.info("Tests completed successfully")


if __name__ == "__main__":
    main()
