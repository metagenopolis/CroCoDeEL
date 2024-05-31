#!/bin/env python3
# -*- coding: utf-8 -*-

import sys
import argparse
import logging
import multiprocessing
from importlib.metadata import version
from crocodeel.execution_description import ExecutionDescription
from crocodeel.search_conta import run_search_conta, run_search_conta_distrib
from crocodeel.plot_conta import run_plot_conta, Defaults as plot_conta_defaults
from crocodeel.easy_wf import run_easy_wf
from crocodeel.test_install import run_test_install


def set_logging() -> None:
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    logging.basicConfig(format='%(asctime)s :: %(levelname)s :: %(message)s')


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
        "create a PDF report in one command.",
    )
    easy_wf_parser.add_argument(
        "-s",
        dest="species_ab_table_fh",
        type=argparse.FileType("r"),
        required=True,
        metavar='SPECIES_ABUNDANCE_TABLE',
        help="Input TSV file giving the species abundance profiles in metagenomic samples",
    )
    easy_wf_parser.add_argument(
       "-c",
        dest="conta_events_fh",
        type=argparse.FileType("w+"),
        required=True,
        metavar='CONTAMINATION_EVENTS_FILE',
        help="Output TSV file listing all contamination events",
    )
    easy_wf_parser.add_argument(
        "-r",
        dest="pdf_report_fh",
        type=argparse.FileType("wb"),
        required=True,
        metavar='PDF_REPORT_FILE',
        help="Output PDF file with scatterplots for all contamination events",
    )
    easy_wf_parser.add_argument(
        "--nproc",
        dest="nproc",
        type=nproc,
        default=multiprocessing.cpu_count(),
        help="Number of parallel processes to search for contaminations "
        "(default: %(default)d)",
    )

    search_conta_distrib_parser = subparsers.add_parser(
        "search_conta", help="Search cross-sample contamination"
    )
    search_conta_distrib_parser.add_argument(
        "-s",
        dest="species_ab_table_fh",
        type=argparse.FileType("r"),
        required=True,
        metavar='SPECIES_ABUNDANCE_TABLE',
        help="Input TSV file giving the species abundance profiles in metagenomic samples.",
    )
    search_conta_distrib_parser.add_argument(
        "-c",
        dest="conta_events_fh",
        type=argparse.FileType("w"),
        required=True,
        metavar='CONTAMINATION_EVENTS_FILE',
        help="Output TSV file listing all contamination events.",
    )
    search_conta_distrib_parser.add_argument(
        "--nproc",
        dest="nproc",
        type=nproc,
        default=multiprocessing.cpu_count(),
        help="Number of parallel processes to search for contaminations "
        "(default: %(default)d)",
    )

    search_conta_distrib_parser = subparsers.add_parser(
        "search_conta_distrib", help="Search cross-sample contamination. "
        "Command used in workflows to distribute computational load across multiple nodes."
    )
    search_conta_distrib_parser.add_argument(
        "-s1",
        dest="species_ab_table_fh",
        type=argparse.FileType("r"),
        required=True,
        metavar='SPECIES_ABUNDANCE_TABLE',
        help="Input TSV file giving the species abundance profiles in metagenomic samples."
        " These samples are potential contamination sources.",
    )
    search_conta_distrib_parser.add_argument(
        "-s2",
        dest="species_ab_table_fh_2",
        type=argparse.FileType("r"),
        required=True,
        metavar='SPECIES_ABUNDANCE_TABLE_2',
        help="Input TSV file giving the species abundance profiles in metagenomic samples."
        " These samples are potential contamination targets.",
    )
    search_conta_distrib_parser.add_argument(
        "-c",
        dest="conta_events_fh",
        type=argparse.FileType("w"),
        required=True,
        metavar='CONTAMINATION_EVENTS_FILE',
        help="Output TSV file listing all contamination events.",
    )
    search_conta_distrib_parser.add_argument(
        "--nproc",
        dest="nproc",
        type=nproc,
        default=multiprocessing.cpu_count(),
        help="Number of parallel processes to search for contaminations "
        "(default: %(default)d)",
    )

    plot_conta_parser = subparsers.add_parser(
        "plot_conta",
        help="Create a PDF report with scatterplots representing "
        "species abundance profiles for each contamination event.",
    )
    plot_conta_parser.add_argument(
        "-s",
        dest="species_ab_table_fh",
        type=argparse.FileType("r"),
        required=True,
        metavar='SPECIES_ABUNDANCE_TABLE',
        help="Input TSV file giving the species abundance profiles in metagenomic samples.",
    )
    plot_conta_parser.add_argument(
        "-c",
        dest="conta_events_fh",
        type=argparse.FileType("r"),
        required=True,
        metavar='CONTAMINATION_EVENTS_FILE',
        help="Input TSV file listing all contaminations events.",
    )
    plot_conta_parser.add_argument(
        "-r",
        dest="pdf_report_fh",
        type=argparse.FileType("wb"),
        required=True,
        metavar='PDF_REPORT_FILE',
        help="Output PDF file with scatterplots for all contamination events.",
    )
    plot_conta_parser.add_argument(
        "--nrow",
        dest="nrow",
        type=int,
        choices=range(plot_conta_defaults.MIN_NROW, plot_conta_defaults.MAX_NROW),
        default=plot_conta_defaults.NROW,
        metavar="NROW",
        help="Number of scatterplots to draw vertically on each page " "(default: %(default)d)",
    )
    plot_conta_parser.add_argument(
        "--ncol",
        dest="ncol",
        type=int,
        choices=range(plot_conta_defaults.MIN_NCOL, plot_conta_defaults.MAX_NCOL),
        default=plot_conta_defaults.NCOL,
        metavar="NCOL",
        help="Number of scatterplots to draw horizontally on each page "
        "(default: %(default)d)",
    )
    plot_conta_parser.add_argument(
        "--no-conta-line",
        dest="no_conta_line",
        action="store_true",
        help="Do not show contamination line in scatterplots.",
    )
    plot_conta_parser.add_argument(
        "--color-conta-species",
        dest="color_conta_species",
        action="store_true",
        help="Use a different color for species introduced by contamination.",
    )

    test_install_parser = subparsers.add_parser(
        "test_install", help="Test if %(prog)s is correctly installed "
        "and generates expected results"
    )
    test_install_parser.add_argument(
        "--keep",
        dest="keep_results",
        action="store_true",
        help="Keep all temporary results files.",
    )

    return parser.parse_args(args=sys.argv[1:] or ["--help"])


def main() -> None:
    set_logging()
    args = get_arguments()

    # Add comment line in output file describing execution context
    if args.command in ("search_conta","easy_wf"):
        exec_desc = ExecutionDescription(args.species_ab_table_fh.name)
        print(exec_desc, file = args.conta_events_fh)

    if args.command == "search_conta_distrib":
        exec_desc = ExecutionDescription(args.species_ab_table_fh.name,
                                         args.species_ab_table_fh_2.name)
        print(exec_desc, file = args.conta_events_fh)

    if args.command == "easy_wf":
        run_easy_wf(vars(args))
    elif args.command == "search_conta":
        run_search_conta(vars(args))
    elif args.command == "search_conta_distrib":
        if args.species_ab_table_fh.name == args.species_ab_table_fh_2.name:
            args.species_ab_table_fh_2.close()
            run_search_conta(vars(args))
        else:
            run_search_conta_distrib(vars(args))
    elif args.command == "plot_conta":
        run_plot_conta(vars(args))
    elif args.command == "test_install":
        run_test_install(args.keep_results)

if __name__ == "__main__":
    main()
