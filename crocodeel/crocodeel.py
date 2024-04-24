#!/bin/env python3
# -*- coding: utf-8 -*-

import sys
import argparse
import logging
import multiprocessing
from utils import load_species_ab_table
from conta_event import ContaminationEventIO
from search_conta import ContaminationSearcherDriver
from plot_conta import ContaminationPlotsReport


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

    subparsers = parser.add_subparsers(
        title="positional arguments",
        help="Select subcommand",
        dest="command",
        required=True,
    )

    search_conta_parser = subparsers.add_parser(
        "search_conta", help="Search cross-sample contamination"
    )
    search_conta_parser.add_argument(
        "-s",
        dest="species_ab_table",
        type=argparse.FileType("r"),
        required=True,
        help="Input TSV file giving the species abundance profiles in metagenomic samples.",
    )
    search_conta_parser.add_argument(
        "--nproc",
        dest="nproc",
        type=nproc,
        default=multiprocessing.cpu_count(),
        help="Number of parallel processes to search for contaminations "
        "(default: %(default)d)",
    )
    search_conta_parser.add_argument(
        "-o",
        dest="output_file",
        type=argparse.FileType("w"),
        required=True,
        help="Output TSV file listing all contaminations events.",
    )

    plot_conta_parser = subparsers.add_parser(
        "plot_conta",
        help="Create a PDF report with scatterplots representing "
        "species abundance profiles for each contamination events.",
    )
    plot_conta_parser.add_argument(
        "-s",
        dest="species_ab_table",
        type=argparse.FileType("r"),
        required=True,
        help="TSV file giving the abundance of species across metagenomic samples.",
    )
    plot_conta_parser.add_argument(
        "-c",
        dest="crocodeel_results",
        type=argparse.FileType("r"),
        required=True,
        help="TSV file generated by CroCoDeEL listing "
        "all detected contamination events.",
    )
    plot_conta_parser.add_argument(
        "--nrow",
        dest="nrow",
        type=int,
        choices=range(1, 11),
        default=4,
        metavar="NROW",
        help="Number of scatterplots to draw vertically on each page "
        "(default: %(default)d)",
    )
    plot_conta_parser.add_argument(
        "--ncol",
        dest="ncol",
        type=int,
        choices=range(1, 11),
        default=4,
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
    plot_conta_parser.add_argument(
        "-o",
        dest="output_file",
        type=argparse.FileType("wb"),
        required=True,
        help="Path to the output PDF file.",
    )

    return parser.parse_args(args=sys.argv[1:] or ["--help"])


def main() -> None:
    set_logging()
    args = get_arguments()
    if args.command == "search_conta":
        species_ab_table = load_species_ab_table(args.species_ab_table)
        args.species_ab_table.close()

        conta_events = ContaminationSearcherDriver(
            species_ab_table,
            nproc=args.nproc
        ).search_contamination()
        contaminated_samples = {conta_event.target for conta_event in conta_events}
        logging.info("%d contamination events detected", len(conta_events))
        logging.info("%d/%d samples contaminated", len(contaminated_samples), species_ab_table.shape[1])

        conta_events.sort(key=lambda e: e.rate, reverse=True)
        ContaminationEventIO.write_tsv(conta_events, args.output_file)
        logging.info("Contamination events sorted by decreasing contamination rate saved in %s", args.output_file.name)
        args.output_file.close()

        logging.warning("Contamination events may be false positives, "
                        "especially when dealing with samples with similar species abundance profiles "
                        "(longitudinal data, animals raised together)")
        logging.warning("Run the plot_conta subcommand to visualize "
                        "and check each reported contamination events")
        

    elif args.command == "plot_conta":
        species_ab_table = load_species_ab_table(args.species_ab_table)
        args.species_ab_table.close()

        conta_events = list(ContaminationEventIO.read_tsv(args.crocodeel_results))
        args.crocodeel_results.close()
        logging.info("%d contamination events loaded", len(conta_events))

        ContaminationPlotsReport(
            species_ab_table=species_ab_table,
            conta_events=conta_events,
            nrow=args.nrow,
            ncol=args.ncol,
            no_contamination_line=args.no_conta_line,
            color_contamination_specific_species=args.color_conta_species,
        ).save_to_pdf(args.output_file)
        logging.info("PDF report saved in %s", args.output_file.name)
        args.output_file.close()


if __name__ == "__main__":
    main()
