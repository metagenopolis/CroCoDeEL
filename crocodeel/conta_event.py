"""Data structures and I/O utilities for contamination events."""

import csv
import logging
from dataclasses import dataclass, field
from typing import Final, Sequence, TextIO

import numpy as np

from crocodeel.exceptions import InputDataError


@dataclass
class ContaminationEvent:
    """Represent a cross-sample contamination event."""

    source: str
    target: str
    rate: float = field(default=0.0)
    probability: float = field(default=0.0)
    conta_line_species: list[str] = field(default_factory=list)


def round_conta_rate(
    rate: float,
    significant_digits: int = 3,
) -> float:
    """
    Round a contamination rate to a specific number of significant digits.

    The contamination rate must be between 0 and 1. Unlike standard
    rounding, this adjusts the decimal precision based on the magnitude
    of the rate to ensure detail is preserved for small values.
    """
    if rate == 0:
        return 0

    magnitude = np.floor(np.log10(rate))
    decimals = int(significant_digits - (magnitude + 1))

    return np.round(rate, decimals)


class ContaminationEventIO:
    """Read and write contamination events in TSV format."""

    REQUIRED_COLUMNS: Final[tuple[str, ...]] = (
        "source",
        "target",
        "rate",
        "probability",
        "contamination_specific_species",
    )

    @staticmethod
    def _validate_columns(
        fieldnames: Sequence[str] | None,
    ) -> None:
        """Validate the columns of a contamination events TSV file."""
        if fieldnames is None:
            raise InputDataError(
                "Contamination events file is empty or has no header."
            )

        missing_columns = [
            column
            for column in ContaminationEventIO.REQUIRED_COLUMNS
            if column not in fieldnames
        ]

        if missing_columns:
            raise InputDataError(
                "Contamination events file is missing required column(s): "
                + ", ".join(missing_columns)
            )

    @staticmethod
    def _parse_float(
        row: dict[str, str],
        column: str,
        row_number: int,
    ) -> float:
        """Parse and validate a finite floating-point value."""
        value = row[column]

        try:
            parsed_value = float(value)
        except (TypeError, ValueError) as error:
            raise InputDataError(
                f"Invalid {column} '{value}' on line {row_number} "
                "of the contamination events file."
            ) from error

        if not np.isfinite(parsed_value):
            raise InputDataError(
                f"Invalid {column} '{value}' on line {row_number} "
                "of the contamination events file: value must be finite."
            )

        return parsed_value

    @staticmethod
    def _validate_event_values(
        source: str,
        target: str,
        rate: float,
        probability: float,
        contamination_specific_species: str,
        row_number: int,
    ) -> None:
        """Validate the values of one contamination event."""
        if not source:
            raise InputDataError(
                f"Missing source sample on line {row_number} "
                "of the contamination events file."
            )

        if not target:
            raise InputDataError(
                f"Missing target sample on line {row_number} "
                "of the contamination events file."
            )

        if not 0 <= rate <= 1:
            raise InputDataError(
                f"Invalid contamination rate '{rate}' on line {row_number}: "
                "value must be between 0 and 1."
            )

        if not 0 <= probability <= 1:
            raise InputDataError(
                f"Invalid contamination probability '{probability}' "
                f"on line {row_number}: value must be between 0 and 1."
            )

        if not contamination_specific_species:
            raise InputDataError(
                "Missing contamination-specific species on line "
                f"{row_number} of the contamination events file."
            )

    @staticmethod
    def read_tsv(fh: TextIO) -> list[ContaminationEvent]:
        """Read and validate contamination events from a TSV file."""
        # Skip comment lines before the header.
        pos = 0

        while True:
            line = fh.readline()

            if not line:
                break

            if not line.startswith("#"):
                break

            pos = fh.tell()

        fh.seek(pos)

        tsv_reader = csv.DictReader(
            fh,
            delimiter="\t",
        )

        ContaminationEventIO._validate_columns(
            tsv_reader.fieldnames,
        )

        conta_events: list[ContaminationEvent] = []

        for row in tsv_reader:
            row_number = tsv_reader.line_num

            if None in row or any(value is None for value in row.values()):
                raise InputDataError(
                    f"Invalid number of fields on line {row_number} "
                    "of the contamination events file."
                )

            source = row["source"]
            target = row["target"]
            contamination_specific_species = row[
                "contamination_specific_species"
            ]

            rate = ContaminationEventIO._parse_float(
                row,
                "rate",
                row_number,
            )
            probability = ContaminationEventIO._parse_float(
                row,
                "probability",
                row_number,
            )

            ContaminationEventIO._validate_event_values(
                source,
                target,
                rate,
                probability,
                contamination_specific_species,
                row_number,
            )

            conta_line_species = contamination_specific_species.split(",")

            conta_events.append(
                ContaminationEvent(
                    source=source,
                    target=target,
                    rate=rate,
                    probability=probability,
                    conta_line_species=conta_line_species,
                )
            )

        logging.info(
            "%d contamination event%s loaded from %s",
            len(conta_events),
            "s" if len(conta_events) > 1 else "",
            fh.name,
        )

        return conta_events

    @staticmethod
    def write_tsv(
        conta_events: list[ContaminationEvent],
        fh: TextIO,
    ) -> None:
        """Write contamination events to a TSV file."""
        print(
            "\t".join(
                [
                    "source",
                    "target",
                    "rate",
                    "probability",
                    "contamination_specific_species",
                ]
            ),
            file=fh,
        )

        for conta_event in conta_events:
            print(
                "\t".join(
                    [
                        conta_event.source,
                        conta_event.target,
                        f"{conta_event.rate:.2e}",
                        str(conta_event.probability),
                        ",".join(conta_event.conta_line_species),
                    ]
                ),
                file=fh,
            )

        logging.info(
            "Contamination events saved in %s",
            fh.name,
        )
