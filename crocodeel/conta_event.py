"""Data structures and I/O utilities for contamination events."""

import csv
import logging
from dataclasses import dataclass, field
from typing import Final, Iterator, Sequence, TextIO

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
    Round a contamination rate expressed as a percentage.

    The contamination rate must be strictly positive. Unlike standard
    rounding, this adjusts the decimal precision based on the magnitude
    of the rate to ensure detail is preserved for small values.
    """
    if rate <= 0:
        raise ValueError("Contamination rate must be strictly positive.")

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

        if not 0 < rate <= 1:
            raise InputDataError(
                f"Invalid contamination rate '{rate}' on line {row_number}: "
                "value must be greater than 0 and at most 1."
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
    def _prepend_line(
        line: str | None,
        fh: TextIO,
    ) -> Iterator[str]:
        """Yield one previously read line followed by the remaining stream."""
        if line is not None:
            yield line
        yield from fh

    @staticmethod
    def _read_header(
        fh: TextIO,
    ) -> tuple[str | None, int]:
        """Read comments before the header and return the header and line number.
        
        Return ``(None, 0)`` when the file is empty or contains only comments.
        """
        for line_number, line in enumerate(fh, start=1):
            if not line.startswith("#"):
                return line, line_number

        return None, 0

    @staticmethod
    def read_tsv(fh: TextIO) -> list[ContaminationEvent]:
        """Read and validate contamination events from a TSV file."""
        # Skip comment lines before the header without requiring a seekable stream.
        header, header_line_number = ContaminationEventIO._read_header(fh)

        tsv_reader = csv.DictReader(
            ContaminationEventIO._prepend_line(header, fh),
            delimiter="\t",
        )

        ContaminationEventIO._validate_columns(
            tsv_reader.fieldnames,
        )

        conta_events: list[ContaminationEvent] = []

        for row in tsv_reader:
            row_number = header_line_number + tsv_reader.line_num - 1

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
        """Write contamination events to TSV format."""
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
