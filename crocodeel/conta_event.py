"""Data structures and I/O utilities for contamination events."""

from dataclasses import dataclass, field
import csv
from typing import TextIO
import logging


@dataclass
class ContaminationEvent:
    """Represent a cross-sample contamination event."""

    source: str
    target: str
    rate: float = field(default=0.0)
    probability: float = field(default=0.0)
    conta_line_species: list[str] = field(default_factory=list)


class ContaminationEventIO:
    """Read and write contamination events in TSV format."""

    @staticmethod
    def read_tsv(fh: TextIO) -> list[ContaminationEvent]:
        """Read contamination events from a TSV file."""
        # Ugly hack to skip comment lines
        pos = 0
        while fh:
            line = fh.readline()
            if not line.startswith("#"):
                break
            pos = fh.tell()
        fh.seek(pos)

        tsv_reader = csv.DictReader(fh, delimiter="\t")

        conta_events = [
            ContaminationEvent(
                row["source"],
                row["target"],
                float(row["rate"]),
                float(row["probability"]),
                (
                    row["contamination_specific_species"].split(",")
                    if row["contamination_specific_species"]
                    else []
                ),
            )
            for row in tsv_reader
        ]

        logging.info(
            "%d contamination event%s loaded from %s",
            len(conta_events),
            "s" if len(conta_events) > 1 else "",
            fh.name,
        )

        return conta_events

    @staticmethod
    def write_tsv(conta_events: list[ContaminationEvent], fh: TextIO) -> None:
        """Write contamination events to a TSV file."""
        # Write header
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

        # Write each event
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
