from dataclasses import dataclass, field
import csv
from typing import TextIO
from pathlib import Path
import logging


@dataclass
class ContaminationEvent:
    source: str
    target: str
    rate: float = field(default=0.0)
    probability: float = field(default=0.0)
    contamination_specific_species: list[str] = field(default_factory=lambda: [])


class ContaminationEventIO:
    @staticmethod
    def read_tsv(fh: TextIO) -> list[ContaminationEvent]:
        # Ugly hack to skip comment lines
        pos = 0
        while fh :
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
                row["contamination_specific_species"].split(","),
            )
            for row in tsv_reader
        ]

        logging.info(
            "%d contamination events loaded from %s",
            len(conta_events),
            Path(fh.name).resolve(),
        )

        return conta_events

    @staticmethod
    def write_tsv(conta_events: list[ContaminationEvent], fh: TextIO) -> None:
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
                        str(conta_event.rate),
                        str(conta_event.probability),
                        ",".join(conta_event.contamination_specific_species),
                    ]
                ),
                file=fh,
            )

        logging.info(
            "Contamination events saved in %s",
            Path(fh.name).resolve(),
        )
