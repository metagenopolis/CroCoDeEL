from dataclasses import dataclass, field
import csv
from typing import TextIO, Iterator


@dataclass
class ContaminationEvent:
    source: str
    target: str
    rate: float = field(default=0.0)
    probability: float = field(default=0.0)
    contamination_specific_species: list[str] = field(default_factory=lambda: [])


class ContaminationEventIO:
    @staticmethod
    def read_tsv(fh: TextIO) -> Iterator[ContaminationEvent]:
        # Ugly hack to skip comment lines
        pos = 0
        while fh :
            line = fh.readline()
            if not line.startswith("#"):
                break
            pos = fh.tell()
        fh.seek(pos)

        tsv_reader = csv.DictReader(fh, delimiter="\t")

        for row in tsv_reader:
            contamination_specific_species = row[
                "contamination_specific_species"
            ].split(",")

            yield ContaminationEvent(
                row["source"],
                row["target"],
                float(row["rate"]),
                float(row["probability"]),
                contamination_specific_species,
            )

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
