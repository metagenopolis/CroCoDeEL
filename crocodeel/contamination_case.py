from dataclasses import dataclass, field
import csv
from typing import TextIO, Iterator


@dataclass
class ContaminationCase:
    source: str
    target: str
    rate: float = field(default=0.0)
    probability: float = field(default=0.0)
    contamination_specific_species: list[str] = field(default_factory=lambda: [])


class ContaminationCaseIO:
    @staticmethod
    def read_tsv(fh: TextIO) -> Iterator[ContaminationCase]:
        tsv_reader = csv.DictReader(fh, delimiter="\t")

        for row in tsv_reader:
            contamination_specific_species = row[
                "contamination_specific_species"
            ].split(",")

            yield ContaminationCase(
                row["source"],
                row["target"],
                float(row["rate"]),
                float(row["probability"]),
                contamination_specific_species,
            )

    @staticmethod
    def write_tsv(contamination_cases: list[ContaminationCase], fh: TextIO) -> None:
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

        # Write each case
        for contamination_case in contamination_cases:
            print(
                "\t".join(
                    [
                        contamination_case.source,
                        contamination_case.target,
                        str(contamination_case.rate),
                        str(contamination_case.probability),
                        ",".join(contamination_case.contamination_specific_species),
                    ]
                ),
                file=fh,
            )
