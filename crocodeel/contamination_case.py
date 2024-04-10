from dataclasses import dataclass
import csv
from typing import TextIO, Iterator


@dataclass
class ContaminationCase:
    source: str
    target: str
    rate: float
    probability: float
    contamination_specific_species: list[int]


class ContaminationCaseIO:
    @staticmethod
    def read_tsv(fh: TextIO) -> Iterator[ContaminationCase]:
        tsv_reader = csv.DictReader(fh, delimiter="\t")

        for row in tsv_reader:
            contamination_specific_species = map(int, row["msp_list"].split(","))
            yield ContaminationCase(
                row["source"],
                row["target"],
                float(row["rate"]),
                float(row["probability"]),
                list(contamination_specific_species)
            )
