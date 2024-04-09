from dataclasses import dataclass
import csv
from typing import TextIO

@dataclass
class ContaminationCase:
    source: str
    target: str
    rate: float
    probability: float
    contamination_specific_species: list[str]

class ContaminationCaseIO:
    @staticmethod
    def read_tsv(fh: TextIO):
        tsv_reader = csv.DictReader(fh, delimiter="\t")

        for row in tsv_reader:
            contamination_specific_species = row["msp_list"].split(",")
            yield ContaminationCase(
                row["source"],
                row["target"],
                float(row["rate"]),
                float(row["probability"]),
                contamination_specific_species,
            )