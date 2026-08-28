"""Store metadata describing a CroCoDeEL execution."""

from datetime import datetime
from getpass import getuser
from importlib.metadata import version
from pathlib import Path
from socket import gethostname
from typing import Optional


class ExecutionDescription:
    """Describe the parameters and environment of a CroCoDeEL execution."""

    def __init__(
        self,
        species_ab_table_fp: Path,
        species_ab_table_2_fp: Optional[Path],
        rf_model_fp: Path,
        filtering_ab_thr_factor: Optional[float],
        probability_cutoff: float,
        rate_cutoff: float,
    ) -> None:
        self.software_version = version("crocodeel")
        self.rf_model_fp = rf_model_fp
        self.hostname = gethostname()
        self.username = getuser()
        self.datetime = datetime.now().replace(microsecond=0).isoformat()
        self.species_ab_table_fp = species_ab_table_fp
        self.species_ab_table_2_fp = species_ab_table_2_fp
        self.filtering_ab_thr_factor = filtering_ab_thr_factor
        self.probability_cutoff = probability_cutoff
        self.rate_cutoff = rate_cutoff

    def __str__(self) -> str:
        fields = [
            ("crocodeel version", self.software_version),
            ("rf_model", self.rf_model_fp),
            ("hostname", self.hostname),
            ("username", self.username),
            ("datetime", self.datetime),
            ("species_ab_table", self.species_ab_table_fp),
        ]

        if self.species_ab_table_2_fp:
            fields.append(("species_ab_table_2", self.species_ab_table_2_fp))

        fields += [
            ("filtering_ab_thr_factor", self.filtering_ab_thr_factor),
            ("probability_cutoff", self.probability_cutoff),
            ("rate_cutoff", self.rate_cutoff),
        ]

        return "# " + " | ".join(f"{label}: {value}" for label, value in fields)
