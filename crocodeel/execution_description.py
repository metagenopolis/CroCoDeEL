from socket import gethostname
from getpass import getuser
from datetime import datetime
from importlib.metadata import version
from pathlib import Path
from typing import Optional


class ExecutionDescription:

    def __init__(
        self,
        species_ab_table_fp: Path,
        species_ab_table_2_fp: Optional[Path],
        rf_model_fp: Path,
        filtering_ab_thr_factor: Optional[float],
        probability_cutoff: float,
        rate_cutoff: float
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
        exec_desc_str = (
            f"# crocodeel version: {self.software_version}"
            f" | rf_model: {self.rf_model_fp}"
            f" | hostname: {self.hostname}"
            f" | username: {self.username}"
            f" | datetime: {self.datetime}"
            f" | species_ab_table: {self.species_ab_table_fp}"
        )

        if self.species_ab_table_2_fp:
            exec_desc_str = (
                exec_desc_str + f" | species_ab_table_2: {self.species_ab_table_2_fp}"
            )

        exec_desc_str = (
            exec_desc_str
            + f" | filtering_ab_thr_factor: {self.filtering_ab_thr_factor}"
        )

        exec_desc_str += (
            f" | probability_cutoff: {self.probability_cutoff}"
            f" | rate_cutoff: {self.rate_cutoff}"
        )

        return exec_desc_str
