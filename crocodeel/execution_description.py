from socket import gethostname
from getpass import getuser
from datetime import datetime
from importlib.metadata import version
from pathlib import Path
from typing import Optional, TextIO
from crocodeel.ressources import RandomForestModel


class ExecutionDescription:

    def __init__(
        self,
        species_ab_table_fh: TextIO,
        species_ab_table_2_fh: Optional[TextIO],
        filtering_ab_thr_factor: Optional[float],
    ) -> None:
        self.software_version = version("crocodeel")
        self.rf_model_version = RandomForestModel.get_version()
        self.hostname = gethostname()
        self.username = getuser()
        self.datetime = datetime.now().replace(microsecond=0).isoformat()
        self.species_ab_table_path = Path(species_ab_table_fh.name).resolve()
        self.species_ab_table_2_path = (
            Path(species_ab_table_2_fh.name).resolve()
            if species_ab_table_2_fh
            else None
        )
        self.filtering_ab_thr_factor = filtering_ab_thr_factor

    def __str__(self) -> str:
        exec_desc_str = (
            f"# crocodeel version: {self.software_version}"
            f" | rf_model_version: {self.rf_model_version}"
            f" | hostname: {self.hostname}"
            f" | username: {self.username}"
            f" | datetime: {self.datetime}"
            f" | species_ab_table: {self.species_ab_table_path}"
        )

        if self.species_ab_table_2_path:
            exec_desc_str = (
                exec_desc_str + f" | species_ab_table_2: {self.species_ab_table_2_path}"
            )

        exec_desc_str = (
            exec_desc_str
            + f" | filtering_ab_thr_factor: {self.filtering_ab_thr_factor}"
        )

        return exec_desc_str
