from socket import gethostname
from getpass import getuser
from datetime import datetime
from importlib.metadata import version
from pathlib import Path
from typing import Optional
from crocodeel.rf_model import RandomForestModel

class ExecutionDescription:
    def __init__(self, species_ab_table_path : str, species_ab_table_2_path : Optional[str] = None) -> None:
        self.software_version = version('crocodeel')
        self.rf_model_version = RandomForestModel.get_version()
        self.hostname = gethostname()
        self.username = getuser()
        self.datetime = datetime.now().replace(microsecond=0).isoformat()
        self.species_ab_table_path = Path(species_ab_table_path).resolve()
        self.species_ab_table_2_path = (
            Path(species_ab_table_2_path).resolve() if species_ab_table_2_path else None
        )

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

        return exec_desc_str
