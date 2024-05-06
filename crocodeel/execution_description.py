from socket import gethostname
from getpass import getuser
from datetime import datetime
from importlib.metadata import version
from pathlib import Path
from crocodeel.rf_model import RandomForestModel

class ExecutionDescription:
    def __init__(self, species_ab_table_path : str) -> None:
        self.software_version = version('crocodeel')
        self.rf_model_version = RandomForestModel.get_version()
        self.hostname = gethostname()
        self.username = getuser()
        self.datetime = datetime.now().replace(microsecond=0).isoformat()
        self.species_ab_table_path = Path(species_ab_table_path).resolve()

    def __str__(self) -> str:
        return (
            f"# crocodeel version: {self.software_version}"
            f" | rf_model_version: {self.rf_model_version}"
            f" | hostname: {self.hostname}"
            f" | username: {self.username}"
            f" | datetime: {self.datetime}"
            f" | species_ab_table: {self.species_ab_table_path}"
        )
