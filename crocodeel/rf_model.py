import importlib.resources
from pathlib import Path
from typing import Final
import warnings
import joblib
from sklearn.exceptions import InconsistentVersionWarning
from sklearn.ensemble import RandomForestClassifier


class RandomForestModel:
    RF_MODEL_FILE: Final[Path] = Path(
        importlib.resources.files().joinpath("models", "crocodeel_rf_Mar2023.joblib")
    )

    @staticmethod
    def load() -> RandomForestClassifier:
        with warnings.catch_warnings():
            warnings.filterwarnings(action="ignore", category=InconsistentVersionWarning)
            return joblib.load(RandomForestModel.RF_MODEL_FILE)

    @staticmethod
    def get_version() -> str:
        return RandomForestModel.RF_MODEL_FILE.stem
