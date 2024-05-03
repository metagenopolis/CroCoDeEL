import importlib.resources
from pathlib import Path
from typing import Final
import joblib
from sklearn.ensemble import RandomForestClassifier


class RandomForestModel:
    RF_MODEL_FILE: Final[Path] = Path(
        importlib.resources.files().joinpath("models", "crocodeel_rf_Mar2023.joblib")
    )

    @staticmethod
    def load() -> RandomForestClassifier:
        return joblib.load(RandomForestModel.RF_MODEL_FILE)
