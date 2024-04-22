# -*- coding: utf-8 -*-

from pathlib import Path
import joblib
from sklearn.ensemble import RandomForestClassifier


class RandomForestModel:
    @staticmethod
    def get_path() -> Path:
        model_path = (
            Path(__file__).resolve().parent / "models" / "crocodeel_last_version.joblib"
        )
        return model_path
    
    @staticmethod
    def load() -> RandomForestClassifier:
        return joblib.load(RandomForestModel.get_path())
