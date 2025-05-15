import importlib.resources
from pathlib import Path
from typing import Final


class TestData:
    SPECIES_ABUNDANCE_TABLE: Final[Path] = Path(
        importlib.resources.files().joinpath("test_data", "mgs_profiles_test.tsv")
    )

    EXPECTED_CONTA_EVENTS_FILE: Final[Path] = Path(
        importlib.resources.files().joinpath(
            "test_data", "results", "contamination_events.tsv"
        )
    )

    EXPECTED_PDF_REPORT_FILE: Final[Path] = Path(
        importlib.resources.files().joinpath(
            "test_data", "results", "contamination_events.pdf"
        )
    )
