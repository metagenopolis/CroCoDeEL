import sys
import logging
import importlib.resources
from pathlib import Path
from typing import Final
from tempfile import NamedTemporaryFile
import filecmp


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


class TestInstall:
    def __init__(self, keep_results: bool) -> None:
        self.keep_results = keep_results
        self.species_ab_table_fh = open(
            TestData.SPECIES_ABUNDANCE_TABLE, "r", encoding="utf-8"
        )
        self.conta_events_fh = NamedTemporaryFile(
            mode="w",
            prefix="contamination_events_",
            suffix=".tsv",
            delete=False,
            delete_on_close=False,
        )
        self.pdf_report_fh = NamedTemporaryFile(
            mode="wb",
            prefix="contamination_events_",
            suffix=".pdf",
            delete=False,
            delete_on_close=False,
        )

    def __del__(self) -> None:
        if not self.species_ab_table_fh.closed:
            self.species_ab_table_fh.close()
        if not self.conta_events_fh.closed:
            self.conta_events_fh.close()
        if not self.conta_events_fh.closed:
            self.pdf_report_fh.close()

        if not self.keep_results:
            Path(self.pdf_report_fh.name).unlink()
            Path(self.conta_events_fh.name).unlink()
            logging.info("Temporary result test files deleted")

    def check_results(self) -> None:
        if filecmp.cmp(TestData.EXPECTED_CONTA_EVENTS_FILE, self.conta_events_fh.name):
            logging.info(
                "All contamination events with expected rates and probabilities found"
            )
        else:
            logging.error("Contamination events found are not those expected")
            sys.exit(1)

        logging.info("Tests completed successfully")
