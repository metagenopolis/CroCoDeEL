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

    EXPECTED_PDF_REPORT_SIZE: Final[int] = EXPECTED_PDF_REPORT_FILE.stat().st_size

    PDF_REPORT_SIZE_TOLERANCE: Final[float] = 0.02

    MIN_PDF_REPORT_SIZE: Final[int] = int(
        (1.0 - PDF_REPORT_SIZE_TOLERANCE) * EXPECTED_PDF_REPORT_SIZE
    )

    MAX_PDF_REPORT_SIZE: Final[int] = int(
        (1.0 + PDF_REPORT_SIZE_TOLERANCE) * EXPECTED_PDF_REPORT_SIZE
    )


class TestInstall:
    def __init__(self, keep_results: bool) -> None:
        logging.info("Running tests on the toy dataset")

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
        if not filecmp.cmp(
            TestData.EXPECTED_CONTA_EVENTS_FILE, self.conta_events_fh.name
        ):
            logging.error("Contamination events found are not those expected")
            sys.exit(1)
        else:
            logging.info(
                "All contamination events with expected rates and probabilities found"
            )

        # Directly comparing PDF contents is unreliable because metadata (e.g., creation date)
        # can cause slight variations. To keep it simple, we verify that the PDF file size
        # is within ±2% of the expected size range.
        pdf_report_size = Path(self.pdf_report_fh.name).stat().st_size
        if pdf_report_size < TestData.MIN_PDF_REPORT_SIZE:
            logging.error(
                "PDF report appears too small: size is %d bytes (expected around %d bytes).",
                pdf_report_size,
                TestData.EXPECTED_PDF_REPORT_SIZE,
            )
            sys.exit(1)
        elif pdf_report_size > TestData.MAX_PDF_REPORT_SIZE:
            logging.error(
                "PDF report appears too large: size is %d bytes (expected around %d bytes).",
                pdf_report_size,
                TestData.EXPECTED_PDF_REPORT_SIZE,
            )
            sys.exit(1)
        else:
            logging.info("PDF report size is within the expected range")

        logging.info("Tests completed successfully")
