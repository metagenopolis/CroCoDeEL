"""End-to-end installation test for CroCoDeEL."""

import filecmp
import importlib.resources
import logging
from pathlib import Path
from tempfile import NamedTemporaryFile
from typing import Final

import pandas as pd

from crocodeel import ab_table_utils
from crocodeel.conta_event import ContaminationEventIO
from crocodeel.exceptions import SelfTestError
from crocodeel.plot_conta import Defaults as plot_conta_defaults
from crocodeel.plot_conta import run_plot_conta
from crocodeel.search_conta import Defaults as search_conta_defaults
from crocodeel.search_conta import run_search_conta


class TestData:
    """Files and expected results used by the self-test suite"""

    SPECIES_ABUNDANCE_TABLE: Final[Path] = Path(
        str(
            importlib.resources.files().joinpath(
                "test_data",
                "mgs_profiles_test.tsv",
            )
        )
    )

    EXPECTED_CONTA_EVENTS_FILE: Final[Path] = Path(
        str(
            importlib.resources.files().joinpath(
                "test_data",
                "results",
                "contamination_events.tsv",
            )
        )
    )

    EXPECTED_PDF_REPORT_FILE: Final[Path] = Path(
        str(
            importlib.resources.files().joinpath(
                "test_data",
                "results",
                "contamination_events.pdf",
            )
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


class SelfTest:
    """Run an end-to-end self-test of CroCoDeEL."""

    def __init__(self, keep_results: bool = False) -> None:
        self.keep_results = keep_results

    def run(self) -> None:
        """Run the complete installation test on the toy dataset."""
        logging.info("Running tests on the toy dataset")

        species_ab_table = self._load_test_data()

        conta_events_fp = self._create_temporary_file(".tsv")
        pdf_report_fp = self._create_temporary_file(".pdf")

        try:
            self._run_search_conta(
                species_ab_table,
                conta_events_fp,
            )

            self._run_plot_conta(
                species_ab_table,
                conta_events_fp,
                pdf_report_fp,
            )

            self._check_results(
                conta_events_fp,
                pdf_report_fp,
            )

        finally:
            if not self.keep_results:
                conta_events_fp.unlink(missing_ok=True)
                pdf_report_fp.unlink(missing_ok=True)
                logging.info("Temporary result test files deleted")

        logging.info("Tests completed successfully")

    @staticmethod
    def _load_test_data() -> pd.DataFrame:
        """Load and normalize the species abundance table used for testing."""
        with TestData.SPECIES_ABUNDANCE_TABLE.open(
            "r",
            encoding="utf8",
        ) as species_ab_table_fh:
            return ab_table_utils.read_filter_normalize(
                species_ab_table_fh,
                filtering_ab_thr_factor=None,
            )

    @staticmethod
    def _create_temporary_file(suffix: str) -> Path:
        """Create an empty temporary file and return its path."""
        mode = "wb" if suffix == ".pdf" else "w"

        with NamedTemporaryFile(
            mode=mode,
            prefix="crocodeel_test_",
            suffix=suffix,
            delete=False,
        ) as fh:
            return Path(fh.name)

    @staticmethod
    def _run_search_conta(
        species_ab_table: pd.DataFrame,
        conta_events_fp: Path,
    ) -> None:
        """Run contamination detection and save the resulting events."""
        with (
            search_conta_defaults.MODEL_FILE.open("rb") as rf_model_fh,
            conta_events_fp.open("w", encoding="utf8") as conta_events_fh,
        ):
            conta_events = run_search_conta(
                species_ab_table,
                None,
                rf_model_fh,
                search_conta_defaults.PROBABILITY_CUTOFF,
                search_conta_defaults.RATE_CUTOFF,
                nproc=1,
            )

            ContaminationEventIO.write_tsv(
                conta_events,
                conta_events_fh,
            )

    @staticmethod
    def _run_plot_conta(
        species_ab_table: pd.DataFrame,
        conta_events_fp: Path,
        pdf_report_fp: Path,
    ) -> None:
        """Read contamination events and generate the PDF report."""
        with (
            conta_events_fp.open("r", encoding="utf8") as conta_events_fh,
            pdf_report_fp.open("wb") as pdf_report_fh,
        ):
            conta_events = ContaminationEventIO.read_tsv(
                conta_events_fh,
            )

            run_plot_conta(
                species_ab_table,
                None,
                conta_events,
                pdf_report_fh,
                plot_conta_defaults.NROW,
                plot_conta_defaults.NCOL,
                no_conta_line=False,
                color_conta_species=False,
            )

    @staticmethod
    def _check_results(
        conta_events_fp: Path,
        pdf_report_fp: Path,
    ) -> None:
        """Check that generated results match the expected results."""
        if not filecmp.cmp(
            TestData.EXPECTED_CONTA_EVENTS_FILE,
            conta_events_fp,
            shallow=False,
        ):
            raise SelfTestError(
                "Contamination events found are not those expected."
            )

        logging.info(
            "All contamination events with expected rates and probabilities found"
        )

        pdf_report_size = pdf_report_fp.stat().st_size

        if pdf_report_size < TestData.MIN_PDF_REPORT_SIZE:
            raise SelfTestError(
                "PDF report appears too small: "
                f"size is {pdf_report_size} bytes "
                f"(expected around {TestData.EXPECTED_PDF_REPORT_SIZE} bytes)."
            )

        if pdf_report_size > TestData.MAX_PDF_REPORT_SIZE:
            raise SelfTestError(
                "PDF report appears too large: "
                f"size is {pdf_report_size} bytes "
                f"(expected around {TestData.EXPECTED_PDF_REPORT_SIZE} bytes)."
            )

        logging.info("PDF report size is within the expected range")
