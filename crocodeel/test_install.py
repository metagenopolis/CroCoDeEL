import importlib.resources
from pathlib import Path
import logging
from tempfile import NamedTemporaryFile
import filecmp
import sys
from typing import Final
from crocodeel.easy_wf import run_easy_wf


class TestData:
    SPECIES_ABUNDANCE_TABLE: Final[Path] = Path(
        importlib.resources.files().joinpath("test_data", "mgs_profiles_test.tsv")
    )

    EXPECTED_CONTA_EVENTS_FILE: Final[Path] = Path(
        importlib.resources.files().joinpath("test_data", "results", "contamination_events.tsv")
    )

    EXPECTED_PDF_REPORT_FILE: Final[Path] = Path(
        importlib.resources.files().joinpath("test_data", "results", "contamination_events.pdf")
    )


def run_test_install(keep_results: bool) -> None:
    with (
        open(TestData.SPECIES_ABUNDANCE_TABLE, "r", encoding="utf-8") as species_ab_table_fh,
        NamedTemporaryFile(
            mode="w+",
            prefix="contamination_events_",
            suffix=".tsv",
            delete=not keep_results,
            delete_on_close=False,
        ) as conta_events_fh,
        NamedTemporaryFile(
            mode="wb",
            prefix="contamination_events_",
            suffix=".pdf",
            delete=not keep_results,
            delete_on_close=False,
        ) as pdf_report_fh,
    ):
        easy_wf_args = {
            "species_ab_table_fh": species_ab_table_fh,
            "conta_events_fh": conta_events_fh,
            "pdf_report_fh": pdf_report_fh,
            "nproc": 1,
        }

        logging.info("Tests on the toy dataset started")
        run_easy_wf(easy_wf_args)

        if filecmp.cmp(TestData.EXPECTED_CONTA_EVENTS_FILE, conta_events_fh.name):
            logging.info("All contamination events with expected rates and probabilities found")
        else:
            logging.error("Contamination events found are not those expected")
            sys.exit(1)

    if not keep_results:
        logging.info("Temporary result files deleted")
    
    logging.info("Tests completed successfully")
