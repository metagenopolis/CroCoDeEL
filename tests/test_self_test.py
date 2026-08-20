"""Unit tests for the CroCoDeEL self-test suite."""

from pathlib import Path
from unittest.mock import MagicMock, patch

import pandas as pd
import pytest

from crocodeel.conta_event import ContaminationEvent
from crocodeel.exceptions import SelfTestError
from crocodeel.plot_conta import Defaults as plot_conta_defaults
from crocodeel.search_conta import Defaults as search_conta_defaults
from crocodeel.self_test import SelfTest, TestData


# ---------------------------------------------------------------------------
# Test data loading
# ---------------------------------------------------------------------------


def test_load_test_data() -> None:
    """Test that the self-test suite loads the expected abundance table."""
    species_ab_table = SelfTest._load_test_data()

    assert species_ab_table.shape == (1990, 14)


# ---------------------------------------------------------------------------
# Temporary file creation
# ---------------------------------------------------------------------------


def test_create_temporary_file() -> None:
    """Test creation of temporary files with the requested suffix."""
    txt_fp = SelfTest._create_temporary_file(".tsv")
    pdf_fp = SelfTest._create_temporary_file(".pdf")

    try:
        with txt_fp.open("w", encoding="utf8") as fh:
            fh.write("test")

        with pdf_fp.open("wb") as fh:
            fh.write(b"test")

        assert txt_fp.read_text(encoding="utf8") == "test"
        assert pdf_fp.read_bytes() == b"test"
    finally:
        # _create_temporary_file() creates files outside pytest's tmp_path,
        # so the test must clean them up explicitly.
        txt_fp.unlink(missing_ok=True)
        pdf_fp.unlink(missing_ok=True)


# ---------------------------------------------------------------------------
# SelfTest.run()
# ---------------------------------------------------------------------------


def test_run_deletes_temporary_results(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """Test that temporary result files are deleted after a successful run."""
    conta_events_fp = tmp_path / "contamination_events.tsv"
    pdf_report_fp = tmp_path / "contamination_report.pdf"

    conta_events_fp.touch()
    pdf_report_fp.touch()

    temporary_files = iter(
        [
            conta_events_fp,
            pdf_report_fp,
        ]
    )

    monkeypatch.setattr(
        SelfTest,
        "_create_temporary_file",
        staticmethod(lambda suffix: next(temporary_files)),
    )
    monkeypatch.setattr(
        SelfTest,
        "_load_test_data",
        staticmethod(lambda: pd.DataFrame()),
    )
    monkeypatch.setattr(
        SelfTest,
        "_run_search_conta",
        staticmethod(
            lambda species_ab_table, conta_events_fp: None,
        ),
    )
    monkeypatch.setattr(
        SelfTest,
        "_run_plot_conta",
        staticmethod(
            lambda species_ab_table, conta_events_fp, pdf_report_fp: None,
        ),
    )
    monkeypatch.setattr(
        SelfTest,
        "_check_results",
        staticmethod(
            lambda conta_events_fp, pdf_report_fp: None,
        ),
    )

    SelfTest(keep_results=False).run()

    assert not conta_events_fp.exists()
    assert not pdf_report_fp.exists()


def test_run_keeps_temporary_results(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """Test that temporary result files are kept when requested."""
    conta_events_fp = tmp_path / "contamination_events.tsv"
    pdf_report_fp = tmp_path / "contamination_report.pdf"

    conta_events_fp.touch()
    pdf_report_fp.touch()

    temporary_files = iter(
        [
            conta_events_fp,
            pdf_report_fp,
        ]
    )

    monkeypatch.setattr(
        SelfTest,
        "_create_temporary_file",
        staticmethod(lambda suffix: next(temporary_files)),
    )
    monkeypatch.setattr(
        SelfTest,
        "_load_test_data",
        staticmethod(lambda: pd.DataFrame()),
    )
    monkeypatch.setattr(
        SelfTest,
        "_run_search_conta",
        staticmethod(
            lambda species_ab_table, conta_events_fp: None,
        ),
    )
    monkeypatch.setattr(
        SelfTest,
        "_run_plot_conta",
        staticmethod(
            lambda species_ab_table, conta_events_fp, pdf_report_fp: None,
        ),
    )
    monkeypatch.setattr(
        SelfTest,
        "_check_results",
        staticmethod(
            lambda conta_events_fp, pdf_report_fp: None,
        ),
    )

    SelfTest(keep_results=True).run()

    assert conta_events_fp.exists()
    assert pdf_report_fp.exists()

    # keep_results=True transfers cleanup responsibility to the caller.
    conta_events_fp.unlink()
    pdf_report_fp.unlink()


def test_run_cleans_up_after_failure(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """Test that temporary result files are deleted when the workflow fails."""
    conta_events_fp = tmp_path / "contamination_events.tsv"
    pdf_report_fp = tmp_path / "contamination_report.pdf"

    conta_events_fp.touch()
    pdf_report_fp.touch()

    temporary_files = iter(
        [
            conta_events_fp,
            pdf_report_fp,
        ]
    )

    def fail_search(species_ab_table, conta_events_fp):
        raise SelfTestError("Self-test failed")

    monkeypatch.setattr(
        SelfTest,
        "_create_temporary_file",
        staticmethod(lambda suffix: next(temporary_files)),
    )
    monkeypatch.setattr(
        SelfTest,
        "_load_test_data",
        staticmethod(lambda: pd.DataFrame()),
    )
    monkeypatch.setattr(
        SelfTest,
        "_run_search_conta",
        staticmethod(fail_search),
    )
    monkeypatch.setattr(
        SelfTest,
        "_run_plot_conta",
        staticmethod(
            lambda species_ab_table, conta_events_fp, pdf_report_fp: None,
        ),
    )
    monkeypatch.setattr(
        SelfTest,
        "_check_results",
        staticmethod(
            lambda conta_events_fp, pdf_report_fp: None,
        ),
    )

    with pytest.raises(
        SelfTestError,
        match="Self-test failed",
    ):
        SelfTest(keep_results=False).run()

    assert not conta_events_fp.exists()
    assert not pdf_report_fp.exists()


# ---------------------------------------------------------------------------
# SelfTest._check_results()
# ---------------------------------------------------------------------------


def test_check_results_accepts_expected_results(tmp_path: Path) -> None:
    """Test that expected contamination results and PDF size are accepted."""
    conta_events_fp = tmp_path / "contamination_events.tsv"
    pdf_report_fp = tmp_path / "contamination_events.pdf"

    conta_events_fp.write_bytes(
        TestData.EXPECTED_CONTA_EVENTS_FILE.read_bytes()
    )
    pdf_report_fp.write_bytes(
        TestData.EXPECTED_PDF_REPORT_FILE.read_bytes()
    )

    SelfTest._check_results(
        conta_events_fp,
        pdf_report_fp,
    )


def test_check_results_rejects_different_contamination_events(
    tmp_path: Path,
) -> None:
    """Test that different contamination events are rejected."""
    conta_events_fp = tmp_path / "contamination_events.tsv"
    pdf_report_fp = tmp_path / "contamination_events.pdf"

    conta_events_fp.write_text(
        "different contamination events\n",
        encoding="utf8",
    )
    pdf_report_fp.write_bytes(
        TestData.EXPECTED_PDF_REPORT_FILE.read_bytes()
    )

    with pytest.raises(
        SelfTestError,
        match="Contamination events found are not those expected",
    ):
        SelfTest._check_results(
            conta_events_fp,
            pdf_report_fp,
        )


def test_check_results_rejects_pdf_that_is_too_small(
    tmp_path: Path,
) -> None:
    """Test that a PDF below the accepted size range is rejected."""
    conta_events_fp = tmp_path / "contamination_events.tsv"
    pdf_report_fp = tmp_path / "contamination_events.pdf"

    conta_events_fp.write_bytes(
        TestData.EXPECTED_CONTA_EVENTS_FILE.read_bytes()
    )
    pdf_report_fp.write_bytes(
        b"0" * (TestData.MIN_PDF_REPORT_SIZE - 1)
    )

    with pytest.raises(
        SelfTestError,
        match="PDF report appears too small",
    ):
        SelfTest._check_results(
            conta_events_fp,
            pdf_report_fp,
        )


def test_check_results_accepts_minimum_pdf_size(
    tmp_path: Path,
) -> None:
    """Test that the minimum accepted PDF size is allowed."""
    conta_events_fp = tmp_path / "contamination_events.tsv"
    pdf_report_fp = tmp_path / "contamination_events.pdf"

    conta_events_fp.write_bytes(
        TestData.EXPECTED_CONTA_EVENTS_FILE.read_bytes()
    )
    pdf_report_fp.write_bytes(
        b"0" * TestData.MIN_PDF_REPORT_SIZE
    )

    SelfTest._check_results(
        conta_events_fp,
        pdf_report_fp,
    )


def test_check_results_rejects_pdf_that_is_too_large(
    tmp_path: Path,
) -> None:
    """Test that a PDF above the accepted size range is rejected."""
    conta_events_fp = tmp_path / "contamination_events.tsv"
    pdf_report_fp = tmp_path / "contamination_events.pdf"

    conta_events_fp.write_bytes(
        TestData.EXPECTED_CONTA_EVENTS_FILE.read_bytes()
    )
    pdf_report_fp.write_bytes(
        b"0" * (TestData.MAX_PDF_REPORT_SIZE + 1)
    )

    with pytest.raises(
        SelfTestError,
        match="PDF report appears too large",
    ):
        SelfTest._check_results(
            conta_events_fp,
            pdf_report_fp,
        )


def test_check_results_accepts_maximum_pdf_size(
    tmp_path: Path,
) -> None:
    """Test that the maximum accepted PDF size is allowed."""
    conta_events_fp = tmp_path / "contamination_events.tsv"
    pdf_report_fp = tmp_path / "contamination_events.pdf"

    conta_events_fp.write_bytes(
        TestData.EXPECTED_CONTA_EVENTS_FILE.read_bytes()
    )
    pdf_report_fp.write_bytes(
        b"0" * TestData.MAX_PDF_REPORT_SIZE
    )

    SelfTest._check_results(
        conta_events_fp,
        pdf_report_fp,
    )


# ---------------------------------------------------------------------------
# SelfTest._run_search_conta() and _run_plot_conta()
# ---------------------------------------------------------------------------


@pytest.fixture
def species_ab_table() -> pd.DataFrame:
    """Return a small abundance table for testing."""
    return pd.DataFrame(
        {
            "sample1": [0.1, 0.01],
            "sample2": [0.2, 0.02],
        },
        index=["species1", "species2"],
    )


def test_run_search_conta(
    species_ab_table: pd.DataFrame,
    tmp_path: Path,
) -> None:
    """Test that _run_search_conta searches and writes contamination events."""
    conta_events_fp = tmp_path / "contamination_events.tsv"

    conta_event = ContaminationEvent(
        source="source",
        target="target",
        rate=0.1,
        probability=0.9,
        conta_line_species=[],
    )
    expected_events = [conta_event]

    rf_model_fh = MagicMock()
    model_fp = MagicMock()
    model_fp.open.return_value.__enter__.return_value = rf_model_fh

    with (
        patch(
            "crocodeel.self_test.search_conta_defaults.MODEL_FILE",
            model_fp,
        ),
        patch(
            "crocodeel.self_test.run_search_conta",
            return_value=expected_events,
        ) as mock_run_search,
        patch(
            "crocodeel.self_test.ContaminationEventIO.write_tsv",
        ) as mock_write_tsv,
    ):
        SelfTest._run_search_conta(
            species_ab_table,
            conta_events_fp,
        )

    model_fp.open.assert_called_once_with("rb")

    mock_run_search.assert_called_once_with(
        species_ab_table,
        None,
        rf_model_fh,
        search_conta_defaults.PROBABILITY_CUTOFF,
        search_conta_defaults.RATE_CUTOFF,
        nproc=1,
    )

    mock_write_tsv.assert_called_once()

    write_args = mock_write_tsv.call_args.args

    assert write_args[0] == expected_events
    assert write_args[1].name == str(conta_events_fp)


def test_run_plot_conta(
    species_ab_table: pd.DataFrame,
    tmp_path: Path,
) -> None:
    """Test that _run_plot_conta reads events and generates the PDF report."""
    conta_events_fp = tmp_path / "contamination_events.tsv"
    pdf_report_fp = tmp_path / "contamination_report.pdf"

    conta_event = ContaminationEvent(
        source="source",
        target="target",
        rate=0.1,
        probability=0.9,
        conta_line_species=[],
    )
    expected_events = [conta_event]

    # The contents are irrelevant because read_tsv() is mocked.
    conta_events_fp.write_text(
        "dummy",
        encoding="utf8",
    )

    with (
        patch(
            "crocodeel.self_test.ContaminationEventIO.read_tsv",
            return_value=expected_events,
        ) as mock_read_tsv,
        patch(
            "crocodeel.self_test.run_plot_conta",
        ) as mock_run_plot,
    ):
        SelfTest._run_plot_conta(
            species_ab_table,
            conta_events_fp,
            pdf_report_fp,
        )

    mock_read_tsv.assert_called_once()

    read_args = mock_read_tsv.call_args.args

    assert read_args[0].name == str(conta_events_fp)

    mock_run_plot.assert_called_once()

    plot_args = mock_run_plot.call_args.args

    assert plot_args[0] is species_ab_table
    assert plot_args[1] is None
    assert plot_args[2] == expected_events
    assert plot_args[3].name == str(pdf_report_fp)
    assert plot_args[4] == plot_conta_defaults.NROW
    assert plot_args[5] == plot_conta_defaults.NCOL

    assert mock_run_plot.call_args.kwargs == {
        "no_conta_line": False,
        "color_conta_species": False,
    }
