"""Unit tests for the CroCoDeEL self-test suite."""

import pytest

from crocodeel.self_test import TestData, SelfTest


def test_load_test_data():
    """Test that the self-test suite loads the expected abundance table."""
    species_ab_table = SelfTest._load_test_data()

    assert species_ab_table.shape == (1990, 14)


def test_check_results_accepts_expected_results(tmp_path):
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


def test_check_results_rejects_different_contamination_events(tmp_path):
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
        RuntimeError,
        match="Contamination events found are not those expected",
    ):
        SelfTest._check_results(
            conta_events_fp,
            pdf_report_fp,
        )


def test_check_results_rejects_pdf_that_is_too_small(tmp_path):
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
        RuntimeError,
        match="PDF report appears too small",
    ):
        SelfTest._check_results(
            conta_events_fp,
            pdf_report_fp,
        )


def test_check_results_rejects_pdf_that_is_too_large(tmp_path):
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
        RuntimeError,
        match="PDF report appears too large",
    ):
        SelfTest._check_results(
            conta_events_fp,
            pdf_report_fp,
        )


def test_check_results_accepts_minimum_pdf_size(tmp_path):
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


def test_check_results_accepts_maximum_pdf_size(tmp_path):
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