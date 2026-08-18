"""Unit tests for contamination plotting."""

from io import BytesIO

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest

from crocodeel.conta_event import ContaminationEvent
from crocodeel.plot_conta import ContaminationPlotsReport


@pytest.fixture
def species_ab_table() -> pd.DataFrame:
    """Return a small log-transformed abundance table for testing."""
    return pd.DataFrame(
        {
            "target": [-np.inf, -2.0, -np.inf],
            "source": [-np.inf, -2.5, -1.0],
        },
        index=["species1", "species2", "species3"],
    )


@pytest.fixture
def sample_species_ab_table() -> pd.DataFrame:
    """Return an abundance table with sample-specific column names."""
    return pd.DataFrame(
        {
            "sample_target": [-np.inf, -2.0],
            "sample_source": [-np.inf, -1.0],
        },
        index=["species1", "species2"],
    )


@pytest.fixture
def conta_event() -> ContaminationEvent:
    """Return a contamination event for testing."""
    return ContaminationEvent(
        source="source",
        target="target",
        rate=0.1,
        probability=0.9,
        conta_line_species=[],
    )


@pytest.fixture
def report(
    species_ab_table: pd.DataFrame,
    conta_event: ContaminationEvent,
) -> ContaminationPlotsReport:
    """Return a contamination plots report for testing."""
    return ContaminationPlotsReport(
        species_ab_table,
        [conta_event],
        nrow=1,
        ncol=1,
        no_conta_line=False,
        color_conta_species=False,
    )


def test_constructor_does_not_modify_species_ab_table(
    species_ab_table: pd.DataFrame,
) -> None:
    """Test that the input abundance table is not modified."""
    original = species_ab_table.copy()

    ContaminationPlotsReport(
        species_ab_table,
        [],
        nrow=1,
        ncol=1,
        no_conta_line=False,
        color_conta_species=False,
    )

    pd.testing.assert_frame_equal(species_ab_table, original)


def test_constructor_computes_pseudo_zero(
    species_ab_table: pd.DataFrame,
) -> None:
    """Test that the pseudo-zero is computed from the minimum finite abundance."""
    report = ContaminationPlotsReport(
        species_ab_table,
        [],
        nrow=1,
        ncol=1,
        no_conta_line=False,
        color_conta_species=False,
    )

    assert report.pseudo_zero == -3


def test_get_species_to_plot_excludes_species_absent_from_both_samples(
    report: ContaminationPlotsReport,
    conta_event: ContaminationEvent,
) -> None:
    """Test that species absent from both samples are excluded."""
    species_to_plot = report._get_species_to_plot(conta_event)

    assert species_to_plot.tolist() == [False, True, True]


def test_get_species_to_plot_includes_species_present_in_either_sample(
    report: ContaminationPlotsReport,
    conta_event: ContaminationEvent,
) -> None:
    """Test that species present in either sample are included."""
    species_to_plot = report._get_species_to_plot(conta_event)

    assert species_to_plot.tolist() == [False, True, True]


def test_add_reference_lines_adds_two_lines(
    report: ContaminationPlotsReport,
    conta_event: ContaminationEvent,
) -> None:
    """Test that both identity and contamination lines are added."""
    fig, ax = plt.subplots()

    try:
        report._add_reference_lines(ax, conta_event)

        assert len(ax.lines) == 2

        identity_line = ax.lines[0]
        contamination_line = ax.lines[1]

        assert identity_line.get_linestyle() == "-"
        assert contamination_line.get_linestyle() == "-"

        for line in (identity_line, contamination_line):
            xdata = line.get_xdata()
            ydata = line.get_ydata()
            assert np.isclose(
                (ydata[1] - ydata[0]) / (xdata[1] - xdata[0]),
                1.0,
            )
    finally:
        plt.close(fig)


def test_add_reference_lines_can_hide_contamination_line(
    species_ab_table: pd.DataFrame,
    conta_event: ContaminationEvent,
) -> None:
    """Test that the contamination line is omitted when requested."""
    report = ContaminationPlotsReport(
        species_ab_table,
        [],
        nrow=1,
        ncol=1,
        no_conta_line=True,
        color_conta_species=False,
    )

    fig, ax = plt.subplots()

    try:
        report._add_reference_lines(ax, conta_event)

        assert len(ax.lines) == 1
    finally:
        plt.close(fig)


def test_create_plot_sets_labels_and_title(
    sample_species_ab_table: pd.DataFrame,
) -> None:
    """Test that a contamination plot has the expected labels and title."""
    conta_event = ContaminationEvent(
        source="sample_source",
        target="sample_target",
        rate=0.1,
        probability=0.9,
        conta_line_species=[],
    )

    report = ContaminationPlotsReport(
        sample_species_ab_table,
        [conta_event],
        nrow=1,
        ncol=1,
        no_conta_line=False,
        color_conta_species=False,
    )

    fig, ax = plt.subplots()

    try:
        report._create_plot(conta_event, ax)

        assert ax.get_xlabel() == "sample_target"
        assert ax.get_ylabel() == "sample_source"
        assert ax.get_title() == "prob = 0.9, rate = 10.0%"
    finally:
        plt.close(fig)


def test_save_to_pdf_creates_non_empty_pdf(
    report: ContaminationPlotsReport,
) -> None:
    """Test that a PDF report is generated."""
    pdf_fh = BytesIO()

    report.save_to_pdf(pdf_fh)

    assert pdf_fh.tell() > 0
    assert pdf_fh.getvalue().startswith(b"%PDF")
