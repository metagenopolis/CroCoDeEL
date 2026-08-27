"""Unit tests for contamination plotting."""

import logging
from io import BytesIO
from unittest.mock import MagicMock, patch

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest

from crocodeel.conta_event import ContaminationEvent
from crocodeel.exceptions import InputDataError
from crocodeel.plot_conta import ContaminationPlotsReport, run_plot_conta

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


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
def conta_event_with_line_species() -> ContaminationEvent:
    """Return a contamination event with a contamination-line species."""
    return ContaminationEvent(
        source="source",
        target="target",
        rate=0.1,
        probability=0.9,
        conta_line_species=["species1"],
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


# ---------------------------------------------------------------------------
# ContaminationPlotsReport constructor
# ---------------------------------------------------------------------------


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

    pd.testing.assert_frame_equal(
        species_ab_table,
        original,
    )


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


def test_constructor_rejects_unknown_source_sample(
    species_ab_table: pd.DataFrame,
    conta_event: ContaminationEvent,
) -> None:
    """Test that an unknown source sample is rejected."""
    conta_event.source = "unknown_source"

    with pytest.raises(
        InputDataError,
        match="unknown source samples: unknown_source",
    ):
        ContaminationPlotsReport(
            species_ab_table,
            [conta_event],
            nrow=1,
            ncol=1,
            no_conta_line=False,
            color_conta_species=False,
        )


def test_constructor_rejects_unknown_target_sample(
    species_ab_table: pd.DataFrame,
    conta_event: ContaminationEvent,
) -> None:
    """Test that an unknown target sample is rejected."""
    conta_event.target = "unknown_target"

    with pytest.raises(
        InputDataError,
        match="unknown target samples: unknown_target",
    ):
        ContaminationPlotsReport(
            species_ab_table,
            [conta_event],
            nrow=1,
            ncol=1,
            no_conta_line=False,
            color_conta_species=False,
        )


def test_constructor_rejects_unknown_source_and_target_samples(
    species_ab_table: pd.DataFrame,
    conta_event: ContaminationEvent,
) -> None:
    """Test that unknown source and target samples are both reported."""
    conta_event.source = "unknown_source"
    conta_event.target = "unknown_target"

    with pytest.raises(
        InputDataError,
        match=(
            "unknown source samples: unknown_source; "
            "unknown target samples: unknown_target"
        ),
    ):
        ContaminationPlotsReport(
            species_ab_table,
            [conta_event],
            nrow=1,
            ncol=1,
            no_conta_line=False,
            color_conta_species=False,
        )


def test_constructor_truncates_unknown_sample_names(
    species_ab_table: pd.DataFrame,
) -> None:
    """Test that long lists of unknown sample names are truncated."""
    conta_events = [
        ContaminationEvent(
            source=f"unknown_source_{i}",
            target="target",
            rate=0.1,
            probability=0.9,
        )
        for i in range(6)
    ]

    with pytest.raises(InputDataError) as exc_info:
        ContaminationPlotsReport(
            species_ab_table,
            conta_events,
            nrow=1,
            ncol=1,
            no_conta_line=False,
            color_conta_species=False,
        )

    message = str(exc_info.value)

    assert "6 sample names (showing 5):" in message
    assert "unknown_source_0" in message
    assert "unknown_source_4" in message
    assert "unknown_source_5" not in message


# ---------------------------------------------------------------------------
# _get_species_to_plot()
# ---------------------------------------------------------------------------


def test_get_species_to_plot(
    report: ContaminationPlotsReport,
    conta_event: ContaminationEvent,
) -> None:
    """Test that species present in either sample are selected."""
    species_to_plot = report._get_species_to_plot(conta_event)

    assert species_to_plot.tolist() == [False, True, True]


# ---------------------------------------------------------------------------
# _set_scatter_colors()
# ---------------------------------------------------------------------------


def test_set_scatter_colors_enabled(
    report: ContaminationPlotsReport,
    conta_event_with_line_species: ContaminationEvent,
) -> None:
    """Test that contamination-line species are highlighted."""
    scatterplot = MagicMock()

    species_to_plot = pd.Series(
        [True, True, True],
        index=report.species_ab_table.index,
    )

    report.color_conta_species = True

    report._set_scatter_colors(
        scatterplot,
        conta_event_with_line_species,
        species_to_plot,
    )

    scatterplot.set_edgecolor.assert_called_once()

    colors = scatterplot.set_edgecolor.call_args.args[0]

    assert colors == [
        "orange",
        "black",
        "black",
    ]


def test_set_scatter_colors_disabled(
    report: ContaminationPlotsReport,
    conta_event_with_line_species: ContaminationEvent,
) -> None:
    """Test that all species are black when coloring is disabled."""
    scatterplot = MagicMock()

    species_to_plot = pd.Series(
        [True, True, True],
        index=report.species_ab_table.index,
    )

    report.color_conta_species = False

    report._set_scatter_colors(
        scatterplot,
        conta_event_with_line_species,
        species_to_plot,
    )

    scatterplot.set_edgecolor.assert_called_once_with(
        "black",
    )


# ---------------------------------------------------------------------------
# _add_reference_lines()
# ---------------------------------------------------------------------------


def test_add_reference_lines_adds_two_lines(
    report: ContaminationPlotsReport,
    conta_event: ContaminationEvent,
) -> None:
    """Test that both identity and contamination lines are added."""
    fig, ax = plt.subplots()

    try:
        report._add_reference_lines(
            ax,
            conta_event,
        )

        assert len(ax.lines) == 2

        identity_line = ax.lines[0]
        contamination_line = ax.lines[1]

        assert identity_line.get_linestyle() == "-"
        assert contamination_line.get_linestyle() == "-"

        # Both reference lines have slope 1 in log-transformed abundance space.
        for line in (identity_line, contamination_line):
            xdata = line.get_xdata()
            ydata = line.get_ydata()

            assert np.isclose(
                (ydata[1] - ydata[0])
                / (xdata[1] - xdata[0]),
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
        report._add_reference_lines(
            ax,
            conta_event,
        )

        assert len(ax.lines) == 1
    finally:
        plt.close(fig)


# ---------------------------------------------------------------------------
# _set_abundance_axis()
# ---------------------------------------------------------------------------


def test_set_abundance_axis(
    report: ContaminationPlotsReport,
) -> None:
    """Test abundance-axis tick locations."""
    fig, ax = plt.subplots()

    try:
        report._set_abundance_axis(ax)

        expected_ticks = np.arange(
            report.pseudo_zero,
            1,
        )

        np.testing.assert_array_equal(
            ax.get_xticks(),
            expected_ticks,
        )
        np.testing.assert_array_equal(
            ax.get_yticks(),
            expected_ticks,
        )
    finally:
        plt.close(fig)


# ---------------------------------------------------------------------------
# _create_plot()
# ---------------------------------------------------------------------------


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
        report._create_plot(
            conta_event,
            ax,
        )

        assert ax.get_xlabel() == "sample_target"
        assert ax.get_ylabel() == "sample_source"
        assert ax.get_title() == "prob = 0.9, rate = 10.0%"
    finally:
        plt.close(fig)


# ---------------------------------------------------------------------------
# _create_page()
# ---------------------------------------------------------------------------


def test_create_page_hides_unused_axes(
    report: ContaminationPlotsReport,
    conta_event: ContaminationEvent,
) -> None:
    """Test that unused axes are hidden when fewer events than panels exist."""
    report.nrow = 2
    report.ncol = 2
    report.conta_events = [conta_event]

    pdf = MagicMock()
    fig = MagicMock()

    axes = np.empty(
        (2, 2),
        dtype=object,
    )
    axes[0, 0] = MagicMock()
    axes[0, 1] = MagicMock()
    axes[1, 0] = MagicMock()
    axes[1, 1] = MagicMock()

    with (
        patch.object(
            report,
            "_create_plot",
        ) as mock_create_plot,
        patch(
            "crocodeel.plot_conta.plt.subplots",
            return_value=(fig, axes),
        ),
        patch(
            "crocodeel.plot_conta.plt.tight_layout",
        ),
        patch(
            "crocodeel.plot_conta.plt.close",
        ),
    ):
        report._create_page(
            page_id=0,
            pdf=pdf,
        )

    mock_create_plot.assert_called_once_with(
        conta_event,
        axes[0, 0],
    )

    # Only the first of four panels contains a contamination event.
    axes[0, 1].axis.assert_called_once_with("off")
    axes[1, 0].axis.assert_called_once_with("off")
    axes[1, 1].axis.assert_called_once_with("off")

    pdf.savefig.assert_called_once_with(fig)


# ---------------------------------------------------------------------------
# save_to_pdf()
# ---------------------------------------------------------------------------


def test_save_to_pdf_creates_non_empty_pdf(
    report: ContaminationPlotsReport,
) -> None:
    """Test that a PDF report is generated."""
    pdf_fh = BytesIO()

    report.save_to_pdf(pdf_fh)

    assert pdf_fh.tell() > 0
    assert pdf_fh.getvalue().startswith(b"%PDF")


# ---------------------------------------------------------------------------
# run_plot_conta()
# ---------------------------------------------------------------------------


def test_run_plot_conta(
    species_ab_table: pd.DataFrame,
    sample_species_ab_table: pd.DataFrame,
    conta_event: ContaminationEvent,
) -> None:
    """Test run_plot_conta() with a second abundance table."""
    pdf_report_fh = MagicMock()
    pdf_report_fh.name = "report.pdf"

    with (
        patch(
            "crocodeel.plot_conta.ContaminationPlotsReport",
        ) as mock_report,
        patch(
            "crocodeel.plot_conta.perf_counter",
            side_effect=[0.0, 1.5],
        ),
    ):
        run_plot_conta(
            species_ab_table=species_ab_table,
            species_ab_table_2=sample_species_ab_table,
            conta_events=[conta_event],
            pdf_report_fh=pdf_report_fh,
            nrow=2,
            ncol=3,
            no_conta_line=True,
            color_conta_species=True,
        )

    mock_report.assert_called_once()

    report_args = mock_report.call_args.args
    report_table = report_args[0]

    missing_species = species_ab_table.index.difference(
        sample_species_ab_table.index,
    )

    assert len(missing_species) > 0

    # Species absent from the second abundance table must be represented
    # as -inf, corresponding to zero abundance after log10 transformation.
    for species in missing_species:
        assert (
            report_table.loc[
                species,
                sample_species_ab_table.columns,
            ]
            == -np.inf
        ).all()

    assert report_args[1:] == (
        [conta_event],
        2,
        3,
        True,
        True,
    )

    mock_report.return_value.save_to_pdf.assert_called_once_with(
        pdf_report_fh,
    )
