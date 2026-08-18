"""Unit tests for the CroCoDeEL command-line interface."""

import argparse

import pytest

from crocodeel.crocodeel import (
    bounded_float_01,
    get_arguments,
    nproc,
    positive_float,
    positive_int,
    readable_file,
    run_easy_workflow,
    writable_file,
)


# ---------------------------------------------------------------------------
# File argument validators
# ---------------------------------------------------------------------------


def test_readable_file(tmp_path):
    """Test that a readable file is accepted."""
    filepath = tmp_path / "input.tsv"
    filepath.write_text("test\n")

    result = readable_file(str(filepath))

    assert result == filepath.resolve()


def test_readable_file_does_not_exist(tmp_path):
    """Test that a missing file is rejected."""
    filepath = tmp_path / "missing.tsv"

    with pytest.raises(argparse.ArgumentTypeError, match="does not exist"):
        readable_file(str(filepath))


def test_readable_file_directory(tmp_path):
    """Test that a directory is rejected."""
    with pytest.raises(
        argparse.ArgumentTypeError,
        match="not a regular file",
    ):
        readable_file(str(tmp_path))


def test_writable_existing_file(tmp_path):
    """Test that an existing writable file is accepted."""
    filepath = tmp_path / "output.tsv"
    filepath.write_text("test\n")

    result = writable_file(str(filepath))

    assert result == filepath.resolve()


def test_writable_file_directory(tmp_path):
    """Test that a directory is rejected as an output file."""
    with pytest.raises(
        argparse.ArgumentTypeError,
        match="is a directory",
    ):
        writable_file(str(tmp_path))


def test_writable_file_missing_parent(tmp_path):
    """Test that a missing parent directory is rejected."""
    filepath = tmp_path / "missing" / "output.tsv"

    with pytest.raises(
        argparse.ArgumentTypeError,
        match="does not exist",
    ):
        writable_file(str(filepath))


def test_writable_file_new_file(tmp_path):
    """Test that a new file in a writable directory is accepted."""
    filepath = tmp_path / "output.tsv"

    result = writable_file(str(filepath))

    assert result == filepath.resolve()


# ---------------------------------------------------------------------------
# Numeric argument validators
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "value, expected",
    [
        ("1", 1.0),
        ("0.1", 0.1),
        ("1.5", 1.5),
        ("1e-3", 0.001),
    ],
)
def test_positive_float_valid(value, expected):
    """Test that positive floating-point values are accepted."""
    assert positive_float(value) == pytest.approx(expected)


@pytest.mark.parametrize(
    "value",
    ["0", "-1", "-0.1"],
)
def test_positive_float_non_positive(value):
    """Test that non-positive values are rejected."""
    with pytest.raises(
        argparse.ArgumentTypeError,
        match="finite number greater than 0",
    ):
        positive_float(value)


@pytest.mark.parametrize(
    "value",
    [
        "abc",
        "1.2.3",
        "",
    ],
)
def test_positive_float_invalid(value):
    """Test that non-numeric values are rejected."""
    with pytest.raises(
        argparse.ArgumentTypeError,
        match=f"{value} is not a number",
    ):
        positive_float(value)


@pytest.mark.parametrize("value", ["nan", "inf", "-inf"])
def test_positive_float_special_values(value):
    """Test that non-finite floating-point values are rejected."""
    with pytest.raises(argparse.ArgumentTypeError):
        positive_float(value)


@pytest.mark.parametrize(
    "value, expected",
    [
        ("1", 1),
        ("10", 10),
        ("999", 999),
    ],
)
def test_positive_int_valid(value, expected):
    """Test that positive integers are accepted."""
    assert positive_int(value) == expected

@pytest.mark.parametrize(
    "value",
    [
        "0",
        "-1",
    ],
)
def test_positive_int_non_positive(value):
    """Test that non-positive integers are rejected."""
    with pytest.raises(
        argparse.ArgumentTypeError,
        match="value must be greater than 0",
    ):
        positive_int(value)


@pytest.mark.parametrize(
    "value",
    [
        "abc",
        "",
        "1.5",
    ],
)
def test_positive_int_invalid(value):
    """Test that non-integer values are rejected."""
    with pytest.raises(
        argparse.ArgumentTypeError,
        match=f"{value} is not an integer",
    ):
        positive_int(value)

@pytest.mark.parametrize(
    "value",
    ["0", "-1", "abc"],
)
def test_nproc_invalid(value):
    """Test that invalid process counts are rejected."""
    with pytest.raises(argparse.ArgumentTypeError):
        nproc(value)


def test_nproc_too_many():
    """Test that more processes than available CPUs are rejected."""
    with pytest.raises(argparse.ArgumentTypeError):
        nproc("999999")


def test_nproc_valid():
    """Test that a valid process count is accepted."""
    assert nproc("1") == 1


@pytest.mark.parametrize(
    "value, expected",
    [
        ("0", 0.0),
        ("1", 1.0),
        ("0.0", 0.0),
        ("1.0", 1.0),
        ("0.5", 0.5),
        ("0.001", 0.001),
        ("0.999", 0.999),
    ],
)
def test_bounded_float_01_accepts_valid_values(value, expected):
    """Test that valid values between 0 and 1 are accepted."""
    assert bounded_float_01(value) == expected


@pytest.mark.parametrize(
    "value",
    [
        "-0.001",
        "-1",
        "1.001",
        "2",
    ],
)
def test_bounded_float_01_rejects_values_outside_range(value):
    """Test that values outside the [0, 1] range are rejected."""
    with pytest.raises(
        argparse.ArgumentTypeError,
        match="value must be a finite float between 0 and 1",
    ):
        bounded_float_01(value)


@pytest.mark.parametrize(
    "value",
    [
        "foo",
        "abc",
        "",
        "not_a_number",
    ],
)
def test_bounded_float_01_rejects_invalid_floats(value):
    """Test that invalid float strings are rejected."""
    with pytest.raises(
        argparse.ArgumentTypeError,
        match=f"{value} is not a valid float",
    ):
        bounded_float_01(value)


@pytest.mark.parametrize(
    "value",
    [
        "nan",
        "NaN",
        "inf",
        "Inf",
        "+inf",
        "-inf",
    ],
)
def test_bounded_float_01_rejects_non_finite_values(value):
    """Test that non-finite float values are rejected."""
    with pytest.raises(
        argparse.ArgumentTypeError,
        match="value must be a finite float between 0 and 1",
    ):
        bounded_float_01(value)


# ---------------------------------------------------------------------------
# Command-line argument parsing
# ---------------------------------------------------------------------------


def test_search_conta_arguments(tmp_path, monkeypatch):
    """Test parsing arguments for the search_conta command."""
    species_file = tmp_path / "species.tsv"
    species_file.write_text("test\n")

    conta_file = tmp_path / "contamination.tsv"

    monkeypatch.setattr(
        "sys.argv",
        [
            "crocodeel",
            "search_conta",
            "-s",
            str(species_file),
            "-c",
            str(conta_file),
        ],
    )

    args = get_arguments()

    assert args.command == "search_conta"
    assert args.species_ab_table_fp == species_file.resolve()
    assert args.conta_events_fp == conta_file.resolve()


def test_plot_conta_arguments(tmp_path, monkeypatch):
    """Test parsing arguments for the plot_conta command."""
    species_file = tmp_path / "species.tsv"
    species_file.write_text("test\n")

    conta_file = tmp_path / "contamination.tsv"
    conta_file.write_text("test\n")

    pdf_file = tmp_path / "report.pdf"

    monkeypatch.setattr(
        "sys.argv",
        [
            "crocodeel",
            "plot_conta",
            "-s",
            str(species_file),
            "-c",
            str(conta_file),
            "-r",
            str(pdf_file),
        ],
    )

    args = get_arguments()

    assert args.command == "plot_conta"
    assert args.species_ab_table_fp == species_file.resolve()
    assert args.conta_events_fp == conta_file.resolve()
    assert args.pdf_report_fp == pdf_file.resolve()


def test_easy_wf_arguments(tmp_path, monkeypatch):
    """Test parsing arguments for the easy workflow."""
    species_file = tmp_path / "species.tsv"
    species_file.write_text("test\n")

    conta_file = tmp_path / "contamination.tsv"
    pdf_file = tmp_path / "report.pdf"

    monkeypatch.setattr(
        "sys.argv",
        [
            "crocodeel",
            "easy_wf",
            "-s",
            str(species_file),
            "-c",
            str(conta_file),
            "-r",
            str(pdf_file),
        ],
    )

    args = get_arguments()

    assert args.command == "easy_wf"
    assert args.species_ab_table_fp == species_file.resolve()
    assert args.conta_events_fp == conta_file.resolve()
    assert args.pdf_report_fp == pdf_file.resolve()


def test_easy_wf_plot_arguments(tmp_path, monkeypatch):
    """Test that easy_wf accepts plotting options."""
    species_file = tmp_path / "species.tsv"
    species_file.write_text("test\n")

    conta_file = tmp_path / "contamination.tsv"
    pdf_file = tmp_path / "report.pdf"

    monkeypatch.setattr(
        "sys.argv",
        [
            "crocodeel",
            "easy_wf",
            "-s",
            str(species_file),
            "-c",
            str(conta_file),
            "-r",
            str(pdf_file),
            "--nrow",
            "2",
            "--ncol",
            "3",
            "--no-conta-line",
            "--color-conta-species",
        ],
    )

    args = get_arguments()

    assert args.nrow == 2
    assert args.ncol == 3
    assert args.no_conta_line is True
    assert args.color_conta_species is True


def test_train_model_does_not_accept_second_abundance_table(
    tmp_path,
    monkeypatch,
):
    """Test that train_model does not accept a second abundance table."""
    species_file = tmp_path / "species.tsv"
    species_file.write_text("test\n")

    species_file_2 = tmp_path / "species2.tsv"
    species_file_2.write_text("test\n")

    model_file = tmp_path / "model.joblib"
    report_file = tmp_path / "report.json"

    monkeypatch.setattr(
        "sys.argv",
        [
            "crocodeel",
            "train_model",
            "-s",
            str(species_file),
            "-s2",
            str(species_file_2),
            "-m",
            str(model_file),
            "-r",
            str(report_file),
        ],
    )

    with pytest.raises(SystemExit):
        get_arguments()


def test_search_conta_requires_contamination_events_file(
    tmp_path,
    monkeypatch,
):
    """Test that search_conta requires an output contamination-events file."""
    species_file = tmp_path / "species.tsv"
    species_file.write_text("test\n")

    monkeypatch.setattr(
        "sys.argv",
        [
            "crocodeel",
            "search_conta",
            "-s",
            str(species_file),
        ],
    )

    with pytest.raises(SystemExit):
        get_arguments()


def test_plot_conta_requires_contamination_events_file(
    tmp_path,
    monkeypatch,
):
    """Test that plot_conta requires an input contamination-events file."""
    species_file = tmp_path / "species.tsv"
    species_file.write_text("test\n")

    pdf_file = tmp_path / "report.pdf"

    monkeypatch.setattr(
        "sys.argv",
        [
            "crocodeel",
            "plot_conta",
            "-s",
            str(species_file),
            "-r",
            str(pdf_file),
        ],
    )

    with pytest.raises(SystemExit):
        get_arguments()


# ---------------------------------------------------------------------------
# Workflow orchestration
# ---------------------------------------------------------------------------


def test_run_easy_workflow_loads_abundance_tables_once(monkeypatch):
    """Test that easy_wf reuses the loaded abundance tables."""
    species_ab_table = object()
    species_ab_table_2 = object()

    load_calls = []
    search_calls = []
    plot_calls = []

    def fake_load_abundance_tables(args):
        load_calls.append(args)
        return species_ab_table, species_ab_table_2

    def fake_run_search(args, table, table_2):
        search_calls.append((table, table_2))
        return ["event"]

    def fake_generate_pdf(args, table, table_2, events):
        plot_calls.append((table, table_2, events))

    monkeypatch.setattr(
        "crocodeel.crocodeel.load_abundance_tables",
        fake_load_abundance_tables,
    )
    monkeypatch.setattr(
        "crocodeel.crocodeel.run_search",
        fake_run_search,
    )
    monkeypatch.setattr(
        "crocodeel.crocodeel.generate_pdf_report",
        fake_generate_pdf,
    )
    monkeypatch.setattr(
        "crocodeel.crocodeel.log_contamination_warnings",
        lambda events: None,
    )

    args = object()

    run_easy_workflow(args)

    assert load_calls == [args]
    assert search_calls == [
        (species_ab_table, species_ab_table_2),
    ]
    assert plot_calls == [
        (species_ab_table, species_ab_table_2, ["event"]),
    ]
