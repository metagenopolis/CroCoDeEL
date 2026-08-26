"""Unit tests for the CroCoDeEL command-line interface."""

import argparse
import logging
from unittest.mock import MagicMock, patch

import pandas as pd
import pytest

from crocodeel.conta_event import ContaminationEvent
from crocodeel.crocodeel import (bounded_float_01, bounded_int,
                                 generate_pdf_report, get_arguments,
                                 get_available_cpu_count,
                                 load_abundance_tables,
                                 load_contamination_events,
                                 log_contamination_warnings, main,
                                 positive_float, positive_int, readable_file,
                                 run_easy_workflow, run_plot_conta_command,
                                 run_search, run_search_conta_command,
                                 run_train_model_command, set_logging,
                                 writable_file)
from crocodeel.exceptions import InputDataError, SelfTestError

# ---------------------------------------------------------------------------
# main()
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "command, function_name",
    [
        ("search_conta", "run_search_conta_command"),
        ("plot_conta", "run_plot_conta_command"),
        ("easy_wf", "run_easy_workflow"),
        ("train_model", "run_train_model_command"),
    ],
)
def test_main_dispatches_command(
    monkeypatch,
    command,
    function_name,
):
    """Test that each CLI command is dispatched with the parsed arguments."""
    monkeypatch.setattr(
        "crocodeel.crocodeel.set_logging",
        lambda: None,
    )

    args = argparse.Namespace(command=command)

    monkeypatch.setattr(
        "crocodeel.crocodeel.get_arguments",
        lambda: args,
    )

    received_args = None

    def fake_command(command_args):
        nonlocal received_args
        received_args = command_args

    monkeypatch.setattr(
        f"crocodeel.crocodeel.{function_name}",
        fake_command,
    )

    assert main() == 0
    assert received_args is args


def test_main_returns_zero_on_success(monkeypatch):
    """Test that main returns zero when the command succeeds."""
    monkeypatch.setattr(
        "crocodeel.crocodeel.set_logging",
        lambda: None,
    )

    monkeypatch.setattr(
        "crocodeel.crocodeel.get_arguments",
        lambda: argparse.Namespace(
            command="self_test",
            keep_results=False,
        ),
    )

    monkeypatch.setattr(
        "crocodeel.crocodeel.SelfTest.run",
        lambda self: None,
    )

    assert main() == 0


@pytest.mark.parametrize(
    "exception",
    [
        InputDataError("Invalid input data."),
        SelfTestError("Self-test failed."),
    ],
)
def test_main_returns_one_on_handled_error(
    monkeypatch,
    caplog,
    exception,
):
    """Test that handled CLI errors result in exit code 1."""
    monkeypatch.setattr(
        "crocodeel.crocodeel.set_logging",
        lambda: None,
    )

    monkeypatch.setattr(
        "crocodeel.crocodeel.get_arguments",
        lambda: argparse.Namespace(
            command="self_test",
            keep_results=False,
        ),
    )

    def raise_error(self):
        raise exception

    monkeypatch.setattr(
        "crocodeel.crocodeel.SelfTest.run",
        raise_error,
    )

    with caplog.at_level(logging.ERROR):
        exit_code = main()

    assert exit_code == 1
    assert str(exception) in caplog.text


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


def test_readable_file_rejects_unreadable_file(tmp_path):
    """Test that an unreadable file is rejected."""
    fp = tmp_path / "unreadable.txt"
    fp.touch()
    fp.chmod(0o000)

    try:
        with pytest.raises(
            argparse.ArgumentTypeError,
            match="is not readable",
        ):
            readable_file(str(fp))
    finally:
        # Restore permissions so pytest can clean up tmp_path.
        fp.chmod(0o644)


def test_readable_file_rejects_directory(tmp_path):
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


def test_writable_file_rejects_directory(tmp_path):
    """Test that a directory is rejected as an output file."""
    with pytest.raises(
        argparse.ArgumentTypeError,
        match="is a directory",
    ):
        writable_file(str(tmp_path))


def test_writable_file_rejects_non_writable_file(tmp_path):
    """Test that a non-writable file is rejected."""
    fp = tmp_path / "output.txt"
    fp.touch()
    fp.chmod(0o444)

    try:
        with pytest.raises(
            argparse.ArgumentTypeError,
            match="is not writable",
        ):
            writable_file(str(fp))
    finally:
        # Restore permissions so pytest can clean up tmp_path.
        fp.chmod(0o644)


def test_writable_file_rejects_missing_parent(tmp_path):
    """Test that a missing parent directory is rejected."""
    filepath = tmp_path / "missing" / "output.tsv"

    with pytest.raises(
        argparse.ArgumentTypeError,
        match="does not exist",
    ):
        writable_file(str(filepath))


def test_writable_file_rejects_non_directory_parent(tmp_path):
    """Test that a non-directory parent path is rejected."""
    parent_path = tmp_path / "parent"
    parent_path.touch()

    output_path = parent_path / "output.txt"

    with pytest.raises(
        argparse.ArgumentTypeError,
        match="is not a directory",
    ):
        writable_file(str(output_path))


def test_writable_file_rejects_non_writable_parent_directory(tmp_path):
    """Test that a non-writable parent directory is rejected."""
    parent_dir = tmp_path / "output"
    parent_dir.mkdir()
    parent_dir.chmod(0o555)

    try:
        output_path = parent_dir / "output.txt"

        with pytest.raises(
            argparse.ArgumentTypeError,
            match="directory .* is not writable",
        ):
            writable_file(str(output_path))
    finally:
        # Restore permissions so pytest can clean up tmp_path.
        parent_dir.chmod(0o755)


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

@pytest.mark.parametrize(
    "value, minimum, maximum, expected",
    [
        ("1", 1, 4, 1),
        ("2", 1, 4, 2),
        ("4", 1, 4, 4),
    ],
)
def test_bounded_int_valid(value, minimum, maximum, expected):
    """Test that bounded_int accepts values within the range."""
    assert bounded_int(value, minimum, maximum) == expected


@pytest.mark.parametrize(
    "value",
    [
        "abc",
        "1.5",
        "0",
        "-1",
        "5",
    ],
)
def test_bounded_int_invalid(value):
    """Test that bounded_int rejects invalid and out-of-range values."""
    with pytest.raises(argparse.ArgumentTypeError):
        bounded_int(value, 1, 4)


def test_get_available_cpu_count_uses_cpu_affinity(monkeypatch):
    """Test that CPU affinity determines the available CPU count."""
    monkeypatch.setattr(
        "crocodeel.crocodeel.os.sched_getaffinity",
        lambda pid: {0, 1, 2, 3},
    )

    assert get_available_cpu_count() == 4


def test_get_available_cpu_count_falls_back_to_cpu_count(monkeypatch):
    """Test the CPU count fallback when CPU affinity is unavailable."""
    monkeypatch.delattr(
        "crocodeel.crocodeel.os.sched_getaffinity",
        raising=False,
    )
    monkeypatch.setattr(
        "crocodeel.crocodeel.multiprocessing.cpu_count",
        lambda: 8,
    )

    assert get_available_cpu_count() == 8


# ---------------------------------------------------------------------------
# Command-line argument parsing
# ---------------------------------------------------------------------------


def test_set_logging(monkeypatch):
    """Test that logging is configured with the expected format and level."""
    mock_basic_config = MagicMock()

    monkeypatch.setattr(
        "crocodeel.crocodeel.logging.basicConfig",
        mock_basic_config,
    )

    set_logging()

    mock_basic_config.assert_called_once_with(
        format="%(asctime)s :: %(levelname)s :: %(message)s",
        level=logging.INFO,
    )


@pytest.mark.parametrize(
    "command, command_args",
    [
        (
            "search_conta",
            [
                "-s",
                "species.tsv",
                "-s2",
                "species2.tsv",
                "-m",
                "model.joblib",
                "--filter-low-ab",
                "20",
                "--probability-cutoff",
                "0.8",
                "--rate-cutoff",
                "0.01",
                "--nproc",
                "1",
                "-c",
                "contamination.tsv",
            ],
        ),
        (
            "plot_conta",
            [
                "-s",
                "species.tsv",
                "-s2",
                "species2.tsv",
                "--filter-low-ab",
                "20",
                "-c",
                "contamination.tsv",
                "-r",
                "report.pdf",
                "--nrow",
                "5",
                "--ncol",
                "6",
                "--no-conta-line",
                "--color-conta-species",
            ],
        ),
        (
            "easy_wf",
            [
                "-s",
                "species.tsv",
                "-s2",
                "species2.tsv",
                "-m",
                "model.joblib",
                "--filter-low-ab",
                "20",
                "-c",
                "contamination.tsv",
                "--probability-cutoff",
                "0.8",
                "--rate-cutoff",
                "0.01",
                "--nproc",
                "1",
                "-r",
                "report.pdf",
                "--nrow",
                "5",
                "--ncol",
                "6",
                "--no-conta-line",
                "--color-conta-species",
            ],
        ),
        (
            "train_model",
            [
                "-s",
                "species.tsv",
                "--filter-low-ab",
                "20",
                "-m",
                "model.joblib",
                "-r",
                "report.json",
                "--test-size",
                "0.25",
                "--ntrees",
                "100",
                "--rng-seed",
                "42",
                "--nproc",
                "1",
            ],
        ),
    ],
)
def test_get_arguments_parses_valid_command(
    tmp_path,
    monkeypatch,
    command,
    command_args,
):
    """Test that each file-based subcommand parses a complete command line."""
    for filename in {
        "species.tsv",
        "species2.tsv",
        "model.joblib",
        "contamination.tsv",
        "report.pdf",
        "report.json",
    }:
        (tmp_path / filename).touch()

    args = [
        "crocodeel",
        command,
        *[
            str(tmp_path / arg)
            if arg in {
                "species.tsv",
                "species2.tsv",
                "model.joblib",
                "contamination.tsv",
                "report.pdf",
                "report.json",
            }
            else arg
            for arg in command_args
        ],
    ]

    monkeypatch.setattr("sys.argv", args)

    parsed_args = get_arguments()

    assert parsed_args.command == command


@pytest.mark.parametrize("command", ["search_conta", "train_model"])
def test_nproc_uses_available_cpu_count_as_upper_bound(
    tmp_path,
    monkeypatch,
    command,
):
    """Test that --nproc uses the available CPU count as its upper bound."""
    available_cpu_count = 4

    monkeypatch.setattr(
        "crocodeel.crocodeel.get_available_cpu_count",
        lambda: available_cpu_count,
    )

    species_file = tmp_path / "species.tsv"
    species_file.write_text("test\n")

    args = [
        "crocodeel",
        command,
        "-s",
        str(species_file),
        "--nproc",
        "4",
    ]

    if command == "search_conta":
        args.extend(["-c", str(tmp_path / "contamination.tsv")])
    else:
        args.extend(
            [
                "-r",
                str(tmp_path / "report.json"),
                "-m",
                str(tmp_path / "model.joblib"),
            ]
        )

    monkeypatch.setattr("sys.argv", args)

    parsed_args = get_arguments()

    assert parsed_args.nproc == 4


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
# Abundance table loading
# ---------------------------------------------------------------------------


def test_load_abundance_tables_without_second_table(
    monkeypatch,
    tmp_path,
) -> None:
    """Test loading a single species abundance table."""
    species_ab_table_fp = tmp_path / "species_abundance.tsv"
    species_ab_table_fp.touch()

    species_ab_table = pd.DataFrame(
        {"sample1": [0.1, 0.2]}
    )

    mock_read = MagicMock(return_value=species_ab_table)

    monkeypatch.setattr(
        "crocodeel.crocodeel.ab_table_utils.read_filter_normalize",
        mock_read,
    )

    args = argparse.Namespace(
        species_ab_table_fp=species_ab_table_fp,
        species_ab_table_2_fp=None,
        filtering_ab_thr_factor=20,
    )

    result, result_2 = load_abundance_tables(args)

    assert result is species_ab_table
    assert result_2 is None

    mock_read.assert_called_once()
    assert mock_read.call_args.args[0].name == str(species_ab_table_fp)
    assert mock_read.call_args.args[1] == 20


def test_load_abundance_tables_with_second_table(
    monkeypatch,
    tmp_path,
) -> None:
    """Test loading and comparing two species abundance tables."""
    species_ab_table_fp = tmp_path / "species_abundance.tsv"
    species_ab_table_2_fp = tmp_path / "species_abundance_2.tsv"

    species_ab_table_fp.touch()
    species_ab_table_2_fp.touch()

    species_ab_table = pd.DataFrame(
        {"sample1": [0.1, 0.2]}
    )
    species_ab_table_2 = pd.DataFrame(
        {"sample2": [0.3, 0.4]}
    )

    mock_read = MagicMock(
        side_effect=[
            species_ab_table,
            species_ab_table_2,
        ]
    )
    mock_compare = MagicMock()

    monkeypatch.setattr(
        "crocodeel.crocodeel.ab_table_utils.read_filter_normalize",
        mock_read,
    )
    monkeypatch.setattr(
        "crocodeel.crocodeel.ab_table_utils.compare_species_names",
        mock_compare,
    )

    args = argparse.Namespace(
        species_ab_table_fp=species_ab_table_fp,
        species_ab_table_2_fp=species_ab_table_2_fp,
        filtering_ab_thr_factor=20,
    )

    result, result_2 = load_abundance_tables(args)

    assert result is species_ab_table
    assert result_2 is species_ab_table_2

    assert mock_read.call_count == 2

    assert mock_read.call_args_list[0].args[0].name == str(
        species_ab_table_fp
    )
    assert mock_read.call_args_list[1].args[0].name == str(
        species_ab_table_2_fp
    )

    assert mock_read.call_args_list[0].args[1] == 20
    assert mock_read.call_args_list[1].args[1] == 20

    mock_compare.assert_called_once_with(
        species_ab_table,
        species_ab_table_2,
    )


def test_load_abundance_tables_same_file(
    monkeypatch,
    tmp_path,
) -> None:
    """Test that the second table is not loaded when both paths are identical."""
    species_ab_table_fp = tmp_path / "species_abundance.tsv"
    species_ab_table_fp.touch()

    species_ab_table = pd.DataFrame(
        {"sample1": [0.1, 0.2]}
    )

    mock_read = MagicMock(return_value=species_ab_table)
    mock_compare = MagicMock()

    monkeypatch.setattr(
        "crocodeel.crocodeel.ab_table_utils.read_filter_normalize",
        mock_read,
    )
    monkeypatch.setattr(
        "crocodeel.crocodeel.ab_table_utils.compare_species_names",
        mock_compare,
    )

    args = argparse.Namespace(
        species_ab_table_fp=species_ab_table_fp,
        species_ab_table_2_fp=species_ab_table_fp,
        filtering_ab_thr_factor=20,
    )

    result, result_2 = load_abundance_tables(args)

    assert result is species_ab_table
    assert result_2 is None

    mock_read.assert_called_once()
    mock_compare.assert_not_called()


# ---------------------------------------------------------------------------
# Command orchestration
# ---------------------------------------------------------------------------


def test_run_search_conta_command(
    monkeypatch,
) -> None:
    """Test orchestration of the search_conta command."""
    args = argparse.Namespace()

    species_ab_table = pd.DataFrame(
        {"sample1": [0.1, 0.2]}
    )
    species_ab_table_2 = pd.DataFrame(
        {"sample2": [0.3, 0.4]}
    )

    conta_events = [
        ContaminationEvent(
            source="sample1",
            target="sample2",
            rate=0.1,
            probability=0.9,
            conta_line_species=["species1"],
        )
    ]

    mock_load_abundance_tables = MagicMock(
        return_value=(species_ab_table, species_ab_table_2)
    )
    mock_run_search = MagicMock(return_value=conta_events)
    mock_log_warnings = MagicMock()

    monkeypatch.setattr(
        "crocodeel.crocodeel.load_abundance_tables",
        mock_load_abundance_tables,
    )
    monkeypatch.setattr(
        "crocodeel.crocodeel.run_search",
        mock_run_search,
    )
    monkeypatch.setattr(
        "crocodeel.crocodeel.log_contamination_warnings",
        mock_log_warnings,
    )

    run_search_conta_command(args)

    mock_load_abundance_tables.assert_called_once_with(args)

    mock_run_search.assert_called_once_with(
        args,
        species_ab_table,
        species_ab_table_2,
    )

    mock_log_warnings.assert_called_once_with(conta_events)


def test_run_plot_conta_command(
    monkeypatch,
) -> None:
    """Test orchestration of the plot_conta command."""
    args = argparse.Namespace()

    species_ab_table = pd.DataFrame(
        {"sample1": [0.1, 0.2]}
    )
    species_ab_table_2 = pd.DataFrame(
        {"sample2": [0.3, 0.4]}
    )

    conta_events = [
        ContaminationEvent(
            source="sample1",
            target="sample2",
            rate=0.1,
            probability=0.9,
            conta_line_species=["species1"],
        )
    ]

    mock_load_abundance_tables = MagicMock(
        return_value=(species_ab_table, species_ab_table_2)
    )
    mock_load_contamination_events = MagicMock(
        return_value=conta_events
    )
    mock_generate_pdf_report = MagicMock()

    monkeypatch.setattr(
        "crocodeel.crocodeel.load_abundance_tables",
        mock_load_abundance_tables,
    )
    monkeypatch.setattr(
        "crocodeel.crocodeel.load_contamination_events",
        mock_load_contamination_events,
    )
    monkeypatch.setattr(
        "crocodeel.crocodeel.generate_pdf_report",
        mock_generate_pdf_report,
    )

    run_plot_conta_command(args)

    mock_load_abundance_tables.assert_called_once_with(args)
    mock_load_contamination_events.assert_called_once_with(args)

    mock_generate_pdf_report.assert_called_once_with(
        args,
        species_ab_table,
        species_ab_table_2,
        conta_events,
    )


def test_run_train_model_command(
    monkeypatch,
    tmp_path,
) -> None:
    """Test orchestration of the train_model command."""
    species_ab_table_fp = tmp_path / "species_abundance.tsv"
    model_fp = tmp_path / "model.joblib"
    json_report_fp = tmp_path / "report.json"

    species_ab_table_fp.write_text(
        "dummy input",
        encoding="utf8",
    )

    args = argparse.Namespace(
        species_ab_table_fp=species_ab_table_fp,
        model_fp=model_fp,
        json_report_fp=json_report_fp,
        filtering_ab_thr_factor=20,
        test_size=0.2,
        ntrees=100,
        rng_seed=42,
        nproc=4,
    )

    species_ab_table = pd.DataFrame(
        {"sample1": [0.1, 0.2]}
    )

    mock_read_filter_normalize = MagicMock(
        return_value=species_ab_table,
    )
    mock_run_train_model = MagicMock()

    monkeypatch.setattr(
        "crocodeel.crocodeel.ab_table_utils.read_filter_normalize",
        mock_read_filter_normalize,
    )
    monkeypatch.setattr(
        "crocodeel.crocodeel.run_train_model",
        mock_run_train_model,
    )

    run_train_model_command(args)

    mock_read_filter_normalize.assert_called_once()

    read_args = mock_read_filter_normalize.call_args.args

    assert read_args[0].name == str(species_ab_table_fp)
    assert read_args[1] == 20

    mock_run_train_model.assert_called_once()

    train_args = mock_run_train_model.call_args.args

    assert train_args[0] is species_ab_table
    assert train_args[1].name == str(model_fp)
    assert train_args[2].name == str(json_report_fp)
    assert train_args[3:] == (
        0.2,
        100,
        42,
        4,
    )

    # The output files are closed by the context managers.
    assert train_args[1].closed
    assert train_args[2].closed


def test_run_search(
    monkeypatch,
    tmp_path,
) -> None:
    """Test orchestration of contamination search and TSV output."""
    species_ab_table = pd.DataFrame(
        {"sample1": [0.1, 0.2]}
    )
    species_ab_table_2 = pd.DataFrame(
        {"sample2": [0.3, 0.4]}
    )

    rf_model_fp = tmp_path / "model.joblib"
    conta_events_fp = tmp_path / "contamination_events.tsv"

    rf_model_fp.write_bytes(b"fake model")

    args = argparse.Namespace(
        species_ab_table_fp=tmp_path / "species_abundance.tsv",
        species_ab_table_2_fp=tmp_path / "species_abundance_2.tsv",
        rf_model_fp=rf_model_fp,
        conta_events_fp=conta_events_fp,
        filtering_ab_thr_factor=20,
        probability_cutoff=0.5,
        rate_cutoff=0.01,
        nproc=2,
    )

    conta_events = [
        ContaminationEvent(
            source="sample1",
            target="sample2",
            rate=0.1,
            probability=0.9,
            conta_line_species=["species1"],
        )
    ]

    mock_run_search_conta = MagicMock(
        return_value=conta_events,
    )
    mock_write_tsv = MagicMock()
    mock_exec_desc = MagicMock(
        return_value="execution description",
    )

    monkeypatch.setattr(
        "crocodeel.crocodeel.run_search_conta",
        mock_run_search_conta,
    )
    monkeypatch.setattr(
        "crocodeel.crocodeel.ContaminationEventIO.write_tsv",
        mock_write_tsv,
    )
    monkeypatch.setattr(
        "crocodeel.crocodeel.ExecutionDescription",
        mock_exec_desc,
    )

    result = run_search(
        args,
        species_ab_table,
        species_ab_table_2,
    )

    assert result == conta_events

    mock_exec_desc.assert_called_once_with(
        args.species_ab_table_fp,
        args.species_ab_table_2_fp,
        args.rf_model_fp,
        args.filtering_ab_thr_factor,
        args.probability_cutoff,
        args.rate_cutoff,
    )

    mock_run_search_conta.assert_called_once()

    search_args = mock_run_search_conta.call_args.args

    assert search_args[0] is species_ab_table
    assert search_args[1] is species_ab_table_2
    assert search_args[3:] == (
        0.5,
        0.01,
        2,
    )

    mock_write_tsv.assert_called_once_with(
        conta_events,
        mock_write_tsv.call_args.args[1],
    )

    with conta_events_fp.open("r", encoding="utf8") as fh:
        assert fh.readline() == "execution description\n"


def test_generate_pdf_report(
    monkeypatch,
    tmp_path,
) -> None:
    """Test orchestration of PDF report generation."""
    pdf_report_fp = tmp_path / "contamination_report.pdf"

    args = argparse.Namespace(
        pdf_report_fp=pdf_report_fp,
        nrow=2,
        ncol=3,
        no_conta_line=False,
        color_conta_species=False,
    )

    species_ab_table = pd.DataFrame(
        {"sample1": [0.1, 0.2]}
    )
    species_ab_table_2 = pd.DataFrame(
        {"sample2": [0.3, 0.4]}
    )

    conta_events = [
        ContaminationEvent(
            source="sample1",
            target="sample2",
            rate=0.1,
            probability=0.9,
            conta_line_species=["species1"],
        )
    ]

    mock_run_plot_conta = MagicMock()

    monkeypatch.setattr(
        "crocodeel.crocodeel.run_plot_conta",
        mock_run_plot_conta,
    )

    generate_pdf_report(
        args,
        species_ab_table,
        species_ab_table_2,
        conta_events,
    )

    mock_run_plot_conta.assert_called_once()

    plot_args = mock_run_plot_conta.call_args.args

    assert plot_args[0] is species_ab_table
    assert plot_args[1] is species_ab_table_2
    assert plot_args[2] is conta_events
    assert plot_args[3].name == str(pdf_report_fp)
    assert plot_args[4:] == (
        2,
        3,
        False,
        False,
    )

    # The output file is closed by the context manager.
    assert plot_args[3].closed


def test_generate_pdf_report_without_events(
    tmp_path,
    caplog,
) -> None:
    """Test that no PDF is generated when no contamination is detected."""
    pdf_report_fp = tmp_path / "report.pdf"
    args = argparse.Namespace(
        pdf_report_fp=pdf_report_fp,
    )

    with (
        patch(
            "crocodeel.crocodeel.run_plot_conta",
        ) as mock_run_plot_conta,
        caplog.at_level(logging.INFO),
    ):
        generate_pdf_report(
            args,
            pd.DataFrame(),
            None,
            [],
        )

    mock_run_plot_conta.assert_not_called()
    assert not pdf_report_fp.exists()
    assert (
        "No PDF report generated."
        in caplog.text
    )


def test_load_contamination_events(
    monkeypatch,
    tmp_path,
) -> None:
    """Test loading contamination events from a TSV file."""
    conta_events_fp = tmp_path / "contamination_events.tsv"
    conta_events_fp.touch()

    conta_events = [
        ContaminationEvent(
            source="sample1",
            target="sample2",
            rate=0.1,
            probability=0.9,
            conta_line_species=["species1"],
        )
    ]

    mock_read_tsv = MagicMock(return_value=conta_events)

    monkeypatch.setattr(
        "crocodeel.crocodeel.ContaminationEventIO.read_tsv",
        mock_read_tsv,
    )

    args = argparse.Namespace(
        conta_events_fp=conta_events_fp,
    )

    result = load_contamination_events(args)

    assert result == conta_events

    mock_read_tsv.assert_called_once()

    file_handle = mock_read_tsv.call_args.args[0]

    assert file_handle.name == str(conta_events_fp)
    assert file_handle.mode == "r"
    assert file_handle.closed


# ---------------------------------------------------------------------------
# Easy workflow
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


# ---------------------------------------------------------------------------
# Contamination warnings
# ---------------------------------------------------------------------------


def test_log_contamination_warnings_does_nothing_for_empty_events(caplog):
    """Test that no warnings are logged when no contamination events are detected."""
    with caplog.at_level(logging.WARNING):
        log_contamination_warnings([])

    assert not caplog.records


def test_log_contamination_warnings_logs_warnings(caplog):
    """Test that warnings are logged when contamination events are detected."""
    conta_events = [
        ContaminationEvent(
            source="sample1",
            target="sample2",
            rate=0.01,
            probability=0.95,
        )
    ]

    with caplog.at_level(logging.WARNING):
        log_contamination_warnings(conta_events)

    assert len(caplog.records) == 4

    messages = [record.message for record in caplog.records]

    assert "decision-support tool" in messages[0]
    assert "false positives" in messages[1]
    assert "manually reviewing the scatterplots" in messages[2]
    assert "CroCoDeEL Interpretation Interface" in messages[3]
