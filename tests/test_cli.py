"""Unit tests for the CroCoDeEL command-line interface."""

import argparse
from pathlib import Path

import pytest

from crocodeel.crocodeel import (
    positive_float,
    nproc,
    readable_file,
    writable_file,
    bounded_float_01,
)


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