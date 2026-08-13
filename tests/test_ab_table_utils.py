"""Unit tests for species abundance table utilities."""

import io
import logging

import pandas as pd
import numpy as np
import pytest

from crocodeel.ab_table_utils import (
    filter_low_ab,
    normalize,
    log_transform,
    read,
    read_filter_normalize,
)


def test_read():
    """Test reading a species abundance table."""
    table = io.StringIO(
        "species_name\tsample1\tsample2\n"
        "species_1\t0.1\t0.2\n"
        "species_2\t0.9\t0.8\n"
    )
    table.name = "species_abundance.tsv"

    result = read(table)

    expected = pd.DataFrame(
        {
            "sample1": [0.1, 0.9],
            "sample2": [0.2, 0.8],
        },
        index=pd.Index(
            ["species_1", "species_2"],
            name="species_name",
        ),
    )

    pd.testing.assert_frame_equal(result, expected)


def test_read_accepts_comments():
    """Test that comment lines are ignored."""
    table = io.StringIO(
        "# Comment\n"
        "species_name\tsample1\n"
        "species_1\t0.1\n"
        "species_2\t0.9\n"
    )
    table.name = "species_abundance.tsv"

    result = read(table)

    assert result.shape == (2, 1)


def test_read_logs_table_dimensions(caplog):
    """Test that table dimensions are logged."""
    table = io.StringIO(
        "species_name\tsample1\tsample2\n"
        "species_1\t0.1\t0.2\n"
        "species_2\t0.9\t0.8\n"
    )
    table.name = "species_abundance.tsv"

    with caplog.at_level(logging.INFO):
        read(table)

    assert "2 species in 2 samples" in caplog.text


def test_read_filter_normalize_converts_integer_species_names_to_strings():
    """Test that integer species names are converted to strings."""
    table = io.StringIO(
        "species_name\tsample1\n"
        "123\t0.1\n"
        "456\t0.9\n"
    )
    table.name = "species_abundance.tsv"

    result = read_filter_normalize(table)

    assert result.index.tolist() == ["123", "456"]
    assert result.index.dtype == object


def test_read_filter_normalize_rejects_invalid_species_names():
    """Test that non-string and non-integer species names are rejected."""
    table = io.StringIO(
        "species_name\tsample1\n"
        "0.5\t0.1\n"
        "0.6\t0.9\n"
    )
    table.name = "species_abundance.tsv"

    with pytest.raises(SystemExit):
        read_filter_normalize(table)


def test_read_filter_normalize_rejects_non_numeric_abundances():
    """Test that non-numeric abundances are rejected."""
    table = io.StringIO(
        "species_name\tsample1\tsample2\n"
        "species_1\t0.1\tfoo\n"
        "species_2\t0.9\tbar\n"
    )
    table.name = "species_abundance.tsv"

    with pytest.raises(SystemExit):
        read_filter_normalize(table)


def test_read_filter_normalize_rejects_negative_abundances():
    """Test that negative abundances are rejected."""
    table = io.StringIO(
        "species_name\tsample1\tsample2\n"
        "species_1\t-0.1\t0.2\n"
        "species_2\t1.1\t0.8\n"
    )
    table.name = "species_abundance.tsv"

    with pytest.raises(SystemExit):
        read_filter_normalize(table)


def test_read_filter_normalize_rejects_empty_samples():
    """Test that samples with no non-zero abundances are rejected."""
    table = io.StringIO(
        "species_name\tsample1\tsample2\n"
        "species_1\t0.0\t0.2\n"
        "species_2\t0.0\t0.8\n"
    )
    table.name = "species_abundance.tsv"

    with pytest.raises(SystemExit):
        read_filter_normalize(table)


def test_normalize():
    """Test normalization to relative abundances."""
    table = pd.DataFrame(
        {
            "sample1": [1.0, 3.0],
            "sample2": [2.0, 2.0],
        },
        index=["species_1", "species_2"],
    )

    result = normalize(table)

    expected = pd.DataFrame(
        {
            "sample1": [0.25, 0.75],
            "sample2": [0.50, 0.50],
        },
        index=["species_1", "species_2"],
    )

    pd.testing.assert_frame_equal(result, expected)


def test_normalize_columns_sum_to_one():
    """Test that normalized samples sum to one."""
    table = pd.DataFrame(
        {
            "sample1": [1.0, 3.0],
            "sample2": [2.0, 2.0],
        }
    )

    result = normalize(table)

    assert result.sum(axis=0)["sample1"] == pytest.approx(1.0)
    assert result.sum(axis=0)["sample2"] == pytest.approx(1.0)


def test_filter_low_ab():
    """Test filtering of low-abundance species."""
    table = pd.DataFrame(
        {
            "sample1": [0.01, 0.10, 0.50, 0.39],
        },
        index=["species_1", "species_2", "species_3", "species_4"],
    )

    result = filter_low_ab(table, filtering_ab_thr_factor=2)

    # Minimum non-zero abundance = 0.01.
    # Threshold = 2 * 0.01 = 0.02.
    # species_1 (0.01) is filtered.
    expected = pd.DataFrame(
        {
            "sample1": [0.0, 0.10, 0.50, 0.39],
        },
        index=["species_1", "species_2", "species_3", "species_4"],
    )

    pd.testing.assert_frame_equal(result, expected)


def test_filter_low_ab_removes_species_at_threshold():
    """Test that abundances equal to the threshold are removed."""
    table = pd.DataFrame(
        {
            "sample1": [0.01, 0.02, 0.03],
        },
        index=["species_1", "species_2", "species_3"],
    )

    result = filter_low_ab(table, filtering_ab_thr_factor=2)

    # Minimum non-zero abundance = 0.01.
    # Threshold = 0.02.
    # The implementation uses <=, so 0.01 and 0.02 are removed.
    assert result.loc["species_1", "sample1"] == 0.0
    assert result.loc["species_2", "sample1"] == 0.0
    assert result.loc["species_3", "sample1"] == 0.03

def test_log_transform():
    """Test log10 transformation of species abundances."""
    table = pd.DataFrame(
        {
            "sample1": [1.0, 0.1, 0.01],
            "sample2": [0.5, 0.25, 0.001],
        },
        index=["species_1", "species_2", "species_3"],
    )

    result = log_transform(table)

    expected = pd.DataFrame(
        {
            "sample1": [0.0, -1.0, -2.0],
            "sample2": [np.log10(0.5), np.log10(0.25), -3.0],
        },
        index=["species_1", "species_2", "species_3"],
    )

    pd.testing.assert_frame_equal(result, expected)


def test_log_transform_converts_zeros_to_negative_infinity():
    """Test that zero abundances are converted to -inf."""
    table = pd.DataFrame(
        {
            "sample1": [1.0, 0.0, 0.1],
            "sample2": [0.0, 0.5, 0.0],
        },
        index=["species_1", "species_2", "species_3"],
    )

    result = log_transform(table)

    assert result.loc["species_1", "sample1"] == 0.0
    assert result.loc["species_2", "sample1"] == -np.inf
    assert result.loc["species_3", "sample1"] == -1.0

    assert result.loc["species_1", "sample2"] == -np.inf
    assert result.loc["species_2", "sample2"] == np.log10(0.5)
    assert result.loc["species_3", "sample2"] == -np.inf


def test_read_filter_normalize_without_filtering():
    """Test reading, normalization, and log transformation without filtering."""
    table = io.StringIO(
        "species_name\tsample1\n"
        "species_1\t1\n"
        "species_2\t3\n"
    )
    table.name = "species_abundance.tsv"

    result = read_filter_normalize(table)

    expected = pd.DataFrame(
        {
            "sample1": [np.log10(0.25), np.log10(0.75)],
        },
        index=pd.Index(
            ["species_1", "species_2"],
            name="species_name",
        ),
    )

    pd.testing.assert_frame_equal(result, expected)

def test_read_filter_normalize_with_filtering():
    """Test filtering, normalization, and log transformation."""
    table = io.StringIO(
        "species_name\tsample1\n"
        "species_1\t1\n"
        "species_2\t3\n"
        "species_3\t10\n"
    )
    table.name = "species_abundance.tsv"

    result = read_filter_normalize(
        table,
        filtering_ab_thr_factor=2,
    )

    # Minimum abundance = 1.
    # Threshold = 2.
    # species_1 is removed; species_2 and species_3 are retained.
    # After filtering, the total abundance is 13.
    # The final output is log10-transformed relative abundance.
    expected = pd.DataFrame(
        {
            "sample1": [
                -np.inf,
                np.log10(3 / 13),
                np.log10(10 / 13),
            ],
        },
        index=pd.Index(
            ["species_1", "species_2", "species_3"],
            name="species_name",
        ),
    )

    pd.testing.assert_frame_equal(result, expected)

def test_read_filter_normalize_rejects_samples_emptied_by_filtering():
    """Test that filtering cannot remove all species from a sample."""
    table = io.StringIO(
        "species_name\tsample1\n"
        "species_1\t1\n"
        "species_2\t2\n"
    )
    table.name = "species_abundance.tsv"

    # Minimum abundance = 1.
    # Threshold = 20.
    # Both species are removed.
    with pytest.raises(SystemExit):
        read_filter_normalize(
            table,
            filtering_ab_thr_factor=20,
        )
