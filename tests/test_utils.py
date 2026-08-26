"""Unit tests for general CroCoDeEL utility functions."""

from crocodeel.utils import format_sample_names


def test_format_sample_names_keeps_all_names() -> None:
    """Test that short lists are formatted without truncation."""
    sample_names = ["sample1", "sample2", "sample3"]

    result = format_sample_names(sample_names)

    assert result == "sample1, sample2, sample3"


def test_format_sample_names_truncates_long_lists() -> None:
    """Test that long lists are truncated and report their total size."""
    sample_names = [
        "sample1",
        "sample2",
        "sample3",
        "sample4",
        "sample5",
        "sample6",
    ]

    result = format_sample_names(sample_names)

    assert result == (
        "6 sample names (showing 5): " "sample1, sample2, sample3, sample4, sample5"
    )


def test_format_sample_names_respects_max_names() -> None:
    """Test that the maximum number of displayed names can be configured."""
    sample_names = ["sample1", "sample2", "sample3", "sample4"]

    result = format_sample_names(sample_names, max_names=2)

    assert result == "4 sample names (showing 2): sample1, sample2"
