"""General utility functions for CroCoDeEL."""

from collections.abc import Sequence


def format_sample_names(
    sample_names: Sequence[str],
    max_names: int = 5,
) -> str:
    """Format sample names, truncating long lists."""

    if len(sample_names) <= max_names:
        return ", ".join(sample_names)

    return f"{len(sample_names)} sample names " f"(showing {max_names}): " + ", ".join(
        sample_names[:max_names]
    )
