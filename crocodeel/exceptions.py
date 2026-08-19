"""Exceptions raised by CroCoDeEL."""


class CrocodeelException(Exception):
    """Base exception for CroCoDeEL errors."""

class InputDataError(CrocodeelException):
    """Raised when input data does not meet CroCoDeEL requirements."""

class SelfTestError(CrocodeelException):
    """Raised when the CroCoDeEL self-test fails."""
