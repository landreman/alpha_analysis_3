"""Ridge-surface computation of the constant-J alpha-particle loss surrogate f.

The scientific and numerical contract is docs/DESIGN.md. Milestone progress is recorded
in docs/STATUS.md.
"""

from alpha_analysis_3.conventions import (
    ABSOLUTE_TOLERANCE,
    LAMBDA_LEVELS,
    MEMORY_TARGET_GB,
    RELATIVE_TOLERANCE,
    RUNTIME_TARGET_SECONDS,
    default_source,
    field_line_direction,
    lambda_from_pitch,
    numerical_error_target,
    pitch_from_lambda,
    source_at_s,
    uniform_source,
)

__all__ = [
    "ABSOLUTE_TOLERANCE",
    "LAMBDA_LEVELS",
    "MEMORY_TARGET_GB",
    "RELATIVE_TOLERANCE",
    "RUNTIME_TARGET_SECONDS",
    "default_source",
    "field_line_direction",
    "lambda_from_pitch",
    "numerical_error_target",
    "pitch_from_lambda",
    "source_at_s",
    "uniform_source",
]

__version__ = "0.1.0"
