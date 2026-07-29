"""Array conversion helpers for benchmark adapters."""

import numpy as np


def to_numpy_uint8(values):
    """Convert a Sage matrix/vector-like object to a NumPy ``uint8`` array.

    Args:
        values: Input sequence to normalize or validate.

    Returns:
        vector: Convert a Sage matrix/vector-like object to a NumPy ``uint8`` array.
    """
    # Keep syndrome normalization separate from the matching solve.
    if isinstance(values, np.ndarray) and values.dtype == np.uint8:
        return values
    if hasattr(values, "nrows") and hasattr(values, "ncols"):
        return np.asarray(
            [
                [int(values[row, col]) for col in range(values.ncols())]
                for row in range(values.nrows())
            ],
            dtype=np.uint8,
        )
    return np.asarray([int(value) for value in values], dtype=np.uint8)
