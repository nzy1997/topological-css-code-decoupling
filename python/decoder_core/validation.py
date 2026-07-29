"""Shared validation and sampling helpers."""

import numpy as np
from sage.all import GF, vector


def sample_bsc_error(n, p):
    """Sample an independent binary symmetric channel error vector.

    Args:
        n: Number of Bernoulli bits or qubits to sample.
        p: Physical error probability.

    Returns:
        vector: Sampled GF(2) error vector.
    """
    samples = np.random.binomial(1, float(p), size=int(n))
    return vector(GF(2), [int(value) for value in samples])


def matrix_to_numpy_uint8(values):
    """Convert a Sage matrix or vector-like object to a NumPy uint8 array.

    Args:
        values: Sage matrix or vector-like object.

    Returns:
        numpy.ndarray: NumPy array with ``uint8`` dtype.
    """
    if hasattr(values, "nrows") and hasattr(values, "ncols"):
        return np.asarray(
            [
                [int(values[row, col]) for col in range(values.ncols())]
                for row in range(values.nrows())
            ],
            dtype=np.uint8,
        )
    return np.asarray([int(value) for value in values], dtype=np.uint8)
