"""Sampling helpers for reproducible benchmark runs."""

import numpy as np
from sage.all import GF, vector


def numpy_rng(seed):
    """Return a NumPy default RNG with a Python-native integer seed.

    Args:
        seed: Random seed for reproducible sampling.

    Returns:
        int: A NumPy default RNG with a Python-native integer seed.
    """
    # Keep syndrome normalization separate from the matching solve.
    return np.random.default_rng(int(seed))


def sample_data_error(num_qubits, p, rng):
    """Sample a reproducible independent data error vector.

    Args:
        num_qubits: Number of data-qubit columns per coarse lattice cell.
        p: Physical error probability.
        rng: Random number generator used for reproducible sampling.

    Returns:
        vector: Sampled GF(2) error vector.

    """
    # Draw reproducible random data before returning the sample.
    return vector(GF(2), [int(value) for value in rng.random(int(num_qubits)) < float(p)])
