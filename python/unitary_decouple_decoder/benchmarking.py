"""Benchmark helpers for unitary-decouple decoder experiments."""

from decoder_core.finite_maps import finite_matrix_from_laurent
from isomorphism.css import antipode_matrix, split_css_blocks


def build_hz_dagger_matrix(decoder):
    """Return finite source ``H_Z^dagger`` for a decoupling decoder.

    Args:
        decoder: Decoupling decoder with source CSS metadata.

    Returns:
        Matrix: Finite source ``H_Z^dagger`` matrix.
    """
    # The benchmark classifier needs the stabilizer column space in the same
    # finite coordinates as decoder corrections.
    hz_input, _hx_input = split_css_blocks(
        decoder.source_epsilon,
        decoder.num_x_checks,
        decoder.num_qubits,
    )
    return finite_matrix_from_laurent(antipode_matrix(hz_input), decoder.shape)
