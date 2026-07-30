"""Benchmark helpers for unitary-decouple-based decoder experiments."""

from decoder_core.finite_maps import finite_matrix_from_laurent
from decoupling.css import dagger_matrix, split_css_blocks


def build_h_z_dagger_finite(decoder):
    """Return finite source ``H_Z^dagger`` for a decoupling decoder.

    Args:
        decoder: Decoupling decoder with source CSS metadata.

    Returns:
        Matrix: Finite source ``H_Z^dagger`` matrix.
    """
    # The benchmark classifier needs the stabilizer column space in the same
    # finite coordinates as decoder corrections.
    h_z_input, _h_x_input = split_css_blocks(
        decoder.input_excitation_map,
        decoder.num_x_checks,
        decoder.num_qubits,
    )
    return finite_matrix_from_laurent(dagger_matrix(h_z_input), decoder.shape)
