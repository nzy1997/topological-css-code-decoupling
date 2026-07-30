"""Run the paper's two color-code examples through the full algorithm."""

import _bootstrap  # noqa: F401

from sage.all import identity_matrix

from decoupling import (
    coarse_grain_to_superlattice,
    dagger_matrix,
    solve_decoupling_unitary,
)
from decoupling.chain_maps.decoupling import verify_qca_decoupling
from decoupling.css import split_css_blocks
from decoupling.paper_examples import COLOR_CODE_488, COLOR_CODE_666


def verify_paper_example(case):
    """Verify one paper example from coarse-graining through QCA decoupling."""
    coarse = coarse_grain_to_superlattice(
        case.excitation_map(),
        *case.superlattice_basis,
    )
    num_x_checks = coarse.nrows() // 2
    num_qubits = coarse.ncols() // 2
    result = solve_decoupling_unitary(
        coarse,
        num_x_checks=num_x_checks,
        num_qubits=num_qubits,
        verify=True,
        compute_psi=True,
    )

    product_x_rank = int(result.diagnostics["product_x_rank"])
    product_z_rank = int(result.diagnostics["product_z_rank"])
    toric_x_rank = result.h_x_tilde.nrows() - product_x_rank
    toric_z_rank = result.h_z_tilde_dagger.ncols() - product_z_rank
    assert toric_x_rank == toric_z_rank
    actual_sectors = (
        num_qubits,
        product_x_rank,
        product_z_rank,
        int(toric_x_rank),
    )
    expected_sectors = (
        case.reference_q,
        case.reference_p_x,
        case.reference_p_z,
        case.reference_t,
    )
    assert actual_sectors == expected_sectors

    h_z, h_x = split_css_blocks(
        result.input_excitation_map,
        num_x_checks,
        num_qubits,
    )
    psi = result.psi
    psi_inverse = result.psi_inverse
    assert psi is not None

    assert result.h_x_tilde * psi.psi_1 == psi.psi_0 * h_x
    assert (
        result.h_z_tilde_dagger * psi.psi_2
        == psi.psi_1 * dagger_matrix(h_z)
    )
    assert (
        h_x * psi_inverse.psi_1_inverse
        == psi_inverse.psi_0_inverse * result.h_x_tilde
    )
    assert (
        dagger_matrix(h_z) * psi_inverse.psi_2_inverse
        == psi_inverse.psi_1_inverse * result.h_z_tilde_dagger
    )

    for forward, backward in (
        (psi.psi_2, psi_inverse.psi_2_inverse),
        (psi.psi_1, psi_inverse.psi_1_inverse),
        (psi.psi_0, psi_inverse.psi_0_inverse),
    ):
        source_identity = identity_matrix(forward.base_ring(), forward.ncols())
        target_identity = identity_matrix(forward.base_ring(), forward.nrows())
        assert backward * forward == source_identity
        assert forward * backward == target_identity

    assert verify_qca_decoupling(result)


for paper_case in (COLOR_CODE_666, COLOR_CODE_488):
    verify_paper_example(paper_case)

print("test_paper_decoupling_end_to_end: ok")
