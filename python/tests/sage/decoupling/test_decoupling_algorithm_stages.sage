"""Verify the intermediate identities in the paper's decoupling algorithm."""

import _bootstrap  # noqa: F401

from decoupling import coarse_grain_to_superlattice, dagger_matrix
from decoupling.chain_maps.decoupling import (
    _PostEliminationMaps,
    _construct_phi_1,
    _construct_phi_1_prime,
    _construct_phi_2,
    _construct_phi_2_prime,
    _solve_homotopy_correction,
    compose_inverse_chain_isomorphisms,
)
from decoupling.chain_maps.elimination import clear_css_mod_j
from decoupling.css import (
    build_standard_blocks,
    evaluate_at_identity,
    split_css_blocks,
)
from decoupling.paper_examples import COLOR_CODE_666


coarse = coarse_grain_to_superlattice(
    COLOR_CODE_666.excitation_map(),
    *COLOR_CODE_666.superlattice_basis,
)
num_x_checks = coarse.nrows() // 2
num_qubits = coarse.ncols() // 2
clearing = clear_css_mod_j(
    coarse,
    num_x_checks=num_x_checks,
    num_qubits=num_qubits,
)

# The clearing routine must retain the exact row and symplectic column maps.
assert (
    clearing.working_matrix
    == clearing.row_map * coarse * clearing.symplectic_column_map
)

h_z, h_x = split_css_blocks(
    clearing.working_matrix,
    num_x_checks,
    num_qubits,
)
h_z_tilde, h_x_tilde = build_standard_blocks(
    clearing.working_matrix,
    num_x_checks,
    num_qubits,
    clearing.product_x_rank,
    clearing.product_z_rank,
)
assert evaluate_at_identity(h_x) == evaluate_at_identity(h_x_tilde)
assert evaluate_at_identity(h_z) == evaluate_at_identity(h_z_tilde)

h_z_dagger = dagger_matrix(h_z)
h_z_tilde_dagger = dagger_matrix(h_z_tilde)
phi_1 = _construct_phi_1(
    h_x,
    h_x_tilde,
    h_z_dagger,
    clearing.product_x_rank,
    clearing.product_z_rank,
    verify=True,
)
assert h_x * phi_1 == h_x_tilde
product_z_start = clearing.product_x_rank
product_z_stop = product_z_start + clearing.product_z_rank
assert (
    phi_1[:, product_z_start:product_z_stop]
    == h_z_dagger[:, : clearing.product_z_rank]
)

phi_2 = _construct_phi_2(
    h_z_dagger,
    phi_1,
    h_z_tilde_dagger,
    clearing.product_z_rank,
    verify=True,
)
assert h_z_dagger * phi_2 == phi_1 * h_z_tilde_dagger

phi_2_prime = _construct_phi_2_prime(
    phi_2,
    clearing.product_z_rank,
)
product_z_rank = clearing.product_z_rank
assert phi_2_prime[:product_z_rank, :] == phi_2[:product_z_rank, :]
assert (
    phi_2_prime[product_z_rank:, :product_z_rank]
    == phi_2[product_z_rank:, :product_z_rank]
)
assert (
    phi_2_prime[product_z_rank:, product_z_rank:]
    == evaluate_at_identity(phi_2[product_z_rank:, product_z_rank:])
)

difference = phi_2_prime - phi_2
eta = _solve_homotopy_correction(
    h_z_tilde_dagger,
    difference,
    product_z_rank,
    verify=True,
)
assert eta * h_z_tilde_dagger == difference

phi_1_prime = _construct_phi_1_prime(phi_1, h_z_dagger, eta)
assert phi_1_prime == phi_1 + h_z_dagger * eta
assert h_z_dagger * phi_2_prime == phi_1_prime * h_z_tilde_dagger
assert phi_1_prime.det().is_unit()
assert phi_2_prime.det().is_unit()

post_maps = _PostEliminationMaps(phi_2_prime, phi_1_prime)
psi_inverse = compose_inverse_chain_isomorphisms(clearing, post_maps)
input_h_z, input_h_x = split_css_blocks(
    coarse,
    num_x_checks,
    num_qubits,
)
assert (
    input_h_x * psi_inverse.psi_1_inverse
    == psi_inverse.psi_0_inverse * h_x_tilde
)
assert (
    dagger_matrix(input_h_z) * psi_inverse.psi_2_inverse
    == psi_inverse.psi_1_inverse * h_z_tilde_dagger
)

print("test_decoupling_algorithm_stages: ok")
