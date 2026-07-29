import _bootstrap  # noqa: F401

from dataclasses import fields, replace

from sage.all import Matrix, block_matrix, zero_matrix

from isomorphism import (
    DecouplingMaps,
    DecouplingResult,
    InverseDecouplingMaps,
    R,
    build_quotient_translation_action,
    choose_smallest_oblique_cell,
    construct_excitation_map,
    decouple_coarse_matrix,
    oblique_coarse_grain,
    periods_from_generators,
    periods_from_translation_action,
    x,
    y,
)
import isomorphism.css as css_module
import isomorphism.chain_maps.decoupling as decoupling_module
from isomorphism.chain_maps.decoupling import (
    _apply_homotopy_correction,
    _left_multiply_binary,
    _phi2_inverse_tail_target,
    _right_multiply_binary,
    compose_total_inverse_maps,
    construct_forward_maps,
    construct_post_elimination_inverse_maps,
    qca_decoupled_excitation_matrix,
    qca_symplectic_matrix,
    stabilizer_redefined_excitation_matrix,
    target_excitation_matrix,
    verify_qca_decoupling,
)
from isomorphism.chain_maps.elimination import clear_css_mod_j
from isomorphism.css import antipode_matrix, evaluate_at_identity, split_css_blocks


assert [item.name for item in fields(DecouplingMaps)] == [
    "phi2",
    "phi1",
    "phi0",
]
assert [item.name for item in fields(InverseDecouplingMaps)] == [
    "phi2_inverse",
    "phi1_inverse",
    "phi0_inverse",
]
assert [item.name for item in fields(DecouplingResult)] == [
    "source_coarse_matrix",
    "hx_standard",
    "hz_standard",
    "inverse_maps",
    "maps",
    "diagnostics",
]

binary = Matrix(R, [[1, 1], [0, 1]])
polynomial = Matrix(R, [[x, 1 + y], [y, 1 + x]])
assert _left_multiply_binary(binary, polynomial) == binary * polynomial
assert _right_multiply_binary(polynomial, binary) == polynomial * binary

epsilon = construct_excitation_map(1 + x, 1 + y)
clearing = clear_css_mod_j(epsilon, num_x_checks=1, num_qubits=2)
post_inverse = construct_post_elimination_inverse_maps(clearing, verify=True)
post_maps = construct_forward_maps(post_inverse, verify=True)

hz_working, hx_working = split_css_blocks(clearing.working_matrix, 1, 2)
hz_standard, hx_standard = css_module.build_standard_blocks(
    clearing.working_matrix,
    clearing.num_x_checks,
    clearing.num_qubits,
    clearing.product_x_rank,
    clearing.product_z_rank,
)
hz_working_dagger = antipode_matrix(hz_working)
hz_standard_dagger = antipode_matrix(hz_standard)
identity = post_inverse.phi1_inverse.parent().identity_matrix()
assert hx_working * post_inverse.phi1_inverse == (
    post_inverse.phi0_inverse * hx_standard
)
assert hz_working_dagger * post_inverse.phi2_inverse == (
    post_inverse.phi1_inverse * hz_standard_dagger
)
assert hx_standard * post_maps.phi1 == post_maps.phi0 * hx_working
assert hz_standard_dagger * post_maps.phi2 == post_maps.phi1 * hz_working_dagger
assert post_inverse.phi1_inverse * post_maps.phi1 == identity
assert post_maps.phi1 * post_inverse.phi1_inverse == identity
post_result = DecouplingResult(
    clearing.working_matrix,
    hx_standard,
    hz_standard,
    post_inverse,
    None,
    {
        "inverse_map_stage": "post_elimination",
        "product_z_rank": clearing.product_z_rank,
    },
)
assert verify_qca_decoupling(post_result)

verification_flags = []
original_solve_right = decoupling_module.solve_right
original_solve_left = decoupling_module.solve_left
original_find_inverse = decoupling_module.find_inverse


def solve_right_spy(A, B, *, verify=True):
    verification_flags.append(("solve_right", verify))
    return original_solve_right(A, B, verify=verify)


def solve_left_spy(A, B, *, verify=True):
    verification_flags.append(("solve_left", verify))
    return original_solve_left(A, B, verify=verify)


def find_inverse_spy(A, *, verify=True):
    verification_flags.append(("find_inverse", verify))
    return original_find_inverse(A, verify=verify)


decoupling_module.solve_right = solve_right_spy
decoupling_module.solve_left = solve_left_spy
decoupling_module.find_inverse = find_inverse_spy
try:
    fast_inverse = construct_post_elimination_inverse_maps(clearing, verify=False)
    fast_maps = construct_forward_maps(fast_inverse, verify=False)
finally:
    decoupling_module.solve_right = original_solve_right
    decoupling_module.solve_left = original_solve_left
    decoupling_module.find_inverse = original_find_inverse

assert verification_flags == [
    ("solve_right", False),
    ("solve_right", False),
    ("solve_left", False),
    ("find_inverse", False),
    ("find_inverse", False),
]
assert fast_inverse.phi2_inverse * fast_maps.phi2 == (
    fast_inverse.phi2_inverse.parent().identity_matrix()
)
assert fast_inverse.phi1_inverse * fast_maps.phi1 == (
    fast_inverse.phi1_inverse.parent().identity_matrix()
)


def reject_total_inverse(_matrix):
    raise AssertionError("composition should reuse the constructed inverse-direction factors")


decoupling_module.find_inverse = reject_total_inverse
try:
    composed_inverse = compose_total_inverse_maps(clearing, post_inverse)
finally:
    decoupling_module.find_inverse = original_find_inverse

assert composed_inverse.phi1_inverse == clearing.xi1 * post_inverse.phi1_inverse

total = decouple_coarse_matrix(epsilon, num_x_checks=1, num_qubits=2)
assert isinstance(total, DecouplingResult)
assert isinstance(total.inverse_maps, InverseDecouplingMaps)
assert isinstance(total.maps, DecouplingMaps)
hz_source, hx_source = split_css_blocks(epsilon, 1, 2)
hz_source_dagger = antipode_matrix(hz_source)
hz_standard_dagger = antipode_matrix(total.hz_standard)
assert total.hx_standard * total.maps.phi1 == total.maps.phi0 * hx_source
assert hz_standard_dagger * total.maps.phi2 == (
    total.maps.phi1 * hz_source_dagger
)
assert hx_source * total.inverse_maps.phi1_inverse == (
    total.inverse_maps.phi0_inverse * total.hx_standard
)
assert hz_source_dagger * total.inverse_maps.phi2_inverse == (
    total.inverse_maps.phi1_inverse
    * antipode_matrix(total.hz_standard)
)
for removed_name in ("input_coarse_matrix", "clearing", "post_inverse_maps"):
    assert not hasattr(total, removed_name)
for map_object in (total.maps, total.inverse_maps):
    for removed_name in ("hx_standard", "hz_standard", "diagnostics"):
        assert not hasattr(map_object, removed_name)

qca_matrix = qca_symplectic_matrix(total.maps, total.inverse_maps)
assert qca_matrix.nrows() == epsilon.ncols()
assert qca_decoupled_excitation_matrix(
    epsilon,
    total.maps,
    total.inverse_maps,
) == epsilon * qca_matrix
assert stabilizer_redefined_excitation_matrix(
    total,
) == target_excitation_matrix(total)
assert verify_qca_decoupling(total)

inverse_calls = []


def reject_forward_inverse(_matrix, *, verify=True):
    inverse_calls.append("find_inverse")
    raise AssertionError("compute_forward_maps=False should skip Laurent inverses")


decoupling_module.find_inverse = reject_forward_inverse
try:
    inverse_only = decouple_coarse_matrix(
        epsilon,
        num_x_checks=1,
        num_qubits=2,
        compute_forward_maps=False,
    )
finally:
    decoupling_module.find_inverse = original_find_inverse

assert inverse_calls == []
assert inverse_only.maps is None
assert isinstance(inverse_only.inverse_maps, InverseDecouplingMaps)
assert verify_qca_decoupling(inverse_only)

hx_standard_swapped = Matrix(
    R,
    [
        [1, 0, 0, 0],
        [0, 0, 1 + x, 1 + y],
    ],
)
hz_standard_swapped = Matrix(
    R,
    [
        [0, 1, 0, 0],
        [0, 0, 1 + y**-1, 1 + x**-1],
    ],
)
hz_source_swapped = Matrix(
    R,
    [
        hz_standard_swapped[1],
        hz_standard_swapped[0],
    ],
)
zero_swapped = zero_matrix(R, 2, 4)
epsilon_swapped = block_matrix(
    R,
    [
        [hx_standard_swapped, zero_swapped],
        [zero_swapped, hz_source_swapped],
    ],
)
decoupling_module.find_inverse = reject_forward_inverse
try:
    inverse_only_swapped = decouple_coarse_matrix(
        epsilon_swapped,
        num_x_checks=2,
        num_qubits=4,
        compute_forward_maps=False,
    )
    swapped_verifies = verify_qca_decoupling(
        inverse_only_swapped,
    )
finally:
    decoupling_module.find_inverse = original_find_inverse

assert inverse_only_swapped.diagnostics["clearing"].xi2 == Matrix(
    R,
    [[0, 1], [1, 0]],
)
assert swapped_verifies
assert inverse_calls == []

inverse_only_swapped_without_diagnostics = replace(
    inverse_only_swapped,
    diagnostics={},
)
try:
    verify_qca_decoupling(
        inverse_only_swapped_without_diagnostics,
        product_z_rank=inverse_only_swapped.diagnostics["product_z_rank"],
    )
except ValueError as exc:
    assert "inverse-only QCA verification requires clearing diagnostics" in str(exc)
    assert "post_phi2_inverse" in str(exc)
    assert "clearing_xi2" in str(exc)
else:
    raise AssertionError(
        "diagnostics-free total inverse maps should be rejected clearly"
    )


class RejectInverseMatrix:
    def __init__(self, matrix):
        self.matrix = matrix

    def inverse(self):
        raise AssertionError("inverse-only composition should not invert the clearing qubit map")

    def nrows(self):
        return self.matrix.nrows()

    def ncols(self):
        return self.matrix.ncols()

    def nonzero_positions_in_row(self, row):
        return self.matrix.nonzero_positions_in_row(row)

    def __getitem__(self, key):
        return self.matrix[key]


wrapped_clearing = replace(clearing, xi1=RejectInverseMatrix(clearing.xi1))
wrapped_total = compose_total_inverse_maps(wrapped_clearing, post_inverse)
assert wrapped_total.phi1_inverse == clearing.xi1 * post_inverse.phi1_inverse

original_antipode_poly = css_module.antipode_poly


def reject_zero_antipode(poly):
    if poly == 0:
        raise AssertionError("antipode_matrix should skip zero entries")
    return original_antipode_poly(poly)


css_module.antipode_poly = reject_zero_antipode
try:
    assert css_module.antipode_matrix(Matrix(R, [[0, x], [y, 0]])) == Matrix(
        R,
        [[0, y**-1], [x**-1, 0]],
    )
finally:
    css_module.antipode_poly = original_antipode_poly

original_evaluate_poly = css_module.evaluate_poly_at_identity


def reject_zero_evaluate(poly):
    if poly == 0:
        raise AssertionError("evaluate_at_identity should skip zero entries")
    return original_evaluate_poly(poly)


css_module.evaluate_poly_at_identity = reject_zero_evaluate
try:
    assert evaluate_at_identity(Matrix(R, [[0, x], [1 + x, 0]])) == Matrix(
        R,
        [[0, 1], [0, 0]],
    )
finally:
    css_module.evaluate_poly_at_identity = original_evaluate_poly

periods_666 = periods_from_generators(
    1 + x + x * y,
    1 + y + x * y,
    max_period=3,
)
cell_666 = choose_smallest_oblique_cell(periods_666.vectors)
coarse_666 = oblique_coarse_grain(
    construct_excitation_map(1 + x + x * y, 1 + y + x * y),
    *cell_666,
)
clearing_666 = clear_css_mod_j(coarse_666, num_x_checks=3, num_qubits=6)
hz_666, hx_666 = split_css_blocks(clearing_666.working_matrix, 3, 6)
hz_standard_666, hx_standard_666 = decoupling_module.build_standard_blocks(
    clearing_666.working_matrix,
    3,
    6,
    clearing_666.product_x_rank,
    clearing_666.product_z_rank,
)
hz_dagger_666 = antipode_matrix(hz_666)
hz_standard_dagger_666 = antipode_matrix(hz_standard_666)
phi1_inverse_666 = decoupling_module._construct_phi1_inverse(
    hx_666,
    hx_standard_666,
    hz_dagger_666,
    clearing_666.product_x_rank,
    clearing_666.product_z_rank,
    verify=False,
)
tail_standard_666 = _phi2_inverse_tail_target(
    phi1_inverse_666,
    hz_standard_dagger_666,
    clearing_666.product_z_rank,
)
assert tail_standard_666 == (
    phi1_inverse_666 * hz_standard_dagger_666
).matrix_from_columns(
    range(
        clearing_666.product_z_rank,
        hz_standard_dagger_666.ncols(),
    )
)
phi2_inverse_666 = decoupling_module._construct_phi2_inverse(
    hz_dagger_666,
    phi1_inverse_666,
    hz_standard_dagger_666,
    clearing_666.product_z_rank,
    verify=False,
)
corrected_phi2_inverse_666 = decoupling_module._phi2_inverse_target(
    phi2_inverse_666,
    clearing_666.product_z_rank,
)
correction_666 = decoupling_module._solve_homotopy_correction(
    hz_standard_dagger_666,
    corrected_phi2_inverse_666 - phi2_inverse_666,
    clearing_666.product_z_rank,
    verify=False,
)
assert _apply_homotopy_correction(
    phi1_inverse_666,
    hz_dagger_666,
    correction_666,
) == phi1_inverse_666 + hz_dagger_666 * correction_666

solver_shapes = []


def solve_right_shape_spy(A, B, *, verify=True):
    solver_shapes.append(("solve_right", A.ncols(), B.ncols()))
    return original_solve_right(A, B, verify=verify)


def solve_left_shape_spy(A, B, *, verify=True):
    solver_shapes.append(("solve_left", B.nrows(), B.ncols()))
    return original_solve_left(A, B, verify=verify)


decoupling_module.solve_right = solve_right_shape_spy
decoupling_module.solve_left = solve_left_shape_spy
try:
    construct_post_elimination_inverse_maps(clearing_666, verify=False)
finally:
    decoupling_module.solve_right = original_solve_right
    decoupling_module.solve_left = original_solve_left

assert solver_shapes[0] == ("solve_right", 6, 5)
assert solver_shapes[1] == ("solve_right", 3, 2)
assert solver_shapes[2] == ("solve_left", 2, 3)

decoupling_module.find_inverse = reject_forward_inverse
try:
    inverse_only_666 = decouple_coarse_matrix(
        coarse_666,
        num_x_checks=3,
        num_qubits=6,
        compute_forward_maps=False,
    )
finally:
    decoupling_module.find_inverse = original_find_inverse

assert verify_qca_decoupling(inverse_only_666)
maps_666 = decouple_coarse_matrix(
    coarse_666,
    num_x_checks=3,
    num_qubits=6,
)
assert maps_666.inverse_maps.phi1_inverse * maps_666.maps.phi1 == (
    maps_666.inverse_maps.phi1_inverse.parent().identity_matrix()
)

check_488 = Matrix(
    R,
    [
        [1 + y, 1 + y, 1 + x, 1 + x],
        [x * y, y, x * y, x],
    ],
)
zero = zero_matrix(R, 2, 4)
epsilon_488 = block_matrix(R, [[check_488, zero], [zero, check_488]])
action_488 = build_quotient_translation_action(epsilon_488, diagnostics=True)
periods_488 = periods_from_translation_action(
    action_488.tx,
    action_488.ty,
    max_period=4,
)
cell_488 = choose_smallest_oblique_cell(periods_488.vectors)
coarse_488 = oblique_coarse_grain(epsilon_488, *cell_488)
maps_488 = decouple_coarse_matrix(
    coarse_488,
    num_x_checks=4,
    num_qubits=8,
)
assert maps_488.inverse_maps.phi1_inverse * maps_488.maps.phi1 == (
    maps_488.inverse_maps.phi1_inverse.parent().identity_matrix()
)

bb_periods = periods_from_generators(
    1 + x + x**-1 * y,
    1 + y + x * y,
    max_period=7,
)
assert bb_periods.square_period == 7

print("test_isomorphism_decoupling: ok")
