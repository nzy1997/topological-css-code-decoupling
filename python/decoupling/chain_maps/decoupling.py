"""Construct the paper's decoupling unitary and chain isomorphisms."""

from dataclasses import dataclass

from sage.all import block_matrix, copy, identity_matrix, zero_matrix

from decoupling.chain_maps.elimination import clear_css_mod_j
from decoupling.chain_maps.result import (
    ChainIsomorphisms,
    DecouplingUnitaryResult,
    InverseChainIsomorphisms,
)
from decoupling.css import (
    dagger_matrix,
    build_standard_blocks,
    evaluate_at_identity,
    split_css_blocks,
)
from decoupling.rings import R
from decoupling.solver import find_inverse, solve_left, solve_right


@dataclass(frozen=True)
class _PostEliminationMaps:
    """Intermediate maps in the paper's post-elimination basis."""

    phi_2_prime: object
    phi_1_prime: object


def _construct_phi_2_prime(phi_2, product_z_rank):
    """Replace the toric lower-right block of ``phi_2`` modulo ``J``.

    Args:
        phi_2: Intermediate degree-two map.
        product_z_rank: Number of Z product rows.

    Returns:
        The paper's corrected degree-two map ``phi_2_prime``.
    """
    standard = copy(phi_2)
    # The algorithm fixes the toric tail modulo J before solving the
    # homotopy correction that restores the exact chain equation.
    standard[product_z_rank:, product_z_rank:] = evaluate_at_identity(
        standard[product_z_rank:, product_z_rank:]
    )
    return standard


def _construct_phi_1(
    h_x,
    h_x_tilde,
    h_z_dagger,
    product_x_rank,
    product_z_rank,
    *,
    verify,
):
    """Solve the non-overwritten columns of the inverse degree-one map.

    Args:
        h_x: Cleared input X-check matrix.
        h_x_tilde: Standard X-check matrix.
        h_z_dagger: Dagger of the cleared Z-check matrix.
        product_x_rank: Number of X product columns.
        product_z_rank: Number of Z product columns.
        verify: If true, verify Laurent solves.

    Returns:
        Laurent matrix ``phi_1`` before homotopy correction.
    """
    phi_1 = zero_matrix(R, h_x.ncols(), h_x_tilde.ncols())
    solved_columns = tuple(range(product_x_rank)) + tuple(
        range(product_x_rank + product_z_rank, h_x_tilde.ncols())
    )
    if solved_columns:
        # These columns are determined by the inverse degree-one chain
        # equation and can be solved as one right-module membership problem.
        solved_targets = h_x_tilde.matrix_from_columns(solved_columns)
        solved_phi_1 = solve_right(h_x, solved_targets, verify=verify)
        for source_col, target_col in enumerate(solved_columns):
            phi_1[:, target_col] = solved_phi_1[:, source_col]
    # Z product columns are overwritten by the stabilizer generators specified
    # by the decoupling-unitary construction.
    phi_1[:, product_x_rank : product_x_rank + product_z_rank] = (
        h_z_dagger[:, :product_z_rank]
    )
    return phi_1


def _solve_homotopy_correction(
    h_z_tilde_dagger,
    difference,
    product_z_rank,
    *,
    verify,
):
    """Solve only the toric rows where the homotopy correction is nonzero.

    Args:
        h_z_tilde_dagger: Dagger of the standard Z-check matrix.
        difference: Required change to the inverse degree-two map.
        product_z_rank: Number of Z product rows.
        verify: If true, verify Laurent solves.

    Returns:
        Homotopy correction matrix.
    """
    correction = zero_matrix(R, difference.nrows(), h_z_tilde_dagger.nrows())
    if product_z_rank == difference.nrows():
        return correction
    # Product rows have no correction.  Only the toric tail contributes.
    tail = difference[product_z_rank:, :]
    correction[product_z_rank:, :] = solve_left(
        h_z_tilde_dagger,
        tail,
        verify=verify,
    )
    return correction


def _construct_phi_2(
    h_z_dagger,
    phi_1,
    h_z_tilde_dagger,
    product_z_rank,
    *,
    verify,
):
    """Solve the toric tail columns of the inverse degree-two map.

    Args:
        h_z_dagger: Dagger of the cleared input Z-check matrix.
        phi_1: Pre-correction degree-one map.
        h_z_tilde_dagger: Dagger of the standard Z-check matrix.
        product_z_rank: Number of Z product rows.
        verify: If true, verify Laurent solves.

    Returns:
        Intermediate degree-two map ``phi_2``.
    """
    phi_2 = zero_matrix(R, h_z_dagger.ncols(), h_z_tilde_dagger.ncols())
    for col in range(product_z_rank):
        phi_2[col, col] = 1
    if product_z_rank == h_z_tilde_dagger.ncols():
        return phi_2
    # The toric tail is solved from the chain equation involving phi_1.
    target_tail = _phi_2_tail_target(
        phi_1,
        h_z_tilde_dagger,
        product_z_rank,
    )
    phi_2_tail = solve_right(h_z_dagger, target_tail, verify=verify)
    for source_col, target_col in enumerate(
        range(product_z_rank, h_z_tilde_dagger.ncols())
    ):
        phi_2[:, target_col] = phi_2_tail[:, source_col]
    return phi_2


def _phi_2_tail_target(phi_1, h_z_tilde_dagger, product_z_rank):
    """Build only the standard tail columns for the ``phi_2`` solve.

    Args:
        phi_1: Pre-correction degree-one map.
        h_z_tilde_dagger: Dagger of the standard Z-check matrix.
        product_z_rank: Number of Z product rows.

    Returns:
        Target matrix used in the ``phi_2`` right solve.
    """
    standard = zero_matrix(
        R,
        phi_1.nrows(),
        h_z_tilde_dagger.ncols() - product_z_rank,
    )
    for target_index, col in enumerate(
        range(product_z_rank, h_z_tilde_dagger.ncols())
    ):
        for row in h_z_tilde_dagger.nonzero_positions_in_column(col):
            # Expand the standard boundary through the inverse degree-one map
            # without forming a dense intermediate product.
            standard[:, target_index] += (
                h_z_tilde_dagger[row, col] * phi_1[:, row]
            )
    return standard


def _construct_phi_1_prime(phi_1, h_z_dagger, eta):
    """Apply ``h_z_dagger * correction`` using only nonzero correction rows.

    Args:
        phi_1: Pre-correction degree-one map.
        h_z_dagger: Dagger of the cleared input Z-check matrix.
        eta: Homotopy correction matrix.

    Returns:
        Corrected inverse degree-one map.
    """
    corrected = copy(phi_1)
    for row in range(eta.nrows()):
        for col in eta.nonzero_positions_in_row(row):
            corrected[:, col] += eta[row, col] * h_z_dagger[:, row]
    return corrected


def _left_multiply_binary(binary, mat):
    """Multiply ``binary * mat`` when ``binary`` has only 0/1 entries.

    Args:
        binary: 0/1 Laurent matrix on the left.
        mat: Laurent matrix on the right.

    Returns:
        Laurent product matrix.

    Raises:
        ValueError: If dimensions are incompatible or ``binary`` is not 0/1.
    """
    if binary.ncols() != mat.nrows():
        raise ValueError("Incompatible matrix dimensions for left multiplication.")
    one = R.one()
    result = zero_matrix(R, binary.nrows(), mat.ncols())
    for row in range(binary.nrows()):
        for col in sorted(binary.nonzero_positions_in_row(row)):
            if binary[row, col] != one:
                raise ValueError("Expected a binary 0/1 matrix.")
            # Binary row maps only select and add rows; this avoids expensive
            # generic Laurent multiplication.
            result[row] += mat[col]
    return result


def _right_multiply_binary(mat, binary):
    """Multiply ``mat * binary`` when ``binary`` has only 0/1 entries.

    Args:
        mat: Laurent matrix on the left.
        binary: 0/1 Laurent matrix on the right.

    Returns:
        Laurent product matrix.

    Raises:
        ValueError: If dimensions are incompatible or ``binary`` is not 0/1.
    """
    if mat.ncols() != binary.nrows():
        raise ValueError("Incompatible matrix dimensions for right multiplication.")
    one = R.one()
    result = zero_matrix(R, mat.nrows(), binary.ncols())
    for col in range(binary.ncols()):
        for row in sorted(binary.nonzero_positions_in_column(col)):
            if binary[row, col] != one:
                raise ValueError("Expected a binary 0/1 matrix.")
            # Binary column maps only select and add columns.
            result[:, col] += mat[:, row]
    return result


def _block_diagonal(left, right):
    """Return the block diagonal matrix with two Laurent blocks.

    Args:
        left: Upper-left Laurent block.
        right: Lower-right Laurent block.

    Returns:
        Block diagonal Laurent matrix.
    """
    upper_zero = zero_matrix(R, left.nrows(), right.ncols())
    lower_zero = zero_matrix(R, right.nrows(), left.ncols())
    return block_matrix(R, [[left, upper_zero], [lower_zero, right]])


def qca_symplectic_matrix(psi, psi_inverse):
    """Return the CSS symplectic matrix for the decoupling QCA.

    Args:
        psi: Forward source-to-standard chain isomorphisms.
        psi_inverse: Inverse standard-to-source chain isomorphisms.

    Returns:
        ``diag(psi_1_inverse, psi_1_dagger)``.
    """
    return _block_diagonal(
        psi_inverse.psi_1_inverse,
        dagger_matrix(psi.psi_1),
    )


def target_excitation_matrix(result):
    """Return the block-diagonal standard excitation matrix."""
    return _block_diagonal(
        result.h_x_tilde,
        dagger_matrix(result.h_z_tilde_dagger),
    )


def qca_decoupled_excitation_matrix(coarse_epsilon, psi, psi_inverse):
    """Apply the explicit forward/inverse QCA factors to the source matrix."""
    return coarse_epsilon * qca_symplectic_matrix(psi, psi_inverse)


def _qca_z_product_from_degree_two_map(
    h_z_tilde,
    degree_two_map,
    product_z_rank,
):
    """Return the QCA Z product from the inverse degree-two certificate.

    Args:
        h_z_tilde: Standard Z-check matrix.
        degree_two_map: Z-side chain-map certificate.
        product_z_rank: Number of Z product rows.

    Returns:
        Z-side QCA product without explicitly constructing ``phi1``.

    Raises:
        ValueError: If ``product_z_rank`` is missing or out of range.
    """
    if product_z_rank is None:
        raise ValueError("product_z_rank is required for inverse-only verification.")
    product_z_rank = int(product_z_rank)
    if not 0 <= product_z_rank <= h_z_tilde.nrows():
        raise ValueError("product_z_rank is outside the Z-check range.")
    toric_rank = h_z_tilde.nrows() - product_z_rank
    if toric_rank == 0:
        return copy(h_z_tilde)

    top_target = h_z_tilde[:product_z_rank, :]
    bottom_target = h_z_tilde[product_z_rank:, :]
    upper_right = degree_two_map[:product_z_rank, product_z_rank:]
    lower_right = degree_two_map[product_z_rank:, product_z_rank:]
    lower_right_dagger_inverse = dagger_matrix(lower_right).inverse()
    coupling = dagger_matrix(upper_right)
    # The triangular certificate gives the same lower block as the explicit
    # forward degree-one map without constructing that Laurent inverse.
    bottom_product = lower_right_dagger_inverse * (bottom_target + coupling * top_target)

    result = zero_matrix(R, h_z_tilde.nrows(), h_z_tilde.ncols())
    result[:product_z_rank, :] = top_target
    result[product_z_rank:, :] = bottom_product
    return result


def _qca_z_product_from_inverse_maps(
    psi_inverse,
    h_z_tilde,
    diagnostics,
    product_z_rank,
):
    """Return the inverse-only QCA Z product with clearing applied safely.

    The triangular certificate belongs to the post-elimination basis. Total
    composition may destroy that block form by left-multiplying ``xi_2``.
    After evaluating the post-elimination certificate, this helper transports
    it through clearing with ``(xi_2_dagger)^-1``. This last inverse is only a
    degree-zero binary basis change, not a Laurent forward decoupling map.

    Args:
        psi_inverse: Total standard-to-input chain isomorphisms.
        h_z_tilde: Standard Z-check matrix.
        diagnostics: Result-level construction diagnostics.
        product_z_rank: Number of Z product rows.

    Returns:
        Z-side QCA product in the source basis.

    Raises:
        ValueError: If total inverse maps lack clearing diagnostics.
    """
    if diagnostics.get("inverse_map_stage") == "post_elimination":
        return _qca_z_product_from_degree_two_map(
            h_z_tilde,
            psi_inverse.psi_2_inverse,
            product_z_rank,
        )

    post_phi_2_prime = diagnostics.get("post_phi_2_prime")
    clearing_xi_2 = diagnostics.get("clearing_xi_2")
    if post_phi_2_prime is None or clearing_xi_2 is None:
        raise ValueError(
            "inverse-only QCA verification requires clearing diagnostics "
            "'post_phi_2_prime' and 'clearing_xi_2' for total inverse maps."
        )

    post_product = _qca_z_product_from_degree_two_map(
        h_z_tilde,
        post_phi_2_prime,
        product_z_rank,
    )
    xi_2_dagger_inverse = dagger_matrix(clearing_xi_2).inverse()
    return _left_multiply_binary(xi_2_dagger_inverse, post_product)


def qca_decoupled_excitation_product(
    result,
    *,
    product_z_rank=None,
):
    """Return the QCA-right-multiplied excitation matrix.

    When forward maps are present, this forms the explicit symplectic product.
    Otherwise the inverse degree-two certificate supplies the Z-side product
    without computing a Laurent inverse.

    Args:
        result: Complete decoupling result.
        product_z_rank: Number of Z product rows, required only when
            it is not recorded in diagnostics.

    Returns:
        Full CSS matrix after applying the QCA column map.
    """
    psi_inverse = result.psi_inverse
    num_x_checks = result.h_x_tilde.nrows()
    num_qubits = result.h_x_tilde.ncols()
    h_z, h_x = split_css_blocks(
        result.input_excitation_map,
        num_x_checks,
        num_qubits,
    )
    x_product = h_x * psi_inverse.psi_1_inverse
    if result.psi is not None:
        z_product = h_z * dagger_matrix(result.psi.psi_1)
    else:
        if product_z_rank is None:
            product_z_rank = result.diagnostics.get("product_z_rank")
        z_product = _qca_z_product_from_inverse_maps(
            psi_inverse,
            dagger_matrix(result.h_z_tilde_dagger),
            result.diagnostics,
            product_z_rank,
        )
    zero_xz = zero_matrix(R, x_product.nrows(), z_product.ncols())
    zero_zx = zero_matrix(R, z_product.nrows(), x_product.ncols())
    return block_matrix(R, [[x_product, zero_xz], [zero_zx, z_product]])


def stabilizer_redefined_excitation_matrix(
    result,
    *,
    product_z_rank=None,
):
    """Apply the QCA and forward stabilizer-generator redefinitions.

    Args:
        result: Complete decoupling result.
        product_z_rank: Optional Z product rank for inverse-only verification.

    Returns:
        Excitation matrix after column QCA and row stabilizer redefinitions.
    """
    psi_inverse = result.psi_inverse
    psi_0 = (
        result.psi.psi_0
        if result.psi is not None
        else psi_inverse.psi_0_inverse.inverse()
    )
    row_redefinition = _block_diagonal(
        psi_0,
        dagger_matrix(psi_inverse.psi_2_inverse),
    )
    return row_redefinition * qca_decoupled_excitation_product(
        result,
        product_z_rank=product_z_rank,
    )


def verify_qca_decoupling(
    result,
    *,
    product_z_rank=None,
):
    """Return whether the QCA-decoupled excitation map equals the standard map.

    Args:
        result: Complete decoupling result.
        product_z_rank: Optional Z product rank for inverse-only verification.

    Returns:
        True when the redefined QCA product equals the standard matrix.
    """
    return stabilizer_redefined_excitation_matrix(
        result,
        product_z_rank=product_z_rank,
    ) == target_excitation_matrix(result)


def construct_post_elimination_maps(clearing, *, verify=False):
    """Construct the paper's post-elimination intermediate maps.

    Set ``verify=True`` to multiply each Laurent-solver witness back into its
    defining equation while debugging. The default matches the benchmarked
    production path and verifies final equations in tests instead.

    Args:
        clearing: Result of finite-field CSS clearing.
        verify: If true, verify every Laurent solve.

    Returns:
        Intermediate ``phi_2_prime`` and ``phi_1_prime``.
    """
    h_z, h_x = split_css_blocks(
        clearing.working_matrix,
        clearing.num_x_checks,
        clearing.num_qubits,
    )
    h_z_tilde, h_x_tilde = build_standard_blocks(
        clearing.working_matrix,
        clearing.num_x_checks,
        clearing.num_qubits,
        clearing.product_x_rank,
        clearing.product_z_rank,
    )
    h_z_dagger = dagger_matrix(h_z)
    h_z_tilde_dagger = dagger_matrix(h_z_tilde)

    phi_1 = _construct_phi_1(
        h_x,
        h_x_tilde,
        h_z_dagger,
        clearing.product_x_rank,
        clearing.product_z_rank,
        verify=verify,
    )
    phi_2 = _construct_phi_2(
        h_z_dagger,
        phi_1,
        h_z_tilde_dagger,
        clearing.product_z_rank,
        verify=verify,
    )
    phi_2_prime = _construct_phi_2_prime(
        phi_2,
        clearing.product_z_rank,
    )
    eta = _solve_homotopy_correction(
        h_z_tilde_dagger,
        phi_2_prime - phi_2,
        clearing.product_z_rank,
        verify=verify,
    )
    phi_1_prime = _construct_phi_1_prime(
        phi_1,
        h_z_dagger,
        eta,
    )
    return _PostEliminationMaps(
        phi_2_prime,
        phi_1_prime,
    )


def compose_inverse_chain_isomorphisms(clearing, post_maps):
    """Compose inverse clearing and post maps into the input basis.

    Args:
        clearing: Finite-field clearing result.
        post_maps: Paper intermediate maps valid on the cleared matrix.

    Returns:
        ``InverseChainIsomorphisms`` valid on the original coarse matrix.
    """
    psi_1_inverse = _left_multiply_binary(
        clearing.xi_1,
        post_maps.phi_1_prime,
    )
    psi_2_inverse = _left_multiply_binary(
        clearing.xi_2,
        post_maps.phi_2_prime,
    )
    return InverseChainIsomorphisms(
        psi_2_inverse,
        psi_1_inverse,
        clearing.xi_0,
    )


def construct_chain_isomorphisms(psi_inverse, *, verify=False):
    """Invert the directly constructed maps only when explicitly requested."""
    return ChainIsomorphisms(
        find_inverse(psi_inverse.psi_2_inverse, verify=verify),
        find_inverse(psi_inverse.psi_1_inverse, verify=verify),
        psi_inverse.psi_0_inverse.inverse(),
    )


def solve_decoupling_unitary(
    input_excitation_map,
    *,
    num_x_checks,
    num_qubits,
    verify=False,
    compute_psi=False,
):
    """Construct the inverse chain isomorphisms of the paper's algorithm.

    Args:
        input_excitation_map: Coarse-grained full CSS excitation map.
        num_x_checks: Number of X-check rows.
        num_qubits: Number of data-qubit columns per CSS half.
        verify: If true, verify Laurent solver witnesses.
        compute_psi: If true, also invert the Laurent maps to construct
            the forward chain isomorphisms.

    Returns:
        Result containing inverse maps always and optional forward maps.
    """
    # Stage 1: choose finite-field pivots and accumulate row/column maps.
    clearing = clear_css_mod_j(
        input_excitation_map,
        num_x_checks=num_x_checks,
        num_qubits=num_qubits,
    )
    # Stage 2: solve the Laurent chain-map equations on the cleared matrix.
    post_maps = construct_post_elimination_maps(
        clearing,
        verify=verify,
    )
    psi_inverse = compose_inverse_chain_isomorphisms(clearing, post_maps)
    h_z_tilde, h_x_tilde = build_standard_blocks(
        clearing.working_matrix,
        clearing.num_x_checks,
        clearing.num_qubits,
        clearing.product_x_rank,
        clearing.product_z_rank,
    )
    psi = (
        construct_chain_isomorphisms(psi_inverse, verify=verify)
        if compute_psi
        else None
    )
    diagnostics = {
        "inverse_map_stage": "total",
        "clearing": clearing,
        "post_elimination_maps": post_maps,
        "product_x_rank": clearing.product_x_rank,
        "product_z_rank": clearing.product_z_rank,
        "post_phi_2_prime": post_maps.phi_2_prime,
        "clearing_xi_2": clearing.xi_2,
    }
    return DecouplingUnitaryResult(
        input_excitation_map,
        h_x_tilde,
        dagger_matrix(h_z_tilde),
        psi_inverse,
        psi,
        diagnostics,
    )
