"""Algorithm SM.2 post maps and total map composition."""

from sage.all import block_matrix, copy, identity_matrix, zero_matrix

from isomorphism.chain_maps.elimination import clear_css_mod_j
from isomorphism.chain_maps.result import (
    DecouplingMaps,
    DecouplingResult,
    InverseDecouplingMaps,
)
from isomorphism.css import (
    antipode_matrix,
    build_standard_blocks,
    evaluate_at_identity,
    split_css_blocks,
)
from isomorphism.rings import R
from isomorphism.solver import find_inverse, solve_left, solve_right


def _phi2_inverse_target(phi2_inverse, product_z_rank):
    """Replace the toric lower-right block by its value mod ``J``.

    Args:
        phi2_inverse: Intermediate phi2_inverse map.
        product_z_rank: Number of Z product rows.

    Returns:
        Copy of ``phi2_inverse`` with the toric tail evaluated at ``x=y=1``.
    """
    standard = copy(phi2_inverse)
    # SM.2 fixes the toric tail modulo J before solving the homotopy
    # correction that restores the exact chain equation.
    standard[product_z_rank:, product_z_rank:] = evaluate_at_identity(
        standard[product_z_rank:, product_z_rank:]
    )
    return standard


def _construct_phi1_inverse(
    HX,
    hx_standard,
    HZ_dagger,
    product_x_rank,
    product_z_rank,
    *,
    verify,
):
    """Solve the non-overwritten columns of the inverse degree-one map.

    Args:
        HX: Cleared input X-check matrix.
        hx_standard: Standard X-check matrix.
        HZ_dagger: Dagger of the cleared Z-check matrix.
        product_x_rank: Number of X product columns.
        product_z_rank: Number of Z product columns.
        verify: If true, verify Laurent solves.

    Returns:
        Laurent matrix ``phi1_inverse`` before homotopy correction.
    """
    phi1_inverse = zero_matrix(R, HX.ncols(), hx_standard.ncols())
    solved_columns = tuple(range(product_x_rank)) + tuple(
        range(product_x_rank + product_z_rank, hx_standard.ncols())
    )
    if solved_columns:
        # These columns are determined by the inverse degree-one chain
        # equation and can be solved as one right-module membership problem.
        solved_targets = hx_standard.matrix_from_columns(solved_columns)
        solved_phi1_inverse = solve_right(HX, solved_targets, verify=verify)
        for source_col, target_col in enumerate(solved_columns):
            phi1_inverse[:, target_col] = solved_phi1_inverse[:, source_col]
    # Z product columns are overwritten by the stabilizer generators specified
    # in Algorithm SM.2.
    phi1_inverse[:, product_x_rank : product_x_rank + product_z_rank] = (
        HZ_dagger[:, :product_z_rank]
    )
    return phi1_inverse


def _solve_homotopy_correction(
    hz_standard_dagger,
    difference,
    product_z_rank,
    *,
    verify,
):
    """Solve only the toric rows where the homotopy correction is nonzero.

    Args:
        hz_standard_dagger: Dagger of the standard Z-check matrix.
        difference: Required change to the inverse degree-two map.
        product_z_rank: Number of Z product rows.
        verify: If true, verify Laurent solves.

    Returns:
        Homotopy correction matrix.
    """
    correction = zero_matrix(R, difference.nrows(), hz_standard_dagger.nrows())
    if product_z_rank == difference.nrows():
        return correction
    # Product rows have no correction.  Only the toric tail contributes.
    tail = difference[product_z_rank:, :]
    correction[product_z_rank:, :] = solve_left(
        hz_standard_dagger,
        tail,
        verify=verify,
    )
    return correction


def _construct_phi2_inverse(
    HZ_dagger,
    phi1_inverse,
    hz_standard_dagger,
    product_z_rank,
    *,
    verify,
):
    """Solve the toric tail columns of the inverse degree-two map.

    Args:
        HZ_dagger: Dagger of the cleared input Z-check matrix.
        phi1_inverse: Pre-correction inverse degree-one map.
        hz_standard_dagger: Dagger of the standard Z-check matrix.
        product_z_rank: Number of Z product rows.
        verify: If true, verify Laurent solves.

    Returns:
        Intermediate phi2_inverse map.
    """
    phi2_inverse = zero_matrix(R, HZ_dagger.ncols(), hz_standard_dagger.ncols())
    for col in range(product_z_rank):
        phi2_inverse[col, col] = 1
    if product_z_rank == hz_standard_dagger.ncols():
        return phi2_inverse
    # The toric tail is solved from the chain equation involving phi1_inverse.
    target_tail = _phi2_inverse_tail_target(
        phi1_inverse,
        hz_standard_dagger,
        product_z_rank,
    )
    phi2_inverse_tail = solve_right(HZ_dagger, target_tail, verify=verify)
    for source_col, target_col in enumerate(
        range(product_z_rank, hz_standard_dagger.ncols())
    ):
        phi2_inverse[:, target_col] = phi2_inverse_tail[:, source_col]
    return phi2_inverse


def _phi2_inverse_tail_target(phi1_inverse, hz_standard_dagger, product_z_rank):
    """Build only the standard tail columns for the phi2_inverse solve.

    Args:
        phi1_inverse: Pre-correction inverse degree-one map.
        hz_standard_dagger: Dagger of the standard Z-check matrix.
        product_z_rank: Number of Z product rows.

    Returns:
        Target matrix used in the phi2_inverse right solve.
    """
    standard = zero_matrix(
        R,
        phi1_inverse.nrows(),
        hz_standard_dagger.ncols() - product_z_rank,
    )
    for target_index, col in enumerate(
        range(product_z_rank, hz_standard_dagger.ncols())
    ):
        for row in hz_standard_dagger.nonzero_positions_in_column(col):
            # Expand the standard boundary through the inverse degree-one map
            # without forming a dense intermediate product.
            standard[:, target_index] += (
                hz_standard_dagger[row, col] * phi1_inverse[:, row]
            )
    return standard


def _apply_homotopy_correction(phi1_inverse, HZ_dagger, correction):
    """Apply ``HZ_dagger * correction`` using only nonzero correction rows.

    Args:
        phi1_inverse: Pre-correction inverse degree-one map.
        HZ_dagger: Dagger of the cleared input Z-check matrix.
        correction: Homotopy correction matrix.

    Returns:
        Corrected inverse degree-one map.
    """
    corrected = copy(phi1_inverse)
    for row in range(correction.nrows()):
        for col in correction.nonzero_positions_in_row(row):
            corrected[:, col] += correction[row, col] * HZ_dagger[:, row]
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


def qca_symplectic_matrix(maps, inverse_maps):
    """Return the CSS symplectic matrix for the decoupling QCA.

    Args:
        maps: Forward source-to-standard maps.
        inverse_maps: Inverse standard-to-source maps.

    Returns:
        ``diag(phi1_inverse, phi1_dagger)``.
    """
    return _block_diagonal(
        inverse_maps.phi1_inverse,
        antipode_matrix(maps.phi1),
    )


def target_excitation_matrix(result):
    """Return the block-diagonal standard excitation matrix."""
    return _block_diagonal(result.hx_standard, result.hz_standard)


def qca_decoupled_excitation_matrix(coarse_epsilon, maps, inverse_maps):
    """Apply the explicit forward/inverse QCA factors to the source matrix."""
    return coarse_epsilon * qca_symplectic_matrix(maps, inverse_maps)


def _qca_z_product_from_phi2_inverse(
    hz_standard,
    phi2_inverse,
    product_z_rank,
):
    """Return the QCA Z product from the inverse degree-two certificate.

    Args:
        hz_standard: Standard Z-check matrix.
        phi2_inverse: Z-side chain map certificate.
        product_z_rank: Number of Z product rows.

    Returns:
        Z-side QCA product without explicitly constructing ``phi1``.

    Raises:
        ValueError: If ``product_z_rank`` is missing or out of range.
    """
    if product_z_rank is None:
        raise ValueError("product_z_rank is required when phi1 is not available.")
    product_z_rank = int(product_z_rank)
    if not 0 <= product_z_rank <= hz_standard.nrows():
        raise ValueError("product_z_rank is outside the Z-check range.")
    toric_rank = hz_standard.nrows() - product_z_rank
    if toric_rank == 0:
        return copy(hz_standard)

    top_target = hz_standard[:product_z_rank, :]
    bottom_target = hz_standard[product_z_rank:, :]
    upper_right = phi2_inverse[:product_z_rank, product_z_rank:]
    lower_right = phi2_inverse[product_z_rank:, product_z_rank:]
    lower_right_dagger_inverse = antipode_matrix(lower_right).inverse()
    coupling = antipode_matrix(upper_right)
    # The triangular certificate gives the same lower block as the explicit
    # forward degree-one map without constructing that Laurent inverse.
    bottom_product = lower_right_dagger_inverse * (bottom_target + coupling * top_target)

    result = zero_matrix(R, hz_standard.nrows(), hz_standard.ncols())
    result[:product_z_rank, :] = top_target
    result[product_z_rank:, :] = bottom_product
    return result


def _qca_z_product_from_inverse_maps(
    inverse_maps,
    hz_standard,
    diagnostics,
    product_z_rank,
):
    """Return the inverse-only QCA Z product with clearing applied safely.

    The triangular certificate belongs to the post-elimination basis. Total
    composition may destroy that block form by left-multiplying ``xi2``.
    After evaluating the post-elimination certificate, this helper transports
    it through clearing with ``(xi2_dagger)^-1``. This last inverse is only a
    degree-zero binary basis change, not a Laurent forward decoupling map.

    Args:
        inverse_maps: Total standard-to-source maps.
        hz_standard: Standard Z-check matrix.
        diagnostics: Result-level construction diagnostics.
        product_z_rank: Number of Z product rows.

    Returns:
        Z-side QCA product in the source basis.

    Raises:
        ValueError: If total inverse maps lack clearing diagnostics.
    """
    if diagnostics.get("inverse_map_stage") == "post_elimination":
        return _qca_z_product_from_phi2_inverse(
            hz_standard,
            inverse_maps.phi2_inverse,
            product_z_rank,
        )

    post_phi2_inverse = diagnostics.get("post_phi2_inverse")
    clearing_xi2 = diagnostics.get("clearing_xi2")
    if post_phi2_inverse is None or clearing_xi2 is None:
        raise ValueError(
            "inverse-only QCA verification requires clearing diagnostics "
            "'post_phi2_inverse' and 'clearing_xi2' for total inverse maps."
        )

    post_product = _qca_z_product_from_phi2_inverse(
        hz_standard,
        post_phi2_inverse,
        product_z_rank,
    )
    xi2_dagger_inverse = antipode_matrix(clearing_xi2).inverse()
    return _left_multiply_binary(xi2_dagger_inverse, post_product)


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
    inverse_maps = result.inverse_maps
    num_x_checks = result.hx_standard.nrows()
    num_qubits = result.hx_standard.ncols()
    HZ, HX = split_css_blocks(
        result.source_coarse_matrix,
        num_x_checks,
        num_qubits,
    )
    x_product = HX * inverse_maps.phi1_inverse
    if result.maps is not None:
        z_product = HZ * antipode_matrix(result.maps.phi1)
    else:
        if product_z_rank is None:
            product_z_rank = result.diagnostics.get("product_z_rank")
        z_product = _qca_z_product_from_inverse_maps(
            inverse_maps,
            result.hz_standard,
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
    inverse_maps = result.inverse_maps
    phi0 = (
        result.maps.phi0
        if result.maps is not None
        else inverse_maps.phi0_inverse.inverse()
    )
    row_redefinition = _block_diagonal(
        phi0,
        antipode_matrix(inverse_maps.phi2_inverse),
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


def construct_post_elimination_inverse_maps(clearing, *, verify=False):
    """Construct inverse-direction maps on the cleared working matrix.

    Set ``verify=True`` to multiply each Laurent-solver witness back into its
    defining equation while debugging. The default matches the benchmarked
    production path and verifies final equations in tests instead.

    Args:
        clearing: Result of finite-field CSS clearing.
        verify: If true, verify every Laurent solve.

    Returns:
        ``InverseDecouplingMaps`` valid for ``clearing.working_matrix``.
    """
    HZ, HX = split_css_blocks(
        clearing.working_matrix,
        clearing.num_x_checks,
        clearing.num_qubits,
    )
    hz_standard, hx_standard = build_standard_blocks(
        clearing.working_matrix,
        clearing.num_x_checks,
        clearing.num_qubits,
        clearing.product_x_rank,
        clearing.product_z_rank,
    )
    HZ_dagger = antipode_matrix(HZ)
    hz_standard_dagger = antipode_matrix(hz_standard)

    phi1_inverse = _construct_phi1_inverse(
        HX,
        hx_standard,
        HZ_dagger,
        clearing.product_x_rank,
        clearing.product_z_rank,
        verify=verify,
    )
    phi2_inverse = _construct_phi2_inverse(
        HZ_dagger,
        phi1_inverse,
        hz_standard_dagger,
        clearing.product_z_rank,
        verify=verify,
    )
    corrected_phi2_inverse = _phi2_inverse_target(
        phi2_inverse,
        clearing.product_z_rank,
    )
    correction = _solve_homotopy_correction(
        hz_standard_dagger,
        corrected_phi2_inverse - phi2_inverse,
        clearing.product_z_rank,
        verify=verify,
    )
    corrected_phi1_inverse = _apply_homotopy_correction(
        phi1_inverse,
        HZ_dagger,
        correction,
    )
    phi0_inverse = identity_matrix(R, HX.nrows())
    return InverseDecouplingMaps(
        corrected_phi2_inverse,
        corrected_phi1_inverse,
        phi0_inverse,
    )


def compose_total_inverse_maps(clearing, post_inverse_maps):
    """Compose inverse clearing and post maps into the input basis.

    Args:
        clearing: Finite-field clearing result.
        post_inverse_maps: Inverse maps valid on the cleared matrix.

    Returns:
        ``InverseDecouplingMaps`` valid on the original coarse matrix.
    """
    phi1_inverse = _left_multiply_binary(
        clearing.xi1,
        post_inverse_maps.phi1_inverse,
    )
    phi2_inverse = _left_multiply_binary(
        clearing.xi2,
        post_inverse_maps.phi2_inverse,
    )
    return InverseDecouplingMaps(
        phi2_inverse,
        phi1_inverse,
        clearing.xi0,
    )


def construct_forward_maps(inverse_maps, *, verify=False):
    """Invert directly constructed maps into the supplement direction."""
    return DecouplingMaps(
        find_inverse(inverse_maps.phi2_inverse, verify=verify),
        find_inverse(inverse_maps.phi1_inverse, verify=verify),
        inverse_maps.phi0_inverse.inverse(),
    )


def decouple_coarse_matrix(
    coarse_epsilon,
    *,
    num_x_checks,
    num_qubits,
    verify=False,
    compute_forward_maps=True,
):
    """Run accumulated clearing plus Algorithm SM.2 post-map construction.

    Args:
        coarse_epsilon: Coarse-grained full CSS excitation matrix.
        num_x_checks: Number of X-check rows.
        num_qubits: Number of data-qubit columns per CSS half.
        verify: If true, verify Laurent solver witnesses.
        compute_forward_maps: If true, compute the source-to-standard maps.

    Returns:
        Result containing inverse maps always and optional forward maps.
    """
    # Stage 1: choose finite-field pivots and accumulate row/column maps.
    clearing = clear_css_mod_j(
        coarse_epsilon,
        num_x_checks=num_x_checks,
        num_qubits=num_qubits,
    )
    # Stage 2: solve the Laurent chain-map equations on the cleared matrix.
    post_inverse_maps = construct_post_elimination_inverse_maps(
        clearing,
        verify=verify,
    )
    inverse_maps = compose_total_inverse_maps(clearing, post_inverse_maps)
    hz_standard, hx_standard = build_standard_blocks(
        clearing.working_matrix,
        clearing.num_x_checks,
        clearing.num_qubits,
        clearing.product_x_rank,
        clearing.product_z_rank,
    )
    maps = (
        construct_forward_maps(inverse_maps, verify=verify)
        if compute_forward_maps
        else None
    )
    diagnostics = {
        "inverse_map_stage": "total",
        "clearing": clearing,
        "post_inverse_maps": post_inverse_maps,
        "product_x_rank": clearing.product_x_rank,
        "product_z_rank": clearing.product_z_rank,
        "post_phi2_inverse": post_inverse_maps.phi2_inverse,
        "clearing_xi2": clearing.xi2,
    }
    return DecouplingResult(
        coarse_epsilon,
        hx_standard,
        hz_standard,
        inverse_maps,
        maps,
        diagnostics,
    )
