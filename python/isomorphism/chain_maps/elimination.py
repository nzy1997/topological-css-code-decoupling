"""Finite-field CSS clearing with accumulated inverse chain maps."""

from sage.all import GF, block_matrix, copy, identity_matrix, zero_matrix

from isomorphism.chain_maps.result import ClearingResult
from isomorphism.css import antipode_matrix, evaluate_at_identity
from isomorphism.rings import R


def _row_add(mat, dst, src):
    """Add one row into another matrix row in place.

    Args:
        mat: Mutable Sage matrix.
        dst: Destination row index.
        src: Source row index.

    Returns:
        None. The matrix is modified in place.
    """
    mat.add_multiple_of_row(dst, src, 1)


def _row_swap(mat, left, right):
    """Swap two matrix rows in place.

    Args:
        mat: Mutable Sage matrix.
        left: First row index.
        right: Second row index.

    Returns:
        None. The matrix is modified in place.
    """
    mat.swap_rows(left, right)


def _row_add_with_map(mat, row_map, dst, src, row_offset):
    """Apply one row addition to data and its accumulated row map.

    Args:
        mat: Working Laurent matrix.
        row_map: Row-map block for the active CSS half.
        dst: Destination row in ``mat``.
        src: Source row in ``mat``.
        row_offset: Offset from full CSS rows to the local row-map block.

    Returns:
        None. Both matrices are modified in place.
    """
    _row_add(mat, dst, src)
    _row_add(row_map, dst - row_offset, src - row_offset)


def _row_swap_with_map(mat, row_map, left, right, row_offset):
    """Apply one row swap to data and its accumulated row map.

    Args:
        mat: Working Laurent matrix.
        row_map: Row-map block for the active CSS half.
        left: First row in ``mat``.
        right: Second row in ``mat``.
        row_offset: Offset from full CSS rows to the local row-map block.

    Returns:
        None. Both matrices are modified in place.
    """
    _row_swap(mat, left, right)
    _row_swap(row_map, left - row_offset, right - row_offset)


def _is_left(col, q):
    """Return whether a CSS column belongs to the left half.

    Args:
        col: Full CSS column index.
        q: Number of columns per CSS half.

    Returns:
        True for the left/X half, false for the right/Z half.
    """
    return col < q


def _col_add(mat, dst, src, q):
    """Apply one symplectic CSS column addition in place.

    Args:
        mat: Mutable full CSS matrix.
        dst: Destination column in one CSS half.
        src: Source column in the same CSS half.
        q: Number of columns per CSS half.

    Returns:
        None. The matrix is modified in place.

    Raises:
        ValueError: If ``dst`` and ``src`` are in different CSS halves.
    """
    if _is_left(dst, q) != _is_left(src, q):
        raise ValueError("CSS clearing cannot mix symplectic halves.")
    mat.add_multiple_of_column(dst, src, 1)
    partner_dst = dst + q if _is_left(dst, q) else dst - q
    partner_src = src + q if _is_left(src, q) else src - q
    # Preserve the CSS symplectic form by applying the contragredient update on
    # partner columns.
    mat.add_multiple_of_column(partner_src, partner_dst, 1)


def _col_swap(mat, left, right, q):
    """Apply one symplectic CSS column swap in place.

    Args:
        mat: Mutable full CSS matrix.
        left: First column in one CSS half.
        right: Second column in the same CSS half.
        q: Number of columns per CSS half.

    Returns:
        None. The matrix is modified in place.

    Raises:
        ValueError: If ``left`` and ``right`` are in different CSS halves.
    """
    if _is_left(left, q) != _is_left(right, q):
        raise ValueError("CSS clearing cannot mix symplectic halves.")
    mat.swap_columns(left, right)
    partner_left = left + q if _is_left(left, q) else left - q
    partner_right = right + q if _is_left(right, q) else right - q
    mat.swap_columns(partner_left, partner_right)


def _col_add_halves(phi, psi, dst, src, q):
    """Apply one accumulated CSS column addition to block-diagonal halves.

    Args:
        phi: Accumulated X-half column map.
        psi: Accumulated Z-half column map.
        dst: Destination full CSS column.
        src: Source full CSS column.
        q: Number of columns per CSS half.

    Returns:
        None. ``phi`` and ``psi`` are modified in place.

    Raises:
        ValueError: If ``dst`` and ``src`` are in different CSS halves.
    """
    if _is_left(dst, q) != _is_left(src, q):
        raise ValueError("CSS clearing cannot mix symplectic halves.")
    if _is_left(dst, q):
        # X-half additions update phi directly and psi contragrediently.
        phi[:, dst] += phi[:, src]
        psi[:, src] += psi[:, dst]
    else:
        local_dst = dst - q
        local_src = src - q
        psi[:, local_dst] += psi[:, local_src]
        phi[:, local_src] += phi[:, local_dst]


def _col_swap_halves(phi, psi, left, right, q):
    """Apply one accumulated CSS column swap to block-diagonal halves.

    Args:
        phi: Accumulated X-half column map.
        psi: Accumulated Z-half column map.
        left: First full CSS column.
        right: Second full CSS column.
        q: Number of columns per CSS half.

    Returns:
        None. ``phi`` and ``psi`` are modified in place.

    Raises:
        ValueError: If ``left`` and ``right`` are in different CSS halves.
    """
    if _is_left(left, q) != _is_left(right, q):
        raise ValueError("CSS clearing cannot mix symplectic halves.")
    local_left = left if _is_left(left, q) else left - q
    local_right = right if _is_left(right, q) else right - q
    phi.swap_columns(local_left, local_right)
    psi.swap_columns(local_left, local_right)


def _col_add_with_map(mat, phi, psi, dst, src, q):
    """Apply one CSS column addition to data and its column map.

    Args:
        mat: Working Laurent matrix.
        phi: Accumulated X-half column map.
        psi: Accumulated Z-half column map.
        dst: Destination full CSS column.
        src: Source full CSS column.
        q: Number of columns per CSS half.

    Returns:
        None. The matrix and accumulated maps are modified in place.
    """
    _col_add(mat, dst, src, q)
    _col_add_halves(phi, psi, dst, src, q)


def _col_swap_with_map(mat, phi, psi, left, right, q):
    """Apply one CSS column swap to data and its column map.

    Args:
        mat: Working Laurent matrix.
        phi: Accumulated X-half column map.
        psi: Accumulated Z-half column map.
        left: First full CSS column.
        right: Second full CSS column.
        q: Number of columns per CSS half.

    Returns:
        None. The matrix and accumulated maps are modified in place.
    """
    _col_swap(mat, left, right, q)
    _col_swap_halves(phi, psi, left, right, q)


def _find_pivot(evaluated, pivot_row, pivot_col, stop_row, stop_col):
    """Find the next pivot in the selected evaluated rectangle.

    Args:
        evaluated: GF(2) matrix evaluated at ``x=y=1``.
        pivot_row: First candidate row.
        pivot_col: First candidate column.
        stop_row: Exclusive row bound.
        stop_col: Exclusive column bound.

    Returns:
        Pair ``(row, col)`` or ``None`` if no pivot exists.
    """
    for row in range(pivot_row, stop_row):
        for col in sorted(evaluated.nonzero_positions_in_row(row)):
            if pivot_col <= col < stop_col:
                return row, col
    return None


def _clear_block(
    evaluated,
    working,
    row_map,
    phi,
    psi,
    *,
    q,
    row_offset,
    pivot_row,
    pivot_col,
    stop_row,
    stop_col,
):
    """Clear one CSS half and return the rank exposed by pivots.

    Args:
        evaluated: GF(2) matrix used for pivot decisions.
        working: Laurent matrix receiving matching row/column operations.
        row_map: Accumulated row map for the active CSS half.
        phi: Accumulated X-half column map.
        psi: Accumulated Z-half column map.
        q: Number of columns per CSS half.
        row_offset: Offset from full CSS rows to local row-map rows.
        pivot_row: Initial pivot row.
        pivot_col: Initial pivot column.
        stop_row: Exclusive row bound for this block.
        stop_col: Exclusive column bound for this block.

    Returns:
        Rank exposed by pivot clearing in the selected block.
    """
    rank = 0
    while pivot_row < stop_row and pivot_col < stop_col:
        pivot = _find_pivot(evaluated, pivot_row, pivot_col, stop_row, stop_col)
        if pivot is None:
            break
        row, col = pivot
        # Move the pivot into place in both the evaluated matrix and the
        # Laurent matrix/map data that must stay algebraically consistent.
        _row_swap(evaluated, pivot_row, row)
        _row_swap_with_map(working, row_map, pivot_row, row, row_offset)
        _col_swap(evaluated, pivot_col, col, q)
        _col_swap_with_map(working, phi, psi, pivot_col, col, q)

        # Clear the pivot column, then the pivot row, using symplectic column
        # operations for the latter.
        for row in sorted(evaluated.nonzero_positions_in_column(pivot_col)):
            if pivot_row <= row < stop_row and row != pivot_row:
                _row_add(evaluated, row, pivot_row)
                _row_add_with_map(working, row_map, row, pivot_row, row_offset)
        for col in sorted(evaluated.nonzero_positions_in_row(pivot_row)):
            if pivot_col <= col < stop_col and col != pivot_col:
                _col_add(evaluated, col, pivot_col, q)
                _col_add_with_map(working, phi, psi, col, pivot_col, q)
        pivot_row += 1
        pivot_col += 1
        rank += 1
    return rank


def clear_css_mod_j(epsilon, *, num_x_checks, num_qubits):
    """Clear evaluated CSS blocks and retain ``(xi2, xi1, xi0)`` maps.

    Args:
        epsilon: Full coarse CSS excitation matrix over the Laurent ring.
        num_x_checks: Number of X-check rows.
        num_qubits: Number of data-qubit columns per CSS half.

    Returns:
        ``ClearingResult`` with the cleared working matrix and accumulated
        row/column maps.
    """
    working = copy(epsilon)
    # Pivot choices are made over the residue field R/(x+1,y+1), represented by
    # evaluating Laurent entries at x=y=1.
    evaluated = evaluate_at_identity(working).change_ring(GF(2))
    x_row_map = identity_matrix(GF(2), num_x_checks)
    z_row_map = identity_matrix(GF(2), epsilon.nrows() - num_x_checks)
    phi = identity_matrix(GF(2), num_qubits)
    psi = identity_matrix(GF(2), num_qubits)

    # First clear the X block on the left half, exposing X product-state rows.
    px = _clear_block(
        evaluated,
        working,
        x_row_map,
        phi,
        psi,
        q=num_qubits,
        row_offset=0,
        pivot_row=0,
        pivot_col=0,
        stop_row=num_x_checks,
        stop_col=num_qubits,
    )
    # Then clear the Z block on the right half, skipping columns already used by
    # X products so target blocks remain disjoint.
    pz = _clear_block(
        evaluated,
        working,
        z_row_map,
        phi,
        psi,
        q=num_qubits,
        row_offset=num_x_checks,
        pivot_row=num_x_checks,
        pivot_col=num_qubits + px,
        stop_row=epsilon.nrows(),
        stop_col=epsilon.ncols(),
    )
    # Lift accumulated GF(2) maps back into the Laurent ring for later chain-map
    # composition.
    x_row_map = x_row_map.change_ring(R)
    z_row_map = z_row_map.change_ring(R)
    phi = phi.change_ring(R)
    psi = psi.change_ring(R)
    row_zero = zero_matrix(R, num_x_checks, epsilon.nrows() - num_x_checks)
    row_map = block_matrix(R, [[x_row_map, row_zero], [row_zero.transpose(), z_row_map]])
    column_zero = zero_matrix(R, num_qubits, num_qubits)
    column_map = block_matrix(R, [[phi, column_zero], [column_zero, psi]])
    xi2 = antipode_matrix(z_row_map)
    xi1 = phi
    xi0 = x_row_map.inverse()
    return ClearingResult(
        epsilon,
        working,
        row_map,
        column_map,
        xi2,
        xi1,
        xi0,
        px,
        pz,
        num_x_checks,
        num_qubits,
    )
