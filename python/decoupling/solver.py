"""Verified Laurent-polynomial matrix equations."""

from sage.all import Matrix, block_matrix, identity_matrix, zero_matrix

from .algebra.singular_backend import (
    SingularBackendError,
    matrix_laurent_to_polynomial,
    matrix_polynomial_to_laurent,
    membership_lift,
)
from .rings import P, PRESENTATION, R


class LaurentSolveError(ValueError):
    """Raised when an exact Laurent matrix equation has no witness."""


def add_inverse_relation_columns(mat):
    """Lift a Laurent matrix and add inverse relations per module row.

    Args:
        mat: Laurent matrix whose columns generate a module.

    Returns:
        Polynomial matrix over ``P`` with extra columns enforcing
        ``a*c=1`` and ``b*d=1`` for every row of the module.
    """
    lifted = matrix_laurent_to_polynomial(mat)
    nrows = mat.nrows()
    relation = zero_matrix(P, nrows, 2 * nrows)
    x_relation, y_relation = PRESENTATION.inverse_relations()
    for row in range(nrows):
        # Singular works over an ordinary polynomial ring, so Laurent inverses
        # are represented by row-local quotient relation columns.
        relation[row, row] = x_relation
        relation[row, row + nrows] = y_relation
    return block_matrix(P, [[lifted, relation]])


def solve_right(A, B, *, verify=True):
    """Return ``X`` with ``A * X == B`` when an exact Laurent solution exists.

    Args:
        A: Left Laurent matrix.
        B: Target Laurent matrix.
        verify: If true, multiply the candidate witness back over the Laurent
            ring.

    Returns:
        Laurent matrix ``X`` satisfying ``A * X == B``.

    Raises:
        LaurentSolveError: If the target is outside the generated Laurent
            module or verification fails.
    """
    generators = add_inverse_relation_columns(A)
    targets = matrix_laurent_to_polynomial(B)
    try:
        # Singular returns polynomial lift coefficients for the augmented
        # generator matrix. Only the first A.ncols rows correspond to the
        # original Laurent matrix columns.
        witness = membership_lift(generators, targets)
    except SingularBackendError as exc:
        raise LaurentSolveError("Right-hand side is not in the generated module.") from exc
    candidate = matrix_polynomial_to_laurent(witness[: A.ncols(), :])
    if verify and A * candidate != B:
        raise LaurentSolveError("Singular witness failed Laurent substitution verification.")
    return candidate


def solve_left(A, B, *, verify=True):
    """Return ``X`` with ``X * A == B`` when an exact Laurent solution exists.

    Args:
        A: Right Laurent matrix.
        B: Target Laurent matrix.
        verify: If true, verify the result by multiplication.

    Returns:
        Laurent matrix ``X`` satisfying ``X * A == B``.

    Raises:
        LaurentSolveError: If no exact left witness exists.
    """
    candidate = solve_right(A.transpose(), B.transpose(), verify=verify).transpose()
    if verify and candidate * A != B:
        raise LaurentSolveError("Left solve failed Laurent substitution verification.")
    return candidate


def _native_inverse(A, *, verify):
    """Return Sage's direct inverse when it stays over the Laurent ring.

    Args:
        A: Square Laurent matrix.
        verify: If true, check both inverse products.

    Returns:
        Laurent inverse matrix, or ``None`` if Sage leaves the Laurent ring.
    """
    try:
        candidate = A.inverse()
    except Exception:  # noqa: BLE001 - any native inversion failure should fall back.
        return None
    if candidate.base_ring() != R:
        try:
            candidate = candidate.change_ring(R)
        except Exception:  # noqa: BLE001 - fraction-field inverses are not Laurent maps.
            return None
    if verify:
        identity = identity_matrix(R, A.nrows())
        if A * candidate != identity or candidate * A != identity:
            return None
    return candidate


def find_inverse(A, *, verify=True):
    """Return the two-sided inverse of a square Laurent matrix.

    Args:
        A: Square Laurent matrix.
        verify: If true, verify the inverse equations.

    Returns:
        Laurent matrix inverse.

    Raises:
        ValueError: If ``A`` is not square.
        LaurentSolveError: If no Laurent inverse can be certified.
    """
    if A.nrows() != A.ncols():
        raise ValueError("Matrix must be square to compute an inverse.")
    identity = identity_matrix(R, A.nrows())
    solve_error = None
    try:
        # The module-membership solve is the main path because it can certify
        # Laurent inverses that Sage's native inverse routes through fractions.
        candidate = solve_right(A, identity, verify=verify)
    except LaurentSolveError as exc:
        solve_error = exc
    else:
        if verify and candidate * A != identity:
            raise LaurentSolveError("Solved right inverse is not a two-sided inverse.")
        return candidate
    # Fall back to Sage only when it can stay inside the Laurent ring.
    candidate = _native_inverse(A, verify=verify)
    if candidate is not None:
        return candidate
    raise solve_error
