"""CSS matrix conventions over F_2[x^{+-1}, y^{+-1}]."""

from sage.all import Matrix, block_matrix, identity_matrix, zero_matrix

from .rings import R, x, y


def antipode_poly(poly):
    """Apply ``x -> x^-1`` and ``y -> y^-1`` to one Laurent polynomial.

    Args:
        poly: Laurent polynomial over ``R``.

    Returns:
        Laurent polynomial after applying the antipode involution.
    """
    return poly(x=x**-1, y=y**-1)


def dagger_matrix(mat):
    """Transpose and apply the Laurent antipode entrywise.

    Args:
        mat: Sage matrix over the Laurent ring.

    Returns:
        Matrix ``mat^dagger`` with shape ``mat.ncols() x mat.nrows()``.
    """
    result = zero_matrix(mat.base_ring(), mat.ncols(), mat.nrows())
    for row, col in mat.nonzero_positions():
        # CSS adjoints combine transpose with reversing lattice translations.
        result[col, row] = antipode_poly(mat[row, col])
    return result


def evaluate_poly_at_identity(poly):
    """Evaluate one Laurent polynomial at ``x=y=1``.

    Args:
        poly: Laurent polynomial over ``R``.

    Returns:
        Element of ``R`` equal to the evaluation at the identity point.
    """
    return poly(x=R.one(), y=R.one())


def evaluate_at_identity(mat):
    """Evaluate a Laurent matrix at ``x=y=1``.

    Args:
        mat: Sage matrix over the Laurent ring.

    Returns:
        Matrix of the same shape whose entries are evaluated at the identity.
    """
    result = zero_matrix(mat.base_ring(), mat.nrows(), mat.ncols())
    for row, col in mat.nonzero_positions():
        result[row, col] = evaluate_poly_at_identity(mat[row, col])
    return result


def build_two_generator_css_excitation_map(f, g):
    """Build the two-generator CSS excitation matrix.

    Args:
        f: First Laurent generator.
        g: Second Laurent generator.

    Returns:
        ``2 x 4`` CSS excitation matrix with ``H_X = [f, g]`` and
        ``H_Z = [g^dagger, f^dagger]`` in the repository block convention.
    """
    return Matrix(
        R,
        [
            [f, g, 0, 0],
            [0, 0, antipode_poly(g), antipode_poly(f)],
        ],
    )


def check_commutation(epsilon, num_qubits):
    """Check the CSS commutation equation.

    Args:
        epsilon: Full CSS excitation matrix.
        num_qubits: Number of qubit columns in one CSS half.

    Returns:
        True when ``epsilon * Lambda * epsilon^dagger`` is zero.
    """
    zero = R.zero()
    identity = identity_matrix(R, num_qubits)
    # Lambda swaps the X and Z column halves in the symplectic CSS pairing.
    symplectic = block_matrix(R, [[zero, identity], [identity, zero]])
    return (epsilon * symplectic * dagger_matrix(epsilon)).is_zero()


def split_css_blocks(epsilon, num_x_checks, num_qubits):
    """Return ``(HZ, HX)`` from a full block CSS excitation matrix.

    Args:
        epsilon: Full CSS excitation matrix.
        num_x_checks: Number of rows in the X-check block.
        num_qubits: Number of columns in one qubit half.

    Returns:
        Pair ``(HZ, HX)`` using the repository convention.
    """
    return (
        epsilon[num_x_checks:, num_qubits:],
        epsilon[:num_x_checks, :num_qubits],
    )


def build_standard_blocks(full_matrix, num_x_checks, num_qubits, px, pz):
    """Return the standard ``(H_Z, H_X)`` blocks after CSS clearing.

    Args:
        full_matrix: Cleared full CSS matrix. Only its shape is used.
        num_x_checks: Number of X-check rows.
        num_qubits: Number of data-qubit columns per CSS half.
        px: Product-state rank exposed in the X block.
        pz: Product-state rank exposed in the Z block.

    Returns:
        Pair of standard Laurent blocks. Product rows are identities; the
        remaining rows are toric checks.
    """
    h_x_tilde = zero_matrix(R, num_x_checks, num_qubits)
    # X product-state rows each pin one data column.
    for index in range(px):
        h_x_tilde[index, index] = 1
    # Non-product X rows become standard toric vertex checks.
    for row in range(px, num_x_checks):
        col = px + pz + 2 * (row - px)
        h_x_tilde[row, col] = 1 + x
        h_x_tilde[row, col + 1] = 1 + y

    num_z_checks = full_matrix.nrows() - num_x_checks
    h_z_tilde = zero_matrix(R, num_z_checks, num_qubits)
    # Z product-state rows start after the X product columns.
    for index in range(pz):
        h_z_tilde[index, px + index] = 1
    # Z toric rows use daggered translations relative to the standard X rows.
    for row in range(pz, num_z_checks):
        col = px + pz + 2 * (row - pz)
        h_z_tilde[row, col] = 1 + y**-1
        h_z_tilde[row, col + 1] = 1 + x**-1
    return h_z_tilde, h_x_tilde
