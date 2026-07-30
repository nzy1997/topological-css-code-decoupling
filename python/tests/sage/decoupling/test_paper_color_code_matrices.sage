"""Verify the literal inverse chain maps printed in the supplement."""

import _bootstrap  # noqa: F401

from sage.all import Matrix

from decoupling import R, dagger_matrix, x, y


def verify_inverse_chain_certificate(
    h_x,
    h_z,
    h_x_tilde,
    h_z_tilde,
    psi_0_inverse,
    psi_1_inverse,
    psi_2_inverse,
):
    """Check the two inverse chain equations for one printed certificate."""
    assert h_x * psi_1_inverse == psi_0_inverse * h_x_tilde
    assert (
        dagger_matrix(h_z) * psi_2_inverse
        == psi_1_inverse * dagger_matrix(h_z_tilde)
    )


h_x_666 = Matrix(
    R,
    [
        [1, 1, x**-1, 1, 1, y**-1],
        [y, 1, 1, x, 1, 1],
        [x * y, y, 1, x * y, x, 1],
    ],
)
h_z_666 = Matrix(
    R,
    [
        [1, x**-1, x**-1 * y**-1, 1, y**-1, x**-1 * y**-1],
        [1, 1, x**-1, 1, 1, y**-1],
        [y, 1, 1, x, 1, 1],
    ],
)
h_x_tilde_666 = Matrix(
    R,
    [
        [1, 0, 0, 0, 0, 0],
        [0, 0, x + 1, y + 1, 0, 0],
        [0, 0, 0, 0, x + 1, y + 1],
    ],
)
h_z_tilde_666 = Matrix(
    R,
    [
        [0, 1, 0, 0, 0, 0],
        [0, 0, 1 + y**-1, 1 + x**-1, 0, 0],
        [0, 0, 0, 0, 1 + y**-1, 1 + x**-1],
    ],
)
psi_0_inverse_666 = Matrix(R, [[1, 0, 0], [1, 1, 0], [1, 0, 1]])
psi_1_inverse_666 = Matrix(
    R,
    [
        [1 + y**-1, 1, 1, 0, 1 + y**-1, 1 + y**-1],
        [
            x * y**-1 + 1 + y**-1,
            x,
            1,
            1,
            x * y**-1 + 1 + y**-1,
            x * y**-1 + y**-1,
        ],
        [x * y**-1, x * y, 0, 0, x * y**-1, x * y**-1],
        [1, 1, 1, 0, 1, 1],
        [1 + y**-1, y, 0, 0, y**-1, 1 + y**-1],
        [x + y, x * y, y, y, x + y, x + y],
    ],
)
psi_2_inverse_666 = Matrix(
    R,
    [
        [1, 1, 1 + x * y**-1],
        [0, y, x + x * y**-1 + y],
        [0, 0, x],
    ],
)

verify_inverse_chain_certificate(
    h_x_666,
    h_z_666,
    h_x_tilde_666,
    h_z_tilde_666,
    psi_0_inverse_666,
    psi_1_inverse_666,
    psi_2_inverse_666,
)
assert psi_0_inverse_666.det() == 1
assert psi_1_inverse_666.det().is_unit()
assert psi_2_inverse_666.det() == x * y


h_x_488 = Matrix(
    R,
    [
        [1, 1, 1, 1, 1, x**-1 * y, 1, x**-1 * y],
        [x, 1, x, 1, y, 1, y, 1],
        [y, 0, 0, 1, y, 0, 0, x**-1 * y],
        [0, y, x, 0, 0, y, y, 0],
    ],
)
h_z_488 = h_x_488
h_x_tilde_488 = Matrix(
    R,
    [
        [1, 0, 0, 0, 0, 0, 0, 0],
        [0, 1, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, x + 1, y + 1, 0, 0],
        [0, 0, 0, 0, 0, 0, x + 1, y + 1],
    ],
)
h_z_tilde_488 = Matrix(
    R,
    [
        [0, 0, 1, 0, 0, 0, 0, 0],
        [0, 0, 0, 1, 0, 0, 0, 0],
        [0, 0, 0, 0, 1 + y**-1, 1 + x**-1, 0, 0],
        [0, 0, 0, 0, 0, 0, 1 + y**-1, 1 + x**-1],
    ],
)
psi_0_inverse_488 = Matrix(
    R,
    [
        [1, 0, 0, 0],
        [1, 0, 1, 0],
        [1, 1, 0, 0],
        [0, 1, 0, 1],
    ],
)
psi_2_inverse_488 = Matrix(
    R,
    [
        [1, 0, y, y],
        [0, 0, 0, x],
        [0, 1, y, 0],
        [0, 0, x * y, 0],
    ],
)
psi_1_inverse_488 = Matrix(
    R,
    [
        [0, 0, 1, y**-1, 1, 0, 1, 0],
        [0, 0, 1, 0, 1, 1, 1, 1],
        [0, 0, 1, 0, 0, 0, 1, 0],
        [1, 0, 1, 1, 0, 0, 1, 1],
        [0, y**-1, 1, y**-1, 1, 0, 1 + y**-1, y**-1],
        [0, 0, x * y**-1, 0, 0, 0, 0, 0],
        [0, y**-1, 1, 0, 1, 1, 1 + y**-1, y**-1],
        [0, 0, x * y**-1, x * y**-1, 0, 0, 0, 0],
    ],
)

verify_inverse_chain_certificate(
    h_x_488,
    h_z_488,
    h_x_tilde_488,
    h_z_tilde_488,
    psi_0_inverse_488,
    psi_1_inverse_488,
    psi_2_inverse_488,
)
assert psi_0_inverse_488.det() == 1
assert psi_1_inverse_488.det() == x**2 * y**-3
assert psi_2_inverse_488.det() == x**2 * y

print("test_paper_color_code_matrices: ok")
