"""Finite periodic projections of Laurent matrices."""

from sage.all import GF, Matrix, vector

from isomorphism import R


def finite_index(component, site_index, shape):
    """Return component-major finite vector index.

    Args:
        component: Component block index in component-major finite vectors.
        site_index: Linear site index inside the finite torus.
        shape: Finite torus shape used for periodic projection.

    Returns:
        int: Component-major vector index.
    """
    return int(component) * shape.size + int(site_index)


def zero_vector(component_count, shape):
    """Return a zero finite vector with one block per component.

    Args:
        component_count: Number of Laurent components represented on the finite torus.
        shape: Finite torus shape used for periodic projection.

    Returns:
        vector: Zero GF(2) vector in component-major order.
    """
    return vector(GF(2), int(component_count) * shape.size)


def finite_matrix_from_laurent(mat, shape):
    """Project a Laurent matrix to a periodic finite GF(2) matrix.

    Args:
        mat: Laurent matrix to project.
        shape: Finite torus shape used for periodic projection.

    Returns:
        Matrix: Periodic finite GF(2) matrix.
    """
    # Use component-major ordering: every Laurent row/column expands into one
    # full torus block.  A monomial x^dx y^dy sends the source site (px, py) to
    # the periodically wrapped target site (px + dx, py + dy).
    out = Matrix(GF(2), mat.nrows() * shape.size, mat.ncols() * shape.size, 0)
    for row in range(mat.nrows()):
        for col in range(mat.ncols()):
            for (dx, dy), coeff in R(mat[row, col]).dict().items():
                if int(coeff) % 2:
                    for px, py in shape.sites():
                        target_site = shape.index(px + int(dx), py + int(dy))
                        source_site = shape.index(px, py)
                        out[
                            finite_index(row, target_site, shape),
                            finite_index(col, source_site, shape),
                        ] += 1
    return out
