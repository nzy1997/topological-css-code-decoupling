"""Coarse-graining maps for Laurent CSS matrices."""

from dataclasses import dataclass, field

from sage.all import QQ, block_matrix, identity_matrix, matrix, vector, zero_matrix

from .css import dagger_matrix
from .rings import R, x, y


def _exact_floor(value):
    """Return the exact integer floor of a Sage rational-like value.

    Args:
        value: Sage rational-like value or Python numeric value.

    Returns:
        Integer floor computed without floating-point conversion.
    """
    if hasattr(value, "floor"):
        return int(value.floor())
    return int(QQ(value).floor())


@dataclass(frozen=True)
class ObliqueCell:
    """Metadata for decomposing points into an oblique coarse cell.

    Attributes:
        u: First integer basis vector of the coarse cell.
        v: Second integer basis vector of the coarse cell.
        basis_matrix: Rational matrix with columns ``u`` and ``v``.
        representatives: Original-lattice representatives inside the
            half-open fundamental parallelogram.
        _representative_slots: Lookup from representative point to block index.
    """

    u: tuple
    v: tuple
    basis_matrix: object
    representatives: tuple
    _representative_slots: dict = field(repr=False)

    @property
    def area(self):
        """Return the number of lattice sites in one oblique cell (the index Lambda).

        Returns:
            Cell index, equal to the number of representatives.
        """
        return len(self.representatives)

    def decompose(self, px, py):
        """Return representative slot and oblique-cell coordinates for a point.

        Args:
            px: Original-lattice x coordinate.
            py: Original-lattice y coordinate.

        Returns:
            Triple ``(slot, cell_x, cell_y)`` such that the point equals
            ``representatives[slot] + cell_x*u + cell_y*v``.
        """
        coeff_u, coeff_v = self.basis_matrix.solve_right(vector(QQ, [px, py]))
        cell_x, cell_y = _exact_floor(coeff_u), _exact_floor(coeff_v)
        # Remove the integer cell part; the remainder lies in the fundamental
        # parallelogram and identifies one representative slot.
        rep = self.basis_matrix * vector(QQ, [coeff_u - cell_x, coeff_v - cell_y])
        representative = (int(rep[0]), int(rep[1]))
        return self._representative_slots[representative], cell_x, cell_y

    def compose(self, slot, cell_x, cell_y):
        """Return the lattice point for a representative slot and cell coordinates.

        Args:
            slot: Representative index.
            cell_x: Coarse coordinate along ``u``.
            cell_y: Coarse coordinate along ``v``.

        Returns:
            Original-lattice integer point.
        """
        rep_x, rep_y = self.representatives[slot]
        return (
            rep_x + cell_x * self.u[0] + cell_y * self.v[0],
            rep_y + cell_x * self.u[1] + cell_y * self.v[1],
        )

    def contains_lattice_vector(self, point):
        """Return whether a vector has integer coordinates in the oblique basis.

        Args:
            point: Two-integer vector in original lattice coordinates.

        Returns:
            True when ``point`` belongs to the coarse lattice spanned by
            ``u`` and ``v``.
        """
        coeff_u, coeff_v = self.basis_matrix.solve_right(vector(QQ, point))
        return coeff_u.denominator() == 1 and coeff_v.denominator() == 1


@dataclass
class PowerCache:
    """Cache powers of one translation matrix, including inverse powers.

    Attributes:
        block: Translation matrix for one positive unit step.
        identity: Identity matrix with the same shape as ``block``.
        powers: Cache from integer powers to matrix powers.
    """

    block: object
    identity: object
    powers: dict = field(default_factory=dict)

    def __post_init__(self):
        """Seed the power cache with the identity matrix."""
        self.powers[0] = self.identity

    def __call__(self, power):
        """Return one cached positive, zero, or inverse matrix power.

        Args:
            power: Integer translation exponent.

        Returns:
            Matrix implementing the requested finite translation.
        """
        if power not in self.powers:
            if power > 0:
                self.powers[power] = self.block**power
            else:
                # The inverse finite translation is represented by the dagger
                # block because Laurent monomial inversion reverses direction.
                self.powers[power] = dagger_matrix(self.block) ** (-power)
        return self.powers[power]


def _cycle_block(var, period):
    """Build one translation block for a square period.

    Args:
        var: Boundary monomial inserted when the translation wraps.
        period: Square cell period.

    Returns:
        ``period x period`` Laurent matrix for one cyclic shift.
    """
    block = zero_matrix(R, period, period)
    block[0, period - 1] = var
    for index in range(period - 1):
        block[index + 1, index] = 1
    return block


def coarse_grain_to_square_superlattice(mat, period):
    """Coarse-grain by the square cell spanned by ``(period, 0)``, ``(0, period)``.

    Args:
        mat: Laurent matrix to rewrite.
        period: Positive square coarse-cell period.

    Returns:
        Block Laurent matrix with each original entry expanded to a
        ``period^2 x period^2`` translation block.

    Raises:
        ValueError: If ``period`` is not positive.
    """
    if period < 1:
        raise ValueError("period must be positive.")
    x_power = PowerCache(_cycle_block(x, period), identity_matrix(R, period))
    y_power = PowerCache(_cycle_block(y, period), identity_matrix(R, period))

    def poly_block(poly):
        """Replace one Laurent entry by a square-cell translation block.

        Args:
            poly: Laurent polynomial entry.

        Returns:
            Matrix block implementing all monomial translations in the square
            cell basis.
        """
        block = zero_matrix(R, period**2, period**2)
        for (px, py), _coeff in poly.dict().items():
            # Square cells factor into independent x/y cyclic translations.
            block += x_power(px).tensor_product(y_power(py))
        return block

    return block_matrix(
        R,
        mat.nrows(),
        mat.ncols(),
        [[poly_block(mat[row, col]) for col in range(mat.ncols())] for row in range(mat.nrows())],
    )


def _cell_basis(u, v):
    """Return integer representatives for one oblique unit cell.

    Args:
        u: First integer cell vector.
        v: Second integer cell vector.

    Returns:
        Pair ``(basis_matrix, representative_slots)``.

    Raises:
        ValueError: If ``u`` and ``v`` are linearly dependent.
    """
    basis_matrix = matrix(QQ, [[u[0], v[0]], [u[1], v[1]]])
    if basis_matrix.det() == 0:
        raise ValueError("u and v must span a two-dimensional cell.")
    corners = [(0, 0), u, v, (u[0] + v[0], u[1] + v[1])]
    representatives = {}
    for point_x in range(
        min(point[0] for point in corners),
        max(point[0] for point in corners) + 1,
    ):
        for point_y in range(
            min(point[1] for point in corners),
            max(point[1] for point in corners) + 1,
        ):
            coeff_u, coeff_v = basis_matrix.solve_right(vector(QQ, [point_x, point_y]))
            # Half-open coordinates avoid double-counting boundary points of
            # the fundamental parallelogram.
            if 0 <= coeff_u < 1 and 0 <= coeff_v < 1:
                representatives[(point_x, point_y)] = len(representatives)
    return basis_matrix, representatives


def oblique_cell(u, v):
    """Return metadata for the oblique cell spanned by independent integer vectors.

    Args:
        u: First integer cell vector.
        v: Second integer cell vector.

    Returns:
        ``ObliqueCell`` instance with representatives and lookup table.

    Raises:
        ValueError: If ``u`` and ``v`` do not span a two-dimensional cell.
    """
    u = (int(u[0]), int(u[1]))
    v = (int(v[0]), int(v[1]))
    basis_matrix, representative_slots = _cell_basis(u, v)
    representatives = [None] * len(representative_slots)
    for representative, slot in representative_slots.items():
        representatives[slot] = representative
    return ObliqueCell(
        u=u,
        v=v,
        basis_matrix=basis_matrix,
        representatives=tuple(representatives),
        _representative_slots=representative_slots,
    )


def _translate_basis_entry(basis_matrix, representatives, point, translation):
    """Translate one representative and reduce it into its oblique cell.

    Args:
        basis_matrix: Rational matrix with columns ``u`` and ``v``.
        representatives: Lookup from representative point to row index.
        point: Representative point as a vector.
        translation: Original-lattice translation vector.

    Returns:
        Pair ``(row, monomial)`` for the translated representative.
    """
    coeff_u, coeff_v = basis_matrix.solve_right(point + translation)
    cell_u, cell_v = _exact_floor(coeff_u), _exact_floor(coeff_v)
    rep = basis_matrix * vector(QQ, [coeff_u - cell_u, coeff_v - cell_v])
    representative = (int(rep[0]), int(rep[1]))
    return representatives[representative], x**cell_u * y**cell_v


def _oblique_translation_blocks(u, v):
    """Return x/y blocks and identity for an oblique cell.

    Args:
        u: First integer cell vector.
        v: Second integer cell vector.

    Returns:
        Triple ``(tx, ty, identity)`` for unit translations in the oblique cell.
    """
    basis_matrix, representatives = _cell_basis(u, v)
    size = len(representatives)
    tx = zero_matrix(R, size, size)
    ty = zero_matrix(R, size, size)
    for representative, col in representatives.items():
        point = vector(QQ, representative)
        # Each representative contributes one nonzero entry per unit
        # translation because translation permutes representatives up to a
        # coarse-cell monomial.
        row_x, elem_x = _translate_basis_entry(
            basis_matrix,
            representatives,
            point,
            vector(QQ, [1, 0]),
        )
        row_y, elem_y = _translate_basis_entry(
            basis_matrix,
            representatives,
            point,
            vector(QQ, [0, 1]),
        )
        tx[row_x, col] = elem_x
        ty[row_y, col] = elem_y
    return tx, ty, identity_matrix(R, size)


def _oblique_monomial_block(cell, px, py):
    """Return the direct block for translating by one Laurent monomial.

    Args:
        cell: Oblique cell metadata.
        px: Original-lattice x exponent.
        py: Original-lattice y exponent.

    Returns:
        ``cell.area x cell.area`` Laurent translation block.

    Examples:
        Let ``u = (2, 0)`` and ``v = (1, 1)``. Then the oblique cell contains the representatives:
            slot 0: (0, 0)
            slot 1: (1, 0)
        Translating by the monomial ``x`` sends these representatives to:
            (0, 0) + (1, 0) = (1, 0) = slot 1
            (1, 0) + (1, 0) = (2, 0) = slot 0 with an additional translation by ``u``.
        Thus the block for ``x`` is:
            [[0, x],
             [1, 0]] 
        (Note the block entry ``x`` now refers to the coarse-grained translation by ``u`` instead of the original-lattice translation by ``x``.)
    """
    block = zero_matrix(R, cell.area, cell.area)
    for col, (rep_x, rep_y) in enumerate(cell.representatives):
        # Translate the representative, reduce back to the representative list,
        # and record the coarse-cell displacement as a Laurent monomial.
        row, cell_x, cell_y = cell.decompose(rep_x + px, rep_y + py)
        block[row, col] = x**cell_x * y**cell_y
    return block


def coarse_grain_to_superlattice(mat, u, v):
    """Coarse-grain by the cell spanned by independent integer vectors.

    Args:
        mat: Laurent matrix to rewrite.
        u: First integer cell vector.
        v: Second integer cell vector.

    Returns:
        Block Laurent matrix with each original entry expanded to an
        ``area x area`` oblique translation block.

    Raises:
        ValueError: If ``u`` and ``v`` are linearly dependent.
    """
    cell = oblique_cell(u, v)
    monomial_blocks = {}

    def poly_block(poly):
        """Replace one Laurent entry by an oblique-cell translation block.

        Args:
            poly: Laurent polynomial entry.

        Returns:
            Matrix block implementing all monomial translations in the oblique
            cell basis.
        """
        block = zero_matrix(R, cell.area, cell.area)
        for (px, py), coeff in poly.dict().items():
            key = (px, py)
            if key not in monomial_blocks:
                # Reuse identical monomial translations across all matrix
                # entries; this is important for large BB cells.
                monomial_blocks[key] = _oblique_monomial_block(cell, px, py)
            block += coeff * monomial_blocks[key]
        return block

    return block_matrix(
        R,
        mat.nrows(),
        mat.ncols(),
        [[poly_block(mat[row, col]) for col in range(mat.ncols())] for row in range(mat.nrows())],
    )
