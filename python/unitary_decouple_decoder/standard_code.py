"""Decoder for standard product-state plus toric-sector CSS codes.

The decoupling pipeline produces standard ``H_X`` matrices whose rows are either
single-qubit product-state checks or one standard toric vertex check.  This
module decodes those rows directly without depending on any code-family-specific
decoder.
"""

from dataclasses import dataclass

from decoder_core import ToricMatchingDecoder
from decoder_core.finite_maps import finite_index, finite_matrix_from_laurent, zero_vector
from isomorphism import R, x, y


@dataclass(frozen=True)
class StandardCodeDecoder:
    """Decode the standard ``H_X`` produced by decoupling."""

    hx_standard: object
    shape: object
    weights: object = None

    def __post_init__(self):
        """Parse product-state and toric rows, then build finite matrices."""
        product_rows, toric_rows = _parse_standard_rows(self.hx_standard)
        toric = (
            ToricMatchingDecoder(self.shape, len(toric_rows), weights=self.weights)
            if toric_rows
            else None
        )
        object.__setattr__(self, "product_rows", tuple(product_rows))
        object.__setattr__(self, "toric_rows", tuple(toric_rows))
        object.__setattr__(self, "toric", toric)
        object.__setattr__(
            self,
            "_product_assignments",
            tuple(
                (
                    finite_index(row, site, self.shape),
                    finite_index(col, site, self.shape),
                )
                for row, col in product_rows
                for site in range(self.shape.size)
            ),
        )
        object.__setattr__(
            self,
            "_toric_syndrome_indices",
            tuple(
                tuple(finite_index(row, site, self.shape) for site in range(self.shape.size))
                for row, _col_h, _col_v in toric_rows
            ),
        )
        if toric is None:
            toric_correction_indices = ()
        else:
            toric_correction_indices = tuple(
                self._standard_correction_index_for_toric_edge(edge)
                for edge in range(toric.correction_size)
            )
        object.__setattr__(self, "_toric_correction_indices", toric_correction_indices)
        object.__setattr__(
            self,
            "hx_standard_finite",
            finite_matrix_from_laurent(self.hx_standard, self.shape),
        )

    @property
    def syndrome_size(self):
        """Return the finite standard syndrome length."""
        return self.hx_standard.nrows() * self.shape.size

    @property
    def correction_size(self):
        """Return the finite standard correction length."""
        return self.hx_standard.ncols() * self.shape.size

    def zero_syndrome(self):
        """Return the zero standard syndrome vector."""
        return zero_vector(self.hx_standard.nrows(), self.shape)

    def zero_correction(self):
        """Return the zero standard correction vector."""
        return zero_vector(self.hx_standard.ncols(), self.shape)

    def decode(self, syndrome):
        """Decode one standard-form syndrome.

        Args:
            syndrome: Standard GF(2) syndrome vector.

        Returns:
            vector: Standard GF(2) correction vector.
        """
        _validate_length(syndrome, self.syndrome_size, "syndrome")
        correction = self.zero_correction()
        # Product-state rows are solved independently by flipping the unique
        # data bit attached to every nonzero syndrome bit.
        for syndrome_index, correction_index in self._product_assignments:
            if syndrome[syndrome_index]:
                correction[correction_index] += 1

        if self.toric is None:
            return correction

        # Toric rows are copied into the stacked toric matching decoder using
        # the same finite site ordering as ``finite_matrix_from_laurent``.
        toric_syndrome = self.toric.zero_syndrome()
        for stack, syndrome_indices in enumerate(self._toric_syndrome_indices):
            vertex_offset = stack * self.toric.vertices_per_stack
            for site, syndrome_index in enumerate(syndrome_indices):
                if syndrome[syndrome_index]:
                    toric_syndrome[vertex_offset + site] = 1

        toric_correction = self.toric.decode(toric_syndrome)
        # Map horizontal and vertical toric edges back to standard-code columns.
        for edge, value in enumerate(toric_correction):
            if not value:
                continue
            correction[self._toric_correction_indices[edge]] += 1
        return correction

    def _standard_correction_index_for_toric_edge(self, edge):
        """Return the standard correction index for one toric edge.

        Args:
            edge: Stacked toric edge index.

        Returns:
            int: Component-major standard correction index.
        """
        stack, orientation, px, py = self.toric.edge_key(edge)
        _row, col_h, col_v = self.toric_rows[stack]
        col = col_h if orientation == "h" else col_v
        return finite_index(col, self.shape.index(px, py), self.shape)


def _parse_standard_rows(hx_standard):
    """Return product and toric row descriptors for a standard matrix.

    Args:
        hx_standard: Standard X-check matrix produced by decoupling.

    Returns:
        tuple: Product-row and toric-row descriptor lists.
    """
    product_rows = []
    toric_rows = []
    one = R.one()
    horizontal = one + x
    vertical = one + y
    for row in range(hx_standard.nrows()):
        entries = [
            (col, R(hx_standard[row, col]))
            for col in range(hx_standard.ncols())
            if hx_standard[row, col] != 0
        ]
        # A product-state row fixes one data bit at every finite site.
        if len(entries) == 1 and entries[0][1] == one:
            product_rows.append((row, entries[0][0]))
            continue
        # A toric row has one horizontal and one vertical incident edge column.
        if len(entries) == 2:
            by_value = {value: col for col, value in entries}
            if horizontal in by_value and vertical in by_value:
                toric_rows.append((row, by_value[horizontal], by_value[vertical]))
                continue
        raise ValueError("Standard H_X row must be a product row or toric row.")
    return product_rows, toric_rows


def _validate_length(values, expected, name):
    """Require an input vector to have the expected finite length.

    Args:
        values: Input sequence to normalize or validate.
        expected: Expected vector length used for validation.
        name: Human-readable field name used in validation errors.
    """
    if len(values) != expected:
        raise ValueError(f"{name} length must be {expected}.")
