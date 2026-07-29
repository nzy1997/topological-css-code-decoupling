"""Check Laurent-map finite projection helpers."""

import _bootstrap  # noqa: F401

from sage.all import GF, Matrix, vector

from decoder_core import TorusShape
from decoder_core.finite_maps import finite_matrix_from_laurent
from isomorphism import R, x, y

shape = TorusShape(2, 2)
mat = Matrix(R, [[1 + x, y]])
finite = finite_matrix_from_laurent(mat, shape)

assert finite.nrows() == shape.size
assert finite.ncols() == 2 * shape.size

source = vector(GF(2), 2 * shape.size)
source[shape.index(0, 0)] = 1
out = finite * source
expected = vector(GF(2), shape.size)
expected[shape.index(0, 0)] = 1
expected[shape.index(1, 0)] = 1
assert out == expected

print("test_decoder_core_finite_maps: ok")
