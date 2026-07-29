import _bootstrap  # noqa: F401

from sage.all import Matrix

from isomorphism import R, x, y
import isomorphism.coarse_graining as coarse_graining_module
from isomorphism.coarse_graining import oblique_coarse_grain, straight_coarse_grain

mat = Matrix(R, [[1 + x, y + x**-1]])
straight = straight_coarse_grain(mat, 2)
assert straight.nrows() == 4
assert straight.ncols() == 8

identity_cell = oblique_coarse_grain(mat, (1, 0), (0, 1))
assert identity_cell == mat

paper_cell = oblique_coarse_grain(Matrix(R, [[1 + x + y]]), (2, 0), (1, 1))
assert paper_cell.nrows() == 2
assert paper_cell.ncols() == 2

original_oblique_translation_blocks = coarse_graining_module._oblique_translation_blocks


def reject_translation_blocks(_u, _v):
    raise AssertionError("oblique_coarse_grain should build monomial blocks directly")


coarse_graining_module._oblique_translation_blocks = reject_translation_blocks
try:
    direct_cell = coarse_graining_module.oblique_coarse_grain(
        Matrix(R, [[1 + x + y]]),
        (2, 0),
        (1, 1),
    )
finally:
    coarse_graining_module._oblique_translation_blocks = original_oblique_translation_blocks

assert direct_cell == paper_cell

print("test_isomorphism_coarse_graining: ok")
