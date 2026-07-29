import _bootstrap  # noqa: F401

from isomorphism import x, y
from isomorphism.periods import (
    choose_smallest_oblique_cell,
    periods_from_generators,
    superlattice_from_generators,
)

trivial = periods_from_generators(1 + x, 1 + y, max_period=1)
assert trivial.square_period == 1
assert (1, 0) in trivial.vectors
assert (0, 1) in trivial.vectors

f_666 = 1 + x + x * y
g_666 = 1 + y + x * y
periods_666 = periods_from_generators(f_666, g_666, max_period=3)
assert periods_666.square_period == 3
assert (3, 0) in periods_666.vectors
assert (2, 1) in periods_666.vectors
assert abs(
    choose_smallest_oblique_cell(periods_666.vectors)[0][0]
    * choose_smallest_oblique_cell(periods_666.vectors)[1][1]
    - choose_smallest_oblique_cell(periods_666.vectors)[0][1]
    * choose_smallest_oblique_cell(periods_666.vectors)[1][0]
) == 3
superlattice_666 = superlattice_from_generators(f_666, g_666, max_period=3)
assert superlattice_666.alpha == 3
assert superlattice_666.gamma == 2
assert superlattice_666.delta == 1
assert superlattice_666.cell == ((3, 0), (2, 1))

bb = periods_from_generators(1 + x + x**-1 * y, 1 + y + x * y, max_period=7)
assert bb.square_period == 7

bb_186 = periods_from_generators(
    1 + x + x**-1 * y**3,
    1 + y + x * y**3,
    max_period=186,
)
assert bb_186.square_period == 186
assert choose_smallest_oblique_cell(bb_186.vectors) == ((14, 32), (16, 10))

bb_217_superlattice = superlattice_from_generators(
    1 + x + x**-1 * y**2,
    1 + y + x**-2 * y**-1,
    max_period=217,
)
assert bb_217_superlattice.alpha == 217
assert bb_217_superlattice.gamma == 67
assert bb_217_superlattice.delta == 1
bb_217 = periods_from_generators(
    1 + x + x**-1 * y**2,
    1 + y + x**-2 * y**-1,
    max_period=217,
)
assert bb_217.square_period == 217
assert len(bb_217.vectors) == 219
assert (217, 0) in bb_217.vectors
assert (67, 1) in bb_217.vectors
assert (0, 217) in bb_217.vectors

try:
    periods_from_generators(1 + x + x * y, 1 + y + x * y, max_period=2)
except ValueError as exc:
    assert "max_period" in str(exc)
else:
    raise AssertionError("Expected bounded ideal period search to fail.")

print("test_isomorphism_periods_ideal: ok")
