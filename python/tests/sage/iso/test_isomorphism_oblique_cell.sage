"""Check oblique cell metadata used by matching-map projection."""

import _bootstrap  # noqa: F401

from isomorphism.coarse_graining import oblique_cell

cell = oblique_cell((1, 2), (2, 1))

assert cell.area == 3
assert len(cell.representatives) == 3

for slot, point in enumerate(cell.representatives):
    recovered_slot, cx, cy = cell.decompose(point[0], point[1])
    assert recovered_slot == slot
    assert (cx, cy) == (0, 0)

slot, cx, cy = cell.decompose(3, 3)
rx, ry = cell.representatives[slot]
assert cell.compose(slot, cx, cy) == (3, 3)
assert (3 - rx, 3 - ry) == (cx * cell.u[0] + cy * cell.v[0], cx * cell.u[1] + cy * cell.v[1])

large = 10**60
slot, cx, cy = cell.decompose(large, large)
assert cell.compose(slot, cx, cy) == (large, large)

slot, cx, cy = cell.decompose(-3, -3)
assert cell.compose(slot, cx, cy) == (-3, -3)

assert cell.contains_lattice_vector((3, 0))
assert cell.contains_lattice_vector((0, 3))
assert not cell.contains_lattice_vector((1, 0))

print("test_isomorphism_oblique_cell: ok")
