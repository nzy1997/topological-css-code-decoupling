"""Check reusable quotient coordinates for Laurent cokernels."""

import _bootstrap  # noqa: F401

from sage.all import Matrix

from isomorphism import R, x, y
from isomorphism.algebra.quotient import quotient_coordinate_system


def assert_value_error(fn, text):
    try:
        fn()
    except ValueError as exc:
        assert text in str(exc)
    else:
        raise AssertionError("Expected ValueError.")


system = quotient_coordinate_system(Matrix(R, [[1 + x + x * y, 1 + y + x * y]]))

assert system.dimension == 2
assert len(system.laurent_basis) == 2
assert system.coordinates([0]).is_zero()
assert system.coordinates([1 + x + x * y]).is_zero()
assert system.coordinates([1 + y + x * y]).is_zero()

origin = system.coordinates([1])
assert origin.length() == 2
assert system.coordinates([x**3]) == origin
assert system.coordinates([y**3]) == origin

rank_two_system = quotient_coordinate_system(
    Matrix(
        R,
        [
            [1 + x + x * y, 1 + y + x * y, 0, 0],
            [0, 0, 1 + x + x * y, 1 + y + x * y],
        ],
    )
)

assert rank_two_system.dimension == 4
assert rank_two_system.coordinates([1, 0]) != rank_two_system.coordinates([0, 1])
assert all(isinstance(entry, tuple) for entry in rank_two_system.laurent_basis)
assert all(len(entry) == 2 for entry in rank_two_system.laurent_basis)
assert all((entry[0] == 0) != (entry[1] == 0) for entry in rank_two_system.laurent_basis)
assert any(entry[0] == 1 and entry[1] == 0 for entry in rank_two_system.laurent_basis)
assert any(entry[0] == 0 and entry[1] == 1 for entry in rank_two_system.laurent_basis)

assert_value_error(
    lambda: quotient_coordinate_system(Matrix(R, [[0, 0]])),
    "finite-dimensional",
)

print("test_isomorphism_quotient_coordinates: ok")
