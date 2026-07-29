import _bootstrap  # noqa: F401

from sage.all import Matrix, identity_matrix

from isomorphism import R, x, y
from isomorphism.css import construct_excitation_map
from isomorphism.periods import (
    build_quotient_translation_action,
    periods_from_translation_action,
    superlattice_from_translation_action,
)

check_488 = Matrix(
    R,
    [
        [1 + y, 1 + y, 1 + x, 1 + x],
        [x * y, y, x * y, x],
    ],
)
action_488 = build_quotient_translation_action(check_488, diagnostics=True)
expected = Matrix([[1, 0], [1, 1]])
assert action_488.tx == expected
assert action_488.ty == expected
assert action_488.diagnostics["quotient_dimension"] == 2
assert action_488.diagnostics["monomial_basis_size"] == 2

periods_488 = periods_from_translation_action(action_488.tx, action_488.ty, max_period=4)
assert periods_488.square_period == 2
assert (2, 0) in periods_488.vectors
assert (1, 1) in periods_488.vectors

superlattice_488 = superlattice_from_translation_action(
    action_488.tx,
    action_488.ty,
    max_period=4,
)
assert superlattice_488.alpha == 2
assert superlattice_488.gamma == 1
assert superlattice_488.delta == 1
assert superlattice_488.cell == ((2, 0), (1, 1))

toric_full = build_quotient_translation_action(
    construct_excitation_map(1 + x, 1 + y),
    diagnostics=True,
)
assert toric_full.tx.nrows() == toric_full.ty.nrows()

try:
    periods_from_translation_action(Matrix([[0, 1], [0, 0]]), identity_matrix(2), max_period=3)
except ValueError as exc:
    assert "max_period" in str(exc)
else:
    raise AssertionError("Expected a bounded non-periodic action to fail.")

print("test_isomorphism_periods_translation: ok")
