import _bootstrap  # noqa: F401

from sage.all import Matrix

from isomorphism import (
    DecouplingMaps,
    DecouplingResult,
    InverseDecouplingMaps,
    ObliqueCell,
    PRESENTATION,
    PeriodData,
    R,
    P,
    SuperlatticeData,
    TranslationAction,
    x,
    y,
)
from isomorphism.css import (
    antipode_matrix,
    check_commutation,
    construct_excitation_map,
    evaluate_at_identity,
)

assert PRESENTATION.laurent_ring is R
assert PRESENTATION.polynomial_ring is P
assert PRESENTATION.inverse_relations() == (1 + PRESENTATION.a * PRESENTATION.c, 1 + PRESENTATION.b * PRESENTATION.d)
assert x.parent() is R
assert y.parent() is R

epsilon = construct_excitation_map(1 + x, 1 + y)
assert epsilon == Matrix(R, [[1 + x, 1 + y, 0, 0], [0, 0, 1 + y**-1, 1 + x**-1]])
assert check_commutation(epsilon, 2)
assert antipode_matrix(Matrix(R, [[x, y**-1]])) == Matrix(R, [[x**-1], [y]])
assert evaluate_at_identity(Matrix(R, [[1 + x, y]]))[0, 1] == 1
assert DecouplingMaps.__name__ == "DecouplingMaps"
assert InverseDecouplingMaps.__name__ == "InverseDecouplingMaps"
assert DecouplingResult.__name__ == "DecouplingResult"
assert PeriodData.__name__ == "PeriodData"
assert SuperlatticeData.__name__ == "SuperlatticeData"
assert TranslationAction.__name__ == "TranslationAction"
assert ObliqueCell.__name__ == "ObliqueCell"

print("test_isomorphism_imports: ok")
