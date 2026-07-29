import _bootstrap  # noqa: F401

from sage.all import FreeModule, Matrix

from isomorphism import P, R, a, b, c, d, x, y
from isomorphism.algebra.singular_backend import (
    laurent_to_polynomial,
    matrix_laurent_to_polynomial,
    matrix_polynomial_to_laurent,
    membership_lift,
    polynomial_to_laurent,
    quotient_dimension,
    quotient_monomial_basis,
    reduce_by_standard_basis,
    standard_basis,
)
from isomorphism.rings import module_element

assert laurent_to_polynomial(1 + x + x**-1 * y**2) == 1 + a + c * b**2
assert polynomial_to_laurent(a * c + b * d + c * b**2) == x**-1 * y**2
assert matrix_laurent_to_polynomial(Matrix(R, [[1 + x, x**-1 * y]])) == Matrix(
    P, [[1 + a, c * b]]
)
assert matrix_polynomial_to_laurent(Matrix(P, [[a * c + b * d + c * b]])) == Matrix(
    R, [[x**-1 * y]]
)

module = FreeModule(P, 1)
generators = [
    module_element(module, [1 + a + a * b]),
    module_element(module, [1 + b + a * b]),
    module_element(module, [1 + a * c]),
    module_element(module, [1 + b * d]),
]
groebner_basis = standard_basis(generators)
assert quotient_dimension(groebner_basis) == 2

monomial_basis = quotient_monomial_basis(groebner_basis)
assert len(monomial_basis) == 2
reduction_input = a**2 * module.basis()[0]
reduced = reduce_by_standard_basis(reduction_input, groebner_basis)
assert reduced == d * module.basis()[0]
assert reduction_input.is_mutable()
assert reduce_by_standard_basis(generators[0], groebner_basis).is_zero()

lift_generators = Matrix(P, [[a, 1 + a]])
lift_targets = Matrix(P, [[a + a**2]])
lift_witness = membership_lift(lift_generators, lift_targets)
assert lift_generators * lift_witness == lift_targets

try:
    membership_lift(Matrix(P, [[a]]), Matrix(P, [[a], [1]]))
except ValueError as exc:
    assert str(exc) == "Singular lift matrices must have the same row count."
else:
    raise AssertionError("membership_lift accepted matrices with different row counts.")

try:
    membership_lift(Matrix(R, [[x]]), Matrix(R, [[x]]))
except ValueError as exc:
    assert str(exc) == "Singular lift generators must be a matrix over P."
else:
    raise AssertionError("membership_lift accepted Laurent-ring generators.")

print("test_isomorphism_singular_backend: ok")
