import _bootstrap  # noqa: F401

from sage.all import Matrix, identity_matrix

from isomorphism import R, x, y
import isomorphism.solver as solver_module
from isomorphism.solver import LaurentSolveError, find_inverse, solve_left, solve_right

A = Matrix(R, [[1, x], [0, 1]])
B = Matrix(R, [[1 + x], [y]])
X = solve_right(A, B)
assert A * X == B

left_A = Matrix(R, [[1, 0], [x, 1]])
left_B = Matrix(R, [[1 + y, x]])
Y = solve_left(left_A, left_B)
assert Y * left_A == left_B

A_inverse = find_inverse(A)
identity = identity_matrix(R, 2)
assert A * A_inverse == identity
assert A_inverse * A == identity

original_solve_right = solver_module.solve_right
original_native_inverse = solver_module._native_inverse


def reject_native_inverse(_A, *, verify=True):
    raise AssertionError("find_inverse should try solve_right before native inversion")


solver_module._native_inverse = reject_native_inverse
try:
    solved_inverse = solver_module.find_inverse(A)
finally:
    solver_module._native_inverse = original_native_inverse

assert A * solved_inverse == identity
assert solved_inverse * A == identity


def failing_solve_right(_A, _B, *, verify=True):
    raise LaurentSolveError("forced solve_right failure")


solver_module.solve_right = failing_solve_right
try:
    fallback_inverse = solver_module.find_inverse(A)
finally:
    solver_module.solve_right = original_solve_right

assert A * fallback_inverse == identity
assert fallback_inverse * A == identity

try:
    solve_right(Matrix(R, [[1 + x]]), Matrix(R, [[1]]))
except LaurentSolveError as exc:
    assert "generated module" in str(exc)
else:
    raise AssertionError("Expected an unsolvable Laurent system to fail.")

print("test_isomorphism_solver: ok")
