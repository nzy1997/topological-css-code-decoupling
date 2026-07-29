import _bootstrap  # noqa: F401

from sage.all import GF, Matrix

from isomorphism import x, y
import isomorphism.chain_maps.elimination as elimination_module
from isomorphism.chain_maps.elimination import clear_css_mod_j
from isomorphism.css import (
    antipode_matrix,
    check_commutation,
    construct_excitation_map,
    evaluate_at_identity,
    split_css_blocks,
)


class NativeAddSpy:
    def __init__(self):
        self.calls = []

    def add_multiple_of_row(self, dst, src, scale):
        self.calls.append(("row", dst, src, scale))

    def add_multiple_of_column(self, dst, src, scale):
        self.calls.append(("column", dst, src, scale))

    def __getitem__(self, key):
        raise AssertionError("row/column helpers should use native in-place adds")


row_spy = NativeAddSpy()
elimination_module._row_add(row_spy, 1, 0)
assert row_spy.calls == [("row", 1, 0, 1)]

column_spy = NativeAddSpy()
elimination_module._col_add(column_spy, 1, 0, q=2)
assert column_spy.calls == [("column", 1, 0, 1), ("column", 2, 3, 1)]

epsilon = construct_excitation_map(1 + x + x * y, 1 + y + x * y)
result = clear_css_mod_j(epsilon, num_x_checks=1, num_qubits=2)

assert result.working_matrix == result.row_map * epsilon * result.symplectic_column_map
assert result.product_x_rank == 1
assert result.product_z_rank == 1
assert evaluate_at_identity(result.working_matrix) == Matrix(
    GF(2),
    [[1, 0, 0, 0], [0, 0, 0, 1]],
)
assert check_commutation(result.working_matrix, 2)

HZ_input, HX_input = split_css_blocks(epsilon, 1, 2)
HZ_working, HX_working = split_css_blocks(result.working_matrix, 1, 2)
assert HX_working == result.xi0.inverse() * HX_input * result.xi1
assert antipode_matrix(HZ_working) == (
    result.xi1.inverse() * antipode_matrix(HZ_input) * result.xi2
)

print("test_isomorphism_elimination: ok")
