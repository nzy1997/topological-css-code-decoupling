# Toric Form Workflow

This is the default ToricBuilder workflow for a two-polynomial presentation over
`GF(2)[x±1, y±1]`.

```julia
using ToricBuilder
using Oscar

Rxy, (x, y) = laurent_polynomial_ring(GF(2), ["x", "y"])
poly_vec = [1 + x + x*y, 1 + y + x*y]

l = find_L([x, y], ideal(Rxy, poly_vec))
A = excitation_matrix(poly_vec)
A_cg = coarse_grain(A, [x, y], [l, l])
result = build_toric_form(poly_vec; show_progress=false)
```

The main outputs are:

- `result.input_matrix`: the coarse-grained input matrix,
- `result.input_blocks.Hz` and `result.input_blocks.Hx`: the input CSS blocks,
- `result.standard_matrix`: the decoupled standard matrix,
- `result.standard_blocks.Hz` and `result.standard_blocks.Hx`: the standard CSS blocks,
- `result.row_transformation`: the row-side basis change,
- `result.phi_1`: the paper's degree-one decoupling map,
- `result.phi_1_inv`: optional, present only when `compute_inverse=true`,
- `result.column_transformation`: optional full symplectic column map, present only when `compute_inverse=true`,
- `result.product_state_num` and `result.toric_num`: decomposition counts.

For a fast end-to-end script, see
[`example/scripts/color_code.jl`](../../example/scripts/color_code.jl).

For lower-level reduction details, use the expert API documentation and source
instead of relying on private helpers.
