# A Degree-Growth Example for Direct Gaussian Elimination

## Summary

This note studies the color-code example

\[
f = 1+x+xy,
\qquad
g = 1+y+xy
\]

over the two-variable Laurent polynomial ring

\[
R = \mathbb F_2[x^{\pm 1},y^{\pm 1}].
\]

After a `3 x 3` coarse graining, the toric-reduction workflow must construct an
exact Laurent-polynomial map `psi_1_inverse`. The relevant public equation is

\[
H_{\mathrm{eff}}\,\psi_1^{-1} = \widetilde H_Z,
\qquad
H_{\mathrm{eff}} = P_Z H_Z,
\]

where `P_Z` is the returned `row_blocks.Hz`. In this example,
`H_eff` is `9 x 18`, `psi_1_inverse` is `18 x 18`, and `Hzt` is `9 x 18`.

The two direct Gaussian implementations tested here have different drawbacks:

1. Divide by polynomial pivots. This gives a valid solve only after extending
   `R` to its fraction field, and the result contains non-unit polynomial
   denominators.
2. Avoid division by cross-multiplying rows. This keeps all intermediate
   entries in `R`, but their degrees and term counts grow rapidly.

For the concrete system above, the division-free route grows from degree `2`
and `2` terms per entry to degree `35` and `96` terms per entry. The exact
Laurent-module witness returned by the current workflow has degree `1` and at
most `3` terms per entry.

These values are observations for the deterministic choices implemented in the
reproduction script: columns are visited from left to right, the first
available nonzero row is selected as the pivot, the division-free route does
not cancel row content, and the fraction-field route sets free variables to
zero. The experiment does not claim that every pivot order, fraction-free
algorithm, or choice of free variables must produce the same values.

## 1. Construction of the Equation

The experiment starts from the excitation matrix associated with `(f,g)` and
coarse-grains both translation directions by three:

```julia
Rxy, (x, y) = laurent_polynomial_ring(GF(2), ["x", "y"])
poly_vector = [1 + x + x*y, 1 + y + x*y]
input_matrix = coarse_grain(
    excitation_matrix(poly_vector),
    [x, y],
    [3, 3],
)
```

The resulting excitation matrix has size `18 x 36`. Its `Z` block has nine
rows and eighteen columns.

The current toric-reduction workflow first performs Gaussian elimination only
on the matrix obtained by evaluating `x=y=1`. Those field operations identify
the product-state and toric sectors without trying to invert arbitrary Laurent
polynomials. The workflow then fixes the standard toric target
`tilde(H)_Z` and constructs an exact Laurent witness by solving a matrix
equation.

For the public result, the equation can be extracted as follows:

```julia
debug = capture_toric_form_debug_matrices(input_matrix; show_progress=false)

input_Hz = debug.input_matrix[1:9, 1:18]
H_eff = debug.row_blocks.Hz * input_Hz
Hzt = debug.standard_blocks.Hz
psi_1_inverse = debug.psi_1_inverse

@assert H_eff * psi_1_inverse == Hzt
```

Thus the comparison below is performed on the exact equation certified by the
current public API, rather than on an artificial matrix.

## 2. Ordinary Gaussian Elimination

Suppose the current pivot is `p` and the entry to eliminate is `a`. Over a
field, the standard row update is

\[
R_i \leftarrow R_i - \frac{a}{p}R_r.
\]

This requires `p` to be invertible. In the Laurent ring `R`, the units are
nonzero monomials times nonzero field constants. A polynomial such as `1+x`
or `1+x+xy` is not a unit. Dividing by such a pivot therefore leaves `R`.

One can formally extend the computation to the fraction field `Frac(R)` and
run ordinary Gauss-Jordan elimination there. The implementation tested here
uses the first available nonzero pivot and sets all free variables to zero.
That produces a valid field solution, but it does not define a
Laurent-polynomial map for this choice. The resulting witness has the following
largest entries:

| quantity | degree | maximum terms |
| --- | ---: | ---: |
| numerator | 7 | 8 |
| denominator | 6 | 10 |

The nontrivial denominator is the obstruction: the returned map is valid over
`Frac(R)`, but it is not the local Laurent-polynomial transformation required
by the toric-reduction construction.

## 3. Division-Free Polynomial Elimination

To remain inside `R`, the naive implementation tested here avoids division by
cross-multiplying. It visits columns from left to right, chooses the first
available nonzero row, and performs no GCD or content cancellation. In
characteristic two, it uses

\[
R_i \leftarrow pR_i + aR_r.
\]

The entry in the pivot column becomes

\[
p a + a p = 0,
\]

so this update eliminates `a` without introducing a denominator.

### 3.1 Forward phase

For the augmented matrix `M = [H_eff | Hzt]`, the forward phase is:

```text
pivot_row = 1
for pivot_column in coefficient columns
    select the first nonzero entry at or below pivot_row
    swap that row into pivot_row
    p = M[pivot_row, pivot_column]

    for target_row below pivot_row
        a = M[target_row, pivot_column]
        M[target_row, :] = p*M[target_row, :] + a*M[pivot_row, :]
    end

    pivot_row += 1
end
```

### 3.2 Backward phase

The backward phase applies the same cross-multiplication from the last pivot
to the first:

```text
for pivot_row from last to first
    p = M[pivot_row, pivot_column]

    for target_row above pivot_row
        a = M[target_row, pivot_column]
        M[target_row, :] = p*M[target_row, :] + a*M[pivot_row, :]
    end
end
```

This resembles a division-free Gauss-Jordan reduction, but there is an
important algebraic caveat. Multiplication of a row by a non-unit is not an
invertible row operation over `R`. The procedure is therefore used here as a
diagnostic model of direct polynomial elimination, not as the exact solver.

### 3.3 Why the degree grows

If `deg` denotes the maximal Laurent monomial degree, one cross-elimination
step obeys the rough bound

\[
\deg(R_i^{\mathrm{new}})
\leq
\max\left(
\deg(p)+\deg(R_i),
\deg(a)+\deg(R_r)
\right).
\]

Every non-unit pivot is multiplied into an entire target row. Later pivots
already contain products created by earlier steps. The backward phase then
multiplies those accumulated pivots into rows that were processed previously.
Consequently, both polynomial degree and the number of monomials can grow much
faster during backward elimination than during the initial forward phase.

## 4. Metrics

For a Laurent polynomial

\[
h = \sum_{(a,b)} c_{a,b}x^a y^b,
\]

the example uses

\[
\deg(h) = \max_{c_{a,b}\neq 0} \left(|a|+|b|\right).
\]

For a matrix, the reported degree is the maximum degree of any entry. The
reported term count is the maximum number of nonzero monomials in any one
entry. These metrics track both spatial range and symbolic expression size.

For a fraction-field element, numerator and denominator statistics are
reported separately.

## 5. Observed Degree Explosion

The full trace for the division-free elimination is:

| phase | pivot row | maximum degree | maximum terms |
| --- | ---: | ---: | ---: |
| input | - | 2 | 2 |
| forward | 1 | 2 | 2 |
| forward | 2 | 2 | 2 |
| forward | 3 | 2 | 2 |
| forward | 4 | 2 | 2 |
| forward | 5 | 2 | 2 |
| forward | 6 | 2 | 3 |
| forward | 7 | 4 | 4 |
| forward | 8 | 8 | 10 |
| forward | 9 | 8 | 10 |
| backward | 9 | 12 | 18 |
| backward | 8 | 22 | 40 |
| backward | 7 | 35 | 96 |
| backward | 6 | 35 | 96 |
| backward | 5 | 35 | 96 |
| backward | 4 | 35 | 96 |
| backward | 3 | 35 | 96 |
| backward | 2 | 35 | 96 |
| backward | 1 | 35 | 96 |

The forward phase remains moderate through the first six pivots. The first
clear growth appears at pivots seven and eight, where the maximum degree rises
from `2` to `8`. The backward phase is the decisive failure point:

\[
(8,10)
\longrightarrow
(12,18)
\longrightarrow
(22,40)
\longrightarrow
(35,96),
\]

where each pair denotes `(maximum degree, maximum terms)`.

Relative to the input, the final largest degree is `17.5` times larger and the
largest term count is `48` times larger. This occurs in a `9 x 18` coefficient
matrix coming from the smallest `3 x 3` coarse-grained instance, not from a
large benchmark.

## 6. Comparison with the Exact Laurent Solve

The three measured routes can be summarized as follows:

| route | result | degree data | term data | acceptable Laurent map? |
| --- | --- | --- | --- | --- |
| ordinary Gauss-Jordan | fraction-field witness | numerator 7, denominator 6 | numerator 8, denominator 10 | no |
| division-free direct elimination | polynomial intermediates | input 2, final 35 | input 2, final 96 | not an invertible Laurent reduction |
| exact Laurent equation solve | `psi_1_inverse` | 1 | 3 | yes |

The current workflow verifies

\[
H_{\mathrm{eff}}\psi_1^{-1} = \widetilde H_Z
\]

exactly in `R`. It obtains a low-degree polynomial witness without introducing
fraction-field denominators and without accumulating the pivot products seen
in direct polynomial elimination.

Conceptually, `solve_laurent_linear` treats the problem as exact module
membership. The implementation prefers Oscar's native right-side solver and
uses a polynomial-module fallback when the native Laurent path is unavailable.
The fallback represents Laurent inverses through quotient relations, computes
module coordinates, maps the witness back to `R`, and verifies the original
matrix equation. See
[`laurent_linear_solve.jl`](../../src/core/laurent_linear_solve.jl) for the
implementation.

## 7. Reproduction

Run the checked example from the repository root:

```bash
julia --project=julia/ToricBuilder \
  julia/ToricBuilder/example/scripts/laurent_gaussian_degree_growth.jl
```

The implementation is in
[`laurent_gaussian_degree_growth.jl`](../../example/scripts/laurent_gaussian_degree_growth.jl).
It verifies the exact equation, records every forward and backward pivot, runs
the fraction-field comparison, and checks the expected metrics.

The degree-growth result in this note comes from applying the two direct
Gaussian alternatives to the exact matrix equation used by the current
workflow. It is a diagnostic example rather than the implementation used by
the main toric-form workflow.

## Conclusion

For the deterministic implementations tested here, direct Gaussian elimination
does not preserve all of the properties required by the toric-reduction map:

- dividing by non-unit pivots produces a fraction-field map rather than a
  Laurent-polynomial map;
- avoiding those divisions by cross-multiplication causes severe symbolic
  growth and does not give invertible Laurent row operations.

This experiment does not rule out a better pivot order, a Bareiss-style
fraction-free method, content cancellation, or a different assignment of free
variables. It does show that a straightforward direct implementation can fail
in either of two concrete ways on the `1+x+xy` example at the smallest
coarse-grained size: the tested field solve introduces nonlocal denominators,
while the tested division-free schedule grows to degree `35` and `96` terms.
The exact Laurent equation solve avoids both outcomes in this instance.
