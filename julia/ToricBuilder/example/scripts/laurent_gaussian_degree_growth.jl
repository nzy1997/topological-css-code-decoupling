using Printf
using ToricBuilder
using Oscar

# See docs/expert/laurent_gaussian_degree_growth.md for the derivation and
# interpretation of the intermediate-expression growth measured here.

function laurent_matrix_stats(A)
    max_degree = 0
    max_terms = 0

    for entry in A
        iszero(entry) && continue
        exponent_vectors = collect(AbstractAlgebra.exponent_vectors(entry))
        entry_degree = maximum(sum(abs, exponents) for exponents in exponent_vectors)
        max_degree = max(max_degree, entry_degree)
        max_terms = max(max_terms, length(exponent_vectors))
    end

    return (; max_degree, max_terms)
end

function fraction_matrix_stats(A)
    numerator_degree = 0
    denominator_degree = 0
    numerator_terms = 0
    denominator_terms = 0

    for entry in A
        numerator_stats = laurent_matrix_stats([numerator(entry)])
        denominator_stats = laurent_matrix_stats([denominator(entry)])
        numerator_degree = max(numerator_degree, numerator_stats.max_degree)
        denominator_degree = max(denominator_degree, denominator_stats.max_degree)
        numerator_terms = max(numerator_terms, numerator_stats.max_terms)
        denominator_terms = max(denominator_terms, denominator_stats.max_terms)
    end

    return (; numerator_degree, denominator_degree, numerator_terms, denominator_terms)
end

function swap_rows!(A, first_row::Int, second_row::Int)
    first_row == second_row && return A
    for column in 1:ncols(A)
        A[first_row, column], A[second_row, column] =
            A[second_row, column], A[first_row, column]
    end
    return A
end

function cross_eliminate_row!(A, target_row::Int, pivot_row::Int, pivot_column::Int; first_column::Int)
    pivot = A[pivot_row, pivot_column]
    coefficient = A[target_row, pivot_column]
    iszero(coefficient) && return A

    # Over GF(2), pivot * coefficient + coefficient * pivot is zero. This
    # avoids division by a non-unit pivot, but multiplies the whole row by it.
    for column in first_column:ncols(A)
        A[target_row, column] =
            pivot * A[target_row, column] + coefficient * A[pivot_row, column]
    end
    return A
end

"""
    direct_polynomial_gaussian_trace(A, B)

Apply naive division-free Gaussian elimination to `[A B]`. Polynomial pivots
are cross-multiplied instead of inverted, so every entry remains in the Laurent
ring. The returned trace records the degree and term-count growth after every
forward and backward pivot.

The schedule visits columns from left to right, picks the first available
nonzero row, and performs no GCD or content cancellation.

This is a diagnostic counterexample, not a replacement linear solver:
cross-multiplication by non-units is not an invertible Laurent row operation.
"""
function direct_polynomial_gaussian_trace(A, B)
    nrows(A) == nrows(B) || throw(DimensionMismatch("A and B must have the same number of rows"))
    base_ring(A) == base_ring(B) || throw(ArgumentError("A and B must use the same base ring"))

    augmented = hcat(copy(A), copy(B))
    coefficient_columns = ncols(A)
    row_count = nrows(A)
    pivot_row = 1
    pivot_columns = Int[]
    forward_trace = NamedTuple[]
    backward_trace = NamedTuple[]

    for pivot_column in 1:coefficient_columns
        pivot_offset = findfirst(!iszero, augmented[pivot_row:row_count, pivot_column])
        isnothing(pivot_offset) && continue

        selected_row = pivot_row + pivot_offset - 1
        swap_rows!(augmented, pivot_row, selected_row)
        for target_row in pivot_row+1:row_count
            cross_eliminate_row!(
                augmented,
                target_row,
                pivot_row,
                pivot_column;
                first_column=pivot_column,
            )
        end

        push!(pivot_columns, pivot_column)
        push!(
            forward_trace,
            (; phase=:forward, pivot=pivot_row, column=pivot_column, laurent_matrix_stats(augmented)...),
        )
        pivot_row += 1
        pivot_row > row_count && break
    end

    for row in length(pivot_columns):-1:1
        pivot_column = pivot_columns[row]
        for target_row in 1:row-1
            cross_eliminate_row!(
                augmented,
                target_row,
                row,
                pivot_column;
                first_column=1,
            )
        end
        push!(
            backward_trace,
            (; phase=:backward, pivot=row, column=pivot_column, laurent_matrix_stats(augmented)...),
        )
    end

    return (;
        pivot_columns,
        input_stats=laurent_matrix_stats(hcat(A, B)),
        forward_trace,
        backward_trace,
        final_stats=laurent_matrix_stats(augmented),
    )
end

"""
    direct_fraction_field_gaussian(A, B)

Apply ordinary Gauss-Jordan elimination after extending the Laurent ring to its
fraction field. This is a valid field solve, but the result generally contains
non-unit polynomial denominators and therefore is not a Laurent-polynomial map.
The returned witness uses the first available nonzero pivot and sets all free
variables to zero.
"""
function direct_fraction_field_gaussian(A, B)
    K = fraction_field(base_ring(A))
    augmented = hcat(map_entries(K, A), map_entries(K, B))
    row_count = nrows(A)
    coefficient_columns = ncols(A)
    pivot_row = 1
    pivot_columns = Int[]

    for pivot_column in 1:coefficient_columns
        pivot_offset = findfirst(!iszero, augmented[pivot_row:row_count, pivot_column])
        isnothing(pivot_offset) && continue

        swap_rows!(augmented, pivot_row, pivot_row + pivot_offset - 1)
        pivot = augmented[pivot_row, pivot_column]
        for column in pivot_column:ncols(augmented)
            augmented[pivot_row, column] = divexact(augmented[pivot_row, column], pivot)
        end

        for target_row in 1:row_count
            target_row == pivot_row && continue
            coefficient = augmented[target_row, pivot_column]
            iszero(coefficient) && continue
            for column in pivot_column:ncols(augmented)
                augmented[target_row, column] -= coefficient * augmented[pivot_row, column]
            end
        end

        push!(pivot_columns, pivot_column)
        pivot_row += 1
        pivot_row > row_count && break
    end

    solution = zero_matrix(K, coefficient_columns, ncols(B))
    for (row, pivot_column) in enumerate(pivot_columns)
        solution[pivot_column, :] = augmented[row, coefficient_columns+1:end]
    end

    map_entries(K, A) * solution == map_entries(K, B) ||
        error("fraction-field Gaussian solution failed verification")
    return (; solution, stats=fraction_matrix_stats(solution))
end

function print_trace(trace, solution_stats)
    println("Metric: max |a|+|b| over monomials x^a*y^b; max terms in one entry")
    @printf("%-12s %8s %8s %8s\n", "stage", "pivot", "degree", "terms")
    @printf("%-12s %8s %8d %8d\n", "input", "-", trace.input_stats.max_degree, trace.input_stats.max_terms)
    for step in trace.forward_trace
        @printf("%-12s %8d %8d %8d\n", "forward", step.pivot, step.max_degree, step.max_terms)
    end
    for step in trace.backward_trace
        @printf("%-12s %8d %8d %8d\n", "backward", step.pivot, step.max_degree, step.max_terms)
    end
    @printf("%-12s %8s %8d %8d\n", "module Q", "-", solution_stats.max_degree, solution_stats.max_terms)
end

function main(; verify_expected::Bool=true)
    Rxy, (x, y) = laurent_polynomial_ring(GF(2), ["x", "y"])
    poly_vector = [1 + x + x*y, 1 + y + x*y]
    input_matrix = coarse_grain(excitation_matrix(poly_vector), [x, y], [3, 3])

    # The current path first fixes the standard target and then solves an
    # equivalent local equation as an exact Laurent-module problem. Its public
    # witness satisfies (row_Hz * input_Hz) * phi_1 = Hzt.
    debug = capture_toric_form_debug_matrices(input_matrix; show_progress=false)
    stabilizer_count = div(nrows(input_matrix), 2)
    qubit_count = div(ncols(input_matrix), 2)
    input_Hz = debug.input_matrix[1:stabilizer_count, 1:qubit_count]
    equation_Hz = debug.row_blocks.Hz * input_Hz
    Hzt = debug.standard_blocks.Hz
    phi_1 = debug.phi_1

    equation_Hz * phi_1 == Hzt ||
        error("current Laurent-module solution failed equation_Hz * phi_1 == Hzt")

    trace = direct_polynomial_gaussian_trace(equation_Hz, Hzt)
    fraction_result = direct_fraction_field_gaussian(equation_Hz, Hzt)
    solution_stats = laurent_matrix_stats(phi_1)
    print_trace(trace, solution_stats)
    println()
    println("Fraction-field Gaussian witness: ", fraction_result.stats)

    if verify_expected
        trace.input_stats == (max_degree=2, max_terms=2) ||
            error("unexpected input metrics: $(trace.input_stats)")
        trace.final_stats == (max_degree=35, max_terms=96) ||
            error("unexpected direct-elimination metrics: $(trace.final_stats)")
        solution_stats == (max_degree=1, max_terms=3) ||
            error("unexpected module-solution metrics: $solution_stats")
        fraction_result.stats == (
            numerator_degree=7,
            denominator_degree=6,
            numerator_terms=8,
            denominator_terms=10,
        ) || error("unexpected fraction-field metrics: $(fraction_result.stats)")
    end

    println()
    println("Verified: equation_Hz * phi_1 == Hzt")
    println("Direct polynomial elimination: degree 2 -> $(trace.final_stats.max_degree), terms 2 -> $(trace.final_stats.max_terms)")
    println("Ordinary field Gaussian: denominator degree $(fraction_result.stats.denominator_degree), denominator terms $(fraction_result.stats.denominator_terms)")
    println("Exact Laurent-module witness: degree $(solution_stats.max_degree), terms $(solution_stats.max_terms)")

    return (; poly_vector, input_matrix, equation_Hz, Hzt, phi_1, trace, fraction_result, solution_stats)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
