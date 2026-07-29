using ToricBuilder
using Oscar

const CODE_LABEL = "BB/Gross [[144,12,12]]"

bit_entry(a) = iszero(a) ? 0 : 1

function print_binary_matrix(name, A)
    println(name, " =")
    for row in 1:nrows(A)
        println("  ", join((bit_entry(A[row, col]) for col in 1:ncols(A)), " "))
    end
end

function print_int_matrix(name, A)
    println(name, " =")
    for row in 1:size(A, 1)
        println("  ", join((A[row, col] for col in 1:size(A, 2)), " "))
    end
end

function matrix_to_bits(A)
    return [bit_entry(A[row, col]) for row in 1:nrows(A), col in 1:ncols(A)]
end

f2_identity(n) = [i == j ? 1 : 0 for i in 1:n, j in 1:n]

function f2_add(A, B)
    rows, cols = size(A)
    return [(A[row, col] + B[row, col]) % 2 for row in 1:rows, col in 1:cols]
end

function f2_mul(A, B)
    rows = size(A, 1)
    inner = size(A, 2)
    cols = size(B, 2)
    return [sum(A[row, k] * B[k, col] for k in 1:inner) % 2 for row in 1:rows, col in 1:cols]
end

function f2_pow(A, exponent::Int)
    result = f2_identity(size(A, 1))
    base = A
    power = exponent
    while power > 0
        if isodd(power)
            result = f2_mul(result, base)
        end
        base = f2_mul(base, base)
        power ÷= 2
    end
    return result
end

function f2_rank(M)
    A = copy(M)
    rows, cols = size(A)
    rank = 0
    pivot_row = 1

    for col in 1:cols
        pivot = findfirst(row -> A[row, col] == 1, pivot_row:rows)
        if isnothing(pivot)
            continue
        end

        pivot += pivot_row - 1
        A[pivot_row, :], A[pivot, :] = copy(A[pivot, :]), copy(A[pivot_row, :])
        for row in 1:rows
            if row != pivot_row && A[row, col] == 1
                A[row, :] = (A[row, :] + A[pivot_row, :]) .% 2
            end
        end

        rank += 1
        pivot_row += 1
        pivot_row > rows && break
    end

    return rank
end

function f2_nullspace(M)
    A = copy(M)
    rows, cols = size(A)
    pivot_row = 1
    pivots = Int[]

    for col in 1:cols
        pivot = findfirst(row -> A[row, col] == 1, pivot_row:rows)
        if isnothing(pivot)
            continue
        end

        pivot += pivot_row - 1
        A[pivot_row, :], A[pivot, :] = copy(A[pivot, :]), copy(A[pivot_row, :])
        for row in 1:rows
            if row != pivot_row && A[row, col] == 1
                A[row, :] = (A[row, :] + A[pivot_row, :]) .% 2
            end
        end

        push!(pivots, col)
        pivot_row += 1
        pivot_row > rows && break
    end

    free_cols = setdiff(collect(1:cols), pivots)
    basis = Vector{Vector{Int}}()
    for free_col in free_cols
        v = zeros(Int, cols)
        v[free_col] = 1
        for (row, pivot_col) in enumerate(pivots)
            if A[row, free_col] == 1
                v[pivot_col] = 1
            end
        end
        push!(basis, v)
    end
    return basis
end

function f2_inverse(A)
    n = size(A, 1)
    augmented = hcat(copy(A), f2_identity(n))
    pivot_row = 1

    for col in 1:n
        pivot = findfirst(row -> augmented[row, col] == 1, pivot_row:n)
        isnothing(pivot) && error("Matrix is not invertible over GF(2)")

        pivot += pivot_row - 1
        augmented[pivot_row, :], augmented[pivot, :] = copy(augmented[pivot, :]), copy(augmented[pivot_row, :])
        for row in 1:n
            if row != pivot_row && augmented[row, col] == 1
                augmented[row, :] = (augmented[row, :] + augmented[pivot_row, :]) .% 2
            end
        end
        pivot_row += 1
    end

    return augmented[:, n+1:2n]
end

function column_rank(columns::Vector{Vector{Int}})
    isempty(columns) && return 0
    return f2_rank(hcat(columns...))
end

function extend_independent_basis!(basis::Vector{Vector{Int}}, candidates; target_dim::Int)
    for candidate in candidates
        if column_rank([basis; [candidate]]) > length(basis)
            push!(basis, candidate)
            length(basis) == target_dim && return basis
        end
    end
    error("Unable to extend independent basis to dimension $target_dim")
end

function unit_vector(n, idx)
    v = zeros(Int, n)
    v[idx] = 1
    return v
end

function standard_unit_basis(n)
    return [unit_vector(n, idx) for idx in 1:n]
end

function change_basis_matrix(A, basis_columns)
    P = hcat(basis_columns...)
    return f2_mul(f2_inverse(P), f2_mul(A, P))
end

function coordinates_in_standard_basis(v, basis_columns)
    return f2_mul(hcat(basis_columns...), reshape(v, :, 1))[:, 1]
end

function vector_label(v, basis_labels)
    terms = [basis_labels[i] for i in eachindex(v) if v[i] == 1]
    return isempty(terms) ? "0" : join(terms, " + ")
end

function build_y_period_basis(Ty)
    dim = size(Ty, 1)
    identity = f2_identity(dim)
    fixed_by_y3 = f2_nullspace(f2_add(f2_pow(Ty, 3), identity))
    fixed_by_y6 = f2_nullspace(f2_add(f2_pow(Ty, 6), identity))

    basis = copy(fixed_by_y3)
    extend_independent_basis!(basis, fixed_by_y6; target_dim=6)
    extend_independent_basis!(basis, standard_unit_basis(dim); target_dim=dim)
    return basis
end

function build_y6_swap_basis(Ty6_adapted)
    basis = [
        unit_vector(8, 7),
        Ty6_adapted[:, 7],
        unit_vector(8, 8),
        Ty6_adapted[:, 8],
    ]
    labels = [
        "s1 = c7",
        "s2 = T_y^6(c7)",
        "s3 = c8",
        "s4 = T_y^6(c8)",
    ]

    for idx in 1:6
        candidate = unit_vector(8, idx)
        if column_rank([basis; [candidate]]) > length(basis)
            push!(basis, candidate)
            push!(labels, "s$(length(basis)) = c$idx")
            length(basis) == 8 && break
        end
    end
    return basis, labels
end

function print_anyon_representation(rep)
    println("anyon_count = ", rep.anyon_count)
    println("anyon_basis =")
    for (idx, basis_element) in enumerate(rep.basis)
        println("  b", idx, " = ", basis_element)
    end
    println()

    println("Matrix convention: column j is the coordinate vector of the translated basis element.")
    print_binary_matrix("T_x", rep.Ax)
    println()
    print_binary_matrix("T_y", rep.Ay)
    println()
    println("commute_check = ", rep.Ax * rep.Ay == rep.Ay * rep.Ax)
end

function print_y_period_basis(period_basis, monomial_labels)
    adapted_labels = [
        "c1(period|3)",
        "c2(period|3)",
        "c3(period|3)",
        "c4(period|3)",
        "c5(period|6)",
        "c6(period|6)",
        "c7(period|12)",
        "c8(period|12)",
    ]

    println()
    println("y-period adapted anyon basis:")
    for (idx, v) in enumerate(period_basis)
        println("  ", adapted_labels[idx], " = ", vector_label(v, monomial_labels))
    end
    return adapted_labels
end

function print_nontrivial_y6_translations(Ty6_adapted, adapted_labels)
    println()
    println("Nontrivial 6-step y translations:")
    for col in 1:size(Ty6_adapted, 2)
        delta = copy(Ty6_adapted[:, col])
        delta[col] = (delta[col] + 1) % 2
        if any(!iszero, delta)
            image = vector_label(Ty6_adapted[:, col], adapted_labels)
            shift = vector_label(delta, adapted_labels)
            println("  T_y^6(", adapted_labels[col], ") = ", image, "  differs by  ", shift)
        end
    end
end

function print_y6_swap_pairs(Ty6_adapted, adapted_labels)
    println()
    println("6-step swap pairs for the two 12-period quotient directions:")
    for col in 7:8
        u = unit_vector(length(adapted_labels), col)
        v = Ty6_adapted[:, col]
        delta = copy(v)
        delta[col] = (delta[col] + 1) % 2
        println("  pair from ", adapted_labels[col], ":")
        println("    u = ", vector_label(u, adapted_labels))
        println("    v = T_y^6(u) = ", vector_label(v, adapted_labels))
        println("    T_y^6(v) = ", vector_label(f2_mul(Ty6_adapted, reshape(v, :, 1))[:, 1], adapted_labels))
        println("    u + v = ", vector_label(delta, adapted_labels), "  (fixed by T_y^6)")
    end
end

function print_swap_basis(Ty6_adapted, period_basis, adapted_labels, monomial_labels)
    swap_basis, swap_labels = build_y6_swap_basis(Ty6_adapted)
    Ty6_swap = change_basis_matrix(Ty6_adapted, swap_basis)

    println()
    println("Swap-pair basis for T_y^6:")
    for (idx, v) in enumerate(swap_basis)
        standard_coords = coordinates_in_standard_basis(v, period_basis)
        println("  ", swap_labels[idx])
        println("    in c-basis: ", vector_label(v, adapted_labels))
        println("    in monomial basis: ", vector_label(standard_coords, monomial_labels))
    end
    println()
    print_int_matrix("T_y^6 in swap-pair basis", Ty6_swap)
end

function main()
    F = GF(2)
    Rxy, (x, y) = laurent_polynomial_ring(F, ["x", "y"])

    # BB/Gross [[144,12,12]] polynomial presentation.
    f = 1 + x + x^-1 * y^3
    g = 1 + y + y^-1 * x^3
    rep = construct_anyon_translation_representation([f, g])

    println("code: ", CODE_LABEL)
    println("f = ", f)
    println("g = ", g)
    println()

    print_anyon_representation(rep)

    Tx = matrix_to_bits(rep.Ax)
    Ty = matrix_to_bits(rep.Ay)
    period_basis = build_y_period_basis(Ty)
    monomial_labels = string.(rep.basis)
    adapted_labels = print_y_period_basis(period_basis, monomial_labels)

    Tx_adapted = change_basis_matrix(Tx, period_basis)
    Ty_adapted = change_basis_matrix(Ty, period_basis)
    Txy_adapted = change_basis_matrix(f2_mul(Tx, Ty), period_basis)
    Ty6_adapted = f2_pow(Ty_adapted, 6)

    println()
    println("Adapted-basis convention: column j is the translated c_j.")
    print_int_matrix("T_x in y-period adapted basis", Tx_adapted)
    println()
    print_int_matrix("T_y in y-period adapted basis", Ty_adapted)
    println()
    print_int_matrix("T_x*T_y in y-period adapted basis", Txy_adapted)
    println()
    print_int_matrix("T_y^6 in y-period adapted basis", Ty6_adapted)

    print_nontrivial_y6_translations(Ty6_adapted, adapted_labels)
    print_y6_swap_pairs(Ty6_adapted, adapted_labels)
    print_swap_basis(Ty6_adapted, period_basis, adapted_labels, monomial_labels)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
