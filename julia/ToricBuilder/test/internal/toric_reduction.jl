using Test
using ToricBuilder
using ToricBuilder: ExcitationMatrix, ExcitationMatrixwithF2Matrix
using Oscar

@testset "gaussian_elimination! API" begin
    @test isdefined(ToricBuilder, :gaussian_elimination!)
end

if isdefined(ToricBuilder, :gaussian_elimination!)
    const gaussian_elimination! = getfield(ToricBuilder, :gaussian_elimination!)

    @testset "gaussian_elimination!" begin
        F = GF(2)
        Rxy, (x, y) = laurent_polynomial_ring(F, ["x", "y"])
        A = excitation_matrix([1 + x + x * y, 1 + y + x * y])
        A4 = coarse_graining(A, [x, y], [3, 3])

        E = ExcitationMatrix(A4)
        E2 = ExcitationMatrixwithF2Matrix(E)
        E2_copy = copy(E2)
        gaussian_elimination!(E2, 1, 1, 18, 1, 1, 9)
        gaussian_elimination!(E2, 26, 19, 36, 10, 10, 18)

        @test E2.m == E2.em.row_transformation * E2_copy.m * E2.em.column_transformation
        @test E2.em.m == E2.em.row_transformation * E2_copy.em.m * E2.em.column_transformation
    end
end

@testset "psi_1_inverse toric correction helpers match dense correction" begin
    F = GF(2)
    Rxy, (x, y) = laurent_polynomial_ring(F, ["x", "y"])

    product_state_num = 2
    toric_num = 2
    stab_num = product_state_num + toric_num
    qubit_num = 2 * product_state_num + 2 * toric_num

    Hxdagger = matrix(Rxy, qubit_num, stab_num, [
        1       x       0       y;
        y       1       x + 1   0;
        x^-1    y^-1    1       x*y;
        0       x       y       1;
        1 + y   0       x^-1    y^2;
        x       1       y^-1    0;
        y       x*y     1       x + y;
        0       y^-1    x^2     1;
    ])

    Adagger = zero_matrix(Rxy, qubit_num, stab_num)
    Adagger[2*product_state_num+1, product_state_num+1] = x + 1
    Adagger[2*product_state_num+2, product_state_num+1] = y^-1 + x
    Adagger[2*product_state_num+3, product_state_num+2] = x*y + 1
    Adagger[2*product_state_num+4, product_state_num+2] = x^-1 + y

    psi_1_inverse_base = matrix(Rxy, qubit_num, qubit_num, [
        1 x 0 y 1 0 x*y 1;
        0 1 y 0 x 1 0 y^-1;
        x^-1 0 1 x*y 0 y 1 0;
        y 1 0 1 x^-1 0 y 1;
        1 0 y^-1 x 1 x*y 0 1;
        x y 1 0 0 1 y^-1 x^-1;
        0 1 x 1 y 0 1 x*y;
        y^-1 0 1 x 1 y 0 1;
    ])

    C = matrix(Rxy, qubit_num, qubit_num, [
        1 0 x 0 1 0 y 0;
        0 1 0 y 0 1 0 x;
        x 0 1 0 y 0 1 0;
        0 y 0 1 0 x 0 1;
        1 0 y 0 1 0 x 0;
        0 1 0 x 0 1 0 y;
        y 0 1 0 x 0 1 0;
        0 x 0 1 0 y 0 1;
    ])

    correction = ToricBuilder._build_psi1_inverse_toric_correction(
        Hxdagger,
        Adagger,
        product_state_num,
        toric_num,
    )

    dense_updated = psi_1_inverse_base + Hxdagger * ToricBuilder._dagger_laurent_matrix(Adagger)
    block_updated = ToricBuilder._apply_psi1_inverse_toric_correction!(copy(psi_1_inverse_base), correction)
    fused = ToricBuilder._compose_psi1_inverse_with_toric_correction(C, psi_1_inverse_base, correction)

    @test block_updated == dense_updated
    @test fused == C * dense_updated
end

@testset "selected column left multiply matches dense composition" begin
    F = GF(2)
    Rxy, (x, y) = laurent_polynomial_ring(F, ["x", "y"])

    C = matrix(Rxy, 4, 4, [
        1       x       0       y;
        0       1       x + 1   0;
        x^-1    0       1       y^-1;
        1 + y   x*y     0       1;
    ])
    M = zero_matrix(Rxy, 4, 6)
    M[:, 1] = matrix(Rxy, 4, 1, [1; x; 0; y])
    M[:, 3] = matrix(Rxy, 4, 1, [x^-1; 0; 1 + y; 0])
    M[:, 5] = matrix(Rxy, 4, 1, [0; y^-1; x*y; 1])
    M[:, 6] = matrix(Rxy, 4, 1, [x + 1; 0; y; x^-1])

    selected = ToricBuilder._left_multiply_column_blocks(
        C,
        M,
        (1:1, 3:3, 5:6),
    )

    @test selected == C * M
    @test iszero(ToricBuilder._left_multiply_column_blocks(C, M, UnitRange{Int}[]))
end

@testset "standard Hxt shape validation" begin
    F = GF(2)
    Rxy, (x, y) = laurent_polynomial_ring(F, ["x", "y"])

    product_state_num = 2
    toric_num = 2
    stab_num = product_state_num + toric_num
    qubit_num = 2 * product_state_num + 2 * toric_num

    Hxt = zero_matrix(Rxy, stab_num, qubit_num)
    for i in 1:product_state_num
        Hxt[i, i + product_state_num] = Rxy(1)
    end
    for i in 1:toric_num
        Hxt[product_state_num + i, 2 * product_state_num + 2 * i - 1] = y^-1 + 1
        Hxt[product_state_num + i, 2 * product_state_num + 2 * i] = x^-1 + 1
    end

    @test isnothing(ToricBuilder._assert_standard_hxt_shape(Hxt, product_state_num, toric_num))

    bad_product = copy(Hxt)
    bad_product[1, product_state_num + 1] = zero(Rxy)
    @test_throws ErrorException ToricBuilder._assert_standard_hxt_shape(
        bad_product,
        product_state_num,
        toric_num,
    )

    bad_extra = copy(Hxt)
    bad_extra[1, 1] = Rxy(1)
    @test_throws ErrorException ToricBuilder._assert_standard_hxt_shape(
        bad_extra,
        product_state_num,
        toric_num,
    )
end

@testset "to_toric_form color code" begin
    F = GF(2)
    Rxy, (x, y) = laurent_polynomial_ring(F, ["x", "y"])
    A = excitation_matrix([1 + x + x*y, 1 + y + x*y])
    A4 = coarse_graining(A, [x, y], [3, 3])

    res = to_toric_form(A4)
    @test check_result(res, A4)
end

@testset "to_toric_form default returns block data without inverse" begin
    F = GF(2)
    Rxy, (x, y) = laurent_polynomial_ring(F, ["x", "y"])
    A = excitation_matrix([1 + x + x*y, 1 + y + x*y])
    A4 = coarse_graining(A, [x, y], [3, 3])

    res = to_toric_form(A4; show_progress=false)

    @test hasproperty(res, :input_matrix)
    @test hasproperty(res, :standard_matrix)
    @test hasproperty(res, :input_blocks)
    @test hasproperty(res, :standard_blocks)
    @test hasproperty(res, :row_blocks)
    @test hasproperty(res, :psi_1_inverse)
    @test res.psi_1 === nothing
    @test !hasproperty(res, Symbol("ph", "i_1"))
    @test !hasproperty(res, Symbol("ph", "i_1_inv"))
    @test res.input_matrix == A4
    @test res.column_transformation === nothing
    @test check_result(res, A4)

    stab_num = size(A4, 1) ÷ 2
    qubit_num = size(A4, 2) ÷ 2
    @test res.input_blocks.Hz == A4[1:stab_num, 1:qubit_num]
    @test res.input_blocks.Hx == A4[stab_num+1:2*stab_num, qubit_num+1:2*qubit_num]
    @test res.standard_matrix == [
        res.standard_blocks.Hz zero_matrix(Rxy, stab_num, qubit_num)
        zero_matrix(Rxy, stab_num, qubit_num) res.standard_blocks.Hx
    ]
end

@testset "to_toric_form bb code L = 3" begin
    F = GF(2)
    Rxy, (x, y) = laurent_polynomial_ring(F, ["x", "y"])
    A = excitation_matrix([1 + x + x^-1, 1 + y + y^-1])
    A4 = coarse_graining(A, [x, y], [3, 3])

    res = to_toric_form(A4)
    @test check_result(res, A4)
end

@testset "to_toric_form enters polynomial correction branch on a real coarse-grained case" begin
    F = GF(2)
    Rxy, (x, y) = laurent_polynomial_ring(F, ["x", "y"])
    A = excitation_matrix([1 + x + x^-1 * y^3, 1 + y + y^-1 * x^3])
    A4 = coarse_graining(A, [x, y], [12, 12])
    res_ref = Ref{Any}(nothing)

    output = mktemp() do path, io
        redirect_stdout(io) do
            res_ref[] = to_toric_form(A4; show_progress=true)
        end
        flush(io)
        close(io)
        @test check_result(res_ref[], A4)
        read(path, String)
    end

    @test occursin("det(B) is a polynomial", output)
end

@testset "to_toric_form inverse mode returns full column transformation" begin
    F = GF(2)
    Rxy, (x, y) = laurent_polynomial_ring(F, ["x", "y"])
    A = excitation_matrix([1 + x + x*y, 1 + y + x*y])
    A4 = coarse_graining(A, [x, y], [3, 3])

    res = to_toric_form(A4; show_progress=false, compute_inverse=true)

    @test !isnothing(res.psi_1)
    @test !isnothing(res.column_transformation)
    @test check_result(res, A4; require_inverse=true)
    identity = identity_matrix(Rxy, size(res.psi_1_inverse, 1))
    @test res.psi_1_inverse * res.psi_1 == identity
    @test res.psi_1 * res.psi_1_inverse == identity
end

@testset "capture_toric_form_debug_matrices" begin
    F = GF(2)
    Rxy, (x, y) = laurent_polynomial_ring(F, ["x", "y"])
    A = excitation_matrix([1 + x + x*y, 1 + y + x*y])
    A4 = coarse_graining(A, [x, y], [3, 3])

    debug = ToricBuilder.capture_toric_form_debug_matrices(A4; show_progress=false)
    res = to_toric_form(A4; show_progress=false)

    @test debug.input_matrix == A4
    @test size(debug.em_matrix) == size(A4)
    @test debug.standard_matrix == res.standard_matrix
    @test debug.standard_blocks == res.standard_blocks
    @test debug.row_transformation == res.row_transformation
    @test debug.row_blocks == res.row_blocks
    @test debug.psi_1_inverse == res.psi_1_inverse
    @test debug.psi_1 === nothing
    @test debug.column_transformation === nothing
    @test debug.product_state_num == res.product_state_num
    @test debug.toric_num == res.toric_num
end

@testset "selected-column Laurent solve helper" begin
    F = GF(2)
    Rxy, (x, y) = laurent_polynomial_ring(F, ["x", "y"])

    A = matrix(Rxy, 3, 3, [
        1 x y;
        0 1 x;
        0 0 1;
    ])
    U = matrix(Rxy, 3, 4, [
        x       y       1       x*y;
        y^-1    x^-1    x + 1   1 + y;
        1       x*y     y^-1    x^-1*y + x;
    ])
    B = A * U

    solved = ToricBuilder._solve_selected_columns(A, B, [1, 3])

    @test size(solved) == size(U)
    @test solved[:, 1:1] == U[:, 1:1]
    @test iszero(solved[:, 2:2])
    @test solved[:, 3:3] == U[:, 3:3]
    @test iszero(solved[:, 4:4])

    empty_solved = ToricBuilder._solve_selected_columns(A, B, Int[])
    @test size(empty_solved) == size(U)
    @test iszero(empty_solved)
end

@testset "initial Pdagger skips known product-state columns" begin
    F = GF(2)
    Rxy, (x, y) = laurent_polynomial_ring(F, ["x", "y"])

    product_state_num = 2
    toric_num = 1
    stab_num = product_state_num + toric_num
    qubit_num = 2 * product_state_num + 2 * toric_num

    Hxdagger = matrix(Rxy, qubit_num, stab_num, [
        1 1 0;
        0 1 1;
        1 0 1;
        0 1 0;
        1 1 1;
        0 0 1;
    ])

    Hxtdagger = matrix(Rxy, qubit_num, stab_num, [
        0 0 0;
        0 0 0;
        1 0 1;
        0 1 1;
        0 0 0;
        0 0 0;
    ])

    expected_Pdagger = matrix(Rxy, stab_num, stab_num, [
        1 0 1;
        0 1 1;
        0 0 0;
    ])

    psi_1_inverse_local = zero_matrix(Rxy, qubit_num, qubit_num)
    psi_1_inverse_local[:, 1] = Hxdagger[:, 1]
    psi_1_inverse_local[:, 2] = Hxdagger[:, 2]
    psi_1_inverse_local[:, product_state_num+1] = Hxdagger[:, 1]
    psi_1_inverse_local[:, product_state_num+2] = Hxdagger[:, 2]

    Pdagger = ToricBuilder._build_initial_pdagger(
        Hxdagger,
        psi_1_inverse_local,
        Hxtdagger,
        product_state_num,
        toric_num,
    )

    @test Pdagger == expected_Pdagger
    @test Hxdagger * Pdagger == psi_1_inverse_local * Hxtdagger
    @test Pdagger[1:product_state_num, 1:product_state_num] ==
          identity_matrix(Rxy, product_state_num)
    @test all(iszero, Pdagger[product_state_num+1:end, 1:product_state_num])

    bad_psi_1_inverse_local = copy(psi_1_inverse_local)
    bad_psi_1_inverse_local[:, product_state_num+1] = zero_matrix(Rxy, qubit_num, 1)
    err = try
        ToricBuilder._build_initial_pdagger(
            Hxdagger,
            bad_psi_1_inverse_local,
            Hxtdagger,
            product_state_num,
            toric_num,
        )
        nothing
    catch e
        e
    end
    @test err isa ErrorException
    @test occursin("product-state columns", sprint(showerror, err))
end

@testset "augmentation ideal telescoping decomposition" begin
    F = GF(2)
    Rxy, (x, y) = laurent_polynomial_ring(F, ["x", "y"])

    samples = [
        zero(Rxy),
        x + 1,
        y^-3 + 1,
        x^2 * y^-1 + x * y + y^-2 + x^-1 * y^2,
    ]

    for d in samples
        a, b = ToricBuilder._decompose_augmentation_entry(d)
        @test (1 + y^-1) * a + (1 + x^-1) * b == d
    end

    err = try
        ToricBuilder._decompose_augmentation_entry(x)
        nothing
    catch e
        e
    end
    @test err isa ErrorException
    @test occursin("augmentation ideal", sprint(showerror, err))
end

@testset "standard Hxt correction solver" begin
    F = GF(2)
    Rxy, (x, y) = laurent_polynomial_ring(F, ["x", "y"])

    product_state_num = 2
    toric_num = 2
    stab_num = product_state_num + toric_num
    qubit_num = 2 * product_state_num + 2 * toric_num

    Hxt = zero_matrix(Rxy, stab_num, qubit_num)
    for i in 1:product_state_num
        Hxt[i, i + product_state_num] = Rxy(1)
    end
    for i in 1:toric_num
        Hxt[product_state_num + i, 2 * product_state_num + 2 * i - 1] = y^-1 + 1
        Hxt[product_state_num + i, 2 * product_state_num + 2 * i] = x^-1 + 1
    end

    P = zero_matrix(Rxy, stab_num, stab_num)
    Pt = zero_matrix(Rxy, stab_num, stab_num)
    P[product_state_num + 1, product_state_num + 1] = x + 1
    P[product_state_num + 1, product_state_num + 2] = x^2 * y^-1 + x * y + y^-2 + x^-1 * y^2
    P[product_state_num + 2, product_state_num + 1] = y^-3 + 1

    Adagger = ToricBuilder._solve_standard_hxt_correction(Hxt, P, Pt, product_state_num, toric_num)
    @test Hxt * Adagger == P + Pt

    bad_outside = copy(P)
    bad_outside[1, 1] = x + 1
    @test_throws ErrorException ToricBuilder._solve_standard_hxt_correction(
        Hxt,
        bad_outside,
        Pt,
        product_state_num,
        toric_num,
    )

    bad_aug = zero_matrix(Rxy, stab_num, stab_num)
    bad_aug[product_state_num + 1, product_state_num + 1] = x
    @test_throws ErrorException ToricBuilder._solve_standard_hxt_correction(
        Hxt,
        bad_aug,
        Pt,
        product_state_num,
        toric_num,
    )
end
