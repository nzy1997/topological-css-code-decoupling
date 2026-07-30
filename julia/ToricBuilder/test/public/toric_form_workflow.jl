using Test
using ToricBuilder
using Oscar

@testset "toric form workflow" begin
    F = GF(2)
    Rxy, (x, y) = laurent_polynomial_ring(F, ["x", "y"])

    poly_vec = [1 + x + x * y, 1 + y + x * y]
    A = excitation_matrix(poly_vec)
    l = find_L([x, y], ideal(poly_vec))
    A_cg = coarse_grain(A, [x, y], [l, l])
    result = build_toric_form(poly_vec; show_progress=false)

    @test size(A_cg) == (2 * l^2, 4 * l^2)
    @test hasproperty(result, :row_transformation)
    @test hasproperty(result, :input_matrix)
    @test hasproperty(result, :standard_matrix)
    @test hasproperty(result, :psi_1_inverse)
    @test result.psi_1 === nothing
    @test !hasproperty(result, Symbol("ph", "i_1"))
    @test !hasproperty(result, Symbol("ph", "i_1_inv"))
    @test result.column_transformation === nothing
    @test check_result(result, result.input_matrix)
    @test result.product_state_num == 1
    @test result.toric_num == 2

    full_result = build_toric_form(poly_vec; show_progress=false, compute_inverse=true)
    @test !isnothing(full_result.psi_1)
    identity = identity_matrix(base_ring(full_result.psi_1_inverse), size(full_result.psi_1_inverse, 1))
    @test full_result.psi_1_inverse * full_result.psi_1 == identity
    @test full_result.psi_1 * full_result.psi_1_inverse == identity
    @test !isnothing(full_result.column_transformation)
    @test check_result(full_result, full_result.input_matrix; require_inverse=true)
    @test full_result.l == l

    prep = ToricBuilder._prepare_transfer_to_toric_form_data(poly_vec)
    @test prep.ideal_triangle === nothing
    @test prep.L_table == make_L_table(poly_vec, prep.l)
end

@testset "excitation matrix construction" begin
    F = GF(2)
    Rxy, (x, y) = laurent_polynomial_ring(F, ["x", "y"])

    A = excitation_matrix([1 + x + x * y, 1 + y + x * y])
    @test A.entries == [
        1 + x + x * y 1 + y + x * y 0 0
        0 0 1 + y^-1 + x^-1 * y^-1 1 + x^-1 + x^-1 * y^-1
    ]
end

@testset "toric form workflow early returns" begin
    F = GF(2)
    Rxy, (x, y) = laurent_polynomial_ring(F, ["x", "y"])
    poly_vec = [1 + x + x * y, 1 + y + x * y]

    short_result = build_toric_form(poly_vec; max_l=2, show_progress=false)
    medium_area_result = build_toric_form(poly_vec; max_area=4, show_progress=false)
    area_limited_result = build_toric_form(poly_vec; max_area=1, show_progress=false)

    @test hasproperty(short_result, :solving_time)
    @test !hasproperty(short_result, :row_transformation)
    @test !hasproperty(short_result, :u_rel)

    @test hasproperty(medium_area_result, :row_transformation)
    @test medium_area_result.area == 3

    @test hasproperty(area_limited_result, :solving_time)
    @test hasproperty(area_limited_result, :area)
    @test hasproperty(area_limited_result, :u_rel)
    @test hasproperty(area_limited_result, :v_rel)
    @test !hasproperty(area_limited_result, :row_transformation)
end

@testset "project to finite code" begin
    F = GF(2)
    Rxy, (x, y) = laurent_polynomial_ring(F, ["x", "y"])
    A = excitation_matrix([1 + x, 1 + y])
    A2 = project_to_finite_code(A, 2, 2)

    @test A2 == [
        1 0 1 0 1 1 0 0 0 0 0 0 0 0 0 0
        0 1 0 1 1 1 0 0 0 0 0 0 0 0 0 0
        1 0 1 0 0 0 1 1 0 0 0 0 0 0 0 0
        0 1 0 1 0 0 1 1 0 0 0 0 0 0 0 0
        0 0 0 0 0 0 0 0 1 1 0 0 1 0 1 0
        0 0 0 0 0 0 0 0 1 1 0 0 0 1 0 1
        0 0 0 0 0 0 0 0 0 0 1 1 1 0 1 0
        0 0 0 0 0 0 0 0 0 0 1 1 0 1 0 1
    ]
end
