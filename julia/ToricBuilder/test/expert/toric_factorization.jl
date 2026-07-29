using ToricBuilder
using Test
using Oscar

@testset "factor_toric_block API" begin
    has_binding = isdefined(ToricBuilder, :factor_toric_block)
    @test has_binding
    @test Base.isexported(ToricBuilder, :factor_toric_block)
    if has_binding
        @test getfield(ToricBuilder, :factor_toric_block) isa Function
    end
end

@testset "toric factorization" begin
    if !isdefined(ToricBuilder, :factor_toric_block)
        @test_broken false
    else
        factor_toric_block = getfield(ToricBuilder, :factor_toric_block)
        F = GF(2)
        Rxy, (x, y) = laurent_polynomial_ring(F, ["x", "y"])
        res = build_toric_form([1 + x + x^-1 * y, 1 + y + x * y]; show_progress=false, compute_inverse=true)
        column_transform_inv, row_transform_inv, a, b = factor_toric_block(3, x, y, Rxy)
        A2 = coarse_graining(excitation_matrix([1 + x, 1 + y]), [x, y], [a, b])
        @test size(column_transform_inv, 1) == size(A2, 2)
        @test size(column_transform_inv, 2) == size(A2, 2)
        @test size(row_transform_inv, 1) == size(A2, 1)
        @test size(row_transform_inv, 2) == size(A2, 1)

        for toric_list in ([2, 1, 1], [4, 0, 0], [0, 1, 3], [0, 0, 1])
            res2 = split_product_state(res, toric_list)
            R = res.standard_matrix.base_ring
            u, v = gens(R)
            @test res2.row_transformation * res.input_matrix * res2.column_transformation == res2.result_matrix
            for i in 1:length(toric_list)
                a, b = res2.toric_cg_size[i]
                A2 = coarse_graining(excitation_matrix([1 + u, 1 + v]), [u, v], [a, b])
                stab_pos = [res2.zstab_pos_list[i]; res2.xstab_pos_list[i]]
                qubit_pos = res2.qubit_pos_list[i]
                qubit_pos_all = [qubit_pos; qubit_pos .+ size(res2.result_matrix, 2) ÷ 2]
                @test res2.result_matrix[stab_pos, qubit_pos_all] == A2
            end
        end

        tuple_res = split_product_state(res, (0, 0, 0))
        @test tuple_res.row_transformation * res.input_matrix * tuple_res.column_transformation == tuple_res.result_matrix

        default_res = build_toric_form([1 + x + x^-1 * y, 1 + y + x * y]; show_progress=false)
        err = try
            split_product_state(default_res, fill(0, default_res.toric_num))
            nothing
        catch e
            e
        end
        @test err isa ArgumentError
        @test occursin("compute_inverse=true", sprint(showerror, err))

        @test_throws ArgumentError split_product_state(res, [-1, 0, 0])
        @test_throws ArgumentError split_product_state(res, [0, 0])
        @test_throws ArgumentError split_product_state(res, [res.product_state_num + 1, 0, 0])
    end
end
