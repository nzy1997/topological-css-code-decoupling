using Test
using ToricBuilder
using ToricBuilder: coarse_graining_with_replace_variable, find_parallelogram_lattice_points
using Oscar

const REPRESENTATION_API_VARIANTS = [
    (
        label="canonical",
        construct_representation=construct_quotient_representation,
        construct_period_table=construct_period_table,
    ),
]

for api in REPRESENTATION_API_VARIANTS
    @testset "representation API ($(api.label))" begin
        @testset "identity case" begin
            F = GF(2)
            Rxy, (x, y) = laurent_polynomial_ring(F, ["x", "y"])
            A = excitation_matrix([1 + x, 1 + y])

            Ax, Ay = api.construct_representation(A)

            I = identity_matrix(base_ring(Ax), nrows(Ax))
            @test Ax == I
            @test Ay == I

            period_table, l = api.construct_period_table(Ax, Ay)
            @test l == 1
            @test Set(period_table) == Set([(1, 0), (0, 1), (1, 1)])
        end

        @testset "z-part case" begin
            F = GF(2)
            Rxy, (x, y) = laurent_polynomial_ring(F, ["x", "y"])
            excitation_map_z = matrix(Rxy, [
                [1 + y, 1 + y, 1 + x, 1 + x],
                [x * y, y, x * y, x],
            ])

            Ax, Ay = api.construct_representation(excitation_map_z)
            I = identity_matrix(base_ring(Ax), nrows(Ax))

            @test Ax * Ay == Ay * Ax
            @test Ax != I
            @test Ay != I

            period_table, l = api.construct_period_table(Ax, Ay)
            @test l == 2
            @test Set(period_table) == Set([(0, 2), (1, 1), (2, 0), (2, 2)])
        end

        @testset "full excitation map case" begin
            F = GF(2)
            Rxy, (x, y) = laurent_polynomial_ring(F, ["x", "y"])
            excitation_map = matrix(Rxy, [
                [1 + y, 1 + y, 1 + x, 1 + x, 0, 0, 0, 0],
                [x * y, y, x * y, x, 0, 0, 0, 0],
                [0, 0, 0, 0, 1 + y, 1 + y, 1 + x, 1 + x],
                [0, 0, 0, 0, x * y, y, x * y, x],
            ])

            Ax, Ay = api.construct_representation(excitation_map)
            I = identity_matrix(base_ring(Ax), nrows(Ax))

            @test size(Ax) == (4, 4)
            @test Ax * Ay == Ay * Ax
            @test Ax != I
            @test Ay != I

            period_table, l = api.construct_period_table(Ax, Ay)
            @test l == 2
            @test Set(period_table) == Set([(0, 2), (1, 1), (2, 0), (2, 2)])
        end

        @testset "pipeline to toric form" begin
            F = GF(2)
            Rxy, (x, y) = laurent_polynomial_ring(F, ["x", "y"])
            excitation_map = matrix(Rxy, [
                [1 + y, 1 + y, 1 + x, 1 + x, 0, 0, 0, 0],
                [x * y, y, x * y, x, 0, 0, 0, 0],
                [0, 0, 0, 0, 1 + y, 1 + y, 1 + x, 1 + x],
                [0, 0, 0, 0, x * y, y, x * y, x],
            ])

            Ax, Ay = api.construct_representation(excitation_map)
            period_table, l = api.construct_period_table(Ax, Ay)

            L_table = zeros(Int, l + 1, l + 1)
            L_table[1, 1] = 1
            for (lx, ly) in period_table
                L_table[lx + 1, ly + 1] = 1
            end

            vertices, area = ToricBuilder.find_minimum_triangle(L_table)
            area *= 2
            @test area > 0

            Ruv, (u, v) = laurent_polynomial_ring(F, ["u", "v"])
            u_rel = x^(vertices[2][1] - 1) * y^(vertices[2][2] - 1)
            v_rel = x^(vertices[3][1] - 1) * y^(vertices[3][2] - 1)
            points = find_parallelogram_lattice_points(vertices...)
            pr = [x^(p[1] - 1) * y^(p[2] - 1) for p in points]

            A_cg = coarse_graining_with_replace_variable(
                excitation_map,
                [u_rel, v_rel],
                Rxy,
                Ruv,
                [x, y],
                [u, v],
                pr,
            )
            res = to_toric_form(A_cg; show_progress=false, compute_inverse=true)
            @test res.row_transformation * A_cg * res.column_transformation == res.standard_matrix
            @test res.toric_num >= 1
        end
    end
end

@testset "anyon translation representation" begin
    F = GF(2)
    Rxy, (x, y) = laurent_polynomial_ring(F, ["x", "y"])
    poly_vec = [1 + x + x*y, 1 + y + x*y]

    rep = construct_anyon_translation_representation(poly_vec)

    @test rep.anyon_count == 2
    @test length(rep.basis) == 2
    @test size(rep.Ax) == (2, 2)
    @test size(rep.Ay) == (2, 2)
    @test rep.Ax * rep.Ay == rep.Ay * rep.Ax

    P = parent(rep.basis[1])
    (xp, yp) = gens(P)
    reduction_ideal = ideal(P, collect(rep.groebner_basis))
    reduction_ordering = getfield(rep.groebner_basis, :ord)

    function coordinate_column(poly, basis, reduction_ideal, reduction_ordering)
        reduced = normal_form(poly, reduction_ideal; ordering=reduction_ordering)
        coords = zero_matrix(F, length(basis), 1)
        for (row, basis_elem) in enumerate(basis)
            basis_coeff = Oscar.coeff(reduced, collect(AbstractAlgebra.exponent_vectors(basis_elem))[1])
            coords[row, 1] = basis_coeff
        end
        return coords
    end

    for (col, basis_elem) in enumerate(rep.basis)
        @test rep.Ax[:, col:col] == coordinate_column(xp * basis_elem, rep.basis, reduction_ideal, reduction_ordering)
        @test rep.Ay[:, col:col] == coordinate_column(yp * basis_elem, rep.basis, reduction_ideal, reduction_ordering)
    end

    laurent_rep = construct_anyon_translation_representation([x^-1 * poly_vec[1], y^-1 * poly_vec[2]])
    @test laurent_rep.anyon_count == rep.anyon_count
    @test size(laurent_rep.Ax) == size(rep.Ax)
    @test size(laurent_rep.Ay) == size(rep.Ay)

    @test_throws ArgumentError construct_anyon_translation_representation([])
    @test_throws ErrorException construct_anyon_translation_representation([1 + x])

    sparse_poly_vec = [1 + x + x^-1*y^3, 1 + y + x^3*y^-1]
    sparse_rep = construct_anyon_translation_representation(sparse_poly_vec)
    @test sparse_rep.anyon_count == 8
    @test size(sparse_rep.Ax) == (8, 8)
    @test size(sparse_rep.Ay) == (8, 8)
    @test sparse_rep.Ax * sparse_rep.Ay == sparse_rep.Ay * sparse_rep.Ax
end
