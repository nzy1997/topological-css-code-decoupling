using ToricBuilder
using Test
using ToricBuilder: ExcitationMatrix, ExcitationMatrixwithF2Matrix, row_switch!
using ToricBuilder: mul_add_row!, row_scale!, column_switch!, mul_add_col!
using Oscar

@testset "ExcitationMatrix Basic" begin
    F = Oscar.GF(2)
    Rxy, (x, y) = laurent_polynomial_ring(F, ["x", "y"])

    A = zero_matrix(Rxy, 2, 4)
    A[1, 1] = 1 + x + x*y
    A[1, 2] = 1 + y + x*y
    A[2, 3] = 1 + y^-1 + x^-1*y^-1
    A[2, 4] = 1 + x^-1 + x^-1*y^-1
    A4 = coarse_graining(A, [x, y], [3, 3])
    E = ExcitationMatrix(A4)
    E2 = ExcitationMatrixwithF2Matrix(E)
    @test E2.em.m == E.m
end

@testset "ExcitationMatrix Operations" begin
    F = Oscar.GF(2)
    Rxy, (x, y) = laurent_polynomial_ring(F, ["x", "y"])
    
    # Create a simple test matrix
    A = zero_matrix(Rxy, 4, 4)
    A[1, 1] = F(1)
    A[1, 2] = x
    A[2, 1] = y
    A[2, 2] = F(1)
    A[3, 3] = F(1)
    A[3, 4] = x + y
    A[4, 3] = x*y
    A[4, 4] = F(1)
    
    @testset "row_switch!" begin
        E = ExcitationMatrix(deepcopy(A))
        Ec = copy(E)
        original_row1 = deepcopy(E.m[1,:])
        original_row2 = deepcopy(E.m[2,:])
        original_trans_row1 = deepcopy(E.row_transformation[1,:])
        original_trans_row2 = deepcopy(E.row_transformation[2,:])
        
        row_switch!(E, 1, 2)
        
        # Check that the rows in the main matrix were swapped
        @test E.m[1,:] == original_row2
        @test E.m[2,:] == original_row1
        
        # Check that the row-transformation matrix was updated in sync
        @test E.row_transformation[1,:] == original_trans_row2
        @test E.row_transformation[2,:] == original_trans_row1
        
        @test E.row_transformation * Ec.m == E.m
        println("✓ row_switch! test passed")
    end
    
    @testset "mul_add_row!" begin
        E = ExcitationMatrix(deepcopy(A))
        original_row1 = deepcopy(E.m[1,:])
        original_row2 = deepcopy(E.m[2,:])
        Ec = copy(E)
        original_trans = deepcopy(E.row_transformation)
        
        # Apply row1 += row2 * x
        mul_add_row!(E, 1, 2, x)
        
        # Check the main matrix
        @test E.m[1,:] == original_row1 + original_row2 .* x
        @test E.m[2,:] == original_row2  # row2 should remain unchanged
        
        # Check the row-transformation matrix
        @test E.row_transformation[1,:] == original_trans[1,:] + original_trans[2,:] .* x
        @test E.row_transformation[2,:] == original_trans[2,:]
        
        @test E.row_transformation * Ec.m == E.m
        println("✓ mul_add_row! test passed")
    end
    
    @testset "row_scale!" begin
        E = ExcitationMatrix(deepcopy(A))
        original_row1 = deepcopy(E.m[1,:])
        original_trans_row1 = deepcopy(E.row_transformation[1,:])
        Ec = copy(E)
        # Apply row1 *= y
        row_scale!(E, 1, y)
        
        # Check the main matrix
        @test E.m[1,:] == original_row1 .* y
        
        # Check the row-transformation matrix
        @test E.row_transformation[1,:] == original_trans_row1 .* y
        
        @test E.row_transformation * Ec.m == E.m
        println("✓ row_scale! test passed")
    end
    
    @testset "column_switch!" begin
        E = ExcitationMatrix(deepcopy(A))
        original_col1 = deepcopy(E.m[:, 1])
        original_col2 = deepcopy(E.m[:, 2])
        original_trans_col1 = deepcopy(E.column_transformation[:, 1])
        original_trans_col2 = deepcopy(E.column_transformation[:, 2])
        Ec = copy(E)
        column_switch!(E, 1, 2)
        
        # Check that the main-matrix columns were swapped, including the paired qubit columns
        @test E.m[:, 1] == original_col2
        @test E.m[:, 2] == original_col1
        
        # Check that the column-transformation matrix was updated in sync
        @test E.column_transformation[:, 1] == original_trans_col2
        @test E.column_transformation[:, 2] == original_trans_col1
        
        @test Ec.m*E.column_transformation == E.m
        println("✓ column_switch! test passed")
    end
    
    @testset "mul_add_col!" begin
        E = ExcitationMatrix(deepcopy(A))
        qubit_num = size(E.m, 2) ÷ 2
        original_col1 = deepcopy(E.m[:, 1])
        original_col2 = deepcopy(E.m[:, 2])
        original_trans = deepcopy(E.column_transformation)
        Ec = copy(E)
        # Apply col1 += col2 * x
        mul_add_col!(E, 1, 2, x)
        
        # Check the main matrix
        @test E.m[:, 1] == original_col1 + original_col2 .* x
        @test E.m[:, 2] == original_col2  # col2 should remain unchanged
        
        # Check the column-transformation matrix
        @test E.column_transformation[:, 1] == original_trans[:, 1] + original_trans[:, 2] .* x
        @test Ec.m*E.column_transformation == E.m
        println("✓ mul_add_col! test passed")
    end
    
    @testset "Transformation Matrix Consistency" begin
        # Test transformation-matrix consistency: R * A * C should equal the transformed matrix
        E = ExcitationMatrix(deepcopy(A))
        original_A = deepcopy(A)
        # Apply a sequence of operations
        row_switch!(E, 1, 2)
        mul_add_row!(E, 1, 2, x)
        column_switch!(E, 1, 2)
        mul_add_col!(E, 1, 2, y)
        
        # Verify that row_transformation * original_A * column_transformation == E.m
        result = E.row_transformation * original_A * E.column_transformation
        @test result == E.m
        println("✓ transformation-matrix consistency test passed")
    end

    @testset "unit multiply-add operations avoid broadcast-scale allocation" begin
        n = 16
        A_dense = zero_matrix(Rxy, n, n)
        for i in 1:n
            for j in 1:n
                A_dense[i, j] = if isodd(i + j)
                    x + y + x*y
                else
                    1 + x^-1*y + x*y^-1
                end
            end
        end

        E_warm = ExcitationMatrix(deepcopy(A_dense))
        mul_add_row!(E_warm, 1, 2, Rxy(1))
        E_warm = ExcitationMatrix(deepcopy(A_dense))
        mul_add_col!(E_warm, 1, 2, Rxy(1))
        GC.gc()

        E_row = ExcitationMatrix(deepcopy(A_dense))
        row_alloc = @allocated mul_add_row!(E_row, 1, 2, Rxy(1))
        @test row_alloc < 55_000
        @test E_row.row_transformation * A_dense == E_row.m

        E_col = ExcitationMatrix(deepcopy(A_dense))
        col_alloc = @allocated mul_add_col!(E_col, 1, 2, Rxy(1))
        @test col_alloc < 110_000
        @test A_dense * E_col.column_transformation == E_col.m
    end
end

@testset "ExcitationMatrixwithF2Matrix Operations" begin
    F = GF(2)
    Rxy, (x, y) = laurent_polynomial_ring(F, ["x", "y"])
    
    # Create a simple test matrix
    A = zero_matrix(Rxy, 4, 4)
    A[1, 1] = F(1)
    A[1, 2] = x
    A[2, 1] = y
    A[2, 2] = F(1)
    A[3, 3] = F(1)
    A[3, 4] = x + y
    A[4, 3] = x*y
    A[4, 4] = F(1)
    
    # Create the expected F2 matrix obtained by setting every variable to 1
    A_F2 = zero_matrix(F, 4, 4)
    A_F2[1, 1] = F(1)
    A_F2[1, 2] = F(1)  # x -> 1
    A_F2[2, 1] = F(1)  # y -> 1
    A_F2[2, 2] = F(1)
    A_F2[3, 3] = F(1)
    A_F2[3, 4] = F(0)  # x + y -> 1 + 1 = 0 in GF(2)
    A_F2[4, 3] = F(1)  # x*y -> 1
    A_F2[4, 4] = F(1)
    
    @testset "Construction and evaluate_all" begin
        E = ExcitationMatrix(deepcopy(A))
        EM = ExcitationMatrixwithF2Matrix(E)
        
        # Verify that m2 is the correct F2 matrix
        @test EM.m == A_F2
        # Verify that m stays unchanged
        @test EM.em.m == A
        
        println("✓ ExcitationMatrixwithF2Matrix construction test passed")
    end
    
    @testset "row_switch!" begin
        E = ExcitationMatrix(deepcopy(A))
        EM = ExcitationMatrixwithF2Matrix(E)
        
        # Store the original data
        original_m_row1 = deepcopy(EM.em.m[1,:])
        original_m_row2 = deepcopy(EM.em.m[2,:])
        original_m2_row1 = deepcopy(EM.m[1,:])
        original_m2_row2 = deepcopy(EM.m[2,:])
        Ec = copy(EM.em.m)
        
        row_switch!(EM, 1, 2)
        
        # Check that the rational-function matrix rows were swapped
        @test EM.em.m[1,:] == original_m_row2
        @test EM.em.m[2,:] == original_m_row1
        
        # Check that the F2 matrix rows were swapped
        @test EM.m[1,:] == original_m2_row2
        @test EM.m[2,:] == original_m2_row1
        
        # Check transformation-matrix consistency
        @test EM.em.row_transformation * Ec == EM.em.m
        
        println("✓ ExcitationMatrixwithF2Matrix row_switch! test passed")
    end
    
    @testset "mul_add_row!" begin
        E = ExcitationMatrix(deepcopy(A))
        EM = ExcitationMatrixwithF2Matrix(E)
        
        original_m_row1 = deepcopy(EM.em.m[1,:])
        original_m_row2 = deepcopy(EM.em.m[2,:])
        original_m2_row1 = deepcopy(EM.m[1,:])
        original_m2_row2 = deepcopy(EM.m[2,:])
        Ec = copy(EM.em.m)
        
        # Apply row1 += row2 * x
        mul_add_row!(EM, 1, 2, F(1))
        
        # Check the rational-function matrix
        @test EM.em.m[1,:] == original_m_row1 + original_m_row2
        @test EM.em.m[2,:] == original_m_row2  # row2 should remain unchanged
        
        # Check the F2 matrix with x evaluated at 1
        @test EM.m[1,:] == original_m2_row1 + original_m2_row2 .* F(1)
        @test EM.m[2,:] == original_m2_row2
        
        # Check transformation-matrix consistency
        @test EM.em.row_transformation * Ec == EM.em.m
        
        println("✓ ExcitationMatrixwithF2Matrix mul_add_row! test passed")
    end
    
    @testset "column_switch!" begin
        E = ExcitationMatrix(deepcopy(A))
        EM = ExcitationMatrixwithF2Matrix(E)
        
        original_m_col1 = deepcopy(EM.em.m[:, 1])
        original_m_col2 = deepcopy(EM.em.m[:, 2])
        original_m2_col1 = deepcopy(EM.m[:, 1])
        original_m2_col2 = deepcopy(EM.m[:, 2])
        Ec = copy(EM.em.m)
        
        column_switch!(EM, 1, 2)
        
        # Check that the rational-function matrix columns were swapped
        @test EM.em.m[:, 1] == original_m_col2
        @test EM.em.m[:, 2] == original_m_col1
        
        # Check that the F2 matrix columns were swapped
        @test EM.m[:, 1] == original_m2_col2
        @test EM.m[:, 2] == original_m2_col1
        
        # Check transformation-matrix consistency
        @test Ec * EM.em.column_transformation == EM.em.m
        
        println("✓ ExcitationMatrixwithF2Matrix column_switch! test passed")
    end
    
    @testset "mul_add_col!" begin
        E = ExcitationMatrix(deepcopy(A))
        EM = ExcitationMatrixwithF2Matrix(E)
        
        original_m_col1 = deepcopy(EM.em.m[:, 1])
        original_m_col2 = deepcopy(EM.em.m[:, 2])
        original_m2_col1 = deepcopy(EM.m[:, 1])
        original_m2_col2 = deepcopy(EM.m[:, 2])
        Ec = copy(EM.em.m)
        
        # Apply col1 += col2 * x
        mul_add_col!(EM, 1, 2, F(1))
        
        # Check the rational-function matrix
        @test EM.em.m[:, 1] == original_m_col1 + original_m_col2
        @test EM.em.m[:, 2] == original_m_col2  # col2 should remain unchanged
        
        # Check the F2 matrix with x evaluated at 1
        @test EM.m[:, 1] == original_m2_col1 + original_m2_col2 .* F(1)
        @test EM.m[:, 2] == original_m2_col2
        
        # Check transformation-matrix consistency
        @test Ec * EM.em.column_transformation == EM.em.m
        
        println("✓ ExcitationMatrixwithF2Matrix mul_add_col! test passed")
    end
    
    @testset "Transformation Consistency" begin
        E = ExcitationMatrix(deepcopy(A))
        EM = ExcitationMatrixwithF2Matrix(E)
        original_A = deepcopy(A)
        original_A_F2 = deepcopy(A_F2)
        
        # Apply a sequence of operations
        row_switch!(EM, 1, 2)
        mul_add_row!(EM, 1, 2, F(1))
        column_switch!(EM, 1, 2)
        mul_add_col!(EM, 1, 2, F(1))
        
        # Verify transformation consistency for the rational-function matrix
        result = EM.em.row_transformation * original_A * EM.em.column_transformation
        @test result == EM.em.m
        
        # Verify transformation consistency for the F2 matrix
        result_F2 = EM.em.row_transformation * original_A_F2 * EM.em.column_transformation
        result_F2 = map(result_F2) do elem
            evaluate(elem, fill(F(1),n_variables(EM.em.m.base_ring)))
        end
        @test result_F2 == EM.m
        
        println("✓ ExcitationMatrixwithF2Matrix transformation consistency test passed")
    end
    
    @testset "m2 consistency with m evaluation" begin
        E = ExcitationMatrix(deepcopy(A))
        EM = ExcitationMatrixwithF2Matrix(E)
        
        # Apply a sequence of operations
        row_switch!(EM, 1, 2)
        mul_add_row!(EM, 3, 1, F(1))
        column_switch!(EM, 2, 1)
        mul_add_col!(EM, 3, 1, F(1))
        
        # Verify that m2 equals the result of evaluating every variable in m at 1
        one_F2 = F(1)
        evaluated_m = map(EM.em.m) do elem
            evaluate(elem, fill(one_F2,n_variables(EM.em.m.base_ring)))
        end
        @test evaluated_m == EM.m
        
        println("✓ m2 consistency with evaluating m test passed")
    end
end

@testset "maximum_deg_term" begin
    F = Oscar.GF(2)
    R, (x, y) = laurent_polynomial_ring(F, ["x", "y"])
    
    # Test 1: Simple case with different degrees
    # x^2 + x^(-1): degree = 2 - (-1) = 3
    # y: degree = 0
    # x: degree = 0  
    # x*y^2 + x*y^(-1): x ∈ [1,1] (range=0), y ∈ [-1,2] (range=3), degree = 3
    A = matrix(R, 2, 2, [
        x^2 + x^(-1),  y,
        x,             x*y^2 + x*y^(-1)
    ])
    
    max_ele,max_degree = ToricBuilder.maximum_deg_term(A)
    # Should be either x^2 + x^(-1) or x*y^2 + x*y^(-1), both have degree 3
    @test max_ele == x^2 + x^(-1) || max_ele == x*y^2 + x*y^(-1)
    
    println("✓ Test 1: Found the maximum-degree element (degree 3)")
    
    # Test 2: Matrix with element of higher degree
    # x^3 + x^(-2): degree = 3 - (-2) = 5
    B = matrix(R, 2, 2, [
        x^3 + x^(-2),  y,
        x,             x*y
    ])
    
    max_elem,max_degree = ToricBuilder.maximum_deg_term(B)
    @test max_elem == x^3 + x^(-2)
    
    println("✓ Test 2: Found the maximum-degree element (degree 5)")
    
    # Test 3: Check degree calculation for multi-variable polynomial
    # x^2*y^3 + x^(-1)*y^(-2): x ∈ [-1,2] (range=3), y ∈ [-2,3] (range=5), degree = 8
    p = x^2 * y^3 + x^(-1) * y^(-2)
    deg = ToricBuilder.polynomial_degree(p)
    @test deg == 8
    
    println("✓ Test 3: Multivariable polynomial degree computed correctly (degree 8)")
end

@testset "polynomial_degree_monomial" begin
    F = Oscar.GF(2)
    R, (x, y) = laurent_polynomial_ring(F, ["x", "y"])

    p = x^2 * y^3 + x^(-1) * y^(-2) + x^(-6)*y^3
    deg = ToricBuilder.polynomial_degree_monomial(p)
    @test deg == 9

    @test ToricBuilder.polynomial_degree_monomial(zero(R)) == 0
    @test ToricBuilder.polynomial_degree_monomial(R(1)) == 0
end

@testset "polynomial_term_count" begin
    F = Oscar.GF(2)
    R, (x, y) = laurent_polynomial_ring(F, ["x", "y"])

    p = x^2 + x^(-1) + y + x^(-6)*y^3
    @test ToricBuilder.polynomial_term_count(p) == 4
    @test ToricBuilder.polynomial_term_count(zero(R)) == 0
    @test ToricBuilder.polynomial_term_count(R(1)) == 1
end

@testset "max_column_monomial_count" begin
    F = Oscar.GF(2)
    R, (x, y) = laurent_polynomial_ring(F, ["x", "y"])

    A = matrix(R, 2, 3, [
        x + y,        0,             x*y + x^-1,
        x + 1,        y^-1,          x^2 + x^(-1) + y + x^(-6)*y^3
    ])

    @test ToricBuilder.max_column_monomial_count(A) == 6
end
