using Test
using Oscar
using ToricBuilder
using ToricBuilder:solve_laurent_linear

@testset "solve_laurent_linear API" begin
    @test isdefined(ToricBuilder, :solve_laurent_linear)
    @test Base.isexported(ToricBuilder, :solve_laurent_linear)
end

@testset "solve_laurent_linear" begin
    F = GF(2)
    R, (x,y,z) = laurent_polynomial_ring(F, ["x","y","z"])
    A = matrix(R, 3, 3, [
    x       y^-1   1;
    z      x^-1   y;
    1      z^-1   x*y^-1
    ])

    u = matrix(R, 3, 1, [
        x*y;
        z^-1;
        1
    ])
    b = A * u

    us = solve_laurent_linear(A, b)
    @test size(us) == size(u)
    @test A * us == b

    # Test case 2: Another exact solution case
    A = matrix(R, 2, 2, [
        x       y;
        z       x*y
    ])
    
    u = matrix(R, 2, 1, [
        y;
        x^-1
    ])
    b = A * u
    
    us = solve_laurent_linear(A, b)
    @test size(us) == size(u)
    @test A * us == b

    # Test case 3: Multiple columns
    A = matrix(R, 3, 3, [
        x       y^-1   1;
        z      x^-1   y;
        1      z^-1   x*y^-1
        ])
    
    u = matrix(R, 3, 2, [
            x*y x*y;
            z^-1 x;
            1 z^-3
        ])
    b = A * u
    us = solve_laurent_linear(A, b)
    @test size(us) == size(u)
    @test A * us == b
    
    # Test case 4: No solution (b not in span of A)
    A = matrix(R, 2, 2, [
        x    y;
        z    x
    ])
    b = matrix(R, 2, 1, [
        x^2 + y;
        z + 1
    ])
    err = try
        solve_laurent_linear(A, b)
        nothing
    catch err
        err
    end
    @test err isa ErrorException
    @test sprint(showerror, err) == "No exact solution exists for A * U = B"
end
