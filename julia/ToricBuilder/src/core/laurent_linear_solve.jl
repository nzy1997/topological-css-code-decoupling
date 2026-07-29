function _setup_module(A, f)
    A_quo = map_entries(f, A)
    RQ = base_ring(A_quo)
    P = base_ring(RQ)
    I = modulus(RQ)
    A_poly = map_entries(x -> P(Oscar.lift(x)), A_quo)

    F = free_module(P, nrows(A_poly))
    n = ncols(A_poly)
    generators = [F(collect(A_poly[:, j])) for j in 1:n]

    m = nrows(A_poly)
    for g in gens(I)
        for i in 1:m
            v = [zero(P) for _ in 1:m]
            v[i] = g
            push!(generators, F(v))
        end
    end

    M, inc = sub(F, generators)
    return (F=F, M=M, inc=inc, n=n, P=P, RQ=RQ)
end

function _solve_single_column(b_col, module_data, f, R)
    (; F, M, n, P, RQ) = module_data

    b_quo = map_entries(f, b_col)
    b_poly = map_entries(x -> P(Oscar.lift(x)), b_quo)
    b_vec = F(collect(b_poly[:, 1]))
    coeffs_vec = Oscar.coordinates(b_vec, M)

    u_coeffs_poly = [coeffs_vec[i] for i in 1:n]
    u_coeffs_quo = [RQ(u_coeffs_poly[i]) for i in 1:n]
    u_laurent = [preimage(f, u_coeffs_quo[i]) for i in 1:n]
    return matrix(R, length(u_laurent), 1, u_laurent)
end

function _validate_laurent_linear_inputs(A, B)
    R = base_ring(A)

    if base_ring(B) != R
        throw(ArgumentError("Matrices A and B must be over the same ring"))
    end

    m, _ = size(A)
    mb, _ = size(B)
    if mb != m
        throw(DimensionMismatch("Number of rows in A ($m) must match number of rows in B ($mb)"))
    end

    return R
end

function _solve_laurent_linear_native(A, B)
    solvable, solution = can_solve_with_solution(A, B; side=:right)
    if solvable
        return solution
    end

    throw(ErrorException("No exact solution exists for A * U = B"))
end

function _is_native_solver_capability_error(err, backtrace)
    if !(err isa MethodError)
        return false
    end

    missing_name = string(err.f)
    if missing_name != "gcdxx" && missing_name != "annihilator"
        return false
    end

    frames = stacktrace(backtrace)
    return any(frame -> startswith(String(frame.func), "can_solve_with_solution"), frames)
end

"""
    solve_laurent_linear(A, B)

Solve the Laurent-polynomial linear system `A * U = B` exactly.

The implementation prefers Oscar's native right-side linear solver. If the
installed Oscar/AbstractAlgebra stack does not support Laurent-polynomial
matrices for that path, it falls back to a module-based exact solver while
preserving the same public API.

# Arguments
- `A`: `m x n` matrix over a Laurent polynomial ring
- `B`: `m x k` matrix over the same Laurent polynomial ring

# Returns
- `U`: `n x k` solution matrix such that `A * U == B`

# Throws
- `ErrorException` if no exact solution exists
- `ArgumentError` or `DimensionMismatch` for invalid inputs
"""
function solve_laurent_linear(A, B)
    _validate_laurent_linear_inputs(A, B)

    try
        return _solve_laurent_linear_native(A, B)
    catch err
        if _is_native_solver_capability_error(err, catch_backtrace())
            return _solve_laurent_linear_legacy(A, B)
        end
        rethrow()
    end
end

function _solve_laurent_linear_legacy(A, B)
    R = _validate_laurent_linear_inputs(A, B)
    f = Oscar._polyringquo(R)
    _, kb = size(B)
    module_data = _setup_module(A, f)

    try
        if kb == 1
            return _solve_single_column(B, module_data, f, R)
        end

        solutions = MatrixElem[]
        for j in 1:kb
            push!(solutions, _solve_single_column(B[:, j:j], module_data, f, R))
        end

        return hcat(solutions...)
    catch err
        if err isa ErrorException && occursin("not liftable to the given generating system", sprint(showerror, err))
            throw(ErrorException("No exact solution exists for A * U = B"))
        end
        rethrow()
    end
end
