function coarse_graining(A, variable_list, enlarge_num::Vector{Int})
    @assert length(variable_list) == length(enlarge_num) "Variable list and expansion dimension list must have the same length"
    
    rows, cols = size(A)
    
    # Calculate expansion factor (product of all enlarge_num)
    expansion = prod(enlarge_num)
    # Create result matrix
    result = zero_matrix(A.base_ring, rows * expansion, cols * expansion)
    
    # Process each entry in the original matrix A
    for i in 1:rows
        for j in 1:cols
            entry = A[i, j]
            
            # Skip zero entries for efficiency
            if iszero(entry)
                continue
            end
            
            # Use enlarge_polynomial to expand this entry into a block matrix
            try
                block = enlarge_polynomial(A.base_ring, entry, variable_list, enlarge_num)
                
                # Place the block in the result matrix
                for bi in 1:expansion
                    for bj in 1:expansion
                        result[(i-1)*expansion + bi, (j-1)*expansion + bj] = block[bi, bj]
                    end
                end
            catch e
                println("Error processing entry at position ($i, $j): ", entry)
                println("Error message: ", e)
                rethrow(e)
            end
        end
    end
    
    return result
end

coarse_grain(args...; kwargs...) = coarse_graining(args...; kwargs...)

"""
    enlarge_monomial(R, variable_list, var_powers::Vector{Int}, enlarge_num::Vector{Int})

Replace each variable in a monomial with a matrix and compute the tensor product.

# Arguments
- `R`: Ring for matrix elements (outermost ring)
- `variable_list`: List of variables, e.g., [x, y]
- `var_powers`: Power of each variable, e.g., [2, 1] represents x^2*y
- `enlarge_num`: Matrix dimension for each variable, e.g., [2, 3]

# Returns
- Result matrix (tensor product of all variable matrices)

# Example
For x^2*y with [x, y], [2, 1], [2, 3]:
- x corresponds to 2, construct matrix [0 x; 1 0], then square it
- y corresponds to 3, construct matrix [0 0 y; 1 0 0; 0 1 0]
- Return the tensor product of the two matrices
"""
function enlarge_monomial(R, variable_list, var_powers::Vector{Int}, enlarge_num::Vector{Int})
    @assert length(variable_list) == length(var_powers) == length(enlarge_num) "All arguments must have the same length"
    
    # Construct corresponding matrix for each variable
    matrices = []
    for (i, var) in enumerate(variable_list)
        n = enlarge_num[i]
        power = var_powers[i]
        
        # Construct base matrix: subdiagonal is 1, top-right corner is the variable
        # [0 0 ... 0 var]
        # [1 0 ... 0 0  ]
        # [0 1 ... 0 0  ]
        # [.............]
        # [0 0 ... 1 0  ]
        base_mat = zero_matrix(R, n, n)
        for j in 1:n-1
            base_mat[j+1, j] = R(1)  # Subdiagonal is 1
        end
        base_mat[1, n] = var  # Top-right corner is the variable
        
        # Calculate matrix power based on exponent
        if power > 0
            mat = base_mat^power
        elseif power == 0
            mat = identity_matrix(R, n)
        else
            # Handle negative powers
            mat = (var^-1 .* base_mat^((n-1)))^(-power)
        end
        
        push!(matrices, mat)
    end
    
    # Compute tensor product (Kronecker product)
    result = matrices[1]
    for i in 2:length(matrices)
        result = kronecker_product(result, matrices[i])
    end
    
    return result
end

"""
    enlarge_polynomial(R, poly, variable_list, enlarge_num::Vector{Int})

Expand each monomial in the polynomial into a matrix and sum them.

# Arguments
- `R`: Ring for matrix elements (outermost ring)
- `poly`: Polynomial (in nested rational function field)
- `variable_list`: List of variables, e.g., [x, y]
- `enlarge_num`: Matrix dimension for each variable, e.g., [2, 3]

# Returns
- Result matrix (weighted sum of all monomial matrices)

# Example
For polynomial x²y + xy² + 1 with enlarge_num=[2, 3]:
- x²y corresponds to one matrix
- xy² corresponds to one matrix
- 1 corresponds to one matrix
- Return the sum of the three matrices
"""
function enlarge_polynomial(R, poly, variable_list, enlarge_num::Vector{Int})
    @assert length(variable_list) == length(enlarge_num) "Arguments must have the same length"
    
    # Calculate dimension of result matrix
    total_dim = prod(enlarge_num)
    
    # Initialize zero matrix
    result = zero_matrix(R, total_dim, total_dim)
    
    # Iterate through all terms of the polynomial (including cases with denominators)
    # For nested rational function fields, recursively extract monomials from each layer
    # Note: Start from the outermost variable (last one), but collect powers correctly
    # Initial denominator offset is all zeros (indicating no denominator)
    denom_offset = zeros(Int, length(variable_list))
    extract_and_add_terms!(result, poly, variable_list, enlarge_num, R, Int[], length(variable_list), denom_offset)
    
    return result
end

"""
    extract_and_add_terms!(result, poly, variable_list, enlarge_num, R, 
                          current_powers, var_idx, denom_offset)

Extract terms from Laurent polynomial and add them to result matrix.

Designed specifically for Laurent polynomial rings created with laurent_polynomial_ring().
Supports native negative powers without requiring fraction notation.

# Supported types
- LaurentMPolyWrap: Multivariate Laurent polynomials (non-recursive extraction)
- LaurentPolyWrap: Univariate Laurent polynomials (recursive extraction)

# Examples
- x^2 + x^(-1) + 1 (native Laurent polynomial notation)
- x*y^(-1) + x^(-1)*y (multivariate with negative powers)
"""
function extract_and_add_terms!(result, poly, variable_list, enlarge_num, R,
                                current_powers::Vector{Int}, var_idx::Int, denom_offset::Vector{Int})
    # For multivariate Laurent polynomials (LaurentMPolyWrap), directly extract all terms
    if poly isa AbstractAlgebra.Generic.LaurentMPolyWrap
        inner_poly = poly.mpoly  # Note: field is 'mpoly' not 'poly' for multivariate
        min_degs = poly.mindegs
        
        # Iterate through all terms
        for (coef, exp_vec) in zip(Oscar.coefficients(inner_poly), exponent_vectors(inner_poly))
            if !iszero(coef)
                # Actual exponents = exponent_vector + minimum_degrees
                actual_exps = [exp_vec[i] + min_degs[i] for i in 1:length(variable_list)]
                
                # Compute matrix for this monomial
                mat = enlarge_monomial(R, variable_list, actual_exps, enlarge_num)
                result .+= coef .* mat
            end
        end
        return
    end
    
    # For univariate Laurent polynomials (LaurentPolyWrap), use recursive approach
    if poly isa AbstractAlgebra.Generic.LaurentPolyWrap
        # Base case: all variables processed
        if var_idx < 1
            if !iszero(poly)
                mat = enlarge_monomial(R, variable_list, current_powers, enlarge_num)
                result .+= poly .* mat
            end
            return
        end
        
        # Extract terms from Laurent polynomial
        inner_pol = poly.poly
        mindeg = poly.mindeg
        
        for d in 0:(length(inner_pol)-1)
            c = coeff(inner_pol, d)
            if !iszero(c)
                actual_exp = d + mindeg
                new_powers = copy(current_powers)
                insert!(new_powers, 1, actual_exp)
                extract_and_add_terms!(result, c, variable_list, enlarge_num, R,
                                     new_powers, var_idx - 1, denom_offset)
            end
        end
        return
    end
    
    # Base case: coefficient is a scalar
    if var_idx < 1
        if !iszero(poly)
            mat = enlarge_monomial(R, variable_list, current_powers, enlarge_num)
            result .+= poly .* mat
        end
        return
    end
    
    # Unsupported type
    error("Only Laurent polynomial rings (created with laurent_polynomial_ring) are supported. Got type: $(typeof(poly))")
end

"""
    kronecker_product(A::MatElem, B::MatElem)

Compute the Kronecker product (tensor product) of two matrices.

If A is an m×n matrix and B is a p×q matrix,
the result is an mp×nq matrix.
"""
function kronecker_product(A::MatElem, B::MatElem)
    m, n = size(A)
    p, q = size(B)
    
    R = base_ring(A)
    result = zero_matrix(R, m*p, n*q)
    result.entries .= kron(A.entries, B.entries)
    return result
end

function find_L(variable_list, I;max_int = 1000)
    variable_num = length(variable_list)
    l_list = Int[]
    for i in 1:variable_num
        for j in 1:max_int
            if (variable_list[i]^j+1) in I
                push!(l_list, j)
                break
            end
        end
    end
    if isempty(l_list)
        return nothing
    else
        return lcm(l_list)
    end
end

function make_L_table(variable_list, I,L)
    variable_num = length(variable_list)
    l_table = zeros(Int, fill(L+1,variable_num)...)
    for ind in CartesianIndices(l_table)
        a = mapreduce(i -> variable_list[i]^(ind.I[i]-1), *, 1:variable_num)
        if a+1 in I
            l_table[ind.I...] =1
        end
    end
    return l_table
end

function _action_order(A; max_int::Int=1000)
    max_int < 1 && return nothing

    nrows(A) == ncols(A) || throw(DimensionMismatch("action matrix must be square"))
    I = identity_matrix(base_ring(A), nrows(A))
    power = I
    for j in 1:max_int
        power = power * A
        power == I && return j
    end
    return nothing
end

function _make_L_table_from_actions(Ax, Ay, L::Int)
    L >= 0 || throw(ArgumentError("L must be nonnegative"))
    base_ring(Ax) == base_ring(Ay) || throw(ArgumentError("Ax and Ay must be over the same base ring"))
    nrows(Ax) == ncols(Ax) || throw(DimensionMismatch("Ax must be square"))
    nrows(Ay) == ncols(Ay) || throw(DimensionMismatch("Ay must be square"))
    size(Ax) == size(Ay) || throw(DimensionMismatch("Ax and Ay must have the same size"))

    x_powers = Matrix{typeof(Ax)}(undef, L + 1, 1)
    y_powers = Matrix{typeof(Ay)}(undef, L + 1, 1)

    I = identity_matrix(base_ring(Ax), nrows(Ax))
    x_powers[1] = I
    y_powers[1] = I
    for i in 1:L
        x_powers[i + 1] = x_powers[i] * Ax
        y_powers[i + 1] = y_powers[i] * Ay
    end

    l_table = zeros(Int, L + 1, L + 1)
    for ix in 0:L
        for iy in 0:L
            if x_powers[ix + 1] * y_powers[iy + 1] == I
                l_table[ix + 1, iy + 1] = 1
            end
        end
    end
    return l_table
end

function _matrix_period_data(poly_vec::Vector; max_int::Int=1000)
    rep = construct_anyon_translation_representation(poly_vec)
    lx = _action_order(rep.Ax; max_int=max_int)
    ly = _action_order(rep.Ay; max_int=max_int)
    (isnothing(lx) || isnothing(ly)) && return (; rep=rep, l=nothing)
    return (; rep=rep, l=lcm(lx, ly))
end

function find_L(poly_vec::Vector; max_int = 1000)
    data = _matrix_period_data(poly_vec; max_int=Int(max_int))
    isnothing(data.l) && return nothing
    data.l > max_int && return nothing
    return data.l
end

function make_L_table(poly_vec::Vector, L)
    data = _matrix_period_data(poly_vec; max_int=max(Int(L), 1))
    return _make_L_table_from_actions(data.rep.Ax, data.rep.Ay, Int(L))
end

function find_minimum_triangle_periods(poly_vec::Vector; max_int = 1000)
    data = _matrix_period_data(poly_vec; max_int=Int(max_int))
    isnothing(data.l) && return nothing
    data.l > max_int && return nothing

    vertices, area = find_minimum_triangle(_make_L_table_from_actions(data.rep.Ax, data.rep.Ay, data.l))
    point1 = (vertices[2][1] - 1, vertices[2][2] - 1)
    point2 = (vertices[3][1] - 1, vertices[3][2] - 1)
    return point1, point2, Int(2 * area)
end


function replace_variable(a, new_var_rels, R_old, R_new, old_vars, new_vars; allowed_remainders=nothing)
    n_old = length(old_vars)
    n_new = length(new_vars)
    
    # Use Oscar's polynomial utilities to extract exponents
    function get_exponents_oscar(poly, vars)
        if iszero(poly)
            return nothing
        end
        
        # For a monomial (coefficient times a product of variable powers)
        # If the polynomial has a single term, read off the exponent vector directly
        exp_vecs = AbstractAlgebra.exponent_vectors(poly)
        exp_vec_array = collect(exp_vecs)
        
        if length(exp_vec_array) != 1
            error("The input polynomial must be a monomial, but it contains $(length(exp_vec_array)) terms")
        end
        
        exps = exp_vec_array[1]
        
        # Build the exponent dictionary
        exp_dict = Dict()
        for (i, var) in enumerate(vars)
            exp_dict[var] = exps[i]
        end
        
        return exp_dict
    end
    
    # Extract the exponents of each monomial
    a_exp = get_exponents_oscar(a, old_vars)
    new_var_exps = [get_exponents_oscar(rel, old_vars) for rel in new_var_rels]
    
    if isnothing(a_exp)
        return R_new(0), R_old(0)
    end
    
    # If an allowed-remainder list is provided, convert it to exponent form
    allowed_remainder_exps = nothing
    if !isnothing(allowed_remainders)
        allowed_remainder_exps = [get_exponents_oscar(rem, old_vars) for rem in allowed_remainders]
    end
    
    # Extract the target exponent vector
    target_exps = [a_exp[var] for var in old_vars]
    
    # Build the coefficient matrix, with one row per new variable and one column per old variable
    coeff_matrix = zeros(Int, n_new, n_old)
    for i in 1:n_new
        for j in 1:n_old
            coeff_matrix[i, j] = new_var_exps[i][old_vars[j]]
        end
    end
    
    # Search for the best coefficient vector
    best_coeffs = nothing
    best_remainder_exps = nothing
    best_score = nothing
    
    # More efficient approach: iterate over allowed remainders and solve the linear system directly
    # This reduces the complexity from O((2*max_range)^n_new) to O(num_remainders)
    
    if isnothing(allowed_remainder_exps)
        # If no allowed-remainder list is provided, fall back to enumeration for backward compatibility
        max_range = maximum(abs.(target_exps)) + 10
        
        function enumerate_coeffs(depth, current_coeffs)
            if depth > n_new
                remainder_exps = copy(target_exps)
                for i in 1:n_new
                    for j in 1:n_old
                        remainder_exps[j] -= coeff_matrix[i, j] * current_coeffs[i]
                    end
                end
                
                if !all(remainder_exps .>= 0)
                    return
                end
                
                coeff_sum = sum(abs.(current_coeffs))
                remainder_deg = sum(abs.(remainder_exps))
                current_score = (remainder_deg, coeff_sum)
                
                if isnothing(best_score) || current_score < best_score
                    best_coeffs = copy(current_coeffs)
                    best_remainder_exps = copy(remainder_exps)
                    best_score = current_score
                end
                return
            end
            
            for coeff in -max_range:max_range
                current_coeffs[depth] = coeff
                enumerate_coeffs(depth + 1, current_coeffs)
            end
        end
        
        enumerate_coeffs(1, zeros(Int, n_new))
    else
        # Use linear algebra: solve one linear system for each allowed remainder
        # We solve M^T * c = target_exps - remainder_exps
        # Here M is coeff_matrix (n_new x n_old), so M^T has shape (n_old x n_new)
        
        # Transpose the coefficient matrix so that M^T[j, i] = M[i, j]
        M_T = zeros(Float64, n_old, n_new)
        for i in 1:n_new
            for j in 1:n_old
                M_T[j, i] = coeff_matrix[i, j]
            end
        end
        
        for allowed_exp in allowed_remainder_exps
            remainder_vec = [allowed_exp[var] for var in old_vars]
            
            # Right-hand side: target_exps - remainder_vec
            rhs = Float64.(target_exps .- remainder_vec)
            
            # Solve the linear system M^T * c = rhs
            # Use least squares to handle overdetermined or underdetermined systems
            try
                if n_new == n_old
                    # Square case: solve directly
                    c_solution = M_T \ rhs
                else
                    # Non-square case: use the pseudoinverse
                    c_solution = pinv(M_T) * rhs
                end
                
                # Check that the solution entries are integers, up to small numerical error
                tolerance = 1e-6
                c_rounded = round.(Int, c_solution)
                
                if all(abs.(c_solution .- c_rounded) .< tolerance)
                    # Verify the solution
                    remainder_exps_check = copy(target_exps)
                    for i in 1:n_new
                        for j in 1:n_old
                            remainder_exps_check[j] -= coeff_matrix[i, j] * c_rounded[i]
                        end
                    end
                    
                    # Confirm that the remainder matches
                    if all(remainder_exps_check .== remainder_vec)
                        # Scoring rule: prefer smaller remainder degree, then smaller coefficient sum
                        coeff_sum = sum(abs.(c_rounded))
                        remainder_deg = sum(abs.(remainder_vec))
                        current_score = (remainder_deg, coeff_sum)
                        
                        if isnothing(best_score) || current_score < best_score
                            best_coeffs = copy(c_rounded)
                            best_remainder_exps = copy(remainder_vec)
                            best_score = current_score
                        end
                    end
                end
            catch e
                # If solving fails, skip this remainder
                continue
            end
        end
    end
    
    # Check whether a valid solution was found
    if isnothing(best_coeffs)
        return (nothing, nothing)
    end
    
    # Construct the result
    new_var_term = R_new(1)
    for i in 1:n_new
        new_var_term *= new_vars[i]^best_coeffs[i]
    end
    
    remainder = R_old(1)
    for i in 1:n_old
        remainder *= old_vars[i]^best_remainder_exps[i]
    end
    
    return (new_var_term, remainder)
end

function coarse_graining_with_replace_variable(A, new_var_rels, R_old, R_new, old_vars, new_vars, allowed_remainders)
    @assert length(old_vars) == length(new_vars) "Variable list and expansion dimension list must have the same length"
    
    num_remainders = length(allowed_remainders)
    rows, cols = size(A)
    
    result = zero_matrix(R_new, rows * num_remainders, cols * num_remainders)
    
    # Process each entry in the original matrix A
    for i in 1:rows
        for j in 1:cols
            entry = A[i, j]
            # @show entry
            # Skip zero entries for efficiency
            if iszero(entry)
                continue
            end
            
            for k in 1:num_remainders
                entry2 = entry * allowed_remainders[k]
                # @show entry2
                for (coef, exp_vec) in zip(Oscar.coefficients(entry2), exponent_vectors(entry2))
                    if !iszero(coef)
                        mono = mapreduce(i -> old_vars[i]^(exp_vec[i]), *, 1:length(old_vars))
                        # @show mono
                        a,remainder = replace_variable(mono, new_var_rels, R_old, R_new, old_vars, new_vars; allowed_remainders=allowed_remainders)
                        # @show a remainder
                        # @show a,remainder
                        re = findfirst(==(remainder), allowed_remainders)
                        # @show re
                        result[(i-1)*num_remainders + k, (j-1)*num_remainders + re] += coef .* a
                    end
                end
            end
            
        end
    end
    
    return result
end

function eliminate_minus_deg(A; eliminate_positive::Bool=false)
    m, n = size(A)
    R = base_ring(A)
    n_vars = Oscar.n_variables(R)
    
    if n_vars == 0
        return A
    end
    
    # Create the result matrix
    result = deepcopy(A)
    
    # Get the variable list
    vars = gens(R)
    
    # Process one row at a time
    for i in 1:m
        # Find the minimum and maximum exponents of each variable across the row
        min_exps = fill(typemax(Int), n_vars)
        max_exps = fill(typemin(Int), n_vars)
        
        for j in 1:n
            elem = A[i, j]
            if !iszero(elem)
                # Get the exponent vectors for this entry
                exp_vecs = AbstractAlgebra.exponent_vectors(elem)
                
                for exps in exp_vecs
                    for k in 1:n_vars
                        min_exps[k] = min(min_exps[k], exps[k])
                        max_exps[k] = max(max_exps[k], exps[k])
                    end
                end
            end
        end
        
        # Construct the compensating monomial
        compensation = R(1)
        
        if eliminate_positive
            # Eliminate positive exponents: if the maximum exponent of x is 2, multiply by x^(-2)
            for k in 1:n_vars
                if max_exps[k] > 0
                    compensation *= vars[k]^(-max_exps[k])
                end
            end
        else
            # Eliminate negative exponents: if the minimum exponent of x is -2, multiply by x^2
            for k in 1:n_vars
                if min_exps[k] < 0
                    compensation *= vars[k]^(-min_exps[k])
                end
            end
        end
        
        # If compensation is needed, multiply the whole row by the compensating monomial
        if compensation != R(1)
            for j in 1:n
                result[i, j] = result[i, j] * compensation
            end
        end
    end
    
    return result
end

"""
    possible_remainders(new_var_rels, R_old, R_new, old_vars, new_vars)

Compute all possible remainders that can appear during variable replacement.

When a monomial in the old variables (x, y, ...) is expressed using new variables (u, v, ...),
the typical form is old_monomial = new_vars^powers * remainder

this function finds all possible remainder monomials that form a complete set of representatives.

# Arguments
- `new_var_rels`: Definitions of the new variables, for example [x*y^2, x^-1*y]
- `R_old`: Laurent polynomial ring in the old variables
- `R_new`: Laurent polynomial ring in the new variables
- `old_vars`: List of old variables
- `new_vars`: List of new variables

# Returns
- An array of possible remainder monomials

# Idea
By Smith normal form theory, the size of the quotient group Z^n / L equals the determinant of the lattice matrix.
We enumerate representatives to obtain a complete and simple remainder set.

# Example
```julia
R, (x, y) = laurent_polynomial_ring(GF(2), ["x", "y"])
Ruv, (u, v) = laurent_polynomial_ring(GF(2), ["u", "v"])
u_rel = x*y^2
v_rel = x^-1*y
rems = possible_remainders([u_rel, v_rel], R, Ruv, [x, y], [u, v])
# Returns [R(1), x, x^-1] or another equivalent representative system
```
"""
function possible_remainders(new_var_rels, R_old, R_new, old_vars, new_vars)
    n_old = length(old_vars)
    n_new = length(new_vars)
    
    # Helper: get the exponent vector of a monomial
    function get_exp_vector(mono, vars)
        if iszero(mono) || mono == R_old(0)
            return zeros(Int, length(vars))
        end
        
        exps = zeros(Int, length(vars))
        
        try
            exp_vecs = AbstractAlgebra.exponent_vectors(mono)
            if !isempty(exp_vecs)
                exp_vec = first(exp_vecs)
                exps = [exp_vec[i] for i in 1:length(vars)]
            end
        catch
            # Fallback: parse the string representation
        end
        
        return exps
    end
    
    # Build the coefficient matrix: M[i, j] is the exponent of old_vars[j] in new_var_rels[i]
    M = zeros(Int, n_new, n_old)
    for i in 1:n_new
        exp_vec = get_exp_vector(new_var_rels[i], old_vars)
        M[i, :] = exp_vec
    end
    
    # Compute the determinant
    if n_new == n_old
        det_val = abs(round(Int, det(M)))
    else
        # If the dimensions do not match, use a heuristic default
        det_val = 10
    end
    
    if det_val == 0
        det_val = 10  # default value when the matrix is rank-deficient
    end
    
    # Strategy: enumerate lattice points in Z^n_old and find representatives modulo M*Z^n_new
    # Simplified version: only consider different exponents of the first variable
    remainders = [R_old(1)]
    
    # Optimization for the two-variable case
    if n_old == 2 && n_new == 2
        # Enumerate exponents of the first variable while keeping the second variable exponent at 0
        for exp_x in 1:(det_val-1)
            mono = old_vars[1]^exp_x
            push!(remainders, mono)
        end
    elseif n_old == 2
        # General case: enumerate monomials in a small range
        tried_exps = Set{Tuple{Int,Int}}()
        push!(tried_exps, (0, 0))
        
        # Enumerate simple monomials
        for exp1 in -det_val:det_val
            for exp2 in -det_val:det_val
                if (exp1, exp2) in tried_exps
                    continue
                end
                
                # Check whether this exponent vector lies in a new equivalence class relative to the existing remainders
                is_new_class = true
                current_vec = [exp1, exp2]
                
                for rem in remainders
                    rem_vec = get_exp_vector(rem, old_vars)
                    diff = current_vec - rem_vec
                    
                    # Check whether diff can be expressed by the columns of M
                    # Simplified test: if diff is an integer linear combination of columns of M, it lies in the same equivalence class
                    try
                        if n_new == 2
                            # Solve M^T * c = diff
                            MT = transpose(M)
                            # Solve over the integers
                            if det(MT) != 0
                                c = MT \ diff
                                if all(isinteger.(c))
                                    is_new_class = false
                                    break
                                end
                            end
                        end
                    catch
                    end
                end
                
                if is_new_class && length(remainders) < det_val
                    mono = old_vars[1]^exp1 * old_vars[2]^exp2
                    push!(remainders, mono)
                    push!(tried_exps, (exp1, exp2))
                end
                
                if length(remainders) >= det_val
                    break
                end
            end
            if length(remainders) >= det_val
                break
            end
        end
    end
    
    # If there are still not enough remainders, add simple monomials
    while length(remainders) < det_val
        idx = length(remainders)
        if idx < length(old_vars) + 1
            push!(remainders, old_vars[idx])
        else
            break
        end
    end
    
    return remainders
end

"""
    find_minimum_triangle(m)

Find the minimum-area triangle in matrix m that includes position (1, 1).

# Arguments
- `m`: Integer matrix

# Returns
- `min_triangle`: The three vertices (p1, p2, p3) of the minimum triangle, with p1 = (1, 1)
- `min_area`: The area of the minimum triangle

# Example
```julia
m = [1 0 0 1;
     0 0 1 0;
     0 1 0 0;
     1 0 0 1]
vertices, area = find_minimum_triangle(m)
```
"""
function find_minimum_triangle(m)
    n_rows, n_cols = size(m)
    
    # Collect all positions with value 1
    one_positions = []
    for i in 1:n_rows
        for j in 1:n_cols
            if m[i, j] == 1
                push!(one_positions, (i, j))
            end
        end
    end
    
    # Fix the first vertex at [1, 1]
    p1 = (1, 1)
    
    # Initialize with the fallback triangle [1, 1], [1, n_cols], [n_rows, 1]
    min_area = abs((1 - 1) * (n_rows - 1) - (n_rows - 1) * (n_cols - 1)) / 2
    min_triangle = (p1, (1, n_cols), (n_rows, 1))
    
    # Compute the cosine of the initial triangle angle for tie-breaking
    v1_init = ((1, n_cols)[1] - p1[1], (1, n_cols)[2] - p1[2])
    v2_init = ((n_rows, 1)[1] - p1[1], (n_rows, 1)[2] - p1[2])
    dot_init = v1_init[1] * v2_init[1] + v1_init[2] * v2_init[2]
    len1_init = sqrt(v1_init[1]^2 + v1_init[2]^2)
    len2_init = sqrt(v2_init[1]^2 + v2_init[2]^2)
    min_cos_angle = dot_init / (len1_init * len2_init)
    
    # Iterate over all possible vertex pairs
    for (i, p2) in enumerate(one_positions)
        if p2 == p1
            continue
        end
        
        for p3 in one_positions[(i+1):end]
            if p3 == p1 || p3 == p2
                continue
            end
            
            # Check whether the three points are collinear and skip zero-area cases
            # Compute the triangle area using the cross-product formula
            # Given p1=(x1, y1), p2=(x2, y2), p3=(x3, y3)
            # area = |x1(y2-y3) + x2(y3-y1) + x3(y1-y2)| / 2
            
            x1, y1 = p1[1], p1[2]
            x2, y2 = p2[1], p2[2]
            x3, y3 = p3[1], p3[2]
            
            area_2x = abs(x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2))
            area = area_2x / 2
            
            # If the area is 0, the points are collinear, so skip them
            if area == 0
                continue
            end
            
            # Compute the angle between the two edges
            # v1 = p2 - p1, v2 = p3 - p1
            v1 = (p2[1] - p1[1], p2[2] - p1[2])
            v2 = (p3[1] - p1[1], p3[2] - p1[2])
            
            # Compute the cosine of the angle: cos(theta) = (v1 * v2) / (|v1| * |v2|)
            dot_product = v1[1] * v2[1] + v1[2] * v2[2]
            len1 = sqrt(v1[1]^2 + v1[2]^2)
            len2 = sqrt(v2[1]^2 + v2[2]^2)
            cos_angle = dot_product / (len1 * len2)
            
            # Update the minimum-area triangle
            # Prefer smaller area; if the area ties, prefer the larger angle (smaller cosine)
            should_update = false
            
            if area < min_area
                should_update = true
            elseif abs(area - min_area) < 1e-10  # Same area, up to floating-point tolerance
                # A larger angle corresponds to a smaller cosine
                if cos_angle < min_cos_angle
                    should_update = true
                end
            end
            
            if should_update
                min_area = area
                min_cos_angle = cos_angle
                min_triangle = (p1, p2, p3)
            end
        end
    end
    
    return min_triangle, min_area
end

"""
    find_parallelogram_lattice_points(p1, p2, p3)

Find all integer lattice points satisfying the selection rules inside the parallelogram defined by three vertices.

# Arguments
- `p1`: First vertex, typically (1, 1)
- `p2`: Second vertex
- `p3`: Third vertex

# Returns
- List of integer lattice points satisfying the conditions

# Parallelogram definition
The parallelogram defined by the three vertices p1, p2, p3:
- Edge vectors: v1 = p2 - p1, v2 = p3 - p1
- Fourth vertex: p4 = p1 + v1 + v2
- Points in the parallelogram: P = p1 + s*v1 + t*v2

# Integer-point selection rules
Select integer points satisfying the following conditions:
- 0 <= s < 1 and 0 <= t < 1 (parameter range)
- Include vertex p1
- Include integer points on edge p1-p2, excluding endpoint p2
- Include integer points on edge p1-p3, excluding endpoint p3
- Exclude integer points on edges p2-p4 and p3-p4
- Exclude vertices p2, p3, and p4
- Include all interior integer points of the parallelogram

# Algorithm
1. Compute the edge vectors v1 and v2
2. For each integer point (x, y) inside the bounding box
3. Solve for the parameters (s, t) using a linear system
4. Check whether 0 <= s < 1 and 0 <= t < 1

# Example
```julia
p1 = (1, 1)
p2 = (1, 3)
p3 = (3, 1)
points = find_parallelogram_lattice_points(p1, p2, p3)
# Returns: [(1,1), (1,2), (2,1), (2,2)]
```
"""
function find_parallelogram_lattice_points(p1, p2, p3)
    # Compute the edge vectors
    v1 = (p2[1] - p1[1], p2[2] - p1[2])
    v2 = (p3[1] - p1[1], p3[2] - p1[2])
    
    # Compute the fourth vertex
    p4 = (p1[1] + v1[1] + v2[1], p1[2] + v1[2] + v2[2])
    
    # Determine the search range (bounding box)
    x_min = min(p1[1], p2[1], p3[1], p4[1])
    x_max = max(p1[1], p2[1], p3[1], p4[1])
    y_min = min(p1[2], p2[2], p3[2], p4[2])
    y_max = max(p1[2], p2[2], p3[2], p4[2])
    
    # Compute det(v1, v2), which is used to test whether a point lies in the parallelogram
    det = v1[1] * v2[2] - v1[2] * v2[1]
    
    if det == 0
        # Degenerate case: the two vectors are collinear
        return [(p1[1], p1[2])]
    end
    
    lattice_points = Tuple{Int64, Int64}[]
    
    # Iterate over all integer points in the bounding box
    for x in x_min:x_max
        for y in y_min:y_max
            # Compute the position relative to p1
            dx = x - p1[1]
            dy = y - p1[2]
            
            # Solve the linear system: dx = s*v1[1] + t*v2[1], dy = s*v1[2] + t*v2[2]
            # Use Cramer's rule
            s = (dx * v2[2] - dy * v2[1]) / det
            t = (v1[1] * dy - v1[2] * dx) / det
            
            # Check whether the point lies in the valid parameter range
            # Condition: 0 <= s < 1 and 0 <= t < 1
            # Use a tolerance to handle floating-point precision
            tolerance = 1e-10
            
            if s >= -tolerance && s < 1 - tolerance && t >= -tolerance && t < 1 - tolerance
                # Verify that this is indeed an integer lattice point, since floating-point arithmetic may introduce error
                if abs(round(x) - x) < tolerance && abs(round(y) - y) < tolerance
                    push!(lattice_points, (Int(round(x)), Int(round(y))))
                end
            end
        end
    end
    
    return lattice_points
end

"""
    replace_variable_inv(a, new_var_rels, R_old, R_new, old_vars, new_vars)

Convert a monomial from the new-variable ring back to the old-variable ring (the inverse of replace_variable).

# Arguments
- `a`: monomial in the new-variable ring `R_new`
- `new_var_rels`: List of new-variable relations expressed in the old variables, for example [x*y^2, x^-1*y]
- `R_old`: Laurent polynomial ring in the old variables
- `R_new`: Laurent polynomial ring in the new variables
- `old_vars`: list of old variables, for example `[x, y]`
- `new_vars`: list of new variables, for example `[u, v]`

# Returns
- monomial in the old-variable ring

# Idea
If u = u_rel and v = v_rel, then u^c_u * v^c_v = u_rel^c_u * v_rel^c_v

# Example
```julia
F = GF(2)
R, (x, y) = laurent_polynomial_ring(F, ["x", "y"])
Ruv, (u, v) = laurent_polynomial_ring(F, ["u", "v"])

u_rel = x*y^2
v_rel = x^-1*y

# Convert u*v back to the old variables
result = replace_variable_inv(u*v, [u_rel, v_rel], R, Ruv, [x, y], [u, v])
# Returns: x*y^2 * x^-1*y = y^3
```
"""
function replace_variable_inv(a, new_var_rels, R_old, R_new, old_vars, new_vars)
    n_new = length(new_vars)
    
    # Helper: extract the exponents of a monomial
    function get_exponents_oscar(poly, vars)
        if iszero(poly)
            return nothing
        end
        
        exp_vecs = AbstractAlgebra.exponent_vectors(poly)
        exp_vec_array = collect(exp_vecs)
        
        if length(exp_vec_array) != 1
            error("The input polynomial must be a monomial, but it contains $(length(exp_vec_array)) terms")
        end
        
        exps = exp_vec_array[1]
        
        # Build the exponent dictionary
        exp_dict = Dict()
        for (i, var) in enumerate(vars)
            exp_dict[var] = exps[i]
        end
        
        return exp_dict
    end
    
    # Extract the exponents of the monomial in the new variables
    a_exp = get_exponents_oscar(a, new_vars)
    
    if isnothing(a_exp)
        return R_old(0)
    end
    
    # Construct the result by replacing each new variable with its old-variable relation
    result = R_old(1)
    for i in 1:n_new
        exp = a_exp[new_vars[i]]
        if exp != 0
            result *= new_var_rels[i]^exp
        end
    end
    
    return result
end

"""
    replace_variable_inv_mat(A, new_var_rels, R_old, R_new, old_vars, new_vars)

Convert a matrix over the new-variable ring back to the old-variable ring by applying replace_variable_inv to each entry.

# Arguments
- `A`: matrix over the new-variable ring `R_new`
- `new_var_rels`: list of new-variable relations expressed in the old variables
- `R_old`: Laurent polynomial ring in the old variables
- `R_new`: Laurent polynomial ring in the new variables
- `old_vars`: List of old variables
- `new_vars`: List of new variables

# Returns
- matrix over the old-variable ring

# Notes
- Matrix entries may be polynomials, i.e. sums of multiple monomials
- Apply replace_variable_inv to each monomial separately and sum the results

# Example
```julia
F = GF(2)
R, (x, y) = laurent_polynomial_ring(F, ["x", "y"])
Ruv, (u, v) = laurent_polynomial_ring(F, ["u", "v"])

u_rel = x*y^2
v_rel = x^-1*y

# Matrix over the new-variable ring
A = zero_matrix(Ruv, 2, 2)
A[1, 1] = u + v
A[1, 2] = u*v
A[2, 1] = u^2
A[2, 2] = 1 + u

# Convert back to the old-variable ring
B = replace_variable_inv_mat(A, [u_rel, v_rel], R, Ruv, [x, y], [u, v])
# B[1,1] = x*y^2 + x^-1*y
# B[1,2] = y^3
# B[2,1] = x^2*y^4
# B[2,2] = 1 + x*y^2
```
"""
function replace_variable_inv_mat(A, new_var_rels, R_old, R_new, old_vars, new_vars)
    rows, cols = size(A)
    
    # Create the result matrix
    result = zero_matrix(R_old, rows, cols)
    
    # Apply replace_variable_inv to each matrix entry
    for i in 1:rows
        for j in 1:cols
            entry = A[i, j]
            
            # Skip zero entries
            if iszero(entry)
                continue
            end
            
            # Handle polynomials by converting each monomial separately
            # Get all terms of the polynomial
            result_poly = R_old(0)
            
            for (coef, exp_vec) in zip(Oscar.coefficients(entry), exponent_vectors(entry))
                # Construct the monomial
                mono = mapreduce(k -> new_vars[k]^(exp_vec[k]), *, 1:length(new_vars); init=R_new(1))
                
                # Convert the monomial
                mono_in_old = replace_variable_inv(mono, new_var_rels, R_old, R_new, old_vars, new_vars)
                
                # Multiply by the coefficient and add it
                result_poly += R_old(coef) * mono_in_old
            end
            
            result[i, j] = result_poly
        end
    end
    
    return result
end

function find_minimum_triangle_ideal(variable_list, I;max_int = 1000)
    @assert length(variable_list) == 2
    l1 = nothing
    for j in 1:max_int
        if (variable_list[1]^j+1) in I
            l1 = j
            break
        end
    end
    isnothing(l1) && return nothing

    for j in 1:max_int
        for i in 0:l1
            if (variable_list[2]^j*variable_list[1]^i+1) in I
                return (i,j),(l1,0), l1*j
            end
        end
    end
    return nothing
end
