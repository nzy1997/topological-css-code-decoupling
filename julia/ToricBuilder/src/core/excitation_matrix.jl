abstract type AbstractExcitationMatrix{MT} end
struct ExcitationMatrix{MT} <: AbstractExcitationMatrix{MT}
    m::MT
    row_transformation::MT
    column_transformation::MT
end

function ExcitationMatrix(A::MT) where MT
    m,n = size(A)
    row_transformation = identity_matrix(A.base_ring,m)
    column_transformation = identity_matrix(A.base_ring,n)
    return ExcitationMatrix{MT}(A,row_transformation,column_transformation)
end
Base.copy(em::ExcitationMatrix) = ExcitationMatrix(copy(em.m),copy(em.row_transformation),copy(em.column_transformation))

struct ExcitationMatrixwithF2Matrix{MT} <: AbstractExcitationMatrix{MT}
    em::ExcitationMatrix{MT}
    m::FqMatrix
end

function ExcitationMatrixwithF2Matrix(em::ExcitationMatrix{MT}) where MT
    F2 = Oscar.GF(2)
    one_F2 = F2(1)
    m = map(em.m) do elem
        evaluate(elem,  fill(one_F2,n_variables(em.m.base_ring))) # Evaluate the inner variables (x, y, ...) at 1
    end
    return ExcitationMatrixwithF2Matrix{MT}(em,m)
end
Base.copy(em::ExcitationMatrixwithF2Matrix) = ExcitationMatrixwithF2Matrix(copy(em.em),deepcopy(em.m))


function _column_switch!(A,col_idx_num1::Int,col_idx_num2::Int)
    qubit_num = size(A,2) ÷ 2
    if col_idx_num1>qubit_num && col_idx_num2>qubit_num
        col_idx_num1 = col_idx_num1 - qubit_num
        col_idx_num2 = col_idx_num2 - qubit_num
    elseif (col_idx_num1>qubit_num && col_idx_num2<=qubit_num) || (col_idx_num1<=qubit_num && col_idx_num2>qubit_num)
        error("col_idx_num1 and col_idx_num2 must be in the same group")
    end
    temp = A[:,col_idx_num1]
    A[:,col_idx_num1] = A[:,col_idx_num2]
    A[:,col_idx_num2] = temp
    temp = A[:,col_idx_num1+qubit_num]
    A[:,col_idx_num1+qubit_num] = A[:,col_idx_num2+qubit_num]
    A[:,col_idx_num2+qubit_num] = temp
    return A
end

function _add_row_unit!(A, row1::Int, row2::Int)
    for col in 1:size(A, 2)
        value = A[row2, col]
        if !iszero(value)
            A[row1, col] += value
        end
    end
    return A
end

function _add_column_unit!(A, col1::Int, col2::Int)
    for row in 1:size(A, 1)
        value = A[row, col2]
        if !iszero(value)
            A[row, col1] += value
        end
    end
    return A
end

function _mul_add_col!(A,col1,col2,a)
    qubit_num = size(A,2) ÷ 2
    if isone(a)
        _add_column_unit!(A, col1, col2)
        _add_column_unit!(A, mod1(col2 + qubit_num, 2 * qubit_num), mod1(col1 + qubit_num, 2 * qubit_num))
    else
        A[:,col1] += A[:,col2] .* a
        φ = hom(parent(a), parent(a), [gens(parent(a))[1]^-1, gens(parent(a))[2]^-1])
        A[:,mod1(col2+qubit_num,2*qubit_num)] += A[:,mod1(col1+qubit_num,2*qubit_num)] .* φ(a)
    end
    return A
end

function _row_switch!(A,row_idx_num1::Int,row_idx_num2::Int)
    temp = A[row_idx_num1,:]
    A[row_idx_num1,:] = A[row_idx_num2,:]
    A[row_idx_num2,:] = temp
    return A
end

function _mul_add_row!(A,row1,row2,a)
    if isone(a)
        _add_row_unit!(A, row1, row2)
    else
        A[row1,:] += A[row2,:] .* a
    end
    return A
end

function _row_scale!(A,row_num,a)
    A[row_num,:] .*= a
    return A
end
function _column_scale!(A,col_num,a)
    qubit_num = size(A,2) ÷ 2
    A[:,col_num] .*= a
    A[:,mod1(col_num+qubit_num,2*qubit_num)] .*= a
    return A
end
# ============================================
# Operations for the ExcitationMatrix type
# ============================================

"""
Perform a row swap on an ExcitationMatrix
Update the row-transformation matrix at the same time
"""
function row_switch!(em::ExcitationMatrix, row_idx_num1::Int, row_idx_num2::Int)
    # Swap rows in the main matrix
    _row_switch!(em.m, row_idx_num1, row_idx_num2)

    _row_switch!(em.row_transformation, row_idx_num1, row_idx_num2)
    return em
end

function column_scale!(em::ExcitationMatrix, col_num, a)
    _column_scale!(em.m, col_num, a)
    _column_scale!(em.column_transformation, col_num, a)
    return em
end

"""
Perform a row multiply-add on an ExcitationMatrix: row1 += row2 * a
Update the row-transformation matrix at the same time
"""
function mul_add_row!(em::ExcitationMatrix, row1::Int, row2::Int, a)
    _mul_add_row!(em.m, row1, row2, a)
    _mul_add_row!(em.row_transformation, row1, row2, a)
    return em
end

"""
Scale a row of an ExcitationMatrix: row_num *= a
Update the row-transformation matrix at the same time
"""
function row_scale!(em::ExcitationMatrix, row_num, a)
    # Apply the operation to the main matrix
    _row_scale!(em.m, row_num, a)
    
    # Apply the same operation to the row-transformation matrix
    _row_scale!(em.row_transformation, row_num, a)
    
    return em
end

"""
Perform a column swap on an ExcitationMatrix, including the paired qubit columns
Update the column-transformation matrix at the same time
"""
function column_switch!(em::ExcitationMatrix, col_idx_num1, col_idx_num2)
    _column_switch!(em.m, col_idx_num1, col_idx_num2)
    _column_switch!(em.column_transformation, col_idx_num1, col_idx_num2)
    return em
end

"""
Perform a column multiply-add on an ExcitationMatrix
col1 += col2 * a
col2+qubit_num += col1+qubit_num * a^-1
Update the column-transformation matrix at the same time
"""
function mul_add_col!(em::ExcitationMatrix, col1, col2, a)
    _mul_add_col!(em.m, col1, col2, a)
    _mul_add_col!(em.column_transformation, col1, col2, a)
    return em
end


function row_switch!(em::ExcitationMatrixwithF2Matrix, row_idx_num1::Int, row_idx_num2::Int)
    # Swap rows in the main matrix
    row_switch!(em.em, row_idx_num1, row_idx_num2)
    _row_switch!(em.m, row_idx_num1, row_idx_num2)
    return em
end

"""
Perform a row multiply-add on an ExcitationMatrix: row1 += row2 * a
Update the row-transformation matrix at the same time
"""
function mul_add_row!(em::ExcitationMatrixwithF2Matrix, row1::Int, row2::Int, a)
    mul_add_row!(em.em, row1, row2, a)
    F2 = Oscar.GF(2)
    one_F2 = F2(1)
    if a isa AbstractAlgebra.Generic.LaurentPolyWrap
        a = evaluate(a.poly, fill(one_F2,n_variables(em.em.m.base_ring)))
    end
    _mul_add_row!(em.m, row1, row2, a)
    return em
end

"""
Perform a column swap on an ExcitationMatrixwithF2Matrix, including the paired qubit columns
Update the column-transformation matrix at the same time
"""
function column_switch!(em::ExcitationMatrixwithF2Matrix, col_idx_num1, col_idx_num2)
    column_switch!(em.em, col_idx_num1, col_idx_num2)
    _column_switch!(em.m, col_idx_num1, col_idx_num2)
    return em
end

"""
Perform a column multiply-add on an ExcitationMatrix
col1 += col2 * a
col2+qubit_num += col1+qubit_num * a^-1
Update the column-transformation matrix at the same time
"""
function mul_add_col!(em::ExcitationMatrixwithF2Matrix, col1, col2, a)
    mul_add_col!(em.em, col1, col2, a)
    F2 = Oscar.GF(2)
    one_F2 = F2(1)
    if a isa AbstractAlgebra.Generic.LaurentPolyWrap
        a = evaluate(a.poly, fill(one_F2,n_variables(em.em.m.base_ring)))
    end
    _mul_add_col!(em.m, col1, col2, a)
    return em
end

function max_length_element(E::ExcitationMatrix)
    qubit_num = size(E.m, 2)
    stab_num = size(E.m, 1)
    max_length = 0
    for i in 1:qubit_num
        for j in 1:stab_num
            if !iszero(E.m[j,i])
                if max_length < length(E.m[j,i])
                    max_length = length(E.m[j,i])
                    @show E.m[j,i]
                end
            end
        end
    end
    return max_length
end

function laurent_conjugate(a_laurent)
    R = parent(a_laurent)
    φ = hom(R, R, [gens(R)[1]^-1, gens(R)[2]^-1])
    return φ(a_laurent)
end

function check_syplectic_form(Em::ExcitationMatrix,row1_idx,row2_idx)
    qubit_num = size(Em.m, 2) ÷ 2
    row1 = Em.m[row1_idx,:]
    row2 = vcat(Em.m[row2_idx,qubit_num+1:2*qubit_num],Em.m[row2_idx,1:qubit_num])
    row2 = laurent_conjugate.(row2)
    return sum(row1 .* row2)
end
check_syplectic_form(Em::ExcitationMatrix) = check_syplectic_form(Em.m)
function check_syplectic_form(A::AbstractAlgebra.Generic.MatSpaceElem)
    qubit_num = size(A, 2) ÷ 2
    stab_num = size(A, 1) ÷ 2
    for i in 1:2*stab_num
        row1 = vcat(A[i,qubit_num+1:2*qubit_num],A[i,1:qubit_num])
        row1 = laurent_conjugate.(row1)
        for j in i+1:2*stab_num
            row2 = A[j,:]
            ip = sum(row1 .* row2)
            if !iszero(ip)
                @show i,j,ip
                return false
            end
        end
    end
    return true
end

"""
    polynomial_degree(p)

Compute the degree of a Laurent polynomial, defined as the sum of the ranges 
of exponents for each variable.

For a Laurent polynomial p, for each variable, we find the minimum and maximum 
exponents appearing in all terms, then sum up (max_exp - min_exp) over all variables.

# Example
```julia
# For p = x^2 + x^(-1), variable x has exponents in [-1, 2], range = 3
# For p = x^2*y + x*y^(-1), x ∈ [1,2] (range=1), y ∈ [-1,1] (range=2), degree = 3
```
"""
function polynomial_degree(p)
    if iszero(p)
        return 0
    end
    
    # Get the number of variables
    R = parent(p)
    n_vars = Oscar.n_variables(R)
    
    if n_vars == 0
        return 0
    end
    
    # Initialize min and max exponents for each variable
    min_exps = fill(typemax(Int), n_vars)
    max_exps = fill(typemin(Int), n_vars)
    
    # Get exponent vectors for Laurent polynomial
    exp_vecs = AbstractAlgebra.exponent_vectors(p)
    
    # Iterate over all exponent vectors
    for exps in exp_vecs
        for i in 1:n_vars
            exp_i = exps[i]
            min_exps[i] = min(min_exps[i], exp_i)
            max_exps[i] = max(max_exps[i], exp_i)
        end
    end
    
    # Compute degree as sum of ranges
    degree = sum(max_exps .- min_exps)
    return degree
end

"""
    polynomial_degree_monomial(p)

Compute the degree of a Laurent polynomial as the maximum monomial degree,
where a monomial degree is the sum of absolute values of exponents.
"""
function polynomial_degree_monomial(p)
    if iszero(p)
        return 0
    end

    R = parent(p)
    n_vars = Oscar.n_variables(R)

    if n_vars == 0
        return 0
    end

    max_degree = 0
    exp_vecs = AbstractAlgebra.exponent_vectors(p)
    for exps in exp_vecs
        deg = sum(abs, exps)
        if deg > max_degree
            max_degree = deg
        end
    end
    return max_degree
end

"""
    polynomial_term_count(p)

Return the number of nonzero terms in a Laurent polynomial.
"""
function polynomial_term_count(p)
    if iszero(p)
        return 0
    end

    R = parent(p)
    n_vars = Oscar.n_variables(R)
    if n_vars == 0
        return 1
    end

    return length(AbstractAlgebra.exponent_vectors(p))
end

"""
    max_column_monomial_count(A)

Compute the maximum monomial count among columns in a polynomial matrix.
Each column count is the sum of `polynomial_term_count` for its entries.
"""
function max_column_monomial_count(A)
    if ncols(A) == 0
        return 0
    end

    max_count = 0
    for j in 1:ncols(A)
        col_count = 0
        for i in 1:nrows(A)
            col_count += polynomial_term_count(A[i, j])
        end
        if col_count > max_count
            max_count = col_count
        end
    end
    return max_count
end

"""
    maximum_deg_term(A)

Find the element with maximum degree in matrix A.

Returns the element in A with the largest degree, where degree is computed 
as the sum of exponent ranges across all variables.

# Example
```julia
R, (x,y) = laurent_polynomial_ring(GF(2), ["x","y"])
A = matrix(R, 2, 2, [x^2 + x^(-1), y, x, x*y^2 + x*y^(-1)])
max_elem = maximum_deg_term(A)  # Returns x^2 + x^(-1) (degree 3)
```
"""
function maximum_deg_term(A)
    max_degree = -1
    max_element = nothing
    
    for i in 1:nrows(A)
        for j in 1:ncols(A)
            elem = A[i, j]
            deg = polynomial_degree_monomial(elem)
            if deg > max_degree
                max_degree = deg
                max_element = elem
            end
        end
    end
    return max_element,max_degree
end
