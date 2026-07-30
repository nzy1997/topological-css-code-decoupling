"""
Perform row elimination on an ExcitationMatrix
Use column start_column to eliminate entries from start_column + 1 through end_column
"""
function row_elemination!(em::AbstractExcitationMatrix, row_num,column_num, start_column, end_column)
    @assert em.m[row_num, column_num] == 1 "Pivot element must be 1"
    
    for i in start_column:end_column
        if em.m[row_num, i] != 0 && i != column_num
            a = em.m[row_num, i]
            mul_add_col!(em, i, column_num, a)
        end
    end
    return em
end

"""
Perform column elimination on an ExcitationMatrix
Use row start_row to eliminate entries from start_row + 1 through end_row
"""
function column_elemination!(em::AbstractExcitationMatrix,row_num,column_num, start_row, end_row)
    @assert em.m[row_num, column_num] == 1 "Pivot element must be 1, got: $(em.m[start_row, row_num])"
    for i in start_row:end_row
        if !iszero(em.m[i, column_num]) && (i != row_num)
            a = em.m[i, column_num]
            mul_add_row!(em, i, row_num, a)
        end
    end
    return em
end

function gaussian_elimination!(em::ExcitationMatrixwithF2Matrix, start_column::Int, min_column::Int, max_column::Int, start_row::Int, min_row::Int,max_row::Int)
    for i in 0:max_row-start_row
        j = findfirst(!iszero, em.m[start_row+i:max_row,start_column+i])
        if isnothing(j)
            j = findfirst(!iszero, em.m[start_row+i,start_column+i:max_column])
            if isnothing(j)
                pos = findfirst(!iszero, em.m[start_row+i:max_row,start_column+i:max_column])
                if isnothing(pos)
                    continue
                end
                row_switch!(em, start_row+i, start_row+i+pos.I[1]-1)
                column_switch!(em, start_column+i, start_column+i+pos.I[2]-1)
            else
                column_switch!(em, start_column+i, start_column+j+i-1)
            end
        elseif j > 1
            row_switch!(em, start_row+i, start_row+j+i-1)
        end
        column_elemination!(em, start_row+i, start_column+i, min_row, max_row)
        row_elemination!(em, start_row+i,start_column+i, min_column, max_column)
    end
end

_dagger_laurent_matrix(A) = transpose(laurent_conjugate.(A))

function _solve_selected_columns(A, B, kept_cols::AbstractVector{<:Integer})
    R = base_ring(A)
    result = zero_matrix(R, ncols(A), ncols(B))
    isempty(kept_cols) && return result

    rhs = zero_matrix(R, nrows(B), length(kept_cols))
    for (kept_idx, original_col) in enumerate(Int.(kept_cols))
        rhs[:, kept_idx:kept_idx] = B[:, original_col:original_col]
    end

    solved = solve_laurent_linear(A, rhs)
    for (kept_idx, original_col) in enumerate(Int.(kept_cols))
        result[:, original_col:original_col] = solved[:, kept_idx:kept_idx]
    end

    return result
end

function _build_initial_pdagger(Hxdagger, psi_1_inverse_local, Hxtdagger, product_state_num::Int, toric_num::Int)
    R = base_ring(Hxdagger)
    stab_num = ncols(Hxdagger)
    Pdagger = zero_matrix(R, stab_num, stab_num)

    if product_state_num > 0
        expected_product_cols = psi_1_inverse_local[:, product_state_num+1:2*product_state_num]
        Hxdagger[:, 1:product_state_num] == expected_product_cols || throw(
            ErrorException("_build_initial_pdagger product-state columns do not match known Hxdagger columns"),
        )

        Pdagger[1:product_state_num, 1:product_state_num] =
            identity_matrix(R, product_state_num)
    end

    if toric_num > 0
        toric_cols = product_state_num + 1:stab_num
        rhs_toric = psi_1_inverse_local * Hxtdagger[:, toric_cols]
        Pdagger[:, toric_cols] = solve_laurent_linear(Hxdagger, rhs_toric)
        Hxdagger * Pdagger[:, toric_cols] == rhs_toric || throw(
            ErrorException("_build_initial_pdagger failed toric columns of Hxdagger * Pdagger == psi_1_inverse_local * Hxtdagger"),
        )
    end

    return Pdagger
end

function _laurent_path_sum(v, k::Integer)
    R = parent(v)
    result = zero(R)

    if k > 0
        for power in 1:k
            result += v^power
        end
    elseif k < 0
        for power in (k + 1):0
            result += v^power
        end
    end

    return result
end

function _decompose_augmentation_entry(d)
    R = parent(d)
    Oscar.n_variables(R) == 2 || throw(
        ErrorException("_decompose_augmentation_entry expects a two-variable Laurent ring"),
    )

    one_values = fill(R(1), Oscar.n_variables(R))
    iszero(evaluate(d, one_values)) || throw(
        ErrorException("correction entry is not in the augmentation ideal"),
    )

    x, y = gens(R)
    a = zero(R)
    b = zero(R)

    for (coef, exp_vec) in zip(Oscar.coefficients(d), AbstractAlgebra.exponent_vectors(d))
        iszero(coef) && continue
        isone(coef) || throw(
            ErrorException("expected GF(2) Laurent coefficients in correction entry"),
        )

        mx = Int(exp_vec[1])
        my = Int(exp_vec[2])
        a += _laurent_path_sum(y, my)
        b += y^my * _laurent_path_sum(x, mx)
    end

    return a, b
end

function _assert_standard_hxt_shape(Hxt, product_state_num::Int, toric_num::Int)
    R = base_ring(Hxt)
    Oscar.n_variables(R) == 2 || throw(
        ErrorException("_assert_standard_hxt_shape expects a two-variable Laurent ring"),
    )

    x, y = gens(R)
    stab_num = product_state_num + toric_num
    qubit_num = 2 * product_state_num + 2 * toric_num
    size(Hxt) == (stab_num, qubit_num) || throw(
        DimensionMismatch("standard Hxt must have size ($stab_num, $qubit_num), got $(size(Hxt))"),
    )

    for row in 1:stab_num
        for col in 1:qubit_num
            expected = zero(R)
            if row <= product_state_num
                if col == row + product_state_num
                    expected = R(1)
                end
            else
                toric_row = row - product_state_num
                if col == 2 * product_state_num + 2 * toric_row - 1
                    expected = y^-1 + 1
                elseif col == 2 * product_state_num + 2 * toric_row
                    expected = x^-1 + 1
                end
            end

            Hxt[row, col] == expected || throw(
                ErrorException("standard Hxt shape mismatch at ($row, $col)"),
            )
        end
    end

    return nothing
end

function _solve_standard_hxt_correction(Hxt, P, Pt, product_state_num::Int, toric_num::Int)
    R = base_ring(Hxt)
    target = P + Pt
    stab_num = product_state_num + toric_num
    qubit_num = ncols(Hxt)

    _assert_standard_hxt_shape(Hxt, product_state_num, toric_num)

    size(target) == (stab_num, stab_num) || throw(
        DimensionMismatch("P + Pt must have size ($stab_num, $stab_num), got $(size(target))"),
    )

    for row in 1:stab_num
        for col in 1:stab_num
            if (row <= product_state_num || col <= product_state_num) && !iszero(target[row, col])
                throw(
                    ErrorException(
                        "standard Hxt correction expected P + Pt to vanish outside lower-right toric block; found nonzero at ($row, $col)",
                    ),
                )
            end
        end
    end

    Adagger = zero_matrix(R, qubit_num, stab_num)
    for toric_row in 1:toric_num
        for toric_col in 1:toric_num
            d = target[product_state_num + toric_row, product_state_num + toric_col]
            a, b = _decompose_augmentation_entry(d)
            a_row = product_state_num + toric_row
            a_col = 2 * product_state_num + 2 * toric_row - 1
            b_col = 2 * product_state_num + 2 * toric_row
            target_col = product_state_num + toric_col

            Hxt[a_row, a_col] * a + Hxt[a_row, b_col] * b == d || throw(
                ErrorException("_solve_standard_hxt_correction failed local toric entry reconstruction"),
            )
            Adagger[a_col, target_col] = a
            Adagger[b_col, target_col] = b
        end
    end

    return Adagger
end

function _build_psi1_inverse_toric_correction(Hxdagger, Adagger, product_state_num::Int, toric_num::Int)
    stab_num = product_state_num + toric_num
    qubit_num = 2 * product_state_num + 2 * toric_num

    size(Hxdagger) == (qubit_num, stab_num) || throw(
        DimensionMismatch("Hxdagger must have size ($qubit_num, $stab_num), got $(size(Hxdagger))"),
    )
    size(Adagger) == (qubit_num, stab_num) || throw(
        DimensionMismatch("Adagger must have size ($qubit_num, $stab_num), got $(size(Adagger))"),
    )

    toric_stab_cols = product_state_num + 1:stab_num
    toric_qubit_cols = 2 * product_state_num + 1:qubit_num

    for row in 1:qubit_num
        row_is_toric = row >= first(toric_qubit_cols) && row <= last(toric_qubit_cols)
        for col in 1:stab_num
            col_is_toric = col >= first(toric_stab_cols) && col <= last(toric_stab_cols)
            if !(row_is_toric && col_is_toric) && !iszero(Adagger[row, col])
                throw(
                    ErrorException(
                        "psi_1_inverse toric correction expected Adagger to vanish outside the toric block; found nonzero at ($row, $col)",
                    ),
                )
            end
        end
    end

    return (
        target_cols=toric_qubit_cols,
        left=Hxdagger[:, toric_stab_cols],
        right=_dagger_laurent_matrix(Adagger[toric_qubit_cols, toric_stab_cols]),
    )
end

function _apply_psi1_inverse_toric_correction!(psi_1_inverse, correction)
    isempty(correction.target_cols) && return psi_1_inverse
    psi_1_inverse[:, correction.target_cols] += correction.left * correction.right
    return psi_1_inverse
end

function _left_multiply_column_blocks(left, right, column_blocks)
    ncols(left) == nrows(right) || throw(
        DimensionMismatch(
            "left column count $(ncols(left)) must match right row count $(nrows(right))",
        ),
    )

    R = base_ring(left)
    result = zero_matrix(R, nrows(left), ncols(right))
    isempty(column_blocks) && return result

    for cols in column_blocks
        isempty(cols) && continue
        first_col = first(cols)
        last_col = last(cols)
        1 <= first_col <= last_col <= ncols(right) || throw(
            BoundsError(right, (:, cols)),
        )

        for row in 1:nrows(left)
            for inner in 1:ncols(left)
                left_value = left[row, inner]
                iszero(left_value) && continue
                for col in cols
                    right_value = right[inner, col]
                    iszero(right_value) && continue
                    result[row, col] += left_value * right_value
                end
            end
        end
    end

    return result
end

function _compose_psi1_inverse_with_toric_correction(column_transformation_block, psi_1_inverse_local, correction)
    psi_1_inverse = _left_multiply_column_blocks(
        column_transformation_block,
        psi_1_inverse_local,
        (1:ncols(psi_1_inverse_local),),
    )
    isnothing(correction) && return psi_1_inverse

    fused_correction = (
        target_cols=correction.target_cols,
        left=_left_multiply_column_blocks(
            column_transformation_block,
            correction.left,
            (1:ncols(correction.left),),
        ),
        right=correction.right,
    )
    return _apply_psi1_inverse_toric_correction!(psi_1_inverse, fused_correction)
end

function _block_diagonal_toric_matrix(Hz, Hx)
    R = Hz.base_ring
    stab_num, qubit_num = size(Hz)
    return [
        Hz zero_matrix(R, stab_num, qubit_num)
        zero_matrix(R, stab_num, qubit_num) Hx
    ]
end

function _toric_form_blocks(input_matrix, standard_matrix, row_transformation)
    stab_num = size(input_matrix, 1) ÷ 2
    qubit_num = size(input_matrix, 2) ÷ 2
    return (
        input_blocks=(
            Hz=input_matrix[1:stab_num, 1:qubit_num],
            Hx=input_matrix[stab_num+1:2*stab_num, qubit_num+1:2*qubit_num],
        ),
        standard_blocks=(
            Hz=standard_matrix[1:stab_num, 1:qubit_num],
            Hx=standard_matrix[stab_num+1:2*stab_num, qubit_num+1:2*qubit_num],
        ),
        row_blocks=(
            Hz=row_transformation[1:stab_num, 1:stab_num],
            Hx=row_transformation[stab_num+1:2*stab_num, stab_num+1:2*stab_num],
        ),
    )
end

function capture_toric_form_debug_matrices(A::AbstractAlgebra.Generic.MatSpaceElem; show_progress::Bool=true, compute_inverse::Bool=false)
    return capture_toric_form_debug_matrices(
        ExcitationMatrix(copy(A));
        show_progress=show_progress,
        compute_inverse=compute_inverse,
    )
end

function capture_toric_form_debug_matrices(E::ExcitationMatrix; show_progress::Bool=true, compute_inverse::Bool=false)
    input_matrix = copy(E.m)
    qubit_num = size(E.m, 2) ÷ 2
    stab_num = size(E.m, 1) ÷ 2

    E2 = ExcitationMatrixwithF2Matrix(copy(E))
    gaussian_elimination!(E2, 1, 1, qubit_num, 1, 1, stab_num)

    product_state_num = findfirst(j -> iszero(E2.m[j, j]), 1:stab_num) - 1
    toric_num = stab_num - product_state_num

    gaussian_elimination!(
        E2,
        qubit_num + product_state_num + 1,
        qubit_num + 1,
        2 * qubit_num,
        stab_num + 1,
        stab_num + 1,
        2 * stab_num,
    )

    em_matrix = copy(E2.em.m)
    res = to_toric_form(copy(input_matrix); show_progress=show_progress, compute_inverse=compute_inverse)

    return (;
        input_matrix=input_matrix,
        em_matrix=em_matrix,
        row_transformation=res.row_transformation,
        row_blocks=res.row_blocks,
        psi_1_inverse=res.psi_1_inverse,
        psi_1=res.psi_1,
        column_transformation=res.column_transformation,
        standard_matrix=res.standard_matrix,
        standard_blocks=res.standard_blocks,
        product_state_num=product_state_num,
        toric_num=toric_num,
    )
end


to_toric_form(A::AbstractAlgebra.Generic.MatSpaceElem; show_progress::Bool=true, compute_inverse::Bool=false) =
    to_toric_form(ExcitationMatrix(copy(A)); show_progress=show_progress, compute_inverse=compute_inverse)
function to_toric_form(E::ExcitationMatrix; show_progress::Bool=true, compute_inverse::Bool=false)
    Rl = E.m[1,1].parent
    (x,y) = gens(Rl)
    qubit_num = size(E.m, 2) ÷ 2
    stab_num = size(E.m, 1) ÷ 2
    input_matrix = copy(E.m)
    # t_start = time()
    if show_progress
        println("--------------------------------")
        println("Guassian elemination on F2 matrix")
        println("--------------------------------")
    end
    E2 = ExcitationMatrixwithF2Matrix(E)
    gaussian_elimination!(E2, 1, 1, qubit_num, 1, 1, stab_num)
    # display(E2.m)
    # return E2
    product_state_num = findfirst(j -> iszero(E2.m[j,j]),1:stab_num)-1
    toric_num = stab_num - product_state_num
    gaussian_elimination!(E2, qubit_num + product_state_num +1, qubit_num +1, 2*qubit_num, stab_num+1, stab_num+1, 2*stab_num)

    # display(E2.m)
    if show_progress
        @show product_state_num
        @show toric_num
    end
    Em = E2.em
    # display(Em.m)
    # return Em

    # @show time() - t_start

    Hzt = zero_matrix(Rl, stab_num, qubit_num)
    for i in 1:product_state_num
        Hzt[i,i] = Rl(1)
    end

    for i in 1:toric_num
        Hzt[product_state_num+i,2*product_state_num + 2*i - 1] = x+1
        Hzt[product_state_num+i,2*product_state_num + 2*i] = y+1
    end

    Hz = Em.m[1:stab_num,1:qubit_num]
    # return Hzt,Hz


    Hxt = zero_matrix(Rl, stab_num, qubit_num)
    for i in 1:product_state_num
        Hxt[i,i+product_state_num] = Rl(1)
    end

    for i in 1:toric_num
        Hxt[product_state_num+i,2*product_state_num + 2*i - 1] = y^-1+1
        Hxt[product_state_num+i,2*product_state_num + 2*i] = x^-1+1
    end
    Hx = Em.m[1+stab_num:2*stab_num,qubit_num+1:2*qubit_num]

    # @show time() - t_start
    # return Hx,Hz,Hxt,Hzt
    if show_progress
        println("--------------------------------")
        println("Solving the equation Hz*psi_1_inverse = Hzt, matrix size: Hz: $(size(Hz)), Hzt: $(size(Hzt))")
        println("--------------------------------")
    end
    kept_psi_1_inverse_cols = Int[]
    append!(kept_psi_1_inverse_cols, 1:product_state_num)
    append!(kept_psi_1_inverse_cols, 2 * product_state_num + 1:qubit_num)
    psi_1_inverse_local = _solve_selected_columns(Hz, Hzt, kept_psi_1_inverse_cols)
    # @show time() - t_start
    # @show Hzt == Hz*psi_1_inverse_local

    psi_1_inverse_local[:, product_state_num+1:2*product_state_num] = _dagger_laurent_matrix(Hx)[:, 1:product_state_num]
    # return Hx,Hz,Hxt,Hzt,psi_1_inverse_local
    # @show time() - t_start
    if show_progress
        println("--------------------------------")
        println("Solving the equation Hxdagger*Pdagger = psi_1_inverse_local*Hxtdagger, matrix size: Hxdagger: $(reverse(size(Hx))), psi_1_inverse_local: $(size(psi_1_inverse_local)), Hxtdagger: $(reverse(size(Hxt)))")
        println("--------------------------------")
    end
    Hxdagger = _dagger_laurent_matrix(Hx)
    Hxtdagger = _dagger_laurent_matrix(Hxt)
    Pdagger = _build_initial_pdagger(Hxdagger, psi_1_inverse_local, Hxtdagger, product_state_num, toric_num)
    P = _dagger_laurent_matrix(Pdagger)
    # @show time() - t_start
    # @show psi_1_inverse_local*transpose(laurent_conjugate.(Hxt)) == transpose(laurent_conjugate.(Hx))*transpose(laurent_conjugate.(P))
    # return P,psi_1_inverse_local

    @assert all(iszero,Pdagger[product_state_num+1:end,1:product_state_num]) "P is not in the correct form"
    pdet = det(P[product_state_num+1:end,product_state_num+1:end])
    psi_1_inverse_correction = nothing
    # @show time() - t_start
    if length(pdet) != 1
        if show_progress
            println("--------------------------------")
            println("det(B) is a polynomial, det(B) = $pdet")
            println("--------------------------------")
        end

        # Pt = map(P) do elem
        #     evaluate(elem,  fill(Rl(1),n_variables(P.base_ring)))
        # end

        Pt = copy(P)
        for i in product_state_num+1:stab_num
            for j in product_state_num+1:stab_num
                Pt[i,j] = evaluate(P[i,j],  fill(Rl(1),n_variables(P.base_ring)))
            end
        end
        # @show time() - t_start
        if show_progress
            println("Solving the equation Hxt*Adagger = P+Pt, matrix size: Hxt: $(size(Hxt)), P: $(size(P))")
            println("--------------------------------")
        end
        Adagger = _solve_standard_hxt_correction(Hxt, P, Pt, product_state_num, toric_num)
        # @show transpose(P+Pt) == Hxt*Atr
        psi_1_inverse_correction = _build_psi1_inverse_toric_correction(Hxdagger, Adagger, product_state_num, toric_num)
        P = Pt
        # @show time() - t_start
    elseif show_progress
        println("--------------------------------")
        println("det(B) is a monomial, det(B) = $pdet")
        println("--------------------------------")
    end

    # @show psi_1_inverse_local*transpose(laurent_conjugate.(Hxt)) == transpose(laurent_conjugate.(Hx))*transpose(laurent_conjugate.(P))
    # @show time() - t_start
    Em.row_transformation[stab_num+1:2*stab_num, stab_num+1:2*stab_num] =
        P * Em.row_transformation[stab_num+1:2*stab_num, stab_num+1:2*stab_num]

    psi_1_inverse = _compose_psi1_inverse_with_toric_correction(
        Em.column_transformation[1:qubit_num, 1:qubit_num],
        psi_1_inverse_local,
        psi_1_inverse_correction,
    )
    psi_1 = nothing
    column_transformation = nothing

    if compute_inverse
        if show_progress
            println("--------------------------------")
            println("Solving the equation psi_1_inverse*psi_1 = I, matrix size: psi_1_inverse: $(size(psi_1_inverse))")
            println("--------------------------------")
        end
        psi_1 = solve_laurent_linear(psi_1_inverse, identity_matrix(Rl, qubit_num))
        column_transformation = copy(Em.column_transformation)
        psi_1_local = psi_1 * Em.column_transformation[1:qubit_num, 1:qubit_num]
        column_transformation[1:qubit_num, 1:qubit_num] = psi_1_inverse
        column_transformation[qubit_num+1:2*qubit_num, qubit_num+1:2*qubit_num] =
            column_transformation[qubit_num+1:2*qubit_num, qubit_num+1:2*qubit_num] *
            _dagger_laurent_matrix(psi_1_local)
    end

    standard_matrix = _block_diagonal_toric_matrix(Hzt, Hxt)
    blocks = _toric_form_blocks(input_matrix, standard_matrix, Em.row_transformation)

    return (;
        input_matrix=input_matrix,
        input_blocks=blocks.input_blocks,
        standard_matrix=standard_matrix,
        standard_blocks=blocks.standard_blocks,
        row_transformation=Em.row_transformation,
        row_blocks=blocks.row_blocks,
        psi_1_inverse=psi_1_inverse,
        psi_1=psi_1,
        column_transformation=column_transformation,
        product_state_num=product_state_num,
        toric_num=toric_num,
    )
end

function check_result(res, A; require_inverse::Bool=false)
    R = A.base_ring
    stab_num = size(A, 1) ÷ 2
    qubit_num = size(A, 2) ÷ 2

    @assert res.input_matrix == A
    @assert res.input_blocks.Hz == A[1:stab_num, 1:qubit_num]
    @assert res.input_blocks.Hx == A[stab_num+1:2*stab_num, qubit_num+1:2*qubit_num]
    @assert res.standard_matrix == _block_diagonal_toric_matrix(res.standard_blocks.Hz, res.standard_blocks.Hx)
    @assert res.row_blocks.Hz * res.input_blocks.Hz * res.psi_1_inverse == res.standard_blocks.Hz
    @assert res.row_blocks.Hx * res.input_blocks.Hx ==
            res.standard_blocks.Hx * _dagger_laurent_matrix(res.psi_1_inverse)

    if require_inverse
        isnothing(res.psi_1) && throw(ArgumentError(
            "Full check requires psi_1; rerun with compute_inverse=true.",
        ))
        isnothing(res.column_transformation) && throw(ArgumentError("Full check requires column_transformation; rerun with compute_inverse=true."))
        identity = identity_matrix(R, size(res.psi_1_inverse, 1))
        @assert res.psi_1_inverse * res.psi_1 == identity
        @assert res.psi_1 * res.psi_1_inverse == identity
        @assert res.row_transformation * A * res.column_transformation == res.standard_matrix
    end

    return true
end
