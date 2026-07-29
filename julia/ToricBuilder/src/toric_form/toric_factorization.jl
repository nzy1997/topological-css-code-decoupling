function factor_toric_block(n::Int, x, y, R)
    n = n+1
    a = 1
    b = n
    for k in 1:isqrt(n)
        if n % k == 0
            a = k
            b = n ÷ k
        end
    end
    A = coarse_graining(excitation_matrix([1+x,1+y]), [x,y],[a,b])
    res2 = to_toric_form(A; show_progress=false, compute_inverse=true)

    column_transform_inv = solve_laurent_linear(res2.column_transformation,identity_matrix(R,size(res2.column_transformation)[1]))
    row_transform_inv = solve_laurent_linear(res2.row_transformation,identity_matrix(R,size(res2.row_transformation)[1]))
    return column_transform_inv, row_transform_inv, a, b
end

function split_product_state(res, size_vec::Tuple{Vararg{Integer}})
    return split_product_state(res, collect(size_vec))
end

function split_product_state(res, size_vec::AbstractVector{<:Integer})
    sizes = Int.(collect(size_vec))
    product_state_num = res.product_state_num
    toric_num = res.toric_num
    if length(sizes) != toric_num
        throw(ArgumentError("size_vec must contain one entry per toric block: got $(length(sizes)), expected $(toric_num)"))
    end
    if any(<(0), sizes)
        throw(ArgumentError("size_vec entries must be nonnegative: got $(sizes)"))
    end
    used_size = sum(sizes)
    if used_size > product_state_num
        throw(ArgumentError("size_vec requests $(used_size) product states, but only $(product_state_num) are available"))
    end

    R = res.standard_matrix.base_ring
    (u,v) = gens(R)
    used_product_state_num = 0
    stab_num,_ = size(res.standard_matrix)
    qubit_pos_list = Vector{Vector{Int}}(undef, toric_num)
    xstab_pos_list = Vector{Vector{Int}}(undef, toric_num)
    zstab_pos_list = Vector{Vector{Int}}(undef, toric_num)
    ab_list = Vector{Tuple{Int,Int}}(undef, toric_num)
    result_mat = copy(res.standard_matrix)
    column_source = _require_column_transformation(res, "split_product_state")
    column_trans = copy(column_source)
    row_trans = copy(res.row_transformation)
    toric_block_cache = Dict{Int,Any}()
    for i in 1:toric_num
        block = get(toric_block_cache, sizes[i], nothing)
        if isnothing(block)
            block = factor_toric_block(sizes[i], u, v, R)
            toric_block_cache[sizes[i]] = block
        end
        column_transform_inv, row_transform_inv, a, b = block
        ab_list[i] = (a, b)
        qubit_pos = [
            collect(used_product_state_num+1:used_product_state_num+sizes[i])...,
            collect(product_state_num+used_product_state_num+1:product_state_num+used_product_state_num+sizes[i])...,
            i*2+2*product_state_num-1,
            2*i+2*product_state_num,
        ]
        qubit_pos2 = [qubit_pos...,(qubit_pos .+ stab_num)...]

        z_stab_pos = [
            collect(used_product_state_num+1:used_product_state_num+sizes[i])...,
            product_state_num+i,
        ]
        x_stab_pos = [
            collect(stab_num ÷2+used_product_state_num+1:stab_num ÷2+used_product_state_num+sizes[i])...,
            product_state_num+i+stab_num ÷2,
        ]
        stab_pos = [z_stab_pos...,x_stab_pos...]

        result_mat[stab_pos,qubit_pos2] = row_transform_inv*result_mat[stab_pos,qubit_pos2]*column_transform_inv

        column_update = identity_matrix(R,size(column_source)[1])
        column_update[qubit_pos2,qubit_pos2] = column_transform_inv
        column_trans = column_trans*column_update

        row_update = identity_matrix(R,size(res.row_transformation)[1])
        row_update[stab_pos,stab_pos] = row_transform_inv
        row_trans = row_update*row_trans
        used_product_state_num += sizes[i]

        qubit_pos_list[i] = qubit_pos
        xstab_pos_list[i] = x_stab_pos
        zstab_pos_list[i] = z_stab_pos
    end
    return (;
        result_matrix=result_mat,
        column_transformation=column_trans,
        row_transformation=row_trans,
        qubit_pos_list=qubit_pos_list,
        xstab_pos_list=xstab_pos_list,
        zstab_pos_list=zstab_pos_list,
        toric_cg_size=ab_list,
        zproduct_state_pos=[(i,i) for i in used_size+1:product_state_num],
        xproduct_state_pos=[
            (i+stab_num ÷ 2 - product_state_num,i)
            for i in product_state_num+used_size+1:product_state_num+product_state_num
        ],
    )
end
