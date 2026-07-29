function excitation_matrix(poly_vec::Vector)
    poly_num = length(poly_vec)
    R = parent(poly_vec[1])
    A = zero_matrix(R, 2, 2 * poly_num)
    for i in 1:poly_num
        A[1, i] = poly_vec[i]
        A[2, 2 * poly_num - i + 1] = laurent_conjugate(poly_vec[i])
    end
    return A
end

function _prepare_transfer_to_toric_form_data(poly_vec::Vector; max_area=150, max_l=200)
    t_start = time()

    Rl = parent(poly_vec[1])
    (x, y) = gens(Rl)

    period_data = _matrix_period_data(poly_vec; max_int=max_l)
    l = period_data.l
    if isnothing(l)
        return (
            ;
            status=:invalid_l,
            l=l,
            solving_time=time() - t_start,
            area=nothing,
            u_rel=nothing,
            v_rel=nothing,
            ideal_triangle=nothing,
        )
    end

    L_table = _make_L_table_from_actions(period_data.rep.Ax, period_data.rep.Ay, l)
    vertices, cg_area = ToricBuilder.find_minimum_triangle(L_table)
    cg_area *= 2

    R_new, (u, v) = laurent_polynomial_ring(Oscar.GF(2), ["u", "v"])
    u_rel = x^(vertices[2][1] - 1) * y^(vertices[2][2] - 1)
    v_rel = x^(vertices[3][1] - 1) * y^(vertices[3][2] - 1)
    if cg_area > max_area
        return (
            ;
            status=:coarse_grained_area_too_large,
            l=l,
            solving_time=time() - t_start,
            area=cg_area,
            u_rel=u_rel,
            v_rel=v_rel,
            ideal_triangle=nothing,
            vertices=vertices,
        )
    end

    points = find_parallelogram_lattice_points(vertices...)
    remainders = [x^(p[1] - 1) * y^(p[2] - 1) for p in points]
    A = coarse_graining_with_replace_variable(
        excitation_matrix(poly_vec),
        [u_rel, v_rel],
        Rl,
        R_new,
        [x, y],
        [u, v],
        remainders,
    )

    return (
        ;
        status=:ok,
        solving_time=time() - t_start,
        l=l,
        area=cg_area,
        u_rel=u_rel,
        v_rel=v_rel,
        ideal_triangle=nothing,
        vertices=vertices,
        points=points,
        remainders=remainders,
        L_table=L_table,
        Rl=Rl,
        R_new=R_new,
        A=A,
    )
end

function build_toric_form(poly_vec::Vector; max_area=150, max_l=200, show_progress::Bool=true, compute_inverse::Bool=false)
    if show_progress
        println("-"^80)
        println("Transfer to toric form:")
        display(poly_vec[1])
        display(poly_vec[2])
        println("-"^80)
    end

    t_start = time()
    prep = _prepare_transfer_to_toric_form_data(poly_vec; max_area=max_area, max_l=max_l)
    prep.status === :no_triangle && return (; solving_time=time() - t_start)
    if prep.status === :ideal_area_too_large
        return (; solving_time=time() - t_start, area=prep.area, v_rel=prep.v_rel, u_rel=prep.u_rel)
    end
    if prep.status === :invalid_l
        return (; l=prep.l, solving_time=time() - t_start)
    end
    if prep.status === :coarse_grained_area_too_large
        return (; l=prep.l, solving_time=time() - t_start, area=prep.area, v_rel=prep.v_rel, u_rel=prep.u_rel)
    end

    if show_progress
        @show prep.ideal_triangle
        @show prep.l
        display(prep.L_table)
        @show prep.area
    end
    local A = prep.A
    local v_rel = prep.v_rel
    local u_rel = prep.u_rel
    local pr = prep.remainders
    local Rl = prep.Rl
    local R_new = prep.R_new
    local x, y = gens(Rl)
    local u, v = gens(R_new)
    if show_progress
        println("We set u = $(u_rel), v = $(v_rel).")
        println("The possible remainders are: $(pr)")
    end

    res = to_toric_form(A; show_progress=show_progress, compute_inverse=compute_inverse)

    phi_1_original = replace_variable_inv_mat(
        res.phi_1,
        [u_rel, v_rel],
        Rl,
        R_new,
        [x, y],
        [u, v],
    )
    max_ele_phi_1, max_degree_phi_1 = maximum_deg_term(phi_1_original)
    max_column_monomial_count_phi_1 = max_column_monomial_count(phi_1_original)

    phi_1_inv_original = nothing
    max_ele_phi_1_inv = nothing
    max_degree_phi_1_inv = nothing
    max_column_monomial_count_phi_1_inv = nothing
    if !isnothing(res.phi_1_inv)
        phi_1_inv_original = replace_variable_inv_mat(
            res.phi_1_inv,
            [u_rel, v_rel],
            Rl,
            R_new,
            [x, y],
            [u, v],
        )
        max_ele_phi_1_inv, max_degree_phi_1_inv = maximum_deg_term(phi_1_inv_original)
        max_column_monomial_count_phi_1_inv = max_column_monomial_count(phi_1_inv_original)
    end
    if show_progress
        @show max_column_monomial_count_phi_1 max_column_monomial_count_phi_1_inv
    end

    return (;
        res...,
        l=prep.l,
        area=prep.area,
        v_rel=v_rel,
        u_rel=u_rel,
        solving_time=time() - t_start,
        A_size=size(A),
        phi_1_size=size(res.phi_1),
        phi_1_original=phi_1_original,
        phi_1_inv_original=phi_1_inv_original,
        max_ele_phi_1=max_ele_phi_1,
        max_degree_phi_1=max_degree_phi_1,
        max_column_monomial_count_phi_1=max_column_monomial_count_phi_1,
        max_ele_phi_1_inv=max_ele_phi_1_inv,
        max_degree_phi_1_inv=max_degree_phi_1_inv,
        max_column_monomial_count_phi_1_inv=max_column_monomial_count_phi_1_inv,
    )
end

function project_to_finite_code(A, l, m)
    F2 = Oscar.GF(2)
    Rl = A.base_ring
    (x, y) = gens(Rl)
    A2 = coarse_graining(A, [x, y], [l, m])
    return [iszero(evaluate(elem, fill(F2(1), n_variables(A.base_ring)))) ? 0 : 1 for elem in A2.entries]
end

function _require_column_transformation(res, function_name::AbstractString)
    isnothing(res.column_transformation) &&
        throw(ArgumentError("$function_name requires column_transformation; rerun build_toric_form(...; compute_inverse=true)."))
    return res.column_transformation
end
