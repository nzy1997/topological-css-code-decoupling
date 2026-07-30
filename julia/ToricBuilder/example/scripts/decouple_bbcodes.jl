using ToricBuilder
using Oscar

F = GF(2)
Rxy, (x, y) = laurent_polynomial_ring(F, ["x", "y"])

ab_list1 = [
    (x * y, x * y),
    (y, x),
    (x^-1, y^-1),
]

ab_list2 = [
    (x * y, x * y),
    (x^-1 * y, x * y),
    (x^2, x^2),
]

ab_list3 = [
    (x * y, x * y),
    # (y, x),
    (x^-1 * y, x * y),
    (x^2, x^2),
    (x^-1, y^-1),
    (x * y, x * y^-1),
    (x^-1, x^3 * y^2),
    (y^-2, x^-2),
    (y^-2, x^2),
    (x^-1 * y, x^-1 * y^-1),
    (x^-2 * y^-1, x^2 * y),
    (x^-1 * y^3, x^3 * y^-1),
    (x^-2, x^-2 * y^2),
    (x^-1 * y^-3, x^3 * y^-1),
    (x^-2 * y, x * y^-2),
    (x^-1 * y^2, x^-2 * y^-1),
    (x^-3 * y, x^3 * y^2),
    (x^-3 * y, x^-5),
    (x^-2 * y, x * y^2),
    (x^-1 * y^-2, x * y^-1),
    (y^2, x^-4 * y),
    (x^-1 * y^2, y^-4),
    (y^-4, x^4),
    (x^-8 * y, x^6 * y^2),
    (x^2 * y^3, x^4 * y),
    (x^2 * y^3, x^2 * y^-2),
    (x^-4, x^-3 * y^2),
    (x^-3 * y, x^-1 * y^-2),
    (x^-3 * y^2, x^-3 * y^-1),
    (x^-2 * y^-5, x^-1 * y^-3),
    (x^-6 * y^-1, x^5),
    (x^3 * y, x^2 * y^-2),
    (x^-2 * y, x^-3 * y^-2),
    (x^-3 * y^-1, x^2 * y^-2),
    (x^-1 * y^-3, y^-6),
    (x^-8 * y^-1, x^5 * y),
    (x * y^-5, x * y^4),
    (x^-1 * y^-1, x^5),
    (x * y^3, x^2 * y^-2),
    (x^-3 * y, x * y^-3),
    (x^-1 * y^-4, x^-3 * y^3),
    (x^-1 * y^-2, x^2 * y^-1),
    (x^3 * y^2, x^-4 * y^4),
    (x^-1 * y^3, x * y^3),
    (x^3 * y^-4, x^-1 * y^-3),
    (x^-3 * y^2, x^-4 * y^-1),
    (x^-6 * y^2, x^2 * y^5),
    (x^-4, x^-1 * y^-3),
    (x^-2 * y^2, x^-1 * y^-2),
    (x^2 * y^2, x^-4 * y),
    (x^-1 * y^3, x^3),
    (x^-3 * y^-2, x^-1 * y^-3),
    (x^3 * y^-3, x^4),
    (x^-4 * y^-3, x^3 * y^-1),
    (x^-2 * y^3, x^2 * y^3),
]

function _latexify_poly(s::AbstractString)
    s = replace(s, r"\^(-?\d+)" => s"^{\1}")
    return replace(s, "*" => " \\cdot ")
end

function fmt_math(v)
    if v === "-" || v === nothing
        return "-"
    end
    return "\$" * _latexify_poly(string(v)) * "\$"
end

function fmt_cell(v)
    if v === "-" || v === nothing
        return "-"
    end
    return string(v)
end

function _markdown_row_from_decoupled_case(a, b, decoupled_case::DecoupledToricCase)
    result = decoupled_case.transfer_result
    l_val = isnothing(result) ? "-" : get(result, :l, "-")
    solving_time = isnothing(result) ? "-" : round(get(result, :solving_time, 0.0), digits=3)

    return "| $(fmt_math(a)) | $(fmt_math(b)) | $(fmt_math(l_val)) | $(fmt_math(_result_value(result, :u_rel))) | $(fmt_math(_result_value(result, :v_rel))) | $(fmt_cell(_result_value(result, :area))) | $(fmt_cell(_result_value(result, :A_size))) | $(fmt_cell(_result_value(result, :product_state_num))) | $(fmt_cell(_result_value(result, :toric_num))) | $(fmt_cell(solving_time)) | $(fmt_math(_result_value(result, :max_ele_psi_1_inverse))) | $(fmt_cell(_result_value(result, :max_degree_psi_1_inverse))) | $(fmt_math(_result_value(result, :max_ele_psi_1))) | $(fmt_cell(_result_value(result, :max_degree_psi_1))) | $(fmt_cell(_result_value(result, :max_column_monomial_count_psi_1_inverse))) | $(fmt_cell(_result_value(result, :max_column_monomial_count_psi_1))) | $(fmt_cell(decoupled_case.case_id)) |"
end

function _result_value(result, key)
    if isnothing(result)
        return "-"
    end
    return get(result, key, "-")
end

function _markdown_row_has_completed_result(cells)
    length(cells) >= 18 || return false
    return all(idx -> strip(cells[idx]) != "-", (8, 9, 10, 12))
end

function _markdown_rows_by_case_id(path)
    rows = Dict{String, String}()
    if isnothing(path) || !isfile(path)
        return rows
    end

    for line in eachline(path)
        startswith(line, "|") || continue
        startswith(line, "|---") && continue
        cells = split(line, "|")
        length(cells) >= 3 || continue
        case_id = strip(cells[end-1])
        if isempty(case_id) || case_id == "Case ID"
            continue
        end
        _markdown_row_has_completed_result(cells) || continue
        rows[case_id] = line
    end

    return rows
end

function _case_id(idx::Int)
    return string("case_", lpad(string(idx), 3, '0'))
end

function _build_case_metadata(idx::Int, a, b)
    return Dict(
        "case_index" => idx,
        "a" => string(a),
        "b" => string(b),
    )
end

function _decoupled_cache_satisfies_request(decoupled_case::DecoupledToricCase; compute_inverse::Bool, capture_debug::Bool)
    capture_debug && isnothing(decoupled_case.debug_result) && return false
    !compute_inverse && return true

    result = decoupled_case.transfer_result
    return !isnothing(result) &&
           hasproperty(result, :psi_1) &&
           !isnothing(result.psi_1) &&
           hasproperty(result, :column_transformation) &&
           !isnothing(result.column_transformation)
end

function _load_or_build_decoupled_case(case_id, poly_vector, cache_path; metadata, overwrite::Bool, show_progress::Bool, compute_inverse::Bool, capture_debug::Bool=false, check_cache::Bool=false, kwargs...)
    overwrite_cache = overwrite
    if isfile(cache_path) && !overwrite
        decoupled_case = load_decoupled_toric_case(cache_path)
        if !check_cache || _decoupled_cache_satisfies_request(decoupled_case; compute_inverse=compute_inverse, capture_debug=capture_debug)
            return decoupled_case
        end
        overwrite_cache = true
    end

    decoupled_case = build_decoupled_toric_case(
        case_id,
        poly_vector;
        metadata=metadata,
        show_progress=show_progress,
        compute_inverse=compute_inverse,
        capture_debug=capture_debug,
        capture_exceptions=true,
        kwargs...,
    )
    save_decoupled_toric_case(cache_path, decoupled_case; overwrite=overwrite_cache)
    return decoupled_case
end

function _warmup_decoupling(poly_vector; show_progress::Bool, compute_inverse::Bool, max_area::Int, max_l::Int)
    redirect_stdout(devnull) do
        redirect_stderr(devnull) do
            build_decoupled_toric_case(
                "__warmup__",
                poly_vector;
                metadata=Dict{String, Any}("purpose" => "timing warmup"),
                show_progress=show_progress,
                compute_inverse=compute_inverse,
                capture_exceptions=true,
                max_area=max_area,
                max_l=max_l,
            )
        end
    end
    return nothing
end

function run_and_save(
    list;
    results_path=joinpath(dirname(dirname(dirname(dirname(@__DIR__)))), "build", "reproduction", "bb_code_decoupling_results.md"),
    cache_dir=joinpath(dirname(dirname(dirname(dirname(@__DIR__)))), "build", "reproduction", "decouple_bbcodes_cache_v3"),
    overwrite::Bool=false,
    show_progress::Bool=false,
    compute_inverse::Bool=false,
    capture_debug::Bool=false,
    check_cache::Bool=false,
    reuse_rows_path=nothing,
    warmup::Bool=true,
    max_area::Int=2000,
    max_l::Int=10000,
)
    mkpath(cache_dir)
    mkpath(dirname(results_path))

    reused_rows = _markdown_rows_by_case_id(reuse_rows_path)

    first_case = iterate(list)
    if warmup && !isnothing(first_case)
        (a, b), _ = first_case
        _warmup_decoupling(
            [1 + x + a, 1 + y + b];
            show_progress=show_progress,
            compute_inverse=compute_inverse,
            max_area=max_area,
            max_l=max_l,
        )
    end

    open(results_path, "w") do io
        println(io, "| \$f=1+x+...\$ | \$g=1+y+...\$ | \$L\$ | \$u\$ | \$v\$ | Area | Matrix size after CG | Product state num | Toric num | Time(s) | Maximum_term_psi_1_inverse | Maximum_degree_psi_1_inverse | Maximum_term_psi_1 | Maximum_degree_psi_1 | Maximum_column_monomial_count_psi_1_inverse | Maximum_column_monomial_count_psi_1 | Case ID |")
        println(io, "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")

        for (idx, (a, b)) in enumerate(list)
            poly_vector = [1 + x + a, 1 + y + b]
            case_id = _case_id(idx)
            cache_path = decoupled_toric_case_path(cache_dir, case_id)
            reused_row = get(reused_rows, case_id, nothing)
            if !overwrite && !check_cache && isfile(cache_path) && !isnothing(reused_row)
                println(io, reused_row)
                flush(io)
                continue
            end

            decoupled_case = _load_or_build_decoupled_case(
                case_id,
                poly_vector,
                cache_path;
                metadata=_build_case_metadata(idx, a, b),
                overwrite=overwrite,
                show_progress=show_progress,
                compute_inverse=compute_inverse,
                capture_debug=capture_debug,
                check_cache=check_cache,
                max_area=max_area,
                max_l=max_l,
            )
            println(io, _markdown_row_from_decoupled_case(a, b, decoupled_case))
            flush(io)
        end
    end
end

function main(; list=ab_list2, kwargs...)
    return run_and_save(list; kwargs...)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
