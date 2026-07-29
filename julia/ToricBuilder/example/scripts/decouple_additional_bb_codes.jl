using ToricBuilder
using Oscar

F = GF(2)
Rxy, (x, y) = laurent_polynomial_ring(F, ["x", "y"])

struct AdditionalBBCodeCase
    case_id::String
    label::String
    poly_vector::Vector
    definition::String
    metadata::Dict{String, Any}
end

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

function fmt_poly_pair(poly_vector)
    return string(fmt_math(poly_vector[1]), "<br>", fmt_math(poly_vector[2]))
end

function _default_definition(poly_vector)
    return "source polynomial pair:<br>" * fmt_poly_pair(poly_vector)
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

function _markdown_row_from_decoupled_case(case::AdditionalBBCodeCase, decoupled_case::DecoupledToricCase)
    result = decoupled_case.transfer_result
    l_val = isnothing(result) ? "-" : get(result, :l, "-")
    solving_time = isnothing(result) ? "-" : round(get(result, :solving_time, 0.0), digits=3)

    return "| $(fmt_cell(case.label)) | $(fmt_cell(case.definition)) | $(fmt_math(l_val)) | $(fmt_math(_result_value(result, :u_rel))) | $(fmt_math(_result_value(result, :v_rel))) | $(fmt_cell(_result_value(result, :area))) | $(fmt_cell(_result_value(result, :A_size))) | $(fmt_cell(_result_value(result, :product_state_num))) | $(fmt_cell(_result_value(result, :toric_num))) | $(fmt_cell(solving_time)) | $(fmt_math(_result_value(result, :max_ele_phi_1))) | $(fmt_cell(_result_value(result, :max_degree_phi_1))) | $(fmt_math(_result_value(result, :max_ele_phi_1_inv))) | $(fmt_cell(_result_value(result, :max_degree_phi_1_inv))) | $(fmt_cell(_result_value(result, :max_column_monomial_count_phi_1))) | $(fmt_cell(_result_value(result, :max_column_monomial_count_phi_1_inv))) | $(fmt_cell(decoupled_case.case_id)) |"
end

function _case_metadata(case::AdditionalBBCodeCase)
    metadata = copy(case.metadata)
    metadata["case_id"] = case.case_id
    metadata["label"] = case.label
    metadata["f"] = string(case.poly_vector[1])
    metadata["g"] = string(case.poly_vector[2])
    metadata["definition"] = case.definition
    return metadata
end

function _table_case(case_id, label, f, g; metadata=Dict{String, Any}(), definition=nothing)
    poly_vector = [f, g]
    case_definition = isnothing(definition) ? _default_definition(poly_vector) : string(definition)
    return AdditionalBBCodeCase(string(case_id), string(label), poly_vector, case_definition, Dict{String, Any}(metadata))
end

function _bunny_case(idx::Int, parameters::AbstractString, l::Int, m::Int, f, g, connectivity::AbstractString)
    case_id = "bunny_" * lpad(string(idx), 3, '0')
    label = "arXiv:2606.22853 $(parameters) l=$(l) m=$(m) $(connectivity)"
    metadata = Dict{String, Any}(
        "source" => "arXiv:2606.22853",
        "source_table" => "code_poly",
        "parameters" => parameters,
        "l" => l,
        "m" => m,
        "connectivity" => connectivity,
    )
    return _table_case(case_id, label, f, g; metadata=metadata)
end

function bunny_code_cases()
    rows = [
        (1, "[[12,2,3]]", 3, 2, 1 + x, x + y, "Hexagon"),
        # Same infinite-plane polynomial pair as [[12,2,3]].
        # ("[[16,2,4]]", 4, 2, 1 + x, x + y, "Hexagon"),
        # ("[[24,2,4]]", 4, 3, 1 + x, x + y, "Hexagon"),
        # ("[[30,2,5]]", 5, 3, 1 + x, x + y, "Hexagon"),
        (5, "[[10,2,3]]", 1, 5, x + y^2, y + y^2, "Square"),
        # Same infinite-plane polynomial pair as [[10,2,3]].
        # ("[[28,2,5]]", 2, 7, x + y^2, y + y^2, "Square"),
        (7, "[[18,4,3]]", 3, 3, x + x * y + x^2 * y, x * y + y^2, "Hexagon"),
        # Same infinite-plane polynomial pair as the previous row.
        # ("[[36,4,4]]", 3, 6, x + x * y + x^2 * y, x * y + y^2, "Hexagon"),
        (9, "[[36,4,4]]", 6, 3, 1 + x + x^2, x + y, "Hexagon"),
        # Same infinite-plane polynomial pair as the previous row.
        # ("[[54,4,5]]", 9, 3, 1 + x + x^2, x + y, "Hexagon"),
        (11, "[[18,4,3]]", 3, 3, 1 + x * y + x^2 * y, x + x * y, "Square"),
        # Same infinite-plane polynomial pair as the previous row.
        # ("[[24,4,4]]", 3, 4, 1 + x * y + x^2 * y, x + x * y, "Square"),
        (13, "[[30,4,5]]", 3, 5, 1 + x + x^2 * y^2, x * y + x * y^2, "Square"),
        (14, "[[12,4,3]]", 3, 2, x + y + x * y + x^2, x * y + x^2 * y, "Hexagon"),
        (15, "[[18,6,3]]", 3, 3, y + x * y + x^2 + x^3, x * y + x^3 * y, "Hexagon"),
        (16, "[[16,4,4]]", 4, 2, x * y + x^2 * y + x^3 + x^3 * y, y^2 + x^3 * y, "Hexagon"),
        (17, "[[12,4,3]]", 6, 1, 1 + x + x^2 + x^3, 1 + x^2, "Square"),
        (18, "[[36,8,4]]", 3, 6, 1 + x + x^2 * y^3, x * y + x * y^2 + x * y^3, "Square"),
        (19, "[[40,4,5]]", 2, 10, x + y, y + y^2 + y^3 + y^4, "Square"),
        (20, "[[70,10,5]]", 5, 7, 1 + y^2 + x + x * y^2, x * y + x * y^2, "Square"),
    ]
    return [_bunny_case(row...) for row in rows]
end

function _sigma(poly)
    R = parent(poly)
    x_var, y_var = gens(R)
    return hom(R, R, [y_var, x_var])(poly)
end

function _reflected_pairing(w, wprime)
    f, g, h = w
    fprime, gprime, hprime = wprime
    return (
        ToricBuilder.laurent_conjugate(f) * _sigma(fprime) +
        ToricBuilder.laurent_conjugate(g) * _sigma(hprime) +
        ToricBuilder.laurent_conjugate(h) * _sigma(gprime)
    )
end

function _sbb_bb_pair(w1, w2)
    f1, g1, h1 = w1
    f2, g2, h2 = w2
    a = _reflected_pairing(w1, w1)
    u = g1 * f1^-1
    v = h1 * f1^-1
    a_bar = ToricBuilder.laurent_conjugate(a)
    return a_bar * (g2 + u * f2), a_bar * (h2 + v * f2)
end

function _sbb_case(
    idx::Int,
    parameters::AbstractString,
    w1,
    w2,
    a1::Tuple{Int, Int},
    a2::Tuple{Int, Int},
)
    f, g = _sbb_bb_pair(w1, w2)
    case_id = "sbb_" * lpad(string(idx), 3, '0')
    label = "arXiv:2605.04151 $(parameters) a1=$(a1) a2=$(a2)"
    metadata = Dict{String, Any}(
        "source" => "arXiv:2605.04151",
        "source_table" => "tab: more examples",
        "parameters" => parameters,
        "a1" => a1,
        "a2" => a2,
        "gauge_triple_1" => string.(w1),
        "gauge_triple_2" => string.(w2),
    )
    definition = string(
        "source SBB gauge triples: `", w1, "`<br>",
        "`", w2, "`<br>",
        "decouple BB pair:<br>",
        fmt_poly_pair([f, g]),
    )
    return _table_case(case_id, label, f, g; metadata=metadata, definition=definition)
end

function sbb_code_cases()
    return [
        _sbb_case(
            1,
            "[[27,6,3]]",
            (y^2, y, 1 + x * y),
            (1 + x * y + x^2, 1, 0),
            (0, 3),
            (3, 0),
        ),
        _sbb_case(
            2,
            "[[60,10,4]]",
            (y^2, x, 1 + x * y),
            (x^2 + x^2 * y^2, y^2 + x * y, 0),
            (0, 5),
            (4, -1),
        ),
        _sbb_case(
            3,
            "[[75,10,5]]",
            (x^2, y^2, x + x^2 * y),
            (1 + y^2, x + y, 0),
            (0, 5),
            (5, 0),
        ),
        _sbb_case(
            4,
            "[[90,12,5]]",
            (y^2, x^2, 1 + x^2 * y^2),
            (1 + x^3 * y, y^2 + x * y, 0),
            (0, 6),
            (5, 2),
        ),
        # Same infinite-plane SBB gauge triples, hence same decouple BB pair, as [[90,12,5]].
        # (
        #     "[[108,12,6]]",
        #     (y^2, x^2, 1 + x^2 * y^2),
        #     (1 + x^3 * y, y^2 + x * y, 0),
        #     (0, 9),
        #     (4, 1),
        # ),
        _sbb_case(
            6,
            "[[126,14,6]]",
            (y^2, x^2, x + x^2 * y),
            (x + x^3, y^2 + x * y, 0),
            (0, 7),
            (6, -1),
        ),
    ]
end

function _cover_bb_case(
    suffix::AbstractString,
    source_table::AbstractString,
    family::AbstractString,
    parameters::AbstractString,
    l::Int,
    m::Int,
    h::Int,
    check_weight::Int,
    f,
    g,
    kd2_over_n::AbstractString,
)
    case_id = "cover_" * string(suffix)
    label = "arXiv:2511.13560 $(parameters) l=$(l) m=$(m) h=$(h) w=$(check_weight)"
    metadata = Dict{String, Any}(
        "source" => "arXiv:2511.13560",
        "source_table" => source_table,
        "family" => family,
        "parameters" => parameters,
        "l" => l,
        "m" => m,
        "h" => h,
        "check_weight" => check_weight,
        "kd2_over_n" => kd2_over_n,
    )
    definition = string(
        "source cover code: ", parameters,
        ", l=", l,
        ", m=", m,
        ", h=", h,
        ", check weight=", check_weight,
        "<br>",
        _default_definition([f, g]),
    )
    return _table_case(case_id, label, f, g; metadata=metadata, definition=definition)
end

function cover_bb_code_cases()
    rows = [
        (
            "w6_k12_001",
            "tab:k-12-codes",
            "weight 6 k=12 fixed-k cover sequence",
            "[[648,12,d<=36]]",
            18,
            18,
            9,
            6,
            x^3 + y^13 + x^12 * y^2,
            y^3 + x^7 * y^12 + x^14 * y^6,
            "<= 24",
        ),
        (
            "w6_k12_002",
            "tab:k-12-codes",
            "weight 6 k=12 fixed-k cover sequence",
            "[[576,12,d<=32]]",
            24,
            12,
            8,
            6,
            x^21 + y^7 + y^2,
            y^9 + x^13 + x^14,
            "<= 21.3",
        ),
        (
            "w6_k12_003",
            "tab:k-12-codes",
            "weight 6 k=12 fixed-k cover sequence",
            "[[504,12,d<=30]]",
            42,
            6,
            7,
            6,
            x^27 + y + y^2,
            y^3 + x^37 + x^38,
            "<= 21.4",
        ),
        (
            "w6_k12_004",
            "tab:k-12-codes",
            "weight 6 k=12 fixed-k cover sequence",
            "[[432,12,d<=24]]",
            18,
            12,
            6,
            6,
            x^15 + y^7 + y^2,
            y^9 + x^7 + x^14,
            "<= 16",
        ),
        (
            "w6_k12_005",
            "tab:k-12-codes",
            "weight 6 k=12 fixed-k cover sequence",
            "[[360,12,24]]",
            30,
            6,
            5,
            6,
            x^9 + y + y^2,
            y^3 + x^25 + x^26,
            "19.2",
        ),
        # [[288,12,18]] is already represented by double_gross_case().
        (
            "w6_k12_007",
            "tab:k-12-codes",
            "weight 6 k=12 fixed-k cover sequence",
            "[[216,12,12]]",
            18,
            6,
            3,
            6,
            x^3 + y + y^2,
            y^3 + x + x^2,
            "8",
        ),
        # The [[144,12,12]] and [[72,12,6]] rows have the same infinite-plane pair as w6_k12_007.
        (
            "w6_k10_001",
            "tab:k-10-codes",
            "weight 6 k=10 fixed-k cover sequence",
            "[[434,10,d<=26]]",
            31,
            7,
            7,
            6,
            1 + x^6 * y + x^27 * y^4,
            1 + x^15 * y^6 + x^24 * y^3,
            "<= 15.6",
        ),
        (
            "w6_k10_002",
            "tab:k-10-codes",
            "weight 6 k=10 fixed-k cover sequence",
            "[[372,10,d<=24]]",
            31,
            6,
            6,
            6,
            1 + x^6 * y^2 + x^27,
            1 + x^15 * y + x^24 * y^4,
            "<= 15.5",
        ),
        (
            "w6_k10_003",
            "tab:k-10-codes",
            "weight 6 k=10 fixed-k cover sequence",
            "[[310,10,d<=22]]",
            31,
            5,
            5,
            6,
            1 + x^6 * y^2 + x^27 * y,
            y^2 + x^15 + x^24,
            "<= 15.6",
        ),
        (
            "w6_k10_004",
            "tab:k-10-codes",
            "weight 6 k=10 fixed-k cover sequence",
            "[[248,10,d<=18]]",
            31,
            4,
            4,
            6,
            1 + x^6 * y + x^27,
            y^2 + x^15 * y^3 + x^24,
            "<= 13.1",
        ),
        (
            "w6_k10_005",
            "tab:k-10-codes",
            "weight 6 k=10 fixed-k cover sequence",
            "[[186,10,14]]",
            31,
            3,
            3,
            6,
            y + x^6 * y^2 + x^27,
            1 + x^15 * y + x^24,
            "10.5",
        ),
        (
            "w6_k10_006",
            "tab:k-10-codes",
            "weight 6 k=10 fixed-k cover sequence",
            "[[124,10,10]]",
            31,
            2,
            2,
            6,
            y + x^6 + x^27,
            1 + x^15 + x^24,
            "8.1",
        ),
        (
            "w6_k10_007",
            "tab:k-10-codes",
            "weight 6 k=10 fixed-k cover sequence",
            "[[62,10,6]]",
            31,
            1,
            1,
            6,
            1 + x^6 + x^27,
            1 + x^15 + x^24,
            "5.8",
        ),
        (
            "w6_k8_001",
            "tab:k-8-codes",
            "weight 6 k=8 fixed-k cover sequence",
            "[[198,8,d<=16]]",
            33,
            3,
            11,
            6,
            x^12 + y + y^2,
            1 + x + x^8,
            "<= 10.3",
        ),
        (
            "w6_k8_002",
            "tab:k-8-codes",
            "weight 6 k=8 fixed-k cover sequence",
            "[[180,8,d<=16]]",
            15,
            6,
            10,
            6,
            x^9 + y + x^6 * y^5,
            x^6 * y^3 + x + x^2,
            "<= 11.4",
        ),
        (
            "w6_k8_003",
            "tab:k-8-codes",
            "weight 6 k=8 fixed-k cover sequence",
            "[[162,8,14]]",
            27,
            3,
            9,
            6,
            1 + y + x^6 * y^2,
            1 + x^25 + x^20,
            "9.7",
        ),
        (
            "w6_k8_004",
            "tab:k-8-codes",
            "weight 6 k=8 fixed-k cover sequence",
            "[[144,8,12]]",
            24,
            3,
            8,
            6,
            1 + y + x^21 * y^2,
            1 + x^22 + x^17,
            "8",
        ),
        (
            "w6_k8_005",
            "tab:k-8-codes",
            "weight 6 k=8 fixed-k cover sequence",
            "[[126,8,10]]",
            21,
            3,
            7,
            6,
            x^9 + y + y^2,
            1 + x + x^8,
            "6.3",
        ),
        # [[108,8,10]] duplicates w6_k12_007.
        (
            "w6_k8_007",
            "tab:k-8-codes",
            "weight 6 k=8 fixed-k cover sequence",
            "[[90,8,10]]",
            15,
            3,
            5,
            6,
            x^9 + y + y^2,
            1 + x^7 + x^2,
            "8.9",
        ),
        (
            "w6_k8_008",
            "tab:k-8-codes",
            "weight 6 k=8 fixed-k cover sequence",
            "[[72,8,8]]",
            12,
            3,
            4,
            6,
            x^9 + y + y^2,
            1 + x^4 + x^11,
            "7.1",
        ),
        (
            "w6_k8_009",
            "tab:k-8-codes",
            "weight 6 k=8 fixed-k cover sequence",
            "[[54,8,6]]",
            9,
            3,
            3,
            6,
            x^3 + y + y^2,
            1 + x + x^2,
            "5.3",
        ),
        # [[36,8,4]] duplicates w6_k8_009.
        (
            "w6_k8_011",
            "tab:k-8-codes",
            "weight 6 k=8 fixed-k cover sequence",
            "[[18,8,2]]",
            3,
            3,
            1,
            6,
            1 + y + y^2,
            1 + x + x^2,
            "1.8",
        ),
        (
            "w6_k6_001",
            "tab:k-6-codes",
            "weight 6 k=6 fixed-k cover sequence",
            "[[154,6,d<=16]]",
            7,
            11,
            11,
            6,
            1 + x^2 * y^3 + x^3 * y^4,
            1 + x^2 * y^8 + x^3 * y^7,
            "<= 10.0",
        ),
        (
            "w6_k6_002",
            "tab:k-6-codes",
            "weight 6 k=6 fixed-k cover sequence",
            "[[140,6,14]]",
            7,
            10,
            10,
            6,
            1 + x^2 * y^5 + x^3 * y^9,
            1 + x^2 * y^6 + x^3 * y^3,
            "8.4",
        ),
        (
            "w6_k6_003",
            "tab:k-6-codes",
            "weight 6 k=6 fixed-k cover sequence",
            "[[126,6,14]]",
            7,
            9,
            9,
            6,
            1 + x^2 * y^5 + x^3 * y,
            1 + x^2 + x^3 * y^2,
            "9.3",
        ),
        (
            "w6_k6_004",
            "tab:k-6-codes",
            "weight 6 k=6 fixed-k cover sequence",
            "[[112,6,12]]",
            7,
            8,
            8,
            6,
            1 + x^2 * y^5 + x^3 * y,
            1 + x^2 * y^6 + x^3 * y^5,
            "7.7",
        ),
        (
            "w6_k6_005",
            "tab:k-6-codes",
            "weight 6 k=6 fixed-k cover sequence",
            "[[98,6,12]]",
            7,
            7,
            7,
            6,
            1 + x^2 * y^2 + x^3 * y,
            1 + x^2 + x^3 * y^2,
            "8.8",
        ),
        (
            "w6_k6_006",
            "tab:k-6-codes",
            "weight 6 k=6 fixed-k cover sequence",
            "[[84,6,10]]",
            7,
            6,
            6,
            6,
            1 + x^2 * y^3 + x^3 * y^2,
            1 + x^2 + x^3 * y,
            "7.1",
        ),
        (
            "w6_k6_007",
            "tab:k-6-codes",
            "weight 6 k=6 fixed-k cover sequence",
            "[[70,6,8]]",
            7,
            5,
            5,
            6,
            y + x^2 * y^4 + x^3 * y,
            y^4 + x^2 + x^3,
            "5.5",
        ),
        (
            "w6_k6_008",
            "tab:k-6-codes",
            "weight 6 k=6 fixed-k cover sequence",
            "[[56,6,8]]",
            7,
            4,
            4,
            6,
            1 + x^2 + x^3 * y^2,
            1 + x^2 * y + x^3,
            "6.9",
        ),
        (
            "w6_k6_009",
            "tab:k-6-codes",
            "weight 6 k=6 fixed-k cover sequence",
            "[[42,6,6]]",
            7,
            3,
            3,
            6,
            1 + x^2 + x^3 * y,
            1 + x^2 + x^3 * y^2,
            "5.1",
        ),
        (
            "w6_k6_010",
            "tab:k-6-codes",
            "weight 6 k=6 fixed-k cover sequence",
            "[[28,6,4]]",
            7,
            2,
            2,
            6,
            y + x^2 + x^3,
            1 + x^2 + x^3,
            "3.4",
        ),
        (
            "w6_k6_011",
            "tab:k-6-codes",
            "weight 6 k=6 fixed-k cover sequence",
            "[[14,6,2]]",
            7,
            1,
            1,
            6,
            1 + x^2 + x^3,
            1 + x^2 + x^3,
            "1.7",
        ),
        (
            "w6_inc_001",
            "tab:increasing-k-codes-wt-6",
            "weight 6 increasing-k cover sequence",
            "[[288,20,6]]",
            12,
            12,
            8,
            6,
            x^9 * y^9 + x^6 * y^7 + y^8,
            1 + x^7 * y^6 + x^2 * y^9,
            "2.5",
        ),
        (
            "w6_inc_002",
            "tab:increasing-k-codes-wt-6",
            "weight 6 increasing-k cover sequence",
            "[[288,16,12]]",
            12,
            12,
            8,
            6,
            x^3 * y^3 + x^6 * y^7 + y^8,
            x^6 + x * y^9 + x^2 * y^9,
            "8",
        ),
        (
            "w6_inc_003",
            "tab:increasing-k-codes-wt-6",
            "weight 6 increasing-k cover sequence",
            "[[252,20,4]]",
            6,
            21,
            7,
            6,
            x^3 * y^18 + y^4 + y^11,
            y^18 + x * y^15 + x^2 * y^9,
            "1.27",
        ),
        (
            "w6_inc_004",
            "tab:increasing-k-codes-wt-6",
            "weight 6 increasing-k cover sequence",
            "[[252,14,12]]",
            6,
            21,
            7,
            6,
            x^3 + y^19 + y^8,
            y^9 + x * y^6 + x^2 * y^15,
            "8",
        ),
        (
            "w8_k14_001",
            "tab:k-14-w-8-codes",
            "weight 8 k=14 fixed-k cover sequence",
            "[[288,14,d<=24]]",
            12,
            12,
            4,
            8,
            x^6 * y^4 + x^11 * y^4 + x^3 * y^6 + x^11 * y^3,
            y^11 + x^2 * y^7 + x^5 * y^11 + x^9 * y^4,
            "<= 28",
        ),
        (
            "w8_k14_002",
            "tab:k-14-w-8-codes",
            "weight 8 k=14 fixed-k cover sequence",
            "[[216,14,d<=20]]",
            18,
            6,
            3,
            8,
            y^4 + x^11 * y^4 + x^3 + x^11 * y^3,
            x^6 * y^5 + x^2 * y + x^5 * y^5 + x^15 * y^4,
            "<= 25.9",
        ),
        (
            "w8_k14_003",
            "tab:k-14-w-8-codes",
            "weight 8 k=14 fixed-k cover sequence",
            "[[144,14,14]]",
            12,
            6,
            2,
            8,
            x^6 * y^4 + x^5 * y^4 + x^3 + x^11 * y^3,
            y^5 + x^8 * y + x^5 * y^5 + x^9 * y^4,
            "19.1",
        ),
        (
            "w8_k14_004",
            "tab:k-14-w-8-codes",
            "weight 8 k=14 fixed-k cover sequence",
            "[[72,14,8]]",
            6,
            6,
            1,
            8,
            y^4 + x^5 * y^4 + x^3 + x^5 * y^3,
            y^5 + x^2 * y + x^5 * y^5 + x^3 * y^4,
            "12.4",
        ),
        (
            "w8_k14_alt_001",
            "tab:k-14-w-8-codes-alt",
            "weight 8 k=14 alternate fixed-k cover sequence",
            "[[256,14,d<=22]]",
            16,
            8,
            4,
            8,
            x^9 * y^3 + y^4 + x^14 + x^3 * y^6,
            x^14 * y + x^4 * y + x^3 + x^13 * y,
            "<= 26.5",
        ),
        (
            "w8_k14_alt_002",
            "tab:k-14-w-8-codes-alt",
            "weight 8 k=14 alternate fixed-k cover sequence",
            "[[192,14,d<=16]]",
            8,
            12,
            3,
            8,
            x * y^3 + y^8 + x^6 * y^4 + x^3 * y^2,
            x^6 * y + x^4 * y^5 + x^3 + x^5 * y,
            "<= 18.7",
        ),
        (
            "w8_k14_alt_003",
            "tab:k-14-w-8-codes-alt",
            "weight 8 k=14 alternate fixed-k cover sequence",
            "[[128,14,12]]",
            8,
            8,
            2,
            8,
            x * y^3 + y^4 + x^6 * y^4 + x^3 * y^6,
            x^6 * y + x^4 * y^5 + x^3 + x^5 * y^5,
            "15.8",
        ),
        (
            "w8_k14_alt_004",
            "tab:k-14-w-8-codes-alt",
            "weight 8 k=14 alternate fixed-k cover sequence",
            "[[64,14,8]]",
            8,
            4,
            1,
            8,
            x * y^3 + 1 + x^6 + x^3 * y^2,
            x^6 * y + x^4 * y + x^3 + x^5 * y,
            "14",
        ),
        (
            "w8_k12_001",
            "tab:k-12-w-8-codes",
            "weight 8 k=12 fixed-k cover sequence",
            "[[160,12,d<=16]]",
            8,
            10,
            5,
            8,
            x * y^3 + y^2 + x^6 * y^8 + x^3 * y^4,
            x^6 * y^7 + x^4 * y + x^3 + x^5 * y^3,
            "<= 19.2",
        ),
        (
            "w8_k12_002",
            "tab:k-12-w-8-codes",
            "weight 8 k=12 fixed-k cover sequence",
            "[[128,12,14]]",
            8,
            8,
            4,
            8,
            x * y^3 + 1 + x^6 * y^6 + x^3 * y^2,
            x^6 * y^7 + x^4 * y^7 + x^3 + x^5 * y,
            "<= 18.4",
        ),
        (
            "w8_k12_003",
            "tab:k-12-w-8-codes",
            "weight 8 k=12 fixed-k cover sequence",
            "[[96,12,10]]",
            8,
            6,
            3,
            8,
            x * y + y^2 + x^6 * y^4 + x^3 * y^4,
            x^6 * y + x^4 * y^5 + x^3 * y^2 + x^5 * y,
            "12.5",
        ),
        (
            "w8_k12_004",
            "tab:k-12-w-8-codes",
            "weight 8 k=12 fixed-k cover sequence",
            "[[64,12,8]]",
            8,
            4,
            2,
            8,
            x * y^3 + 1 + x^6 * y^2 + x^3,
            x^6 * y + x^4 * y + x^3 * y^2 + x^5 * y,
            "12",
        ),
        (
            "w8_k12_005",
            "tab:k-12-w-8-codes",
            "weight 8 k=12 fixed-k cover sequence",
            "[[32,12,4]]",
            8,
            2,
            1,
            8,
            x * y + 1 + x^6 + x^3,
            x^6 * y + x^4 * y + x^3 + x^5 * y,
            "6",
        ),
        (
            "w8_k10_001",
            "tab:k-10-w-8-codes",
            "weight 8 k=10 fixed-k cover sequence",
            "[[168,10,d<=18]]",
            6,
            14,
            7,
            8,
            y^12 + x^5 * y^8 + x^3 * y^4 + x^5 * y^5,
            y^7 + x^2 * y^7 + x^5 * y^5 + x^3 * y^6,
            "<= 19.3",
        ),
        (
            "w8_k10_002",
            "tab:k-10-w-8-codes",
            "weight 8 k=10 fixed-k cover sequence",
            "[[144,10,d<=16]]",
            6,
            12,
            6,
            8,
            1 + x^5 * y^4 + x^3 * y^4 + x^5 * y^3,
            y^5 + x^2 * y^9 + x^5 * y^11 + x^3 * y^4,
            "<= 17.8",
        ),
        (
            "w8_k10_003",
            "tab:k-10-w-8-codes",
            "weight 8 k=10 fixed-k cover sequence",
            "[[120,10,14]]",
            6,
            10,
            5,
            8,
            y^6 + x^5 + x^3 * y^2 + x^5 * y^3,
            y + x^2 * y + x^5 * y^9 + x^3 * y^2,
            "<= 16.3",
        ),
        (
            "w8_k10_004",
            "tab:k-10-w-8-codes",
            "weight 8 k=10 fixed-k cover sequence",
            "[[96,10,12]]",
            12,
            4,
            4,
            8,
            y^2 + x^11 + x^9 + x^5 * y,
            y^3 + x^2 * y^3 + x^5 * y^3 + x^3 * y^2,
            "15",
        ),
        (
            "w8_k10_005",
            "tab:k-10-w-8-codes",
            "weight 8 k=10 fixed-k cover sequence",
            "[[72,10,8]]",
            6,
            6,
            3,
            8,
            y^4 + x^5 + x^3 * y^4 + x^5 * y,
            y^5 + x^2 * y^5 + x^5 * y + x^3 * y^2,
            "8.9",
        ),
        (
            "w8_k10_006",
            "tab:k-10-w-8-codes",
            "weight 8 k=10 fixed-k cover sequence",
            "[[48,10,6]]",
            6,
            4,
            2,
            8,
            y^2 + x^5 + x^3 + x^5 * y,
            y^3 + x^2 * y^3 + x^5 * y^3 + x^3 * y^2,
            "7.5",
        ),
        (
            "w8_k10_007",
            "tab:k-10-w-8-codes",
            "weight 8 k=10 fixed-k cover sequence",
            "[[24,10,4]]",
            6,
            2,
            1,
            8,
            1 + x^5 + x^3 + x^5 * y,
            y + x^2 * y + x^5 * y + x^3,
            "6.7",
        ),
        (
            "w8_k8_001",
            "tab:k-8-w-8-codes",
            "weight 8 k=8 fixed-k cover sequence",
            "[[140,8,d<=16]]",
            7,
            10,
            5,
            8,
            x^4 * y^6 + y^5 + x^5 * y^5 + x^3 * y^9,
            x^5 * y^3 + x^3 * y + x^4 * y + y^8,
            "<= 14.6",
        ),
        (
            "w8_k8_002",
            "tab:k-8-w-8-codes",
            "weight 8 k=8 fixed-k cover sequence",
            "[[112,8,14]]",
            7,
            8,
            4,
            8,
            x^4 * y^4 + y^3 + x^5 * y^7 + x^3 * y^7,
            x^5 * y^7 + x^3 * y^5 + x^4 * y^7 + 1,
            "<= 14",
        ),
        (
            "w8_k8_003",
            "tab:k-8-w-8-codes",
            "weight 8 k=8 fixed-k cover sequence",
            "[[84,8,10]]",
            7,
            6,
            3,
            8,
            x^4 * y^4 + y^5 + x^5 * y^3 + x^3 * y^5,
            x^5 * y + x^3 * y^3 + x^4 * y^3 + y^2,
            "9.5",
        ),
        (
            "w8_k8_004",
            "tab:k-8-w-8-codes",
            "weight 8 k=8 fixed-k cover sequence",
            "[[56,8,8]]",
            7,
            4,
            2,
            8,
            x^4 * y^2 + y + x^5 * y + x^3 * y^3,
            x^5 * y^3 + x^3 * y^3 + x^4 * y^3 + y^2,
            "9.1",
        ),
        (
            "w8_k8_005",
            "tab:k-8-w-8-codes",
            "weight 8 k=8 fixed-k cover sequence",
            "[[28,8,4]]",
            7,
            2,
            1,
            8,
            x^4 + y + x^5 * y + x^3 * y,
            x^5 * y + x^3 * y + x^4 * y + 1,
            "4.6",
        ),
        (
            "w8_inc_001",
            "tab:increasing-k-codes-wt-8",
            "weight 8 increasing-k cover sequence",
            "[[288,36,d<=8]]",
            12,
            12,
            8,
            8,
            x^9 * y^3 + x^8 * y^4 + x^10 + x^7 * y^10,
            x^6 * y^9 + x^4 * y + x^11 + x^9 * y,
            "<= 8",
        ),
        (
            "w8_inc_002",
            "tab:increasing-k-codes-wt-8",
            "weight 8 increasing-k cover sequence",
            "[[288,24,d<=14]]",
            12,
            12,
            8,
            8,
            x^9 * y^11 + x^4 * y^4 + x^10 + x^3 * y^10,
            x^10 * y^5 + x^8 * y^5 + x^11 * y^8 + x * y^5,
            "<= 16.3",
        ),
        (
            "w8_inc_003",
            "tab:increasing-k-codes-wt-8",
            "weight 8 increasing-k cover sequence",
            "[[288,20,d<=22]]",
            12,
            12,
            8,
            8,
            x^5 * y^3 + x^4 * y^8 + x^2 * y^4 + x^7 * y^6,
            x^10 * y + x^8 * y + x^11 * y^4 + x^5 * y^5,
            "<= 33.6",
        ),
        (
            "w8_inc_004",
            "tab:increasing-k-codes-wt-8",
            "weight 8 increasing-k cover sequence",
            "[[224,56,4]]",
            4,
            28,
            7,
            8,
            x * y^7 + y^4 + x^2 * y^4 + x^3 * y^14,
            x^2 * y^25 + y^25 + x^3 + x * y^21,
            "4",
        ),
        (
            "w8_inc_005",
            "tab:increasing-k-codes-wt-8",
            "weight 8 increasing-k cover sequence",
            "[[224,32,8]]",
            4,
            28,
            7,
            8,
            x * y^23 + y^16 + x^2 * y^4 + x^3 * y^18,
            x^2 * y^9 + y^21 + x^3 + x * y^9,
            "9.1",
        ),
        (
            "w8_inc_006",
            "tab:increasing-k-codes-wt-8",
            "weight 8 increasing-k cover sequence",
            "[[224,26,d<=14]]",
            28,
            4,
            7,
            8,
            x^13 * y^3 + x^16 + x^10 + x^11 * y^2,
            x^26 * y + x^16 * y + x^3 + x * y,
            "<= 22.8",
        ),
        # The [[224,20,d<=14]] row has a missing term in the arXiv source polynomial A.
        (
            "w8_inc_008",
            "tab:increasing-k-codes-wt-8",
            "weight 8 increasing-k cover sequence",
            "[[192,30,8]]",
            24,
            4,
            6,
            8,
            x^21 * y^3 + x^4 + x^10 + x^15 * y^2,
            x^22 * y + x^4 * y + x^3 + x^21 * y,
            "10",
        ),
        (
            "w8_inc_009",
            "tab:increasing-k-codes-wt-8",
            "weight 8 increasing-k cover sequence",
            "[[192,26,12]]",
            24,
            4,
            6,
            8,
            x^9 * y^3 + x^12 + x^14 + x^23 * y^2,
            x^10 * y + x^20 * y + x^11 + x * y,
            "19.5",
        ),
        (
            "w8_inc_010",
            "tab:increasing-k-codes-wt-8",
            "weight 8 increasing-k cover sequence",
            "[[192,20,d<=16]]",
            24,
            4,
            6,
            8,
            x^9 * y^3 + x^4 + x^6 + x^19 * y^2,
            x^22 * y + x^8 * y + x^7 + x^5 * y,
            "<= 26.7",
        ),
        (
            "w8_inc_011",
            "tab:increasing-k-codes-wt-8",
            "weight 8 increasing-k cover sequence",
            "[[192,18,d<=16]]",
            24,
            4,
            6,
            8,
            x^21 * y^3 + x^20 + x^2 + x^23 * y^2,
            x^18 * y + y + x^11 + x^17 * y,
            "<= 24",
        ),
        (
            "w8_inc_012",
            "tab:increasing-k-codes-wt-8",
            "weight 8 increasing-k cover sequence",
            "[[160,40,4]]",
            20,
            4,
            5,
            8,
            x^17 * y^3 + x^16 + x^6 + x^7 * y^2,
            x^10 * y + y + x^11 + x * y,
            "4",
        ),
        (
            "w8_inc_013",
            "tab:increasing-k-codes-wt-8",
            "weight 8 increasing-k cover sequence",
            "[[160,24,8]]",
            20,
            4,
            5,
            8,
            x * y^3 + x^16 + x^14 + x^19 * y^2,
            x^6 * y + x^8 * y + x^3 + x * y,
            "9.6",
        ),
        (
            "w8_inc_014",
            "tab:increasing-k-codes-wt-8",
            "weight 8 increasing-k cover sequence",
            "[[160,16,d<=14]]",
            20,
            4,
            5,
            8,
            x^17 * y^3 + x^12 + x^10 + x^15 * y^2,
            x^6 * y + x^12 * y + x^7 + x * y,
            "<= 19.6",
        ),
        (
            "w8_inc_015",
            "tab:increasing-k-codes-wt-8",
            "weight 8 increasing-k cover sequence",
            "[[96,24,4]]",
            12,
            4,
            3,
            8,
            x * y^3 + 1 + x^6 + x^7 * y^2,
            x^2 * y + x^8 * y + x^3 + x^9 * y,
            "4",
        ),
        (
            "w8_inc_016",
            "tab:increasing-k-codes-wt-8",
            "weight 8 increasing-k cover sequence",
            "[[96,20,8]]",
            12,
            4,
            3,
            8,
            x * y^3 + x^4 + x^2 + x^11 * y^2,
            x^6 * y + x^8 * y + x^11 + x^9 * y,
            "13.3",
        ),
        (
            "w8_inc_017",
            "tab:increasing-k-codes-wt-8",
            "weight 8 increasing-k cover sequence",
            "[[96,16,4]]",
            12,
            4,
            3,
            8,
            x^5 * y^3 + x^4 + x^10 + x^11 * y^2,
            x^6 * y + y + x^3 + x^9 * y,
            "2.7",
        ),
        (
            "w8_inc_018",
            "tab:increasing-k-codes-wt-8",
            "weight 8 increasing-k cover sequence",
            "[[96,12,8]]",
            12,
            4,
            3,
            8,
            x^9 * y^3 + 1 + x^6 + x^3 * y^2,
            x^10 * y + x^8 * y + x^7 + x^5 * y,
            "8",
        ),
        # [[64,14,8]] duplicates w8_k14_alt_004.
        (
            "w8_inc_020",
            "tab:increasing-k-codes-wt-8",
            "weight 8 increasing-k cover sequence",
            "[[32,8,4]]",
            4,
            4,
            1,
            8,
            x * y^3 + 1 + x^2 + x^3 * y^2,
            x^2 * y + y + x^3 + x * y,
            "4",
        ),
    ]
    return [_cover_bb_case(row...) for row in rows]
end

function double_gross_case()
    metadata = Dict{String, Any}(
        "source" => "Bravyi et al., Nature 627, 778-782 (2024)",
        "doi" => "10.1038/s41586-024-07107-7",
        "arxiv" => "2308.07915",
        "source_code" => "[[144,12,12]] Gross code",
        "parameters" => "[[288,12,18]]",
        "name" => "double_gross",
    )
    return _table_case(
        "double_gross",
        "Doubled Gross code [[288,12,18]] (Bravyi et al., 2024)",
        x^3 + y^2 + y^7,
        y^3 + x + x^2;
        metadata=metadata,
    )
end

function all_additional_bb_code_cases()
    return vcat(bunny_code_cases(), sbb_code_cases(), cover_bb_code_cases(), [double_gross_case()])
end

function _decoupled_cache_satisfies_request(decoupled_case::DecoupledToricCase; compute_inverse::Bool, capture_debug::Bool)
    capture_debug && isnothing(decoupled_case.debug_result) && return false
    !compute_inverse && return true

    result = decoupled_case.transfer_result
    return !isnothing(result) &&
           hasproperty(result, :phi_1_inv) &&
           !isnothing(result.phi_1_inv) &&
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
    cases=all_additional_bb_code_cases();
    results_path=joinpath(dirname(dirname(dirname(dirname(@__DIR__)))), "build", "reproduction", "additional_bb_code_decoupling_results.md"),
    cache_dir=joinpath(dirname(dirname(dirname(dirname(@__DIR__)))), "build", "reproduction", "additional_bb_code_decoupling_cache_v2"),
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

    first_case = iterate(cases)
    if warmup && !isnothing(first_case)
        case, _ = first_case
        _warmup_decoupling(
            case.poly_vector;
            show_progress=show_progress,
            compute_inverse=compute_inverse,
            max_area=max_area,
            max_l=max_l,
        )
    end

    open(results_path, "w") do io
        println(io, "| Code | Definition / decouple input | \$L\$ | \$u\$ | \$v\$ | Area | Matrix size after CG | Product state num | Toric num | Time(s) | Maximum_term_phi_1 | Maximum_degree_phi_1 | Maximum_term_phi_1_inv | Maximum_degree_phi_1_inv | Maximum_column_monomial_count_phi_1 | Maximum_column_monomial_count_phi_1_inv | Case ID |")
        println(io, "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")

        for case in cases
            cache_path = decoupled_toric_case_path(cache_dir, case.case_id)
            reused_row = get(reused_rows, case.case_id, nothing)
            if !overwrite && !check_cache && isfile(cache_path) && !isnothing(reused_row)
                println(io, reused_row)
                flush(io)
                continue
            end

            decoupled_case = _load_or_build_decoupled_case(
                case.case_id,
                case.poly_vector,
                cache_path;
                metadata=_case_metadata(case),
                overwrite=overwrite,
                show_progress=show_progress,
                compute_inverse=compute_inverse,
                capture_debug=capture_debug,
                check_cache=check_cache,
                max_area=max_area,
                max_l=max_l,
            )
            println(io, _markdown_row_from_decoupled_case(case, decoupled_case))
            flush(io)
        end
    end
end

function main(; cases=all_additional_bb_code_cases(), kwargs...)
    return run_and_save(cases; kwargs...)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
