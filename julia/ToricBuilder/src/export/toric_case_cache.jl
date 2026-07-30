const DECOUPLED_TORIC_CASE_FORMAT_VERSION = 3
const REQUIRED_DECOUPLED_TRANSFER_RESULT_FIELDS = (
    :input_matrix,
    :input_blocks,
    :standard_matrix,
    :standard_blocks,
    :row_transformation,
    :row_blocks,
    :psi_1_inverse,
)

const V2_TRANSFER_RESULT_FIELD_RENAMES = Dict(
    :phi_1 => :psi_1_inverse,
    :phi_1_inv => :psi_1,
    :phi_1_size => :psi_1_inverse_size,
    :phi_1_original => :psi_1_inverse_original,
    :phi_1_inv_original => :psi_1_original,
    :max_ele_phi_1 => :max_ele_psi_1_inverse,
    :max_degree_phi_1 => :max_degree_psi_1_inverse,
    :max_column_monomial_count_phi_1 => :max_column_monomial_count_psi_1_inverse,
    :max_ele_phi_1_inv => :max_ele_psi_1,
    :max_degree_phi_1_inv => :max_degree_psi_1,
    :max_column_monomial_count_phi_1_inv => :max_column_monomial_count_psi_1,
)

const V2_DEBUG_RESULT_FIELD_RENAMES = Dict(
    :phi_1 => :psi_1_inverse,
    :phi_1_inv => :psi_1,
)

const REQUIRED_V1_FULL_TRANSFER_RESULT_FIELDS = (
    :mat_after_coarse_graining,
    :result_matrix,
    :row_transformation,
    :column_transformation,
)

struct DecoupledToricCase
    format_version::Int
    case_id::String
    status::Symbol
    metadata::Dict{String, Any}
    poly_vec
    transfer_result
    debug_result
    created_at::String
    runtime_info::Dict{String, Any}
end

decoupled_toric_case_path(dir::AbstractString, case_id::AbstractString) =
    joinpath(dir, string(case_id, ".jls"))

function build_decoupled_toric_case(
    case_id,
    poly_vec;
    metadata=Dict(),
    show_progress::Bool=false,
    capture_exceptions::Bool=false,
    compute_inverse::Bool=false,
    capture_debug::Bool=false,
    kwargs...,
)
    normalized_metadata = _normalize_cached_toric_case_metadata(metadata)
    poly_vec_snapshot = copy(poly_vec)
    transfer_result = nothing
    debug_result = nothing

    try
        transfer_result = build_toric_form(
            poly_vec_snapshot;
            show_progress=show_progress,
            compute_inverse=compute_inverse,
            kwargs...,
        )
    catch err
        _rethrow_unless_capturable(err, capture_exceptions)
        return _failed_decoupled_toric_case(
            case_id,
            poly_vec_snapshot,
            normalized_metadata,
            nothing,
            "build_toric_form failed: $(sprint(showerror, err))";
            exception_type=string(nameof(typeof(err))),
        )
    end

    missing_fields = _missing_decoupled_transfer_result_fields(transfer_result)
    if !isempty(missing_fields)
        return _failed_decoupled_toric_case(
            case_id,
            poly_vec_snapshot,
            normalized_metadata,
            transfer_result,
            "build_toric_form did not produce reusable decoupled matrices; missing fields: $(join(missing_fields, ", "))",
        )
    end

    if capture_debug
        try
            debug_result = capture_toric_form_debug_matrices(
                transfer_result.input_matrix;
                show_progress=show_progress,
                compute_inverse=compute_inverse,
            )
        catch err
            _rethrow_unless_capturable(err, capture_exceptions)
            return _failed_decoupled_toric_case(
                case_id,
                poly_vec_snapshot,
                normalized_metadata,
                transfer_result,
                "capture_toric_form_debug_matrices failed: $(sprint(showerror, err))";
                exception_type=string(nameof(typeof(err))),
            )
        end
    end

    return DecoupledToricCase(
        DECOUPLED_TORIC_CASE_FORMAT_VERSION,
        string(case_id),
        :ok,
        normalized_metadata,
        poly_vec_snapshot,
        transfer_result,
        debug_result,
        _cached_toric_case_timestamp(),
        _cached_toric_case_runtime_info(),
    )
end

function save_decoupled_toric_case(path::AbstractString, decoupled_case::DecoupledToricCase; overwrite::Bool=false)
    if ispath(path) && !overwrite
        throw(ArgumentError("Refusing to overwrite existing decoupled cache file: $path"))
    end
    dir = dirname(path)
    mkpath(dir)
    tmp_path, io = mktemp(dir; cleanup=false)
    close(io)
    try
        Oscar.save(tmp_path, _decoupled_toric_case_payload(decoupled_case))
        if overwrite && ispath(path)
            rm(path; force=true)
        end
        Base.Filesystem.rename(tmp_path, path)
    catch
        ispath(tmp_path) && rm(tmp_path; force=true)
        rethrow()
    end
    return path
end

function load_decoupled_toric_case(path::AbstractString)
    payload = Oscar.load(path)
    return _migrate_decoupled_toric_case_payload(payload)
end

function list_decoupled_toric_cases(dir::AbstractString)
    if !isdir(dir)
        return String[]
    end
    paths = String[]
    for entry in readdir(dir)
        path = joinpath(dir, entry)
        if isfile(path) && endswith(entry, ".jls")
            push!(paths, path)
        end
    end
    sort!(paths)
    return paths
end

function migrate_cached_toric_cases_to_decoupled(
    src_dir::AbstractString,
    dst_dir::AbstractString;
    overwrite::Bool=false,
)
    paths = list_decoupled_toric_cases(src_dir)
    mkpath(dst_dir)
    written = String[]
    for src_path in paths
        payload = Oscar.load(src_path)
        decoupled_case = _migrate_decoupled_toric_case_payload(payload)
        dst_path = decoupled_toric_case_path(dst_dir, decoupled_case.case_id)
        save_decoupled_toric_case(dst_path, decoupled_case; overwrite=overwrite)
        push!(written, dst_path)
    end
    return written
end

function _migrate_decoupled_toric_case_payload(payload)
    version = Int(_payload_field(payload, :format_version))
    version == 1 && return _migrate_v1_payload_to_decoupled(payload)
    version == 2 && return _migrate_v2_payload_to_decoupled(payload)
    version == 3 && return _decoupled_toric_case_from_payload(payload)
    throw(ArgumentError("Unsupported DecoupledToricCase format version $version"))
end

function _migrate_v1_payload_to_decoupled(payload)
    version = Int(_payload_field(payload, :format_version))
    version == 1 || throw(ArgumentError("Expected v1 CachedToricCase payload, got format version $version"))
    transfer_result = _restore_serialization_safe_value(_payload_field_or(payload, :transfer_result, nothing))
    migrated_transfer_result = if isnothing(transfer_result)
        nothing
    elseif _has_v1_full_transfer_result(transfer_result)
        _migrate_v1_transfer_result(transfer_result)
    else
        _migrate_v1_partial_transfer_result(transfer_result)
    end
    debug_result = _restore_serialization_safe_value(_payload_field_or(payload, :debug_result, nothing))
    migrated_debug_result = isnothing(debug_result) ? nothing : _migrate_v1_debug_result(debug_result)
    status = _payload_field(payload, :status)
    return DecoupledToricCase(
        DECOUPLED_TORIC_CASE_FORMAT_VERSION,
        string(_payload_field(payload, :case_id)),
        status isa Symbol ? status : Symbol(status),
        _normalize_cached_toric_case_metadata(_payload_field(payload, :metadata)),
        _payload_field(payload, :poly_vec),
        migrated_transfer_result,
        migrated_debug_result,
        string(_payload_field(payload, :created_at)),
        _normalize_cached_toric_case_metadata(_payload_field(payload, :runtime_info)),
    )
end

function _migrate_v2_payload_to_decoupled(payload)
    version = Int(_payload_field(payload, :format_version))
    version == 2 || throw(ArgumentError("Expected v2 DecoupledToricCase payload, got format version $version"))
    transfer_result = _restore_serialization_safe_value(_payload_field_or(payload, :transfer_result, nothing))
    debug_result = _restore_serialization_safe_value(_payload_field_or(payload, :debug_result, nothing))
    canonical_payload = Dict{String, Any}(
        "format_version" => DECOUPLED_TORIC_CASE_FORMAT_VERSION,
        "case_id" => _payload_field(payload, :case_id),
        "status" => _payload_field(payload, :status),
        "metadata" => _payload_field(payload, :metadata),
        "poly_vec" => _payload_field(payload, :poly_vec),
        "created_at" => _payload_field(payload, :created_at),
        "runtime_info" => _payload_field(payload, :runtime_info),
    )
    if !isnothing(transfer_result)
        canonical_payload["transfer_result"] = _serialization_safe_value(
            _rename_namedtuple_fields(transfer_result, V2_TRANSFER_RESULT_FIELD_RENAMES),
        )
    end
    if !isnothing(debug_result)
        canonical_payload["debug_result"] = _serialization_safe_value(
            _rename_namedtuple_fields(debug_result, V2_DEBUG_RESULT_FIELD_RENAMES),
        )
    end
    return _decoupled_toric_case_from_payload(canonical_payload)
end

function _has_v1_full_transfer_result(old)
    for field in REQUIRED_V1_FULL_TRANSFER_RESULT_FIELDS
        if !hasproperty(old, field) || isnothing(getproperty(old, field))
            return false
        end
    end
    return true
end

function _migrate_v1_transfer_result(old)
    input_matrix = old.mat_after_coarse_graining
    standard_matrix = old.result_matrix
    row_transformation = old.row_transformation
    column_transformation = old.column_transformation
    stab_num = size(input_matrix, 1) ÷ 2
    qubit_num = size(input_matrix, 2) ÷ 2
    psi_1_inverse = column_transformation[1:qubit_num, 1:qubit_num]
    psi_1 = _dagger_laurent_matrix(column_transformation[qubit_num+1:2*qubit_num, qubit_num+1:2*qubit_num])
    return (;
        input_matrix=input_matrix,
        input_blocks=(
            Hz=input_matrix[1:stab_num, 1:qubit_num],
            Hx=input_matrix[stab_num+1:2*stab_num, qubit_num+1:2*qubit_num],
        ),
        standard_matrix=standard_matrix,
        standard_blocks=(
            Hz=standard_matrix[1:stab_num, 1:qubit_num],
            Hx=standard_matrix[stab_num+1:2*stab_num, qubit_num+1:2*qubit_num],
        ),
        row_transformation=row_transformation,
        row_blocks=(
            Hz=row_transformation[1:stab_num, 1:stab_num],
            Hx=row_transformation[stab_num+1:2*stab_num, stab_num+1:2*stab_num],
        ),
        psi_1_inverse=psi_1_inverse,
        psi_1=psi_1,
        column_transformation=column_transformation,
        product_state_num=old.product_state_num,
        toric_num=old.toric_num,
        l=get(old, :l, nothing),
        area=get(old, :area, nothing),
        u_rel=get(old, :u_rel, nothing),
        v_rel=get(old, :v_rel, nothing),
        solving_time=get(old, :solving_time, nothing),
        A_size=get(old, :A_size, size(input_matrix)),
        psi_1_inverse_size=size(psi_1_inverse),
        max_ele_psi_1_inverse=get(old, :max_eleQ, nothing),
        max_degree_psi_1_inverse=get(old, :max_degreeQ, nothing),
        max_column_monomial_count_psi_1_inverse=get(old, :max_column_monomial_countQ, nothing),
        max_ele_psi_1=get(old, :max_eleQinv, nothing),
        max_degree_psi_1=get(old, :max_degreeQinv, nothing),
        max_column_monomial_count_psi_1=get(old, :max_column_monomial_countQinv, nothing),
    )
end

function _migrate_v1_partial_transfer_result(old)
    return (;
        input_matrix=nothing,
        input_blocks=nothing,
        standard_matrix=nothing,
        standard_blocks=nothing,
        row_transformation=nothing,
        row_blocks=nothing,
        psi_1_inverse=nothing,
        psi_1=nothing,
        column_transformation=nothing,
        product_state_num=get(old, :product_state_num, nothing),
        toric_num=get(old, :toric_num, nothing),
        l=get(old, :l, nothing),
        area=get(old, :area, nothing),
        u_rel=get(old, :u_rel, nothing),
        v_rel=get(old, :v_rel, nothing),
        solving_time=get(old, :solving_time, nothing),
        A_size=get(old, :A_size, nothing),
        psi_1_inverse_size=get(old, :Q_size, nothing),
        max_ele_psi_1_inverse=get(old, :max_eleQ, nothing),
        max_degree_psi_1_inverse=get(old, :max_degreeQ, nothing),
        max_column_monomial_count_psi_1_inverse=get(old, :max_column_monomial_countQ, nothing),
        max_ele_psi_1=get(old, :max_eleQinv, nothing),
        max_degree_psi_1=get(old, :max_degreeQinv, nothing),
        max_column_monomial_count_psi_1=get(old, :max_column_monomial_countQinv, nothing),
    )
end

function _migrate_v1_debug_result(old)
    input_matrix = old.input_matrix
    standard_matrix = old.result_matrix
    row_transformation = old.row_transformation
    column_transformation = old.column_transformation
    stab_num = size(input_matrix, 1) ÷ 2
    qubit_num = size(input_matrix, 2) ÷ 2
    psi_1_inverse = isnothing(column_transformation) ? nothing : column_transformation[1:qubit_num, 1:qubit_num]
    psi_1 = isnothing(column_transformation) ? nothing : _dagger_laurent_matrix(column_transformation[qubit_num+1:2*qubit_num, qubit_num+1:2*qubit_num])
    return (;
        input_matrix=input_matrix,
        em_matrix=old.em_matrix,
        row_transformation=row_transformation,
        row_blocks=(
            Hz=row_transformation[1:stab_num, 1:stab_num],
            Hx=row_transformation[stab_num+1:2*stab_num, stab_num+1:2*stab_num],
        ),
        psi_1_inverse=psi_1_inverse,
        psi_1=psi_1,
        column_transformation=column_transformation,
        standard_matrix=standard_matrix,
        standard_blocks=(
            Hz=standard_matrix[1:stab_num, 1:qubit_num],
            Hx=standard_matrix[stab_num+1:2*stab_num, qubit_num+1:2*qubit_num],
        ),
        product_state_num=old.product_state_num,
        toric_num=old.toric_num,
    )
end

function _failed_decoupled_toric_case(
    case_id,
    poly_vec_snapshot,
    metadata,
    transfer_result,
    reason::AbstractString;
    exception_type::Union{Nothing, String}=nothing,
)
    failed_metadata = copy(metadata)
    failed_metadata["failure_reason"] = reason
    if !isnothing(exception_type)
        failed_metadata["exception_type"] = exception_type
    end
    return DecoupledToricCase(
        DECOUPLED_TORIC_CASE_FORMAT_VERSION,
        string(case_id),
        :failed,
        failed_metadata,
        poly_vec_snapshot,
        transfer_result,
        nothing,
        _cached_toric_case_timestamp(),
        _cached_toric_case_runtime_info(),
    )
end

function _decoupled_toric_case_payload(decoupled_case::DecoupledToricCase)
    payload = Dict{String, Any}(
        "format_version" => DECOUPLED_TORIC_CASE_FORMAT_VERSION,
        "case_id" => decoupled_case.case_id,
        "status" => decoupled_case.status,
        "metadata" => _serialization_safe_value(decoupled_case.metadata),
        "poly_vec" => decoupled_case.poly_vec,
        "created_at" => decoupled_case.created_at,
        "runtime_info" => _serialization_safe_value(decoupled_case.runtime_info),
    )
    if !isnothing(decoupled_case.transfer_result)
        payload["transfer_result"] = _serialization_safe_value(decoupled_case.transfer_result)
    end
    if !isnothing(decoupled_case.debug_result)
        payload["debug_result"] = _serialization_safe_value(decoupled_case.debug_result)
    end
    return payload
end

function _serialization_safe_value(value)
    if isnothing(value)
        return nothing
    elseif value isa AbstractDict
        sanitized = _empty_serialization_dict(value)
        for (key, entry) in pairs(value)
            sanitized_key = _serialization_safe_value(key)
            sanitized_entry = _serialization_safe_value(entry)
            if !isnothing(sanitized_key) && !isnothing(sanitized_entry)
                sanitized[sanitized_key] = sanitized_entry
            end
        end
        return sanitized
    elseif value isa NamedTuple
        return _serialization_safe_named_tuple(value)
    elseif value isa Tuple
        return tuple((_serialization_safe_value(entry) for entry in value)...)
    elseif value isa AbstractArray
        return map(_serialization_safe_value, value)
    end
    return value
end

function _serialization_safe_named_tuple(value::NamedTuple)
    names = propertynames(value)
    values = Any[]
    fields = Dict{String, Any}()
    omitted_any = false

    for name in names
        sanitized_entry = _serialization_safe_value(getproperty(value, name))
        if !isnothing(sanitized_entry)
            push!(values, sanitized_entry)
            fields[string(name)] = sanitized_entry
        else
            omitted_any = true
        end
    end

    if !omitted_any
        return NamedTuple{names}(Tuple(values))
    end

    return Dict{String, Any}(
        "__cache_type__" => "namedtuple",
        "__field_order__" => string.(names),
        "__fields__" => fields,
    )
end

function _rename_namedtuple_fields(value::NamedTuple, renames::AbstractDict{Symbol, Symbol})
    names = Symbol[]
    values = Any[]
    seen = Set{Symbol}()
    for name in propertynames(value)
        renamed = get(renames, name, name)
        renamed in seen && throw(ArgumentError(
            "Cannot rename cache field $name to $renamed because that field already exists",
        ))
        push!(seen, renamed)
        push!(names, renamed)
        push!(values, getproperty(value, name))
    end
    return NamedTuple{Tuple(names)}(Tuple(values))
end

function _decoupled_toric_case_from_payload(payload)
    required_fields = (:format_version, :case_id, :status, :metadata, :poly_vec, :created_at, :runtime_info)
    missing_fields = [string(field) for field in required_fields if !_payload_has_field(payload, field)]
    if !isempty(missing_fields)
        throw(ArgumentError("Cache payload is missing required fields: $(join(missing_fields, ", "))"))
    end
    metadata = _normalize_cached_toric_case_metadata(_payload_field(payload, :metadata))
    runtime_info = _normalize_cached_toric_case_metadata(_payload_field(payload, :runtime_info))
    status = _payload_field(payload, :status)
    normalized_status = status isa Symbol ? status : Symbol(status)
    transfer_result = _restore_serialization_safe_value(_payload_field_or(payload, :transfer_result, nothing))
    debug_result = _restore_serialization_safe_value(_payload_field_or(payload, :debug_result, nothing))
    return DecoupledToricCase(
        Int(_payload_field(payload, :format_version)),
        string(_payload_field(payload, :case_id)),
        normalized_status,
        metadata,
        _payload_field(payload, :poly_vec),
        transfer_result,
        debug_result,
        string(_payload_field(payload, :created_at)),
        runtime_info,
    )
end

function _restore_serialization_safe_value(value)
    if isnothing(value)
        return nothing
    elseif value isa AbstractDict && _is_namedtuple_wrapper(value)
        field_order = _payload_field(value, Symbol("__field_order__"))
        fields = _payload_field(value, Symbol("__fields__"))
        names = Symbol[]
        restored_values = Any[]
        for field_name in field_order
            symbol_name = Symbol(field_name)
            push!(names, symbol_name)
            field_value = _payload_field_or(fields, symbol_name, nothing)
            push!(restored_values, _restore_serialization_safe_value(field_value))
        end
        return NamedTuple{Tuple(names)}(Tuple(restored_values))
    elseif value isa AbstractDict
        restored = _empty_serialization_dict(value)
        for (key, entry) in pairs(value)
            restored[_restore_serialization_safe_value(key)] = _restore_serialization_safe_value(entry)
        end
        return restored
    elseif value isa NamedTuple
        names = propertynames(value)
        restored_values = map(name -> _restore_serialization_safe_value(getproperty(value, name)), names)
        return NamedTuple{names}(Tuple(restored_values))
    elseif value isa Tuple
        return tuple((_restore_serialization_safe_value(entry) for entry in value)...)
    elseif value isa AbstractArray
        return map(_restore_serialization_safe_value, value)
    end
    return value
end

function _is_namedtuple_wrapper(value::AbstractDict)
    return _payload_field_or(value, Symbol("__cache_type__"), nothing) == "namedtuple"
end

function _empty_serialization_dict(value::AbstractDict)
    if keytype(value) === String
        return Dict{String, Any}()
    elseif keytype(value) === Symbol
        return Dict{Symbol, Any}()
    end
    return Dict{Any, Any}()
end

function _normalize_cached_toric_case_metadata(metadata)
    normalized = Dict{String, Any}()
    if isnothing(metadata)
        return normalized
    end

    for (key, value) in pairs(metadata)
        normalized[string(key)] = value
    end
    return normalized
end

function _missing_decoupled_transfer_result_fields(transfer_result)
    missing = String[]
    if isnothing(transfer_result)
        return [string(field) for field in REQUIRED_DECOUPLED_TRANSFER_RESULT_FIELDS]
    end
    for field in REQUIRED_DECOUPLED_TRANSFER_RESULT_FIELDS
        if !hasproperty(transfer_result, field) || isnothing(getproperty(transfer_result, field))
            push!(missing, string(field))
        end
    end
    return missing
end

function _rethrow_unless_capturable(err, capture_exceptions::Bool)
    if !capture_exceptions || err isa InterruptException
        rethrow()
    end
    return nothing
end

_cached_toric_case_timestamp() = string(round(Int, time()))

function _cached_toric_case_runtime_info()
    return Dict{String, Any}(
        "julia_version" => string(VERSION),
    )
end

_payload_has_field(payload::NamedTuple, field::Symbol) = hasproperty(payload, field)
_payload_has_field(payload::AbstractDict, field::Symbol) = haskey(payload, field) || haskey(payload, string(field))
_payload_has_field(payload, field::Symbol) = false

_payload_field(payload::NamedTuple, field::Symbol) = getproperty(payload, field)
function _payload_field(payload::AbstractDict, field::Symbol)
    if haskey(payload, field)
        return payload[field]
    elseif haskey(payload, string(field))
        return payload[string(field)]
    end
    throw(KeyError(field))
end

function _payload_field_or(payload, field::Symbol, default)
    if _payload_has_field(payload, field)
        return _payload_field(payload, field)
    end
    return default
end
