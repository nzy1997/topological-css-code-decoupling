using Test
using ToricBuilder
using ToricBuilder:
    DECOUPLED_TORIC_CASE_FORMAT_VERSION,
    build_decoupled_toric_case,
    decoupled_toric_case_path,
    load_decoupled_toric_case,
    migrate_cached_toric_cases_to_decoupled,
    save_decoupled_toric_case
using Oscar

@testset "DecoupledToricCase round-trip default" begin
    F = GF(2)
    Rxy, (x, y) = laurent_polynomial_ring(F, ["x", "y"])
    poly_vec = [1 + x + x*y, 1 + y + x*y]
    original_first = poly_vec[1]

    mktempdir() do tmpdir
        path = joinpath(tmpdir, "xy_case.jls")
        decoupled = build_decoupled_toric_case("xy_case", poly_vec; show_progress=false)
        @test decoupled.status == :ok
        @test decoupled.format_version == DECOUPLED_TORIC_CASE_FORMAT_VERSION
        @test decoupled.debug_result === nothing
        @test DECOUPLED_TORIC_CASE_FORMAT_VERSION == 3
        @test hasproperty(decoupled.transfer_result, :psi_1_inverse)
        @test decoupled.transfer_result.psi_1 === nothing
        @test !hasproperty(decoupled.transfer_result, :phi_1)
        @test !hasproperty(decoupled.transfer_result, :phi_1_inv)
        @test decoupled.transfer_result.column_transformation === nothing
        @test hasproperty(decoupled.transfer_result, :input_blocks)
        @test hasproperty(decoupled.transfer_result, :standard_blocks)
        poly_vec[1] = 1 + x
        @test decoupled.poly_vec[1] == original_first

        save_decoupled_toric_case(path, decoupled)
        raw_payload = Oscar.load(path)
        raw_transfer = ToricBuilder._restore_serialization_safe_value(
            ToricBuilder._payload_field(raw_payload, :transfer_result),
        )
        @test hasproperty(raw_transfer, :psi_1_inverse)
        @test !hasproperty(raw_transfer, :phi_1)
        @test !hasproperty(raw_transfer, :phi_1_inv)
        loaded = load_decoupled_toric_case(path)

        @test loaded.case_id == "xy_case"
        @test loaded.status == :ok
        @test loaded.debug_result === nothing
        @test loaded.transfer_result.input_matrix == decoupled.transfer_result.input_matrix
        @test loaded.transfer_result.standard_matrix == decoupled.transfer_result.standard_matrix
        @test loaded.transfer_result.input_blocks.Hx == decoupled.transfer_result.input_blocks.Hx
        @test loaded.transfer_result.input_blocks.Hz == decoupled.transfer_result.input_blocks.Hz
        @test loaded.transfer_result.standard_blocks.Hz == decoupled.transfer_result.standard_blocks.Hz
        @test loaded.transfer_result.standard_blocks.Hx == decoupled.transfer_result.standard_blocks.Hx
    end
end

@testset "DecoupledToricCase debug capture is opt-in" begin
    F = GF(2)
    Rxy, (x, y) = laurent_polynomial_ring(F, ["x", "y"])
    poly_vec = [1 + x + x*y, 1 + y + x*y]

    decoupled = build_decoupled_toric_case(
        "xy_debug",
        poly_vec;
        show_progress=false,
        capture_debug=true,
    )
    @test decoupled.status == :ok
    @test !isnothing(decoupled.debug_result)
    @test decoupled.debug_result.input_matrix == decoupled.transfer_result.input_matrix
    @test decoupled.debug_result.standard_matrix == decoupled.transfer_result.standard_matrix
    @test decoupled.debug_result.psi_1_inverse == decoupled.transfer_result.psi_1_inverse
end

@testset "DecoupledToricCase inverse mode round-trip" begin
    F = GF(2)
    Rxy, (x, y) = laurent_polynomial_ring(F, ["x", "y"])
    poly_vec = [1 + x + x*y, 1 + y + x*y]

    mktempdir() do tmpdir
        path = joinpath(tmpdir, "xy_inverse.jls")
        decoupled = build_decoupled_toric_case("xy_inverse", poly_vec; show_progress=false, compute_inverse=true)
        @test !isnothing(decoupled.transfer_result.psi_1)
        @test !isnothing(decoupled.transfer_result.column_transformation)
        save_decoupled_toric_case(path, decoupled)
        loaded = load_decoupled_toric_case(path)
        loaded_ring = base_ring(loaded.transfer_result.psi_1_inverse)
        @test loaded.transfer_result.psi_1_inverse * loaded.transfer_result.psi_1 ==
              identity_matrix(loaded_ring, size(loaded.transfer_result.psi_1_inverse, 1))
    end
end

@testset "load_decoupled_toric_case rejects unsupported version" begin
    mktempdir() do tmpdir
        path = joinpath(tmpdir, "legacy_case.jls")
        payload = Dict{String, Any}(
            "format_version" => 4,
            "case_id" => "legacy_case",
            "status" => :ok,
            "metadata" => Dict{String, String}(),
            "poly_vec" => String[],
            "created_at" => "0",
            "runtime_info" => Dict{String, String}(),
        )
        Oscar.save(path, payload)

        err = try
            load_decoupled_toric_case(path)
            nothing
        catch caught
            caught
        end
        @test err isa ArgumentError
        @test occursin("Unsupported DecoupledToricCase format version 4", sprint(showerror, err))
    end
end

@testset "v1 cache migration extracts canonical psi maps" begin
    F = GF(2)
    Rxy, (x, y) = laurent_polynomial_ring(F, ["x", "y"])
    old = build_toric_form([1 + x + x*y, 1 + y + x*y]; show_progress=false, compute_inverse=true)
    q = size(old.column_transformation, 1) ÷ 2
    expected_psi_1_inverse = old.column_transformation[1:q, 1:q]
    expected_psi_1 = ToricBuilder._dagger_laurent_matrix(old.column_transformation[q+1:2*q, q+1:2*q])

    mktempdir() do tmpdir
        src_dir = joinpath(tmpdir, "v1")
        dst_dir = joinpath(tmpdir, "v3")
        mkpath(src_dir)
        v1_path = joinpath(src_dir, "case_001.jls")
        Oscar.save(v1_path, Dict{String, Any}(
            "format_version" => 1,
            "case_id" => "case_001",
            "status" => :ok,
            "metadata" => Dict{String, Any}("source" => "synthetic"),
            "poly_vec" => [1 + x + x*y, 1 + y + x*y],
            "transfer_result" => (
                mat_after_coarse_graining=old.input_matrix,
                result_matrix=old.standard_matrix,
                row_transformation=old.row_transformation,
                column_transformation=old.column_transformation,
                product_state_num=old.product_state_num,
                toric_num=old.toric_num,
                l=old.l,
                area=old.area,
                u_rel=old.u_rel,
                v_rel=old.v_rel,
                solving_time=old.solving_time,
                A_size=old.A_size,
                Q_size=size(old.column_transformation),
                max_eleQ=old.max_ele_psi_1_inverse,
                max_degreeQ=old.max_degree_psi_1_inverse,
                max_eleQinv=old.max_ele_psi_1,
                max_degreeQinv=old.max_degree_psi_1,
                max_column_monomial_countQ=old.max_column_monomial_count_psi_1_inverse,
                max_column_monomial_countQinv=old.max_column_monomial_count_psi_1,
            ),
            "created_at" => "0",
            "runtime_info" => Dict{String, Any}("julia_version" => string(VERSION)),
        ))

        written = migrate_cached_toric_cases_to_decoupled(src_dir, dst_dir)
        @test length(written) == 1
        migrated = load_decoupled_toric_case(joinpath(dst_dir, "case_001.jls"))
        @test migrated.format_version == 3
        @test migrated.transfer_result.psi_1_inverse == expected_psi_1_inverse
        @test migrated.transfer_result.psi_1 == expected_psi_1
        @test migrated.transfer_result.column_transformation == old.column_transformation
        @test migrated.transfer_result.input_blocks.Hz == old.input_blocks.Hz
        @test migrated.transfer_result.standard_blocks.Hx == old.standard_blocks.Hx
        @test migrated.transfer_result.max_ele_psi_1_inverse == old.max_ele_psi_1_inverse
        @test migrated.transfer_result.max_ele_psi_1 == old.max_ele_psi_1
    end
end

@testset "load_decoupled_toric_case migrates v2 fields in memory" begin
    F = GF(2)
    Rxy, (x, y) = laurent_polynomial_ring(F, ["x", "y"])
    current = build_decoupled_toric_case(
        "v2_case",
        [1 + x + x*y, 1 + y + x*y];
        show_progress=false,
        compute_inverse=true,
        capture_debug=true,
    )
    v3_fields = (
        :psi_1_inverse,
        :psi_1,
        :psi_1_inverse_size,
        :psi_1_inverse_original,
        :psi_1_original,
        :max_ele_psi_1_inverse,
        :max_degree_psi_1_inverse,
        :max_column_monomial_count_psi_1_inverse,
        :max_ele_psi_1,
        :max_degree_psi_1,
        :max_column_monomial_count_psi_1,
    )
    retained_fields = Tuple(name for name in propertynames(current.transfer_result) if !(name in v3_fields))
    retained_transfer = NamedTuple{retained_fields}(
        Tuple(getproperty(current.transfer_result, name) for name in retained_fields),
    )
    legacy_v2_transfer = (;
        retained_transfer...,
        phi_1=current.transfer_result.psi_1_inverse,
        phi_1_inv=current.transfer_result.psi_1,
        phi_1_size=current.transfer_result.psi_1_inverse_size,
        phi_1_original=current.transfer_result.psi_1_inverse_original,
        phi_1_inv_original=current.transfer_result.psi_1_original,
        max_ele_phi_1=current.transfer_result.max_ele_psi_1_inverse,
        max_degree_phi_1=current.transfer_result.max_degree_psi_1_inverse,
        max_column_monomial_count_phi_1=current.transfer_result.max_column_monomial_count_psi_1_inverse,
        max_ele_phi_1_inv=current.transfer_result.max_ele_psi_1,
        max_degree_phi_1_inv=current.transfer_result.max_degree_psi_1,
        max_column_monomial_count_phi_1_inv=current.transfer_result.max_column_monomial_count_psi_1,
    )
    debug_fields = (:psi_1_inverse, :psi_1)
    retained_debug_fields = Tuple(name for name in propertynames(current.debug_result) if !(name in debug_fields))
    retained_debug = NamedTuple{retained_debug_fields}(
        Tuple(getproperty(current.debug_result, name) for name in retained_debug_fields),
    )
    legacy_v2_debug = (;
        retained_debug...,
        phi_1=current.debug_result.psi_1_inverse,
        phi_1_inv=current.debug_result.psi_1,
    )

    mktempdir() do tmpdir
        path = joinpath(tmpdir, "v2_case.jls")
        payload = ToricBuilder._decoupled_toric_case_payload(current)
        payload["format_version"] = 2
        payload["transfer_result"] = ToricBuilder._serialization_safe_value(legacy_v2_transfer)
        payload["debug_result"] = ToricBuilder._serialization_safe_value(legacy_v2_debug)
        Oscar.save(path, payload)

        migrated = load_decoupled_toric_case(path)
        @test migrated.format_version == 3
        @test migrated.transfer_result.psi_1_inverse == current.transfer_result.psi_1_inverse
        @test migrated.transfer_result.psi_1 == current.transfer_result.psi_1
        @test migrated.transfer_result.psi_1_inverse_size == current.transfer_result.psi_1_inverse_size
        @test migrated.transfer_result.max_degree_psi_1_inverse == current.transfer_result.max_degree_psi_1_inverse
        @test migrated.transfer_result.max_degree_psi_1 == current.transfer_result.max_degree_psi_1
        @test !hasproperty(migrated.transfer_result, :phi_1)
        @test !hasproperty(migrated.transfer_result, :phi_1_inv)
        @test migrated.debug_result.psi_1_inverse == current.debug_result.psi_1_inverse
        @test !hasproperty(migrated.debug_result, :phi_1)
        @test !hasproperty(migrated.debug_result, :phi_1_inv)
        @test Int(ToricBuilder._payload_field(Oscar.load(path), :format_version)) == 2
    end
end

@testset "v1 failed partial cache migration preserves metrics" begin
    F = GF(2)
    Rxy, (x, y) = laurent_polynomial_ring(F, ["x", "y"])

    mktempdir() do tmpdir
        src_dir = joinpath(tmpdir, "v1")
        dst_dir = joinpath(tmpdir, "v3")
        mkpath(src_dir)
        Oscar.save(joinpath(src_dir, "case_failed.jls"), Dict{String, Any}(
            "format_version" => 1,
            "case_id" => "case_failed",
            "status" => :failed,
            "metadata" => Dict{String, Any}("source" => "partial legacy"),
            "poly_vec" => [1 + x, 1 + y],
            "transfer_result" => (
                solving_time=1.25,
                area=34,
                u_rel=x,
                v_rel=y,
            ),
            "created_at" => "0",
            "runtime_info" => Dict{String, Any}("julia_version" => string(VERSION)),
        ))

        written = migrate_cached_toric_cases_to_decoupled(src_dir, dst_dir)
        @test length(written) == 1
        migrated = load_decoupled_toric_case(joinpath(dst_dir, "case_failed.jls"))
        @test migrated.status == :failed
        @test migrated.transfer_result.solving_time == 1.25
        @test migrated.transfer_result.area == 34
        @test migrated.transfer_result.u_rel == x
        @test migrated.transfer_result.v_rel == y
        @test migrated.transfer_result.input_matrix === nothing
        @test migrated.transfer_result.psi_1_inverse === nothing
        @test migrated.transfer_result.psi_1 === nothing
    end
end

@testset "load_decoupled_toric_case migrates v1 payload in memory" begin
    mktempdir() do tmpdir
        path = joinpath(tmpdir, "old.jls")
        Oscar.save(path, Dict{String, Any}(
            "format_version" => 1,
            "case_id" => "old",
            "status" => :failed,
            "metadata" => Dict{String, Any}(),
            "poly_vec" => String[],
            "created_at" => "0",
            "runtime_info" => Dict{String, Any}(),
        ))
        migrated = load_decoupled_toric_case(path)
        @test migrated.format_version == 3
        @test migrated.case_id == "old"
    end
end

@testset "DecoupledToricCase failed status" begin
    F = GF(2)
    Rxy, (x, y) = laurent_polynomial_ring(F, ["x", "y"])
    poly_vec = [1 + x + x*y, 1 + y + x*y]

    decoupled = build_decoupled_toric_case("too_small_area", poly_vec; max_area=0, show_progress=false)
    @test decoupled.status == :failed
    @test decoupled.debug_result === nothing
    @test haskey(decoupled.metadata, "failure_reason")

    mktempdir() do tmpdir
        path = joinpath(tmpdir, "too_small_area.jls")
        save_decoupled_toric_case(path, decoupled)
        loaded = load_decoupled_toric_case(path)
        @test loaded.case_id == "too_small_area"
        @test loaded.status == :failed
        @test loaded.transfer_result.area == decoupled.transfer_result.area
    end
end

@testset "DecoupledToricCase overwrite protection and path" begin
    F = GF(2)
    Rxy, (x, y) = laurent_polynomial_ring(F, ["x", "y"])
    poly_vec = [1 + x + x*y, 1 + y + x*y]

    mktempdir() do tmpdir
        path = decoupled_toric_case_path(tmpdir, "xy_case")
        decoupled = build_decoupled_toric_case("xy_case", poly_vec; show_progress=false)
        save_decoupled_toric_case(path, decoupled)
        @test_throws ArgumentError save_decoupled_toric_case(path, decoupled)
        @test endswith(path, joinpath("xy_case.jls"))
    end
end

@testset "DecoupledToricCase exception capture is opt-in" begin
    bad_poly_vec = [1, 2]

    @test_throws MethodError build_decoupled_toric_case("bad_case", bad_poly_vec; show_progress=false)

    decoupled = build_decoupled_toric_case(
        "bad_case",
        bad_poly_vec;
        show_progress=false,
        capture_exceptions=true,
    )
    @test decoupled.status == :failed
    @test decoupled.debug_result === nothing
    @test decoupled.metadata["exception_type"] == "MethodError"
end
