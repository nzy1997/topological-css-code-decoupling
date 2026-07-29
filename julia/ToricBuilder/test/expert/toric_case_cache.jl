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
        @test decoupled.transfer_result.phi_1_inv === nothing
        @test decoupled.transfer_result.column_transformation === nothing
        @test hasproperty(decoupled.transfer_result, :input_blocks)
        @test hasproperty(decoupled.transfer_result, :standard_blocks)
        poly_vec[1] = 1 + x
        @test decoupled.poly_vec[1] == original_first

        save_decoupled_toric_case(path, decoupled)
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
    @test decoupled.debug_result.phi_1 == decoupled.transfer_result.phi_1
end

@testset "DecoupledToricCase inverse mode round-trip" begin
    F = GF(2)
    Rxy, (x, y) = laurent_polynomial_ring(F, ["x", "y"])
    poly_vec = [1 + x + x*y, 1 + y + x*y]

    mktempdir() do tmpdir
        path = joinpath(tmpdir, "xy_inverse.jls")
        decoupled = build_decoupled_toric_case("xy_inverse", poly_vec; show_progress=false, compute_inverse=true)
        @test !isnothing(decoupled.transfer_result.phi_1_inv)
        @test !isnothing(decoupled.transfer_result.column_transformation)
        save_decoupled_toric_case(path, decoupled)
        loaded = load_decoupled_toric_case(path)
        loaded_ring = base_ring(loaded.transfer_result.phi_1)
        @test loaded.transfer_result.phi_1 * loaded.transfer_result.phi_1_inv ==
              identity_matrix(loaded_ring, size(loaded.transfer_result.phi_1, 1))
    end
end

@testset "load_decoupled_toric_case rejects non-v2 payload with migration message" begin
    mktempdir() do tmpdir
        path = joinpath(tmpdir, "legacy_case.jls")
        payload = Dict{String, Any}(
            "format_version" => 1,
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
        @test occursin("migrate_cached_toric_cases_to_decoupled", sprint(showerror, err))
    end
end

@testset "v1 cache migration extracts phi_1 and phi_1_inv" begin
    F = GF(2)
    Rxy, (x, y) = laurent_polynomial_ring(F, ["x", "y"])
    old = build_toric_form([1 + x + x*y, 1 + y + x*y]; show_progress=false, compute_inverse=true)
    q = size(old.column_transformation, 1) ÷ 2
    expected_phi_1 = old.column_transformation[1:q, 1:q]
    expected_phi_1_inv = ToricBuilder._dagger_laurent_matrix(old.column_transformation[q+1:2*q, q+1:2*q])

    mktempdir() do tmpdir
        src_dir = joinpath(tmpdir, "v1")
        dst_dir = joinpath(tmpdir, "v2")
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
                max_eleQ=old.max_ele_phi_1,
                max_degreeQ=old.max_degree_phi_1,
                max_eleQinv=old.max_ele_phi_1_inv,
                max_degreeQinv=old.max_degree_phi_1_inv,
                max_column_monomial_countQ=old.max_column_monomial_count_phi_1,
                max_column_monomial_countQinv=old.max_column_monomial_count_phi_1_inv,
            ),
            "created_at" => "0",
            "runtime_info" => Dict{String, Any}("julia_version" => string(VERSION)),
        ))

        written = migrate_cached_toric_cases_to_decoupled(src_dir, dst_dir)
        @test length(written) == 1
        migrated = load_decoupled_toric_case(joinpath(dst_dir, "case_001.jls"))
        @test migrated.transfer_result.phi_1 == expected_phi_1
        @test migrated.transfer_result.phi_1_inv == expected_phi_1_inv
        @test migrated.transfer_result.column_transformation == old.column_transformation
        @test migrated.transfer_result.input_blocks.Hz == old.input_blocks.Hz
        @test migrated.transfer_result.standard_blocks.Hx == old.standard_blocks.Hx
        @test migrated.transfer_result.max_ele_phi_1 == old.max_ele_phi_1
        @test migrated.transfer_result.max_ele_phi_1_inv == old.max_ele_phi_1_inv
    end
end

@testset "v1 failed partial cache migration preserves metrics" begin
    F = GF(2)
    Rxy, (x, y) = laurent_polynomial_ring(F, ["x", "y"])

    mktempdir() do tmpdir
        src_dir = joinpath(tmpdir, "v1")
        dst_dir = joinpath(tmpdir, "v2")
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
        @test migrated.transfer_result.phi_1 === nothing
        @test migrated.transfer_result.phi_1_inv === nothing
    end
end

@testset "load_decoupled_toric_case rejects v1 payload with migration message" begin
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
        err = try
            load_decoupled_toric_case(path)
            nothing
        catch caught
            caught
        end
        @test err isa ArgumentError
        @test occursin("migrate_cached_toric_cases_to_decoupled", sprint(showerror, err))
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
