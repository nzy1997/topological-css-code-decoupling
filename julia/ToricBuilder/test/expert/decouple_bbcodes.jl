using Test
using ToricBuilder

@testset "run_and_save writes case id column to markdown" begin
    script_path = joinpath(dirname(dirname(@__DIR__)), "example", "scripts", "decouple_bbcodes.jl")
    source = read(script_path, String)
    source = replace(source, r"\nrun_and_save\(ab_list3\)\s*$" => "\n")
    @test occursin("\"decouple_bbcodes_cache_v3\"", source)
    @test occursin("max_area::Int=2000", source)
    @test occursin("warmup::Bool=true", source)
    @test occursin("capture_debug::Bool=false", source)
    @test occursin("check_cache::Bool=false", source)
    @test occursin("reuse_rows_path=nothing", source)

    harness = Module(:TestingCasesHarness)
    Core.eval(harness, :(using ToricBuilder))
    Core.eval(harness, :(using Oscar))
    Base.include_string(harness, source, script_path)

    stub_case = ToricBuilder.DecoupledToricCase(
        ToricBuilder.DECOUPLED_TORIC_CASE_FORMAT_VERSION,
        "case_001",
        :ok,
        Dict{String, Any}(),
        nothing,
        (
            l=3,
            u_rel="x",
            v_rel="y",
            area=4,
            A_size=5,
            product_state_num=6,
            toric_num=7,
            solving_time=1.2345,
            max_ele_psi_1_inverse="x*y",
            max_degree_psi_1_inverse=8,
            max_ele_psi_1=nothing,
            max_degree_psi_1=nothing,
            max_column_monomial_count_psi_1_inverse=10,
            max_column_monomial_count_psi_1=nothing,
        ),
        nothing,
        "0",
        Dict{String, Any}(),
    )

    no_inverse_case = ToricBuilder.DecoupledToricCase(
        ToricBuilder.DECOUPLED_TORIC_CASE_FORMAT_VERSION,
        "case_001",
        :ok,
        Dict{String, Any}(),
        nothing,
        (psi_1=nothing, column_transformation=nothing),
        nothing,
        "0",
        Dict{String, Any}(),
    )
    inverse_case = ToricBuilder.DecoupledToricCase(
        ToricBuilder.DECOUPLED_TORIC_CASE_FORMAT_VERSION,
        "case_001",
        :ok,
        Dict{String, Any}(),
        nothing,
        (psi_1=:psi_1, column_transformation=:column_transformation),
        nothing,
        "0",
        Dict{String, Any}(),
    )
    cache_satisfies_request = getfield(harness, :_decoupled_cache_satisfies_request)
    @test cache_satisfies_request(no_inverse_case; compute_inverse=false, capture_debug=false)
    @test !cache_satisfies_request(no_inverse_case; compute_inverse=false, capture_debug=true)
    @test !cache_satisfies_request(no_inverse_case; compute_inverse=true, capture_debug=false)
    @test cache_satisfies_request(inverse_case; compute_inverse=true, capture_debug=false)

    mktempdir() do tmpdir
        cache_path = joinpath(tmpdir, "case_001.jls")
        poly_vector = [1 + getfield(harness, :x) + getfield(harness, :x) * getfield(harness, :y), 1 + getfield(harness, :y) + getfield(harness, :x) * getfield(harness, :y)]
        default_case = ToricBuilder.build_decoupled_toric_case("case_001", poly_vector; show_progress=false)
        @test default_case.transfer_result.psi_1 === nothing
        ToricBuilder.save_decoupled_toric_case(cache_path, default_case)

        reused_case = getfield(harness, :_load_or_build_decoupled_case)(
            "case_001",
            poly_vector,
            cache_path;
            metadata=Dict{String, Any}(),
            overwrite=false,
            show_progress=false,
            compute_inverse=true,
        )

        @test reused_case.transfer_result.psi_1 === nothing
        @test ToricBuilder.load_decoupled_toric_case(cache_path).transfer_result.psi_1 === nothing

        refreshed_case = getfield(harness, :_load_or_build_decoupled_case)(
            "case_001",
            poly_vector,
            cache_path;
            metadata=Dict{String, Any}(),
            overwrite=false,
            show_progress=false,
            compute_inverse=true,
            check_cache=true,
        )

        @test !isnothing(refreshed_case.transfer_result.psi_1)
        @test !isnothing(ToricBuilder.load_decoupled_toric_case(cache_path).transfer_result.psi_1)
    end

    Core.eval(harness, quote
        warmup_calls = Ref(0)
        warmup_show_progress = Ref{Union{Nothing, Bool}}(nothing)
        load_or_build_calls = Ref(0)
        function _warmup_decoupling(poly_vector::Vector; show_progress::Bool, compute_inverse::Bool, max_area::Int, max_l::Int)
            warmup_calls[] += 1
            warmup_show_progress[] = show_progress
            return nothing
        end

        function _load_or_build_decoupled_case(case_id::AbstractString, poly_vector::Vector, cache_path::AbstractString; metadata, overwrite::Bool, show_progress::Bool, compute_inverse::Bool, capture_debug::Bool=false, check_cache::Bool=false, kwargs...)
            load_or_build_calls[] += 1
            return $stub_case
        end
    end)

    mktempdir() do tmpdir
        results_path = joinpath(tmpdir, "results.md")
        cache_dir = joinpath(tmpdir, "cache")
        getfield(harness, :run_and_save)([(getfield(harness, :x), getfield(harness, :y))]; results_path=results_path, cache_dir=cache_dir)

        lines = readlines(results_path)

        @test getfield(harness, :warmup_calls)[] == 1
        @test getfield(harness, :warmup_show_progress)[] === false
        @test length(lines) == 3
        @test occursin("Case ID", lines[1])
        @test occursin("Maximum_term_psi_1_inverse", lines[1])
        @test occursin("Maximum_term_psi_1", lines[1])
        @test !occursin("Maximum_term_Q", lines[1])
        @test endswith(lines[2], "|---|")
        @test occursin("| case_001 |", lines[3])
        @test occursin("| \$x \\cdot y\$ | 8 | - | - | 10 | - |", lines[end])
    end

    mktempdir() do tmpdir
        results_path = joinpath(tmpdir, "results.md")
        cache_dir = joinpath(tmpdir, "cache")
        mkpath(cache_dir)
        write(joinpath(cache_dir, "case_001.jls"), "not a loadable cache")
        reuse_rows_path = joinpath(tmpdir, "previous.md")
        reused_row = "| old a | old b | old l | old u | old v | old area | old A | old product | old toric | old time | old phi | old deg | - | - | - | - | case_001 |"
        write(
            reuse_rows_path,
            "| h | h | h | h | h | h | h | h | h | h | h | h | h | h | h | h | Case ID |\n|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|\n$reused_row\n",
        )

        getfield(harness, :warmup_calls)[] = 0
        getfield(harness, :load_or_build_calls)[] = 0
        getfield(harness, :run_and_save)(
            [(getfield(harness, :x), getfield(harness, :y))];
            results_path=results_path,
            cache_dir=cache_dir,
            reuse_rows_path=reuse_rows_path,
            warmup=false,
        )

        lines = readlines(results_path)
        @test lines[3] == reused_row
        @test getfield(harness, :load_or_build_calls)[] == 0
        @test getfield(harness, :warmup_calls)[] == 0
    end

    mktempdir() do tmpdir
        results_path = joinpath(tmpdir, "results.md")
        cache_dir = joinpath(tmpdir, "cache")
        mkpath(cache_dir)
        write(joinpath(cache_dir, "case_001.jls"), "old failed cache")
        reuse_rows_path = joinpath(tmpdir, "previous.md")
        incomplete_row = "| old a | old b | old l | old u | old v | old area | - | - | - | old time | - | - | - | - | - | - | case_001 |"
        write(
            reuse_rows_path,
            "| h | h | h | h | h | h | h | h | h | h | h | h | h | h | h | h | Case ID |\n|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|\n$incomplete_row\n",
        )

        getfield(harness, :load_or_build_calls)[] = 0
        getfield(harness, :run_and_save)(
            [(getfield(harness, :x), getfield(harness, :y))];
            results_path=results_path,
            cache_dir=cache_dir,
            reuse_rows_path=reuse_rows_path,
            warmup=false,
        )

        lines = readlines(results_path)
        @test lines[3] != incomplete_row
        @test getfield(harness, :load_or_build_calls)[] == 1
    end

    mktempdir() do tmpdir
        results_path = joinpath(tmpdir, "results.md")
        cache_dir = joinpath(tmpdir, "cache")
        getfield(harness, :warmup_calls)[] = 0
        getfield(harness, :warmup_show_progress)[] = nothing
        getfield(harness, :run_and_save)(
            [(getfield(harness, :x), getfield(harness, :y))];
            results_path=results_path,
            cache_dir=cache_dir,
            show_progress=true,
        )

        @test getfield(harness, :warmup_calls)[] == 1
        @test getfield(harness, :warmup_show_progress)[] === true
    end

    mktempdir() do tmpdir
        results_path = joinpath(tmpdir, "results.md")
        cache_dir = joinpath(tmpdir, "cache")
        getfield(harness, :warmup_calls)[] = 0
        getfield(harness, :run_and_save)(
            [(getfield(harness, :x), getfield(harness, :y))];
            results_path=results_path,
            cache_dir=cache_dir,
            warmup=false,
        )

        @test getfield(harness, :warmup_calls)[] == 0
    end
end
