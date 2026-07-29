using Test
using ToricBuilder

@testset "additional BB code decoupling script inventory and markdown" begin
    script_path = joinpath(dirname(dirname(@__DIR__)), "example", "scripts", "decouple_additional_bb_codes.jl")
    source = read(script_path, String)
    @test occursin("\"additional_bb_code_decoupling_cache_v2\"", source)
    @test occursin("warmup::Bool=true", source)
    @test occursin("capture_debug::Bool=false", source)
    @test occursin("check_cache::Bool=false", source)
    @test occursin("reuse_rows_path=nothing", source)

    harness = Module(:AdditionalBBCodesHarness)
    Core.eval(harness, :(using ToricBuilder))
    Core.eval(harness, :(using Oscar))
    Base.include_string(harness, source, script_path)

    cases = getfield(harness, :all_additional_bb_code_cases)()
    case_ids = getproperty.(cases, :case_id)
    hx = getfield(harness, :x)
    hy = getfield(harness, :y)

    @test occursin("max_area::Int=2000", source)
    @test length(cases) == 99
    polynomial_pairs = [(string(case.poly_vector[1]), string(case.poly_vector[2])) for case in cases]
    @test length(unique(polynomial_pairs)) == length(polynomial_pairs)
    @test "bunny_001" in case_ids
    @test "bunny_005" in case_ids
    @test !("bunny_002" in case_ids)
    @test "sbb_003" in case_ids
    @test "sbb_006" in case_ids
    @test !("sbb_005" in case_ids)
    @test "cover_w6_k12_001" in case_ids
    @test "cover_w8_k14_alt_004" in case_ids
    @test "cover_w8_inc_020" in case_ids
    @test !("gross" in case_ids)
    @test "double_gross" in case_ids
    @test count(case -> get(case.metadata, "source", nothing) == "arXiv:2511.13560", cases) == 80

    bunny = cases[findfirst(==("bunny_001"), case_ids)]
    @test bunny.poly_vector[1] == 1 + hx
    @test bunny.poly_vector[2] == hx + hy
    @test bunny.metadata["source"] == "arXiv:2606.22853"

    sbb = cases[findfirst(==("sbb_003"), case_ids)]
    @test sbb.poly_vector[1] == hx^3 * hy^-2 + hx^2 * hy^-1 + hy^2 + 1
    @test sbb.poly_vector[2] == hx^2 * hy + hx^2 * hy^-1 + hx + hx * hy^-2
    @test sbb.metadata["source"] == "arXiv:2605.04151"

    cover_w6 = cases[findfirst(==("cover_w6_k12_001"), case_ids)]
    @test cover_w6.poly_vector[1] == hx^3 + hy^13 + hx^12 * hy^2
    @test cover_w6.poly_vector[2] == hy^3 + hx^7 * hy^12 + hx^14 * hy^6
    @test cover_w6.metadata["source"] == "arXiv:2511.13560"
    @test cover_w6.metadata["source_table"] == "tab:k-12-codes"

    cover_w8 = cases[findfirst(==("cover_w8_k14_alt_004"), case_ids)]
    @test cover_w8.poly_vector[1] == hx * hy^3 + 1 + hx^6 + hx^3 * hy^2
    @test cover_w8.poly_vector[2] == hx^6 * hy + hx^4 * hy + hx^3 + hx^5 * hy
    @test cover_w8.metadata["check_weight"] == 8

    double_gross = cases[findfirst(==("double_gross"), case_ids)]
    @test double_gross.poly_vector[1] == hx^3 + hy^2 + hy^7
    @test double_gross.poly_vector[2] == hy^3 + hx + hx^2

    stub_case = ToricBuilder.DecoupledToricCase(
        ToricBuilder.DECOUPLED_TORIC_CASE_FORMAT_VERSION,
        "bunny_001",
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
            max_ele_phi_1="x*y",
            max_degree_phi_1=8,
            max_ele_phi_1_inv=nothing,
            max_degree_phi_1_inv=nothing,
            max_column_monomial_count_phi_1=10,
            max_column_monomial_count_phi_1_inv=nothing,
        ),
        nothing,
        "0",
        Dict{String, Any}(),
    )

    no_inverse_case = ToricBuilder.DecoupledToricCase(
        ToricBuilder.DECOUPLED_TORIC_CASE_FORMAT_VERSION,
        "bunny_001",
        :ok,
        Dict{String, Any}(),
        nothing,
        (phi_1_inv=nothing, column_transformation=nothing),
        nothing,
        "0",
        Dict{String, Any}(),
    )
    inverse_case = ToricBuilder.DecoupledToricCase(
        ToricBuilder.DECOUPLED_TORIC_CASE_FORMAT_VERSION,
        "bunny_001",
        :ok,
        Dict{String, Any}(),
        nothing,
        (phi_1_inv=:phi_1_inv, column_transformation=:column_transformation),
        nothing,
        "0",
        Dict{String, Any}(),
    )
    cache_satisfies_request = getfield(harness, :_decoupled_cache_satisfies_request)
    @test cache_satisfies_request(no_inverse_case; compute_inverse=false, capture_debug=false)
    @test !cache_satisfies_request(no_inverse_case; compute_inverse=false, capture_debug=true)
    @test !cache_satisfies_request(no_inverse_case; compute_inverse=true, capture_debug=false)
    @test cache_satisfies_request(inverse_case; compute_inverse=true, capture_debug=false)

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
        getfield(harness, :run_and_save)([bunny, sbb]; results_path=results_path, cache_dir=cache_dir)

        lines = readlines(results_path)

        @test getfield(harness, :warmup_calls)[] == 1
        @test getfield(harness, :warmup_show_progress)[] === false
        @test length(lines) == 4
        @test startswith(lines[1], "| Code | Definition / decouple input | \$L\$ |")
        @test occursin("Maximum_term_phi_1", lines[1])
        @test occursin("Maximum_term_phi_1_inv", lines[1])
        @test !occursin("Maximum_term_Q", lines[1])
        @test occursin("bunny_001", lines[3])
        @test occursin("<br>", lines[3])
        @test occursin("| 4 | 5 | 6 | 7 |", lines[3])
        @test occursin("| \$x \\cdot y\$ | 8 | - | - | 10 | - |", lines[end])
        @test occursin("source SBB gauge triples", lines[4])
        @test occursin("(x^2, y^2, x^2*y + x)", lines[4])
        @test occursin("decouple BB pair", lines[4])
    end

    mktempdir() do tmpdir
        results_path = joinpath(tmpdir, "results.md")
        cache_dir = joinpath(tmpdir, "cache")
        mkpath(cache_dir)
        write(joinpath(cache_dir, "bunny_001.jls"), "not a loadable cache")
        reuse_rows_path = joinpath(tmpdir, "previous.md")
        reused_row = "| old code | old input | old l | old u | old v | old area | old A | old product | old toric | old time | old phi | old deg | - | - | - | - | bunny_001 |"
        write(
            reuse_rows_path,
            "| Code | Definition / decouple input | h | h | h | h | h | h | h | h | h | h | h | h | h | h | Case ID |\n|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|\n$reused_row\n",
        )

        getfield(harness, :warmup_calls)[] = 0
        getfield(harness, :load_or_build_calls)[] = 0
        getfield(harness, :run_and_save)(
            [bunny];
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
        write(joinpath(cache_dir, "bunny_001.jls"), "old failed cache")
        reuse_rows_path = joinpath(tmpdir, "previous.md")
        incomplete_row = "| old code | old input | old l | old u | old v | old area | - | - | - | old time | - | - | - | - | - | - | bunny_001 |"
        write(
            reuse_rows_path,
            "| Code | Definition / decouple input | h | h | h | h | h | h | h | h | h | h | h | h | h | h | Case ID |\n|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|\n$incomplete_row\n",
        )

        getfield(harness, :load_or_build_calls)[] = 0
        getfield(harness, :run_and_save)(
            [bunny];
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
            [bunny, sbb];
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
            [bunny, sbb];
            results_path=results_path,
            cache_dir=cache_dir,
            warmup=false,
        )

        @test getfield(harness, :warmup_calls)[] == 0
    end
end
