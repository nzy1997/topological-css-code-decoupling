using Test
using ToricBuilder

const TORICBUILDER_DIR = dirname(dirname(@__DIR__))
const REPO_ROOT = dirname(dirname(TORICBUILDER_DIR))
const EXAMPLE_DIR = joinpath(TORICBUILDER_DIR, "example")
const SCRIPT_DIR = joinpath(EXAMPLE_DIR, "scripts")

const RETAINED_SCRIPT_FILES = (
    "color_code.jl",
    "decouple_bbcodes.jl",
    "decouple_additional_bb_codes.jl",
    "gross_1441212_anyon_period_analysis.jl",
    "laurent_gaussian_degree_growth.jl",
    "plot_area_comparison.jl",
    "plot_decoding_result_from_data.jl",
    "transported_cnot_support.jl",
)

function tracked_example_files()
    Sys.which("git") === nothing && return nothing
    ispath(joinpath(REPO_ROOT, ".git")) || return nothing
    output = read(`git -C $REPO_ROOT ls-files -- julia/ToricBuilder/example`, String)
    prefix = "julia/ToricBuilder/"
    return [
        startswith(path, prefix) ? path[length(prefix)+1:end] : path
        for path in filter(!isempty, split(output, '\n'))
    ]
end

@testset "script examples define guarded entrypoints" begin
    for filename in RETAINED_SCRIPT_FILES
        path = joinpath(SCRIPT_DIR, filename)
        @test isfile(path)
        source = isfile(path) ? read(path, String) : ""
        @test occursin("abspath(PROGRAM_FILE) == @__FILE__", source)
    end
end

@testset "color-code example prints a reproducible summary" begin
    path = joinpath(SCRIPT_DIR, "color_code.jl")
    script_mod = Module(:ColorCodeExample)
    Base.include(script_mod, path)

    output_text, result = mktemp() do _, output
        result = redirect_stdout(output) do
            Core.eval(script_mod, :(main()))
        end
        flush(output)
        seekstart(output)
        return read(output, String), result
    end

    @test result.l == 3
    @test size(result.A_cg) == (18, 36)
    @test result.result.product_state_num == 1
    @test result.result.toric_num == 2
    @test output_text == join(
        (
            "Color-code toric-form reproduction",
            "period L: 3",
            "coarse-grained matrix size: 18 x 36",
            "decomposition: 1 product-state sector + 2 toric sectors",
            "verification: ok",
        ),
        '\n',
    ) * "\n"
end

@testset "exported example files match the release allowlist" begin
    tracked = tracked_example_files()
    if isnothing(tracked)
        @test true
    else
        expected = Set(joinpath("example", "scripts", filename) for filename in RETAINED_SCRIPT_FILES)
        @test Set(tracked) == expected
    end
end

@testset "Laurent Gaussian degree-growth example remains reproducible" begin
    path = joinpath(SCRIPT_DIR, "laurent_gaussian_degree_growth.jl")
    script_mod = Module(:LaurentGaussianDegreeGrowthExample)
    Base.include(script_mod, path)

    output_text, result = mktemp() do _, output
        result = redirect_stdout(output) do
            Core.eval(script_mod, :(main()))
        end
        flush(output)
        seekstart(output)
        return read(output, String), result
    end

    @test result.equation_Hz * result.phi_1 == result.Hzt
    @test result.trace.input_stats == (max_degree=2, max_terms=2)
    @test result.trace.final_stats == (max_degree=35, max_terms=96)
    @test result.fraction_result.stats == (
        numerator_degree=7,
        denominator_degree=6,
        numerator_terms=8,
        denominator_terms=10,
    )
    @test result.solution_stats == (max_degree=1, max_terms=3)
    @test occursin("Direct polynomial elimination: degree 2 -> 35, terms 2 -> 96", output_text)
end

@testset "example top-level contains scripts only" begin
    @test isdir(SCRIPT_DIR)
    @test !ispath(joinpath(EXAMPLE_DIR, "generated"))
    @test !ispath(joinpath(EXAMPLE_DIR, "reference"))
    @test !ispath(joinpath(EXAMPLE_DIR, "testing_case_cache"))
end
