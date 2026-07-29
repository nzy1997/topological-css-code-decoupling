using Test

const RELEASE_TORICBUILDER_DIR = dirname(dirname(@__DIR__))
const RELEASE_REPO_ROOT = dirname(dirname(RELEASE_TORICBUILDER_DIR))
const RELEASE_SCRIPT_DIR = joinpath(RELEASE_TORICBUILDER_DIR, "example", "scripts")
const RELEASE_RESULTS_DIR = joinpath(RELEASE_TORICBUILDER_DIR, "results")

const RELEASE_RETAINED_SCRIPT_FILES = (
    "color_code.jl",
    "decouple_bbcodes.jl",
    "decouple_additional_bb_codes.jl",
    "gross_1441212_anyon_period_analysis.jl",
    "laurent_gaussian_degree_growth.jl",
    "plot_area_comparison.jl",
    "plot_decoding_result_from_data.jl",
    "transported_cnot_support.jl",
)

const RELEASE_RETAINED_RESULT_FILES = (
    "additional_bb_code_decoupling_results.md",
    "bb_code_decoupling_results.md",
    "decoding_benchmark.json",
    "transported_cnot_support.json",
)

function release_tracked_files(paths...)
    Sys.which("git") === nothing && return nothing
    ispath(joinpath(RELEASE_REPO_ROOT, ".git")) || return nothing
    rooted_paths = [joinpath("julia", "ToricBuilder", path) for path in paths]
    output = read(`git -C $RELEASE_REPO_ROOT ls-files -- $rooted_paths`, String)
    prefix = "julia/ToricBuilder/"
    return [
        startswith(path, prefix) ? path[length(prefix)+1:end] : path
        for path in filter(!isempty, split(output, '\n'))
    ]
end

@testset "release example and result surfaces are stable" begin
    for filename in RELEASE_RETAINED_SCRIPT_FILES
        @test isfile(joinpath(RELEASE_SCRIPT_DIR, filename))
    end
    for filename in RELEASE_RETAINED_RESULT_FILES
        @test isfile(joinpath(RELEASE_RESULTS_DIR, filename))
    end

    tracked_scripts = release_tracked_files("example")
    if isnothing(tracked_scripts)
        @test true
    else
        expected = Set(joinpath("example", "scripts", filename) for filename in RELEASE_RETAINED_SCRIPT_FILES)
        @test Set(tracked_scripts) == expected
    end

    tracked_results = release_tracked_files("results")
    if isnothing(tracked_results)
        @test true
    else
        expected = Set(joinpath("results", filename) for filename in RELEASE_RETAINED_RESULT_FILES)
        @test Set(tracked_results) == expected
    end
end

@testset "example scripts use guarded entrypoints" begin
    for filename in RELEASE_RETAINED_SCRIPT_FILES
        source = read(joinpath(RELEASE_SCRIPT_DIR, filename), String)
        @test occursin("abspath(PROGRAM_FILE) == @__FILE__", source)
    end
end

@testset "Julia metadata declares the reproducibility environment" begin
    project = read(joinpath(RELEASE_TORICBUILDER_DIR, "Project.toml"), String)
    @test occursin("Zhongyi Ni", project)
    @test occursin("MIT", project)
    for dependency in ("CairoMakie", "DataFrames", "GLM", "LaTeXStrings")
        @test occursin(dependency, project)
    end
    @test !occursin("DelimitedFiles", project)
end

@testset "tracked tree excludes generated clutter" begin
    tracked = release_tracked_files("example", "results", "test", "bench")
    if isnothing(tracked)
        @test true
    else
        @test !any(endswith(".jls"), tracked)
        @test !any(path -> occursin("debug_matrix_outputs", path), tracked)
        @test !any(path -> occursin("testing_case_cache", path), tracked)
        @test !any(path -> startswith(path, "bench/"), tracked)
    end
end
