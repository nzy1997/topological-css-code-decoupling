using Test
using JSON

const TORICBUILDER_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const RESULTS_DIR = joinpath(TORICBUILDER_ROOT, "results")

const AREA_PLOT_TEST_MODULE = Module(:AreaPlotPublishedResults)
Base.include(AREA_PLOT_TEST_MODULE, joinpath(TORICBUILDER_ROOT, "example", "scripts", "plot_area_comparison.jl"))

function markdown_data_row_count(path)
    return count(eachline(path)) do line
        startswith(line, "|") && !startswith(line, "|---") && !occursin("Case ID", line)
    end
end

@testset "published BB scaling data" begin
    main_path = joinpath(RESULTS_DIR, "bb_code_decoupling_results.md")
    additional_path = joinpath(RESULTS_DIR, "additional_bb_code_decoupling_results.md")

    @test markdown_data_row_count(main_path) == 54
    @test markdown_data_row_count(additional_path) == 99

    main_data = AREA_PLOT_TEST_MODULE.parse_markdown_table(main_path)
    additional_data = AREA_PLOT_TEST_MODULE.parse_markdown_table(additional_path)
    combined = AREA_PLOT_TEST_MODULE.parse_markdown_tables([main_path, additional_path])

    @test length(main_data.areas) == 31
    @test length(additional_data.areas) == 29
    @test length(combined.qs) == 60
    @test combined.qs == 2 .* combined.areas

    time_data = AREA_PLOT_TEST_MODULE.positive_pairs(combined.qs, combined.times)
    degree_data = AREA_PLOT_TEST_MODULE.positive_pairs(combined.qs, combined.degrees)
    @test length(time_data.xs) == 60
    @test length(degree_data.xs) == 59

    time_slope, _ = AREA_PLOT_TEST_MODULE.fit_power_law(time_data.xs, time_data.ys)
    degree_slope, _ = AREA_PLOT_TEST_MODULE.fit_power_law(degree_data.xs, degree_data.ys)
    @test round(time_slope; digits=2) == 1.86
    @test round(degree_slope; digits=2) == 0.65

    plot_source = read(joinpath(TORICBUILDER_ROOT, "example", "scripts", "plot_area_comparison.jl"), String)
    @test occursin(raw"\deg(\psi_1^{-1})", plot_source)

    additional_source = read(joinpath(TORICBUILDER_ROOT, "example", "scripts", "decouple_additional_bb_codes.jl"), String)
    additional_table = read(additional_path, String)
    @test occursin("10.1038/s41586-024-07107-7", additional_source)
    @test occursin("Doubled Gross code [[288,12,18]] (Bravyi et al., 2024)", additional_table)
end

@testset "decoder benchmark provenance" begin
    path = joinpath(RESULTS_DIR, "decoding_benchmark.json")
    payload = JSON.parsefile(path)
    metadata = payload["metadata"]

    @test metadata["schema_version"] == 1
    @test metadata["code_family"] == "bivariate bicycle"
    @test metadata["check_polynomials"] == ["1 + x + x^-1*y", "1 + y + x*y"]
    @test metadata["toric_code_distances"] == [4, 6, 8, 10]
    @test occursin("Apple M4", metadata["runtime_environment"])
    @test !isempty(metadata["source_pipeline"])
    @test payload["distances"] == metadata["toric_code_distances"]

    for distance in payload["distances"]
        key = string(distance)
        @test haskey(payload["logical_error"][key], "unitary_decouple")
        @test haskey(payload["decoding_time"][key], "unitary_decouple")
        obsolete_series = "t" * "toric"
        @test !haskey(payload["logical_error"][key], obsolete_series)
        @test !haskey(payload["decoding_time"][key], obsolete_series)
    end
end
