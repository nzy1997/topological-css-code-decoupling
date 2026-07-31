using Test
using JSON

const TORICBUILDER_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const RESULTS_DIR = joinpath(TORICBUILDER_ROOT, "results")

const AREA_PLOT_TEST_MODULE = Module(:AreaPlotPublishedResults)
Base.include(AREA_PLOT_TEST_MODULE, joinpath(TORICBUILDER_ROOT, "example", "scripts", "plot_area_comparison.jl"))
const DECODER_PLOT_TEST_MODULE = Module(:DecoderPlotPublishedResults)
Base.include(
    DECODER_PLOT_TEST_MODULE,
    joinpath(TORICBUILDER_ROOT, "example", "scripts", "plot_decoding_result_from_data.jl"),
)

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

    @test payload["schema_version"] == 1
    @test payload["dataset_class"] == "archived"
    @test payload["seed_status"] == "not_recorded"
    @test payload["distances"] == [4, 6, 8, 10]
    @test payload["statistics"]["interval"] == "profile likelihood"
    @test payload["statistics"]["likelihood_ratio"] == 1000.0
    @test payload["provenance"]["environment"]["status"] == "not_recorded"
    @test payload["provenance"]["methods"]["unitary"]["matching_weights"] == "unweighted"
    @test payload["provenance"]["methods"]["bposd"]["prior"] == "matched_per_point"
    @test payload["provenance"]["timing"]["included"] == ["decoder call"]

    for distance in payload["distances"]
        key = string(distance)
        @test haskey(payload["logical_error"][key], "unitary")
        @test haskey(payload["decoding_time"][key], "unitary")
        @test length(payload["logical_error"][key]["bposd"]["pvec"]) == 17
        @test length(payload["logical_error"][key]["unitary"]["pvec"]) == 10
        obsolete_series = "t" * "toric"
        @test !haskey(payload["logical_error"][key], obsolete_series)
        @test !haskey(payload["decoding_time"][key], obsolete_series)
    end

    unitary = payload["logical_error"]["4"]["unitary"]
    @test unitary["nsim"][1] == 1_088_545
    @test unitary["error_count"][1] == 2_040
    @test unitary["interval"]["low"][1] ≈ 0.0017241678059600049
    @test unitary["interval"]["high"][1] ≈ 0.0020323834596710486
    @test payload["decoding_time"]["4"]["unitary"]["time_res"][1] ==
        5.317783355712891e-6

    zero_failure = Dict(
        "interval" => Dict(
            "estimate" => [0.0],
            "low" => [0.0],
            "high" => [6.907752893026142e-7],
            "likelihood_ratio" => 1000.0,
        ),
    )
    plotted = DECODER_PLOT_TEST_MODULE.logical_plot_series(zero_failure)
    @test plotted.center == [6.907752893026142e-10]
    @test plotted.lower == [0.0]
    @test plotted.upper ≈ [6.900845140133117e-7]
    @test zero_failure["interval"]["high"][1] >
        DECODER_PLOT_TEST_MODULE.LOGICAL_ERROR_PANEL_LIMITS.y[1]
end
