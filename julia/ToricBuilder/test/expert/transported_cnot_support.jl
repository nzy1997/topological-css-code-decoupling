using Test
using JSON
using ToricBuilder

const TORICBUILDER_ROOT_CNOT = normpath(joinpath(@__DIR__, "..", ".."))
const CNOT_RESULTS_PATH = joinpath(TORICBUILDER_ROOT_CNOT, "results", "transported_cnot_support.json")

include(joinpath(TORICBUILDER_ROOT_CNOT, "example", "scripts", "transported_cnot_support.jl"))

function measurement_tuples(case_payload)
    return [
        (
            control=Int(row["control"]),
            target=Int(row["target"]),
            support=Int(row["support"]),
        )
        for row in case_payload["measurements"]
    ]
end

@testset "transported CNOT published results" begin
    payload = JSON.parsefile(CNOT_RESULTS_PATH)
    @test payload["metadata"]["schema_version"] == 1
    @test payload["metadata"]["transport_formula"] == "Q = psi_1_inverse * E * psi_1"
    @test length(payload["cases"]) == 6
    @test sum(length(case["measurements"]) for case in payload["cases"]) == 40

    expected = Dict(case["case_id"] => measurement_tuples(case) for case in payload["cases"])
    specs = transported_cnot_case_specs()
    @test [spec.case_id for spec in specs] == ["case_001", "case_002", "case_003", "case_004", "case_005", "case_006"]

    for spec in specs
        result = build_toric_form(spec.poly_vector; show_progress=false, compute_inverse=true)
        @test result.toric_num == spec.toric_num
        @test transported_cnot_measurements(result) == expected[spec.case_id]
    end
end
