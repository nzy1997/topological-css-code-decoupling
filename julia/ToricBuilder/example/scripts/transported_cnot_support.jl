using ToricBuilder
using ToricBuilder: max_column_monomial_count
using Oscar
using JSON

const CNOT_TORICBUILDER_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const CNOT_REPOSITORY_ROOT = normpath(joinpath(CNOT_TORICBUILDER_ROOT, "..", ".."))
const CNOT_OUTPUT_PATH = joinpath(CNOT_REPOSITORY_ROOT, "build", "reproduction", "transported_cnot_support.json")

function transported_cnot_case_specs()
    F = GF(2)
    _, (x, y) = laurent_polynomial_ring(F, ["x", "y"])
    return [
        (; case_id="case_001", f="1 + x + x*y", g="1 + y + x*y", poly_vector=[1 + x + x*y, 1 + y + x*y], toric_num=2),
        (; case_id="case_002", f="1 + x + x^-1*y", g="1 + y + x*y", poly_vector=[1 + x + x^-1*y, 1 + y + x*y], toric_num=3),
        (; case_id="case_003", f="1 + x + x^2", g="1 + y + x^2", poly_vector=[1 + x + x^2, 1 + y + x^2], toric_num=2),
        (; case_id="case_004", f="1 + x + x^-1", g="1 + y + y^-1", poly_vector=[1 + x + x^-1, 1 + y + y^-1], toric_num=4),
        (; case_id="case_005", f="1 + x + x*y", g="1 + y + x*y^-1", poly_vector=[1 + x + x*y, 1 + y + x*y^-1], toric_num=3),
        (; case_id="case_006", f="1 + x + x^-1", g="1 + y + x^3*y^2", poly_vector=[1 + x + x^-1, 1 + y + x^3*y^2], toric_num=4),
    ]
end

function toric_pair_cnot_matrix(result, target::Int, control::Int)
    1 <= target <= result.toric_num || throw(ArgumentError("target must index a toric-code sector"))
    1 <= control <= result.toric_num || throw(ArgumentError("control must index a toric-code sector"))
    target != control || throw(ArgumentError("control and target must be distinct"))
    isnothing(result.phi_1_inv) && throw(ArgumentError("transported CNOTs require compute_inverse=true"))

    qsize = size(result.phi_1, 1)
    elementary = identity_matrix(base_ring(result.phi_1), qsize)
    offset = 2 * result.product_state_num
    for component in 0:1
        target_column = offset + 2 * target - 1 + component
        control_column = offset + 2 * control - 1 + component
        elementary[:, target_column] += elementary[:, control_column]
    end
    return elementary
end

function transported_cnot_measurements(result)
    measurements = NamedTuple[]
    for target in 1:result.toric_num
        for control in 1:result.toric_num
            target == control && continue
            elementary = toric_pair_cnot_matrix(result, target, control)
            transported = result.phi_1 * elementary * result.phi_1_inv
            push!(
                measurements,
                (;
                    control,
                    target,
                    support=max_column_monomial_count(transported),
                ),
            )
        end
    end
    return measurements
end

function transported_cnot_payload()
    cases = Dict{String, Any}[]
    for spec in transported_cnot_case_specs()
        result = build_toric_form(spec.poly_vector; show_progress=false, compute_inverse=true)
        result.toric_num == spec.toric_num || error("unexpected toric sector count for $(spec.case_id)")
        measurements = [
            Dict("control" => row.control, "target" => row.target, "support" => row.support)
            for row in transported_cnot_measurements(result)
        ]
        push!(
            cases,
            Dict(
                "case_id" => spec.case_id,
                "f" => spec.f,
                "g" => spec.g,
                "toric_sectors" => spec.toric_num,
                "measurements" => measurements,
            ),
        )
    end

    return Dict(
        "metadata" => Dict(
            "schema_version" => 1,
            "transport_formula" => "Q = phi_1 * E * phi_1_inv",
            "support_definition" => "maximum total Laurent-monomial count in a column of Q",
            "source" => "Supplementary table: support spreading of transported inter-copy CNOTs",
            "software" => "Julia 1.12.5 and ToricBuilder.jl",
        ),
        "cases" => cases,
    )
end

function save_transported_cnot_support(output_path::AbstractString=CNOT_OUTPUT_PATH)
    payload = transported_cnot_payload()
    mkpath(dirname(output_path))
    open(output_path, "w") do io
        JSON.print(io, payload, 2)
        println(io)
    end
    return output_path
end

if abspath(PROGRAM_FILE) == @__FILE__
    output_path = save_transported_cnot_support()
    println("Saved transported-CNOT support data to $(output_path)")
end
