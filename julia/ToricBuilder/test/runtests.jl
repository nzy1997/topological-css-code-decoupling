using ToricBuilder
using Test

@testset "public" begin
    include("public/api_surface.jl")
    include("public/examples.jl")
    include("public/docs_hygiene.jl")
    include("public/toric_form_workflow.jl")
    include("public/published_results.jl")
end

@testset "expert" begin
    include("expert/laurent_linear_solve.jl")
    include("expert/excitation_matrix.jl")
    include("expert/quotient_representation.jl")
    include("expert/toric_factorization.jl")
    include("expert/toric_case_cache.jl")
    include("expert/decouple_bbcodes.jl")
    include("expert/decouple_additional_bb_codes.jl")
    include("expert/transported_cnot_support.jl")
end

@testset "internal" begin
    include("internal/toric_reduction.jl")
    include("internal/coarse_graining.jl")
end
