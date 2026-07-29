module ToricBuilder

using AbstractAlgebra
using Oscar
using LinearAlgebra
import JSON

export excitation_matrix
export coarse_graining, coarse_grain, find_L, make_L_table
export to_toric_form, check_result
export build_toric_form
export project_to_finite_code, split_product_state
export factor_toric_block
export solve_laurent_linear
export capture_toric_form_debug_matrices
export construct_quotient_representation, construct_period_table
export construct_anyon_translation_representation
export DecoupledToricCase, DECOUPLED_TORIC_CASE_FORMAT_VERSION
export build_decoupled_toric_case, save_decoupled_toric_case, load_decoupled_toric_case
export list_decoupled_toric_cases, decoupled_toric_case_path
export migrate_cached_toric_cases_to_decoupled

include("core/excitation_matrix.jl")
include("toric_form/coarse_graining.jl")
include("toric_form/toric_reduction.jl")
include("core/laurent_linear_solve.jl")
include("toric_form/toric_workflow.jl")
include("toric_form/toric_factorization.jl")
include("export/toric_case_cache.jl")
include("core/quotient_representation.jl")
end
