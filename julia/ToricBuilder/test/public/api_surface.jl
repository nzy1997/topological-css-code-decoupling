using Test
using ToricBuilder

@testset "public API surface" begin
    mainline_names = (
        :excitation_matrix,
        :coarse_grain,
        :find_L,
        :build_toric_form,
    )

    expert_names = (
        :coarse_graining,
        :make_L_table,
        :to_toric_form,
        :check_result,
        :solve_laurent_linear,
        :construct_quotient_representation,
        :construct_anyon_translation_representation,
        :construct_period_table,
        :factor_toric_block,
        :project_to_finite_code,
        :split_product_state,
        :capture_toric_form_debug_matrices,
        :DecoupledToricCase,
        :DECOUPLED_TORIC_CASE_FORMAT_VERSION,
        :build_decoupled_toric_case,
        :save_decoupled_toric_case,
        :load_decoupled_toric_case,
        :list_decoupled_toric_cases,
        :decoupled_toric_case_path,
        :migrate_cached_toric_cases_to_decoupled,
    )

    for name in mainline_names
        @test isdefined(ToricBuilder, name)
        @test Base.isexported(ToricBuilder, name)
    end

    for name in expert_names
        @test isdefined(ToricBuilder, name)
        @test Base.isexported(ToricBuilder, name)
    end

    for function_name in (
        :excitation_matrix,
        :coarse_grain,
        :find_L,
        :build_toric_form,
        :to_toric_form,
        :check_result,
        :solve_laurent_linear,
        :construct_quotient_representation,
        :construct_anyon_translation_representation,
        :construct_period_table,
        :factor_toric_block,
        :split_product_state,
        :project_to_finite_code,
        :capture_toric_form_debug_matrices,
        :build_decoupled_toric_case,
        :save_decoupled_toric_case,
        :load_decoupled_toric_case,
        :list_decoupled_toric_cases,
        :decoupled_toric_case_path,
        :migrate_cached_toric_cases_to_decoupled,
    )
        @test getfield(ToricBuilder, function_name) isa Function
    end

    for legacy_name in (
        :transfer_to_toric_form,
        :to_code,
        :guassian_elemination!,
        :construct_representation_by_mat,
        :construct_period_table_by_representation,
        :cg_toric,
        :CachedToricCase,
        :CACHED_TORIC_CASE_FORMAT_VERSION,
        :build_cached_toric_case,
        :save_cached_toric_case,
        :load_cached_toric_case,
        :list_cached_toric_cases,
        :cached_toric_case_path,
    )
        @test !isdefined(ToricBuilder, legacy_name)
        @test !Base.isexported(ToricBuilder, legacy_name)
    end
end
