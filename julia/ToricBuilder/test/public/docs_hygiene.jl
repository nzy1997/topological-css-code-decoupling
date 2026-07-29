using Test

const DOCS_TORICBUILDER_ROOT = dirname(dirname(@__DIR__))
const DOCS_REPO_ROOT = dirname(dirname(DOCS_TORICBUILDER_ROOT))

const RETAINED_DOC_PATHS = (
    joinpath("docs", "workflows", "toric_form_workflow.md"),
    joinpath("docs", "expert", "laurent_gaussian_degree_growth.md"),
)

function tracked_doc_files()
    Sys.which("git") === nothing && return nothing
    ispath(joinpath(DOCS_REPO_ROOT, ".git")) || return nothing
    output = read(`git -C $DOCS_REPO_ROOT ls-files -- julia/ToricBuilder/docs`, String)
    prefix = "julia/ToricBuilder/"
    return [
        startswith(path, prefix) ? path[length(prefix)+1:end] : path
        for path in filter(!isempty, split(output, '\n'))
    ]
end

@testset "README documents the retained API and examples" begin
    readme_path = joinpath(DOCS_TORICBUILDER_ROOT, "README.md")
    @test isfile(readme_path)
    readme = isfile(readme_path) ? read(readme_path, String) : ""

    for heading in ("## Mainline Workflow", "## API Layers", "## Documentation Map", "## Examples")
        @test occursin(heading, readme)
    end
    for name in ("excitation_matrix", "coarse_grain", "find_L", "build_toric_form")
        @test occursin("`$name`", readme)
    end

    @test occursin("docs/workflows/toric_form_workflow.md", readme)
    @test occursin("docs/expert/laurent_gaussian_degree_growth.md", readme)
    @test occursin("example/scripts/color_code.jl", readme)
    @test !occursin("matching", lowercase(readme))
    @test !occursin("decode bundle", lowercase(readme))
end

@testset "purpose-specific docs match the allowlist" begin
    for path in RETAINED_DOC_PATHS
        @test isfile(joinpath(DOCS_TORICBUILDER_ROOT, path))
    end
    tracked = tracked_doc_files()
    isnothing(tracked) ? (@test true) : (@test Set(tracked) == Set(RETAINED_DOC_PATHS))
end

@testset "public docs avoid internal helper names" begin
    for path in (joinpath("README.md"), joinpath("docs", "workflows", "toric_form_workflow.md"))
        doc = read(joinpath(DOCS_TORICBUILDER_ROOT, path), String)
        @test !occursin("_prepare_transfer_to_toric_form_data", doc)
        @test !occursin("_project_laurent_matrix_to_f2", doc)
        @test !occursin("_expand_binary_matrix_domain", doc)
        @test !occursin("_expand_binary_matrix_both_sides", doc)
    end
end
