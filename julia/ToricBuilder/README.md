# ToricBuilder.jl

ToricBuilder is a Julia/Oscar research package for constructing toric-form
data for two-dimensional translation-invariant CSS codes over `GF(2)`. It uses
Laurent-polynomial presentations and exact algebraic operations.

## Mainline Workflow

```julia
using ToricBuilder
using Oscar

Rxy, (x, y) = laurent_polynomial_ring(GF(2), ["x", "y"])
poly_vec = [1 + x + x*y, 1 + y + x*y]

l = find_L([x, y], ideal(Rxy, poly_vec))
A = excitation_matrix(poly_vec)
A_cg = coarse_grain(A, [x, y], [l, l])
toric = build_toric_form(poly_vec; show_progress=false, compute_inverse=true)

@assert check_result(toric, toric.input_matrix; require_inverse=true)
```

The mainline entry points are:

- `excitation_matrix`: construct the Laurent excitation matrix;
- `coarse_grain`: project to a finite periodic cell;
- `find_L`: compute the period bound used by the workflow;
- `build_toric_form`: construct the standard form and basis changes.

See [`docs/workflows/toric_form_workflow.md`](docs/workflows/toric_form_workflow.md)
for a longer walkthrough.

## API Layers

### Mainline

- `excitation_matrix`
- `coarse_grain`
- `find_L`
- `build_toric_form`

### Expert

- `coarse_graining` and `make_L_table`
- `to_toric_form` and `check_result`
- `solve_laurent_linear`
- `construct_quotient_representation`
- `construct_anyon_translation_representation`
- `construct_period_table`
- `factor_toric_block`, `project_to_finite_code`, and `split_product_state`
- `capture_toric_form_debug_matrices`
- `DecoupledToricCase` and the associated save/load/list helpers

Underscored functions remain implementation details.

## Documentation Map

- [`docs/workflows/toric_form_workflow.md`](docs/workflows/toric_form_workflow.md):
  toric-form construction and verification.
- [`docs/expert/laurent_gaussian_degree_growth.md`](docs/expert/laurent_gaussian_degree_growth.md):
  an instrumented comparison between exact Laurent solving and two naive
  Gaussian-elimination alternatives.
- [`../../docs/reproduction.md`](../../docs/reproduction.md): paper-facing
  result tables, figures, and numerical checks.

## Examples

Runnable scripts live under [`example/scripts/`](example/scripts/):

- [`color_code.jl`](example/scripts/color_code.jl): a small toric-form example;
- [`decouple_bbcodes.jl`](example/scripts/decouple_bbcodes.jl): the main BB set;
- [`decouple_additional_bb_codes.jl`](example/scripts/decouple_additional_bb_codes.jl):
  the additional published BB set;
- [`plot_area_comparison.jl`](example/scripts/plot_area_comparison.jl): combined
  time/locality scaling plot;
- [`plot_decoding_result_from_data.jl`](example/scripts/plot_decoding_result_from_data.jl):
  decoder benchmark figure;
- [`transported_cnot_support.jl`](example/scripts/transported_cnot_support.jl):
  all transported inter-copy CNOT support values;
- [`gross_1441212_anyon_period_analysis.jl`](example/scripts/gross_1441212_anyon_period_analysis.jl):
  Gross-code period analysis;
- [`laurent_gaussian_degree_growth.jl`](example/scripts/laurent_gaussian_degree_growth.jl):
  symbolic degree-growth diagnostic.

Tracked numerical inputs live in [`results/`](results/). Every generated file
is written under the repository-level `build/reproduction/` directory.

## Testing

From the repository root:

```bash
julia --project=julia/ToricBuilder -e 'using Pkg; Pkg.test()'
```

The full test suite includes a fresh reconstruction of the six BB codes in the
transported-CNOT table.

## Dependencies

Direct dependencies are declared in [`Project.toml`](Project.toml), including
Oscar, AbstractAlgebra, JSON, CairoMakie, DataFrames, GLM, and LaTeXStrings.
