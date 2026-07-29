# topological-css-code-decoupling

Research software and reproducibility data for the article *Decoupling 2D translation-invariant topological CSS codes*.

The repository contains two complementary implementations:

- [`python/`](python/) is the SageMath/Python paper companion for period
  search, coarse graining, decoupling maps, worked notebooks, BB-instance
  validation, and the unitary-decoupling decoder.
- [`julia/ToricBuilder/`](julia/ToricBuilder/) is the Julia/Oscar
  implementation used for toric-form construction, the BB-code result tables,
  scaling analysis, and transported-CNOT calculations.

They implement the same mathematical program but are not API-equivalent
packages.

## Setup

Python/SageMath:

```bash
cd python
sage -pip install -e '.[test]'
```

Julia/Oscar:

```bash
julia --project=julia/ToricBuilder -e 'using Pkg; Pkg.instantiate()'
```

## Verification

Run the repository checks with:

```bash
make verify
make test-julia
make test-python
```

The quick paper-companion reproduction writes only to the ignored
`build/reproduction/` directory:

```bash
make reproduce
```

See [`docs/reproduction.md`](docs/reproduction.md) for individual commands,
expected outputs, full-validation options, and the numerical environment.

## Published Data

The tracked Julia result inputs are under
[`julia/ToricBuilder/results/`](julia/ToricBuilder/results/):

- two complete BB-code decoupling tables;
- the decoder benchmark data used by the article figure;
- support spreading for all 40 transported CNOTs reported in the supplement.

Generated plots, caches, and executed notebooks are intentionally not tracked.

## Authorship

- Zhongyi Ni: Julia/Oscar implementation.
- Mingxin He: Python/SageMath implementation.

Other article authors contributed mathematical and implementation ideas but
did not write the code in this repository.

## Citation And License

Use [`CITATION.cff`](CITATION.cff) for the software and preferred article
citation metadata. The repository is distributed under the MIT license; see
[`LICENSE`](LICENSE).
