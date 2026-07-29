# Sage Decoder: decoupling implementation

This directory contains the SageMath implementation of the decoupling algorithm
for two-dimensional translational CSS codes and the unitary decoupling decoder
built on top of it. The public packages are:

- `isomorphism`: period search, coarse graining, and decoupling maps;
- `unitary_decouple_decoder`: the unitary decoupling decoder;
- `decoder_core`: shared finite-torus, matching, validation, and sampling
  utilities.

## Requirements and installation

Use Python 3.10 or newer inside a SageMath environment. From this directory,
install the packages and test dependencies with:

```bash
sage -pip install -e '.[test]'
```

For an already built wheel, use:

```bash
sage -pip install sage_decoder_decoupling-0.1.0-py3-none-any.whl
```

The base dependencies are NumPy and PyMatching. The optional `benchmark`
dependency group adds `ldpc` and Matplotlib for the reproduction scripts.

## Quick start: run the supplement algorithm

The supplement uses a source CSS complex and a standard decoupled complex:

```text
source degree 2 --hz_source_dagger--> source degree 1 --hx_source--> source degree 0
```

The forward maps go from the source complex to the standard complex and obey

```text
hx_standard * phi1 = phi0 * hx_source
hz_standard_dagger * phi2 = phi1 * hz_source_dagger
```

The elimination algorithm directly constructs the inverse direction:

```text
hx_source * phi1_inverse = phi0_inverse * hx_standard
hz_source_dagger * phi2_inverse
    = phi1_inverse * hz_standard_dagger
```

The following Sage code runs the algorithm on the smallest two-generator
example and checks both directions:

```python
from isomorphism import construct_excitation_map, decouple_coarse_matrix, x, y
from isomorphism.css import antipode_matrix, split_css_blocks

epsilon = construct_excitation_map(1 + x, 1 + y)
result = decouple_coarse_matrix(
    epsilon,
    num_x_checks=1,
    num_qubits=2,
    verify=True,
)

hz_source, hx_source = split_css_blocks(epsilon, 1, 2)
hz_source_dagger = antipode_matrix(hz_source)
hz_standard_dagger = antipode_matrix(result.hz_standard)

assert result.hx_standard * result.maps.phi1 == result.maps.phi0 * hx_source
assert (
    hz_standard_dagger * result.maps.phi2
    == result.maps.phi1 * hz_source_dagger
)
assert (
    hx_source * result.inverse_maps.phi1_inverse
    == result.inverse_maps.phi0_inverse * result.hx_standard
)
assert (
    hz_source_dagger * result.inverse_maps.phi2_inverse
    == result.inverse_maps.phi1_inverse * hz_standard_dagger
)
```

`decouple_coarse_matrix(...)` returns a `DecouplingResult` with exactly these
fields:

```text
source_coarse_matrix
hx_standard
hz_standard
inverse_maps
maps
diagnostics
```

`inverse_maps` always contains `phi2_inverse`, `phi1_inverse`, and
`phi0_inverse`. By default, `maps` also contains the forward maps `phi2`,
`phi1`, and `phi0`. Use `compute_forward_maps=False` when only the directly
constructed inverse maps are needed:

```python
inverse_only = decouple_coarse_matrix(
    epsilon,
    num_x_checks=1,
    num_qubits=2,
    compute_forward_maps=False,
)
assert inverse_only.maps is None
```

Computing forward maps requires Laurent matrix inversions. The `verify=True`
option adds solver-witness and chain-map checks, so both options can be
expensive for large coarse matrices.

For a code in the original translation cell, first determine a square or
oblique period, then apply `straight_coarse_grain(...)` or
`oblique_coarse_grain(...)`. The notebooks below show this complete workflow.

## Quick start: unitary decoupling decoder

`DecouplingDecoder.from_bb(...)` accepts two BB Laurent generators. It finds a
coarse cell, constructs the inverse decoupling maps, projects them to an
`Lx`-by-`Ly` finite torus, and decodes through the standard complex:

```python
from sage.all import GF, vector

from isomorphism import x, y
from unitary_decouple_decoder import DecouplingDecoder

decoder = DecouplingDecoder.from_bb(
    1 + x + x * y,
    1 + y + x * y,
    2,
    2,
)
assert decoder.check_chain_relation()

error = vector(GF(2), decoder.correction_size)
error[0] = 1
syndrome = decoder.syndrome(error)
correction = decoder.decode(syndrome, verify=True)

assert decoder.syndrome(correction) == syndrome
```

For an already coarse-grained CSS excitation matrix, construct the same public
decoder with:

```python
decoder = DecouplingDecoder.from_coarse_css(
    coarse_epsilon,
    num_x_checks=num_x_checks,
    num_qubits=num_qubits,
    Lx=Lx,
    Ly=Ly,
)
```

All finite vectors use component-major ordering. The decoder data flow is:

```text
standard_syndrome = phi0_finite * source_syndrome
standard_correction = standard_decoder.decode(standard_syndrome)
source_correction = phi1_inverse_finite * standard_correction
```

This is why syndrome transport uses the forward degree-zero map while
correction lift-back uses the inverse degree-one map.

## Worked notebooks

The detailed examples are Sage Jupyter notebooks included with this
paper-companion repository:

- [6.6.6 color code](results/isomorphism/666_color_code.ipynb)
- [4.8.8 color code](results/isomorphism/488_color_code.ipynb)
- [Bivariate bicycle code instances](results/isomorphism/bb_code_instances.ipynb)

Start Jupyter through Sage from this directory so the local packages are
available:

```bash
sage -n jupyter
```

The 6.6.6 and BB notebooks start from two Laurent generators. The 4.8.8
notebook demonstrates the full-matrix API and an explicitly chosen oblique
coarse cell.

## Tests and reproduction

Run pure-Python tests from this directory:

```bash
for test in tests/python/test_*.py; do sage -python "$test" || exit 1; done
```

Run Sage tests from each test directory so `_bootstrap.py` is importable:

```bash
cd tests/sage/core
for test in test_*.sage; do sage "$test" || exit 1; done

cd ../iso
for test in test_*.sage; do sage "$test" || exit 1; done

cd ../unitary
for test in test_*.sage; do sage "$test" || exit 1; done
```

`scripts/iso/validate_bb_instances.sage` validates the BB instances.
`scripts/unitary/` contains the baseline decoder comparison entry points.
Small generated summaries are under `results/`.

## Authorship and license

Python/SageMath code in this directory was written by Mingxin He. The release
repository is licensed under the MIT license; see the repository root
`LICENSE` and `CITATION.cff`.
