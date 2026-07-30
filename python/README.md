# Decoupling 2D topological CSS codes

This directory contains the SageMath implementation of the decoupling
algorithm for two-dimensional translational topological CSS codes and the
unitary-decouple-based decoder built from its chain isomorphisms.

The public packages are:

- `decoupling`: translation representations, anyon-preserving superlattices,
  coarse graining, standard complexes, and chain isomorphisms;
- `unitary_decouple_based_decoder`: finite-torus decoding through the standard
  product-state and toric-code complex;
- `decoder_core`: shared finite-torus, matching, sampling, and benchmark
  utilities.

## Requirements and installation

Use Python 3.10 or newer inside a SageMath environment:

```bash
sage -pip install -e .
```

NumPy and PyMatching are base dependencies. The optional `benchmark` extra
adds `ldpc` and Matplotlib:

```bash
sage -pip install -e '.[benchmark]'
```

## Decoupling-unitary API

For an input CSS chain

```text
degree 2 --H_Z^dagger--> degree 1 --H_X--> degree 0
```

the paper's forward chain isomorphisms obey

```text
H_X_tilde * psi_1 = psi_0 * H_X
H_Z_tilde_dagger * psi_2 = psi_1 * H_Z_dagger
```

The algorithm directly constructs the inverse direction:

```text
H_X * psi_1_inverse = psi_0_inverse * H_X_tilde
H_Z_dagger * psi_2_inverse
    = psi_1_inverse * H_Z_tilde_dagger
```

```python
from decoupling import (
    build_two_generator_css_excitation_map,
    dagger_matrix,
    solve_decoupling_unitary,
    x,
    y,
)
from decoupling.css import split_css_blocks

input_excitation_map = build_two_generator_css_excitation_map(
    1 + x,
    1 + y,
)
result = solve_decoupling_unitary(
    input_excitation_map,
    num_x_checks=1,
    num_qubits=2,
    compute_psi=True,
)

h_z, h_x = split_css_blocks(input_excitation_map, 1, 2)
assert result.h_x_tilde * result.psi.psi_1 == result.psi.psi_0 * h_x
assert (
    result.h_z_tilde_dagger * result.psi.psi_2
    == result.psi.psi_1 * dagger_matrix(h_z)
)
assert (
    h_x * result.psi_inverse.psi_1_inverse
    == result.psi_inverse.psi_0_inverse * result.h_x_tilde
)
assert (
    dagger_matrix(h_z) * result.psi_inverse.psi_2_inverse
    == result.psi_inverse.psi_1_inverse * result.h_z_tilde_dagger
)
```

`DecouplingUnitaryResult` exposes exactly:

```text
input_excitation_map
h_x_tilde
h_z_tilde_dagger
psi_inverse
psi
diagnostics
```

`psi_inverse` always contains `psi_2_inverse`, `psi_1_inverse`, and
`psi_0_inverse`. The default `compute_psi=False` avoids full Laurent matrix
inversions. Set `compute_psi=True` only when `psi_2`, `psi_1`, and `psi_0` are
needed explicitly.

The internal construction follows the paper's stages:
`xi_2`, `xi_1`, `xi_0`, then `phi_1`, `phi_2`, `phi_2_prime`, `eta`, and
`phi_1_prime`. Only their total compositions are called `psi_inverse`.

## Unitary-decouple-based decoder

```python
from sage.all import GF, vector

from decoupling import x, y
from unitary_decouple_based_decoder import UnitaryDecoupleBasedDecoder

decoder = UnitaryDecoupleBasedDecoder.from_bb(
    1 + x + x * y,
    1 + y + x * y,
    2,
    2,
)
error = vector(GF(2), decoder.correction_size)
error[0] = 1
syndrome = decoder.syndrome(error)
correction = decoder.decode(syndrome, verify=True)
assert decoder.syndrome(correction) == syndrome
```

Finite vectors use component-major ordering. Syndrome transport uses
`psi_0_finite`; correction lift-back uses `psi_1_inverse_finite`.

## Paper reproductions

The three decoupling runners use the paper-specified superlattice bases and
fail with a nonzero exit status if any structural, chain, inverse, or
symplectic identity fails:

```bash
sage scripts/decoupling/reproduce_666_color_code.sage
sage scripts/decoupling/reproduce_488_color_code.sage
sage scripts/decoupling/reproduce_bb_codes.sage
```

The shared BB catalog contains all 30 rows from the
`f=1+x+a, g=1+y+b` table and all 29 rows from the arbitrary-polynomial table.
The default publication run validates the 56 rows whose paper runtimes are
shorter than `benchmark-19`, plus `benchmark-19` itself as the cutoff
reference. Only the two longer rows, `benchmark-29` and `table-17`, are skipped
by default; either remains available through an explicit row selector. The
runner supports stable selectors, atomic checkpoints, resume, retries, and
optional positive timeouts:

```bash
sage scripts/decoupling/reproduce_bb_codes.sage \
  --rows benchmark-01,table-03

sage scripts/decoupling/reproduce_bb_codes.sage \
  --resume --retry-failed --retries 1
```

Selected-row runs use `bb_codes.selection.*` filenames and never overwrite the
canonical 57-row cutoff-qualified result. In-progress resume data are kept in
a private atomic checkpoint; canonical JSON, CSV, and Markdown files are
published only after every selected row passes.

Canonical JSON, CSV, and Markdown outputs are written under
`results/decoupling/`. They record actual local runtime and
`deg(psi_1_inverse)`; paper timing is intentionally omitted. Chain
isomorphisms are non-unique, so validation uses mathematical identities rather
than entry-by-entry agreement with printed matrices.

The explanatory notebooks consume the complete shared catalog and canonical
validated-subset results:

- [6.6.6 color code](results/decoupling/666_color_code.ipynb)
- [4.8.8 color code](results/decoupling/488_color_code.ipynb)
- [BB-code catalog and validated subset](results/decoupling/bb_code_instances.ipynb)

The decoder comparison uses
`H_X=(1+x+x^-1*y, 1+y+x*y)`, distances 4 and 6, error probabilities
0.01, 0.03, 0.05, 0.07, and 0.10, 100,000 shots per point, and seed 20260607:

```bash
sage scripts/unitary_decouple_based_decoder/benchmark_paper_bb_family.sage
```

Its CSV, Markdown, and plot are written under
`results/unitary_decouple_based_decoder/`.

## Authorship and license

The Python/SageMath implementation was written by Mingxin He. This public
repository is distributed under the MIT license; see the repository root
`LICENSE` and `CITATION.cff`.
