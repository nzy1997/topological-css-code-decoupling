# Translational CSS-code isomorphism

The `isomorphism` package finds translation periods, coarse-grains Laurent CSS
matrices, and constructs a chain isomorphism from a coarse source code to its
standard product-state and toric-code decomposition.

## Convention

The package follows the supplement convention. The supplement's theorem draws
the maps from the source chain to the standard chain even though its algorithm
constructs the inverse direction before returning the forward maps.

For a source chain

```text
source degree 2 --hz_source_dagger--> source degree 1 --hx_source--> source degree 0
```

and its standard chain, the public forward maps satisfy

```text
hx_standard * phi1 == phi0 * hx_source
hz_standard_dagger * phi2 == phi1 * hz_source_dagger
```

The Laurent solvers construct the inverse direction first:

```text
hx_source * phi1_inverse == phi0_inverse * hx_standard
hz_source_dagger * phi2_inverse == phi1_inverse * hz_standard_dagger
```

This explicit naming prevents syndrome transport from being confused with
correction lift-back. A source syndrome is transported by `phi0`; a standard
correction is lifted back by `phi1_inverse`.

## Main API

```python
from isomorphism import (
    DecouplingMaps,
    DecouplingResult,
    InverseDecouplingMaps,
    construct_excitation_map,
    decouple_coarse_matrix,
    x,
    y,
)

epsilon = construct_excitation_map(1 + x, 1 + y)
result = decouple_coarse_matrix(
    epsilon,
    num_x_checks=1,
    num_qubits=2,
)
```

`DecouplingResult` contains:

- `source_coarse_matrix`: the source excitation matrix;
- `hx_standard` and `hz_standard`: the standard CSS blocks;
- `inverse_maps`: inverse maps composed into the source basis;
- `maps`: forward maps, or `None` when forward-map construction is disabled;
- `diagnostics`: advanced clearing, rank, and inverse-only QCA state.

`InverseDecouplingMaps` contains `phi2_inverse`, `phi1_inverse`,
and `phi0_inverse`.

`DecouplingMaps` contains only `phi2`, `phi1`, and `phi0`.

When only the directly constructed direction is needed, skip the two Laurent
inversions:

```python
result = decouple_coarse_matrix(
    epsilon,
    num_x_checks=1,
    num_qubits=2,
    compute_forward_maps=False,
)
assert result.maps is None
```

Inverse maps are always present.

## Clearing maps

`clear_css_mod_j(...)` records three inverse-direction basis maps:

- `xi2` on degree two;
- `xi1` on degree one;
- `xi0` on degree zero.

It also retains the complete row and symplectic column transformations for
debugging the clearing stage.

## QCA verification

The verification helpers make both directions explicit:

```python
from isomorphism.chain_maps.decoupling import (
    qca_decoupled_excitation_matrix,
    qca_symplectic_matrix,
    verify_qca_decoupling,
)

symplectic = qca_symplectic_matrix(result.maps, result.inverse_maps)
decoupled = qca_decoupled_excitation_matrix(
    result.source_coarse_matrix,
    result.maps,
    result.inverse_maps,
)
assert verify_qca_decoupling(result)
```

Inverse-only verification uses the post-elimination triangular
`phi2_inverse` certificate, then transports its Z-side product through the
clearing map. This requires inverting the degree-zero binary basis change
`xi2`, but it does not construct any Laurent forward decoupling map.

## Public result metadata

The package root exports the reusable search and geometry result types:

- `PeriodData`;
- `SuperlatticeData`;
- `TranslationAction`;
- `ObliqueCell`.

It also exports `DecouplingMaps`, `InverseDecouplingMaps`, and
`DecouplingResult`.

## Verification

Run tests from their directory so `_bootstrap.py` is importable:

```bash
cd tests/sage/iso
for test in test_*.sage; do
  sage "$test" || exit 1
done
```
