# Translational CSS-code decoupling

The `decoupling` package computes anyon-preserving superlattices,
coarse-grains Laurent CSS excitation maps, and constructs the chain
isomorphisms to the paper's standard product-state and toric-code complex.

## Chain-map convention

For an input chain

```text
degree 2 --H_Z^dagger--> degree 1 --H_X--> degree 0
```

the public maps satisfy

```text
H_X_tilde * psi_1 = psi_0 * H_X
H_Z_tilde_dagger * psi_2 = psi_1 * H_Z_dagger

H_X * psi_1_inverse = psi_0_inverse * H_X_tilde
H_Z_dagger * psi_2_inverse
    = psi_1_inverse * H_Z_tilde_dagger
```

The algorithm directly constructs `psi_inverse`. Computing `psi` requires
Laurent matrix inversion and is opt-in.

## Main API

```python
from decoupling import (
    build_two_generator_css_excitation_map,
    solve_decoupling_unitary,
    x,
    y,
)

input_excitation_map = build_two_generator_css_excitation_map(1 + x, 1 + y)
result = solve_decoupling_unitary(
    input_excitation_map,
    num_x_checks=1,
    num_qubits=2,
)

assert result.psi is None
assert result.psi_inverse.psi_1_inverse.nrows() == 2
```

`DecouplingUnitaryResult` contains:

- `input_excitation_map`;
- `h_x_tilde`;
- `h_z_tilde_dagger`;
- `psi_inverse`;
- optional `psi`;
- `diagnostics`.

`InverseChainIsomorphisms` contains `psi_2_inverse`, `psi_1_inverse`, and
`psi_0_inverse`. `ChainIsomorphisms` contains `psi_2`, `psi_1`, and `psi_0`.

To request explicit forward maps:

```python
result = solve_decoupling_unitary(
    input_excitation_map,
    num_x_checks=1,
    num_qubits=2,
    compute_psi=True,
)
```

## Construction stages

Finite-field clearing records `xi_2`, `xi_1`, and `xi_0`. The Laurent stage
constructs `phi_1`, `phi_2`, `phi_2_prime`, the homotopy `eta`, and
`phi_1_prime`. Their total inverse-direction compositions are exposed only as
`psi_2_inverse`, `psi_1_inverse`, and `psi_0_inverse`.

Inverse-only QCA validation uses the triangular `phi_2_prime` certificate and
the degree-zero `xi_2` basis change. It does not construct a Laurent forward
map.

## Shared paper examples

`decoupling.paper_examples` is the single source of truth for:

- the 6.6.6 color code;
- the 4.8.8 color code;
- all 30 BB benchmark-family rows;
- all 29 arbitrary-polynomial BB rows.

The result notebooks and reproduction scripts import these definitions rather
than duplicating polynomials or reference values.
