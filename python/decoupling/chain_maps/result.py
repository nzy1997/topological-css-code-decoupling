"""Result objects for the decoupling stages."""

from dataclasses import dataclass, field


@dataclass(frozen=True)
class ClearingResult:
    """Finite-field CSS clearing with accumulated row and column maps.

    Attributes:
        input_matrix: Original coarse CSS excitation matrix.
        working_matrix: Cleared Laurent matrix after mod-J pivot choices.
        row_map: Block-diagonal row redefinition accumulated during clearing.
        symplectic_column_map: Full CSS symplectic column map.
        xi_2: Degree-two inverse map accumulated during clearing.
        xi_1: Degree-one inverse map accumulated during clearing.
        xi_0: Degree-zero inverse map accumulated during clearing.
        product_x_rank: Number of X product rows found by clearing.
        product_z_rank: Number of Z product rows found by clearing.
        num_x_checks: Number of X-check rows.
        num_qubits: Number of data-qubit columns per CSS half.
        diagnostics: Optional algorithm diagnostics.
    """

    input_matrix: object
    working_matrix: object
    row_map: object
    symplectic_column_map: object
    xi_2: object
    xi_1: object
    xi_0: object
    product_x_rank: int
    product_z_rank: int
    num_x_checks: int
    num_qubits: int
    diagnostics: dict = field(default_factory=dict)


@dataclass(frozen=True)
class ChainIsomorphisms:
    """The theorem's chain isomorphisms from input to standard form.

    Attributes:
        psi_2: Degree-two input-to-standard isomorphism.
        psi_1: Degree-one input-to-standard isomorphism.
        psi_0: Degree-zero input-to-standard isomorphism.
    """

    psi_2: object
    psi_1: object
    psi_0: object


@dataclass(frozen=True)
class InverseChainIsomorphisms:
    """The inverse chain isomorphisms constructed directly by the algorithm.

    Attributes:
        psi_2_inverse: Degree-two standard-to-input isomorphism.
        psi_1_inverse: Degree-one standard-to-input isomorphism.
        psi_0_inverse: Degree-zero standard-to-input isomorphism.
    """

    psi_2_inverse: object
    psi_1_inverse: object
    psi_0_inverse: object


@dataclass(frozen=True)
class DecouplingUnitaryResult:
    """Paper-aligned output of the decoupling-unitary algorithm.

    Attributes:
        input_excitation_map: Coarse CSS excitation map passed to the algorithm.
        h_x_tilde: Standard-complex ``\\tilde H_X`` matrix.
        h_z_tilde_dagger: Standard-complex ``\\tilde H_Z^\\dagger`` matrix.
        psi_inverse: Directly constructed inverse chain isomorphisms.
        psi: Forward chain isomorphisms, or ``None`` when not requested.
        diagnostics: Advanced construction and verification state.
    """

    input_excitation_map: object
    h_x_tilde: object
    h_z_tilde_dagger: object
    psi_inverse: InverseChainIsomorphisms
    psi: object
    diagnostics: dict
