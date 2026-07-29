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
        xi2: Degree-two inverse map accumulated during clearing.
        xi1: Degree-one inverse map accumulated during clearing.
        xi0: Degree-zero inverse map accumulated during clearing.
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
    xi2: object
    xi1: object
    xi0: object
    product_x_rank: int
    product_z_rank: int
    num_x_checks: int
    num_qubits: int
    diagnostics: dict = field(default_factory=dict)


@dataclass(frozen=True)
class DecouplingMaps:
    """Forward maps from the source chain to the standard chain.

    Attributes:
        phi2: Degree-two source-to-standard map.
        phi1: Degree-one source-to-standard map.
        phi0: Degree-zero source-to-standard map.
    """

    phi2: object
    phi1: object
    phi0: object


@dataclass(frozen=True)
class InverseDecouplingMaps:
    """Directly constructed maps from the standard chain to the source chain.

    Attributes:
        phi2_inverse: Degree-two standard-to-source map.
        phi1_inverse: Degree-one standard-to-source map.
        phi0_inverse: Degree-zero standard-to-source map.
    """

    phi2_inverse: object
    phi1_inverse: object
    phi0_inverse: object


@dataclass(frozen=True)
class DecouplingResult:
    """Standards, maps, and diagnostics from coarse-matrix decoupling.

    Attributes:
        source_coarse_matrix: Coarse CSS matrix passed to the algorithm.
        hx_standard: Standard ``H_X`` block.
        hz_standard: Standard ``H_Z`` block.
        inverse_maps: Inverse maps composed back to the input basis.
        maps: Forward maps, or ``None`` when their computation was disabled.
        diagnostics: Advanced construction and verification state.
    """

    source_coarse_matrix: object
    hx_standard: object
    hz_standard: object
    inverse_maps: InverseDecouplingMaps
    maps: object
    diagnostics: dict
