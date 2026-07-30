"""Public entrypoints for translational CSS-code decoupling algorithms."""

from .chain_maps.decoupling import solve_decoupling_unitary
from .chain_maps.result import (
    ChainIsomorphisms,
    DecouplingUnitaryResult,
    InverseChainIsomorphisms,
)
from .coarse_graining import (
    ObliqueCell,
    oblique_cell,
    coarse_grain_to_superlattice,
    coarse_grain_to_square_superlattice,
)
from .css import build_two_generator_css_excitation_map, dagger_matrix
from .periods import (
    compute_translation_representation,
    choose_minimal_area_superlattice_basis,
    PeriodData,
    find_anyon_preserving_superlattice_from_generators,
    find_anyon_preserving_superlattice_from_translation_representation,
    SuperlatticeData,
    superlattice_from_generators,
    superlattice_from_translation_action,
    TranslationAction,
)
from .rings import PRESENTATION, P, R, a, b, c, d, x, y
from .solver import LaurentSolveError, find_inverse, solve_left, solve_right

__all__ = [
    "LaurentSolveError",
    "ChainIsomorphisms",
    "DecouplingUnitaryResult",
    "InverseChainIsomorphisms",
    "ObliqueCell",
    "PRESENTATION",
    "P",
    "PeriodData",
    "R",
    "SuperlatticeData",
    "TranslationAction",
    "a",
    "b",
    "compute_translation_representation",
    "c",
    "choose_minimal_area_superlattice_basis",
    "build_two_generator_css_excitation_map",
    "dagger_matrix",
    "d",
    "solve_decoupling_unitary",
    "find_inverse",
    "oblique_cell",
    "coarse_grain_to_superlattice",
    "find_anyon_preserving_superlattice_from_generators",
    "find_anyon_preserving_superlattice_from_translation_representation",
    "solve_left",
    "solve_right",
    "coarse_grain_to_square_superlattice",
    "superlattice_from_generators",
    "superlattice_from_translation_action",
    "x",
    "y",
]
