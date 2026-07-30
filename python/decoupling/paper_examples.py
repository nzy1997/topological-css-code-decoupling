"""Authoritative code instances reported in the paper's two BB-code tables.

Only structural reference data are stored here.  Runtime and chain-map degree
are implementation outputs and are therefore deliberately absent.
"""

from dataclasses import dataclass

from sage.all import Matrix, block_matrix, zero_matrix

from decoupling.css import build_two_generator_css_excitation_map
from decoupling.rings import R, x, y


Exponent = tuple[int, int]
Basis = tuple[Exponent, Exponent]


def _polynomial(exponents):
    """Build a Laurent polynomial over GF(2) from exponent pairs."""
    value = R.zero()
    for x_power, y_power in exponents:
        value += x**x_power * y**y_power
    return value


@dataclass(frozen=True)
class BBCodeInstance:
    """One complete row from a BB-code table in the paper."""

    collection: str
    row: int
    f_exponents: tuple[Exponent, ...]
    g_exponents: tuple[Exponent, ...]
    reference_square_period: int
    superlattice_basis: Basis
    reference_q: int
    reference_p_x: int
    reference_p_z: int
    reference_t: int

    @property
    def identifier(self):
        """Return the stable row identifier used by scripts and result files."""
        return f"{self.collection}-{self.row:02d}"

    def polynomials(self):
        """Return the complete Laurent pair ``(f, g)``."""
        return _polynomial(self.f_exponents), _polynomial(self.g_exponents)

    def as_reference_record(self):
        """Return JSON-serializable structural reference data."""
        f, g = self.polynomials()
        return {
            "id": self.identifier,
            "collection": self.collection,
            "row": self.row,
            "f": str(f),
            "g": str(g),
            "reference_square_period": self.reference_square_period,
            "superlattice_basis": [list(vector) for vector in self.superlattice_basis],
            "reference_q": self.reference_q,
            "reference_p_x": self.reference_p_x,
            "reference_p_z": self.reference_p_z,
            "reference_t": self.reference_t,
        }


@dataclass(frozen=True)
class ColorCodeInstance:
    """One color-code example specified explicitly in the paper."""

    label: str
    reference_square_period: int
    superlattice_basis: Basis
    reference_q: int
    reference_p_x: int
    reference_p_z: int
    reference_t: int

    def excitation_map(self):
        """Return the paper's full CSS excitation map."""
        if self.label == "666":
            return build_two_generator_css_excitation_map(
                1 + x + x * y,
                1 + y + x * y,
            )
        if self.label == "488":
            check = Matrix(
                R,
                [
                    [1 + y, 1 + y, 1 + x, 1 + x],
                    [x * y, y, x * y, x],
                ],
            )
            zero = zero_matrix(R, 2, 4)
            return block_matrix(R, [[check, zero], [zero, check]])
        raise ValueError(f"Unknown color-code label: {self.label}")


COLOR_CODE_666 = ColorCodeInstance(
    "666",
    3,
    ((3, 0), (2, 1)),
    6,
    1,
    1,
    2,
)

COLOR_CODE_488 = ColorCodeInstance(
    "488",
    2,
    ((2, 0), (1, 1)),
    8,
    2,
    2,
    2,
)


def _benchmark_case(row, a, b, period, basis, q, p, t):
    """Build a row with ``f=1+x+a`` and ``g=1+y+b``."""
    return BBCodeInstance(
        "benchmark",
        row,
        ((0, 0), (1, 0), a),
        ((0, 0), (0, 1), b),
        period,
        basis,
        q,
        p,
        p,
        t,
    )


def _table_case(row, f, g, period, basis, q, p_x, p_z, t):
    """Build a row from the paper's arbitrary-polynomial table."""
    return BBCodeInstance(
        "table",
        row,
        tuple(f),
        tuple(g),
        period,
        basis,
        q,
        p_x,
        p_z,
        t,
    )


BENCHMARK_BB_CODE_INSTANCES = (
    _benchmark_case(1, (1, 1), (1, 1), 3, ((2, 1), (1, 2)), 6, 1, 2),
    _benchmark_case(2, (-1, 1), (1, 1), 7, ((3, 1), (2, 3)), 14, 4, 3),
    _benchmark_case(3, (2, 0), (2, 0), 3, ((2, 1), (1, 2)), 6, 1, 2),
    _benchmark_case(4, (-1, 0), (0, -1), 3, ((3, 0), (0, 3)), 18, 5, 4),
    _benchmark_case(5, (1, 1), (1, -1), 7, ((3, 2), (1, 3)), 14, 4, 3),
    _benchmark_case(6, (-1, 0), (3, 2), 3, ((3, 0), (0, 3)), 18, 5, 4),
    _benchmark_case(7, (0, -2), (-2, 0), 21, ((5, 2), (2, 5)), 42, 16, 5),
    _benchmark_case(8, (0, -2), (2, 0), 15, ((3, 3), (0, 5)), 30, 11, 4),
    _benchmark_case(9, (-1, 1), (-1, -1), 31, ((5, 3), (3, 8)), 62, 26, 5),
    _benchmark_case(10, (-2, -1), (2, 1), 3, ((1, 1), (0, 3)), 6, 1, 2),
    _benchmark_case(11, (-1, 3), (3, -1), 12, ((12, 0), (0, 12)), 288, 136, 8),
    _benchmark_case(12, (-2, 0), (-2, 2), 7, ((7, 0), (0, 7)), 98, 43, 6),
    _benchmark_case(13, (-2, 1), (1, -2), 63, ((15, 6), (6, 15)), 378, 181, 8),
    _benchmark_case(14, (-1, 2), (-2, -1), 217, ((19, 10), (3, 13)), 434, 209, 8),
    _benchmark_case(15, (-3, 1), (-5, 0), 21, ((1, 1), (0, 21)), 42, 16, 5),
    _benchmark_case(16, (-2, 1), (1, 2), 105, ((9, 1), (3, 12)), 210, 98, 7),
    _benchmark_case(17, (-1, -2), (1, -1), 63, ((7, 3), (0, 9)), 126, 57, 6),
    _benchmark_case(18, (0, 2), (-4, 1), 73, ((14, 5), (5, 7)), 146, 64, 9),
    _benchmark_case(19, (-1, 2), (0, -4), 105, ((35, 7), (0, 21)), 1470, 725, 10),
    _benchmark_case(20, (0, -4), (4, 0), 255, ((15, 13), (0, 17)), 510, 239, 16),
    _benchmark_case(21, (-4, 0), (-3, 2), 21, ((21, 0), (0, 21)), 882, 431, 10),
    _benchmark_case(22, (-2, -5), (-1, -3), 42, ((6, 6), (2, 16)), 168, 77, 7),
    _benchmark_case(23, (-8, -1), (5, 1), 21, ((4, 1), (3, 6)), 42, 16, 5),
    _benchmark_case(24, (1, -5), (1, 4), 186, ((18, 42), (10, 44)), 744, 363, 9),
    _benchmark_case(25, (-1, -1), (5, 0), 105, ((15, 6), (10, 11)), 210, 98, 7),
    _benchmark_case(26, (1, 3), (2, -2), 63, ((35, 7), (7, 14)), 882, 432, 9),
    _benchmark_case(27, (-1, -2), (2, -1), 217, ((13, 3), (10, 19)), 434, 209, 8),
    _benchmark_case(28, (-1, 3), (1, 3), 186, ((16, 10), (14, 32)), 744, 363, 9),
    _benchmark_case(29, (2, 2), (-4, 1), 889, ((101, 5), (4, 9)), 1778, 879, 10),
    _benchmark_case(30, (-1, 3), (3, 0), 217, ((17, 9), (8, 17)), 434, 209, 8),
)


TABLE_BB_CODE_INSTANCES = (
    _table_case(1, ((1, 0), (0, 0)), ((1, 0), (0, 1)), 1, ((1, 0), (0, 1)), 2, 0, 0, 1),
    _table_case(2, ((1, 0), (0, 2)), ((0, 2), (0, 1)), 1, ((1, 0), (0, 1)), 2, 0, 0, 1),
    _table_case(3, ((2, 1), (1, 1), (1, 0)), ((1, 1), (0, 2)), 3, ((2, 1), (1, 2)), 6, 1, 1, 2),
    _table_case(4, ((2, 0), (1, 0), (0, 0)), ((1, 0), (0, 1)), 3, ((2, 1), (1, 2)), 6, 1, 1, 2),
    _table_case(5, ((2, 1), (1, 1), (0, 0)), ((1, 1), (1, 0)), 3, ((3, 0), (0, 1)), 6, 1, 1, 2),
    _table_case(6, ((2, 2), (1, 0), (0, 0)), ((1, 2), (1, 1)), 3, ((3, 0), (0, 1)), 6, 1, 1, 2),
    _table_case(7, ((3, 1), (3, 0), (2, 1), (1, 1)), ((3, 1), (0, 2)), 4, ((1, 1), (0, 4)), 8, 1, 1, 3),
    _table_case(8, ((2, 3), (1, 0), (0, 0)), ((1, 3), (1, 2), (1, 1)), 3, ((3, 0), (0, 3)), 18, 5, 5, 4),
    _table_case(9, ((1, 0), (0, 1)), ((0, 4), (0, 3), (0, 2), (0, 1)), 4, ((2, 2), (1, 3)), 8, 1, 1, 3),
    _table_case(10, ((0, 1), (-1, 2), (-2, 2), (-2, 1)), ((1, 1), (0, 2), (0, 0), (-2, 0)), 3, ((1, 1), (0, 3)), 6, 0, 0, 3),
    _table_case(11, ((1, 2), (1, 0), (-1, 3), (-2, 4)), ((1, 3), (1, 1), (0, 2), (0, 0)), 40, ((2, 2), (0, 40)), 160, 71, 71, 9),
    _table_case(12, ((3, -2), (2, -1), (0, 2), (0, 0)), ((2, 1), (2, -1), (1, 0), (1, -2)), 40, ((2, 2), (0, 40)), 160, 71, 71, 9),
    _table_case(13, ((3, 1), (0, 0), (-1, 3), (-2, 4)), ((3, 3), (1, 1), (0, 2), (-2, 0)), 24, ((4, 4), (0, 24)), 192, 76, 76, 20),
    _table_case(14, ((3, 0), (1, 0), (-1, 3), (-2, 4)), ((3, 1), (2, 0), (1, 1), (0, 0)), 56, ((2, 2), (0, 56)), 224, 101, 101, 11),
    _table_case(15, ((3, 0), (0, 2), (0, 1)), ((2, 0), (1, 0), (0, 3)), 12, ((12, 0), (0, 12)), 288, 136, 136, 8),
    _table_case(16, ((27, 0), (6, 0), (0, 1)), ((24, 0), (15, 0), (0, 0)), 93, ((93, 0), (0, 1)), 186, 78, 78, 15),
    _table_case(17, ((9, 0), (0, 2), (0, 1)), ((8, 0), (1, 0), (0, 0)), 63, ((36, 27), (27, 36)), 1134, 551, 551, 16),
    _table_case(18, ((3, 0), (0, 2), (0, 1)), ((2, 0), (1, 0), (0, 0)), 3, ((3, 0), (0, 3)), 18, 5, 5, 4),
    _table_case(19, ((3, 1), (2, 2), (0, 0)), ((3, 2), (2, 0), (0, 0)), 105, ((21, 3), (0, 5)), 210, 98, 98, 7),
    _table_case(20, ((3, 2), (2, 3), (0, 0)), ((3, 1), (2, 0), (0, 0)), 217, ((14, 23), (7, 27)), 434, 209, 209, 8),
    _table_case(21, ((3, 2), (2, 0), (0, 0)), ((3, 0), (2, 1), (0, 0)), 7, ((7, 0), (0, 1)), 14, 4, 4, 3),
    _table_case(22, ((3, 1), (2, 0), (0, 0)), ((3, 2), (2, 0), (0, 0)), 7, ((7, 0), (0, 1)), 14, 4, 4, 3),
    _table_case(23, ((3, 0), (2, 0), (0, 1)), ((3, 0), (2, 0), (0, 0)), 7, ((7, 0), (0, 1)), 14, 4, 4, 3),
    _table_case(24, ((6, 2), (3, 0), (1, 3), (0, 0)), ((6, 1), (5, 1), (4, 1), (3, 2)), 120, ((8, 6), (0, 30)), 480, 230, 230, 10),
    _table_case(25, ((6, 0), (3, 0), (1, 1), (0, 0)), ((6, 1), (5, 1), (4, 1), (3, 0)), 8, ((8, 0), (0, 4)), 64, 25, 25, 7),
    _table_case(26, ((5, 1), (5, 0), (3, 0), (0, 0)), ((5, 1), (3, 0), (2, 1), (0, 1)), 186, ((18, 2), (12, 22)), 744, 362, 362, 10),
    _table_case(27, ((5, 3), (4, 4), (3, 5), (0, 5)), ((5, 1), (4, 3), (3, 3), (0, 2)), 56, ((28, 2), (0, 4)), 224, 103, 103, 9),
    _table_case(28, ((5, 1), (4, 0), (3, 1), (0, 1)), ((5, 1), (4, 1), (3, 1), (0, 0)), 56, ((28, 2), (0, 4)), 224, 103, 103, 9),
    _table_case(29, ((3, 0), (0, 7), (0, 2)), ((2, 0), (1, 0), (0, 3)), 126, ((12, 6), (6, 66)), 1512, 744, 744, 12),
)


BB_CODE_INSTANCES = BENCHMARK_BB_CODE_INSTANCES + TABLE_BB_CODE_INSTANCES


def select_bb_code_instances(selectors=None):
    """Select rows by stable identifier, failing on unknown identifiers."""
    if not selectors:
        return BB_CODE_INSTANCES
    requested = set(selectors)
    by_identifier = {case.identifier: case for case in BB_CODE_INSTANCES}
    unknown = sorted(requested - set(by_identifier))
    if unknown:
        raise ValueError(f"Unknown BB-code row identifiers: {', '.join(unknown)}")
    return tuple(case for case in BB_CODE_INSTANCES if case.identifier in requested)


__all__ = [
    "BBCodeInstance",
    "BB_CODE_INSTANCES",
    "BENCHMARK_BB_CODE_INSTANCES",
    "COLOR_CODE_488",
    "COLOR_CODE_666",
    "ColorCodeInstance",
    "TABLE_BB_CODE_INSTANCES",
    "select_bb_code_instances",
]
