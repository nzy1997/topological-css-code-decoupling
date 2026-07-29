"""Shared BB Table I case inventory."""

from dataclasses import dataclass


@dataclass(frozen=True)
class BBRow:
    """One BB Table I row represented by monomial exponent pairs."""

    row: int
    a_exponents: tuple[int, int]
    b_exponents: tuple[int, int]
    paper_period: int


BB_ROWS = (
    BBRow(1, (1, 1), (1, 1), 3),
    BBRow(2, (-1, 1), (1, 1), 7),
    BBRow(3, (2, 0), (2, 0), 3),
    BBRow(4, (-1, 0), (0, -1), 3),
    BBRow(5, (1, 1), (1, -1), 7),
    BBRow(6, (-1, 0), (3, 2), 3),
    BBRow(7, (0, -2), (-2, 0), 21),
    BBRow(8, (0, -2), (2, 0), 15),
    BBRow(9, (-1, 1), (-1, -1), 31),
    BBRow(10, (-2, -1), (2, 1), 3),
    BBRow(11, (-1, 3), (3, -1), 12),
    BBRow(12, (-2, 0), (-2, 2), 7),
    BBRow(13, (-2, 1), (1, -2), 63),
    BBRow(14, (-1, 2), (-2, -1), 217),
    BBRow(15, (-3, 1), (-5, 0), 21),
    BBRow(16, (-2, 1), (1, 2), 105),
    BBRow(17, (-1, -2), (1, -1), 63),
    BBRow(18, (0, 2), (-4, 1), 73),
    BBRow(19, (0, -4), (4, 0), 255),
    BBRow(20, (-4, 0), (-3, 2), 21),
    BBRow(21, (-2, -5), (-1, -3), 42),
    BBRow(22, (-8, -1), (5, 1), 21),
    BBRow(23, (1, -5), (1, 4), 186),
    BBRow(24, (-1, -1), (5, 0), 105),
    BBRow(25, (1, 3), (2, -2), 63),
    BBRow(26, (-1, -2), (2, -1), 217),
    BBRow(27, (-1, 3), (1, 3), 186),
    BBRow(28, (-1, 3), (3, 0), 217),
)

BB_ROWS_LTE_63 = tuple(row for row in BB_ROWS if row.paper_period <= 63)
