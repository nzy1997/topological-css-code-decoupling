"""Shared Laurent and lifted polynomial ring presentation."""

from dataclasses import dataclass

from sage.all import GF, LaurentPolynomialRing, PolynomialRing
from sage.modules.free_module_element import FreeModuleElement_generic_dense as module_element


@dataclass(frozen=True)
class LaurentPresentation:
    """Ring objects and inverse-variable mapping used by Singular lifts.

    Attributes:
        laurent_ring: Sage Laurent ring ``GF(2)[x^{+-1}, y^{+-1}]``.
        polynomial_ring: Ordinary polynomial ring used by Singular.
        x: Laurent ``x`` generator.
        y: Laurent ``y`` generator.
        a: Lifted positive ``x`` generator.
        b: Lifted positive ``y`` generator.
        c: Lifted inverse ``x^-1`` generator.
        d: Lifted inverse ``y^-1`` generator.
    """

    laurent_ring: object
    polynomial_ring: object
    x: object
    y: object
    a: object
    b: object
    c: object
    d: object

    def inverse_relations(self):
        """Return the quotient relations presenting the Laurent ring.

        Returns:
            Pair of polynomial relations ``1 + a*c`` and ``1 + b*d``. Over
            ``GF(2)`` these impose ``a*c = 1`` and ``b*d = 1``.
        """
        return (1 + self.a * self.c, 1 + self.b * self.d)


R = LaurentPolynomialRing(GF(2), "x, y")
x, y = R.gens()
P = PolynomialRing(GF(2), "a, b, c, d", 4)
a, b, c, d = P.gens()

PRESENTATION = LaurentPresentation(R, P, x, y, a, b, c, d)
