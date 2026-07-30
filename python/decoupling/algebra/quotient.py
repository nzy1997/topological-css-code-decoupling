"""Finite quotient-module bases for anyon-preserving translations."""

from dataclasses import dataclass, field

from sage.all import FreeModule, GF, vector

from decoupling.algebra.singular_backend import (
    laurent_to_polynomial,
    polynomial_to_laurent,
    quotient_dimension,
    quotient_monomial_basis,
    reduce_by_standard_basis,
    standard_basis,
)
from decoupling.rings import P, R, PRESENTATION, module_element


@dataclass(frozen=True)
class QuotientModuleBasis:
    """Standard basis plus Singular monomial basis for a finite quotient.

    Attributes:
        module: Free module over the lifted polynomial ring.
        standard_basis: Singular standard basis for the quotient relations.
        monomial_basis: Singular monomial basis of the finite quotient.
        diagnostics: Optional size and dimension diagnostics.
    """

    module: object
    standard_basis: object
    monomial_basis: tuple
    diagnostics: dict = field(default_factory=dict)

    def normal_form(self, element):
        """Reduce one module element to quotient normal form.

        Args:
            element: Module element over the lifted polynomial ring.

        Returns:
            Normal form modulo ``standard_basis``.
        """
        return reduce_by_standard_basis(element, self.standard_basis)


def _monomial_key(module_element):
    """Use the one-term module vector itself as the basis lookup key.

    Args:
        module_element: Singular/Sage module monomial.

    Returns:
        Immutable module element usable as a dictionary key.
    """
    if hasattr(module_element, "set_immutable"):
        module_element.set_immutable()
    return module_element


def _laurent_module_entry(module, module_monomial):
    """Project one module monomial to a Laurent row-entry tuple.

    Args:
        module: Free module containing ``module_monomial``.
        module_monomial: One-term module basis monomial.

    Returns:
        Tuple of Laurent entries representing the monomial in module-row
        coordinates.
    """
    entries = [R.zero()] * len(module.basis())
    entries[module_monomial.leading_support()] = polynomial_to_laurent(
        module_monomial.leading_coefficient()
    )
    return tuple(entries)


@dataclass(frozen=True)
class QuotientCoordinateSystem:
    """Coordinate vectors for a finite Laurent quotient module.

    Attributes:
        quotient: Underlying quotient basis object.
        laurent_basis: Laurent row-entry representatives for each coordinate.
        _index: Lookup from monomial normal-form term to coordinate index.
    """

    quotient: QuotientModuleBasis
    laurent_basis: tuple
    _index: dict

    @property
    def dimension(self):
        """Return the quotient dimension over GF(2).

        Returns:
            Number of quotient basis monomials.
        """
        return len(self.laurent_basis)

    def normal_form(self, entries):
        """Reduce a Laurent module element represented by row entries.

        Args:
            entries: Iterable of Laurent entries, one per module row.

        Returns:
            Lifted polynomial normal form in the quotient.

        Raises:
            ValueError: If the number of entries does not match the module rank.
        """
        raw = tuple(entries)
        if len(raw) != len(self.quotient.module.basis()):
            raise ValueError("Entry count must match the quotient module rank.")
        # Lift Laurent row entries to the polynomial presentation before
        # reducing by the quotient standard basis.
        lifted = module_element(
            self.quotient.module,
            [laurent_to_polynomial(R(entry)) for entry in raw],
        )
        return self.quotient.normal_form(lifted)

    def coordinates(self, entries):
        """Return GF(2) coordinates of one Laurent module element.

        Args:
            entries: Iterable of Laurent entries, one per module row.

        Returns:
            GF(2) vector in the quotient monomial basis.

        Raises:
            ValueError: If the reduced normal form contains a term outside the
                cached quotient basis.
        """
        coords = vector(GF(2), self.dimension)
        reduced = self.normal_form(entries)
        for support, coeff in reduced.monomial_coefficients().items():
            while not coeff.is_zero():
                # Singular groups terms by module support; peel leading
                # monomials one at a time to fill binary coordinates.
                term = coeff.leading_monomial() * self.quotient.module.basis()[support]
                key = _monomial_key(term)
                if key not in self._index:
                    raise ValueError("Normal form uses a monomial outside the quotient basis.")
                coords[self._index[key]] += 1
                coeff -= coeff.leading_monomial()
        return coords


def _lift_column_module(mat):
    """Return ``P``-module generators for ``im(mat)`` plus inverse relations.

    Args:
        mat: Laurent matrix whose columns generate a submodule.

    Returns:
        Pair ``(module, generators)`` in the lifted polynomial presentation.
    """
    module = FreeModule(P, mat.nrows())
    generators = []
    for col in range(mat.ncols()):
        # Each Laurent matrix column becomes one module generator over P.
        generators.append(
            module_element(
                module,
                [laurent_to_polynomial(entry) for entry in mat.column(col)],
            )
        )
    x_relation, y_relation = PRESENTATION.inverse_relations()
    for basis_element in module.basis():
        # Add inverse relations in every module component so normal forms are
        # computed modulo the Laurent presentation.
        generators.append(x_relation * basis_element)
        generators.append(y_relation * basis_element)
    return module, generators


def quotient_module_basis(mat, *, diagnostics=False):
    """Return the finite quotient basis for ``P^m / lifted im(mat)``.

    Args:
        mat: Laurent matrix defining a cokernel quotient.
        diagnostics: If true, include generator and basis sizes.

    Returns:
        ``QuotientModuleBasis`` object.

    Raises:
        ValueError: If the quotient is not finite-dimensional.
    """
    module, generators = _lift_column_module(mat)
    groebner_basis = standard_basis(generators)
    dimension = quotient_dimension(groebner_basis)
    if dimension < 0:
        raise ValueError("Translation quotient is not finite-dimensional.")
    monomial_basis = quotient_monomial_basis(groebner_basis)
    return QuotientModuleBasis(
        module,
        groebner_basis,
        monomial_basis,
        {
            "generator_count": len(generators),
            "standard_basis_size": len(groebner_basis),
            "quotient_dimension": dimension,
            "monomial_basis_size": len(monomial_basis),
        }
        if diagnostics
        else {},
    )


def quotient_coordinate_system(mat, *, diagnostics=False):
    """Return finite quotient coordinates for the cokernel of ``mat``.

    Args:
        mat: Laurent matrix defining a finite cokernel quotient.
        diagnostics: If true, include quotient construction diagnostics.

    Returns:
        ``QuotientCoordinateSystem`` for normal-form coordinate conversion.
    """
    quotient = quotient_module_basis(mat, diagnostics=diagnostics)
    index = {
        _monomial_key(basis_element): row
        for row, basis_element in enumerate(quotient.monomial_basis)
    }
    # Store Laurent representatives so callers can inspect or reconstruct the
    # quotient basis without dealing with Singular module monomials.
    laurent_basis = tuple(
        _laurent_module_entry(quotient.module, basis_element)
        for basis_element in quotient.monomial_basis
    )
    return QuotientCoordinateSystem(quotient, laurent_basis, index)
