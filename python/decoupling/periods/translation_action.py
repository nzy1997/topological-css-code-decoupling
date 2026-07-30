"""Quotient representations of lattice translations."""

from dataclasses import dataclass, field

from sage.all import GF, Matrix

from decoupling.algebra.quotient import quotient_module_basis
from decoupling.rings import a, b


@dataclass(frozen=True)
class TranslationAction:
    """Matrices for multiplication by x and y on a quotient basis.

    Attributes:
        tx: Binary matrix for multiplication by lifted ``x``.
        ty: Binary matrix for multiplication by lifted ``y``.
        basis: Singular monomial basis for the finite quotient.
        diagnostics: Optional quotient-construction diagnostics.
    """

    tx: object
    ty: object
    basis: tuple
    diagnostics: dict = field(default_factory=dict)


def _monomial_key(module_element):
    """Use the one-term module vector itself as the basis lookup key.

    Args:
        module_element: Singular/Sage module monomial.

    Returns:
        Immutable module element suitable as a dictionary key.
    """
    if hasattr(module_element, "set_immutable"):
        module_element.set_immutable()
    return module_element


def _action_matrix(quotient, variable):
    """Reduce multiplication by one lifted translation variable.

    Args:
        quotient: Finite quotient module basis.
        variable: Lifted polynomial variable ``a`` or ``b``.

    Returns:
        GF(2) matrix for multiplication by ``variable`` on the quotient basis.

    Raises:
        ValueError: If a normal form references a monomial outside the basis.
    """
    index = {
        _monomial_key(basis_element): row
        for row, basis_element in enumerate(quotient.monomial_basis)
    }
    action = Matrix(GF(2), len(index), len(index), 0)
    for col, basis_element in enumerate(quotient.monomial_basis):
        # Multiply one basis monomial, reduce modulo the quotient, and expand
        # the reduced normal form in the same monomial basis.
        reduced = quotient.normal_form(variable * basis_element)
        for support, coeff in reduced.monomial_coefficients().items():
            while not coeff.is_zero():
                term = coeff.leading_monomial() * quotient.module.basis()[support]
                key = _monomial_key(term)
                if key not in index:
                    raise ValueError("Normal form uses a monomial outside the quotient basis.")
                action[index[key], col] = 1
                coeff -= coeff.leading_monomial()
    return action


def compute_translation_representation(mat, *, diagnostics=False):
    """Build ``T_x`` and ``T_y`` for the finite quotient associated with ``mat``.

    Args:
        mat: Laurent matrix whose cokernel quotient carries translation action.
        diagnostics: If true, include quotient basis diagnostics.

    Returns:
        ``TranslationAction`` containing action matrices and quotient basis.
    """
    quotient = quotient_module_basis(mat, diagnostics=diagnostics)
    return TranslationAction(
        _action_matrix(quotient, a),
        _action_matrix(quotient, b),
        quotient.monomial_basis,
        quotient.diagnostics,
    )
