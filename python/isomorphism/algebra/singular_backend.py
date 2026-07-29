"""Narrow wrappers around Singular primitives used by isomorphism."""

from copy import copy

from sage.all import Sequence
from sage.libs.singular.function_factory import ff

from isomorphism.rings import P, R, a, b, c, d, x, y

_kbase = ff.kbase
_lift = ff.lift
_reduce = ff.reduce
_std = ff.std
_vdim = ff.vdim


class SingularBackendError(RuntimeError):
    """Wrap Singular failures with the algorithm context that triggered them."""


def laurent_to_polynomial(poly):
    """Lift one Laurent polynomial to ``F_2[a,b,c,d]``.

    Args:
        poly: Laurent polynomial over ``R``.

    Returns:
        Polynomial where positive exponents use ``a,b`` and negative
        exponents use ``c,d``.
    """
    lifted = P.zero()
    for (px, py), coeff in poly.dict().items():
        # Ordinary polynomial rings cannot store negative powers, so x^-k and
        # y^-k are represented by independent inverse variables c^k and d^k.
        lifted += coeff * (a**px if px >= 0 else c**(-px)) * (b**py if py >= 0 else d**(-py))
    return lifted


def polynomial_to_laurent(poly):
    """Project one lifted polynomial back to the Laurent ring.

    Args:
        poly: Polynomial over ``P``.

    Returns:
        Laurent polynomial after identifying ``a=x``, ``b=y``, ``c=x^-1``,
        and ``d=y^-1``.
    """
    projected = R.zero()
    for (pa, pb, pc, pd), coeff in poly.dict().items():
        projected += coeff * x ** (pa - pc) * y ** (pb - pd)
    return projected


def matrix_laurent_to_polynomial(mat):
    """Lift every entry of a Laurent matrix.

    Args:
        mat: Sage matrix over ``R``.

    Returns:
        Matrix of the same shape over ``P``.
    """
    return mat.apply_map(laurent_to_polynomial)


def matrix_polynomial_to_laurent(mat):
    """Project every entry of a lifted matrix.

    Args:
        mat: Sage matrix over ``P``.

    Returns:
        Matrix of the same shape over ``R``.
    """
    return mat.apply_map(polynomial_to_laurent)


def _immutable_sequence(elements):
    """Make standard-basis generators hashable for Singular attributes.

    Args:
        elements: Iterable of Singular/Sage module elements.

    Returns:
        Immutable Sage ``Sequence``.
    """
    frozen = []
    for element in elements:
        if hasattr(element, "set_immutable"):
            element.set_immutable()
        frozen.append(element)
    return Sequence(frozen, immutable=True)


def standard_basis(generators):
    """Return an immutable Singular standard basis for generators.

    Args:
        generators: Module generators over ``P``.

    Returns:
        Immutable standard basis sequence.

    Raises:
        SingularBackendError: If Singular ``std`` fails.
    """
    try:
        return _immutable_sequence(_std(Sequence(generators)))
    except RuntimeError as exc:
        raise SingularBackendError("Singular std failed for module generators.") from exc


def _standard_attributes(groebner_basis):
    """Tell Singular that the immutable Sage sequence is already standard.

    Args:
        groebner_basis: Immutable standard basis sequence.

    Returns:
        Attribute dictionary accepted by Sage's Singular wrappers.
    """
    return {groebner_basis: {"isSB": 1}}


def reduce_by_standard_basis(element, groebner_basis):
    """Compute the normal form of one element by a standard basis.

    Args:
        element: Polynomial/module element to reduce.
        groebner_basis: Immutable Singular standard basis.

    Returns:
        Reduced normal form.

    Raises:
        SingularBackendError: If Singular ``reduce`` fails.
    """
    reduced_element = copy(element)
    if hasattr(reduced_element, "set_immutable"):
        reduced_element.set_immutable()
    try:
        # Passing isSB avoids recomputing the standard basis inside Singular.
        return _reduce(reduced_element, groebner_basis, attributes=_standard_attributes(groebner_basis))
    except RuntimeError as exc:
        raise SingularBackendError("Singular reduce failed for a standard basis.") from exc


def quotient_dimension(groebner_basis):
    """Return Singular ``vdim`` for a quotient by a standard basis.

    Args:
        groebner_basis: Immutable Singular standard basis.

    Returns:
        Integer quotient dimension. Singular returns negative values for
        infinite-dimensional quotients.

    Raises:
        SingularBackendError: If Singular ``vdim`` fails.
    """
    try:
        return int(_vdim(groebner_basis, attributes=_standard_attributes(groebner_basis)))
    except RuntimeError as exc:
        raise SingularBackendError("Singular vdim failed for a standard basis.") from exc


def quotient_monomial_basis(groebner_basis):
    """Return Singular ``kbase`` monomials for a finite quotient module.

    Args:
        groebner_basis: Immutable Singular standard basis.

    Returns:
        Tuple of monomial basis elements.

    Raises:
        ValueError: If the quotient is infinite-dimensional.
        SingularBackendError: If Singular ``kbase`` fails.
    """
    dimension = quotient_dimension(groebner_basis)
    if dimension < 0:
        raise ValueError("The quotient module is not finite-dimensional.")
    try:
        basis = _kbase(groebner_basis, attributes=_standard_attributes(groebner_basis))
    except RuntimeError as exc:
        raise SingularBackendError("Singular kbase failed for a finite quotient.") from exc
    return tuple(basis)


def _validate_lift_matrix(mat, name):
    """Require one polynomial matrix accepted by the lift wrapper.

    Args:
        mat: Candidate Sage matrix.
        name: Human-readable input name.

    Returns:
        None. The function validates inputs only.

    Raises:
        TypeError: If ``mat`` is not matrix-like.
        ValueError: If ``mat`` is not over ``P``.
    """
    if not all(hasattr(mat, method) for method in ("base_ring", "nrows", "ncols")):
        raise TypeError(f"Singular lift {name} must be a matrix over P.")
    if mat.base_ring() != P:
        raise ValueError(f"Singular lift {name} must be a matrix over P.")


def membership_lift(generators, targets):
    """Return lift coefficients for target columns of polynomial matrices.

    Args:
        generators: Polynomial matrix whose columns generate a module.
        targets: Polynomial matrix whose columns should be lifted into that
            generated module.

    Returns:
        Polynomial coefficient matrix ``witness`` satisfying
        ``generators * witness == targets`` in Singular's module sense.

    Raises:
        TypeError: If either input is not matrix-like.
        ValueError: If rings or row counts are incompatible.
        SingularBackendError: If Singular ``lift`` fails.
    """
    _validate_lift_matrix(generators, "generators")
    _validate_lift_matrix(targets, "targets")
    if generators.nrows() != targets.nrows():
        raise ValueError("Singular lift matrices must have the same row count.")
    try:
        return _lift(generators, targets)
    except RuntimeError as exc:
        raise SingularBackendError("Singular lift failed for polynomial matrix inputs.") from exc
