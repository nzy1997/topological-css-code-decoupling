"""Residual logical-failure classification helpers."""

from dataclasses import dataclass


def classify_attempt(hx_matrix, classifier, true_error, syndrome, correction):
    """Return ``(logical_failure, decode_failure)`` for one decoder output.

    Args:
        hx_matrix: Finite X-check matrix used for syndrome validation.
        classifier: Logical-failure classifier for zero-syndrome residuals.
        true_error: Sampled physical error used as the benchmark ground truth.
        syndrome: GF(2) syndrome vector to decode or classify.
        correction: Candidate GF(2) correction vector.

    Returns:
        object: ``(logical_failure, decode_failure)`` for one decoder output.

    """
    # Keep syndrome normalization separate from the matching solve.
    if hx_matrix * correction != syndrome:
        return True, True
    residual = true_error + correction
    return classifier.is_logical_failure(residual), False


def classify_verified_attempt(
    hx_matrix,
    classifier,
    true_error,
    syndrome,
    correction,
    *,
    decoder_name="decoder",
):
    """Classify a decoder output that is expected to always match syndrome.

    Args:
        hx_matrix: Finite X-check matrix used for syndrome validation.
        classifier: Logical-failure classifier for zero-syndrome residuals.
        true_error: Sampled physical error used as the benchmark ground truth.
        syndrome: GF(2) syndrome vector to decode or classify.
        correction: Candidate GF(2) correction vector.
        decoder_name: Human-readable decoder name for diagnostics.

    Returns:
        vector: Classify a decoder output that is expected to always match syndrome.

    """
    # Keep syndrome normalization separate from the matching solve.
    if hx_matrix * correction != syndrome:
        raise ValueError(f"{decoder_name} returned a correction with mismatched syndrome.")
    residual = true_error + correction
    return classifier.is_logical_failure(residual), False


@dataclass
class LogicalFailureClassifier:
    """Classify residual vectors modulo a stabilizer column space."""

    hx_matrix: object
    stabilizer_matrix: object

    def __post_init__(self):
        """Cache checks cutting out ``im(stabilizer_matrix)``.

        Args:
            None.

        Returns:
            None: This function mutates local state or performs validation only.
        """
        # A residual is a stabilizer exactly when it is orthogonal to every
        # vector in the right kernel of the stabilizer-column span.  Caching the
        # basis makes repeated Monte Carlo classification cheap and consistent.
        self._stabilizer_checks = (
            self.stabilizer_matrix.transpose().right_kernel().basis_matrix()
        )

    def has_zero_syndrome(self, residual):
        """Return whether ``residual`` has zero measured syndrome.

        Args:
            residual: Zero-syndrome residual after adding error and correction.

        Returns:
            bool: True when the requested invariant or condition holds.
        """
        # Keep syndrome normalization separate from the matching solve.
        return (self.hx_matrix * residual).is_zero()

    def is_stabilizer_residual(self, residual):
        """Return whether ``residual`` lies in the stabilizer column space.

        Args:
            residual: Zero-syndrome residual after adding error and correction.

        Returns:
            bool: True when the requested invariant or condition holds.
        """
        # Keep syndrome normalization separate from the matching solve.
        return (self._stabilizer_checks * residual).is_zero()

    def is_logical_failure(self, residual):
        """Return whether a zero-syndrome residual is a nontrivial logical.

        Args:
            residual: Zero-syndrome residual after adding error and correction.

        Returns:
            bool: True when the requested invariant or condition holds.
        """
        # Keep syndrome normalization separate from the matching solve.
        if not self.has_zero_syndrome(residual):
            return True
        return not self.is_stabilizer_residual(residual)
