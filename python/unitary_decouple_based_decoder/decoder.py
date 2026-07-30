"""Finite decoders built from the paper's decoupling chain maps."""

from sage.all import GF, vector

from decoder_core import TorusShape
from decoder_core.finite_maps import finite_matrix_from_laurent, zero_vector
from decoder_core.logical import LogicalFailureClassifier, classify_attempt
from decoder_core.validation import sample_bsc_error
from decoupling import (
    R,
    dagger_matrix,
    choose_minimal_area_superlattice_basis,
    build_two_generator_css_excitation_map,
    solve_decoupling_unitary,
    coarse_grain_to_superlattice,
    find_anyon_preserving_superlattice_from_generators,
)
from decoupling.css import split_css_blocks

from .standard_code import StandardComplexDecoder


class UnitaryDecoupleBasedDecoder:
    """Decode source-code X errors through a standard decoupled code."""

    def __init__(
        self,
        input_excitation_map,
        *,
        num_x_checks,
        num_qubits,
        Lx,
        Ly,
        verify_maps=False,
        weights=None,
    ):
        """Construct finite inverse-direction maps from a coarse CSS matrix.

        Use :meth:`from_bb` or :meth:`from_coarse_css` for public
        construction.

        Args:
            input_excitation_map: Coarse CSS excitation map over the Laurent ring.
            num_x_checks: Number of source X-check rows per coarse cell.
            num_qubits: Number of source data-qubit columns per coarse cell.
            Lx: Number of coarse periods in the first torus direction.
            Ly: Number of coarse periods in the second torus direction.
            verify_maps: Whether to verify Laurent solver witnesses.
            weights: Optional matching weights for standard-code toric edges.
        """
        self.shape = TorusShape(int(Lx), int(Ly))
        self.input_excitation_map = input_excitation_map
        self.num_x_checks = int(num_x_checks)
        self.num_qubits = int(num_qubits)
        if self.input_excitation_map.nrows() != 2 * self.num_x_checks:
            raise ValueError(
                "input_excitation_map row count must equal 2 * num_x_checks."
            )
        if self.input_excitation_map.ncols() != 2 * self.num_qubits:
            raise ValueError(
                "input_excitation_map column count must equal 2 * num_qubits."
            )

        self.h_z_input, self.h_x_input = split_css_blocks(
            self.input_excitation_map,
            self.num_x_checks,
            self.num_qubits,
        )
        self.decoupling_result = solve_decoupling_unitary(
            self.input_excitation_map,
            num_x_checks=self.num_x_checks,
            num_qubits=self.num_qubits,
            verify=bool(verify_maps),
            compute_psi=False,
        )
        self.psi_inverse = self.decoupling_result.psi_inverse
        self.h_x_tilde = self.decoupling_result.h_x_tilde

        # All finite vectors use component-major ordering.  The paper's
        # inverse maps go from the standard chain back to the source chain.
        self.h_x_input_finite = finite_matrix_from_laurent(
            self.h_x_input,
            self.shape,
        )
        self._h_z_dagger_input_finite = None
        self._logical_failure_classifier = None
        self.psi_0_inverse_finite = finite_matrix_from_laurent(
            self.psi_inverse.psi_0_inverse,
            self.shape,
        )
        self.psi_0_finite = self.psi_0_inverse_finite.inverse()
        self.psi_1_inverse_finite = finite_matrix_from_laurent(
            self.psi_inverse.psi_1_inverse,
            self.shape,
        )
        self.standard_decoder = StandardComplexDecoder(
            self.h_x_tilde,
            self.shape,
            weights=weights,
        )
        self.h_x_tilde_finite = self.standard_decoder.h_x_tilde_finite

    @classmethod
    def from_bb(
        cls,
        f,
        g,
        Lx,
        Ly,
        *,
        max_period=None,
        verify_maps=False,
        weights=None,
    ):
        """Construct a decoder for a BB code.

        Args:
            f: First BB-code Laurent generator.
            g: Second BB-code Laurent generator.
            Lx: Number of coarse periods in the first torus direction.
            Ly: Number of coarse periods in the second torus direction.
            max_period: Optional search budget for period detection.
            verify_maps: Whether to verify Laurent solver witnesses.
            weights: Optional matching weights for standard-code toric edges.

        Returns:
            UnitaryDecoupleBasedDecoder: Decoder for the coarse-grained BB code.
        """
        f = R(f)
        g = R(g)
        shape = TorusShape(int(Lx), int(Ly))
        period_limit = (
            max(512, shape.Lx, shape.Ly)
            if max_period is None
            else int(max_period)
        )
        bb_excitation_map = build_two_generator_css_excitation_map(f, g)
        periods = find_anyon_preserving_superlattice_from_generators(f, g, max_period=period_limit)
        superlattice_basis = choose_minimal_area_superlattice_basis(
            periods.period_vectors
        )
        input_excitation_map = coarse_grain_to_superlattice(
            bb_excitation_map,
            *superlattice_basis,
        )
        decoder = cls.from_coarse_css(
            input_excitation_map,
            num_x_checks=input_excitation_map.nrows() // 2,
            num_qubits=input_excitation_map.ncols() // 2,
            Lx=shape.Lx,
            Ly=shape.Ly,
            verify_maps=verify_maps,
            weights=weights,
        )
        decoder.f = f
        decoder.g = g
        decoder.max_period = period_limit
        decoder.bb_excitation_map = bb_excitation_map
        decoder.anyon_preserving_superlattice = periods
        decoder.superlattice_basis = superlattice_basis
        return decoder

    @classmethod
    def from_coarse_css(
        cls,
        input_excitation_map,
        *,
        num_x_checks,
        num_qubits,
        Lx,
        Ly,
        verify_maps=False,
        weights=None,
    ):
        """Construct a decoder from a pre-coarse-grained CSS matrix.

        Args:
            input_excitation_map: Coarse CSS excitation map over the Laurent ring.
            num_x_checks: Number of source X-check rows per coarse cell.
            num_qubits: Number of source data-qubit columns per coarse cell.
            Lx: Number of coarse periods in the first torus direction.
            Ly: Number of coarse periods in the second torus direction.
            verify_maps: Whether to verify Laurent solver witnesses.
            weights: Optional matching weights for standard-code toric edges.

        Returns:
            UnitaryDecoupleBasedDecoder: Decoder for the supplied coarse CSS matrix.
        """
        return cls(
            input_excitation_map,
            num_x_checks=num_x_checks,
            num_qubits=num_qubits,
            Lx=Lx,
            Ly=Ly,
            verify_maps=verify_maps,
            weights=weights,
        )

    @property
    def syndrome_size(self):
        """Return the finite source X-syndrome length."""
        return self.num_x_checks * self.shape.size

    @property
    def correction_size(self):
        """Return the finite source X-error length."""
        return self.num_qubits * self.shape.size

    @property
    def h_z_dagger_input_finite(self):
        """Return the cached finite source ``H_Z^dagger`` matrix."""
        if self._h_z_dagger_input_finite is None:
            self._h_z_dagger_input_finite = finite_matrix_from_laurent(
                dagger_matrix(self.h_z_input),
                self.shape,
            )
        return self._h_z_dagger_input_finite

    @property
    def logical_failure_classifier(self):
        """Return the cached source-code logical-failure classifier."""
        if self._logical_failure_classifier is None:
            self._logical_failure_classifier = LogicalFailureClassifier(
                self.h_x_input_finite,
                self.h_z_dagger_input_finite,
            )
        return self._logical_failure_classifier

    def zero_syndrome(self):
        """Return the zero finite source syndrome."""
        return zero_vector(self.num_x_checks, self.shape)

    def zero_correction(self):
        """Return the zero finite source correction."""
        return zero_vector(self.num_qubits, self.shape)

    def syndrome_vector(self, syndrome):
        """Normalize a source syndrome-like object to a GF(2) vector."""
        return _binary_vector(syndrome, self.syndrome_size, "syndrome")

    def correction_vector(self, correction):
        """Normalize a source correction-like object to a GF(2) vector."""
        return _binary_vector(correction, self.correction_size, "correction")

    def syndrome(self, correction):
        """Return the source X syndrome of a finite X-error vector."""
        return self.h_x_input_finite * self.correction_vector(correction)

    def sample_error(self, p):
        """Sample an independent BSC X-error."""
        return sample_bsc_error(self.correction_size, p)

    def decode(self, syndrome, *, verify=False):
        """Decode one finite source X syndrome.

        Args:
            syndrome: Source GF(2) syndrome vector.
            verify: Whether to verify the returned correction's syndrome.

        Returns:
            vector: Source GF(2) correction vector.
        """
        source_syndrome = self.syndrome_vector(syndrome)
        standard_syndrome = self.psi_0_finite * source_syndrome
        standard_correction = self.standard_decoder.decode(standard_syndrome)
        source_correction = self.psi_1_inverse_finite * standard_correction
        if verify and self.h_x_input_finite * source_correction != source_syndrome:
            raise ValueError("Decoded correction does not match the source syndrome.")
        return source_correction

    def check_chain_relation(self):
        """Return whether the finite inverse chain relation holds."""
        return (
            self.h_x_input_finite * self.psi_1_inverse_finite
            == self.psi_0_inverse_finite * self.h_x_tilde_finite
        )

    def decode_error(self, error, *, verify=True):
        """Decode and classify one supplied finite source X error.

        Args:
            error: Source GF(2) X-error vector.
            verify: Whether to verify the returned correction's syndrome.

        Returns:
            dict: Error, syndrome, correction, physical residual, classification,
            and success fields.
        """
        source_error = self.correction_vector(error)
        syndrome = self.syndrome(source_error)
        correction = self.decode(syndrome, verify=verify)
        residual = source_error + correction
        logical_failure, decode_failure = classify_attempt(
            self.h_x_input_finite,
            self.logical_failure_classifier,
            source_error,
            syndrome,
            correction,
        )
        return {
            "error": source_error,
            "syndrome": syndrome,
            "correction": correction,
            "residual": residual,
            "residual_syndrome": self.syndrome(residual),
            "logical_failure": logical_failure,
            "decode_failure": decode_failure,
            "success": not logical_failure and not decode_failure,
        }

    def decode_sample(self, p, *, verify=True):
        """Sample and decode one BSC error.

        Args:
            p: Physical error probability.
            verify: Whether to verify the returned correction's syndrome.

        Returns:
            dict: Error, syndrome, correction, physical residual,
            classification, and success fields.
        """
        return self.decode_error(self.sample_error(p), verify=verify)


def _binary_vector(values, expected, name):
    """Return a length-checked GF(2) vector."""
    if len(values) != expected:
        raise ValueError(f"{name} length must be {expected}.")
    if _is_gf2_vector(values):
        return values
    normalized = []
    for value in values:
        integer = int(value)
        if integer not in (0, 1):
            raise ValueError(f"{name} entries must be binary.")
        normalized.append(integer)
    return vector(GF(2), normalized)


def _is_gf2_vector(values):
    """Return whether ``values`` is already a Sage GF(2) vector."""
    try:
        return values.base_ring() == GF(2)
    except AttributeError:
        return False
