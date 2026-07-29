"""Finite decoders built from supplement-oriented decoupling maps."""

from sage.all import GF, vector

from decoder_core import TorusShape
from decoder_core.finite_maps import finite_matrix_from_laurent, zero_vector
from decoder_core.validation import sample_bsc_error
from isomorphism import (
    R,
    choose_smallest_oblique_cell,
    construct_excitation_map,
    decouple_coarse_matrix,
    oblique_coarse_grain,
    periods_from_generators,
)
from isomorphism.css import split_css_blocks

from .standard_code import StandardCodeDecoder


class DecouplingDecoder:
    """Decode source-code X errors through a standard decoupled code."""

    def __init__(
        self,
        source_epsilon,
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
            source_epsilon: Coarse CSS excitation matrix over the Laurent ring.
            num_x_checks: Number of source X-check rows per coarse cell.
            num_qubits: Number of source data-qubit columns per coarse cell.
            Lx: Number of coarse periods in the first torus direction.
            Ly: Number of coarse periods in the second torus direction.
            verify_maps: Whether to verify Laurent solver witnesses.
            weights: Optional matching weights for standard-code toric edges.
        """
        self.shape = TorusShape(int(Lx), int(Ly))
        self.source_epsilon = source_epsilon
        self.num_x_checks = int(num_x_checks)
        self.num_qubits = int(num_qubits)
        if self.source_epsilon.nrows() != 2 * self.num_x_checks:
            raise ValueError("source_epsilon row count must equal 2 * num_x_checks.")
        if self.source_epsilon.ncols() != 2 * self.num_qubits:
            raise ValueError("source_epsilon column count must equal 2 * num_qubits.")

        _hz_source, self.hx_source = split_css_blocks(
            self.source_epsilon,
            self.num_x_checks,
            self.num_qubits,
        )
        self.decoupling_result = decouple_coarse_matrix(
            self.source_epsilon,
            num_x_checks=self.num_x_checks,
            num_qubits=self.num_qubits,
            verify=bool(verify_maps),
            compute_forward_maps=False,
        )
        self.inverse_maps = self.decoupling_result.inverse_maps
        self.hx_standard = self.decoupling_result.hx_standard

        # All finite vectors use component-major ordering.  The supplement's
        # inverse maps go from the standard chain back to the source chain.
        self.hx_source_finite = finite_matrix_from_laurent(
            self.hx_source,
            self.shape,
        )
        self.phi0_inverse_finite = finite_matrix_from_laurent(
            self.inverse_maps.phi0_inverse,
            self.shape,
        )
        self.phi0_finite = self.phi0_inverse_finite.inverse()
        self.phi1_inverse_finite = finite_matrix_from_laurent(
            self.inverse_maps.phi1_inverse,
            self.shape,
        )
        self.standard_decoder = StandardCodeDecoder(
            self.hx_standard,
            self.shape,
            weights=weights,
        )
        self.hx_standard_finite = self.standard_decoder.hx_standard_finite

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
            DecouplingDecoder: Decoder for the coarse-grained BB code.
        """
        f = R(f)
        g = R(g)
        shape = TorusShape(int(Lx), int(Ly))
        period_limit = (
            max(512, shape.Lx, shape.Ly)
            if max_period is None
            else int(max_period)
        )
        bb_epsilon = construct_excitation_map(f, g)
        periods = periods_from_generators(f, g, max_period=period_limit)
        cell = choose_smallest_oblique_cell(periods.vectors)
        source_epsilon = oblique_coarse_grain(bb_epsilon, *cell)
        decoder = cls.from_coarse_css(
            source_epsilon,
            num_x_checks=source_epsilon.nrows() // 2,
            num_qubits=source_epsilon.ncols() // 2,
            Lx=shape.Lx,
            Ly=shape.Ly,
            verify_maps=verify_maps,
            weights=weights,
        )
        decoder.f = f
        decoder.g = g
        decoder.max_period = period_limit
        decoder.bb_epsilon = bb_epsilon
        decoder.periods = periods
        decoder.cell = cell
        return decoder

    @classmethod
    def from_coarse_css(
        cls,
        coarse_epsilon,
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
            coarse_epsilon: Coarse CSS excitation matrix over the Laurent ring.
            num_x_checks: Number of source X-check rows per coarse cell.
            num_qubits: Number of source data-qubit columns per coarse cell.
            Lx: Number of coarse periods in the first torus direction.
            Ly: Number of coarse periods in the second torus direction.
            verify_maps: Whether to verify Laurent solver witnesses.
            weights: Optional matching weights for standard-code toric edges.

        Returns:
            DecouplingDecoder: Decoder for the supplied coarse CSS matrix.
        """
        return cls(
            coarse_epsilon,
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
        return self.hx_source_finite * self.correction_vector(correction)

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
        standard_syndrome = self.phi0_finite * source_syndrome
        standard_correction = self.standard_decoder.decode(standard_syndrome)
        source_correction = self.phi1_inverse_finite * standard_correction
        if verify and self.hx_source_finite * source_correction != source_syndrome:
            raise ValueError("Decoded correction does not match the source syndrome.")
        return source_correction

    def check_chain_relation(self):
        """Return whether the finite inverse chain relation holds."""
        return (
            self.hx_source_finite * self.phi1_inverse_finite
            == self.phi0_inverse_finite * self.hx_standard_finite
        )

    def decode_sample(self, p, *, verify=True):
        """Sample and decode one BSC error.

        Args:
            p: Physical error probability.
            verify: Whether to verify the returned correction's syndrome.

        Returns:
            dict: Error, syndrome, correction, residual, and success fields.
        """
        error = self.sample_error(p)
        syndrome = self.syndrome(error)
        correction = self.decode(syndrome, verify=verify)
        return {
            "error": error,
            "syndrome": syndrome,
            "correction": correction,
            "residual": self.syndrome(error + correction),
            "success": self.syndrome(correction) == syndrome,
        }


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
